"""Pure-function, runtime-contract, and real-data tests for sequence numbering."""

from __future__ import annotations

import copy
import csv
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.sequence_numbering import (  # noqa: E402
    ANARCII_PARAMETERS,
    ANARCII_VERSION,
    REQUIRED_INPUT_FIELDS,
    InputSequence,
    SequenceNumberingError,
    audit_numbering_result,
    build_numbering_audits,
    imgt_region,
    sha256_bytes,
    validate_sequence_rows,
)
from antibody_optimization.sequence_numbering_artifacts import (  # noqa: E402
    INPUT_TABLE_NAME,
    POSITION_FIELDS,
    SEQUENCE_REVIEW_FIELDS,
    numbering_position_rows,
    validate_manifest_document,
)
from antibody_optimization import sequence_numbering_runtime  # noqa: E402


RECORDS = ROOT / "docs" / "result_artifacts" / "nb_expression" / INPUT_TABLE_NAME
INPUT_MANIFEST = ROOT / "docs" / "result_artifacts" / "nb_expression" / "manifest.json"
SCRIPT = ROOT / "scripts" / "input_baseline" / "build_sequence_review.py"
FIXED_GENERATED_AT = "2026-08-03T15:30:00+08:00"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def synthetic_record(
    sequence: str = "MAQACDGS", *, sample_uid: str = "TEST_ONLY__VHH"
) -> InputSequence:
    return InputSequence(
        sample_uid=sample_uid,
        provider_code="TEST_ONLY",
        source_sample_id="VHH",
        sequence_raw=sequence,
        sequence_length_aa=len(sequence),
        sequence_sha256=sha256_bytes(sequence.encode("ascii")),
        source_sequence_scope="unknown",
        source_vhh_region_sequence="DO_NOT_WRITE_BACK",
        source_sequence_review_flags="test_only_fixture",
    )


def synthetic_success_result() -> dict[str, object]:
    return {
        "numbering": [
            ((1, " "), "Q"),
            ((2, " "), "-"),
            ((3, " "), "A"),
            ((3, "A"), "C"),
            ((4, " "), "D"),
        ],
        "chain_type": "H",
        "score": 20.0,
        "query_start": 2,
        "query_end": 5,
        "error": None,
        "scheme": "imgt",
    }


class SequenceNumberingPureFunctionTests(unittest.TestCase):
    def test_validate_sequence_rows_checks_literal_length_and_hash_without_mutation(self) -> None:
        sequence = "QVQL"
        row = {
            "sample_uid": "TEST_ONLY__1",
            "provider_code": "TEST_ONLY",
            "source_sample_id": "1",
            "sequence_raw": sequence,
            "sequence_length_aa": str(len(sequence)),
            "sequence_sha256": sha256_bytes(sequence.encode("ascii")),
            "sequence_scope": "unknown",
            "vhh_region_sequence": "SOURCE_SENTINEL",
            "sequence_review_flags": "",
        }
        before = copy.deepcopy(row)
        records = validate_sequence_rows(
            REQUIRED_INPUT_FIELDS, [row], expected_count=1
        )
        self.assertEqual(records[0].source_vhh_region_sequence, "SOURCE_SENTINEL")
        self.assertEqual(row, before)

        bad_length = {**row, "sequence_length_aa": "5"}
        with self.assertRaisesRegex(SequenceNumberingError, "length mismatch"):
            validate_sequence_rows(REQUIRED_INPUT_FIELDS, [bad_length], expected_count=1)
        bad_hash = {**row, "sequence_sha256": "0" * 64}
        with self.assertRaisesRegex(SequenceNumberingError, "SHA-256 mismatch"):
            validate_sequence_rows(REQUIRED_INPUT_FIELDS, [bad_hash], expected_count=1)

    def test_manifest_validation_binds_csv_bytes_count_and_header(self) -> None:
        fieldnames = list(REQUIRED_INPUT_FIELDS)
        manifest = {
            "record_count": 1,
            "outputs": {
                INPUT_TABLE_NAME: {"sha256": "a" * 64, "size_bytes": 123}
            },
            "tables": {INPUT_TABLE_NAME: fieldnames},
            "sequence_integrity": {"per_record_sha256": True},
            "dataset": "test_only",
            "parser_version": "test_only",
        }
        validated = validate_manifest_document(
            manifest,
            records_size_bytes=123,
            records_sha256="a" * 64,
            csv_fieldnames=fieldnames,
            expected_count=1,
        )
        self.assertTrue(validated["per_record_sha256"])
        bad = copy.deepcopy(manifest)
        bad["outputs"][INPUT_TABLE_NAME]["sha256"] = "b" * 64
        with self.assertRaisesRegex(SequenceNumberingError, "SHA-256"):
            validate_manifest_document(
                bad,
                records_size_bytes=123,
                records_sha256="a" * 64,
                csv_fieldnames=fieldnames,
                expected_count=1,
            )

    def test_success_mapping_uses_inclusive_bounds_and_preserves_gaps_insertions(self) -> None:
        audit = audit_numbering_result(
            synthetic_record(), synthetic_success_result()
        )
        self.assertEqual(audit.numbering_status, "pass")
        self.assertEqual(audit.sequence_scope_status, "provisional_numbered_domain")
        self.assertEqual(audit.provisional_numbered_span_sequence, "QACD")
        self.assertEqual(
            audit.provisional_numbered_span_sha256,
            sha256_bytes(b"QACD"),
        )
        self.assertEqual(
            [position.sequence_index_0based for position in audit.positions],
            [2, None, 3, 4, 5],
        )
        self.assertEqual([position.label for position in audit.positions], ["1", "2", "3", "3A", "4"])
        rows = numbering_position_rows([audit])
        self.assertEqual(
            [row["sequence_index_1based"] for row in rows],
            [3, "", 4, 5, 6],
        )
        self.assertNotIn("vhh_region_sequence", SEQUENCE_REVIEW_FIELDS)
        self.assertEqual(tuple(rows[0]), POSITION_FIELDS)

    def test_reconstruction_mismatch_and_multidomain_like_keys_are_rejected(self) -> None:
        raw = synthetic_success_result()
        raw["query_end"] = 6
        with self.assertRaisesRegex(SequenceNumberingError, "reconstruct"):
            audit_numbering_result(synthetic_record(), raw)

        record = synthetic_record()
        with self.assertRaisesRegex(SequenceNumberingError, "identities"):
            build_numbering_audits(
                [record],
                {
                    record.sample_uid: synthetic_success_result(),
                    f"{record.sample_uid}-1": synthetic_success_result(),
                },
            )

    def test_failed_result_has_no_invented_span_or_positions(self) -> None:
        audit = audit_numbering_result(
            synthetic_record(),
            {
                "numbering": None,
                "chain_type": "F",
                "score": 13.0,
                "query_start": None,
                "query_end": None,
                "error": "Score less than cut off.",
                "scheme": "imgt",
            },
        )
        self.assertEqual(audit.numbering_status, "failed")
        self.assertEqual(audit.sequence_scope_status, "unresolved")
        self.assertEqual(audit.provisional_numbered_span_sha256, "")
        self.assertEqual(audit.positions, ())

    def test_imgt_region_boundaries_are_explicit(self) -> None:
        expected = {
            1: "FR1",
            26: "FR1",
            27: "CDR1",
            38: "CDR1",
            39: "FR2",
            55: "FR2",
            56: "CDR2",
            65: "CDR2",
            66: "FR3",
            104: "FR3",
            105: "CDR3",
            117: "CDR3",
            118: "FR4",
            128: "FR4",
        }
        self.assertEqual({position: imgt_region(position) for position in expected}, expected)
        with self.assertRaises(SequenceNumberingError):
            imgt_region(129)


class SequenceNumberingRuntimeTests(unittest.TestCase):
    def test_runtime_uses_only_the_pinned_anarcii_configuration(self) -> None:
        record = synthetic_record(sequence="Q", sample_uid="TEST_ONLY__RUNTIME")
        calls: dict[str, object] = {}

        class FakeAnarcii:
            def __init__(self, **kwargs: object) -> None:
                calls["constructor"] = kwargs

            def number(self, values: object, *, scfv: bool) -> object:
                calls["number_values"] = values
                calls["scfv"] = scfv
                return None

            def to_scheme(self, *, scheme: str) -> dict[str, object]:
                calls["scheme"] = scheme
                return {
                    record.sample_uid: {
                        "numbering": [((1, " "), "Q")],
                        "chain_type": "H",
                        "score": 20.0,
                        "query_start": 0,
                        "query_end": 0,
                        "error": None,
                        "scheme": "imgt",
                    }
                }

        fake_module = types.ModuleType("anarcii")
        fake_module.Anarcii = FakeAnarcii
        with mock.patch.object(
            sequence_numbering_runtime.importlib_metadata,
            "version",
            return_value=ANARCII_VERSION,
        ), mock.patch.dict(sys.modules, {"anarcii": fake_module}):
            audits = sequence_numbering_runtime.run_anarcii_numbering([record])

        self.assertEqual(
            calls["constructor"],
            {
                "seq_type": "antibody",
                "mode": "accuracy",
                "cpu": True,
                "ncpu": 1,
                "batch_size": 8,
            },
        )
        self.assertIs(calls["scfv"], False)
        self.assertEqual(calls["scheme"], "imgt")
        self.assertEqual(audits[0].numbering_status, "pass")
        self.assertEqual(
            ANARCII_PARAMETERS,
            {
                "seq_type": "antibody",
                "mode": "accuracy",
                "scheme": "imgt",
                "cpu": True,
                "ncpu": 1,
                "batch_size": 8,
                "scfv": False,
            },
        )

    def test_runtime_rejects_any_other_anarcii_version(self) -> None:
        with mock.patch.object(
            sequence_numbering_runtime.importlib_metadata,
            "version",
            return_value="2.0.9",
        ):
            with self.assertRaisesRegex(RuntimeError, "2.0.8 is required"):
                sequence_numbering_runtime.run_anarcii_numbering([synthetic_record()])


class SequenceNumberingRealDataTests(unittest.TestCase):
    def run_command(self, *arguments: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, *arguments],
            cwd=ROOT,
            capture_output=True,
            text=True,
            encoding="utf-8",
            env={**os.environ, "PYTHONUTF8": "1"},
            check=False,
        )

    def arguments(self, output_dir: Path, run_summary: Path) -> list[str]:
        return [
            str(SCRIPT),
            "--records",
            str(RECORDS),
            "--input-manifest",
            str(INPUT_MANIFEST),
            "--output-dir",
            str(output_dir),
            "--run-summary",
            str(run_summary),
            "--generated-at",
            FIXED_GENERATED_AT,
        ]

    def test_real_47_sequence_run_acceptance_bom_and_overwrite_safety(self) -> None:
        records_before = digest(RECORDS)
        input_manifest_before = digest(INPUT_MANIFEST)
        with tempfile.TemporaryDirectory(
            prefix=".sequence-numbering-test-", dir=ROOT
        ) as temporary_dir:
            base = Path(temporary_dir)
            output_dir = base / "artifacts"
            run_summary_path = base / "sequence_numbering.json"
            arguments = self.arguments(output_dir, run_summary_path)

            first = self.run_command(*arguments)
            self.assertEqual(first.returncode, 0, first.stderr)
            review_path = output_dir / "sequence_numbering_review.csv"
            positions_path = output_dir / "sequence_numbering_positions.csv"
            artifact_manifest_path = output_dir / "sequence_numbering_manifest.json"
            self.assertEqual(review_path.read_bytes()[:3], b"\xef\xbb\xbf")
            self.assertEqual(positions_path.read_bytes()[:3], b"\xef\xbb\xbf")

            with review_path.open("r", encoding="utf-8-sig", newline="") as handle:
                reader = csv.DictReader(handle)
                self.assertEqual(tuple(reader.fieldnames or ()), SEQUENCE_REVIEW_FIELDS)
                review_rows = list(reader)
            self.assertEqual(len(review_rows), 47)
            statuses = {status: 0 for status in ("pass", "failed")}
            for row in review_rows:
                statuses[row["numbering_status"]] += 1
            self.assertEqual(statuses, {"pass": 46, "failed": 1})

            by_uid = {row["sample_uid"]: row for row in review_rows}
            nb252 = by_uid["LTT__Nb252"]
            self.assertEqual(nb252["numbering_status"], "pass")
            self.assertEqual(nb252["sequence_scope_status"], "provisional_numbered_domain")
            self.assertEqual(nb252["chain_type"], "H")
            self.assertEqual(nb252["scheme"], "imgt")
            self.assertEqual(nb252["query_start_0based_inclusive"], "0")
            self.assertEqual(nb252["query_end_0based_inclusive"], "125")
            self.assertEqual(nb252["numbered_non_gap_count"], "126")
            self.assertEqual(nb252["first_numbered_imgt_position"], "1")
            self.assertEqual(nb252["last_numbered_imgt_position"], "128")
            self.assertEqual(nb252["unnumbered_c_sequence"], "GS")
            self.assertTrue(nb252["score"])

            wcc_4_11 = by_uid["WCC__4-11"]
            self.assertEqual(wcc_4_11["chain_type"], "L")
            self.assertEqual(wcc_4_11["last_numbered_imgt_position"], "117")
            self.assertIn(
                "anarcii_chain_type_L", wcc_4_11["numbering_review_flags"]
            )

            wcc_4_1 = by_uid["WCC__4-1"]
            self.assertEqual(wcc_4_1["query_start_0based_inclusive"], "2")
            self.assertEqual(wcc_4_1["query_end_0based_inclusive"], "88")
            self.assertEqual(wcc_4_1["numbered_non_gap_count"], "87")
            self.assertEqual(wcc_4_1["last_numbered_imgt_position"], "102")

            failure = by_uid["WCC__4-28"]
            self.assertEqual(failure["numbering_status"], "failed")
            self.assertEqual(failure["sequence_scope_status"], "unresolved")
            self.assertEqual(failure["error"], "Score less than cut off.")
            self.assertEqual(failure["provisional_numbered_span_sha256"], "")

            with positions_path.open("r", encoding="utf-8-sig", newline="") as handle:
                reader = csv.DictReader(handle)
                self.assertEqual(tuple(reader.fieldnames or ()), POSITION_FIELDS)
                position_rows = list(reader)
            self.assertTrue(position_rows)
            self.assertEqual({row["numbering_status"] for row in position_rows}, {"pass"})
            self.assertFalse(any(row["sample_uid"] == "WCC__4-28" for row in position_rows))
            nb252_non_gap = [
                row
                for row in position_rows
                if row["sample_uid"] == "LTT__Nb252" and row["is_gap"] == "false"
            ]
            self.assertEqual(len(nb252_non_gap), 126)
            self.assertEqual(
                [int(row["sequence_index_1based"]) for row in nb252_non_gap],
                list(range(1, 127)),
            )
            self.assertEqual(nb252_non_gap[-1]["numbering_position_label"], "128")

            artifact_manifest = json.loads(
                artifact_manifest_path.read_text(encoding="utf-8")
            )
            run_summary = json.loads(run_summary_path.read_text(encoding="utf-8"))
            for document in (artifact_manifest, run_summary):
                self.assertEqual(len(document["samples"]), 47)
                sample = document["samples"][0]
                self.assertTrue(
                    {
                        "sample_uid",
                        "numbering_status",
                        "sequence_scope_status",
                        "provisional_numbered_span_sha256",
                    }.issubset(sample)
                )
                self.assertEqual(
                    document["tool"]["parameters"], ANARCII_PARAMETERS
                )
                self.assertEqual(
                    document["acceptance"]["score_assertion_policy"],
                    "scores recorded but no exact score asserted",
                )

            hashes_before_refusal = {
                path.name: digest(path)
                for path in (review_path, positions_path, artifact_manifest_path, run_summary_path)
            }
            refused = self.run_command(*arguments)
            self.assertNotEqual(refused.returncode, 0)
            self.assertIn("Refusing to overwrite", refused.stderr)
            self.assertEqual(
                {
                    path.name: digest(path)
                    for path in (review_path, positions_path, artifact_manifest_path, run_summary_path)
                },
                hashes_before_refusal,
            )

            overwritten = self.run_command(*arguments, "--overwrite")
            self.assertEqual(overwritten.returncode, 0, overwritten.stderr)

        self.assertEqual(digest(RECORDS), records_before)
        self.assertEqual(digest(INPUT_MANIFEST), input_manifest_before)


if __name__ == "__main__":
    unittest.main()
