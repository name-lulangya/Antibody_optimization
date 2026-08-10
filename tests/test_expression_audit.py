"""Tests for the conservative expression-data comparability audit."""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.expression_audit import (  # noqa: E402
    ASSAY_METADATA_FIELDS,
    EXPECTED_SEMANTICS_COUNTS,
    EXPECTED_SOURCE_COUNTS,
    SAMPLE_REVIEW_FIELDS,
    VIEW_FIELDS,
    ExpressionAuditError,
    build_assay_metadata_rows,
    build_allowed_use_manifest,
    build_comparability_view,
    build_sample_review_rows,
    load_and_validate_inputs,
    load_cross_provider_confirmation,
    sha256_file,
    validate_written_audit,
    write_csv,
)


ARTIFACTS = ROOT / "docs" / "result_artifacts" / "nb_expression"
RECORDS = ARTIFACTS / "nb_expression_records.csv"
CONTEXT = ARTIFACTS / "assay_context.csv"
MANIFEST = ARTIFACTS / "manifest.json"
CONFIRMATION = (
    ROOT
    / "docs"
    / "result_artifacts"
    / "input_baseline"
    / "reviews"
    / "expression_cross_provider_confirmation.json"
)
SCRIPT = ROOT / "scripts" / "input_baseline" / "build_expression_audit.py"
FIXED_TIME = "2026-08-03T18:00:00+08:00"


class ExpressionAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.inputs = load_and_validate_inputs(RECORDS, CONTEXT, MANIFEST)
        cls.metadata = build_assay_metadata_rows(
            cls.inputs.assay_contexts, generated_at=FIXED_TIME
        )
        cls.samples = build_sample_review_rows(cls.inputs)
        cls.view = build_comparability_view(cls.samples, cls.metadata)

    def test_real_inputs_have_frozen_counts_and_semantics(self) -> None:
        self.assertEqual(
            Counter(row["provider_code"] for row in self.inputs.records),
            Counter(EXPECTED_SOURCE_COUNTS),
        )
        self.assertEqual(
            Counter(row["observation_semantics"] for row in self.inputs.records),
            Counter(EXPECTED_SEMANTICS_COUNTS),
        )
        self.assertEqual(len({row["sample_uid"] for row in self.inputs.records}), 47)
        self.assertTrue(all(row["sequence_scope"] == "unknown" for row in self.inputs.records))

    def test_metadata_does_not_promote_same_system_claim(self) -> None:
        indexed = {
            (row["provider_code"], row["field_name"]): row for row in self.metadata
        }
        for provider in EXPECTED_SOURCE_COUNTS:
            claim = indexed[(provider, "system_equivalence_claim")]
            self.assertEqual(claim["field_value"], "same_system")
            self.assertEqual(claim["evidence_status"], "user_provided")
            self.assertEqual(claim["reviewed_at"], FIXED_TIME)
            self.assertEqual(
                indexed[(provider, "cross_provider_protocol_equivalence")][
                    "evidence_status"
                ],
                "unknown_not_reported",
            )
        self.assertEqual(indexed[("WCC", "medium")]["field_value"], "TB")
        self.assertEqual(
            indexed[("WCC", "yield_stage")]["field_value"], "post_purification"
        )
        for provider in ("LTT", "LLJ"):
            self.assertEqual(indexed[(provider, "medium")]["field_value"], "")
            self.assertEqual(indexed[(provider, "yield_stage")]["field_value"], "")

    def test_sample_decisions_are_conservative_and_semantics_specific(self) -> None:
        self.assertEqual(len(self.samples), 47)
        self.assertTrue(
            all(row["cross_assay_pooling_status"] == "blocked" for row in self.samples)
        )
        self.assertTrue(all(row["nb252_transfer_status"] == "blocked" for row in self.samples))
        self.assertTrue(all(row["review_version"] == "1.1.0" for row in self.samples))
        individual = [row for row in self.samples if row["provider_code"] in {"LTT", "WCC"}]
        llj = [row for row in self.samples if row["provider_code"] == "LLJ"]
        self.assertEqual(len(individual), 31)
        self.assertTrue(
            all(row["within_assay_numeric_use_status"] == "conditional" for row in individual)
        )
        self.assertTrue(
            all(row["within_assay_ordinal_use_status"] == "conditional" for row in llj)
        )
        self.assertTrue(all(row["point_estimate_mg"] == "" for row in llj))

    def test_allowed_use_manifest_keeps_release_gates_blocked(self) -> None:
        manifest = build_allowed_use_manifest(
            inputs=self.inputs,
            metadata_rows=self.metadata,
            sample_rows=self.samples,
            view_rows=self.view,
            generated_at=FIXED_TIME,
            output_hashes={},
        )
        self.assertEqual(manifest["gates"]["expression_audit_gate"], "pass")
        self.assertEqual(manifest["gates"]["cross_assay_pooling_gate"], "blocked")
        self.assertEqual(manifest["gates"]["nb252_transfer_gate"], "blocked")
        self.assertEqual(
            manifest["sequence_audit_summary"]["numbering_status_counts"],
            {"pending": 47},
        )

    def test_collaborator_confirmation_passes_semantics_aware_pooling_only(self) -> None:
        confirmation = load_cross_provider_confirmation(CONFIRMATION)
        metadata = build_assay_metadata_rows(
            self.inputs.assay_contexts,
            generated_at=FIXED_TIME,
            comparability_confirmation=confirmation,
        )
        samples = build_sample_review_rows(
            self.inputs, comparability_confirmation=confirmation
        )
        view = build_comparability_view(samples, metadata)
        manifest = build_allowed_use_manifest(
            inputs=self.inputs,
            metadata_rows=metadata,
            sample_rows=samples,
            view_rows=view,
            generated_at=FIXED_TIME,
            output_hashes={},
            comparability_confirmation=confirmation,
        )

        self.assertTrue(
            all(row["cross_assay_pooling_status"] == "pass" for row in samples)
        )
        self.assertTrue(all(row["nb252_transfer_status"] == "blocked" for row in samples))
        self.assertEqual(manifest["gates"]["cross_assay_pooling_gate"], "pass")
        self.assertEqual(manifest["gates"]["nb252_transfer_gate"], "blocked")
        self.assertEqual(
            manifest["comparability_confirmation"]["status"], "confirmed"
        )
        self.assertIn("pooled_continuous_yield_model", manifest["blocked_uses"])
        self.assertNotIn("cross_assay_pooling", manifest["blocked_uses"])

    def test_optional_sequence_summary_is_joined_but_not_promoted(self) -> None:
        sample = self.inputs.records[0]
        with tempfile.TemporaryDirectory() as temporary_dir:
            path = Path(temporary_dir) / "sequence_summary.json"
            path.write_text(
                json.dumps(
                    {
                        "samples": [
                            {
                                "sample_uid": sample["sample_uid"],
                                "numbering_status": "pass",
                                "sequence_scope_status": "provisional_numbered_domain",
                                "provisional_numbered_span_sha256": "a" * 64,
                            }
                        ]
                    }
                ),
                encoding="utf-8",
            )
            inputs = load_and_validate_inputs(RECORDS, CONTEXT, MANIFEST, path)
        rows = build_sample_review_rows(inputs)
        joined = next(row for row in rows if row["sample_uid"] == sample["sample_uid"])
        missing = next(row for row in rows if row["sample_uid"] != sample["sample_uid"])
        self.assertEqual(joined["numbering_status"], "pass")
        self.assertEqual(joined["provisional_numbered_span_sha256"], "a" * 64)
        self.assertEqual(joined["construct_comparability_status"], "pending")
        self.assertEqual(joined["nb252_transfer_status"], "blocked")
        self.assertEqual(missing["numbering_status"], "pending")
        self.assertEqual(missing["sequence_scope_status"], "unknown_not_reported")

    def test_written_csvs_have_bom_fixed_columns_and_validated_counts(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            paths = {
                "assay_metadata": root / "assay_metadata_review.csv",
                "sample_review": root / "sample_comparability_review.csv",
                "view": root / "expression_comparability_view.csv",
            }
            write_csv(paths["assay_metadata"], ASSAY_METADATA_FIELDS, self.metadata)
            write_csv(paths["sample_review"], SAMPLE_REVIEW_FIELDS, self.samples)
            write_csv(paths["view"], VIEW_FIELDS, self.view)
            report = validate_written_audit(paths)
            self.assertEqual(report["status"], "pass")
            self.assertEqual(report["sample_review_row_count"], 47)
            for path in paths.values():
                self.assertTrue(path.read_bytes().startswith(b"\xef\xbb\xbf"))
                with path.open("r", encoding="utf-8-sig", newline="") as handle:
                    self.assertIsNotNone(csv.DictReader(handle).fieldnames)

    def test_manifest_hash_check_rejects_modified_records(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            altered = Path(temporary_dir) / "nb_expression_records.csv"
            shutil.copyfile(RECORDS, altered)
            with altered.open("ab") as handle:
                handle.write(b"\n")
            with self.assertRaisesRegex(ExpressionAuditError, "Input hash mismatch"):
                load_and_validate_inputs(altered, CONTEXT, MANIFEST)

    def test_cli_is_transactional_and_refuses_default_overwrite(self) -> None:
        input_hashes = {
            path: sha256_file(path) for path in (RECORDS, CONTEXT, MANIFEST)
        }
        with tempfile.TemporaryDirectory(prefix=".expression-audit-test-", dir=ROOT) as temporary_dir:
            root = Path(temporary_dir)
            output = root / "artifacts"
            summary = root / "summary.json"
            command = [
                sys.executable,
                str(SCRIPT),
                "--records",
                str(RECORDS),
                "--assay-context",
                str(CONTEXT),
                "--manifest",
                str(MANIFEST),
                "--output-dir",
                str(output),
                "--run-summary",
                str(summary),
                "--generated-at",
                FIXED_TIME,
            ]
            result = subprocess.run(
                command,
                cwd=ROOT,
                capture_output=True,
                text=True,
                encoding="utf-8",
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(summary.is_file())
            run_summary = json.loads(summary.read_text(encoding="utf-8"))
            self.assertEqual(run_summary["validation"]["sample_review_row_count"], 47)
            rerun = subprocess.run(
                command,
                cwd=ROOT,
                capture_output=True,
                text=True,
                encoding="utf-8",
                check=False,
            )
            self.assertNotEqual(rerun.returncode, 0)
            self.assertIn("Refusing to overwrite", rerun.stderr)
        self.assertEqual(
            {path: sha256_file(path) for path in input_hashes}, input_hashes
        )


if __name__ == "__main__":
    unittest.main()
