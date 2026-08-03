from __future__ import annotations

import csv
import json
import sys
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from antibody_optimization.input_integrity import (  # noqa: E402
    EXPECTED_CXS_SHA256,
    InputIntegrityError,
    build_input_freeze_manifest,
)


CXS = PROJECT_ROOT / "Nb252-optimization.cxs"
DOCX = PROJECT_ROOT / "nb序列及产量（1L）.docx"
EXPRESSION_MANIFEST_PATH = (
    PROJECT_ROOT / "docs/result_artifacts/nb_expression/manifest.json"
)
RECORDS_PATH = (
    PROJECT_ROOT / "docs/result_artifacts/nb_expression/nb_expression_records.csv"
)
SEQUENCE_MANIFEST_PATH = (
    PROJECT_ROOT
    / "docs/result_artifacts/input_baseline/sequence/sequence_numbering_manifest.json"
)
AUDIT_MANIFEST_PATH = (
    PROJECT_ROOT
    / "docs/result_artifacts/input_baseline/expression/allowed_use_manifest.json"
)
NUMBERING_REVIEW_PATH = (
    PROJECT_ROOT
    / "docs/result_artifacts/input_baseline/sequence/sequence_numbering_review.csv"
)
NUMBERING_POSITIONS_PATH = (
    PROJECT_ROOT
    / "docs/result_artifacts/input_baseline/sequence/sequence_numbering_positions.csv"
)
SAMPLE_COMPARABILITY_PATH = (
    PROJECT_ROOT
    / "docs/result_artifacts/input_baseline/expression/sample_comparability_review.csv"
)


def read_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    assert isinstance(value, dict)
    return value


def read_records() -> list[dict[str, str]]:
    with RECORDS_PATH.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


class InputIntegrityTests(unittest.TestCase):
    def build(self, records: list[dict[str, str]]) -> dict[str, object]:
        return build_input_freeze_manifest(
            source_cxs=CXS,
            source_docx=DOCX,
            expression_manifest_path=EXPRESSION_MANIFEST_PATH,
            expression_manifest=read_json(EXPRESSION_MANIFEST_PATH),
            expression_records_path=RECORDS_PATH,
            expression_records=records,
            sequence_manifest_path=SEQUENCE_MANIFEST_PATH,
            sequence_manifest=read_json(SEQUENCE_MANIFEST_PATH),
            numbering_review_path=NUMBERING_REVIEW_PATH,
            numbering_positions_path=NUMBERING_POSITIONS_PATH,
            expression_audit_manifest_path=AUDIT_MANIFEST_PATH,
            expression_audit_manifest=read_json(AUDIT_MANIFEST_PATH),
            sample_comparability_path=SAMPLE_COMPARABILITY_PATH,
            generated_at="2026-08-03T20:00:00+08:00",
            python_version="3.11.15",
            gemmi_version="0.7.5",
            expected_cxs_sha256=EXPECTED_CXS_SHA256.upper(),
        )

    def test_real_sources_and_all_47_literal_sequences_freeze(self) -> None:
        manifest = self.build(read_records())
        self.assertEqual(manifest["status"], "pass")
        self.assertEqual(
            manifest["sources"]["collaborator_cxs"]["sha256"],
            EXPECTED_CXS_SHA256,
        )
        self.assertEqual(manifest["sequence_identity"]["record_count"], 47)
        samples = manifest["sequence_identity"]["per_sample"]
        self.assertEqual(len(samples), 47)
        nb252 = next(row for row in samples if row["sample_uid"] == "LTT__Nb252")
        self.assertEqual(nb252["sequence_length_aa"], 128)
        self.assertEqual(
            nb252["sequence_sha256"],
            "df5b83ddde8a3486383c12afe45e22af6a358f507eab5503d5dbd4430710288d",
        )

    def test_literal_sequence_change_is_rejected(self) -> None:
        records = read_records()
        records[0] = dict(records[0])
        records[0]["sequence_raw"] = "A" + records[0]["sequence_raw"][1:]
        with self.assertRaisesRegex(InputIntegrityError, "Sequence SHA-256 mismatch"):
            self.build(records)

    def test_wrong_tool_contract_is_rejected(self) -> None:
        with self.assertRaisesRegex(InputIntegrityError, "Expected Gemmi"):
            build_input_freeze_manifest(
                source_cxs=CXS,
                source_docx=DOCX,
                expression_manifest_path=EXPRESSION_MANIFEST_PATH,
                expression_manifest=read_json(EXPRESSION_MANIFEST_PATH),
                expression_records_path=RECORDS_PATH,
                expression_records=read_records(),
                sequence_manifest_path=SEQUENCE_MANIFEST_PATH,
                sequence_manifest=read_json(SEQUENCE_MANIFEST_PATH),
                numbering_review_path=NUMBERING_REVIEW_PATH,
                numbering_positions_path=NUMBERING_POSITIONS_PATH,
                expression_audit_manifest_path=AUDIT_MANIFEST_PATH,
                expression_audit_manifest=read_json(AUDIT_MANIFEST_PATH),
                sample_comparability_path=SAMPLE_COMPARABILITY_PATH,
                generated_at="2026-08-03T20:00:00+08:00",
                python_version="3.11.15",
                gemmi_version="0.7.4",
            )


if __name__ == "__main__":
    unittest.main()
