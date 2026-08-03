"""Regression tests against the unchanged collaborator source document."""

from __future__ import annotations

import shutil
import sys
import tempfile
import unittest
from collections import Counter
from decimal import Decimal
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.nb_expression import (  # noqa: E402
    parse_expression_docx,
    validate_records,
)
from antibody_optimization.nb_expression_artifacts import (  # noqa: E402
    assay_context_rows,
    render_qc_svg,
    validate_written_outputs,
    write_assay_context_csv,
    write_fasta,
    write_qc_plot_data,
    write_raw_transcription_csv,
    write_samples_csv,
    write_wide_records_csv,
    write_yield_observations_csv,
)


SOURCE = ROOT / "nb序列及产量（1L）.docx"
SOURCE_SHA256 = "a6e4022f0978fbd70a0e04dc78f479140ab6f55caaa90b467fb77a62eb5db5d1"
EXPECTED_IDS = {
    "LTT": [
        "Nb290",
        "Nb270",
        "Nb314",
        "Nb5",
        "Nb257",
        "Nb232",
        "Nb294",
        "Nb256",
        "Nb233",
        "Nb302",
        "Nb253",
        "Nb21",
        "Nb252",
        "Q1",
        "P26",
        "Q17",
        "P27",
        "S53",
        "S31",
        "S4",
        "S19",
        "S73",
        "S26",
    ],
    "WCC": ["4-7", "4-34", "4-1", "4-28", "4-11", "4-40", "4-36", "4-42"],
    "LLJ": [
        "1G7",
        "2A4",
        "2C8",
        "1E2",
        "2E11",
        "1H7",
        "2C6",
        "2E2",
        "1G11",
        "1H6",
        "1F12",
        "2H8",
        "1B11",
        "2F12",
        "2G7",
        "1G4",
    ],
}


class ExpressionDocxTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.records, cls.paragraphs = parse_expression_docx(
            SOURCE, document_title_culture_volume_l=Decimal("1")
        )

    def test_source_contract_and_all_clone_ids(self) -> None:
        report = validate_records(
            self.records,
            self.paragraphs,
            expected_source_sha256=SOURCE_SHA256,
        )
        self.assertEqual(report["status"], "pass")
        self.assertEqual(len(self.records), 47)
        self.assertEqual(len(self.paragraphs), 100)
        for provider, expected in EXPECTED_IDS.items():
            observed = [
                record.source_sample_id
                for record in self.records
                if record.provider_code == provider
            ]
            self.assertEqual(observed, expected)

    def test_every_sequence_has_unique_verified_hash(self) -> None:
        self.assertEqual(len({record.sequence_raw for record in self.records}), 47)
        self.assertEqual(len({record.sequence_sha256 for record in self.records}), 47)
        for record in self.records:
            source_sequence = self.paragraphs[
                record.source_sequence_paragraph_index - 1
            ].raw_text
            self.assertEqual(record.sequence_raw, source_sequence)
            self.assertEqual(record.sequence_length_aa, len(source_sequence))

    def test_yield_semantics_remain_separate(self) -> None:
        semantics = Counter(record.observation_semantics for record in self.records)
        self.assertEqual(
            semantics,
            {
                "individual_approximate": 31,
                "group_lower_bound": 9,
                "group_approximate": 7,
            },
        )
        llj = [record for record in self.records if record.provider_code == "LLJ"]
        self.assertTrue(all(record.point_estimate_mg is None for record in llj))
        lower_bound = [r for r in llj if r.observation_semantics == "group_lower_bound"]
        self.assertEqual(len(lower_bound), 9)
        self.assertTrue(all(r.lower_bound_mg == Decimal("20") for r in lower_bound))
        target = next(record for record in self.records if record.sample_uid == "LTT__Nb252")
        self.assertEqual(target.point_estimate_mg, Decimal("0.5"))
        self.assertEqual(target.sequence_sha256[:12], "df5b83ddde8a")

    def test_all_written_sequences_round_trip(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            output = Path(temporary_dir)
            paths = {
                "samples": output / "samples.csv",
                "yields": output / "yield_observations.csv",
                "assay_context": output / "assay_context.csv",
                "wide": output / "nb_expression_records.csv",
                "raw_transcription": output / "raw_transcription.csv",
                "fasta": output / "nb_expression_sequences.fasta",
            }
            write_samples_csv(paths["samples"], self.records)
            write_yield_observations_csv(paths["yields"], self.records)
            write_assay_context_csv(paths["assay_context"], self.records)
            write_wide_records_csv(paths["wide"], self.records)
            write_raw_transcription_csv(paths["raw_transcription"], self.paragraphs)
            write_fasta(paths["fasta"], self.records)
            report = validate_written_outputs(paths, self.records, self.paragraphs)
            self.assertEqual(report["status"], "pass")
            self.assertEqual(report["csv_sequence_mismatch_count"], 0)
            self.assertEqual(report["fasta_sequence_mismatch_count"], 0)

    def test_qc_plot_is_reproducible_from_compact_data(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            output = Path(temporary_dir)
            plot_data = output / "qc_plot_data.csv"
            figure = output / "qc.svg"
            write_qc_plot_data(plot_data, self.records)
            render_qc_svg(plot_data, figure, SOURCE_SHA256)
            svg = figure.read_text(encoding="utf-8")
            self.assertIn("Records by source section", svg)
            self.assertIn("Yield-value semantics", svg)
            self.assertIn(SOURCE_SHA256, svg)

    def test_raw_paragraph_text_is_preserved_separately(self) -> None:
        whitespace_indices = [
            paragraph.nonempty_index
            for paragraph in self.paragraphs
            if paragraph.raw_text != paragraph.parse_text
        ]
        self.assertEqual(whitespace_indices, [40, 42])
        sequence_indices = {record.source_sequence_paragraph_index for record in self.records}
        self.assertTrue(sequence_indices.isdisjoint(whitespace_indices))

    def test_explicit_title_volume_survives_source_rename(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            renamed = Path(temporary_dir) / "renamed_source.docx"
            shutil.copyfile(SOURCE, renamed)
            renamed_records, _ = parse_expression_docx(
                renamed, document_title_culture_volume_l=Decimal("1")
            )
        self.assertEqual(
            [record.culture_volume_l for record in renamed_records],
            [record.culture_volume_l for record in self.records],
        )
        self.assertTrue(
            all(
                record.volume_evidence == "document_title"
                for record in renamed_records
                if record.provider_code in {"LTT", "LLJ"}
            )
        )
        unknown_records, _ = parse_expression_docx(
            SOURCE, document_title_culture_volume_l=None
        )
        unknown_contexts = {
            row["provider_code"]: row for row in assay_context_rows(unknown_records)
        }
        for provider in ("LTT", "LLJ"):
            self.assertEqual(unknown_contexts[provider]["culture_volume_l"], "")
            self.assertEqual(unknown_contexts[provider]["volume_evidence"], "unavailable")
        self.assertEqual(unknown_contexts["WCC"]["culture_volume_l"], "1")
        self.assertEqual(unknown_contexts["WCC"]["volume_evidence"], "entry_text")


if __name__ == "__main__":
    unittest.main()
