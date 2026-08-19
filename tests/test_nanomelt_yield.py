from __future__ import annotations

import csv
from types import SimpleNamespace
from pathlib import Path

import pytest
import antibody_optimization.nanomelt_yield as nanomelt_yield

from antibody_optimization.nanomelt_yield import (
    NanoMeltYieldError,
    analyze_nanomelt_associations,
    build_nanomelt_validation_inputs,
    normalize_nanomelt_scores,
    verify_anarci_runtime,
    verify_required_openmm_platforms,
)
from antibody_optimization.nanomelt_yield_plot import render_nanomelt_yield_figure


ROOT = Path(__file__).resolve().parents[1]
TEST_ONLY_NANOMELT_NOT_SCORED = {"WCC__4-1", "WCC__4-28", "WCC__4-11", "WCC__4-42"}


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _real_samples() -> list[dict[str, object]]:
    return build_nanomelt_validation_inputs(
        _csv(ROOT / "docs/result_artifacts/nb_expression/nb_expression_records.csv"),
        _csv(ROOT / "docs/result_artifacts/input_baseline/sequence/sequence_numbering_review.csv"),
        _csv(ROOT / "docs/result_artifacts/input_baseline/sequence/sequence_numbering_positions.csv"),
    )["sample_rows"]


def _raw_scores(samples, *, correlated: bool = False, omit=frozenset()):
    rows = []
    for index, sample in enumerate(samples):
        if sample["sample_uid"] in omit:
            continue
        sequence = str(sample["sequence_raw"])
        scored = sequence[:-2] if sample["sample_uid"] == "LTT__Nb252" else sequence
        if correlated and sample["observation_semantics"] == "individual_approximate":
            tm = 50.0 + float(sample["numeric_yield_value"])
        elif correlated:
            tm = 50.0 + 5.0 * int(sample["llj_ordinal_level"])
        else:
            tm = 55.0 + index / 10.0
        rows.append({"ID": sample["sample_uid"], "Aligned Sequence": scored, "Sequence": scored, "NanoMelt Tm (C)": tm})
    return rows


def test_real_plan_preserves_47_samples_and_existing_clusters() -> None:
    samples = _real_samples()
    assert len(samples) == 47
    assert len({row["sample_uid"] for row in samples}) == 47
    assert len({row["sequence_cluster_90"] for row in samples}) == 40


def test_anarci_runtime_validation_ignores_stale_distribution_metadata(tmp_path: Path) -> None:
    package = tmp_path / "lib" / "python3.10" / "site-packages" / "anarci"
    hmm_dir = package / "dat" / "HMMs"
    hmm_dir.mkdir(parents=True)
    module_file = package / "__init__.py"
    module_file.write_text("", encoding="utf-8")
    for suffix in ("", ".h3f", ".h3i", ".h3m", ".h3p"):
        (hmm_dir / f"ALL.hmm{suffix}").write_text("test-only", encoding="utf-8")
    conda_meta = tmp_path / "conda-meta"
    conda_meta.mkdir()
    (conda_meta / "anarci-2024.05.21-test.json").write_text(
        '{"name":"anarci","version":"2024.05.21"}', encoding="utf-8"
    )
    module = SimpleNamespace(
        __file__=str(module_file),
        anarci=lambda: None,
        run_anarci=lambda: None,
        validate_sequence=lambda: None,
        scheme_short_to_long={},
    )
    result = verify_anarci_runtime(module, tmp_path, expected_conda_version="2024.05.21")
    assert result["conda_package_version"] == "2024.05.21"
    assert result["hmm_pressed_index_count"] == 4


def test_openmm_platform_gate_allows_compute_node_extras() -> None:
    assert verify_required_openmm_platforms(
        ["Reference", "CPU", "OpenCL"], ["Reference", "CPU"]
    ) == ["Reference", "CPU", "OpenCL"]
    with pytest.raises(NanoMeltYieldError, match="lacks required platforms"):
        verify_required_openmm_platforms(["Reference", "OpenCL"], ["Reference", "CPU"])


def test_normalization_records_exact_nb252_domain_trimming() -> None:
    samples = _real_samples()
    normalized = normalize_nanomelt_scores(samples, _raw_scores(samples))
    nb252 = next(row for row in normalized if row["sample_uid"] == "LTT__Nb252")
    assert len(nb252["sequence_raw"]) == 128
    assert nb252["scored_length_aa"] == 126
    assert nb252["trimmed_n_terminal"] == ""
    assert nb252["trimmed_c_terminal"] == "GS"


def test_normalization_preserves_47_identities_with_four_nanomelt_not_scored() -> None:
    samples = _real_samples()
    normalized = normalize_nanomelt_scores(
        samples,
        _raw_scores(samples, omit=TEST_ONLY_NANOMELT_NOT_SCORED),
        expected_pass_count=43,
    )
    assert len(normalized) == 47
    assert sum(row["scoring_status"] == "pass" for row in normalized) == 43
    excluded = {row["sample_uid"] for row in normalized if row["scoring_status"] == "nanomelt_not_scored"}
    assert excluded == TEST_ONLY_NANOMELT_NOT_SCORED


def test_normalization_rejects_nonunique_or_mismatched_domains() -> None:
    samples = _real_samples()
    raw = _raw_scores(samples)
    raw[0]["Sequence"] = "AAAA"
    raw[0]["Aligned Sequence"] = "AAAA"
    with pytest.raises(NanoMeltYieldError, match="unique raw-sequence segment"):
        normalize_nanomelt_scores(samples, raw)


def test_association_reports_cluster_cv_llj_and_nb252_influence(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(nanomelt_yield, "RESAMPLING_REPLICATES", 200)
    samples = _real_samples()
    scores = normalize_nanomelt_scores(
        samples,
        _raw_scores(samples, correlated=True, omit=TEST_ONLY_NANOMELT_NOT_SCORED),
        expected_pass_count=43,
    )
    result = analyze_nanomelt_associations(samples, scores)
    primary = result["primary"]
    assert primary["numeric_n"] == 27
    assert primary["llj_ordinal_n"] == 16
    assert primary["not_scored_count"] == 4
    assert primary["stratified_spearman_rho"] > 0
    assert primary["without_nb252_stratified_spearman_rho"] > 0
    assert len(result["leave_one_out_rows"]) == 27
    assert {row["model"] for row in result["cv_rows"]} == {
        "provider_only",
        "provider_plus_nanomelt_tm",
    }
    assert result["evidence_level"] == "weak_ranking_evidence"
    assert len(result["classification_rows"]) == 2
    assert len(result["classification_prediction_rows"]) == 54
    assert all(float(row["roc_auc"]) > 0.9 for row in result["classification_rows"])


def test_plot_renders_exact_compact_analysis(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(nanomelt_yield, "RESAMPLING_REPLICATES", 200)
    samples = _real_samples()
    result = analyze_nanomelt_associations(
        samples,
        normalize_nanomelt_scores(
            samples,
            _raw_scores(samples, correlated=True, omit=TEST_ONLY_NANOMELT_NOT_SCORED),
            expected_pass_count=43,
        ),
    )
    png, svg = tmp_path / "result.png", tmp_path / "result.svg"
    render_nanomelt_yield_figure(
        result["sample_rows"],
        result["primary"],
        result["cv_rows"],
        result["leave_one_out_rows"],
        result["classification_rows"],
        result["classification_prediction_rows"],
        png_path=png,
        svg_path=svg,
    )
    assert png.stat().st_size > 10_000
    assert svg.stat().st_size > 1_000
