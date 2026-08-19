"""Leakage-controlled binary validation for small reported-yield datasets."""

from __future__ import annotations

import math
from typing import Mapping, Sequence

import numpy as np


class YieldClassificationError(ValueError):
    """Raised when the classification validation contract cannot be evaluated."""


def nested_yield_classification(
    rows: Sequence[Mapping[str, object]],
    feature: str,
    *,
    outer_scheme: str,
) -> dict[str, object]:
    """Evaluate a higher-is-better score with all thresholds fitted inside each outer fold.

    The outcome is high versus low yield relative to the median of the matching
    provider in the outer training set. The score threshold maximizes training
    MCC, with balanced accuracy and then the higher threshold as deterministic
    tie-breakers. ``outer_scheme`` is either leave-one-sample-out or
    leave-one-sequence-cluster-out.
    """

    if len(rows) < 8:
        raise YieldClassificationError("At least eight numeric observations are required")
    if outer_scheme not in {"leave_one_out", "leave_one_cluster_out"}:
        raise YieldClassificationError(f"Unsupported outer scheme: {outer_scheme}")
    providers = np.asarray([str(row["provider_code"]) for row in rows])
    values = np.asarray([float(row[feature]) for row in rows], dtype=float)
    yields = np.asarray([float(row["numeric_yield_value"]) for row in rows], dtype=float)
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(yields)):
        raise YieldClassificationError("Classification values must be finite")
    groups = (
        np.asarray([str(row["sample_uid"]) for row in rows])
        if outer_scheme == "leave_one_out"
        else np.asarray([str(row["sequence_cluster_90"]) for row in rows])
    )

    predictions: list[dict[str, object]] = []
    for fold_id, held_group in enumerate(dict.fromkeys(groups.tolist()), 1):
        held = groups == held_group
        train = ~held
        yield_thresholds = _provider_medians(yields[train], providers[train])
        train_labels = np.asarray(
            [int(value >= yield_thresholds[provider]) for value, provider in zip(yields[train], providers[train], strict=True)]
        )
        if len(set(train_labels.tolist())) != 2:
            raise YieldClassificationError("An outer training fold has only one yield class")
        score_threshold = _best_mcc_threshold(values[train], train_labels)
        for index in np.flatnonzero(held):
            provider = providers[index]
            if provider not in yield_thresholds:
                raise YieldClassificationError(f"No training yield threshold for provider {provider}")
            label = int(yields[index] >= yield_thresholds[provider])
            prediction = int(values[index] >= score_threshold)
            predictions.append(
                {
                    "outer_scheme": outer_scheme,
                    "fold_id": fold_id,
                    "held_group": str(held_group),
                    "sample_uid": str(rows[index]["sample_uid"]),
                    "provider_code": provider,
                    "numeric_yield_value": yields[index],
                    "provider_training_yield_threshold": yield_thresholds[provider],
                    "feature": feature,
                    "feature_value": values[index],
                    "training_score_threshold": score_threshold,
                    "observed_high_yield": label,
                    "predicted_high_yield": prediction,
                }
            )
    if len(predictions) != len(rows):
        raise YieldClassificationError("Outer-fold predictions do not cover all numeric rows")
    return {
        "summary": _summarize_predictions(predictions),
        "prediction_rows": predictions,
    }


def fixed_yield_nested_classification(
    rows: Sequence[Mapping[str, object]],
    feature: str,
    *,
    outer_scheme: str,
    yield_threshold: float,
) -> dict[str, object]:
    """Evaluate a higher-is-better score against one fixed numeric-yield cutoff.

    The yield label is fixed before resampling.  Within each outer fold, only
    the score cutoff is fitted by maximizing training MCC, with balanced
    accuracy and then the higher cutoff used as deterministic tie-breakers.
    """

    values, yields, groups = _fixed_yield_inputs(rows, feature, outer_scheme, yield_threshold)
    labels = (yields >= yield_threshold).astype(int)
    predictions: list[dict[str, object]] = []
    for fold_id, held_group in enumerate(dict.fromkeys(groups.tolist()), 1):
        held = groups == held_group
        train = ~held
        if len(set(labels[train].tolist())) != 2:
            raise YieldClassificationError("An outer training fold has only one yield class")
        score_threshold = _best_mcc_threshold(values[train], labels[train])
        predictions.extend(
            _fixed_prediction_row(
                rows[index], feature, values[index], yields[index], labels[index],
                score_threshold, yield_threshold, outer_scheme, fold_id, str(held_group),
            )
            for index in np.flatnonzero(held)
        )
    if len(predictions) != len(rows):
        raise YieldClassificationError("Outer-fold predictions do not cover all numeric rows")
    summary = _summarize_predictions(predictions)
    summary.update({"fixed_yield_threshold": float(yield_threshold), "threshold_fit_scope": "outer_training_fold"})
    return {"summary": summary, "prediction_rows": predictions}


def fixed_yield_apparent_classification(
    rows: Sequence[Mapping[str, object]],
    feature: str,
    *,
    yield_threshold: float,
) -> dict[str, object]:
    """Fit and describe one score cutoff on all rows without claiming test performance."""

    values, yields, _ = _fixed_yield_inputs(rows, feature, "leave_one_out", yield_threshold)
    labels = (yields >= yield_threshold).astype(int)
    if len(set(labels.tolist())) != 2:
        raise YieldClassificationError("The fixed yield threshold produces only one class")
    score_threshold = _best_mcc_threshold(values, labels)
    predictions = [
        _fixed_prediction_row(
            row, feature, value, observed_yield, label, score_threshold,
            yield_threshold, "apparent_full_sample", 1, "all_numeric_rows",
        )
        for row, value, observed_yield, label in zip(rows, values, yields, labels, strict=True)
    ]
    summary = _summarize_predictions(predictions)
    summary.update({"fixed_yield_threshold": float(yield_threshold), "threshold_fit_scope": "all_numeric_rows_apparent"})
    return {"summary": summary, "prediction_rows": predictions}


def _fixed_yield_inputs(
    rows: Sequence[Mapping[str, object]], feature: str, outer_scheme: str, yield_threshold: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if len(rows) < 8:
        raise YieldClassificationError("At least eight numeric observations are required")
    if outer_scheme not in {"leave_one_out", "leave_one_cluster_out"}:
        raise YieldClassificationError(f"Unsupported outer scheme: {outer_scheme}")
    if not math.isfinite(yield_threshold):
        raise YieldClassificationError("The fixed yield threshold must be finite")
    values = np.asarray([float(row[feature]) for row in rows], dtype=float)
    yields = np.asarray([float(row["numeric_yield_value"]) for row in rows], dtype=float)
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(yields)):
        raise YieldClassificationError("Classification values must be finite")
    groups = (
        np.asarray([str(row["sample_uid"]) for row in rows])
        if outer_scheme == "leave_one_out"
        else np.asarray([str(row["sequence_cluster_90"]) for row in rows])
    )
    return values, yields, groups


def _fixed_prediction_row(
    row: Mapping[str, object],
    feature: str,
    feature_value: float,
    observed_yield: float,
    observed_label: int,
    score_threshold: float,
    yield_threshold: float,
    outer_scheme: str,
    fold_id: int,
    held_group: str,
) -> dict[str, object]:
    return {
        "outer_scheme": outer_scheme,
        "fold_id": fold_id,
        "held_group": held_group,
        "sample_uid": str(row["sample_uid"]),
        "provider_code": str(row["provider_code"]),
        "numeric_yield_value": float(observed_yield),
        "fixed_yield_threshold": float(yield_threshold),
        "feature": feature,
        "feature_value": float(feature_value),
        "training_score_threshold": float(score_threshold),
        "observed_high_yield": int(observed_label),
        "predicted_high_yield": int(feature_value >= score_threshold),
    }


def _provider_medians(yields: np.ndarray, providers: np.ndarray) -> dict[str, float]:
    result = {}
    for provider in sorted(set(providers.tolist())):
        selected = yields[providers == provider]
        if len(selected) < 3:
            raise YieldClassificationError(f"Insufficient training observations for provider {provider}")
        result[provider] = float(np.median(selected))
    return result


def _best_mcc_threshold(values: np.ndarray, labels: np.ndarray) -> float:
    unique = np.unique(values)
    candidates = [float(np.nextafter(unique[0], -np.inf))]
    candidates.extend(float((left + right) / 2) for left, right in zip(unique[:-1], unique[1:], strict=True))
    candidates.append(float(np.nextafter(unique[-1], np.inf)))
    ranked = []
    for threshold in candidates:
        metrics = _binary_metrics(labels, (values >= threshold).astype(int))
        ranked.append((metrics["mcc"], metrics["balanced_accuracy"], threshold))
    return max(ranked)[2]


def _summarize_predictions(rows: Sequence[Mapping[str, object]]) -> dict[str, object]:
    labels = np.asarray([int(row["observed_high_yield"]) for row in rows])
    predicted = np.asarray([int(row["predicted_high_yield"]) for row in rows])
    scores = np.asarray([float(row["feature_value"]) for row in rows])
    thresholds = np.asarray([float(row["training_score_threshold"]) for row in rows])
    metrics = _binary_metrics(labels, predicted)
    metrics.update(
        {
            "n": len(rows),
            "positive_count": int(labels.sum()),
            "prevalence": float(labels.mean()),
            "roc_auc": _roc_auc(labels, scores),
            "pr_auc_average_precision": _average_precision(labels, scores),
            "score_threshold_median": float(np.median(thresholds)),
            "score_threshold_q1": float(np.quantile(thresholds, 0.25)),
            "score_threshold_q3": float(np.quantile(thresholds, 0.75)),
            "score_threshold_min": float(thresholds.min()),
            "score_threshold_max": float(thresholds.max()),
        }
    )
    return metrics


def _binary_metrics(labels: np.ndarray, predictions: np.ndarray) -> dict[str, float | int]:
    tp = int(np.sum((labels == 1) & (predictions == 1)))
    tn = int(np.sum((labels == 0) & (predictions == 0)))
    fp = int(np.sum((labels == 0) & (predictions == 1)))
    fn = int(np.sum((labels == 1) & (predictions == 0)))
    sensitivity = tp / (tp + fn) if tp + fn else 0.0
    specificity = tn / (tn + fp) if tn + fp else 0.0
    denominator = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return {
        "tp": tp,
        "tn": tn,
        "fp": fp,
        "fn": fn,
        "mcc": (tp * tn - fp * fn) / denominator if denominator else 0.0,
        "balanced_accuracy": (sensitivity + specificity) / 2,
        "sensitivity": sensitivity,
        "specificity": specificity,
    }


def _roc_auc(labels: np.ndarray, scores: np.ndarray) -> float:
    positives = scores[labels == 1]
    negatives = scores[labels == 0]
    if not len(positives) or not len(negatives):
        return 0.5
    wins = sum(float(pos > neg) + 0.5 * float(pos == neg) for pos in positives for neg in negatives)
    return wins / (len(positives) * len(negatives))


def _average_precision(labels: np.ndarray, scores: np.ndarray) -> float:
    positives = int(labels.sum())
    if positives == 0:
        return 0.0
    order = np.argsort(-scores, kind="stable")
    ranked = labels[order]
    precision = np.cumsum(ranked) / np.arange(1, len(ranked) + 1)
    return float(np.sum(precision * ranked) / positives)
