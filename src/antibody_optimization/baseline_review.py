"""Validate the unified chain, orange, and construct baseline review.

This module owns the stable, machine-readable contract shared by the structure
baseline and temporary-interface entry points.  It deliberately does not infer
chain roles, decide which saved-session orange class represents the collaborator
annotation, or infer an authoritative construct boundary.  Passing records must
bind each of those claims to explicit review evidence.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime
from pathlib import Path
from typing import Mapping, Sequence

from .structure_inventory import ChainSelector


REVIEW_SCHEMA_VERSION = 1
REVIEW_TYPE = "nb252_structure_orange_and_construct_baseline_review"
BUILTIN_ORANGE_RGB = (255, 165, 0)
ROLE_FIELDS = (
    "source_model_name",
    "auth_asym_id",
    "label_asym_id",
    "entity_id",
    "confirmed_role",
    "confirmation_status",
    "confirmed_by",
    "confirmed_at",
    "confirmation_note",
)
ORANGE_CHANNELS = (
    "ribbon",
    "atom",
    "surface_vertex_or_uniform_patch",
)
STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
_HISTOGRAM_CHANNELS = (
    ("atom_rgba_histogram_json", "atom"),
    ("surface_rgba_histogram_json", "surface_vertex_or_uniform_patch"),
)


class BaselineReviewError(ValueError):
    """Raised when a baseline review or its bound evidence is inconsistent."""


def role_key(row: Mapping[str, object]) -> tuple[str, str, str, str]:
    """Return the complete inventory identity used to bind a reviewed chain."""

    return (
        str(row.get("source_model_name", "")),
        str(row.get("auth_asym_id", "")),
        str(row.get("label_asym_id", "")),
        str(row.get("entity_id", "")),
    )


def chain_identity_sha256(inventory: Sequence[Mapping[str, object]]) -> str:
    """Hash the immutable identity/count fields of a structure chain inventory."""

    fields = (
        "source_model_name",
        "auth_asym_id",
        "label_asym_id",
        "entity_id",
        "entity_type",
        "observed_sequence_sha256",
        "residue_count",
        "atom_count",
    )
    identities = [
        {field: row.get(field, "") for field in fields}
        for row in sorted(inventory, key=role_key)
    ]
    payload = json.dumps(
        identities, sort_keys=True, separators=(",", ":")
    ).encode()
    return hashlib.sha256(payload).hexdigest()


def color_channels_for_rgba(
    row: Mapping[str, object], rgba: tuple[int, int, int, int]
) -> list[str]:
    """Return saved color channels containing one exact RGBA class.

    Ribbon colors and positive-count atom/surface histograms are considered.
    Selection state is intentionally absent from this contract.
    """

    channels: list[str] = []
    ribbon = _rgba_tuple(row.get("ribbon_rgba", ""), allow_empty=True)
    if ribbon == rgba:
        channels.append("ribbon")
    for field, channel in _HISTOGRAM_CHANNELS:
        histogram = _color_histogram(row.get(field, "{}"), field)
        if any(color == rgba and count > 0 for color, count in histogram.items()):
            channels.append(channel)
    return channels


def build_review_template(
    *,
    inventory_rows: Sequence[Mapping[str, object]],
    color_rows: Sequence[Mapping[str, object]],
    source_binding: Mapping[str, str],
    reported_sequence: str | None = None,
    reported_sequence_sha256: str | None = None,
) -> dict[str, object]:
    """Build the pending unified review without assigning scientific roles.

    Every exported RGBA class whose RGB is exactly ``(255, 165, 0)`` is listed
    separately.  Alpha is retained exactly, so transparent orange cannot be
    silently collapsed into the opaque built-in class.
    """

    chain_reviews = [
        {
            "source_model_name": row["source_model_name"],
            "auth_asym_id": row["auth_asym_id"],
            "label_asym_id": row["label_asym_id"],
            "entity_id": row["entity_id"],
            "confirmed_role": "",
            "confirmation_status": "pending_user_review",
            "confirmed_by": "",
            "confirmed_at": "",
            "confirmation_note": "",
            "provisional_role_suggestion": _provisional_role_suggestion(row),
            "suggestion_evidence": _suggestion_evidence(row),
        }
        for row in sorted(inventory_rows, key=role_key)
    ]

    candidates: list[dict[str, object]] = []
    candidate_residues: set[tuple[str, ...]] = set()
    for row in color_rows:
        residue_identity = (
            str(row.get("model_name", "")),
            str(row.get("chimerax_chain_id", "")),
            str(row.get("mmcif_chain_id", "")),
            str(row.get("auth_seq_id", "")),
            _clean_code(row.get("insertion_code", "")),
            str(row.get("residue_name", "")),
        )
        for rgba, channels in _exact_orange_classes(row):
            candidate_residues.add(residue_identity)
            candidates.append(
                {
                    "model_name": residue_identity[0],
                    "chimerax_chain_id": residue_identity[1],
                    "mmcif_chain_id": residue_identity[2],
                    "auth_seq_id": residue_identity[3],
                    "insertion_code": residue_identity[4],
                    "residue_name": residue_identity[5],
                    "exact_rgba": list(rgba),
                    "candidate_channels": channels,
                    "status": "candidate_unconfirmed",
                }
            )
    candidates.sort(
        key=lambda row: (
            str(row["model_name"]),
            str(row["chimerax_chain_id"]),
            str(row["mmcif_chain_id"]),
            str(row["auth_seq_id"]),
            str(row["insertion_code"]),
            str(row["residue_name"]),
            tuple(row["exact_rgba"]),
        )
    )

    return {
        "schema_version": REVIEW_SCHEMA_VERSION,
        "review_type": REVIEW_TYPE,
        "status": "pending_user_review",
        "source_binding": dict(source_binding),
        "instructions": {
            "single_review_contract": (
                "Confirm every chain role and the collaborator orange annotation "
                "in this one file, then save it as baseline_review.json."
            ),
            "confirmed_role_examples": ["Nb252_VHH", "NK2R", "NKA", "other"],
            "orange_channels": list(ORANGE_CHANNELS),
            "orange_caution": (
                "Exact saved-session color is candidate evidence only; cite the "
                "collaborator annotation or an explicit visual review."
            ),
        },
        "chain_reviews": chain_reviews,
        "orange_annotation_review": {
            "status": "pending_user_review",
            "confirmed_rgb": [],
            "confirmed_rgba": [],
            "confirmed_channels": [],
            "confirmed_by": "",
            "confirmed_at": "",
            "confirmation_note": "",
            "evidence": "",
            "candidate_exact_builtin_orange_residue_count": len(
                candidate_residues
            ),
            "candidates": candidates,
        },
        "authoritative_construct_review": {
            "status": "pending_collaborator_confirmation",
            "reported_sequence_sha256": reported_sequence_sha256 or "",
            "reported_sequence_length_aa": (
                len(reported_sequence) if reported_sequence is not None else ""
            ),
            "authoritative_sequence": "",
            "authoritative_sequence_sha256": "",
            "reported_start_1based_inclusive": "",
            "reported_end_1based_inclusive": "",
            "construct_scope": "",
            "terminal_gs_interpretation": "",
            "confirmed_by": "",
            "confirmed_at": "",
            "confirmation_note": "",
            "evidence": "",
            "instructions": (
                "Only set status=confirmed from explicit construct evidence. The exact "
                "authoritative sequence must be a contiguous literal span of the reported "
                "sequence, with 1-based inclusive boundaries and an explicit disposition "
                "for the reported terminal GS."
            ),
        },
    }


def load_authoritative_construct_confirmation(
    path: Path,
    *,
    reported_sequence: str,
    reported_sequence_sha256: str,
) -> dict[str, object] | None:
    """Return a strict construct confirmation, or ``None`` while still pending.

    A reviewer name alone is deliberately insufficient.  A passing record is
    bound to the exact reported sequence hash, an exact authoritative
    subsequence and inclusive boundaries, construct scope, terminal-GS
    interpretation, timestamp, and evidence.
    """

    review = _read_json_object(path)
    raw = review.get("authoritative_construct_review")
    if raw is None:
        return None
    if not isinstance(raw, dict):
        raise BaselineReviewError(
            "authoritative_construct_review must be an object"
        )
    status = raw.get("status")
    if status in {"pending_collaborator_confirmation", "pending_user_review"}:
        return None
    if status != "confirmed":
        raise BaselineReviewError(
            "authoritative_construct_review status is unsupported"
        )
    calculated_reported_hash = hashlib.sha256(
        reported_sequence.encode("ascii")
    ).hexdigest()
    if calculated_reported_hash != reported_sequence_sha256.lower():
        raise BaselineReviewError("reported sequence hash is internally inconsistent")
    if raw.get("reported_sequence_sha256") != reported_sequence_sha256:
        raise BaselineReviewError(
            "construct review is not bound to the current reported sequence"
        )
    authoritative_sequence = raw.get("authoritative_sequence")
    if (
        not isinstance(authoritative_sequence, str)
        or not authoritative_sequence
        or set(authoritative_sequence) - STANDARD_AMINO_ACIDS
    ):
        raise BaselineReviewError(
            "authoritative construct sequence must use standard uppercase amino acids"
        )
    authoritative_hash = hashlib.sha256(
        authoritative_sequence.encode("ascii")
    ).hexdigest()
    if raw.get("authoritative_sequence_sha256") != authoritative_hash:
        raise BaselineReviewError(
            "authoritative construct sequence SHA-256 is absent or incorrect"
        )
    try:
        start = int(raw.get("reported_start_1based_inclusive"))
        end = int(raw.get("reported_end_1based_inclusive"))
    except (TypeError, ValueError) as exc:
        raise BaselineReviewError(
            "authoritative construct boundaries must be integer 1-based positions"
        ) from exc
    if not 1 <= start <= end <= len(reported_sequence):
        raise BaselineReviewError("authoritative construct boundaries are out of range")
    if reported_sequence[start - 1 : end] != authoritative_sequence:
        raise BaselineReviewError(
            "authoritative construct sequence does not equal its reported-sequence span"
        )
    scope = raw.get("construct_scope")
    if scope not in {"mature_vhh", "full_expression_construct"}:
        raise BaselineReviewError(
            "construct_scope must be mature_vhh or full_expression_construct"
        )
    gs_interpretation = raw.get("terminal_gs_interpretation")
    allowed_gs = {
        "included_in_authoritative_construct",
        "excluded_expression_flank",
        "not_present_in_reported_sequence",
    }
    if gs_interpretation not in allowed_gs:
        raise BaselineReviewError(
            "terminal_gs_interpretation is absent or unsupported"
        )
    if reported_sequence.endswith("GS"):
        if end == len(reported_sequence) and gs_interpretation != (
            "included_in_authoritative_construct"
        ):
            raise BaselineReviewError(
                "a construct ending at the reported C terminus must explicitly include GS"
            )
        if end == len(reported_sequence) - 2 and gs_interpretation != (
            "excluded_expression_flank"
        ):
            raise BaselineReviewError(
                "a construct ending before terminal GS must mark it as excluded flank"
            )
        if end not in {len(reported_sequence), len(reported_sequence) - 2}:
            raise BaselineReviewError(
                "reported terminal GS disposition is unresolved for these boundaries"
            )
    elif gs_interpretation != "not_present_in_reported_sequence":
        raise BaselineReviewError(
            "terminal GS cannot be interpreted when absent from the reported sequence"
        )
    required_text = (
        "confirmed_by",
        "confirmed_at",
        "confirmation_note",
        "evidence",
    )
    for field in required_text:
        if not isinstance(raw.get(field), str) or not raw[field].strip():
            raise BaselineReviewError(
                f"authoritative construct review field {field!r} is required"
            )
    try:
        confirmed_at = datetime.fromisoformat(str(raw["confirmed_at"]))
    except ValueError as exc:
        raise BaselineReviewError(
            "authoritative construct confirmed_at must be ISO-8601"
        ) from exc
    if confirmed_at.tzinfo is None:
        raise BaselineReviewError(
            "authoritative construct confirmed_at must include a UTC offset"
        )
    return {
        "status": "confirmed",
        "reported_sequence_sha256": reported_sequence_sha256,
        "reported_start_1based_inclusive": start,
        "reported_end_1based_inclusive": end,
        "authoritative_sequence": authoritative_sequence,
        "authoritative_sequence_sha256": authoritative_hash,
        "length_aa": len(authoritative_sequence),
        "construct_scope": scope,
        "terminal_gs_interpretation": gs_interpretation,
        **{field: str(raw[field]) for field in required_text},
    }


def confirmed_review_components(
    review: Mapping[str, object], *, orange_semantics: str | None = None
) -> tuple[list[dict[str, str]], dict[str, object]]:
    """Validate and normalize one confirmed unified review object.

    This checks the review schema itself.  Call :func:`load_confirmed_review`
    when the chain set, source hashes, and selected color class must also be
    bound to concrete current exports.
    """

    if (
        review.get("schema_version") != REVIEW_SCHEMA_VERSION
        or review.get("review_type") != REVIEW_TYPE
    ):
        raise BaselineReviewError("confirmed review has the wrong schema/review_type")
    if review.get("status") != "confirmed":
        raise BaselineReviewError("baseline_review.json status must be confirmed")
    if not isinstance(review.get("source_binding"), dict):
        raise BaselineReviewError("baseline review source_binding must be an object")

    raw_rows = review.get("chain_reviews")
    if not isinstance(raw_rows, list):
        raise BaselineReviewError("baseline review chain_reviews must be a list")
    roles: list[dict[str, str]] = []
    seen: set[tuple[str, str, str, str]] = set()
    for raw in raw_rows:
        if not isinstance(raw, dict):
            raise BaselineReviewError("baseline review chain row is not an object")
        row = {
            str(key): "" if value is None else str(value)
            for key, value in raw.items()
        }
        missing = [field for field in ROLE_FIELDS if field not in row]
        if missing:
            raise BaselineReviewError(
                f"baseline review chain row lacks fields: {missing}"
            )
        key = role_key(row)
        if key in seen:
            raise BaselineReviewError(f"duplicate baseline review chain row: {key}")
        required = (
            "confirmed_role",
            "confirmed_by",
            "confirmed_at",
            "confirmation_note",
        )
        if row["confirmation_status"] != "confirmed" or any(
            not row[field] for field in required
        ):
            raise BaselineReviewError(
                f"chain review row is not fully confirmed: {key}"
            )
        seen.add(key)
        roles.append(row)

    orange = review.get("orange_annotation_review")
    if not isinstance(orange, dict) or orange.get("status") != "confirmed":
        raise BaselineReviewError("orange_annotation_review must be confirmed")
    rgb = _integer_channel_tuple(
        orange.get("confirmed_rgb"), 3, "confirmed_rgb"
    )
    rgba = _integer_channel_tuple(
        orange.get("confirmed_rgba"), 4, "confirmed_rgba"
    )
    if rgba[:3] != rgb:
        raise BaselineReviewError("confirmed RGB and RGBA disagree")
    channels = orange.get("confirmed_channels")
    allowed = set(ORANGE_CHANNELS)
    if (
        not isinstance(channels, list)
        or not channels
        or not all(isinstance(channel, str) for channel in channels)
    ):
        raise BaselineReviewError(
            "confirmed orange channels are empty or unsupported"
        )
    channel_set = set(channels)
    if not channel_set <= allowed:
        raise BaselineReviewError(
            "confirmed orange channels are empty or unsupported"
        )
    for field in ("confirmed_by", "confirmed_at", "confirmation_note", "evidence"):
        if not isinstance(orange.get(field), str) or not orange[field].strip():
            raise BaselineReviewError(f"orange review field {field!r} is required")

    normalized: dict[str, object] = {
        "status": "confirmed",
        "confirmed_rgb": list(rgb),
        "confirmed_rgba": list(rgba),
        "confirmed_channels": sorted(channel_set),
        **{
            field: orange[field]
            for field in (
                "confirmed_by",
                "confirmed_at",
                "confirmation_note",
                "evidence",
            )
        },
    }
    if orange_semantics is not None:
        normalized["semantics"] = orange_semantics
    return roles, normalized


def load_confirmed_review(
    path: Path,
    *,
    inventory_rows: Sequence[Mapping[str, object]],
    color_rows: Sequence[Mapping[str, object]],
    expected_binding: Mapping[str, str],
) -> tuple[
    dict[tuple[str, str, str, str], dict[str, str]], dict[str, object]
]:
    """Load a confirmed review and bind it to current chains and color export."""

    review = _read_json_object(path)
    roles_in_order, orange = confirmed_review_components(review)
    if review.get("source_binding") != dict(expected_binding):
        raise BaselineReviewError(
            "baseline review source_binding does not match current inputs"
        )
    roles = {role_key(row): row for row in roles_in_order}
    if set(roles) != {role_key(row) for row in inventory_rows}:
        raise BaselineReviewError(
            "baseline review chain set does not match inventory"
        )
    rgba = tuple(orange["confirmed_rgba"])
    channels = set(orange["confirmed_channels"])
    if not any(
        set(color_channels_for_rgba(row, rgba)) & channels for row in color_rows
    ):
        raise BaselineReviewError(
            "confirmed orange RGBA/channel is absent from export"
        )
    return roles, orange


def validated_required_roles(
    roles: Mapping[tuple[str, str, str, str], Mapping[str, str]],
    *,
    reference_model: str,
    af3_model: str,
) -> dict[tuple[str, str], list[ChainSelector]]:
    """Validate the roles needed for Nb252 mapping and return chain selectors."""

    result: dict[tuple[str, str], list[ChainSelector]] = {}
    for (model, auth, label, _entity_id), row in roles.items():
        selector = ChainSelector(model, auth, label)
        result.setdefault((model, row["confirmed_role"]), []).append(selector)
    for model in (reference_model, af3_model):
        if len(result.get((model, "Nb252_VHH"), [])) != 1:
            raise BaselineReviewError(
                f"{model} must have exactly one confirmed Nb252_VHH chain"
            )
    if not result.get((reference_model, "NK2R")):
        raise BaselineReviewError(
            f"{reference_model} requires at least one confirmed NK2R chain"
        )
    return result


def interface_selectors_from_roles(
    rows: Sequence[Mapping[str, str]], *, experimental_model: str
) -> tuple[ChainSelector, list[ChainSelector]]:
    """Return the reviewed VHH and receptor selectors for one experiment model."""

    confirmed = [
        row
        for row in rows
        if row.get("source_model_name") == experimental_model
        and row.get("confirmation_status") == "confirmed"
    ]
    vhh_rows = [row for row in confirmed if row.get("confirmed_role") == "Nb252_VHH"]
    receptor_rows = [row for row in confirmed if row.get("confirmed_role") == "NK2R"]
    if len(vhh_rows) != 1 or not receptor_rows:
        raise BaselineReviewError(
            "experimental model needs one confirmed Nb252_VHH and at least one NK2R"
        )

    def make_selector(row: Mapping[str, str]) -> ChainSelector:
        return ChainSelector(
            experimental_model, row["auth_asym_id"], row["label_asym_id"]
        )

    return make_selector(vhh_rows[0]), [make_selector(row) for row in receptor_rows]


def confirmed_orange_residue_channels(
    rows: Sequence[Mapping[str, object]],
    selector: ChainSelector,
    rgba: tuple[int, int, int, int],
    confirmed_channels: set[str],
) -> dict[tuple[str, str, int, str, str], list[str]]:
    """Map exact reviewed RGBA/channel evidence to residues of one reviewed chain."""

    if not confirmed_channels or not confirmed_channels <= set(ORANGE_CHANNELS):
        raise BaselineReviewError("confirmed orange channels are empty or unsupported")
    result: dict[tuple[str, str, int, str, str], list[str]] = {}
    for row in rows:
        if (
            row.get("model_name") != selector.model_name
            or row.get("chimerax_chain_id") != selector.auth_asym_id
            or row.get("mmcif_chain_id") != selector.label_asym_id
        ):
            continue
        channels = [
            channel
            for channel in color_channels_for_rgba(row, rgba)
            if channel in confirmed_channels
        ]
        if not channels:
            continue
        try:
            auth_seq_id = int(str(row["auth_seq_id"]))
        except (KeyError, TypeError, ValueError) as exc:
            raise BaselineReviewError(
                "orange color row has an invalid auth_seq_id"
            ) from exc
        key = (
            selector.auth_asym_id,
            selector.label_asym_id,
            auth_seq_id,
            _clean_code(row.get("insertion_code", "")),
            str(row.get("residue_name", "")),
        )
        if key in result:
            raise BaselineReviewError(f"duplicate orange residue color row: {key}")
        result[key] = channels
    return result


def selector_dict(selector: ChainSelector) -> dict[str, str]:
    """Serialize a chain selector using the shared manifest field names."""

    return {
        "model_name": selector.model_name,
        "auth_asym_id": selector.auth_asym_id,
        "label_asym_id": selector.label_asym_id,
    }


def _exact_orange_classes(
    row: Mapping[str, object],
) -> list[tuple[tuple[int, int, int, int], list[str]]]:
    colors: set[tuple[int, int, int, int]] = set()
    ribbon = _rgba_tuple(row.get("ribbon_rgba", ""), allow_empty=True)
    if ribbon and ribbon[:3] == BUILTIN_ORANGE_RGB:
        colors.add(ribbon)
    for field, _channel in _HISTOGRAM_CHANNELS:
        for rgba, count in _color_histogram(row.get(field, "{}"), field).items():
            if count > 0 and rgba[:3] == BUILTIN_ORANGE_RGB:
                colors.add(rgba)
    return [
        (rgba, color_channels_for_rgba(row, rgba)) for rgba in sorted(colors)
    ]


def _provisional_role_suggestion(row: Mapping[str, object]) -> str:
    exact_sequence = any(
        bool(row.get(field))
        for field in (
            "exact_source_polymer_match_to_reported",
            "exact_entity_match_to_authoritative",
            "exact_observed_match_to_authoritative",
        )
    )
    if exact_sequence:
        return "Nb252_VHH_candidate"
    description = str(row.get("source_entity_description", "")).lower()
    if any(term in description for term in ("nanobody", "vhh", "antibody")):
        return "VHH_candidate_from_source_entity_description"
    return "unassigned"


def _suggestion_evidence(row: Mapping[str, object]) -> str:
    """Serialize only explicit sequence/entity evidence; never chain order."""

    evidence = {
        "exact_source_polymer_match_to_reported": bool(
            row.get("exact_source_polymer_match_to_reported")
        ),
        "exact_analysis_entity_match_to_reported": bool(
            row.get("exact_entity_match_to_authoritative")
        ),
        "exact_observed_match_to_reported": bool(
            row.get("exact_observed_match_to_authoritative")
        ),
        "source_entity_description": str(
            row.get("source_entity_description", "")
        ),
        "source_polymer_sequence_provenance": str(
            row.get("source_polymer_sequence_provenance", "absent")
        ),
    }
    return json.dumps(evidence, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _color_histogram(
    value: object, field: str
) -> dict[tuple[int, int, int, int], int]:
    try:
        histogram = json.loads(str(value))
    except json.JSONDecodeError as exc:
        raise BaselineReviewError(f"invalid color histogram {field}") from exc
    if not isinstance(histogram, dict):
        raise BaselineReviewError(f"color histogram {field} is not an object")
    normalized: dict[tuple[int, int, int, int], int] = {}
    for key, raw_count in histogram.items():
        rgba = _rgba_tuple(key)
        try:
            count = int(raw_count)
        except (TypeError, ValueError) as exc:
            raise BaselineReviewError(
                f"color histogram {field} has a non-integer count"
            ) from exc
        normalized[rgba] = normalized.get(rgba, 0) + count
    return normalized


def _rgba_tuple(
    value: object, *, allow_empty: bool = False
) -> tuple[int, int, int, int] | tuple[()]:
    text = str(value).strip()
    if allow_empty and not text:
        return ()
    try:
        result = tuple(int(part.strip()) for part in text.split(","))
    except ValueError as exc:
        raise BaselineReviewError(f"invalid RGBA value {value!r}") from exc
    if len(result) != 4 or any(channel < 0 or channel > 255 for channel in result):
        raise BaselineReviewError(f"invalid RGBA value {value!r}")
    return result  # type: ignore[return-value]


def _integer_channel_tuple(
    value: object, length: int, label: str
) -> tuple[int, ...]:
    if not isinstance(value, list) or len(value) != length:
        raise BaselineReviewError(
            f"{label} must contain {length} integer channels"
        )
    try:
        result = tuple(int(item) for item in value)
    except (TypeError, ValueError) as exc:
        raise BaselineReviewError(f"{label} contains a non-integer channel") from exc
    if any(item < 0 or item > 255 for item in result):
        raise BaselineReviewError(f"{label} channel is outside 0..255")
    return result


def _read_json_object(path: Path) -> dict[str, object]:
    lexical = Path(path).expanduser().absolute()
    if lexical.is_symlink() or not lexical.is_file():
        raise BaselineReviewError(
            f"input must be a regular non-symlink file: {lexical}"
        )
    try:
        value = json.loads(lexical.resolve(strict=True).read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise BaselineReviewError(f"invalid JSON: {lexical}") from exc
    if not isinstance(value, dict):
        raise BaselineReviewError(f"JSON top-level is not an object: {lexical}")
    return value


def _clean_code(value: object) -> str:
    text = str(value or "")
    return "" if text in {" ", "\x00", ".", "?"} else text
