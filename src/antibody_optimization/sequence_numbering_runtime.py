"""Pinned ANARCII runtime adapter for the sequence-numbering baseline.

This module is the only executable bridge to ANARCII.  Semantic normalization
and validation remain in :mod:`antibody_optimization.sequence_numbering`, so
tests can exercise them without loading a model.
"""

from __future__ import annotations

from collections.abc import Sequence
from importlib import metadata as importlib_metadata

from .sequence_numbering import (
    ANARCII_VERSION,
    InputSequence,
    NumberingAudit,
    SequenceNumberingError,
    build_numbering_audits,
)


def run_anarcii_numbering(
    records: Sequence[InputSequence],
) -> tuple[NumberingAudit, ...]:
    """Run only ANARCII 2.0.8 antibody/accuracy/IMGT on one CPU.

    The fixed invocation uses ``cpu=True``, ``ncpu=1``, ``batch_size=8``, and
    ``scfv=False``.  No fallback version, device, sequence type, model mode, or
    numbering scheme is selected automatically.
    """

    installed_version = importlib_metadata.version("anarcii")
    if installed_version != ANARCII_VERSION:
        raise RuntimeError(
            f"ANARCII {ANARCII_VERSION} is required; found {installed_version}"
        )
    from anarcii import Anarcii

    runner = Anarcii(
        seq_type="antibody",
        mode="accuracy",
        cpu=True,
        ncpu=1,
        batch_size=8,
    )
    runner.number(
        {record.sample_uid: record.sequence_raw for record in records},
        scfv=False,
    )
    raw_results = runner.to_scheme(scheme="imgt")
    if not isinstance(raw_results, dict):
        raise SequenceNumberingError(
            "ANARCII returned a serialized or non-mapping result for this 47-sequence run"
        )
    return build_numbering_audits(records, raw_results)
