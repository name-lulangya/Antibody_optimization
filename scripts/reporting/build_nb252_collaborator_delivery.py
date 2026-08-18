"""Build the collaborator-facing Nb252 report and sequence delivery folder."""

from __future__ import annotations

import argparse
import csv
import shutil
import tempfile
from collections import Counter
from pathlib import Path


CATEGORY_LABELS = {
    "affinity_focused_single": "亲和力导向单突",
    "property_focused_single": "性质导向单突",
    "balanced_combination": "平衡组合",
    "affinity_supported_double": "亲和力支持双突",
}


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--final-candidates", required=True, type=Path)
    parser.add_argument("--wt-samples", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _validate_and_sort_candidates(rows: list[dict[str, str]], wt: str) -> list[dict[str, str]]:
    if len(rows) != 30:
        raise ValueError(f"Expected 30 final candidates, found {len(rows)}")
    rows = sorted(rows, key=lambda row: int(row["final_panel_order"]))
    if [int(row["final_panel_order"]) for row in rows] != list(range(1, 31)):
        raise ValueError("Final panel order must be exactly 1-30")

    sequences = set()
    for row in rows:
        sequence = row["sequence"]
        if len(sequence) != 128 or not sequence.endswith("SSGS"):
            raise ValueError(f"Invalid sequence boundary for {row['mutation_set']}")
        if sequence in sequences:
            raise ValueError(f"Duplicate final sequence for {row['mutation_set']}")
        sequences.add(sequence)
        if row["panel_category"] not in CATEGORY_LABELS:
            raise ValueError(f"Unknown panel category: {row['panel_category']}")
        expected = list(wt)
        for mutation in row["mutation_set"].split(";"):
            source, target = mutation[0], mutation[-1]
            position = int(mutation[1:-1])
            if wt[position - 1] != source:
                raise ValueError(f"WT mismatch for {mutation}")
            expected[position - 1] = target
        if "".join(expected) != sequence:
            raise ValueError(f"Mutation label and sequence disagree for {row['mutation_set']}")
    return rows


def _wrap_fasta(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[index : index + width] for index in range(0, len(sequence), width))


def _write_delivery(staging: Path, rows: list[dict[str, str]], wt: str, report: Path) -> None:
    shutil.copy2(report, staging / report.name)

    fasta_lines: list[str] = []
    catalog_lines = [
        "# Nb252最终30条候选序列",
        "",
        "> 以下序列均为计算优先、建议实验测试的候选，尚未经过实验验证。WT不计入30条候选。",
        "",
        "| 序号 | 样品名称 | 突变 | 类别 | 长度（aa） | 完整氨基酸序列 |",
        "|---:|---|---|---|---:|---|",
    ]
    for order, row in enumerate(rows, start=1):
        sample_name = f"Nb252-C{order:02d}"
        category = CATEGORY_LABELS[row["panel_category"]]
        fasta_lines.extend(
            [
                f">{sample_name} mutation={row['mutation_set']} category={row['panel_category']}",
                _wrap_fasta(row["sequence"]),
            ]
        )
        catalog_lines.append(
            f"| {order} | {sample_name} | {row['mutation_set']} | {category} | 128 | `{row['sequence']}` |"
        )

    (staging / "Nb252_final_30_sequences.fasta").write_text(
        "\n".join(fasta_lines) + "\n", encoding="utf-8", newline="\n"
    )
    (staging / "Nb252_final_30_sequence_catalog.md").write_text(
        "\n".join(catalog_lines) + "\n", encoding="utf-8", newline="\n"
    )
    (staging / "Nb252_WT_control.fasta").write_text(
        f">Nb252-WT parent_control\n{_wrap_fasta(wt)}\n", encoding="utf-8", newline="\n"
    )

    counts = Counter(row["panel_category"] for row in rows)
    readme = f"""# Nb252合作者交付包

本文件夹用于阶段性沟通和后续实验订购，包含项目报告、最终30条候选序列和WT母本对照。

## 文件说明

- `{report.name}`：项目阶段报告，介绍数据、计算流程、筛选逻辑和主要结论。
- `Nb252_final_30_sequences.fasta`：30条候选的标准FASTA；样品名为`Nb252-C01`至`Nb252-C30`。
- `Nb252_final_30_sequence_catalog.md`：便于人工查看和复制的候选序列表。
- `Nb252_WT_control.fasta`：完整128-aa Nb252母本，建议作为独立实验对照，不计入30条候选。

## 候选构成

- 亲和力导向单突：{counts['affinity_focused_single']}条
- 性质导向单突：{counts['property_focused_single']}条
- 平衡组合：{counts['balanced_combination']}条
- 亲和力支持双突：{counts['affinity_supported_double']}条

## 使用说明

- 突变编号以完整128-aa Nb252母本序列的1-based位置为准。
- 所有候选均保留母本末端`SSGS`，序列长度均为128 aa。
- 候选类别表示计算设计侧重点，不代表已经实验测得亲和力、表达量或稳定性改善。
- 所有30条序列均建议在相同构建体、表达体系和检测条件下与WT平行验证。
"""
    (staging / "README.md").write_text(readme, encoding="utf-8", newline="\n")


def main() -> int:
    args = _parse_args()
    for path in (args.final_candidates, args.wt_samples, args.report):
        if not path.is_file():
            raise FileNotFoundError(path)
    if args.output_dir.exists():
        raise FileExistsError(f"Output directory already exists: {args.output_dir}")

    wt_rows = [row for row in _read_csv(args.wt_samples) if row["sample_uid"] == "LTT__Nb252"]
    if len(wt_rows) != 1:
        raise ValueError(f"Expected one LTT__Nb252 row, found {len(wt_rows)}")
    wt = wt_rows[0]["sequence_raw"]
    if len(wt) != 128 or not wt.endswith("SSGS"):
        raise ValueError("Authoritative Nb252 WT must be 128 aa and end in SSGS")

    rows = _validate_and_sort_candidates(_read_csv(args.final_candidates), wt)
    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{args.output_dir.name}-stage-", dir=args.output_dir.parent))
    try:
        _write_delivery(staging, rows, wt, args.report)
        args.output_dir.mkdir()
        for staged_file in staging.iterdir():
            shutil.copy2(staged_file, args.output_dir / staged_file.name)
        shutil.rmtree(staging)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        shutil.rmtree(args.output_dir, ignore_errors=True)
        raise
    print(f"Created collaborator delivery with {len(rows)} candidates: {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
