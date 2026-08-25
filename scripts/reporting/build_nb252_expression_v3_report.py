from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from antibody_optimization.v3_report_data import load_v3_report_data  # noqa: E402
from antibody_optimization.v3_report_document import (  # noqa: E402
    build_v3_report_document,
    write_v3_report_manifest,
)
from antibody_optimization.v3_report_figures import render_v3_report_figures  # noqa: E402


DEFAULT_TEMPLATE = (
    ROOT
    / "docs/result_artifacts/weekly_report_result/"
    "report_2026_W34_nb252_expression_route/"
    "Nb252_BL21_expression_optimization_project_report.docx"
)
DEFAULT_OUTPUT_DIR = (
    ROOT
    / "docs/result_artifacts/weekly_report_result/"
    "Nb252_V3_expression_report"
)
DEFAULT_DOCX = DEFAULT_OUTPUT_DIR / "Nb252_BL21_expression_optimization_V3_project_report.docx"
DEFAULT_PDF = DEFAULT_OUTPUT_DIR / "Nb252_BL21_expression_optimization_V3_project_report.pdf"
DEFAULT_MANIFEST = DEFAULT_OUTPUT_DIR / "Nb252_BL21_expression_optimization_V3_report_manifest.json"
DEFAULT_FIGURES = DEFAULT_OUTPUT_DIR / "figures"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build the V3 Nb252 expression project report without modifying V2 artifacts."
    )
    parser.add_argument("--template", type=Path, default=DEFAULT_TEMPLATE)
    parser.add_argument("--output", type=Path, default=DEFAULT_DOCX)
    parser.add_argument("--figures", type=Path, default=DEFAULT_FIGURES)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument(
        "--bind_pdf",
        type=Path,
        default=None,
        help="Bind an already-rendered PDF to the report manifest.",
    )
    parser.add_argument(
        "--finalize_only",
        action="store_true",
        help="Do not rebuild DOCX/figures; bind existing DOCX and optional PDF into the manifest.",
    )
    return parser


def main() -> int:
    args = _parser().parse_args()
    data = load_v3_report_data(ROOT)
    template_hash_before = _sha256(args.template)
    if args.finalize_only:
        if not args.output.is_file():
            raise FileNotFoundError(f"Existing V3 report is missing: {args.output}")
        metadata = {
            "output_docx": args.output.as_posix(),
            "sha256": _sha256(args.output),
            "candidate_count": 30,
            "single_count": 15,
            "double_count": 15,
            "antifold_role": data["antifold_policy"],
        }
        figure_metadata = {
            path.stem: {"path": path.as_posix(), "sha256": _sha256(path)}
            for path in sorted(args.figures.glob("v3_*.png"))
        }
    else:
        figure_result = render_v3_report_figures(ROOT, args.figures)
        figure_paths = {
            key: Path(value[0])
            for key, value in figure_result.items()
            if key != "source_data_csv"
        }
        metadata = build_v3_report_document(
            template_docx=args.template,
            output_docx=args.output,
            data=data,
            report_figures=figure_paths,
        )
        figure_metadata = {
            key: {
                "png": Path(value[0]).as_posix(),
                "png_sha256": _sha256(Path(value[0])),
                "svg": Path(value[1]).as_posix(),
                "svg_sha256": _sha256(Path(value[1])),
            }
            for key, value in figure_result.items()
            if key != "source_data_csv"
        }
        source_csv = Path(figure_result["source_data_csv"])
        figure_metadata["source_data_csv"] = {
            "path": source_csv.as_posix(),
            "sha256": _sha256(source_csv),
        }
    if _sha256(args.template) != template_hash_before:
        raise RuntimeError("Historical V2 template changed during V3 report generation")
    write_v3_report_manifest(
        path=args.manifest,
        document_metadata=metadata,
        pdf_path=args.bind_pdf,
        template_docx=args.template,
        source_artifacts=data["source_artifacts"],
        figure_metadata=figure_metadata,
        repository_root=ROOT,
    )
    print(f"V3 report DOCX: {args.output}")
    print(f"V3 report manifest: {args.manifest}")
    if args.bind_pdf:
        print(f"Bound PDF: {args.bind_pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
