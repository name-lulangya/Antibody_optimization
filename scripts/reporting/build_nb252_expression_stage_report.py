from __future__ import annotations

import csv
import json
import shutil
from collections import Counter
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager
from docx import Document
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Inches, Pt, RGBColor


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs/result_artifacts/weekly_report_result/report_2026_W34_nb252_expression_route"
FIG = OUT / "figures"
TABLE = OUT / "tables"
DELIVERY = OUT / "delivery"

SELECTION_DIR = ROOT / "docs/result_artifacts/candidate_design/expression_single_mutant_trial_selection_v2_20260820"
CONSERVATION_DIR = ROOT / "docs/result_artifacts/input_baseline/vhh_conservation_consensus_v2_20260819"
LOGO_DIR = ROOT / "docs/result_artifacts/input_baseline/vhh_sequence_logos_20260818"
LANDSCAPE_DIR = ROOT / "docs/result_artifacts/candidate_design/expression_single_mutant_landscape_v1_20260820"

REPORT = OUT / "Nb252_BL21_expression_optimization_stage_report_2026_W31_W34.docx"
GUIDE = OUT / "Nb252_predictor_and_selection_guide.docx"

BLACK = "000000"
GRAY = "555555"
LIGHT = "E8E8E8"
MID = "B7B7B7"
FONT = "Microsoft YaHei"


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv_bom(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _set_run_font(run, size: float | None = None, bold: bool | None = None) -> None:
    run.font.name = FONT
    run._element.get_or_add_rPr().rFonts.set(qn("w:eastAsia"), FONT)
    run._element.get_or_add_rPr().rFonts.set(qn("w:ascii"), FONT)
    run._element.get_or_add_rPr().rFonts.set(qn("w:hAnsi"), FONT)
    run.font.color.rgb = RGBColor(0, 0, 0)
    if size is not None:
        run.font.size = Pt(size)
    if bold is not None:
        run.bold = bold


def _style_document(doc: Document) -> None:
    section = doc.sections[0]
    section.top_margin = Cm(1.9)
    section.bottom_margin = Cm(1.7)
    section.left_margin = Cm(2.1)
    section.right_margin = Cm(2.1)
    for name, size, bold in [
        ("Normal", 10.5, False),
        ("Title", 21, True),
        ("Subtitle", 11.5, False),
        ("Heading 1", 15, True),
        ("Heading 2", 12.5, True),
        ("Caption", 9, False),
    ]:
        style = doc.styles[name]
        style.font.name = FONT
        style._element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
        style.font.color.rgb = RGBColor(0, 0, 0)
        style.font.size = Pt(size)
        style.font.bold = bold
    normal = doc.styles["Normal"].paragraph_format
    normal.line_spacing = 1.25
    normal.space_after = Pt(5)
    doc.styles["Heading 1"].paragraph_format.space_before = Pt(4)
    doc.styles["Heading 1"].paragraph_format.space_after = Pt(7)
    doc.styles["Heading 2"].paragraph_format.space_before = Pt(3)
    doc.styles["Heading 2"].paragraph_format.space_after = Pt(5)


def _set_cell_text(cell, text: object, *, bold: bool = False, size: float = 8.5) -> None:
    cell.text = ""
    p = cell.paragraphs[0]
    p.paragraph_format.space_after = Pt(0)
    run = p.add_run(str(text))
    _set_run_font(run, size=size, bold=bold)
    cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER


def _shade(cell, color: str) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = tc_pr.find(qn("w:shd"))
    if shd is None:
        shd = OxmlElement("w:shd")
        tc_pr.append(shd)
    shd.set(qn("w:fill"), color)


def _add_table(doc: Document, headers: list[str], rows: list[list[object]], widths: list[float] | None = None, size: float = 8.5):
    table = doc.add_table(rows=1, cols=len(headers))
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = False
    for idx, header in enumerate(headers):
        _set_cell_text(table.rows[0].cells[idx], header, bold=True, size=size)
        _shade(table.rows[0].cells[idx], LIGHT)
    for row in rows:
        cells = table.add_row().cells
        for idx, value in enumerate(row):
            _set_cell_text(cells[idx], value, size=size)
    if widths:
        for row in table.rows:
            for idx, width in enumerate(widths):
                row.cells[idx].width = Inches(width)
    return table


def _add_bullet(doc: Document, text: str, level: int = 0) -> None:
    p = doc.add_paragraph(style="List Bullet" if level == 0 else "List Bullet 2")
    p.paragraph_format.space_after = Pt(3)
    run = p.add_run(text)
    _set_run_font(run, size=10.2)


def _add_caption(doc: Document, text: str) -> None:
    p = doc.add_paragraph(style="Caption")
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_before = Pt(3)
    p.paragraph_format.space_after = Pt(5)
    run = p.add_run(text)
    _set_run_font(run, size=8.8)


def _add_figure(doc: Document, path: Path, caption: str, width: float = 6.6) -> None:
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_after = Pt(0)
    p.add_run().add_picture(str(path), width=Inches(width))
    _add_caption(doc, caption)


def _fmt(value: str, digits: int = 3) -> str:
    return f"{float(value):.{digits}f}"


def _copy_source_figures() -> dict[str, Path]:
    FIG.mkdir(parents=True, exist_ok=True)
    sources = {
        "identity": ROOT / "docs/result_artifacts/input_baseline/summary/input_baseline_qc.png",
        "constraint_tracks": CONSERVATION_DIR / "nb252_conservation_constraint_tracks.png",
        "natural_logo": LOGO_DIR / "nb252_neighbor_sequence_logo_with_regions.png",
        "project_logo": LOGO_DIR / "project_expression_vhh_sequence_logo_with_regions.png",
        "landscape": LANDSCAPE_DIR / "expression_single_mutant_landscape.png",
        "scatter": LANDSCAPE_DIR / "expression_single_mutant_scatter.png",
        "selection": SELECTION_DIR / "expression_single_mutant_selection.png",
    }
    copied: dict[str, Path] = {}
    for key, source in sources.items():
        if not source.is_file():
            raise FileNotFoundError(source)
        target = FIG / f"source_{key}.png"
        shutil.copy2(source, target)
        copied[key] = target
    return copied


def _configure_plot_font() -> None:
    candidates = [Path(r"C:\Windows\Fonts\msyh.ttc"), Path(r"C:\Windows\Fonts\simhei.ttf")]
    font_path = next((p for p in candidates if p.exists()), None)
    if font_path:
        font_manager.fontManager.addfont(str(font_path))
        plt.rcParams["font.family"] = font_manager.FontProperties(fname=str(font_path)).get_name()
    plt.rcParams["axes.unicode_minus"] = False


def _make_route_figure() -> Path:
    _configure_plot_font()
    path = FIG / "figure_current_route.png"
    fig, ax = plt.subplots(figsize=(11.5, 2.6), dpi=220)
    ax.axis("off")
    labels = [
        "权威128-aa亲本\n+ 24界面位点冻结",
        "天然VHH保守性\n+ 80位置不可常规突变",
        "847条约束单突\n+ 单突且保留SSGS",
        "NetSolP / NanoMelt\n+ AntiFold统一评价",
        "幅度分档\n+ 风险与位点多样性",
        "30条计算试选\n+ 等待BL21实验",
    ]
    xs = [0.08, 0.25, 0.42, 0.59, 0.76, 0.93]
    for idx, (x, label) in enumerate(zip(xs, labels)):
        ax.text(x, 0.52, label, ha="center", va="center", fontsize=9.5,
                bbox=dict(boxstyle="round,pad=0.55", facecolor="white", edgecolor="black", linewidth=1.1))
        if idx < len(xs) - 1:
            ax.annotate("", xy=(xs[idx + 1] - 0.075, 0.52), xytext=(x + 0.075, 0.52),
                        arrowprops=dict(arrowstyle="->", color="black", lw=1.2))
    ax.text(0.5, 0.06, "当前路线不使用Rosetta/Flex ddG，不显式优化亲和力；组合突变等待单突实验结果后再讨论。",
            ha="center", va="center", fontsize=9.2)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def _make_funnel_figure() -> Path:
    _configure_plot_font()
    path = FIG / "figure_selection_funnel.png"
    labels = ["约束单突空间", "幅度短名单", "规则合格", "计算试选"]
    values = [847, 40, 37, 30]
    fig, ax = plt.subplots(figsize=(7.5, 3.5), dpi=220)
    y = list(range(len(labels)))
    ax.barh(y, values, color=["#8C8C8C", "#686868", "#4A4A4A", "#222222"])
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("候选数")
    ax.set_title("847条单突到30条计算试选的筛选漏斗")
    for yi, value in zip(y, values):
        ax.text(value + 8, yi, str(value), va="center", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlim(0, 910)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def _trial_reason(row: dict[str, str]) -> str:
    tier = row["selection_tier"]
    mapping = {
        "A_multi_family": "至少两个预测家族达到中等/强改善，且无中等以上恶化",
        "B_single_family_strong": "一个预测家族强改善，其余指标未出现明确恶化",
        "C_single_family_moderate": "一个预测家族中等改善，其余指标整体可接受",
        "D_controlled_tradeoff": "存在单项中等恶化，但有独立中等改善并保留机制多样性",
        "E_stable_word_exploratory": "稳定词新增假设对照；四项性质仅弱/近中性证据",
    }
    return mapping[tier]


def _short_mutation_label(row: dict[str, str]) -> str:
    return row["mutation_reported_label"].split()[-1]


def _candidate_summary_rows(rows: list[dict[str, str]]) -> list[list[object]]:
    output = []
    for row in rows:
        source = "实验复合物" if row["antifold_selection_source"] == "experimental_complex_context" else "AF3单体回退"
        output.append([
            row["trial_selection_order"], _short_mutation_label(row), row["region"],
            row["selection_tier"].split("_", 1)[0], source, _trial_reason(row),
        ])
    return output


def _candidate_metric_rows(rows: list[dict[str, str]]) -> list[list[object]]:
    output = []
    for row in rows:
        output.append([
            _short_mutation_label(row), row["region"],
            _fmt(row["netsolp_delta_usability_vs_current_wt"]),
            _fmt(row["netsolp_delta_solubility_vs_current_wt"]),
            _fmt(row["nanomelt_delta_predicted_apparent_tm_c_vs_current_wt"], 2),
            _fmt(row["antifold_selection_delta_log_probability"]),
            row["stable_word_effect"], row["selection_tier"].split("_", 1)[0],
        ])
    return output


def _add_page_title(doc: Document, title: str, subtitle: str | None = None) -> None:
    p = doc.add_paragraph(style="Title")
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run(title)
    _set_run_font(r, size=21, bold=True)
    if subtitle:
        p2 = doc.add_paragraph(style="Subtitle")
        p2.alignment = WD_ALIGN_PARAGRAPH.CENTER
        r2 = p2.add_run(subtitle)
        _set_run_font(r2, size=11.5)


def _add_footer(doc: Document) -> None:
    for section in doc.sections:
        footer = section.footer
        p = footer.paragraphs[0]
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        run = p.add_run("Nb252 BL21表达优化阶段报告｜计算结果仅用于实验优先级排序")
        _set_run_font(run, size=8)


def _build_report(trial: list[dict[str, str]], reserve: list[dict[str, str]], figures: dict[str, Path]) -> None:
    doc = Document()
    _style_document(doc)
    _add_page_title(doc, "Nb252纳米抗体BL21表达优化阶段报告", "报告期：2026-W31–W34｜面向导师与合作者")
    doc.add_paragraph()
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("当前版本聚焦表达量优化：冻结实验界面、避开天然保守位点、仅交付单突假设。")
    _set_run_font(r, size=13, bold=True)
    doc.add_paragraph()
    _add_figure(doc, figures["route"], "图1  当前项目路线。亲和力/Rosetta路线已经退出当前候选生成与排序。", 6.8)
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("版本日期：2026-08-20")
    _set_run_font(r, size=10)

    doc.add_page_break()
    doc.add_heading("执行摘要", level=1)
    _add_bullet(doc, "目标：在不突变实验界面、关键保守位点、二硫键Cys22/Cys95及末端SSGS的前提下，提出30条Nb252单突序列用于BL21表达实验。")
    _add_bullet(doc, "设计空间：天然VHH保守性和结构保护规则联合冻结80个reported positions；47个常规可扫位置产生846条非Cys单突，另加Q5V共识回变，共847条。")
    _add_bullet(doc, "工具：NetSolP提供S/U，NanoMelt提供预测表观Tm，AntiFold提供结构条件下的逐位点氨基酸兼容性。三者均是预测信号，不等同于实测产量。")
    _add_bullet(doc, "筛选：原始小数保留用于追溯，但候选选择只使用预定义幅度档位。847条中40条进入幅度短名单，排除2条序列风险和2条多重中等恶化后得到36条规则候选；另加入T99F稳定词假设对照，形成37条可选池，最终试选30条、候补7条。")
    _add_bullet(doc, "当前结论：30条仍是计算试选而非最终实验释放。F30（9条）、Q1（5条）和T27（4条）存在明显位点集中，需由项目负责人确认是否接受该探索密度。")
    _add_bullet(doc, "实验建议：与Nb252 WT同批表达，至少3个生物重复；以mg/L产量为主终点，同时记录可溶/总表达比例及纯化回收。单突数据返回后再决定组合突变。")
    doc.add_heading("当前试选面板概览", level=2)
    _add_table(doc, ["项目", "结果"], [
        ["候选形式", "30条单氨基酸替换；均为128 aa并保留末端SSGS"],
        ["候选构成", "26条严格核心 + 3条受控权衡 + 1条稳定词探索（T99F）"],
        ["位点覆盖", "13个reported positions；CDR 19条，FR 11条"],
        ["AntiFold证据", "24条实验复合物视图；6条实验缺失坐标位置使用独立AF3 VHH-only回退"],
        ["释放状态", "待用户/导师复核，尚未称为最终实验面板"],
    ], widths=[1.45, 5.2], size=9)

    doc.add_page_break()
    doc.add_heading("1. 项目目标与路线调整", level=1)
    doc.add_paragraph("项目最初同时讨论亲和力、稳定性和表达量，但与导师讨论后，当前目标收敛为Nb252在BL21体系中的表达量优化。此前亲和力扫描、Rosetta/Flex ddG和双突组合工作只保留为历史技术探索，不再参与本轮候选的生成、过滤、排序或解释。")
    doc.add_heading("1.1 当前设计原则", level=2)
    for text in [
        "实验结构中定义的全部VHH界面残基冻结，不再以预测亲和力改善为理由开放。",
        "不把FR与CDR机械划分为表达区和结合区；只要不属于硬冻结集合，FR/CDR均可被评价。",
        "最终实验阶段先测试单突；任何组合突变均等待单突实验结果。",
        "工具信号以互补证据和风险约束使用，不以简单加权总分或微小小数差异驱动选择。",
    ]:
        _add_bullet(doc, text)
    _add_figure(doc, figures["route"], "图2  当前有效路线及阶段边界。", 6.8)

    doc.add_page_break()
    doc.add_heading("2. 权威输入与结构保护边界", level=1)
    doc.add_paragraph("设计亲本为合作者确认的128-aa Nb252完整reported sequence。末端SSGS不是linker，是构建体组成部分，reported positions 125–128在所有候选中保持不变。实验复合物NK2R–Nb252是结合姿态和界面保护的主要结构证据；AF3 VHH仅补充实验未建模区的预测视图。")
    _add_table(doc, ["输入/证据", "本项目用途", "不能据此声称"], [
        ["Nb252 128-aa亲本", "定义WT、候选序列与末端SSGS", "不能自行裁剪或改写构建体"],
        ["实验NK2R–Nb252复合物", "界面保护、实验坐标可评价范围", "界面残基不是能量热点排名"],
        ["AF3 VHH预测", "实验缺失区的独立结构兼容性参考", "不能替代实验CDR3或受体接触"],
        ["47条BL21产量序列", "验证预测指标是否与产量有关", "不能训练高容量、可外推的表达模型"],
    ], widths=[1.55, 2.5, 2.55], size=8.7)
    doc.add_heading("2.1 不可突变集合", level=2)
    _add_bullet(doc, "24个实验界面残基：由严格<4.0 Å聚合物重原子中心距离定义复现。")
    _add_bullet(doc, "54个Nb252自身等于天然优势残基的高置信保守位点。")
    _add_bullet(doc, "Cys22/Cys95、末端SSGS以及其他合同规定的不可变位置。")
    _add_bullet(doc, "Q5不称为亲本保守位点：天然共识为V，因此只开放Q5V共识回变。")
    _add_figure(doc, figures["identity"], "图3  序列、实验坐标覆盖、AF3坐标覆盖及实验界面位置基线。灰色代表实验坐标不可评价，不等同于序列缺失。", 6.8)

    doc.add_page_break()
    doc.add_heading("3. 天然VHH保守性与可变位置", level=1)
    doc.add_paragraph("保守性参考采用TNP论文公开的4,059条非冗余天然VHH序列。经ANARCII/IMGT编号与完整性审核，4,057条通过；先按完整IMGT域90% identity去冗余得到3,784个有效簇，再按framework-only identity和coverage定义Nb252邻域。邻域包含1,564条序列、1,532个有效簇。按簇等权，避免大量近重复序列放大某一克隆家族。")
    doc.add_heading("3.1 高置信保守位点的定义", level=2)
    _add_bullet(doc, "邻域优势残基频率≥0.90、位置覆盖率≥0.80、有效簇数≥50。")
    _add_bullet(doc, "全局与邻域优势残基一致，且全局优势频率≥0.80。")
    _add_bullet(doc, "只有当Nb252自身残基等于该共同优势残基时，才因天然保守性完全冻结。")
    _add_bullet(doc, "满足统计保守性但Nb252偏离共识的位置不任意开放；仅允许共识回变。")
    _add_figure(doc, figures["constraint_tracks"], "图4  Nb252天然保守性与设计约束轨道。硬保守、界面、二硫键和末端构成设计禁区。", 6.8)

    doc.add_page_break()
    doc.add_heading("3.2 序列Logo提供的两类参照", level=1)
    doc.add_paragraph("天然邻域Logo反映经过簇权重校正的Nb252邻近VHH序列规律；项目47条表达序列Logo反映现有样本集合内部的残基分布。后者不能替代天然保守性合同，因为样本来源、构建体范围和序列完整性并不完全一致，且样本量较小。")
    _add_figure(doc, figures["natural_logo"], "图5  Nb252天然邻域序列Logo，按IMGT FR/CDR区域标注。", 6.7)
    _add_figure(doc, figures["project_logo"], "图6  项目47条表达序列中的可编号H链Logo，仅用于描述项目样本序列空间。", 6.7)

    doc.add_page_break()
    doc.add_heading("4. 847条约束单突空间", level=1)
    doc.add_paragraph("约束合并后，128个reported positions中80个不可常规突变。其余47个位置对除Cys外的19种替换进行完整枚举，共846条；加上Q5V共识回变，总计847条。每条候选均为唯一单突、长度128 aa、保留末端SSGS、不引入新Cys。")
    _add_table(doc, ["步骤", "数量", "说明"], [
        ["reported positions", 128, "权威亲本序列索引"],
        ["联合冻结位置", 80, "界面、天然保守、二硫键、末端及合同保护"],
        ["常规可扫描位置", 47, "每个位置19种非Cys替换"],
        ["常规单突", 846, "47×19，全部为单突"],
        ["共识回变", 1, "Q5V，仅此一条"],
        ["完整评价空间", 847, "进入统一性质计算，不是最终候选数"],
    ], widths=[1.7, 1.0, 3.9], size=9)
    doc.add_heading("4.1 为什么不先按FR/CDR切割", level=2)
    doc.add_paragraph("表达量可能同时受表面电荷、疏水性、局部构象稳定和折叠兼容性影响，这些因素跨越FR/CDR边界。因此本轮只使用结构界面和保守性作为硬约束；其余可变位置统一进入同一性质评价框架。")
    _add_figure(doc, figures["funnel"], "图7  完整约束空间到计算试选的数量漏斗。", 6.3)

    doc.add_page_break()
    doc.add_heading("5. 工具选择与证据角色", level=1)
    doc.add_paragraph("47条序列的BL21产量数据用于判断工具指标是否具备项目内解释价值。验证同时查看连续关联和预设方向下的样本外分类表现；没有任何工具被证明可直接预测未测Nb252突变体的mg/L产量。")
    _add_table(doc, ["工具/指标", "含义与方向", "当前证据角色"], [
        ["NetSolP S", "预测solubility；越高通常越有利", "有限产量关联；作为可溶性软偏好"],
        ["NetSolP U", "模型定义的usability；0–1，越高通常越有利", "与S同属一个模型家族，不重复计票"],
        ["NanoMelt Tm", "预测表观熔解温度（°C）；越高通常代表更稳定", "作为稳定性约束，不称为产量排序器"],
        ["AntiFold ΔlogP", "结构条件下突变残基相对WT的对数概率变化；越高越兼容", "实验复合物视图优先；实验缺失位置用独立AF3 VHH-only回退"],
        ["稳定词", "简并字母表下新增连续片段", "未验证为产量指标，仅保留一个探索对照T99F"],
    ], widths=[1.25, 2.55, 2.8], size=8.6)
    doc.add_heading("5.1 已测试但未进入当前排序的工具", level=2)
    doc.add_paragraph("RP3Net与PLM_Sol已在47条数据上完成验证，但未达到预先设定的独立、稳定、样本外信息要求；TNP和nanoBERT也不在当前精简工具集。它们不参与本轮847条候选的排序或候选解释。")
    doc.add_heading("5.2 幅度分档而非小数排序", level=2)
    _add_bullet(doc, "保留所有原始S/U/Tm/AntiFold数值，但先映射为强有利、中等有利、弱/近中性、中等不利、强不利。")
    _add_bullet(doc, "同一档内的小数差异不用于决定先后，避免模型噪声和微小波动主导候选。")
    _add_bullet(doc, "NetSolP S与U计作一个预测家族；AntiFold与NanoMelt分别代表结构兼容性与稳定性。")

    doc.add_page_break()
    doc.add_heading("6. 847条候选的统一性质评价", level=1)
    doc.add_paragraph("NetSolP和NanoMelt覆盖全部847条。AntiFold中，721条候选所在位置在实验VHH复合物坐标中可评价，优先使用实验复合物视图；126条位于实验结构缺失区，保持实验视图not_evaluable，并单独使用AF3 VHH-only结果作为回退。两种来源在图表和候选清单中明确区分。")
    _add_figure(doc, figures["landscape"], "图8  847条约束单突的四指标热图。AntiFold面板中126条实验坐标不可评价项以AF3 VHH-only结果补充并保留来源标记。", 6.9)
    doc.add_paragraph("热图显示，大多数单突的NetSolP S/U和NanoMelt Tm变化很小。当前流程因此不把极小变化解释为真实改善；只有达到预定义中等或强幅度时才形成主要入选证据。")

    doc.add_page_break()
    doc.add_heading("6.1 指标之间的关系", level=1)
    doc.add_paragraph("NetSolP S/U、NanoMelt Tm和AntiFold ΔlogP描述不同方面，并不要求同步变化。散点图用于识别明显权衡、孤立极端值和实验坐标缺失位置，而不是拟合一条统一最优线。")
    _add_figure(doc, figures["scatter"], "图9  四指标关键组合散点图。点形区分AntiFold实验复合物视图与AF3-only回退；星号标记稳定词新增。", 6.8)
    doc.add_heading("6.2 稳定词作为探索性软偏好", level=2)
    doc.add_paragraph("1,336条稳定词按固定12符号简并字母表进行大小写敏感、允许重叠的精确连续子串匹配。847条单突中22条新增24个命中，825条不变，没有减少项。但在47条产量数据中，稳定词特征未显示稳定正相关或可靠分类能力，因此不作为硬筛选器，只保留T99F用于实验检验该假设。")

    doc.add_page_break()
    doc.add_heading("7. 从847条到30条：幅度分档筛选", level=1)
    doc.add_paragraph("筛选先看幅度和风险，再看位点与机制多样性。847条中40条满足“无明显恶化且至少一个预测家族达到中等改善”。其中2条引入明确序列风险，2条同时出现两个中等不利指标，均不进入可释放池。剩余36条由26条严格核心和10条受控权衡组成；另将不在40条短名单中的T99F作为稳定词探索例外加入，形成37条可供试选，最终选择30条、保留7条候补。")
    _add_table(doc, ["阶段", "数量", "判定"], [
        ["完整约束空间", 847, "全部满足序列与硬约束"],
        ["幅度短名单", 40, "至少一个中等/强改善，且无明显恶化"],
        ["规则合格", 36, "排除2条序列风险和2条多重中等不利"],
        ["探索例外", 1, "T99F新增稳定词；仅弱/近中性性质证据"],
        ["计算试选", 30, "26严格核心 + 3受控权衡 + T99F"],
        ["候补", 7, "受控权衡备选"],
    ], widths=[1.4, 0.8, 4.5], size=9)
    _add_figure(doc, figures["selection"], "图10  当前30条计算试选：漏斗、位点分布和四指标幅度档位。空心三角表示AntiFold采用AF3-only回退；金色星号表示T99F稳定词探索候选。", 6.9)

    doc.add_page_break()
    doc.add_heading("8. 当前30条计算试选", level=1)
    doc.add_paragraph("30条覆盖13个位置，分为A–E五类。A/B/C为严格核心，D为受控权衡，E为稳定词探索。下表给出简要入选逻辑；完整序列与原始预测值见交付FASTA、CSV和附录。")
    for start in range(0, 30, 10):
        if start:
            doc.add_page_break()
            doc.add_heading(f"8. 当前30条计算试选（续：{start + 1}–{min(start + 10, 30)}）", level=1)
        _add_table(doc, ["序", "单突", "区段", "层", "AntiFold来源", "主要入选理由"],
                   _candidate_summary_rows(trial[start:start + 10]),
                   widths=[0.38, 0.65, 0.55, 0.35, 1.0, 3.5], size=7.8)

    doc.add_page_break()
    doc.add_heading("9. 位点集中、风险与不确定性", level=1)
    _add_table(doc, ["观察", "当前解释", "需要的决策/实验"], [
        ["F30占9条", "该位置在NetSolP/AntiFold中出现较多可接受替换，规则筛选自然集中", "确认是否接受同一位点的广泛替换探索；若实验容量有限，可保留化学性质差异最大的子集"],
        ["Q1占5条", "N端替换主要由单一预测家族支持", "关注起始端加工及表达构建体上下文，实验中保持同一载体和标签"],
        ["T27占4条、F29占2条", "实验结构缺失，AntiFold使用AF3-only回退", "解释时降低结构证据等级；不可称为实验复合物支持"],
        ["T99F稳定词探索", "新增长度5稳定词，但四项性质无中等改善", "作为机制对照，不与A–C层等价"],
        ["工具与产量相关有限", "预测可减少明显风险，但不能保证mg/L提高", "WT同批、多重复、保留全部阴性/中性结果"],
    ], widths=[1.25, 2.55, 2.8], size=8.5)
    doc.add_heading("9.1 风险控制", level=2)
    _add_bullet(doc, "所有30条均未触碰界面、硬保守位置、Cys22/Cys95和末端SSGS，也未引入新Cys。")
    _add_bullet(doc, "筛选检查新Pro、局部疏水/电荷变化、N-糖基化、脱酰胺/异构化和氧化易感残基等序列风险。")
    _add_bullet(doc, "受控权衡候选最多允许单项中等不利；同时出现两个中等不利的候选不进入面板。")
    _add_bullet(doc, "AntiFold回退只补充结构兼容性参考，不把实验缺失区转化为实验可评价区。")

    doc.add_page_break()
    doc.add_heading("10. 下一步实验与决策计划", level=1)
    doc.add_heading("10.1 面板释放前复核", level=2)
    _add_bullet(doc, "确认30条中F30/Q1/T27的位点集中度是否符合本轮探索目的。")
    _add_bullet(doc, "确认3条受控权衡和T99F探索对照是否占用实验名额；必要时从7条候补中一对一替换。")
    _add_bullet(doc, "冻结最终FASTA版本后，不再因同档内小数差异调整顺序。")
    doc.add_heading("10.2 BL21表达实验", level=2)
    _add_bullet(doc, "WT与30条单突使用相同载体、信号肽、标签、诱导、培养、裂解、纯化和定量流程。")
    _add_bullet(doc, "建议至少3个独立生物重复；主终点为纯化后mg/L产量，并记录总表达、可溶上清、沉淀与纯化回收。")
    _add_bullet(doc, "预先定义成功标准，例如相对WT提高幅度、重复一致性和不低于WT的可溶比例。")
    _add_bullet(doc, "保留全部阴性、中性和失败结果，用于下一轮局部规则更新。")
    doc.add_heading("10.3 组合突变的阶段门", level=2)
    doc.add_paragraph("本轮不生成或提交组合突变。只有当两个或多个单突在相同实验条件下分别显示可重复改善，且无明显共同风险时，才进入组合设计；组合后需要重新评价NetSolP、NanoMelt和AntiFold，并重新检查界面/保守/末端合同。")
    _add_table(doc, ["阶段门", "当前状态"], [
        ["结构与界面身份", "通过"],
        ["天然VHH保守性合同", "通过"],
        ["847条单突性质矩阵", "通过"],
        ["30条计算试选", "已生成，待人工复核"],
        ["最终实验面板释放", "未通过"],
        ["组合突变设计", "阻断，等待单突实验"],
    ], widths=[2.6, 4.0], size=9)

    doc.add_page_break()
    doc.add_heading("附录A 30条试选的原始变化值", level=1)
    doc.add_paragraph("下表保留原始模型变化值用于追溯。候选选择采用幅度档位，表内小数不用于同档内精细排序。AntiFold ΔlogP为突变残基相对WT的对数概率变化；正值为预测兼容性提高。")
    for start in range(0, 30, 10):
        if start:
            doc.add_page_break()
            doc.add_heading(f"附录A（续：{start + 1}–{min(start + 10, 30)}）", level=1)
        _add_table(doc, ["单突", "区段", "ΔU", "ΔS", "ΔTm °C", "AntiFold ΔlogP", "稳定词", "层"],
                   _candidate_metric_rows(trial[start:start + 10]),
                   widths=[0.65, 0.55, 0.55, 0.55, 0.65, 1.0, 0.8, 0.35], size=7.7)

    doc.add_page_break()
    doc.add_heading("附录B 7条候补", level=1)
    _add_table(doc, ["单突", "区段", "ΔU", "ΔS", "ΔTm °C", "AntiFold ΔlogP", "候补理由"], [
        [_short_mutation_label(row), row["region"], _fmt(row["netsolp_delta_usability_vs_current_wt"]),
         _fmt(row["netsolp_delta_solubility_vs_current_wt"]), _fmt(row["nanomelt_delta_predicted_apparent_tm_c_vs_current_wt"], 2),
         _fmt(row["antifold_selection_delta_log_probability"]), "受控权衡；用于替换位点集中或探索名额"]
        for row in reserve
    ], widths=[0.65, 0.55, 0.55, 0.55, 0.65, 1.0, 2.2], size=7.8)
    doc.add_heading("附录C 术语边界", level=1)
    _add_bullet(doc, "计算试选：通过当前计算规则、推荐进入人工复核的实验假设。")
    _add_bullet(doc, "最终实验面板：经项目负责人确认并冻结用于实际构建/表达的序列集合。当前尚未达到该状态。")
    _add_bullet(doc, "预测改善：模型输出方向有利，不等同于实测产量、溶解度、Tm或亲和力提高。")
    _add_bullet(doc, "AF3-only回退：实验坐标缺失位置的独立预测结构参考，不是实验复合物证据。")
    _add_bullet(doc, "稳定词探索：用户指定的可解释软偏好，尚未被47条产量数据验证。")

    _add_footer(doc)
    doc.core_properties.title = "Nb252纳米抗体BL21表达优化阶段报告"
    doc.core_properties.subject = "2026-W31–W34当前表达优化路线"
    doc.core_properties.author = "Antibody_optimization project"
    doc.core_properties.keywords = "Nb252, BL21, expression, NetSolP, NanoMelt, AntiFold"
    doc.save(REPORT)


def _build_guide() -> None:
    doc = Document()
    _style_document(doc)
    section = doc.sections[0]
    section.top_margin = Cm(1.25)
    section.bottom_margin = Cm(1.2)
    _add_page_title(doc, "Nb252候选指标与筛选说明", "一页速览｜用于阅读30条试选清单")
    _add_table(doc, ["指标", "方向", "在本项目中的正确解释"], [
        ["NetSolP U", "越高通常越有利", "0–1的usability模型输出；与S同属一个工具家族，不双重计票"],
        ["NetSolP S", "越高通常越有利", "预测solubility；与BL21产量仅有限相关，是软偏好"],
        ["NanoMelt Tm", "越高通常越稳定", "预测表观Tm（°C）；作为稳定性约束，不是mg/L预测"],
        ["AntiFold ΔlogP", "正值更兼容", "突变残基相对WT的结构条件对数概率变化；实验复合物优先，缺失区用AF3-only回退"],
        ["稳定词", "新增仅作探索", "简并片段新增；未验证为产量指标，仅T99F作为对照"],
    ], widths=[1.1, 1.1, 4.45], size=8.2)
    doc.add_heading("筛选逻辑", level=2)
    _add_bullet(doc, "先硬排除：24界面位置、天然硬保守位置、Cys22/Cys95、末端SSGS及其他不可变位置。")
    _add_bullet(doc, "再按幅度档位：强/中等有利、弱/近中性、中等/强不利；同档内小数不排序。")
    _add_bullet(doc, "A–C层为严格核心；D层为单项中等不利的受控权衡；E层T99F是稳定词假设对照。")
    _add_bullet(doc, "空心三角：AntiFold使用AF3 VHH-only回退；金色星号：稳定词探索候选。")
    _add_bullet(doc, "当前30条是计算试选，需人工确认位点集中后才能冻结为实验面板。")
    doc.add_heading("实验判读", level=2)
    p = doc.add_paragraph("所有预测均不是实验结果。建议WT与30条单突同批、至少3个生物重复，以纯化后mg/L为主终点，同时记录总表达、可溶比例和纯化回收；单突结果返回前不组合。")
    for run in p.runs:
        _set_run_font(run, size=9.2)
    doc.save(GUIDE)


def _write_delivery(trial: list[dict[str, str]], reserve: list[dict[str, str]]) -> None:
    DELIVERY.mkdir(parents=True, exist_ok=True)
    shutil.copy2(SELECTION_DIR / "expression_single_mutant_trial30.fasta", DELIVERY / "Nb252_trial30_single_mutants.fasta")
    shutil.copy2(SELECTION_DIR / "expression_single_mutant_reserve.fasta", DELIVERY / "Nb252_reserve7_single_mutants.fasta")
    fields = [
        "trial_selection_order", "mutation_reported_label", "reported_sequence_index_1based", "region",
        "imgt_position_label", "selection_tier", "trial_selection_reason", "sequence",
        "netsolp_predicted_usability", "netsolp_delta_usability_vs_current_wt",
        "netsolp_predicted_solubility", "netsolp_delta_solubility_vs_current_wt",
        "nanomelt_predicted_apparent_tm_c", "nanomelt_delta_predicted_apparent_tm_c_vs_current_wt",
        "antifold_selection_source", "antifold_selection_delta_log_probability", "stable_word_effect",
        "final_experimental_panel_released",
    ]
    _write_csv_bom(DELIVERY / "Nb252_trial30_single_mutants.csv", trial, fields)
    _write_csv_bom(DELIVERY / "Nb252_reserve7_single_mutants.csv", reserve, fields)
    readme = (
        "# Nb252 BL21表达优化阶段交付包\n\n"
        "本目录面向导师和合作者，包含当前表达优化路线的阶段报告、指标说明、30条计算试选和7条候补。\n\n"
        "- `Nb252_BL21_expression_optimization_stage_report_2026_W31_W34.docx/pdf`：完整阶段报告。\n"
        "- `Nb252_predictor_and_selection_guide.docx/pdf`：一页指标与筛选说明。\n"
        "- `Nb252_trial30_single_mutants.fasta/csv`：30条计算试选；尚未等同于最终实验释放。\n"
        "- `Nb252_reserve7_single_mutants.fasta/csv`：7条候补。\n\n"
        "所有候选均为128-aa单突并保留末端SSGS。预测值不是实测表达量、溶解度、Tm或亲和力。\n"
    )
    (DELIVERY / "README.md").write_text(readme, encoding="utf-8")


def main() -> int:
    if OUT.exists():
        raise FileExistsError(f"Output directory already exists: {OUT}")
    FIG.mkdir(parents=True)
    TABLE.mkdir(parents=True)
    trial = _read_csv(SELECTION_DIR / "expression_single_mutant_trial30.csv")
    reserve = _read_csv(SELECTION_DIR / "expression_single_mutant_reserve.csv")
    if len(trial) != 30 or len(reserve) != 7:
        raise ValueError(f"Unexpected selection counts: trial={len(trial)}, reserve={len(reserve)}")
    if any(row["final_experimental_panel_released"].lower() != "false" for row in trial):
        raise ValueError("Trial rows must remain unreleased")

    figures = _copy_source_figures()
    figures["route"] = _make_route_figure()
    figures["funnel"] = _make_funnel_figure()
    _build_report(trial, reserve, figures)
    _build_guide()
    _write_delivery(trial, reserve)

    position_counts = Counter(row["reported_sequence_index_1based"] for row in trial)
    manifest = {
        "schema_version": 1,
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "reporting_period": "2026-W31--W34",
        "active_route": "BL21_expression_only_single_mutants",
        "trial_count": 30,
        "reserve_count": 7,
        "unique_trial_positions": len(position_counts),
        "trial_position_counts": dict(sorted(position_counts.items(), key=lambda item: int(item[0]))),
        "final_experimental_panel_released": False,
        "report_docx": REPORT.name,
        "guide_docx": GUIDE.name,
        "delivery_directory": "delivery",
        "source_selection_release": "expression_single_mutant_trial_selection_v2_20260820",
        "excluded_current_route_methods": ["Rosetta", "PyRosetta", "Flex ddG", "affinity ranking", "double mutants"],
    }
    (OUT / "report_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    shutil.copy2(REPORT, DELIVERY / REPORT.name)
    shutil.copy2(GUIDE, DELIVERY / GUIDE.name)
    print(REPORT)
    print(GUIDE)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
