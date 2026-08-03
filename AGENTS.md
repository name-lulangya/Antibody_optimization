# Agent Rules for Antibody_optimization

## Project Goal

This project optimizes the experimentally identified NK2R nanobody Nb252 for multiple properties, including binding affinity, stability, and expression, while preserving scientific interpretability, structural provenance, and experimental testability.

## Project Scientific Context

- `Nb252-optimization.cxs` is the current structure session supplied by the collaborator. It contains three models:
  - `NK2R-252.pdb`: the experimentally resolved binding conformation of NK2R in complex with Nb252.
  - `NK2R-NKA.pdb`: the NK2R complex with its ligand NKA.
  - `fold_2r_252_nomg_model_0.cif`: the AlphaFold 3 (AF3) prediction of the Nb252 VHH structure.
- These three names currently identify models inside the ChimeraX session; they are not standalone PDB/mmCIF files in the repository unless they are explicitly exported and verified.
- Parts of the nanobody are not built in the experimental `NK2R-252.pdb` structure. Missing coordinates must not be interpreted as deleted residues or proof that the corresponding sequence is absent.
- The AF3 VHH prediction is provided to supplement the structural view of those unbuilt regions. According to the collaborator/user's visual ChimeraX inspection, alignment to the experimentally resolved nanobody is good outside CDR3, while CDR3 differs; this observation has not yet been quantitatively verified. It is a starting point for analysis, not a license to treat the two structures as interchangeable.
- The experimental NK2R–Nb252 complex is the primary structural evidence for the observed binding pose and interface. The AF3 model is a prediction and must remain labeled as such. The NK2R–NKA complex provides ligand-binding context; it is not an Nb252 complex or direct evidence for an Nb252 mutation effect.

## Non-Negotiable Rules

- Read this file and `docs/codex-handoff.md` before starting substantive work.
- Treat the repository as the primary code, documentation, Codex, and report-editing workspace. Store large datasets, checkpoints, exhaustive prediction tables, and generated experiment outputs only in locations explicitly established by project documentation or the user.
- The configured Git remote `origin` is `git@github.com:name-lulangya/Antibody_optimization.git`.
- The project is planned to be checked out under the remote-server parent directory `/homes/Tianlab/luly25/`. Do not assume the final repository root until the server checkout has been created and verified.
- The local Windows project Conda environment is named `ab_optim` and currently resolves to `D:\miniconda\envs\ab_optim`; it was created with `conda create -n ab_optim python=3.11`. Install all locally compatible project packages and tools into this environment by default; do not use `base` or unrelated local environments for project work.
- The project's primary remote conda environment is `/data/software/env/luly25/ab_optim`.
- The local named environment `ab_optim` and the remote path `/data/software/env/luly25/ab_optim` are distinct environments. Do not assume that packages installed in one are present in the other. Before a local installation, confirm Windows and Python 3.11 compatibility and dependency/licence requirements; do not force AF3, PyRosetta, nanoBERT, Linux-only binaries, licensed binaries, or incompatible CUDA/PyTorch stacks into the local environment when the documented server or tool-specific environment is the valid execution route.
- PyRosetta at `/data/software/env/luly25/multi_ligand` and nanoBERT at `/data/software/env/luly25/vhh-lm` are separate tool-specific environments, not the project environment. Activate them only for their corresponding tool stages, use the pinned versions, model revision, commands, and offline cache settings recorded in `docs/codex-handoff.md` and the current weekly history, and do not assume their packages are available in `ab_optim` or silently substitute another environment. The remote login alias is still not established.
- The remote Slurm resource policy is established in the Long-Task Script Rules below and matches the reference project. Do not deviate from that policy unless the user explicitly changes it for this project or a specific task.
- Do not run heavy training, exhaustive sequence design, full-library structure prediction, large docking/refinement campaigns, or long inference in an unconfirmed environment. First estimate resource needs and agree on an appropriate execution route.
- Do not edit generated data or results unless the task explicitly asks for it.
- Do not overwrite or revert user changes. Check `git status --short` before edits and before any commit.
- Preserve computational and scientific provenance. Any substantive implementation, design workflow change, structural preparation step, or result-producing workflow change must be documented in the current weekly history file.
- When planning or coding, if required information is missing, ambiguous, or raises a material question, stop and ask the user directly rather than proceeding silently on assumptions. This applies especially when different scientific choices could materially change residue mappings, candidate selection, structural conclusions, or experimental priorities; do not silently choose a consequential chain mapping, numbering convention, CDR definition, alignment, scoring protocol, or assay interpretation.
- When an implementation route is changed, replace the old route in the active script or workflow and remove obsolete branches, CLI options, compatibility code, and automatic fallbacks. Do not retain the superseded executable route as a fallback unless the user explicitly requests it. Historical experiment scripts may remain as provenance; active code may retain a brief rationale comment, not commented-out executable blocks.
- Never substitute mock, synthetic, fabricated, or placeholder sequences, structures, scores, assays, or metadata for required project inputs merely to make code, checks, smoke runs, integration runs, or result-producing workflows succeed. Deliberately synthetic fixtures are allowed only in isolated unit tests, must be explicitly named as test-only, and must never be selected automatically by runtime code.
- Never report a computational prediction as an experimental observation. Clearly label predicted, modeled, inferred, and experimentally measured values in code outputs, tables, figures, history, handoff documents, and conversation.

## Code Style

- Prefer small, conservative Python changes that follow nearby script style.
- Keep existing experiment scripts stable unless the task is explicitly to modernize or align them.
- Prefer explicit CLI arguments or config values for reusable scripts. Avoid new hardcoded paths, chain IDs, residue ranges, numbering offsets, CDR boundaries, or scoring thresholds.
- Use structured readers and writers for PDB/mmCIF, FASTA, CSV, parquet, HDF5, JSON, and pickle files instead of ad hoc string parsing when practical.
- Preserve identifiers that matter for structural traceability, including model name, source file, chain ID, residue name, author and label numbering where applicable, insertion code, alternate location, and atom selection.
- Keep multiprocessing and distributed code rerunnable and explicit about inputs, outputs, partial files, and overwrite behavior. Add restart or resume state only when the Long-Task Script Rules below require it.
- Add comments only for non-obvious logic, especially around residue mapping, structural alignment, CDR semantics, missing coordinates, scoring aggregation, candidate filtering, or distributed sharding.

## Script Scope and Complexity Rules

- Keep experiment entry points focused on orchestration: parse arguments, connect explicit inputs and outputs, call reusable utilities, and record run metadata. Move stable algorithms, shared validation, generic summaries, and lightweight I/O into `src/antibody_optimization/` when they have a clear reusable contract.
- Give each script one primary workflow responsibility, such as structure preparation, sequence/numbering reconciliation, interface analysis, candidate generation, property prediction, multi-objective ranking, statistical analysis, or artifact rendering. Split independent stages instead of concentrating computation, plotting, reporting, and unrelated design routes in one entry point.
- Treat script length as a review signal rather than a hard design target. When a Python script approaches or exceeds roughly 500 lines, inspect it for duplicated helpers, unused code, an overlong `main()`, or multiple responsibilities. Scripts above roughly 800 lines should normally be split, or the reason for keeping one cohesive file must be recorded in the current weekly history.
- Do not split code mechanically to satisfy a line target. Extract only logic with stable inputs and outputs that can be tested independently or reused; keep experiment-specific paths, mutation sets, parameter grids, and scientific orchestration in the experiment entry point.
- Experiment scripts must not import implementation code from other experiment scripts. Promote shared behavior to `src/antibody_optimization/`, add focused tests, and update `src/README.md` before reusing it across workflows.
- After a substantial script change, check for unused functions, duplicated helpers, obsolete branches, and unclear stage boundaries. Record material responsibility splits or reusable-code extractions in the current weekly history.

## Incremental Refactoring Rules

- Refactor gradually by extraction, not by moving old experiment scripts in bulk.
- Keep old scripts in place unless the user explicitly asks to remove or relocate them.
- Before adding a script, duplicating helper code, or extracting a utility, check `src/README.md` and search `src/antibody_optimization/` with `rg` for reusable functions.
- If a suitable utility already exists, prefer importing or extending it over copying logic into an experiment script.
- Add new reusable, stable logic under `src/antibody_optimization/` when duplicate logic appears across experiment folders.
- New scripts should prefer importing validated utilities from `src/antibody_optimization/` instead of copying logic from older experiment folders.
- Old scripts should be connected to `src` utilities only when fixing a bug, unifying current semantics, or performing a deliberate migration with small-sample verification.
- Keep `src/README.md` current. When adding, removing, renaming, or materially changing a `src/antibody_optimization/` utility, update the index in the same task.
- New `src` utilities must document:
  - the code logic and algorithmic assumptions;
  - all inputs and return values;
  - the intended scope of use and what is deliberately out of scope.
- Start with pure functions and low-risk helpers before extracting I/O-heavy, structure-modifying, or distributed logic.
- Do not change scientific semantics during a mechanical extraction. If semantics change, document it as a separate intentional behavior change.
- Prefer focused unit tests for every new reusable utility.
- Existing scripts should adopt `src` utilities only after the extracted behavior is verified against the old implementation on a small, real sample or an explicitly test-only fixture.

## Structural Data and Optimization Principles

### Structure identity and provenance

- Preserve the original collaborator-provided `Nb252-optimization.cxs` and its source models. Do not overwrite original coordinates when cleaning, aligning, renumbering, completing, mutating, minimizing, docking, or exporting structures.
- Every derived structure must have traceable provenance. Record at least the source file and model, software and version, command or script and parameters, chains and residues affected, alignment transform when relevant, modeling/refinement method, output path, and timestamp.
- Keep experimental structures, AF3 predictions, completed/hybrid models, minimized structures, docked poses, and mutation models visibly distinct in filenames, manifests, metadata, tables, and figures.
- Do not silently combine experimental and predicted coordinates. If a completed or hybrid model is created, record exactly which residues and atoms came from each source, how boundaries were joined, what optimization was performed, and which downstream analyses used the hybrid instead of the original experimental structure.
- Missing residues, unresolved side chains, alternate conformers, nonstandard residues, protonation states, glycans, ions, waters, and other hetero components must be handled explicitly. Do not delete or rebuild them without recording the decision and its consequence.
- Never transfer a contact, confidence, secondary-structure, CDR, or interface annotation from one model to another solely because the structures appear aligned. Transfer requires an explicit residue mapping and validation of the mapped sequence and coordinates.

### Chain identity, residue numbering, and sequence mapping

- Before residue-level analysis or mutation design, identify and record the chain IDs and molecular roles in every source model. Do not infer chain roles from chain order or visual position alone.
- Reconcile the Nb252 sequence across the experimental structure, the AF3 model, and any collaborator-provided or assay sequence. Record gaps, unresolved coordinates, terminal tags, construct boundaries, engineered mutations, modified residues, and sequence discrepancies.
- Retain both source numbering and any standardized antibody numbering. Mapping tables must include model/source, chain ID, source residue number, insertion code, residue identity, standardized position if used, and sequence index.
- Never describe a mutation with an ambiguous bare residue number. Include the molecule or chain, wild-type residue, position and numbering scheme, and mutant residue, for example in a documented form equivalent to `Nb252 chain X, IMGT 105, Y105F`.
- Validate wild-type residue identity against the authoritative sequence and the relevant structure before generating or interpreting a mutant.
- Do not renumber structures in place. Write a new derived file and a reversible mapping table.

### CDR definitions and structural alignment

- State the antibody numbering and CDR definition used for every CDR-dependent analysis. Acceptable schemes may include IMGT, Kabat, Chothia, AHo, or a tool-specific definition, but the chosen scheme, tool, version, input sequence, and resulting boundaries must be recorded.
- Do not assume that “CDR3” has identical residue boundaries under different schemes. When comparing reports or tools, reconcile the definitions explicitly.
- For every reported structural alignment, record the reference and mobile models, chain selections, residue/atom selections, sequence-alignment method if any, fitting algorithm, outlier-rejection rule, number of fitted residues/atoms, RMSD or other metric, software/version, and saved transform or reproducible command.
- Because the experimental and AF3 Nb252 structures differ in CDR3, do not fit on CDR3 by default when the goal is to compare framework orientation. Conversely, do not exclude CDR3 silently when the scientific question concerns CDR3 conformation or the binding interface.
- Report whether an RMSD is calculated over backbone, Cα, all heavy atoms, framework only, CDRs, interface residues, or another explicit atom set. RMSD values without selection provenance are not comparable.

### Interface interpretation and mutation design

- Use `NK2R-252.pdb` as the primary evidence for the experimentally observed NK2R–Nb252 binding pose and interface, subject to its unresolved regions and experimental-model limitations.
- In the `NK2R-252` model inside `Nb252-optimization.cxs`, the collaborator colored putative interface VHH residues orange and described them as being within 4 Å of the partner, but the distance definition and exact residue list have not been documented. Treat this as an unverified collaborator annotation: identify and remap the residues, reproduce the selection with an explicit atom/distance rule, and use extra caution before mutating this region.
- Use `fold_2r_252_nomg_model_0.cif` to examine the predicted VHH fold and unbuilt regions, while preserving AF3 confidence information when available. AF3 coordinates do not experimentally establish the CDR3 conformation or receptor contacts.
- Use `NK2R-NKA.pdb` to analyze the NKA-bound context, overlap, steric relationships, or possible mechanism only after NK2R chains and residue mappings are verified. Structural proximity or overlap alone does not prove competition, agonism, antagonism, or functional effect.
- Define interface contacts, buried surface area, hydrogen bonds, salt bridges, clashes, and energetic terms with explicit software, versions, parameters, protonation assumptions, and distance/angle cutoffs.
- Candidate mutation generation must preserve a traceable path from structural or sequence rationale to the exact candidate sequence. Record parental sequence, mutation set, numbering scheme, source structure/model, design method/version, constraints, and filtering decisions.
- Avoid treating a single structure, single predicted pose, or single scoring function as ground truth. When feasible, assess sensitivity to unresolved coordinates, alternative CDR3 conformations, model preparation, and scoring method.
- Flag mutations that may affect paratope geometry, framework packing, disulfide bonds, conserved antibody positions, charge patches, aggregation propensity, chemical liabilities, protease sensitivity, glycosylation motifs, or expression. Do not automatically filter them without documenting the rule and rationale.
- Preserve diversity in proposed experimental panels when multiple mechanistic hypotheses remain plausible. Do not collapse the design to a single top-scoring sequence without reporting uncertainty and alternatives.

### Multi-objective optimization and experimental validation

- Affinity, stability, expression, specificity, developability, and functional activity are distinct objectives. Do not claim global improvement from a gain in one property while ignoring measured or predicted tradeoffs in others.
- Keep raw property values, units, assay conditions, construct format, replicate counts, uncertainty, and data provenance. Do not combine measurements from different assay formats or conditions as directly comparable without an explicit normalization and caveat.
- Computational affinity or energy scores are model-specific ranking signals, not measured \(K_D\), \(k_{\mathrm{on}}\), \(k_{\mathrm{off}}\), potency, or functional activity. Do not convert them into physical assay values without a validated calibration.
- Multi-objective rankings must record objective definitions, direction, transformations, missing-value handling, constraints, weights or Pareto method, tie-breaking, and software/version. Never infer or reuse a universal weighting formula.
- Separate hard constraints from soft preferences. Preserve per-objective values in outputs so a composite rank never hides a severe regression.
- Compare candidates only under aligned computational and experimental protocols: same parental construct, numbering/mapping, structure-preparation route, evaluated residues, model versions, seeds or ensembles, assay format, controls, thresholds, and metric implementation. Disclose unavoidable differences and do not present them as a fully controlled superiority comparison.
- Candidate recommendations are hypotheses for experimental testing unless supported by actual experiments. Clearly distinguish “recommended for testing,” “computationally prioritized,” and “experimentally validated.”
- Experimental validation plans should include the appropriate parent/control constructs, replicate strategy, assay units and conditions, and predefined decision criteria. Do not claim optimization success until the relevant experimental measurements support it.
- Preserve negative, neutral, and failed experimental outcomes with the same provenance as positive results; do not retrospectively omit them from model assessment.

### Numerical result integrity

- Every numerical result reported in conversation, history, handoff, run summaries, tables, captions, or plots, and every value used to create a figure, must be read from or recomputed from concrete project data/result files or authoritative project records such as machine-readable run summaries, manifests, assay exports, and tracked artifacts. Never use conversational memory, manually transcribed narrative values, or an earlier caption as the sole data source.
- Before reporting or plotting derived results, reopen and verify source paths, schemas, sequence/structure identifiers, sample counts, units, assay conditions, replicate and aggregation semantics, model versions, seeds, thresholds, and filtering. Retain enough provenance to identify those sources.
- If a required source is unavailable, state that the value is unverified and do not present or plot it as an observed result.

## Git Rules

- Run `git status --short` before edits, before commits, and before final reporting.
- Do not commit or push unless the user explicitly asks.
- Do not create a standalone history entry only for a commit or push operation.
- If substantive changes are not yet documented, update the relevant current weekly history entry before committing.
- If the changes are already documented and the user only asks to commit or push them, update the existing relevant history entry's Git sync status section, usually titled `Git sync status`, instead of adding a new entry.
- Before committing documented changes, update the relevant existing history entry's Git sync status section with `ready to commit`, the intended commit message, and the target remote/branch.
- Keep code changes and matching documentation/history updates in the same commit when practical.
- Report the commit hash and push result in the final response.
- Never use destructive Git commands such as `git reset --hard` or checkout-based reverts unless explicitly requested.

## Long-Task Script Rules

- Before writing or materially changing a script, estimate expected runtime and identify both whether it is likely to exceed one hour and whether it is likely to exceed five hours.
- Require a pre-run `--check_only` and a remote real-data smoke run only for scripts expected to exceed five hours. For scripts expected to finish within five hours, use static checks and a clear run command; do not add or require a smoke stage by default unless the user explicitly requests one.
- For tasks expected to finish within five hours, do not implement checkpoint resume, partial-output recovery, or per-stage restart state by default. If such a task is interrupted, fix the underlying issue and rerun from the beginning with explicit overwrite behavior.
- Add resume or partial-output recovery only when the expected runtime exceeds five hours, the workflow contains an independently expensive completed stage that must be preserved, or the user explicitly requests it. Keep the recovery contract narrow and documented.
- Do not add speculative fallback paths or defensive recovery state machines merely to make a short task appear more robust. Prefer one explicit execution route so it remains clear which logic produced the result.
- Scripts expected to run longer than one hour must include `tqdm` progress bars around primary long-running loops when feasible. If `tqdm` is unsuitable, emit regular progress logs with comparable visibility.
- Scripts expected to run longer than one hour must be paired with a Slurm submission file and should run through Slurm on the remote Linux server.
- Slurm scripts must use `#SBATCH --partition=batch` by default. Do not use `#SBATCH --partition=gpu` unless the user explicitly requests that partition for a specific task.
- Slurm scripts must request at least one GPU because the current environment cannot request CPU-only jobs separately. Use the project resource convention that one GPU is paired with 12 CPUs. Do not explicitly add `#SBATCH --mem` by default; only specify memory when the user asks for it or a task has a clear non-default memory requirement.
- Multi-GPU jobs must run only on `n1` or `n2` and must not use `n3`. For single-node multi-GPU wrappers, request one node and add `#SBATCH --exclude=n3`; do not list both `n1,n2` with `--nodelist`, because that can request both nodes while the current `torchrun` launch is single-node.
- Before launching a long or resource-intensive task, confirm the remote login route, conda environment, actual resource availability, and input/output locations from concrete project documentation or the user. Do not invent hostnames or environment activation commands, and do not deviate from the established Slurm policy silently.
- Long-running scripts must expose clear input/output paths and explicit overwrite behavior; expose resume behavior only when the runtime and recovery rules above require it.
- Record enough metadata to reproduce a run: command, arguments, input structures/sequences, model and database versions, checkpoints where relevant, mapping/CDR/alignment definitions, scoring parameters, output directory, software environment, and timestamp.
- Prefer resumable partial outputs for large library prediction, structure generation, docking, refinement, or design jobs only when they are expected to exceed five hours or contain an independently expensive reusable stage.
- Emit compact progress and final statistics files for large jobs.
- Avoid interactive prompts in scripts intended for remote batch execution.

## Run Summary Handoff Rules

- Existing scripts do not need to be retrofitted only to add git-tracked summaries. Future result-producing scripts should write a lightweight git-tracked run summary unless the task explicitly says not to.
- Apply these rules to old scripts only when they are otherwise being updated.
- Full experiment outputs should remain under established ignored result/output directories. Git-tracked run summaries are the lightweight handoff layer for making essential result facts visible without tracking the full outputs.
- Prefer one summary folder per experiment under `docs/run_summaries/`, for example `docs/run_summaries/candidate_design/<experiment>/`.
- Each major script or experiment stage should write one small summary file in that experiment summary folder when the run finishes.
- Summary files may use Markdown, CSV, TSV, JSON, or another lightweight text format. Choose the format that makes the recorded statistics easiest to inspect.
- Summary content should record result facts only: run time, script/command, key parameters, source structure and sequence identifiers, chain/numbering/CDR/alignment definitions, input/output paths, model/checkpoint/database versions, simple statistics, property values, compact tables, and failure messages if any.
- Run summaries may include moderately detailed result statistics when they remain readable and git-trackable, such as elapsed time, sequence or candidate counts, pass/fail counts, score distributions, interface counts, coverage, cache/shard statistics, and compact per-objective tables.
- Prefer compact fields, short Markdown tables, or small CSV/TSV/JSON summaries over dumping large result files into Git-tracked documentation.
- Do not use run summaries for interpretation, next-step planning, narrative analysis, or daily history. They are only an intermediate transfer record for observed or computed results.
- When synchronized summaries are folded into `docs/history/YYYY-Www.md`, copy only the key information needed for project history, not every table or low-level metric.
- If a synchronized summary omits useful but available small metadata, inspect the concrete source files before concluding that the information is unavailable.
- Result-producing scripts must not automatically commit or push summary files. Commit, push, branch creation, and merge/cherry-pick remain explicit user-controlled steps.

## Git-Tracked Result Artifact Rules

- In addition to run summaries, result-producing workflows must write small Git-tracked result artifacts for important results and may add further compact artifacts when they materially improve review, plotting, or agent analysis.
- Every experiment's important result must be presented to the user in at least one Git-tracked figure under `docs/result_artifacts/`; it must not remain implicit only in CSV/JSON/parquet files, run summaries, or ignored outputs. Retain the figure, executable plotting script, and exact compact source data so the result can be independently reviewed and reproduced.
- Store artifacts under `docs/result_artifacts/`, mirroring the experiment hierarchy where practical, for example `docs/result_artifacts/candidate_design/<route>/`.
- Store weekly or presentation-ready report artifact packages under `docs/result_artifacts/weekly_report_result/`, separate from experiment-source artifact directories. Use subdirectories such as `report_<YYYY_Www>_<topic>/` and keep curated final figures, compact tables, manifests, and README files there.
- Local analysis or report-generation runs should write small final tables and figures only under `docs/result_artifacts/`. Do not create or update additional local ignored result mirrors unless the user explicitly asks.
- Large compute runs may write full outputs under established ignored result directories and also write or copy compact Git-trackable artifacts into `docs/result_artifacts/` for later review.
- Good artifact candidates include compact candidate tables, per-objective comparison tables, residue/interface summaries, mapping tables, assay summaries, small JSON/TSV manifests, and final analysis figures in `.png`, `.pdf`, or `.svg`.
- Every important decision-facing, report-facing, or publication-facing figure must retain a Git-tracked executable plotting script and the exact compact data table consumed by that script. Also retain provenance identifying upstream real-data files or authoritative records and key filtering, aggregation, uncertainty, and plotting parameters. The figure must be reproducible without conversational context or manually copied narrative values.
- Structure images must identify the source model or derived-structure manifest, chains and molecular roles, residue numbering/CDR scheme where relevant, aligned models and fitting selection where relevant, and the visualization script/session or reproducible commands.
- If an important figure depends on raw inputs too large for Git, keep those raw inputs in an established ignored result directory and track the extraction/aggregation script, exact compact plotted table, and a manifest containing upstream paths and, when practical, hashes or stable dataset/run identifiers. Do not reconstruct plot data manually from prose, captions, or remembered values.
- Analysis figures and result plots should use a clean scientific style: restrained color palettes, readable axis labels, consistent panel sizing, clear legends, no title/legend/label overlap, and 600 dpi export by default unless the user explicitly requests another resolution.
- Every figure containing error bars, intervals, or uncertainty bands must identify their exact meaning inside the figure through an axis label, legend, or concise figure note, for example `mean +/- SD`, `mean +/- SE`, `95% CI`, or `median with P10-P90`.
- Every nonstandard transformed or standardized axis must state its unit or definition inside the figure when omission could cause misinterpretation, for example `Δscore relative to Nb252`, `log10(x + 1)`, or a normalized multi-objective value.
- For multi-panel figures with a colorbar, reserve an explicit colorbar axis or sufficient right-side margin instead of letting the colorbar share or subtract space from the main panels; verify the saved image so the colorbar does not overlap or visually squeeze the rightmost plot.
- Tables intended for weekly reports or Word documents should use a three-line table style: top rule, header midrule, bottom rule, and no vertical grid. Keep columns and font sizes compact enough for a normal Word page width; split overly wide tables instead of making text unreadable.
- Do not track large generated outputs as artifacts: exhaustive candidate libraries, per-residue or per-pose full prediction tables, embeddings, trajectories, large structure ensembles, HDF5 files, parquet shards, large pickle files, checkpoints, or bulky intermediates should remain in ignored result/output directories.
- As a default size guideline, prefer artifacts under 5-10 MB. For larger but important results, create a compact artifact containing only the needed columns, rows, metrics, representative structures, or final figure.
- Result artifacts should contain observed or computed results, not narrative interpretation. Put interpretation and key conclusions in history or handoff documents.
- When a run produces small useful tables or final figures, copy or write them into `docs/result_artifacts/` before synchronization. Scripts must still not automatically commit or push them.
- History files may link to Git-tracked result artifact figures when the figure helps interpret a recorded result. Prefer links for most figures and embed images only for key final plots that directly support the entry's main conclusion.
- Keep embedded history figures rare and compact; as a default, use at most one embedded figure per history entry and keep exploratory, diagnostic, or repeated plots as linked artifacts instead.

## Weekly Audit Rules

- A systematic audit of the immediately preceding ISO week's code, structures, sequences, mappings, results, figures, and documentation is mandatory before that week's Word report is formally written or finalized.
- Before starting a weekly report task, verify that the corresponding audit report exists under `docs/weekly_audits/YYYY-Www.md` and is marked complete.
- If the user asks to write or finalize a weekly report before the required audit is complete, do not silently proceed. Remind the user to issue an explicit audit instruction, for example: `执行上一 ISO 周的项目周审计，系统检查代码、结构与序列映射、结果、图表和文档，并按风险等级生成审计报告。`
- Audit the prior week's changed code and tests, source and derived structures, sequence/numbering/CDR mappings, result-producing workflows, run summaries, tracked artifacts, important ignored outputs when available, figure source data and plotting scripts, history, handoff, and report claims.
- Check both engineering and scientific validity, including rule compliance, input/output contracts, chain and residue identity, missing-coordinate handling, experimental/predicted model separation, alignment and CDR semantics, design provenance, assay and unit alignment, seeds and model versions, threshold sources, leakage risk, aggregation logic, numerical consistency, figure reproducibility, and documentation consistency.
- Base every audit finding on concrete files or records. Cite relevant paths, rows, fields, commits, or line references. Mark unavailable inputs as `unverified`; never infer audit conclusions from conversational memory.
- Classify findings as `Critical`, `High`, `Medium`, `Low`, or `Observation`. Each finding must state the evidence, affected scope, likely impact, required action, and current status.
- Store each completed audit in `docs/weekly_audits/YYYY-Www.md`. Record the audited date range and commit range, inputs inspected, checks performed, findings by risk, unverified items, unresolved issues, and whether the weekly report may be finalized.
- Unresolved `Critical` or `High` findings block weekly-report finalization until they are corrected and rechecked. Lower-risk limitations must be disclosed when they affect report interpretation.
- Weekly audits are report-first by default. Do not modify experimental results, rewrite artifacts, launch compute jobs, or apply broad fixes during the audit unless the user explicitly approves remediation work.
- Automated inventory or validation scripts may support the audit, but they do not replace manual review of scientific semantics, structure/sequence provenance, experimental/predicted evidence separation, result provenance, and figure-to-data consistency.

## Documentation Maintenance Rules

- Do not create a separate README for every experiment or batch-script folder. Record experiment-specific inputs, outputs, commands, assumptions, and implementation details in the current weekly history file by default. Add a directory-level README only when the directory exposes a stable reusable interface, contains multiple long-lived workflows needing a shared entry point, or the user explicitly requests one.
- When a history entry introduces or materially changes a script, document that script's main inputs and outputs in the same entry. Include relevant source structures, sequences, mapping tables, assay datasets, model/checkpoint/database versions, and important CLI inputs, plus generated files/directories, formats, and whether they are ignored full results or Git-tracked summaries/artifacts. For multi-script pipelines, distinguish each stage's input/output contract.
- Weekly report drafting should follow the formatting style of the previous weekly report when available, including font family, font size, paragraph spacing, figure/table placement, and caption style.
- Weekly reports should describe the workflow as well as the results. Include relevant structure preparation, mapping, candidate design, filtering, prediction, experimental testing, scoring, and analysis processes when needed to understand the conclusion.
- Every experiment or conclusion included in a weekly report should be supported by a figure or table whenever possible. Figures must have clear captions, and tables must have clear notes/captions so the report remains readable without oral explanation.
- Before formally writing or editing a Word weekly report, first provide a proposed report outline and wait for the user's approval.
- Weekly report content does not need to follow strict chronological order. Organize by scientific story and user-specified scope; if the user assigns an experiment to a later report, defer it even if it appears in the current week's history.
- Whenever completed or active work is deliberately excluded from the current weekly report, record the deferred items and intended future reporting destination in `docs/codex-handoff.md` in the same task. Keep them visible until incorporated into a later report or explicitly removed from scope.
- `AGENTS.md` contains project-level rules and must not become a daily log.
- `docs/codex-handoff.md` is a replace-in-place current-state snapshot for new sessions, not an append-only history log. Keep its fixed sections.
- Maintain a `Last updated: YYYY-MM-DD HH:mm:ss` field near the top of `docs/codex-handoff.md` so staleness is visible.
- Use `docs/codex-handoff.md` only for current project status, active experiments, recent implementation state, decision-relevant cautions, and immediate next steps. Do not put current experiment status in `AGENTS.md`.
- Update `docs/codex-handoff.md` whenever a substantive implementation, workflow, experiment status, run-summary result, caution, or next step changes the current state a new agent needs to know. Every update must remove or rewrite superseded facts in the same task.
- After an experiment completes, remove its execution command and pending status from `Suggested Next Steps`; replace them with the observed result and actual unresolved follow-up, if any.
- When an implementation or scientific route changes, remove the superseded route from `Current Project Status` and `Suggested Next Steps`. Preserve historical details only in weekly history or run summaries.
- Before finishing any substantive implementation or result-sync task, audit the handoff for contradictions, especially an experiment being marked both completed and pending.
- Keep `Recent Changes` limited to the latest 5-8 milestones that affect current decisions. Move older completed work to `docs/history/` and retain only a link when needed.
- Keep `Suggested Next Steps` limited to at most 5 concrete, currently actionable items ordered by priority. Remove completed, abandoned, blocked-without-action, and speculative long-term tasks.
- Do not duplicate stable rules from `AGENTS.md`, detailed metrics already stored in run summaries, or chronological experiment logs in the handoff.
- Keep `docs/codex-handoff.md` below 120 lines and approximately 20 KB. If an update exceeds either limit, compact or remove stale content before completing the task.
- Before committing a handoff update, verify that its active route, latest completed result, cautions, and next steps agree with the current weekly history and available run summaries.
- When the user says run results or summaries have been synchronized to the workspace, the default task is to read the synchronized summary/result metadata, fold key facts into current weekly history, update `docs/codex-handoff.md` if project state changed, and analyze the run outcome.
- Result-sync analysis should state what succeeded or failed, the most important metrics or counts, caveats and missing statistics, and whether any follow-up source-file inspection or script change is needed. Do not only restate the raw summary.
- When a user-provided table, synchronized output excerpt, or local result inspection leads to a substantive interpretation that changes current understanding, record the key conclusion in `docs/history/YYYY-Www.md` and update `docs/codex-handoff.md` if it affects current state or next steps. This applies even when no new script ran and no run summary exists.
- `docs/history/YYYY-Www.md` is the detailed weekly log. Use Chinese. Each entry must include a timestamp formatted as `YYYY-MM-DD HH:mm:ss`.
- Keep entries inside each `docs/history/YYYY-Www.md` in chronological ascending order: oldest at the top, newest appended at the bottom.
- Every `docs/history/YYYY-Www.md` file must end with this unique append marker: `<!-- HISTORY_APPEND_MARKER: insert new entries immediately above this line -->`.
- New history entries must be inserted immediately above the append marker. Do not anchor history edits on repeated template text such as `Git sync status`.
- After editing a history file, verify all `## YYYY-MM-DD HH:mm:ss` headings remain in chronological ascending order.
- If history heading order verification fails, fix it before committing or sending the final response.
- `docs/weekly_history.md` is the long-term weekly summary index. Use Chinese. Do not duplicate every detailed history entry there.
- Before starting substantive work or before writing history:
  - Check the current system date and ISO week.
  - If still in the current week, insert detailed records immediately above the append marker in `docs/history/YYYY-Www.md`.
  - If a new week has started, summarize the previous week into `docs/weekly_history.md`, create the new weekly history file with the append marker at the end, and continue there.
