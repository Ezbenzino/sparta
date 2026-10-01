# PROMPT — 交给 Codex / AI Coding Agent 优化 SPARTA 投稿项目

> 用法：在项目根目录 `D:\sparta` 启动 Codex，把下面整段（从 `=== BEGIN ===` 到 `=== END ===`）粘给它。

=== BEGIN ===

You are working in the git repo at the current working directory (`D:\sparta`).
This is a **finished, already-submission-ready computational-biology manuscript** — NOT a greenfield project.
Your job is to **polish it to the standards of a computational biology journal**, NOT to redesign the science, NOT to change any result number, and NOT to add new experiments.

## 1. Project identity (read this first, do not change these facts)

- **Project**: SPARTA — "Two transport operators, one substrate: a graph model of T-cell migration and antibody penetration barriers in tumour tissue from spatial transcriptomics".
- **Target journal**: *Computational Biology and Chemistry* (Elsevier, subscription route, no APC; IF ≈ 2.6–2.9, Q3). The manuscript already exists at `docs/manuscript_cbc_draft.md` and a generated Word file `docs/manuscript_cbc_draft.docx`.
- **Author**: single author, Yize Li, Hangzhou Medical College, `lllyz630031258@gmail.com`. No funding. No co-authors. Do NOT add any author affiliation, funding line, or consortium.
- **Public code**: already on GitHub `https://github.com/Ezbenzino/sparta` (public, main branch, tag `v2.0.0`) and archived at Zenodo DOI `https://doi.org/10.5281/zenodo.23086431`. MIT license. `CITATION.cff` and `LICENSE` are filled.

## 2. Hard constraints — you MUST NOT violate these

1. **Do not change any reported number or conclusion.** The following are locked results, already in the paper and in the code output; if a re-run disagrees, FIX THE CODE to reproduce them, do NOT "correct" the paper:
   - 19 sections / 7 patients; cSCC 15 sections / 6 patients (GSE144239, two platform generations), metastatic melanoma 4 sections / 1 patient PtB (GSE250636). CSCC13 was rejected (UMI 289.5 < 300) and stays rejected.
   - Partial Spearman ρ positive in 18/19 sections; median ≈ +0.217.
   - Size-exclusion channel accounts for median 97.5% of B_mAb variance at β=3 (full β-sweep dependence already reported).
   - After removing the shared matrix input: coupling stays significant in 12/15 cSCC sections, median retention 81%; melanoma 3/4 reported as unresolved (one patient).
   - Single-section runtime < 0.2 s on CPU, no learned parameters.
   - 60–65% of edges exclude IgG (median mesh 4.36–4.68 nm).
2. **Do not invent data, datasets, or citations.** If a number in the paper cannot be reproduced from the existing `results/validation/*.json` and `data/ledger.csv`, flag it to me — do not fabricate a fix.
3. **Do not touch the experimental conclusion.** The paper deliberately reports its own negative/weak results (melanoma stratum unresolved, rearrangement control confounded with depth, antigen channel at detection floor). Keep them.
4. Python interpreter is `.venv\Scripts\python.exe` on Windows. `node` is available for `scripts/generate_docx.js`. Do not reinstall or upgrade packages in a way that breaks existing `results/` outputs.
5. After every change, re-run what you changed and show the command + output. Keep git commits small and conventional.

## 3. What I want you to actually do (in this priority order)

### P0 — Reproducibility & engineering hygiene (journal reviewers will check this)
- Add/verify `requirements.txt` (or `environment.yml`) that pins the exact versions used: numpy, scipy, networkx, pandas, matplotlib, scikit-learn, statsmodels, PIL/Pillow. Derive versions from the running `.venv`, do not guess.
- Add a top-level `README.md` (English) with: what SPARTA does in 3 bullets, the locked result summary, one-command reproduction steps (`python scripts/run_all.py` if such an orchestrator exists — if not, write a thin `run_all.py` that calls the existing `run_01..run_16` scripts in order and is idempotent), expected runtime, and how to cite (the Zenodo DOI).
- Audit every script under `scripts/` for: fixed random seeds (all permutation tests must use a fixed seed, e.g. `numpy.random.default_rng(20260828)`), hard-coded absolute paths (`d:\sparta` → use `pathlib` relative to repo root), and `plt.savefig` DPI/`bbox_inches`. Make figures regenerate at 300 DPI.
- Verify the `.gitignore` actually excludes `node_modules/`, `.venv/`, large `docs/*.docx` (they are ~8 MB each), raw data caches, and `__pycache__`. The repo currently ships 119 files / ~1.1 MB — keep it that way.

### P1 — Correctness audit of the two core operators (this is the science reviewers will attack)
- Read `scripts/` implementing (a) the source–sink minimum-cut cellular barrier and (b) the screened-Poisson diffusion–absorption antibody barrier with size exclusion.
- Confirm against the equations stated in `docs/manuscript_cbc_draft.md` §2.4:
  - min-cut: source = endothelial-rich nodes, sink = malignant/immune-poor nodes, capacity decreasing in matrix resistance, barrier = 1/max-flow, preflow-push on the largest connected component.
  - screened Poisson: vessel nodes = unit concentration, edge conductance implements mesh ξ = ξ0·exp(−β·x), edge is size-excluded when ξ < r (r=5.5 nm) with a conductance floor, antigen as absorption sink, steady state solved as a linear Helmholtz system.
- For each, add a **unit test** under `tests/` (pytest) that checks a tiny hand-computable graph: e.g. two-node graph, three-node chain, fully-excluded vs fully-open mesh. The test must pass and must document the expected value. This is what a reviewer means by "implementation correctness".
- If you find a real discrepancy between code and paper, STOP and report it to me with file:line — do not silently patch the paper.

### P2 — Figure & manuscript polish (visual/formatting, not science)
- Fix the known cosmetic issue: in `results/figures/fig3_size_scan.png` panel (b), the per-section numeric labels above the bars overlap. Adjust rotation/offset/font size so they are legible; keep the same data. Regenerate via the existing figure script.
- Unify figure style: font family, font size, axis/spine color, and the blue (B_cell) / orange (B_mAb) palette across fig1–fig5 so they look like one paper.
- Confirm `results/figures/graphical_abstract_portrait.png` (1452×3432, vertical) is the one embedded; the old horizontal version must not be referenced anywhere.
- Run `node scripts/generate_docx.js` after any figure change and confirm `docs/manuscript_cbc_draft.docx` regenerates without errors and still contains: 14 references, author "Yize Li", the GitHub URL and Zenodo DOI, "The author received no specific funding", and the portrait graphical abstract.

## 4. Workflow rules

1. **First turn**: run a read-only audit — list the repo tree, read `docs/manuscript_cbc_draft.md`, `README*`, `requirements*`, and the core operator scripts. Then output a short plan (P0/P1/P2 checklist) and wait for my go-ahead before editing.
2. Make changes in small commits; after each P0/P1 item, re-run the relevant script/test and paste the tail of the output as evidence.
3. If something in the existing pipeline is already correct, say so and move on — do not refactor for the sake of it.
4. At the end, give me: (a) a list of files changed, (b) test/run output proving reproduction, (c) any discrepancy you found but did not fix, and (d) the exact git commands I should run to push.

Do NOT touch: `docs/manuscript_cbc_draft_zh.md`, the Zenodo/GitHub releases, `data/` raw ledger values, or the manuscript's Results/Discussion prose. If you think prose must change, write the proposed edit as a suggestion and let me decide.

=== END ===

## 给你（用户）的补充说明（不粘进 Codex）

- 这个提示词把"能改什么/不能改什么"锁死了，防止 Codex 自作主张改实验数字或加实验——那是你最该防的。
- P0（复现性/README/种子）是 CBC 审稿人最常挑、也最容易补的，优先让它做。
- P1（两个算子的单元测试）是真正能提升论文可信度的部分：加几个手算小图的 pytest，审稿人会觉得实现严谨。
- 如果 Codex 在 P1 报告"代码和论文公式对不上"，把它的 `file:line` 发给我，我帮你判断是改代码还是改论文——这一步不要让它自己拍板。
