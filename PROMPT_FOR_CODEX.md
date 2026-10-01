# PROMPT — 交给 Codex / AI Coding Agent 完成 SPARTA 投稿前最后工程收尾

> 用法：在项目根目录 `D:\sparta` 启动 Codex，把下面整段（从 `=== BEGIN ===` 到 `=== END ===`）粘给它。

=== BEGIN ===

You are working in the git repo at the current working directory (`D:\sparta`).
This is a **finished, submission-ready computational-biology manuscript** targeting *Computational Biology and Chemistry* (Elsevier, hybrid, subscription route = no APC). Do NOT redesign the science, do NOT change any result number, do NOT add new experiments. The paper prose is frozen in `docs/manuscript_cbc_draft.md`; your job is engineering reproducibility and Word export.

## 1. Locked facts — do not change

- Single author: **Yize Li**, Hangzhou Medical College, `lllyz630031258@gmail.com`. No funding, no co-authors.
- GitHub: `https://github.com/Ezbenzino/sparta`; Zenodo DOI: `https://doi.org/10.5281/zenodo.23086431`; MIT license.
- Cohort: 19 sections / 7 patients (15 cSCC / 6 patients from GSE144239 across two platform generations; 4 melanoma deposits / 1 patient PtB from GSE250636). CSCC13 rejected and stays rejected.
- Locked numbers: partial ρ positive in 18/19 (median +0.217); size-exclusion channel median 97.5% of B_mAb variance at β=3; shared-input removal leaves coupling significant in 12/15 cSCC (median 81% retention); single section <0.2 s CPU.
- **New results already in the manuscript** (do not re-run, do not "fix"):
  - Null crosslink permutation: median real 97.6% vs null 97.6% (permutation confirms the 97.5% is a rank-normalisation property, not a data signal — the paper now reports this honestly).
  - In silico stromal intervention: 30% reduction in ecm/crosslink/caf lowers B_cell by median 60.8% and B_mAb by 49.5%, in 19/19 sections; 50% reduction → 80.2% / 71.3%.
  - Naive baseline: graph B_cell field partial ρ vs B_mAb median +0.217 vs naive 0.5·(ECM+CAF) at +0.166; graph stronger in 14/19.
- Python: `.venv\Scripts\python.exe` (Windows). Node is available. Do not reinstall/upgrade packages.

## 2. Current state — what is already done (verify, do not redo)

- `requirements.txt`, `environment.yml`, `README.md`, `scripts/run_all.py`, `tests/test_barrier.py` exist from a prior session. **First task**: run `pytest tests/ -v` and paste the output. If tests fail, fix them. If `run_all.py` does not call `run_22_null_crosslink.py`, `run_23_stromal_intervention.py`, `run_24_naive_baseline.py`, add them as optional tail stages (they are fast, read-only on the existing h5ad/graph artefacts).
- Fig 3b bar labels were already de-overlapped (per-bar numbers removed; median dashed line kept). Do not re-add per-bar labels.
- Manuscript markdown `docs/manuscript_cbc_draft.md` already contains the new paragraphs. Do not edit the prose.

## 3. Remaining work — do these in order

### P0 (highest priority): sync `scripts/generate_docx.js` with the current manuscript
The Word export script `scripts/generate_docx.js` is hard-coded with an older version of the prose. It now lags `docs/manuscript_cbc_draft.md` in these specific places — update the hard-coded strings to match:

1. **Abstract**: after "...no learned parameters," add the stromal-intervention sentence: "In silico reduction of matrix density and crosslinking by 30% lowers both barriers by a median 61% and 50% respectively, in 19 of 19 sections across both tumour types."
2. **§1 Introduction**: before "Our contributions are:", add the positioning paragraph about methodological contribution (not clinical biomarker).
3. **§3.2**: after the β-sweep paragraph, add the permutation control paragraph (null crosslink: median real 97.6% vs permuted 97.6%, 95% null interval 96.8–98.4%, 2/6 sections p<0.05).
4. **§3.5**: after the spatial-rearrangement paragraph, add the "In silico stromal co-targeting" paragraph (30% reduction → B_cell −60.8%, B_mAb −49.5%; 50% → −80.2%/−71.3%; both barriers drop in 19/19).
5. **§3.6**: before "## 4. Discussion", add the "Why a graph at all?" paragraph (naive 0.5·(ECM+CAF) baseline: median partial ρ +0.166 vs graph +0.217; graph stronger in 14/19).
6. **§4 Limitations**: add the single-author sentence ("This is a single-author work: the analysis was not independently checked by a second analyst; the code is released for community verification.").
7. **§5 Conclusions**: before "For matrix-rich tumours the model predicts", add the quantitative stromal-intervention sentence.
8. **References**: verify refs [3] (Ospina/spatialGE), [4] (Singhal/BANKSY), [6] (Ren/SpaceFlow), [13] (Riaz), [14] (Jain) match the markdown.
9. After editing, run `node scripts/generate_docx.js` and confirm `docs/manuscript_cbc_draft.docx` regenerates with no error. Open it mentally (or via `unzip -p word/document.xml | grep -c`) to confirm the new strings are present.

### P1: Make the repo reviewer-proof
- Run `pytest tests/ -v`. If `tests/test_barrier.py` only covers compute_b_cell/compute_b_mab on synthetic graphs, add a test that the size-exclusion floor is non-negative and that a two-node fully-excluded graph yields B_mAb = conductance-floor (hand-computable).
- Add a short `tests/test_new_experiments.py` that loads `results/validation/null_crosslink_check.json`, `stromal_intervention.json`, `naive_baseline.json` and asserts: (a) all 19 slides appear in stromal_intervention, (b) median B_cell drop at 30% is between 50% and 70%, (c) naive_baseline n_slides == 19. These are smoke tests, not science.
- Audit `scripts/run_22_null_crosslink.py`, `run_23_stromal_intervention.py`, `run_24_naive_baseline.py`: they must use fixed seeds (`numpy.random.default_rng(20261001)`), no hard-coded `d:\sparta` paths, and write only to `results/validation/`.
- Update `README.md` to mention the three new scripts and the new result files.

### P2: Final verification
- Run `python scripts/run_all.py` (or at least `run_22`, `run_23`, `run_24`) and confirm they reproduce the JSON files without errors.
- Confirm `.gitignore` excludes `docs/manuscript_cbc_draft.docx` (it is ~3.7 MB and must not be committed).
- `git add -A && git commit -m "Sync docx generator with new experiments; add smoke tests" && git push origin main`.

## 4. Hard do-nots

- Do NOT edit `docs/manuscript_cbc_draft.md` prose (it is frozen).
- Do NOT change any number in the paper, the JSON outputs, or the figure data.
- Do NOT add new experiments, new citations beyond what is already in the markdown, or new figures.
- Do NOT touch `data/ledger.csv`, the Zenodo DOI, the GitHub releases, or the LICENSE/CITATION.cff.
- If you find a real bug where code output disagrees with the paper, STOP and report it with file:line — do not patch the paper.

## 5. Deliverable

At the end, paste: (a) list of files changed, (b) `pytest` summary output, (c) the `node scripts/generate_docx.js` output confirming the docx regenerated, (d) any discrepancy you found but did not fix, (e) the final git log one-liner.

=== END ===

## 给你（用户）的补充说明（不粘进 Codex）

- Codex 之前已经做了一部分 P0-C（requirements/README/tests 框架），这次主要是**同步 generate_docx.js**——这是机械但繁琐的活，让它干。
- 三个新实验脚本（run_22/23/24）已经在本地跑过、结果 json 已生成、论文 markdown 已更新。Codex 不需要重跑科学实验，只需要：(1) 确认 JSON 能复现，(2) 把论文新段落同步进 Word 导出脚本，(3) 加冒烟测试。
- 投稿前你只需要做：在 Editorial Manager 上传 docx + graphical abstract + cover letter + 3 个建议审稿人。
