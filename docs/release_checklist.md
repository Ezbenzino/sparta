# CMPB analysis archive and software release checklist

Use this checklist to make the submitted manuscript traceable to one immutable software-and-results snapshot. A DOI must identify the actual deposited artifact; never copy a draft DOI or mint one before the archive exists.

## Freeze contents

- [ ] Confirm the exact source commit/tag and record the working-tree diff used for that release.
- [ ] Include package source, `pyproject.toml`, dependency specification, configs, analysis scripts, plotting scripts, tests, README/run instructions, and citation metadata.
- [ ] Include the per-section JSON outputs and figures needed to reproduce every manuscript number and panel. Keep raw GEO data and private credentials out of the archive; document public accessions and download steps.
- [ ] Include the validation outputs `s2_matched_selection.json`, `spatial_null_check.json`, `metric_connectivity_sensitivity.json`, `benchmark_tools.json`, and the figure inputs used by the CMPB draft.
- [ ] Include a manifest with file checksums, software/runtime versions, config paths, random seeds, and generation commands.
- [ ] Run the documented workflow from a clean environment as far as public inputs and available compute permit; list any expensive analyses not rerun and identify the archived outputs they depend on.
- [ ] Verify the archive opens, required files are present, and no raw/private data or local absolute paths are included unintentionally.

## Publish and cite

- [ ] Review all staged files before creating or pushing a public release. The working tree currently contains many pre-existing uncommitted project changes; do not bulk-stage them without reviewing the diff.
- [ ] Create a versioned GitHub release that matches the archived source and results.
- [ ] Deposit the exact release archive in Zenodo (or an equivalent public repository), verify the landing page and files, then record the minted DOI in `CITATION.cff`, the manuscript data-availability statement, and README.
- [ ] Validate `CITATION.cff` and ensure package version, tag, release date, repository URL and DOI refer to the same snapshot.
- [ ] Rebuild the manuscript DOCX, figures and Highlights from the frozen source. Record checksums so the uploaded files can be matched to the archive.

## Final CMPB submission checks

- [ ] Verify live CMPB Guide for Authors instructions and Editorial Manager file fields on the submission date.
- [ ] Confirm author names, affiliations, corresponding-author details, CRediT roles, ethics statement, funding, competing interests, data statement and AI declaration.
- [ ] Ensure all results in the manuscript agree with the archived JSON; check all figure legends, citations and references.
- [ ] Upload the manuscript, separate Highlights and figure/supplement files in the requested editable formats and resolutions.
- [ ] Confirm the manuscript is not under consideration elsewhere and that all authors approved the submitted version.
