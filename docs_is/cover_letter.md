Yize Li
Hangzhou Medical College, Hangzhou, Zhejiang, China
lllyz630031258@gmail.com

{{date}}

Prof. Dongqing Wei and Prof. Peng Li
Co-Editors-in-Chief, *Interdisciplinary Sciences: Computational Life Sciences*

Dear Professor Wei and Professor Li,

Please consider the enclosed manuscript, "SPARTA: graph transport operators and structural null models for immune-cell and IgG-transport barriers in spatial omics", for publication as an original research article in *Interdisciplinary Sciences: Computational Life Sciences*.

Many questions about tumours are questions about transport: can effector T cells that leave a vessel reach tumour nests, and can an antibody reach the tumour core before it is excluded or bound? SPARTA poses both problems on one spatial graph, as a source–sink minimum cut for migrating cells and a screened diffusion–absorption equation for an IgG-sized molecule, and gives equal weight to the statistics needed to interpret them: a calibrated surrogate test, patient-level inference and structural null models that measure how much coupling the model produces by itself.

- In {{syn_n_tissues}} simulated tissues with planted barriers and an agent-based ground truth, the minimum cut tracked cell arrival (Spearman ρ = {{syn_sparta_bcell_rho}}) and was the only summary that separated closed capsules from capsules with a 5% gap.
- In multiplexed protein images of {{cx_n_cores}} colorectal-cancer cores from {{cx_n_patients}} patients, the cut computed from the stromal scaffold alone, with all immune cells withheld, predicted measured CD8^+^ T-cell depletion from tumour cores (patient-level ρ = {{cx_rho_pat}}) under a protocol written before the analysis; it did so no better than peritumoural matrix density, which the manuscript reports.
- Graph-spectral surrogates held the false-positive rate near 5% on all {{cal30_n_graphs}} real section graphs, whereas a spot-level test reached {{cal_pooled_naive}}.
- In 19 transcriptomics sections from seven patients the two fields co-varied in every patient (pooled ρ = {{pl_mean}}), and the association replicated in an independent cohort of four melanoma patients ({{rep_pl_mean}}); structural nulls attribute a median {{dec_share_constr_pct}}% and {{rep_share_constr_pct}}% of the coupling, respectively, to shared inputs and geometry.
- Simple stromal-density, distance and neighbourhood summaries reproduced only part of either operator; in the images, a geometric summary added after the analysis predicted depletion best, but the cut kept information beyond composition and geometry.
- In-model intervention maps show a focal cellular barrier (top 1% of spots: a median {{iv_top1_cell_median}}% of a full breach) and a distributed antibody barrier ({{iv_top1_mab_median}}%).

The work combines graph algorithms, transport physics and spatial statistics to address a question from tumour biology, which I believe suits the journal's interdisciplinary scope.

The manuscript is prepared for double-blind review: identifying information is on the separate title page, and the code, derived node tables and result files are provided as an anonymised archive (Online Resource 3). All data are public (GEO GSE144239 and GSE250636; 10x Genomics; Thrane et al. 2018; Schürch et al. 2020), and every analysis runs on a laptop CPU.

This manuscript has not been published and is not under consideration elsewhere; I am its sole author, approve the submitted version, have no competing interests and received no funding for this work.

Thank you for considering this submission.

Sincerely,

Yize Li
