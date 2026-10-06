Yize Li
Hangzhou Medical College, Hangzhou, Zhejiang, China
lllyz630031258@gmail.com

{{date}}

Prof. Dongqing Wei and Prof. Peng Li
Co-Editors-in-Chief, *Interdisciplinary Sciences: Computational Life Sciences*

Dear Professor Wei and Professor Li,

Please consider the enclosed manuscript, "SPARTA: graph structural operators and null models for characterizing spatial tissue architecture in spatial omics", for publication as an original research article in *Interdisciplinary Sciences: Computational Life Sciences*.

Spatial omics can map cell states and matrix programmes, but composition alone does not describe how those programmes are connected. SPARTA formulates two graph structural operators on one spatial graph: a source–sink minimum cut for separation between vessel-proximal immune positions and tumour-nest positions, and a screened diffusion–absorption field for a molecular-sized in silico probe. The manuscript gives equal weight to the statistics needed to interpret these operators: a calibrated surrogate test, patient-level inference and structural null models that measure how much coupling the model produces by itself.

- In {{syn_n_tissues}} simulated tissues with planted barriers and an agent-based ground truth, the minimum cut tracked cell arrival (Spearman ρ = {{syn_sparta_bcell_rho}}) and was the only summary that separated closed capsules from capsules with a 5% gap.
- In multiplexed protein images of {{cx_n_cores}} colorectal-cancer cores from {{cx_n_patients}} patients, the cut computed from the stromal scaffold alone, with all immune cells withheld, was associated with measured CD8^+^ T-cell depletion from tumour cores (patient-level ρ = {{cx_rho_pat}}) under a protocol written before the analysis; this is structural validation against held-out cell positions rather than a prognostic test.
- Graph-spectral surrogates held the false-positive rate near 5% on all {{cal30_n_graphs}} real section graphs, whereas a spot-level test reached {{cal_pooled_naive}}.
- In 19 transcriptomics sections from eight patients the two fields co-varied in every patient (pooled ρ = {{pl_mean}}), and the association replicated in an independent cohort of four melanoma patients ({{rep_pl_mean}}); the explicitly labelled primary-plus-replication analysis covered 12 patients.
- A 2026 public extension cohort of 50 sections from six GEO series, spanning Visium and Slide-seqV2, reproduced the direction in 48 sections (pooled ρ = {{ext2026_pl_mean}}); patient relationships were unavailable for 11 cSCC sections and were handled as conservative section-level units. Structural nulls and shared-ECM ablation again reduced but did not eliminate the association.
- Simple stromal-density, distance and neighbourhood summaries reproduced only part of either operator; in the images, a geometric summary added after the analysis showed the strongest single association, while the cut retained information beyond composition and geometry.
- In-model perturbation maps identify node sets that affect the computational outputs: the minimum cut is relatively node-sensitive, whereas the diffusion field is less focal.

The work combines graph algorithms, spatial statistics and structural validation, which I believe suits the journal's interdisciplinary scope.

The manuscript is prepared for double-blind review: identifying information is on the separate title page, and the code, derived node tables and result files are provided as an anonymised archive (Online Resource 3). All data are public (GEO GSE144239, GSE250636, GSE200278, GSE289745, GSE300445, GSE316760, GSE320041 and GSE321832; 10x Genomics; Thrane et al. 2018; Schürch et al. 2020), and every analysis runs on a laptop CPU.

This manuscript has not been published and is not under consideration elsewhere; I am its sole author, approve the submitted version, have no competing interests and received no funding for this work.

Thank you for considering this submission.

Sincerely,

Yize Li
