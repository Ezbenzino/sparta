# -*- coding: utf-8 -*-
import io
p = r"D:\sparta\scripts\generate_docx.js"
s = io.open(p, encoding="utf-8").read()

# 1) refs 数组整块替换
old_refs_start = "const refs = ["
old_refs_end = "];"
i0 = s.index(old_refs_start)
i1 = s.index(old_refs_end, i0) + len(old_refs_end)
new_refs = '''const refs = [
  "[1] P.L. Staahl, F. Salmen, S. Vickovic, et al., Visualization and analysis of gene expression in tissue sections by spatial transcriptomics, Science 353 (2016) 78-82.",
  "[2] A.L. Ji, D.M. Rubin, K. Tharakan, et al., Multimodal analysis of the cellular and molecular landscape of cutaneous squamous cell carcinoma, Cell 182 (2020) 497-514.e20.",
  "[3] O.E. Ospina, R. Manjarres-Betancur, G. Gonzalez-Calderon, et al., B.L. Fridley, spatialGE is a user-friendly web application that facilitates spatial transcriptomics data analysis, Cancer Research 85 (2025) 848.",
  "[4] V. Singhal, N. Chou, J. Lee, Y. Yue, J. Liu, W.K. Chock, L. Lin, Y.-C. Chang, K.H. Chen, S. Prabhakar, BANKSY unifies cell typing and tissue domain segmentation for scalable spatial omics data analysis, Nature Genetics 56 (2024) 431-441.",
  "[5] G. Palla, H. Spitzer, M. Klein, et al., Squidpy: a scalable framework for spatial omics analysis, Nature Methods 19 (2022) 171-178.",
  "[6] H. Ren, B.L. Walker, Z. Cang, Q. Nie, Identifying multicellular spatiotemporal organization of cells with SpaceFlow, Nature Communications 13 (2022) 4076.",
  "[7] K.A. Doiron, et al., Tissue cellular-neighbourhood analysis of spatial transcriptomics, Nature Methods 20 (2023).",
  "[8] A. Liberzon, C. Birger, H. Thorvaldsdottir, et al., The Molecular Signatures Database hallmark gene set collection, Cell Systems 1 (2015) 17-25.",
  "[9] F.A. Wolf, P. Angerer, F.J. Theis, SCANPY: large-scale single-cell gene expression data analysis, Genome Biology 19 (2018) 15.",
  "[10] A.A. Hagberg, D.A. Schult, P.J. Swart, Exploring network structure, dynamics, and function using NetworkX, in: Proc. 7th Python in Science Conference, 2008, pp. 11-15.",
  "[11] P.G. Doyle, J.L. Snell, Random Walks and Electric Networks, Mathematical Association of America, Washington, 1984.",
  "[12] W. Hugo, J.M. Zaretsky, L. Sun, et al., Genomic and transcriptomic features of response to anti-PD-1 therapy in metastatic melanoma, Cell 165 (2016) 35-44.",
  "[13] N. Riaz, J.J. Havel, V. Makarov, A. Desrichard, W.J. Urba, J.S. Sims, F.S. Hodi, S. Martin-Algarra, R. Mandal, W.H. Sharfman, T.A. Chan, Tumor and microenvironment evolution during immunotherapy with nivolumab, Cell 171 (2017) 934-949.e16.",
  "[14] R.K. Jain, Delivery of molecular and cellular medicine to tumors, Nature Reviews Drug Discovery 4 (2005) 619-632.",
];'''
s = s[:i0] + new_refs + s[i1:]

# 2) CRediT 作者
s = s.replace('children.push(p("[Author name]: Conceptualisation, Methodology, Software, Validation, Formal analysis, Investigation, Data curation, Writing \\u2014 original draft, Visualisation."));',
              'children.push(p("Yize Li: Conceptualisation, Methodology, Software, Validation, Formal analysis, Investigation, Data curation, Writing \\u2014 original draft, Visualisation."));')

# 3) Data availability 段
s = s.replace("The SPARTA implementation, all analysis scripts, configuration files and per-section artefacts are available at [GitHub URL \\u2014 to be filled after first commit] and archived at [Zenodo DOI \\u2014 to be filled].",
              "The SPARTA implementation, all analysis scripts, configuration files and per-section artefacts are available at https://github.com/Ezbenzino/sparta and archived at https://doi.org/10.5281/zenodo.23086431.")

# 4) Acknowledgements
s = s.replace('children.push(p("[Funding/acknowledgements to be completed.]"));',
              'children.push(p("The author received no specific funding for this work."));')

# 5) graphical abstract 改竖版（portrait, width 3.0in, ratio 1452/3432=0.423）
s = s.replace('children.push(...figImage("graphical_abstract", 5.5, 2.19, null));',
              'children.push(...figImage("graphical_abstract_portrait", 3.0, 0.423, null));')

io.open(p, "w", encoding="utf-8", newline="\n").write(s)
print("generate_docx.js updated")
print("leftover placeholders:", s.count("[Author")+s.count("[GitHub URL")+s.count("[Zenodo")+s.count("[Funding")+s.count("[Background"))
