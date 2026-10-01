/**
 * generate_docx.js — CBC 投稿完全体初稿 Word 文档生成
 *
 * 从 manuscript_cbc_draft.md 的正文文本 + results/figures/*.png + JSON 产物
 * 生成带插图的 .docx 文件。
 */
const fs = require("fs");
const path = require("path");
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell, ImageRun,
  Header, Footer, AlignmentType, LevelFormat, HeadingLevel, BorderStyle,
  WidthType, ShadingType, VerticalAlign, PageNumber, PageBreak,
} = require("docx");

const ROOT = "d:\\sparta";
const FIG = path.join(ROOT, "results", "figures");
const VAL = path.join(ROOT, "results", "validation");

// ── 样式常量 ──
const FONT = { ascii: "Times New Roman", hAnsi: "Times New Roman", eastAsia: "Microsoft YaHei" };
const SZ_BODY = 24;       // 12pt
const SZ_H1 = 28;         // 14pt
const SZ_H2 = 26;         // 13pt
const SZ_TITLE = 32;      // 16pt
const SZ_SMALL = 20;      // 10pt
const SZ_CAP = 22;        // 11pt
const COLOR_HEAD = "1a1a1a";
const COLOR_MUTED = "666666";
const COLOR_CELL = "2a78d6";
const COLOR_MAB = "eb6834";

const border = { style: BorderStyle.SINGLE, size: 1, color: "BBBBBB" };
const borders = { top: border, bottom: border, left: border, right: border };
const noBorders = {
  top: { style: BorderStyle.NONE, size: 0 },
  bottom: { style: BorderStyle.NONE, size: 0 },
  left: { style: BorderStyle.NONE, size: 0 },
  right: { style: BorderStyle.NONE, size: 0 },
};

// ── 辅助 ──
function p(text, opts = {}) {
  return new Paragraph({
    alignment: opts.align || AlignmentType.JUSTIFIED,
    spacing: { after: opts.after ?? 120, line: opts.line ?? 276 },
    children: [new TextRun({ text, font: FONT, size: opts.size || SZ_BODY, bold: opts.bold, italics: opts.italics, color: opts.color })],
  });
}
function pRuns(runs, opts = {}) {
  return new Paragraph({
    alignment: opts.align || AlignmentType.JUSTIFIED,
    spacing: { after: opts.after ?? 120, line: opts.line ?? 276 },
    children: runs,
  });
}
function h1(text) {
  return new Paragraph({
    heading: HeadingLevel.HEADING_1,
    spacing: { before: 360, after: 180 },
    children: [new TextRun({ text, font: FONT, size: SZ_H1, bold: true, color: COLOR_HEAD })],
  });
}
function h2(text) {
  return new Paragraph({
    heading: HeadingLevel.HEADING_2,
    spacing: { before: 240, after: 120 },
    children: [new TextRun({ text, font: FONT, size: SZ_H2, bold: true, color: COLOR_HEAD })],
  });
}
function emptyP() {
  return new Paragraph({ children: [new TextRun("")] });
}
function bullet(text) {
  return new Paragraph({
    numbering: { reference: "bullets", level: 0 },
    spacing: { after: 60, line: 276 },
    children: [new TextRun({ text, font: FONT, size: SZ_BODY })],
  });
}
function numItem(text) {
  return new Paragraph({
    numbering: { reference: "numbers", level: 0 },
    spacing: { after: 60, line: 276 },
    children: [new TextRun({ text, font: FONT, size: SZ_BODY })],
  });
}
function figImage(filename, widthIn, ratio, caption) {
  const h = Math.round(widthIn / ratio);
  const wPt = Math.round(widthIn * 72);
  const hPt = Math.round(h * 72);
  const imgPath = path.join(FIG, filename + ".png");
  const data = fs.readFileSync(imgPath);
  const children = [
    new Paragraph({
      alignment: AlignmentType.CENTER,
      spacing: { before: 200, after: 80 },
      children: [new ImageRun({
        type: "png",
        data,
        transformation: { width: wPt, height: hPt },
        altText: { title: filename, description: caption || filename, name: filename },
      })],
    }),
  ];
  if (caption) {
    children.push(new Paragraph({
      alignment: AlignmentType.CENTER,
      spacing: { after: 200 },
      children: [new TextRun({ text: caption, font: FONT, size: SZ_CAP, color: COLOR_MUTED })],
    }));
  }
  return children;
}
function cell(text, opts = {}) {
  return new TableCell({
    borders,
    width: { size: opts.w || 2000, type: WidthType.DXA },
    shading: opts.fill ? { fill: opts.fill, type: ShadingType.CLEAR } : undefined,
    margins: { top: 40, bottom: 40, left: 80, right: 80 },
    verticalAlign: VerticalAlign.CENTER,
    children: [new Paragraph({
      alignment: opts.align || AlignmentType.LEFT,
      spacing: { after: 0 },
      children: [new TextRun({ text, font: FONT, size: opts.size || SZ_SMALL, bold: opts.bold })],
    })],
  });
}
function headerCell(text, w) {
  return cell(text, { w, fill: "D5E8F0", bold: true, align: AlignmentType.CENTER });
}

// ── 载入 JSON 产物 ──
function loadJSON(p) {
  return JSON.parse(fs.readFileSync(p, "utf-8"));
}

// ════════════════════════════════════════════════════════════════
// 构建 Table 1 (队列概况)
function buildTable1() {
  const ledger = fs.readFileSync(path.join(ROOT, "data", "ledger.csv"), "utf-8");
  const lines = ledger.trim().split("\n");
  const header = lines[0].split(",");
  const rows = lines.slice(1).map(l => {
    const vals = {};
    header.forEach((h, i) => { vals[h] = l.split(",")[i]; });
    return vals;
  }).filter(r => r.status === "ingested" || r.status === "rejected");

  const cols = [
    { key: "slide_id", w: 800, label: "Section" },
    { key: "patient", w: 900, label: "Patient" },
    { key: "cancer_type", w: 700, label: "Type" },
    { key: "platform", w: 700, label: "Platform" },
    { key: "n_spots", w: 600, label: "Spots" },
    { key: "median_umi", w: 700, label: "Median UMI" },
    { key: "status", w: 800, label: "Status" },
  ];
  const totalW = cols.reduce((s, c) => s + c.w, 0);

  const headerRow = new TableRow({
    cantSplit: true,
    children: cols.map(c => headerCell(c.label, c.w)),
  });

  const dataRows = rows.map(r => new TableRow({
    cantSplit: true,
    children: cols.map(c => {
      const v = r[c.key] || "";
      const isRejected = r.status === "rejected";
      return cell(v, { w: c.w, size: 18, align: AlignmentType.CENTER, bold: isRejected });
    }),
  }));

  return new Table({
    width: { size: 100, type: WidthType.PERCENTAGE },
    columnWidths: cols.map(c => c.w),
    rows: [headerRow, ...dataRows],
  });
}

// ════════════════════════════════════════════════════════════════
// 构建 Table 2 (R3b 分层耦合表)
function buildTable2() {
  const chk = loadJSON(path.join(VAL, "shared_ecm_check.json"));
  const slides = chk.per_slide;
  // 计算分层统计
  const strata = [
    { name: "cSCC, Visium", filter: s => s.startsWith("CSCC") && parseInt(s.slice(4)) <= 4 },
    { name: "cSCC, 1st-gen ST", filter: s => s.startsWith("CSCC") && parseInt(s.slice(4)) >= 5 },
    { name: "cSCC combined", filter: s => s.startsWith("CSCC") },
    { name: "Melanoma (one patient)", filter: s => s.startsWith("MEL") },
  ];
  const pmap = {
    CSCC01:"P4",CSCC02:"P4",CSCC03:"P6",CSCC04:"P6",
    CSCC05:"P2",CSCC06:"P2",CSCC07:"P2",CSCC08:"P5",CSCC09:"P5",CSCC10:"P5",
    CSCC11:"P9",CSCC12:"P9",CSCC14:"P10",CSCC15:"P10",CSCC16:"P10",
    MEL01:"PtB",MEL02:"PtB",MEL03:"PtB",MEL04:"PtB",
  };
  const results = strata.map(st => {
    const members = Object.keys(slides).filter(st.filter);
    const n = members.length;
    const patients = new Set(members.map(s => pmap[s])).size;
    const rhos = members.map(s => slides[s]["切断共享ECM"].rho_partial);
    const median = rhos.sort((a, b) => a - b)[Math.floor(rhos.length / 2)];
    const sigCount = members.filter(s => slides[s]["切断共享ECM"].p_partial < 0.05 && slides[s]["切断共享ECM"].rho_partial > 0).length;
    return { name: st.name, sections: n, patients, rho: median.toFixed(3), sig: sigCount + " / " + n };
  });

  const cols = [
    { w: 2400, label: "Stratum" },
    { w: 1400, label: "Sections / patients" },
    { w: 1400, label: "ρ after removal" },
    { w: 1400, label: "Significant" },
    { w: 2760, label: "Patients majority-significant" },
  ];

  const rows = [
    new TableRow({ cantSplit: true, children: cols.map(c => headerCell(c.label, c.w)) }),
    ...results.map(r => new TableRow({
      cantSplit: true,
      children: [
        cell(r.name, { w: 2400 }),
        cell(r.sections + " / " + r.patients, { w: 1400, align: AlignmentType.CENTER }),
        cell("+" + r.rho, { w: 1400, align: AlignmentType.CENTER }),
        cell(r.sig, { w: 1400, align: AlignmentType.CENTER }),
        cell(r.patients > 1 ? "see text" : "—", { w: 2760, align: AlignmentType.CENTER }),
      ],
    })),
  ];

  return new Table({
    width: { size: 100, type: WidthType.PERCENTAGE },
    columnWidths: cols.map(c => c.w),
    rows,
  });
}

// ════════════════════════════════════════════════════════════════
// 构建 Table 3 (运行时对比)
function buildTable3() {
  const rt = loadJSON(path.join(VAL, "runtime_benchmark.json"));
  const per = rt.per_slide;
  const slides = Object.keys(per).sort();
  const cols = [
    { w: 900, label: "Section" },
    { w: 700, label: "Platform" },
    { w: 600, label: "Nodes" },
    { w: 1200, label: "SPARTA (s)" },
    { w: 1200, label: "BANKSY-style (s)" },
    { w: 1200, label: "Squidpy (s)" },
    { w: 900, label: "Peak MEM (MB)" },
  ];
  const rows = [
    new TableRow({ cantSplit: true, children: cols.map(c => headerCell(c.label, c.w)) }),
    ...slides.map(s => {
      const d = per[s];
      return new TableRow({
        cantSplit: true,
        children: [
          cell(s, { w: 900, size: 18, align: AlignmentType.CENTER }),
          cell(d.platform === "visium" ? "Visium" : "1st-gen ST", { w: 700, size: 18, align: AlignmentType.CENTER }),
          cell(String(d.n_nodes), { w: 600, size: 18, align: AlignmentType.CENTER }),
          cell(d.t_sparta_core_s.toFixed(3), { w: 1200, size: 18, align: AlignmentType.CENTER }),
          cell(d.t_banksy_style_s.toFixed(3), { w: 1200, size: 18, align: AlignmentType.CENTER }),
          cell(d.t_squidpy_s.toFixed(3), { w: 1200, size: 18, align: AlignmentType.CENTER }),
          cell(d.peak_mem_mb.toFixed(1), { w: 900, size: 18, align: AlignmentType.CENTER }),
        ],
      });
    }),
    // 汇总行
    new TableRow({
      cantSplit: true,
      children: [
        cell("Median (19 sections)", { w: 900, size: 18, bold: true }),
        cell("—", { w: 700, size: 18, align: AlignmentType.CENTER }),
        cell("370–2673", { w: 600, size: 18, align: AlignmentType.CENTER }),
        cell("0.025", { w: 1200, size: 18, align: AlignmentType.CENTER, bold: true }),
        cell("0.271", { w: 1200, size: 18, align: AlignmentType.CENTER }),
        cell("4.091", { w: 1200, size: 18, align: AlignmentType.CENTER }),
        cell("≤ 15.9", { w: 900, size: 18, align: AlignmentType.CENTER }),
      ],
    }),
  ];

  return new Table({
    width: { size: 100, type: WidthType.PERCENTAGE },
    columnWidths: cols.map(c => c.w),
    rows,
  });
}

// ════════════════════════════════════════════════════════════════
// 主文档构建
const children = [];

// ── 标题 ──
children.push(new Paragraph({
  alignment: AlignmentType.CENTER,
  spacing: { before: 0, after: 240, line: 300 },
  children: [new TextRun({
    text: "Two transport operators, one substrate: a graph model of T-cell migration and antibody penetration barriers in tumour tissue from spatial transcriptomics",
    font: FONT, size: SZ_TITLE, bold: true, color: COLOR_HEAD,
  })],
}));

// ── Highlights ──
children.push(new Paragraph({
  spacing: { before: 200, after: 120 },
  children: [new TextRun({ text: "Highlights", font: FONT, size: SZ_H2, bold: true, color: COLOR_HEAD })],
}));
const highlights = [
  "Graph operators turn spatial transcriptomics into transport barrier quantities.",
  "Minimum cut yields the T-cell barrier; a screened Poisson field yields the IgG barrier.",
  "The antibody barrier is percolation-limited: 60\u201365% of matrix edges exclude IgG.",
  "Both barriers share a matrix origin and co-move in 18 of 19 tumour sections.",
  "Deterministic CPU implementation analyses one section in under 0.2 s.",
];
highlights.forEach(h => children.push(bullet(h)));

// ── Graphical abstract ──
children.push(new Paragraph({
  spacing: { before: 300, after: 120 },
  children: [new TextRun({ text: "Graphical abstract", font: FONT, size: SZ_H2, bold: true, color: COLOR_HEAD })],
}));
children.push(...figImage("graphical_abstract_portrait", 3.0, 0.423, null));

children.push(new Paragraph({ children: [new PageBreak()] }));

// ── Abstract ──
children.push(new Paragraph({
  spacing: { before: 200, after: 120 },
  children: [new TextRun({ text: "Abstract", font: FONT, size: SZ_H2, bold: true, color: COLOR_HEAD })],
}));
children.push(p("Monoclonal antibodies and cytotoxic T cells must both cross the tumour extracellular matrix to reach their targets, yet they differ by three orders of magnitude in size, and the analytical tools of spatial transcriptomics return labels rather than transport quantities. We introduce two deterministic graph operators that share one spatial graph and differ only in edge-weight semantics. T-cell migration is modelled as a source\u2013sink minimum cut between immune-entry and tumour-core compartments, returning a barrier strength together with its blockade geometry. Antibody transport is modelled as a screened Poisson diffusion\u2013absorption field whose conductances implement size exclusion against the IgG hydrodynamic radius (5.5 nm) and whose sinks represent target-antigen binding. Applied to 19 sections from 7 patients spanning two cutaneous tumour types and two spatial-platform generations, the model shows that the antibody barrier is percolation-limited \u2014 the effective mesh size falls below 5.5 nm on 60\u201365% of edges \u2014 and that size exclusion, not antigen availability or bulk matrix density, accounts for a median 97.5% of barrier variance (at the default scale parameter \u03b2 = 3; full dependence on \u03b2 in Section 3.2). Ablating the matrix channel shared by both operators leaves their spatial coupling positive and significant in 12 of 15 squamous sections (median retention 81%), indicating that the two delivery problems are mathematically distinct but physically co-localised in the tissue. A molecular size scan on the same tissue confirms the distinction: the identical graph is nearly transparent to a 0.5 nm solute and strongly obstructive to an IgG. In silico reduction of matrix density and crosslinking by 30% lowers both barriers by a median 61% and 50% respectively, in 19 of 19 sections across both tumour types. The full analysis of one section runs in under 0.2 s on CPU with no learned parameters, making the model usable as a reproducible component of translational pipelines."));

children.push(pRuns([
  new TextRun({ text: "Keywords: ", font: FONT, size: SZ_BODY, bold: true }),
  new TextRun({ text: "spatial transcriptomics; graph algorithms; minimum cut; screened Poisson equation; antibody transport; tumour microenvironment", font: FONT, size: SZ_BODY, italics: true }),
], { after: 200 }));

// ════ 1. Introduction ════
children.push(h1("1. Introduction"));

children.push(p("Solid tumours are difficult places to deliver things to. A therapeutic IgG is a 150 kDa globular protein with a hydrodynamic radius near 5.5 nm; a cytotoxic T cell is an object of roughly 10 \u00b5m. Both must arrive from the vasculature, through interstitial matrix whose collagen density and crosslinking vary strongly over distances comparable to a few spot diameters of a spatial assay, before either antigen binding or tumour-cell killing can occur. The physical side of this problem has a long modelling tradition \u2014 compartmental pharmacokinetics, Krogh-cylinder penetration models, and reaction\u2013diffusion descriptions of macromolecule transport through matrix \u2014 but those models are formulated at the tissue or organ level with effective parameters, and they are not driven by molecular measurements of a specific patient section. The molecular side now produces such measurements: spatial transcriptomics resolves, at 10\u2013100 \u00b5m pitch, the local abundance of collagen and crosslinking-associated transcripts, of stromal and immune signatures, and of the therapeutic target itself."));

children.push(p("What spatial transcriptomics has lacked is an operator vocabulary that turns those measurements into transport quantities. The prevailing analyses return labels: spatial-domain segmentations and their embeddings, neighbourhood enrichment scores, cell-type maps. These answer what is next to what; they do not return a conductance, a capacity, or a barrier with a source and a sink. Even where graph-cut machinery has entered the field, it has entered as a clustering loss: the tissue-cellular-neighbourhood method of Doiron et al. [7] and the Spatial-RGCN domain identifier [8] optimise min-cut objectives to produce assignments, and the cut is discarded once labels are obtained. The transport question \u2014 which cross-section of the matrix limits flux from a specified source to a specified sink, and by how much \u2014 is not asked."));

children.push(p("We ask it, and we ask it for two delivery modalities at once, because the two have a clinically important relationship. If the T-cell barrier and the antibody barrier were spatially dissociable, matrix-directed intervention could be timed or dosed to relieve one modality selectively; if they are co-located because both are expressions of the same matrix substrate, then matrix normalisation improves both classes of delivery simultaneously and neither is reachable independently. This is a quantitative question about transport in measured tissue, which is precisely what a graph-transport model can answer and a label-based analysis cannot."));

children.push(p("The model \u2014 we refer to the implementation as SPARTA \u2014 is built on one spatial graph per section and two operators. For cell migration, immune-entry nodes (endothelial-rich) are the source and tumour-core nodes (malignant-rich, immune-poor) are the sink; edge capacities decrease with the local matrix score, and the barrier is the reciprocal of the maximum source\u2013sink flow, with the minimising cut returned as an explicit blockade geometry. For antibody transport, the same graph carries a screened Poisson diffusion\u2013absorption operator: conductances implement size exclusion through an effective mesh size that contracts with the local crosslinking score, and absorption represents target-antigen binding. The two operators share no functional form, but they share one input \u2014 the matrix \u2014 and the paper's central experiment removes that shared input entirely and asks whether the coupling survives."));

children.push(p("Three design decisions follow from the computational setting. Parameters are partitioned into physically anchored values (the IgG radius, the oxygen diffusion limit), qualitative scale parameters that are reported with their full sensitivity dependence, and sensitivity parameters scanned on grids (Section 2.4); no parameter was fitted to any outcome, clinical or otherwise. Every inferential statistic is evaluated with the patient, not the section, as the independent unit \u2014 19 sections come from 7 patients, and counts are reported at both levels throughout. And every claim the operators generate is tested by counterfactual computation on the same tissue: spatial rearrangement of the resistance field, removal of the blockade as a contiguous arc versus scattered fragments, and a molecular size sweep from small molecule to large antibody."));

children.push(p("Applied to 19 sections (15 primary cutaneous squamous cell carcinomas from 6 patients, 4 metastatic melanoma deposits from 1 patient; two platform generations of spatial transcriptomics four years apart), the model yields four findings. The antibody barrier is percolation-limited: at the default scale parameters the effective mesh size is below the IgG radius on 60\u201365% of edges, and size exclusion accounts for a median 97.5% (at \u03b2 = 3) of antibody-barrier variance. The two barriers are mathematically distinct \u2014 the same tissue graph is nearly transparent to a 0.5 nm solute and strongly obstructive to an IgG \u2014 but physically coupled: after controlling for distance from vasculature, the per-spot barrier fields are positively correlated in 18/19 sections, and the dissociation of the two barriers is rarer than chance in 18/19. Removing the shared matrix input entirely leaves the coupling significant in 12/15 squamous sections with a median 81% retention, so the coupling is substantially tissue-borne rather than an artefact of shared edge weights. And the blockade behaves as a connected structure rather than a pile of resistant spots: opening a contiguous gap in the cut band leaves ~15% more residual barrier than removing the same material at scattered positions, with a dose\u2013response over removal fraction in 15/19 sections."));

children.push(p("We frame this as a methodological contribution on public data, not as a clinical biomarker: no parameter was set using clinical outcome labels, and the biological claims are bounded throughout by the patient counts reported. The deliverable is a well-defined, deterministic, reproducible transport operator that turns spatial measurements into barrier quantities, plus a testable prediction about what stromal modification would do to both delivery modalities simultaneously."));

children.push(p("Our contributions are:"));
children.push(numItem("A pair of interpretable, deterministic graph-transport operators for spatial molecular data \u2014 a source\u2013sink minimum-cut barrier for cell migration and a screened Poisson diffusion\u2013absorption barrier for macromolecule transport \u2014 connected by an explicit parameter-identity discipline (anchored / scale / sensitivity);"));
children.push(numItem("A counterfactual validation framework (rearrangement, blockade continuity, molecular size sweep, shared-input ablation) applicable to any operator claiming to measure transport in tissue;"));
children.push(numItem("An application across two platform generations quantifying the coupling of the two delivery barriers, with a direct implication for the design of matrix-directed combination therapies and for the choice of antibody- vs cell-based modalities in matrix-rich tumours;"));
children.push(numItem("A reference implementation that analyses a full section in 0.02\u20130.19 s on CPU with no learned parameters."));

// ════ 2. Materials and methods ════
children.push(h1("2. Materials and methods"));

children.push(h2("2.1. Data and admission"));
children.push(p("Nineteen sections from two public cohorts were analysed: fifteen treatment-naive primary cutaneous squamous cell carcinomas (CSCC01\u2013CSCC16 minus CSCC13) from GSE144239 [2] and four extracranial metastatic deposits of cutaneous melanoma (MEL01\u2013MEL04) from GSE250636 [3]. The four melanoma sections are separate deposits (sternum, cecal nodule, chest wall, ribcage) of a single patient; the second patient in GSE250636 contributes only leptomeningeal deposits, a different anatomical compartment, and was not included (Section 4.4). Patient assignment follows the GEO sample metadata and is recorded with replicate structure in the analysis ledger; all counts below are reported at both section and patient level, and n = 19 is not treated as a sample size."));
children.push(p("Sections passed seven pre-declared admission criteria (tumour-content fraction, minimum spots, minimum median UMI, detected genes, endothelial signal, treatment status and site) recorded before analysis. Two cohort-level threshold relaxations for the 2016- and 2020-generation platforms, with dates and justifications, were fixed before ingestion; no section was admitted by a run-time override. One section (CSCC13) failed the relaxed median-UMI threshold (289.5 vs 300) and was excluded; its admission record is retained. Per-section outcomes are in Table 1."));
children.push(p("The squamous cohort deliberately spans two platform generations: four 2020-generation Visium sections (2 patients) and eleven 2016-generation first-generation ST sections on a staggered 200 \u00b5m array (4 patients, three technical replicates each). Agreement across that platform shift is used as a robustness argument throughout."));

children.push(h2("2.2. Spatial graphs"));
children.push(p("Coordinates were rescaled so that the median nearest-neighbour distance equals the platform's nominal pitch (100 \u00b5m Visium, 200 \u00b5m first-generation ST). The rescaling matters: the first-generation array is staggered (array indices satisfy x + y \u2261 0 mod 2), so nearest neighbours lie on the diagonal and one index unit corresponds to 141.4 \u00b5m, not 200 \u00b5m; assuming otherwise inflates every micrometre quantity by \u221a2. Graphs were built by radius adjacency (radius 150 \u00b5m Visium, 300 \u00b5m first-generation ST; the latter connects both 200 \u00b5m and 282.8 \u00b5m neighbours; mean degree 4.9\u20137.5 after quality control). Spots below 500 UMI were removed, which fragments the shallow sections; we report the fragmentation rather than assume connectivity: analysed graphs have 1\u201354 components with 90.0\u2013100% of nodes in the largest, and 8/19 sections have some source or sink nodes outside the largest component, which biases the cellular barrier slightly upward there. The operator counts and reports the affected nodes per section."));

children.push(h2("2.3. Signature scores"));
children.push(p("Inputs to both operators are within-section signature scores (rank-normalised within each section): endothelial and T/NK signatures defining source and sink compartments; a malignant signature defining the tumour core; CAF (cancer-associated fibroblast), core-matrisome and crosslinking-associated signatures defining resistance and conductance; a two-gene target-antigen signature for anti-PD-1 agents (CD274, PDCD1LG2); and hypoxia (MSigDB v7.1 HALLMARK_HYPOXIA [8]), proliferation and efflux signatures for the metabolic barrier. Scores are not deconvolved proportions; consequences are discussed in Section 4.4."));

children.push(h2("2.4. The two operators and parameter identity"));
children.push(pRuns([
  new TextRun({ text: "Cell migration \u2014 source\u2013sink minimum cut. ", font: FONT, size: SZ_BODY, bold: true }),
  new TextRun({ text: "Endothelial-rich nodes form the source (immune entry), malignant-rich immune-poor nodes the sink; quantile thresholds (scan class below) define membership, with a border fallback when a compartment is empty. Edge capacities are decreasing functions of the local matrix resistance (shared core-matrisome and CAF scores). The section-level barrier B_cell is the reciprocal of the maximum source\u2013sink flow, computed exactly (preflow-push) on the largest component; the minimising cut set is returned as a blockade band with width distribution (6.7\u201359.1% of spots across sections). The per-spot field used in correlation analyses is the shortest-path accumulated migration cost from immune-entry nodes; the two quantities are reported for their respective purposes and are not interchangeable (Section 4.4).", font: FONT, size: SZ_BODY }),
]));
children.push(pRuns([
  new TextRun({ text: "Antibody transport \u2014 screened Poisson diffusion\u2013absorption. ", font: FONT, size: SZ_BODY, bold: true }),
  new TextRun({ text: "Vessel nodes hold unit concentration; edge conductance decreases with the local matrix score and implements size exclusion through an effective mesh size \u03be = \u03be\u2080\u00b7exp(\u2212\u03b2\u00b7x), where x is the within-section rank-normalised crosslinking score and r = 5.5 nm the IgG hydrodynamic radius: an edge whose \u03be falls below r contributes only a conductance floor. Target-antigen expression acts as a distributed absorption sink. The steady state solves a screened Poisson (Helmholtz-type) system on the graph; the barrier B_mAb is the mean concentration deficit over the tumour-nest core, and the per-spot deficit field is the correlate used in Section 3.4. With absorption removed the operator reduces to the harmonic (effective-resistance) limit [11]. Variance of B_mAb is decomposed by single-channel ablation into a matrix channel (shared with the cut), a size-exclusion channel (crosslinking) and an antigen channel.", font: FONT, size: SZ_BODY }),
]));
children.push(pRuns([
  new TextRun({ text: "Parameter identity. ", font: FONT, size: SZ_BODY, bold: true }),
  new TextRun({ text: "Parameters fall into three classes. Physically anchored: IgG hydrodynamic radius r = 5.5 nm; oxygen diffusion limit d\u2080 = 130 \u00b5m. These take literature values and are never fitted. Qualitative scale parameters: \u03be\u2080 = 20 nm and \u03b2 = 3 set the steepness of the exclusion law. They are not mesh-size measurements: because x is rank-normalised, the model fixes the mesh range to [\u03be\u2080e^\u2212\u03b2, \u03be\u2080] and the complete-exclusion rank to ln(\u03be\u2080/r)/\u03b2 = 0.430 identically in every section \u2014 a constancy verified in the data (median effective mesh size 4.36\u20134.68 nm; excluded edge fraction 60.0\u201365.2% across two tumour types, two platform generations and a 29-fold depth span). Every claim depending on \u03be\u2080 or \u03b2 is reported with its dependence over a \u03b2 grid (Section 3.2). Sensitivity parameters: the four compartment quantile thresholds, the graph radius, the effective binding rate, and the metabolic weights, each scanned on a grid. No parameter of any class was set using clinical response labels.", font: FONT, size: SZ_BODY }),
]));

children.push(h2("2.5. Statistics"));
children.push(p("Permutation tests use 500 permutations per section and mode (identical across sections; empirical resolution 1/501). Multiple testing is controlled by Benjamini\u2013Hochberg within each test family (one test per section; the rearrangement modes and removal fractions form separate families). Effect sizes are reported as ratios alongside z-scores, because the permutation nulls have small variance and z-scores overstate practical magnitude. Discordant fractions are compared against their chance expectation (1 \u2212 q)\u00b2 = 6.25% at the quartile threshold. All correlations are Spearman; partial correlations residualise both fields on distance to the nearest vessel (quadratic, in rank space), where that distance is the weighted shortest path in micrometres \u2014 not a hop count: the hop approximation is biased in opposite directions on the two lattices (median \u00d71.109 Visium, \u00d70.894 first-generation ST), and its effect is reported as a sensitivity analysis wherever it changes anything."));

// ════ 3. Results ════
children.push(h1("3. Results"));

children.push(h2("3.1. Cohort and graphs"));
children.push(p("Table 1 summarises the cohort: 19 sections, 7 patients, 2 tumour types, 2 platform generations; 370\u20132 673 nodes per section after quality control; median UMI per spot spanning 567\u201316 686 (a 29-fold span used deliberately as a robustness axis). Graph construction differs between platforms only in geometry (Section 2.2), so operator outputs are comparable across the platform shift."));

// Table 1
children.push(new Paragraph({
  spacing: { before: 200, after: 80 },
  children: [new TextRun({ text: "Table 1. Sections analysed.", font: FONT, size: SZ_CAP, bold: true, color: COLOR_MUTED })],
}));
children.push(buildTable1());
children.push(emptyP());

// Fig 1
children.push(...figImage("fig1_framework", 5.5, 1.45,
  "Figure 1. Framework: one spatial graph, two edge-weight semantics, two transport operators (B_cell and B_mAb) sharing a matrix substrate."));

children.push(h2("3.2. The antibody barrier is percolation-limited"));
children.push(p("Single-channel ablation decomposes B_mAb variance into the matrix channel shared with the cut, the size-exclusion channel, and the antigen channel. Across the 19 sections, size exclusion accounts for 95.9\u201398.9% (median 97.5%) and the shared matrix channel for 0.8\u20132.8% (median 1.8%); the two inputs are not collinear in any section (|r| = 0.010\u20130.512). On the 15 sections where the two-gene target-antigen signature is detectable at all, the antigen channel accounts for 0.47\u20133.06% (median 0.96%); on the four shallowest first-generation sections neither gene passes detection, which is missing data and is reported as such rather than as zero."));
children.push(p("The mechanism is measured directly. At default parameters the effective mesh size falls below the IgG radius on the median edge of every section (median \u03be = 4.36\u20134.68 nm; 60.0\u201365.2% of edges size-excluded, median 62.7%), and only 0.0\u20137.8% of spots (median 1.0%) saturate the numerical conductance floor. The antibody barrier is therefore not a smooth gradient but a largely percolation-limited field: transport proceeds through the minority of edges whose local mesh remains open to a 5.5 nm solute."));
children.push(p("This dominance is a property of the scale parameter and we say so. Sweeping \u03b2 from 0.25 to 12 at \u03bb = 3 on five sections spanning both tumour types and both platform generations, the size-exclusion share of B_mAb variance moves from 0.2\u20130.3% (\u03b2 = 0.25) through 50.3\u201375.0% (\u03b2 = 1.5) to 95.9\u201397.7% (\u03b2 = 3) and 81.0\u201396.0% (\u03b2 = 12), with the five sections agreeing within a few points at every value. The defensible claim is not an empirical one about skin cancer but a conditional one about the model: given an exclusion law steep enough that mid-range crosslinking already contracts the mesh below the IgG radius \u2014 which at \u03b2 = 3 it does on ~63% of edges \u2014 the barrier is dominated by that axis rather than by bulk matrix density or antigen availability. The antigen channel is weak under every parameter setting examined (0.3\u201317.9% over the full grid); with a two-gene signature this is a detection floor, not evidence that binding sites are unimportant (Section 4.4)."));

children.push(p("Permutation control for the rank-normalised input. Because the crosslinking score is rank-normalised within each section, the fraction of edges below the IgG radius is partly fixed by construction. To separate model artefact from data signal, we permuted the crosslinking score 50 times within each of six representative sections (spanning both tumour types and both platform generations) and re-ran the variance decomposition on each permutation. The size-exclusion share of B_mAb variance is essentially unchanged: median 97.6% under the real crosslinking field versus 97.6% under permutation (95% null interval 96.8\u201398.4%); only 2 of 6 sections show a permutation p < 0.05, and in those the effect size is less than one percentage point. We therefore report the 97.5% figure explicitly as a property of the operator given \u03b2 = 3, not as an empirical measurement of crosslinking architecture \u2014 a framing already stated in Section 2.4 and confirmed here rather than softened after review."));

children.push(pRuns([
  new TextRun({ text: "Graph-radius sensitivity. ", font: FONT, size: SZ_BODY, bold: true }),
  new TextRun({ text: "Because the spatial graph is constructed by a hard radius threshold, we tested whether the main conclusions depend on that choice. We re-graphed all 19 sections at four Visium radii (100, 150, 200, 250 \u00b5m) and three first-generation-ST radii (200, 300, 400 \u00b5m), re-defined source/sink/vessel sets at each radius, and recomputed both operators (65 configurations total). The median partial coupling between the B_cell field and B_mAb is +0.225 across configurations (96.9% positive, versus +0.217 at default). The median B_cell reduction after 30% stromal scaling is 59.6% and B_mAb reduction 49.5% (versus 60.8% and 49.5% at default). Both barriers fall simultaneously in 89.2% of configurations. The coupling, the stromal co-targeting prediction, and the direction of the size-exclusion result are therefore not artifacts of a single graph radius.", font: FONT, size: SZ_BODY }),
]));

// Fig 2
children.push(...figImage("fig2_driver_decomposition", 6.0, 2.65,
  "Figure 2. The antibody barrier decomposes onto one axis \u2014 whose dominance is a modelling choice, not a measurement. (a) Three-channel share of B_mAb variance (log axis). (b) Crosslinking share as a function of \u03b2."));

children.push(h2("3.3. Molecular size separates the two transport problems"));
children.push(p("Sweeping the hydrodynamic radius from 0.5 nm to 10 nm at fixed tissue structure raises the mean barrier in the tumour-nest core monotonically in every section (MEL01: 1.48 \u2192 13.20; CSCC03: 2.86 \u2192 21.56; CSCC04: 2.40 \u2192 19.29 across eight radii, no saturation at the numerical ceiling). The identical graph is nearly transparent to a small molecule and strongly obstructive to an IgG-sized molecule. This is the cleanest demonstration that the two operators pose different physical problems even though they are driven by the same substrate \u2014 a distinction no cell-migration model produces, and one that is a deterministic consequence of the operator rather than a statistical test."));

// Fig 3
children.push(...figImage("fig3_size_scan", 6.0, 2.48,
  "Figure 3. Molecular size scan and percolation statistics. (a) Mean B_mAb in tumour core vs hydrodynamic radius for all 19 sections. (b) Percentage of edges with mesh size below IgG radius."));

children.push(h2("3.4. The two barriers do not dissociate"));
children.push(p("The central question is whether the two barriers are spatially separable in real tissue. They are not. After residualising both per-spot fields on distance to the nearest vessel, they remain positively correlated in 18/19 sections and significantly so in 17 (median partial Spearman \u03c1 = +0.217; all seven patients positive at patient level: +0.141 to +0.323). Within the squamous cohort the two platform generations agree (Visium median +0.270, first-generation ST +0.207). The dissociation zone \u2014 lowest quartile of the cell barrier, highest quartile of the antibody barrier \u2014 occupies a median 3.7% of spots against a 6.25% chance expectation, below chance in 18/19 sections: discordance is rarer than random, a second independent expression of the coupling."));

children.push(p("Part of this association is guaranteed by construction \u2014 the core-matrisome score enters both edge-capacity and conductance \u2014 so the decisive experiment removes the shared input entirely, recomputing the cellular barrier from the fibroblast signature alone and the antibody barrier from crosslinking and antigen alone, such that the two operators share no input variable. Across 19 sections the median partial correlation falls from +0.217 to +0.152; 17/19 remain positive, 15/19 significantly. Stratified:"));

// Table 2
children.push(new Paragraph({
  spacing: { before: 200, after: 80 },
  children: [new TextRun({ text: "Table 2. R3b stratified coupling after shared-input removal.", font: FONT, size: SZ_CAP, bold: true, color: COLOR_MUTED })],
}));
children.push(buildTable2());
children.push(emptyP());

children.push(p("The squamous coupling survives on both platform generations, under both the weighted and the hop-count geometric control (5/6 patients either way), and is essentially independent of sequencing depth (\u03c1 = +0.19) and section size (\u03c1 = +0.03). We read it as organised obstruction: a desmoplastic front in which the fibroblast band, its deposited collagen, and the crosslinking that contracts the mesh are one physical structure, so the cellular and macromolecular obstructions coincide because they are the same wall. The melanoma stratum cannot be decided by this design \u2014 it is one patient, and the only stratum whose classification flips with the distance convention \u2014 and we report it as unresolved rather than as a contrast. An earlier eight-section version of this analysis asserted a squamous-versus-melanoma split on the strength of 4/4 versus 1/4; the squamous side does not survive extension to 15 sections (12/15), and we withdraw the contrast."));
children.push(p("The negative half of the conclusion is robust to everything tested: the dissociation-zone fraction stays below chance in 18/19 before the shared-input removal and 19/19 after it, under both distance conventions, across both platform generations, independent of depth and size. The two barriers do not dissociate."));

// Fig 4
children.push(...figImage("fig4_coupling_forest", 6.0, 2.88,
  "Figure 4. The two barriers do not dissociate. (a) Partial \u03c1 per section. (b) Forest plot: before vs after shared-input removal. (c) Dissociation zone vs chance."));

children.push(h2("3.5. Counterfactual experiments"));
children.push(pRuns([
  new TextRun({ text: "Blockade continuity. ", font: FONT, size: SZ_BODY, bold: true }),
  new TextRun({ text: "Removing resistance material from a contiguous arc of the min-cut band, compared with removing the same material from scattered positions within the same cut set (same material, same amount, only arrangement differs), leaves the residual barrier higher by a factor of 1.038\u20131.390 (median 1.147) at the pre-specified 20% removal fraction; 16/19 sections are significant after per-section correction, five of seven patients have every section significant, and no patient has none. The effect shows a dose\u2013response over removal fraction (residual barrier at 30% > 5% in 15/19 sections) \u2014 the falsifiable prediction of a connected blockade and not the behaviour of uncorrelated noise. We regard the dose\u2013response as the stronger evidence and the single-level significances as secondary. The effect is real but modest (real bands are thick and redundant; a synthetic one-spot-wide ring gives 7.45\u00d7 under the same experiment), and we do not rest the central argument on it.", font: FONT, size: SZ_BODY }),
]));
children.push(pRuns([
  new TextRun({ text: "Spatial rearrangement. ", font: FONT, size: SZ_BODY, bold: true }),
  new TextRun({ text: "Permuting spot positions while keeping composition fixed moves the section-level cellular barrier above its permutation null by 0.89\u00d7\u20134.49\u00d7, significant in 12/19 sections, organised by patient (three patients fully significant, two fully non-significant). We report this as a positive control, not as evidence: the effect ratio correlates with sequencing depth (\u03c1 = +0.72) and graph connectivity (\u03c1 = +0.46), which are themselves correlated, and the non-significant patients sit at the shallow, fragmented end of the range. Composition-based measures cannot move under this null at all \u2014 their permutation nulls are degenerate \u2014 which is the point of difference: a barrier is a statement about arrangement, and only an arrangement-sensitive quantity can be evidence about it.", font: FONT, size: SZ_BODY }),
]));

children.push(pRuns([
  new TextRun({ text: "In silico stromal co-targeting. ", font: FONT, size: SZ_BODY, bold: true }),
  new TextRun({ text: "The coupling implies a directly testable prediction: if the two barriers share a matrix substrate, reducing matrix density and crosslinking in silico should lower both barriers at once. We tested this by scaling the ecm, caf and crosslinking scores by (1 \u2212 reduction) at reduction levels of 20%, 30% and 50%, and recomputing both operators. At a 30% reduction, the section-level minimum-cut barrier falls by a median 60.8% (cSCC 61.1%, melanoma 58.1%) and the mean tumour-core B_mAb by 49.5% (cSCC 49.5%, melanoma 49.6%); both barriers fall in 19 of 19 sections. At 50% reduction the drops are 80.2% and 71.3% respectively. The two modalities respond in the same direction and with similar magnitude across both tumour types, which is the quantitative basis for the paper\u2019s clinical claim: stromal-directed intervention (LOX/TGF-\u03b2 inhibition, anti-fibrotic combinations) is predicted to improve cellular and macromolecular delivery simultaneously rather than relieving one at the expense of the other. This is a model-based prediction, not a measured treatment response; it is stated as such.", font: FONT, size: SZ_BODY }),
]));

// Fig 5
children.push(...figImage("fig5_counterfactual", 6.0, 2.48,
  "Figure 5. Counterfactual experiments. (a) S1 spatial rearrangement by patient. (b) S2 blockade continuity dose\u2013response."));

children.push(h2("3.6. The barrier is not a rewrite of density, and its relation to domain tools"));
children.push(p("The per-spot cellular field correlates only moderately with the CAF signature (\u03c1 = 0.313\u20130.554, median 0.458) and essentially not at all with the T/NK signature (median \u22120.018), and the commonly used immune-exclusion proxy \u2014 mean distance from T-cell-rich spots to the tumour core \u2014 changes sign across replicate sections of the same tumour type (permutation z from \u22122.4 to +5.9, negative in 7 sections), which is one motivation for defining the barrier as a transport quantity rather than a distance summary."));
children.push(p("Relation to domain segmentation follows the same logic. A BANKSY-style segmentation [4] recovers part of the structure (min-cut edges are enriched on domain boundaries by 1.10\u20131.87\u00d7, median 1.48\u00d7, over the domain-count grid k = 4\u201312), but only 4\u201318% of domain-boundary edges lie on the min cut: a segmentation draws a boundary wherever expression changes, while the transport formulation selects the one cross-section that limits flux between a specified source and sink and returns its capacity. Four properties follow that a domain boundary does not have: a capacity, an orientation, a counterfactual handle, and modality transfer. Where cuts have appeared in spatial pipelines before \u2014 as clustering losses in tissue-neighbourhood assignment [7] and domain identification [8] \u2014 the cut optimises an assignment and is discarded; here it is the quantity of interest. Computational cost is part of the comparison: on the same graphs, the full five-stage SPARTA analysis runs in 0.02\u20130.19 s per section (entire cohort 0.9 s, peak memory < 16 MB, commodity CPU), against 0.16\u20131.20 s for BANKSY-style segmentation and 3.8\u20135.0 s for neighbourhood enrichment [5] (both reference implementations run within this repository for controlled comparison; graph-learning domain methods were not run and are not claimed as comparisons)."));

children.push(pRuns([
  new TextRun({ text: "Why a graph at all? ", font: FONT, size: SZ_BODY, bold: true }),
  new TextRun({ text: "A natural reviewer question is whether a simple node-level matrix score (0.5\u00b7(ECM + CAF)) already captures the barrier, making the min-cut and screened-Poisson machinery unnecessary. We compared the per-spot graph B_cell field against this naive score on all 19 sections, using partial Spearman correlation with B_mAb after residualising on vessel distance. The graph field yields a median partial \u03c1 of +0.217 versus +0.166 for the naive score; the graph field is stronger in 14 of 19 sections and reaches p < 0.05 in 17 of 19 versus 15 of 19 for the naive score. The advantage is modest rather than transformative (median \u0394\u03c1 = +0.048) and is concentrated in sections where source\u2013sink geometry is non-trivial (e.g. CSCC11/CSCC12, where \u0394\u03c1 = +0.14/+0.17). We read this honestly: a matrix score captures most of the barrier signal because the matrix is the dominant substrate, but the graph operator adds the spatial arrangement \u2014 the capacity of the narrowest cut and the screened-diffusion field \u2014 that a node score cannot represent.", font: FONT, size: SZ_BODY }),
]));

// Table 3
children.push(new Paragraph({
  spacing: { before: 200, after: 80 },
  children: [new TextRun({ text: "Table 3. Runtime comparison: SPARTA vs BANKSY-style segmentation vs Squidpy neighbourhood enrichment (19 sections).", font: FONT, size: SZ_CAP, bold: true, color: COLOR_MUTED })],
}));
children.push(buildTable3());
children.push(emptyP());

// ════ 4. Discussion ════
children.push(pRuns([
  new TextRun({ text: "Alignment with measured cell distributions. ", font: FONT, size: SZ_BODY, bold: true }),
  new TextRun({ text: "Beyond internal consistency, we asked whether the per-spot barrier fields align with measured cell distributions in the same tissue. After residualising on vessel distance, the B_cell field shows a weak but directionally consistent negative partial correlation with the T/NK signature (median \u03c1 = \u22120.022, significant in 6/19 sections) and with the CD8 T-cell signature (median \u03c1 = \u22120.013, significant in 5/19): barrier-high spots tend to contain fewer T cells, though the effect is small because T-cell localisation is also driven by antigen availability and inflammatory cues not modelled here. The B_mAb field shows a weak positive partial correlation with proliferation (median \u03c1 = +0.046, significant in 7/19), consistent with antibody-blocked nests retaining proliferating cells. We do not overstate these alignments: they are weak, they point in the expected direction, and they provide an independent (if modest) check that the operators are not purely mathematical constructs.", font: FONT, size: SZ_BODY }),
]));

children.push(h1("4. Discussion"));

children.push(pRuns([
  new TextRun({ text: "What the model says about delivering biologics. ", font: FONT, size: SZ_BODY, bold: true }),
  new TextRun({ text: "The dominant term in the antibody barrier is not how much matrix there is, but whether the local mesh is open to the molecule: transport is percolation-limited, proceeding through the minority of edges whose effective mesh exceeds the hydrodynamic radius. Two corollaries follow. First, bulk measurements of matrix abundance are a poor proxy for penetration \u2014 two tissues with identical average matrix but different mesh topology can differ substantially in delivered fraction, which is consistent with the instability we observe when scoring bulk RNA-seq cohorts with the same signature logic (Section 4.4). Second, the model provides a per-section, parameter-swept way to ask what an intervention buys: the molecular size scan is the degenerate case, and replacing the target set (ERBB2, TACSTD2, NECTIN4) converts the antibody problem into an antibody\u2013drug-conjugate problem with a measurable binding-site channel, at no change to the operator.", font: FONT, size: SZ_BODY }),
]));
children.push(pRuns([
  new TextRun({ text: "One substrate, two modalities. ", font: FONT, size: SZ_BODY, bold: true }),
  new TextRun({ text: "After removing the only input the two operators share, their coupling retains a median 81% of its magnitude in the squamous cohort and remains significant in 5/6 patients. The practical reading for therapy design is that in matrix-rich cutaneous tumours the two delivery problems are not independently addressable: matrix normalisation is predicted to relieve both simultaneously, and modality selection (antibody versus cell) cannot route around the matrix. The corollary for measurement is that any assay meant to predict delivery must resolve arrangement, not only abundance \u2014 and the dissociation zone, the region where one modality is delivered and the other is not, occupies less than the chance fraction in every section examined.", font: FONT, size: SZ_BODY }),
]));
children.push(pRuns([
  new TextRun({ text: "Relation to continuum transport models. ", font: FONT, size: SZ_BODY, bold: true }),
  new TextRun({ text: "Macromolecule penetration has a mature continuum tradition \u2014 Krogh cylinders, reaction\u2013diffusion and pharmacokinetic models \u2014 whose parameters are effective and whose geometry is idealised. The graph formulation is complementary: it inherits the geometry of a measured section, defines barrier quantities on that geometry directly, and sacrifices physical continuity (edges are 100\u2013300 \u00b5m) for data fidelity. The size-exclusion channel is where the two views meet, since it carries the only nanometre-scale quantity (the hydrodynamic radius) in an otherwise micrometre-scale model.", font: FONT, size: SZ_BODY }),
]));
children.push(pRuns([
  new TextRun({ text: "Limitations. ", font: FONT, size: SZ_BODY, bold: true }),
  new TextRun({ text: "This is a single-author work: the analysis was not independently checked by a second analyst, and the code is released for community verification rather than having passed internal peer review. The study's principal limitation is sample size and its asymmetry: six squamous patients across two platform generations against one melanoma patient; all inferential statistics are patient-level, and no conclusion is drawn from the cross-cohort comparison. Sequencing depth spans 29-fold and is partly confounded with patient; results that track depth (the rearrangement control) are presented as controls rather than evidence. Scale parameters \u03be\u2080 and \u03b2 are not mesh measurements \u2014 the rank-normalised input fixes the exclusion fraction by construction \u2014 so the 97.5% share is a statement about the model given \u03b2, with the full dependence reported. Inputs are signature scores rather than deconvolved proportions; most acutely, the two-gene anti-PD-1 target set sits at the detection floor, so the weak antigen channel is a resolution limit, not a biological finding. The per-spot cellular field is a shortest-path accumulated cost, which permits cost-free detours and reads closer to effective depth than to blockade; the section-level cut is the topological object. No public dataset combines spatial transcriptomics with checkpoint-response labels, so the bulk concept check (one cohort significant in the expected direction, one not) is motivational, and its instability is itself the finding we report: bulk abundance discards the arrangement that makes matrix a barrier. Finally, the melanoma cohort lacks treatment annotation, so post-treatment alteration of barrier structure cannot be excluded there.", font: FONT, size: SZ_BODY }),
]));
children.push(pRuns([
  new TextRun({ text: "Outlook. ", font: FONT, size: SZ_BODY, bold: true }),
  new TextRun({ text: "Three extensions require no modification of the operators: target-set replacement for ADC-style problems; deconvolution-based inputs for a heterogeneous graph; and application to imaging-based platforms at single-cell resolution, where the 10 \u00b5m cell and the 5.5 nm molecule are separated by more than two orders of magnitude and the two barriers might yet prove separable at a resolution spot-based assays cannot reach.", font: FONT, size: SZ_BODY }),
]));

// ════ 5. Conclusions ════
children.push(h1("5. Conclusions"));
children.push(p("We formulated cell and antibody delivery in tumour tissue as two transport operators on one spatial-transcriptomic graph: a source\u2013sink minimum cut and a screened Poisson diffusion\u2013absorption field with size exclusion. The model is deterministic, interpretable, parameter-honest, and fast enough (sub-second per section on CPU) to serve as a component of translational analysis pipelines. On 19 sections from 7 patients across two platform generations it establishes three quantitative facts: the antibody barrier is percolation-limited by size exclusion; the two barriers are mathematically distinct yet physically co-located, with a coupling that survives removal of every shared input in the squamous cohort; and the blockade behaves as a connected structure. Quantitatively, a 30% in silico reduction in matrix density and crosslinking lowers the minimum-cut barrier by 61% and the antibody barrier by 50% across all 19 sections. For matrix-rich tumours the model predicts that matrix-directed intervention improves both classes of delivery at once \u2014 and that neither quantity is measurable without resolving spatial arrangement."));

// ════ 后置件 ════
children.push(h1("CRediT authorship contribution statement"));
children.push(p("Yize Li: Conceptualisation, Methodology, Software, Validation, Formal analysis, Investigation, Data curation, Writing \u2014 original draft, Visualisation."));

children.push(h1("Declaration of competing interest"));
children.push(p("The authors declare no competing interests."));

children.push(h1("Data availability"));
children.push(p("All analysed data are publicly available: spatial transcriptomic cohorts from the Gene Expression Omnibus (accessions GSE144239 [2] and GSE250636 [3]); the hypoxia signature from MSigDB v7.1 [8]. The SPARTA implementation, all analysis scripts, configuration files and per-section artefacts are available at https://github.com/Ezbenzino/sparta and archived at https://doi.org/10.5281/zenodo.23086431."));

children.push(h1("Acknowledgements"));
children.push(p("The author received no specific funding for this work."));

// ════ References ════
children.push(h1("References"));
const refs = [
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
];
refs.forEach(r => children.push(new Paragraph({
  spacing: { after: 80, line: 240 },
  alignment: AlignmentType.JUSTIFIED,
  children: [new TextRun({ text: r, font: FONT, size: SZ_SMALL })],
})));

// ════ 构建文档 ════
const doc = new Document({
  styles: {
    default: {
      document: {
        run: { font: FONT, size: SZ_BODY },
      },
    },
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: SZ_H1, bold: true, font: FONT, color: COLOR_HEAD },
        paragraph: { spacing: { before: 360, after: 180 }, outlineLevel: 0, keepNext: false, keepLines: false } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: SZ_H2, bold: true, font: FONT, color: COLOR_HEAD },
        paragraph: { spacing: { before: 240, after: 120 }, outlineLevel: 1, keepNext: false, keepLines: false } },
    ],
  },
  numbering: {
    config: [
      { reference: "bullets",
        levels: [{ level: 0, format: LevelFormat.BULLET, text: "\u2022", alignment: AlignmentType.LEFT,
          style: { paragraph: { indent: { left: 720, hanging: 360 } } } }] },
      { reference: "numbers",
        levels: [{ level: 0, format: LevelFormat.DECIMAL, text: "%1.", alignment: AlignmentType.LEFT,
          style: { paragraph: { indent: { left: 720, hanging: 360 } } } }] },
    ],
  },
  sections: [{
    properties: {
      page: {
        size: { width: 12240, height: 15840 },
        margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 },
      },
    },
    headers: {
      default: new Header({ children: [new Paragraph({
        alignment: AlignmentType.RIGHT,
        children: [new TextRun({ text: "SPARTA \u2014 CBC submission draft", font: FONT, size: 18, color: "999999" })],
      })],
    }),
    },
    footers: {
      default: new Footer({ children: [new Paragraph({
        alignment: AlignmentType.CENTER,
        children: [new TextRun({ text: "Page ", font: FONT, size: 18, color: "999999" }),
                   new TextRun({ children: [PageNumber.CURRENT], font: FONT, size: 18, color: "999999" })],
      })],
    }),
    },
    children,
  }],
});

const outPath = path.join(ROOT, "docs", "manuscript_cbc_draft.docx");
Packer.toBuffer(doc).then(buffer => {
  fs.writeFileSync(outPath, buffer);
  console.log("Done: " + outPath);
});
