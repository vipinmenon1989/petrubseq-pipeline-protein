# petrubseq-pipeline-protein

A reproducible pipeline for **Perturb-seq** and **Perturb-CITE-seq** experiments.

Start with RNA counts and CRISPR guide information; add antibody-derived tag
(ADT) protein measurements if you have them. The pipeline checks and aligns
your inputs, assigns perturbations to cells, runs RNA and protein QC and
preprocessing, and writes one processed AnnData object together with a
self-contained HTML report. Optional analyses then ask which transcriptional
states the cells fall into and which perturbations favour them, how strongly each
cell shows its perturbation's response, which genes respond together and, when
protein is measured, which surface proteins changed.

```bash
petrubseq-protein run --config config/demo_papalexi.yaml
```

## What can I use it for?

| experiment | what you provide | supported |
|---|---|---|
| Perturb-seq | RNA + guide counts (or guide calls) | yes; protein is optional |
| Perturb-CITE-seq / ECCITE-seq | RNA + guides + ADT counts or normalized protein | yes |
| perturbations already assigned upstream | RNA + a per-cell guide list or a guide/perturbation column in your metadata | yes |
| several 10x lanes | one matrix per lane plus an optional per-lane sample sheet | yes |
| input formats | dense CSV/TSV, 10x MTX, 10x H5, `.h5ad`; one combined matrix or one per modality | yes |

Protein is never required for the core pipeline or for the RNA-based analyses.
The protein analyses simply skip themselves when there is no protein.

## What questions does the pipeline answer?

The core pipeline answers the first question for every run. The others come
from optional analyses: the cell-state analyses are switched on with
`analysis.clustering.enabled: true`, the perturbation-response, program and
protein analyses with `analysis.perturbation_effects.enabled: true`.
Every perturbation comparison uses **perturbed cells** (cells with exactly one
targeting guide) and **control cells** (cells with exactly one non-targeting
guide). Cells with ambiguous, multiple or no guide calls are kept in the object
but left out of these comparisons.

| group | question | analysis |
|---|---|---|
| core | Did the experiment and guide assignment work? | QC and perturbation assignment |
| cell states | What transcriptional states exist? | Leiden clustering |
| | Does a perturbation preferentially occupy one of those states? | perturbation × cluster enrichment |
| | Does a perturbation concentrate locally in expression space, without cluster boundaries? | lochNESS |
| perturbation response | How strongly does each cell show its perturbation's response? | PS |
| programs | Which perturbations behave similarly? Which downstream genes respond together? | perturbation modules, gene programs |
| multimodal | Which proteins change? Do RNA and protein phenotypes go together? | protein effects, RNA–protein concordance |

### Did the experiment and guide assignment work?

RNA QC (genes and UMIs per cell, mitochondrial and ribosomal fractions),
protein QC when ADTs are present (ADT depth, antibodies detected, isotype
background), and guide QC (guide UMIs per cell, guides detected, top-versus-
second guide). Each cell gets a perturbation class. Every filtering step is
recorded, and QC figures are drawn before and after filtering so you can see
what was removed.

### What transcriptional states exist?

**Leiden clustering** divides the cells into discrete states from their RNA
profiles alone, using the same neighbour graph as the UMAP (guide labels play no
part). Clusters are numbered, not named, and their number depends on the chosen
resolution. The report shows each cluster's size and make-up: which perturbation
classes, samples and lanes it contains, and its library depth. A cluster made of
one sample, one lane, low-depth cells or ambiguous guide calls is flagged.
Method: Traag et al. 2019.

### Does a perturbation preferentially occupy one of those states?

**Perturbation × cluster enrichment.** For each target and cluster, compare the
fraction of the target's cells in that cluster with the fraction of control cells
there (Fisher's exact test, false-discovery-rate corrected). The result reads
directly: "65 % of JAK2 cells sit in cluster 4, versus no control cells". It
also reports how many of the target's guides agree, and can be stratified by
sample or lane so that a cluster that merely differs between samples does not
pass as a perturbation effect.

### Does a perturbation concentrate locally, without cluster boundaries?

**lochNESS.** For each cell, look at its nearest neighbours in expression
space and ask whether cells carrying perturbation X are over-represented
there compared with their overall frequency. Averaged over X's own cells, this
tells you whether X's cells cluster together in a state of their own. It
measures *where cells end up*, not how strongly the target gene was knocked
down: a perturbation can change a few genes strongly without moving cells to a
new state, and the reverse.
Method: Huang et al. 2023.

Cluster enrichment and lochNESS look at the same question from two sides and are
never combined:

| | cluster enrichment | lochNESS |
|---|---|---|
| needs discrete clusters | yes; depends on the resolution | no |
| reads as | "X is over-represented in cluster k" | "X's cells sit near each other" |
| can see | preference for a defined state | concentration inside one cluster or across two |

A perturbation can show both, either or neither. The report puts the strongest
cluster result and the lochNESS summary side by side for each target.

### How strongly does each cell show its perturbation's response?

**Perturbation-response score (PS).** A perturbation changes a set of
downstream genes. PS asks, cell by cell, how much of that change a perturbed
cell actually shows: near 0 it looks like a control, near 1 it shows the full
response. Combined with the targeted gene's own expression, it separates cells
where the knockdown worked from escapers and non-responders.

PS is normalized separately for each target, so compare PS between cells of the
same target, not between targets. To compare response strength across targets,
use `auc_vs_control` and `d_vs_control` in `ps_targets.csv`, which say how well
perturbed cells separate from controls.
Method: Song et al. 2025.

### Which perturbations produce similar responses?

**Perturbation modules.** Perturbations whose response profiles across genes
are correlated are grouped into modules (`M1`, `M2`, ...).

### Which downstream genes respond together?

**Gene programs.** Downstream genes whose responses move together across
perturbations are grouped into programs (`P1`, `P2`, ...). Each cell also gets
an activity score for every program.

Four terms are kept distinct:

| term | meaning |
|---|---|
| perturbation target | the gene a guide is designed to act on |
| response gene | a gene whose expression changes downstream of a perturbation |
| perturbation module | a group of perturbations with similar responses |
| gene program | a group of response genes that change together |

Modules and programs come from the same perturbation × gene table of fold
changes, clustered along its two axes. They are numbered, never given
biological names; the tables list their members so you can interpret them.
Design: Zhou et al. 2023.

### Which measured proteins changed?

*Requires protein/ADT measurements.*

**Protein effects.** For each perturbation and each measured protein, perturbed
cells are compared with controls on the normalized protein values. The table
reports an effect size, a p-value and a false-discovery rate. It also says
whether the direction holds in each sample and for each guide against the same
target.

### Do transcriptional and protein phenotypes agree?

*Requires protein/ADT measurements.*

**RNA–protein concordance** asks three questions:

* **PS ↔ protein.** Among the cells of one perturbation, do cells that respond
  more strongly at the RNA level have more, or less, of a protein?
* **lochNESS ↔ protein.** Across perturbations, do the ones that move cells
  into a distinct state also change a protein more? lochNESS describes a whole
  perturbed population, so this is asked across targets, not cell by cell.
* **Gene program ↔ protein.** Does a program's activity track a protein from
  cell to cell? And do perturbations that change a program also change the
  protein?

These are correlations. They show which RNA and protein signals move
together. They do not show that a program causes or mediates a protein change.

## Quick start

```bash
git clone https://github.com/vipinmenon1989/petrubseq-pipeline-protein.git
cd petrubseq-pipeline-protein

conda env create -f environment.yml
conda activate petrubseq-protein
pip install -e .

petrubseq-protein run --config config/demo_papalexi.yaml
```

The repository already contains a deterministic 1,800-cell subset of the
Papalexi et al. ECCITE-seq screen (`demo/data/papalexi_eccite/`, about 20 MB),
so nothing needs to be downloaded. The run takes about a minute and under
1 GB of memory.

Results are written next to the repository, in
`../results/demo/papalexi_eccite/`. Set `PETRUBSEQ_RESULTS_ROOT` to write them
elsewhere. Open `report.html` there first; it is a single file with every
figure embedded.

## What does a run produce?

1. **A processed AnnData object** (`.h5ad`): normalized RNA, count layer, protein
   matrices, guide counts, perturbation labels, QC columns, embeddings and, when
   enabled, the per-cell analysis results.
2. **`report.html`**: one self-contained report covering inputs, QC before and
   after filtering, perturbation assignment, embeddings, the optional analyses
   and provenance.
3. **CSV tables** for QC, filtering, perturbation coverage and every analysis.
4. **Figures**, also embedded in the report.
5. **The resolved configuration, a run log and a manifest** recording inputs,
   software versions and timings.

```
results/<run>/
├── report.html
├── report.md                      same content as Markdown, figures linked
├── processed/<name>.h5ad
├── tables/                        QC and perturbation tables
│   ├── cell_states/               clusters, their composition, perturbation × cluster enrichment
│   └── perturbation_effects/      PS, lochNESS, programs, modules, protein effects, concordance
├── figures/                       grouped by report section
└── logs/                          run.log, resolved_config.yaml, run_manifest.json
```

The slot-by-slot contents are in `docs/PROCESSED_OBJECT.md`, and every table
and figure is listed in `docs/OUTPUTS.md`.

## Demo dataset

Papalexi E., Mimitou E.P., Butler A.W., et al. *Characterizing the molecular
regulation of inhibitory immune checkpoints with multimodal single-cell
screens.* Nat Genet 53, 322–331 (2021). GEO GSE153056.

The bundled subset comes from the paper's pooled ECCITE-seq screen
(GSM4633614–GSM4633618) in IFN-γ stimulated THP-1 cells.

| | bundled subset |
|---|---|
| cells | 1,800 of 20,729; hashing singlets sampled per hashed sample × guide target, seed 0 |
| RNA | 18,649 genes, raw UMI counts |
| protein | 4 ADTs: CD86, PD-L1 (`PDL1`), PD-L2 (`PDL2`), TIM-3 (`CD366`); no isotype controls |
| guides | 111 sgRNAs: 25 gene targets with 3–4 guides each, 9 non-targeting (`NT`), 1 eGFP |
| design | 3 hashed samples over 8 lanes; every sample IFN-γ stimulated |

This is a demonstration subset, not a reanalysis of the publication. With 1,800
cells, 21–102 cells per target and four proteins, its statistics are
illustrative. Because every sample is stimulated, there is no stimulated-versus-
control comparison. Provenance is in `demo/data/papalexi_eccite/README.md`.

### What should I see in the demo?

These observations come from the small bundled subset. They are listed to show
what the outputs look like, not as findings.

* IFNGR1 and JAK2 show the strongest lochNESS self-enrichment, followed by
  SMAD4 and IFNGR2: their cells sit together in their own region of expression
  space.
* IFNGR1, IFNGR2, JAK2, STAT1, IRF1 and POU2F2 form one perturbation module.
  IFNGR1, IFNGR2, JAK2 and STAT1 lower program P1 the most (mean log2FC about
  −0.3); P1 includes CXCL9, CXCL10, SOCS1, WARS and HLA-DRA.
* JAK2, IFNGR1 and IFNGR2 cells have lower PD-L1 protein than controls. The
  direction holds in all three samples and for all four guides of each target.
* Leiden (resolution 1.0) finds 10 clusters. Cluster 4 contains no control cells
  at all; 176 of its 185 single-guide targeting cells carry JAK2, IFNGR1, IFNGR2 or STAT1
  guides, so those four targets are enriched there (every guide of each target
  agrees, and the result holds when stratified by sample). The same targets have
  the strongest lochNESS scores, so here the two views agree.
* Two clusters are flagged for review rather than interpreted: cluster 2 is 85 %
  one hashed sample, and cluster 7 has a median library size 1.5 times the rest.
* Within IFNGR1, IFNGR2 and STAT1 cells, a higher PS goes with lower PD-L1
  (Spearman ρ between −0.38 and −0.58).
* Cells expressing program P1 more also tend to have more PD-L1 (ρ 0.40). That
  association is much weaker among the control cells alone (ρ 0.13), so it
  comes mostly from differences between perturbations.

## Inputs

| data | needed |
|---|---|
| RNA counts, or log-normalized RNA | always |
| guide counts, a per-cell guide list, or a guide/perturbation column in the metadata | for perturbation assignment and every perturbation analysis |
| ADT / protein, raw counts or normalized values | optional |
| per-cell metadata (condition, sample, library size, ...) | optional; required by default, see below |
| per-lane sample sheet, antibody annotation table, precomputed embedding | optional |

Supported layouts: dense CSV/TSV (either orientation, several files that split
the cells), 10x MTX directories, 10x H5 files and `.h5ad` files. A combined
feature matrix (Gene Expression / Antibody Capture / CRISPR Guide Capture) is
read once and split by feature type. Separate matrices per modality work too,
and several lanes can be combined. `config/examples/` has one ready-to-edit
config per layout.

A minimal config for a combined 10x matrix with RNA, ADT and guides:

```yaml
dataset:
  name: my_screen
  input_dir: ${DATA_ROOT}/my_screen
inputs:
  multiplexed:
    format: mtx
    path: filtered_feature_bc_matrix
    feature_types: {rna: [Gene Expression], protein: [Antibody Capture], guide: [CRISPR Guide Capture]}
    var_names: name
  metadata: {file: cells.csv, format: auto}
columns:
  cell_id: barcode
  condition: condition
perturbation:
  guide_target_regex: "^(?P<target>.+?)_(?P<index>\\d+)$"
  control_classes: {non_targeting: ["^NTC$"]}
analysis:
  clustering:
    enabled: true
  perturbation_effects:
    enabled: true
output:
  dir: ${RESULTS_ROOT}/my_screen
```

The defaults assume protein and a metadata table are present. For RNA and
guides only, with no protein and no metadata, set:

```yaml
inputs:
  protein: {required: false}
  metadata: {required: false}
alignment:
  required: [rna]
```

(and leave out the `protein` feature type and the `metadata` and `columns`
entries above).

## Guide assignment

There are two routes, chosen with `perturbation.assignment.source`:

* **Provided assignments.** A per-cell guide list (`inputs.guide_assignments`)
  or a guide/perturbation column in the metadata.
* **Called from guide counts.** The `dominant` rule (default) gives a cell its
  top guide when that guide has at least 3 UMIs and more than twice the
  runner-up. Otherwise the cell is `ambiguous`, or `unassigned` when no guide
  reaches 3 UMIs. The `threshold` rule assigns every guide with at least
  3 UMIs, for high-MOI screens with several guides per cell.

The default `auto` uses whichever source exists. If both exist, the run stops
and asks you to set `provided` or `guide_counts`. Either way both are kept, and
the report shows how well they agree.

Each cell ends up in one class:

| class | meaning |
|---|---|
| `single_targeting` | one guide, against a gene |
| `single_control` | one control guide (e.g. non-targeting) |
| `multi_targeting` | several guides, all targeting |
| `multi_control` | several guides, all controls |
| `mixed_control_targeting` | targeting and control guides together |
| `ambiguous` | guide counts present but no dominant guide (dominant rule only) |
| `unassigned` | no guide |

Guides are mapped to targets with a regular expression or an explicit table.
Controls are recognised by patterns you set (`perturbation.control_classes`).

## RNA and protein processing

**RNA.** Raw counts or log-normalized input are detected automatically.
Raw counts are normalized once, after filtering (library-size normalization to
10,000, then log1p). Log-normalized input is left as it is. Then come up to
3,000 highly variable genes, PCA (50 components), a 15-neighbour graph and UMAP.

**Protein (optional).** Raw ADT counts are kept. The default representation
is a per-protein CLR across cells. Normalized protein values supplied with
the data are preserved, and never treated as counts. Isotype controls, when
present, are used for background QC and left out of the protein PCA and UMAP.

## Analysis methods

Perturbed cells are single-guide targeting cells of one target. Controls are
single-guide non-targeting cells. Full definitions, parameters and limitations
are in `docs/CELL_STATES.md` and `docs/PERTURBATION_EFFECTS.md`.

**Leiden clustering** (Traag et al. 2019). `scanpy.tl.leiden` with the igraph
backend on the existing RNA neighbour graph (15 neighbours on the RNA PCA),
resolution 1.0, 2 iterations, fixed seed. All QC-passing cells are clustered,
since that is the population the graph was built on.

**Perturbation × cluster enrichment.** For each target (≥ 10 single-guide cells)
and cluster (≥ 20 cells), a two-sided Fisher exact test compares the target's
cells in and out of the cluster with the controls' cells in and out of it.
The table reports:

* the odds ratio as Fisher returns it (0 or ∞ when a count is zero), and a
  finite Haldane-corrected log2 odds ratio for ranking and the heatmap;
* a Benjamini–Hochberg FDR over every tested target × cluster pair of the run;
* a low-power flag when the cluster holds fewer than 10 control cells;
* for each guide with ≥ 5 cells, whether it falls on the same side of the
  control fraction as the target.

An optional Cochran–Mantel–Haenszel test across a sample or lane column
(`stratify_by`) is reported beside the Fisher result, not instead of it.

**PS** (Song et al. 2025; implementation follows PS_python). For each target:

1. Take the 100 genes with the highest t-statistic in perturbed versus control
   cells (the most increased genes).
2. Estimate each gene's shift under the perturbation by least squares.
3. Project every cell onto that signature.
4. Centre the scores on the controls, clip them to [0, 3] and divide by the
   maximum, giving [0, 1].

Targets whose gene is expressed in under 1 % of control cells are skipped, because
knockdown cannot be measured there. The raw projection is fitted and scored on the
same cells, so even a perturbation without effect separates slightly from the
controls.

**lochNESS** (Huang et al. 2023; implementation ported from pertTF). For every
cell and target: `lochNESS = (fraction of the cell's k nearest neighbours that
carry the target) / (the target's overall fraction) − 1`, in PCA space.
k = 300 by default, capped at 10 % of the cells. 0 means chance; above 0 means
enriched. A label-permutation null gives a z-score and FDR per target.

Because the reference frequency counts all cells, a perturbation without effect
can still score above 0 wherever the control cells sit.
`delta_own_vs_control` removes that: it is close to 0 for a target whose cells
sit with the controls.

**Gene programs and perturbation modules** (design of Zhou et al. 2023):

1. For every target with ≥ 20 cells, compute log2 fold changes against controls
   for a panel of response genes: each target's top significant genes, detected
   in ≥ 5 % of cells.
2. Cluster genes on 1 − Pearson correlation into programs (4 by default).
3. Cluster targets on 1 − Spearman correlation into modules (9 by default).

A program's per-cell activity is the mean z-scored expression of its genes.

**Protein effects.** On the CLR scale, the difference of means is a log fold
change. The table also gives Cohen's d, a Mann–Whitney U p-value, a
Benjamini–Hochberg FDR over all target × protein tests, and sign consistency
per sample and per guide.

**RNA–protein concordance.** All three questions use Spearman correlations:

| association | computed over | multiple-testing correction |
|---|---|---|
| PS ↔ protein | the perturbed cells of one target | FDR across all within-target tests; a pooled all-targets value is reported as descriptive only |
| lochNESS ↔ protein | targets | FDR across proteins |
| program ↔ protein | single-guide and control cells, and targets | FDR across program × protein pairs |

Small groups and a small protein panel limit what these can show.

## Configuration

```bash
petrubseq-protein init-config my.yaml             # every key with its default
petrubseq-protein validate --config my.yaml       # check the config and that every input exists
petrubseq-protein audit --config my.yaml          # inspect the inputs before running
petrubseq-protein run --config my.yaml
```

`${DATA_ROOT}` and `${RESULTS_ROOT}` default to `../data` and `../results` next
to the repository. Override them with `PETRUBSEQ_DATA_ROOT` and
`PETRUBSEQ_RESULTS_ROOT`. Unknown keys are rejected.

Every default and the reason for it is in `docs/DEFAULTS.md`; example configs
are in `config/examples/`. Dataset-specific details belong in the YAML, for
example a target label that differs from its gene symbol
(`analysis.perturbation_effects.target_gene_map: {PDL1: CD274}` in the demo).

## Reproducibility and validation

* `pytest` runs 114 tests on synthetic data. They cover:
  * every input format, including checks that each format produces the same
    processed object;
  * guide calling;
  * QC and filtering;
  * the reports;
  * each analysis, including clustering and cluster enrichment, on data with a
    known answer.
* The bundled demo is prepared deterministically; re-preparing it gives
  byte-identical files. The processed demo object passes the contract validator
  (`scripts/validate_processed.py`: 73 checks, 0 failures, 1 expected warning
  because the panel has no isotype controls).
* PS, lochNESS and the gene-program / module calculations reproduce the group's
  reference implementation exactly on shared inputs
  (`scripts/compare_reference_perturbation.py`).
* Every change is compared against a frozen processed object of the Frangieh et
  al. SCP1064 Perturb-CITE-seq dataset (218,331 cells) and its smoke subset,
  slot by slot, and accepted only with zero unexplained numerical differences
  (`docs/REGRESSION.md`).

## Advanced documentation

| document | contents |
|---|---|
| `docs/DEFAULTS.md` | every parameter, its default and the reasoning |
| `docs/OUTPUTS.md` | run directory, tables, figures, report sections |
| `docs/PROCESSED_OBJECT.md` | the processed `.h5ad`, slot by slot |
| `docs/CELL_STATES.md` | Leiden clustering, cluster QC flags and perturbation × cluster enrichment in detail |
| `docs/PERTURBATION_EFFECTS.md` | PS, lochNESS, programs, modules, protein effects, concordance in detail |
| `docs/HPC.md` | running on a SLURM cluster |
| `docs/REGRESSION.md` | the regression protocol |
| `docs/datasets/` | audits of the SCP1064 and Papalexi source data |
| `docs/reference/` | how the analyses relate to the reference implementation |

## Rebuilding the demo

This is optional; the Quick start uses the bundled files. To rebuild the subset
from the public accession:

```bash
python demo/fetch_papalexi_data.py --dest ${DATA_ROOT}/ECCITE-seq      # downloads from GEO; skips files already present
python demo/prepare_papalexi_demo.py --source ${DATA_ROOT}/ECCITE-seq --output /tmp/papalexi_eccite
python demo/check_demo_reproducibility.py /tmp/papalexi_eccite demo/data/papalexi_eccite   # identical to the bundle
```

## Citation

Methods:

* Song B., Liu D., Dai W., et al. Decoding heterogeneous single-cell perturbation
  responses. *Nat Cell Biol* 27, 493–504 (2025). (PS)
* Huang X., Henck J., Qiu C., et al. Single-cell, whole-embryo phenotyping of
  mammalian developmental disorders. *Nature* 623, 772–781 (2023). (lochNESS)
* Traag V.A., Waltman L., van Eck N.J. From Louvain to Leiden: guaranteeing
  well-connected communities. *Sci Rep* 9, 5233 (2019). (Leiden clustering)
* Zhou P., Shi H., Huang H., et al. Single-cell CRISPR screens in vivo map T cell
  fate regulomes in cancer. *Nature* 624, 154–163 (2023). (gene programs and
  perturbation modules)

Data:

* Papalexi E., Mimitou E.P., Butler A.W., et al. Characterizing the molecular
  regulation of inhibitory immune checkpoints with multimodal single-cell
  screens. *Nat Genet* 53, 322–331 (2021). GEO GSE153056. (bundled demo)
* Frangieh C.J., Melms J.C., Thakore P.I., et al. Multimodal pooled
  Perturb-CITE-seq screens in patient models define mechanisms of cancer immune
  evasion. *Nat Genet* 53, 332–341 (2021). Single Cell Portal SCP1064.
  (regression dataset)

Engineering conventions follow weili-lab/perturbseq-pipeline.

## License

MIT; see `LICENSE`. The bundled demo data are derived from GEO GSE153056;
please cite the original publication when using them.
