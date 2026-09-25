# Papalexi et al. 2021 ECCITE-seq (GSE153056) — local data audit

Audited 2026-09-25 from the files under `${DATA_ROOT}/ECCITE-seq/`
(`<project>/data/ECCITE-seq/`, immutable
source data; nothing was modified). Only headers, first columns, the small
matrices (ADT / HTO / GDO) and the tiny barcode tables were read on the login
node; the full cDNA matrix is only touched by the SLURM preparation job.

Citation: Papalexi E., Mimitou E.P., Butler A.W., et al. *Characterizing the
molecular regulation of inhibitory immune checkpoints with multimodal
single-cell screens.* Nat Genet 53, 322–331 (2021). GEO accession
**GSE153056** (samples GSM4633605–GSM4633618).

## 1. Inventory (23 files, 256 MB)

`Papalexi.tar` (128 MB) contains exactly the 22 GEO supplementary files listed
below (same names and sizes); it is redundant and not used.

| GSM | file | size | experiment | modality |
|---|---|---|---|---|
| GSM4633605 | `GSM4633605_CITE_cDNA_counts.tsv.gz` | 22.8 MB | CITE-seq (no guides) | RNA counts |
| GSM4633606 | `GSM4633606_CITE_ADT_{Barcodes.csv,counts.tsv}.gz` | 170 B / 104 KB | CITE-seq | 9 ADTs (PDL1, CD86, IDO1, TIM3, PDL2, HVEM, CD47, BTLA, CD66a) |
| GSM4633607 | `GSM4633607_CITE_HTO_{Barcodes.csv,counts.tsv}.gz` | 123 B / 73 KB | CITE-seq | 6 HTOs (HTO22–28) |
| GSM4633608 | `GSM4633608_ECCITE_Arrayed_cDNA_counts.tsv.gz` | 39.9 MB | ECCITE arrayed | RNA counts |
| GSM4633609 | `GSM4633609_ECCITE_Arrayed_ADT_{Barcodes.csv,counts.tsv}.gz` | 124 B / 151 KB | ECCITE arrayed | 4 ADTs |
| GSM4633610/11 | `GSM4633610_ECCITE_Arrayed_HTO_counts.tsv.gz`, `..._HTO1_Barcodes`, `GSM4633611_..._HTO2_Barcodes` | 187 KB / 187 B / 194 B | ECCITE arrayed | HTOs named by arrayed guide (HTO_ETV7g1, HTO_NT5, ...) |
| GSM4633612/13 | `GSM4633612_ECCITE_Arrayed_GDO_Barcodes.csv.gz`, `..._GO_CITE03_counts.tsv.gz`, `GSM4633613_..._GO_LENTI_counts.tsv.gz` | 270 B / 170 KB / 105 KB | ECCITE arrayed | 16 guide barcodes, two guide-capture libraries |
| GSM4633614 | `GSM4633614_ECCITE_cDNA_counts.tsv.gz` | 64.1 MB | **ECCITE pooled screen** | RNA counts |
| GSM4633615 | `GSM4633615_ECCITE_ADT_{Barcodes.csv,counts.tsv}.gz` | 116 B / 196 KB | ECCITE pooled | 4 ADTs |
| GSM4633616/17 | `GSM4633616_ECCITE_HTO_counts.tsv.gz`, `GSM4633616_ECCITE_HTO_1-4_Barcodes.csv.gz`, `GSM4633617_ECCITE_HTO_5-8_Barcodes.csv.gz` | 170 KB / 164 B / 160 B | ECCITE pooled | 12 HTO rows; two 8-HTO barcode sets |
| GSM4633618 | `GSM4633618_ECCITE_GDO_{Barcodes.csv,counts.tsv}.gz` | 1.2 KB / 338 KB | ECCITE pooled | 112 guide barcodes / 111 guide count rows |

Three experiments are present: the CITE-seq control experiment (no CRISPR),
the arrayed ECCITE experiment (guides delivered per well, identified by
HTO) and the **pooled ECCITE-seq screen**. The pooled screen is the only one
with RNA + ADT + a per-cell guide-count matrix, so it is the Stage D demo.

## 2. Pooled ECCITE-seq screen (GSM4633614–GSM4633618)

| SOURCE FILE | MODALITY | FORMAT | DIMENSIONS | STATE | CELL IDENTIFIER | FEATURE IDENTIFIER | ROLE | STAGE-C ADAPTER | NOTES |
|---|---|---|---|---|---|---|---|---|---|
| `GSM4633614_ECCITE_cDNA_counts.tsv.gz` | RNA | dense TSV, features x cells, quoted header, first column = gene symbol | 18,649 genes x 20,729 cells | raw UMI counts (integers) | `l<lane>_<16nt barcode>` (lanes l1–l8, globally unique, no `-1` suffix) | gene symbol (unique, 13 `MT-` genes, 101 `RPS/RPL`) | `inputs.rna` | `dense_csv` (or, after demo preparation, `multiplexed` MTX) | 64 MB gzip; ~1.5 GB dense; read only on SLURM |
| `GSM4633615_ECCITE_ADT_counts.tsv.gz` | protein (ADT) | dense TSV, features x cells | 4 x 20,729 | raw UMI counts | same as RNA (identical 20,729 barcodes) | `CD86`, `PDL1`, `PDL2`, `CD366` | `inputs.protein_counts` | `dense_csv` / MTX `Antibody Capture` | **no isotype controls** in the panel; median 342 ADT UMIs / cell |
| `GSM4633615_ECCITE_ADT_Barcodes.csv.gz` | protein feature table | 2-column CSV (antibody barcode, name), no header | 4 rows | – | – | antibody name | feature metadata | `protein.feature_table` (after preparation) | gives only the name; gene symbols added in the demo table (`CD274`, `PDCD1LG2`, `HAVCR2`) |
| `GSM4633618_ECCITE_GDO_counts.tsv.gz` | guide (sgRNA) counts | dense TSV, features x cells | 111 x 20,729 | raw UMI counts | same 20,729 barcodes | guide ID `<TARGET>g<n>` (e.g. `PDL1g1`, `NTg5`, `eGFPg1`) | `inputs.guide_counts` | `dense_csv` / MTX `CRISPR Guide Capture` | median 217 guide UMIs / cell; dominant rule (>= 3 UMIs, > 2x runner-up): 19,681 assigned, 1,048 ambiguous, 0 unassigned |
| `GSM4633618_ECCITE_GDO_Barcodes.csv.gz` | guide library | 2-column CSV (protospacer, guide ID) | 112 rows | – | – | guide ID | guide -> target mapping | `perturbation.guide_target_regex` / `guide_target_table` | `NTg6` is in the barcode file but absent from the count matrix; 25 gene targets x 3–4 guides, 9 `NT` non-targeting guides (`NTg1-5,7-10`), `eGFPg1` (transgene control, 0 dominant cells) |
| `GSM4633616_ECCITE_HTO_counts.tsv.gz` | cell hashing (HTO) | dense TSV, features x cells | 12 x 20,729 | raw UMI counts | same 20,729 barcodes | `rep{1-4}-{tx,ctrl}`, `PDL1g{1,2}-{tx,ctrl}` | sample demultiplexing | none (HTO demultiplexing is done by the demo preparation script, not the pipeline) | see section 3 |
| `GSM4633616_ECCITE_HTO_1-4_Barcodes.csv.gz`, `GSM4633617_ECCITE_HTO_5-8_Barcodes.csv.gz` | HTO barcode sets | 2-column CSV | 8 rows each | – | – | HTO sequence -> label | lane-group specific HTO naming | – | the **same 8 sequences** carry different labels in the two files (`AGGACCATCCAA` = `rep3_tx` in lanes 1–4, `rep1_tx` in lanes 5–8) |

Relationships: all four count matrices share the identical, ordered set of
20,729 cell barcodes (pairwise overlap 20,729 / 20,729), so alignment is exact
and no barcode reconciliation is needed. The `l<k>_` prefix encodes the 10x
lane (l3 2,390, l4 2,582, l6 3,049, l7 2,976, l8 3,098 cells; the full per-lane
table is in `demo_manifest.json`); it is already globally unique, so the
Stage C multi-lane suffixing is *not* needed (the lane is recorded as
metadata instead). Every matrix is dense, features x cells, integer.

## 3. Sample / condition metadata (derived, with a caveat)

There is no per-cell metadata file. Sample identity must come from the HTO
matrix. Per lane, the HTO counts are:

| lane group | HTO rows with real signal (as labelled in `GSM4633616_..._HTO_counts`) | other rows |
|---|---|---|
| l1–l4 (median 52 HTO UMIs / cell) | `rep1-tx` (95 % of cells by argmax); `rep2-tx` at ~9 UMIs / cell, i.e. ambient | `rep1-ctrl`, `rep2-ctrl`, `PDL1g*` at <= 1 % |
| l5–l8 (median 135 HTO UMIs / cell) | `rep3-tx` and `rep4-tx` (two pooled samples, ~45 % / 40 %) | `rep3-ctrl`, `rep4-ctrl`, `PDL1g*` at ~1 % |

Caveats recorded for the demo:

* The count-matrix row labels contradict the per-lane barcode files: lanes
  1–4 signal sits in rows named `rep1-*` although the lanes 1–4 barcode file
  names those sequences `rep3_*`, and vice versa. **Replicate identity
  (rep1 vs rep3 …) is therefore ambiguous**; the demo keeps the deposited
  label verbatim as `hto_sample` and does not rename it.
* Under either naming, every HTO with real signal is a `tx` (IFN-γ
  stimulated) sample; the `ctrl` HTOs are at background level. **All cells of
  the pooled screen are IFN-γ stimulated**; `stimulation` is constant
  (`IFNg`). The screen's contrast is perturbed vs non-targeting cells within
  the stimulated population, not stimulated vs control.
* HTO demultiplexing in the demo preparation uses a simple deterministic
  rule (top HTO >= 10 UMIs and > 3 x the runner-up = singlet; otherwise
  `doublet_or_negative`). Only singlets enter the demo subset; the rule and
  counts are written to the manifest.

## 4. Decisions for the demo (Stage D)

* Experiment: pooled ECCITE-seq screen (section 2), one coherent
  experiment; the CITE-seq and arrayed experiments are not combined with it.
* Layout written by `demo/prepare_papalexi_demo.py`: one 10x-style
  `filtered_feature_bc_matrix/` (Gene Expression + Antibody Capture + CRISPR
  Guide Capture) read by the Stage C `inputs.multiplexed` MTX adapter,
  `cells.csv` (lane, HTO call, sample, stimulation), `antibodies.csv`
  (`protein.feature_table`), `guide_targets.csv`
  (`perturbation.guide_target_table`), `demo_manifest.json`.
* Guide -> target: `^(?P<target>.+?)g\d+$` (and the explicit table);
  controls: `non_targeting` = `^NT$`, `transgene` = `^eGFP$`.
* Protein: raw ADT counts -> CLR (`protein.normalization: auto`); no
  isotype controls exist, so isotype-based flags are inactive (documented
  warning); 4 antibodies -> `protein.n_pcs: 3`, `umap.min_features: 3`.
* Guide assignment: only guide counts exist -> `perturbation.assignment.source:
  auto` resolves to `guide_counts`, dominant rule (reference defaults).
* Genericity gate (D1): every change above is **A. configuration** or
  **B. demo preparation**. No generic pipeline deficiency (C) was found for
  this dataset; nothing dataset-specific (D) was added to the package.

## 5. Derived demo (`${DATA_ROOT}/demo/papalexi_eccite/`, SLURM job 20846967)

| item | value |
|---|---|
| source cells / genes / antibodies / guides | 20,729 / 18,649 / 4 / 111; cells per lane l1–l8: 2,185 / 2,271 / 2,390 / 2,582 / 2,178 / 3,049 / 2,976 / 3,098 |
| HTO singlets (top >= 10 UMIs, > 3x runner-up) | 17,473 (3,256 doublet/negative); samples `rep1-tx` 6,630, `rep3-tx` 5,576, `rep4-tx` 5,267 |
| preliminary guide call (dominant rule) | 19,681 assigned, 1,048 ambiguous, 0 unassigned |
| selection | seed 0, 1,800 cells, 81 strata (HTO sample x target incl. `NT` and `ambiguous`), floor 8 per stratum (capped by stratum size; smallest stratum 5) |
| selected composition | `rep1-tx` 684 / `rep3-tx` 571 / `rep4-tx` 545; `NT` 199, `ambiguous` 85, 25 gene targets 21–102 cells each |
| derived matrix | 18,764 features x 1,800 cells (18,649 GEX + 4 ADT + 111 guides), 6,347,995 non-zeros, integer |
| reproducibility | three independent preparations byte-identical (`demo/check_demo_reproducibility.py`) |
| source integrity | sizes + mtimes unchanged before/after; sha256 of every source file in `demo_manifest.json` |
