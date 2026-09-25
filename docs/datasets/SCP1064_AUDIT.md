# SCP1064 dataset audit

Local copy audited: `<project>/SCP1064` (since moved, layout unchanged, to `<project>/data/SCP1064/`, i.e. `${DATA_ROOT}/SCP1064`)
(Single Cell Portal download of Frangieh et al. 2021, *Multimodal pooled
Perturb-CITE-seq screens in patient models define mechanisms of cancer immune
evasion*, Nat Genet 53:332). Audited 2026-09-24 with shell streaming tools
(`zcat`, `cut`, `awk`, `sort`, `comm`) and `petrubseq-protein audit --config
config/scp1064.yaml` (which samples 2,000 RNA cells and loads both protein
matrices fully). **The local files are the source of truth**; nothing below is
taken from the portal description.

Reproduce the machine-readable part with:

```bash
petrubseq-protein audit --config config/scp1064.yaml --out ../results/SCP1064/audit/scp1064_audit.json
```

## 1. File inventory

Total on disk: **7.7 GB** (16 files). All matrices are plain-text CSV; the
large ones are gzip-compressed. All cell-bearing files use the same identifier
scheme `CELL_1 … CELL_218331` (no 10x barcodes, no `-1` suffix, no sample
prefix).

| path | type / role | gzip | size | rows x cols (data) | orientation | status |
|---|---|---|---|---|---|---|
| `other/RNA_expression.csv.gz` | RNA expression, **all cells** | yes | 3.5 GB (≈37 GB text) | 23,712 genes x 218,331 cells | genes x cells, first column `GENE` | complete; ln(TPM+1) |
| `expression/RNA_expression_subset1a.csv.gz` | RNA chunk (cells 1–27,291) | yes | 504 MB | 23,712 x 27,291 | genes x cells | complete |
| `expression/RNA_expression_subset1b_alpha.csv.gz` | RNA chunk (27,292–40,936) | yes | 257 MB | 23,712 x 13,645 | genes x cells | complete |
| `expression/RNA_expression_subset1b_beta.csv.gz` | RNA chunk (40,937–54,582) | yes | 248 MB | 23,712 x 13,646 | genes x cells | complete |
| `expression/RNA_expression_subset2.csv.gz` | RNA chunk (54,583–109,164) | yes | 850 MB | 23,712 x 54,582 | genes x cells | complete |
| `expression/RNA_expression_subset3.csv.gz` | RNA chunk (109,165–163,746) | yes | 874 MB | 23,712 x 54,582 | genes x cells | complete |
| `expression/RNA_expression_subset4.csv.gz` | RNA chunk (163,747–218,331) | yes | 939 MB | 23,712 x 54,585 | genes x cells | complete |
| `documentation/Protein_expression.csv.gz` | **normalized** protein (ADT) | yes | 5.5 MB (37 MB text) | 20 antibodies x 218,331 cells | proteins x cells, first column `Protein` | complete; derived |
| `other/raw_CITE_expression.csv.gz` | **raw** ADT UMI counts | yes | 3.7 MB (15 MB text) | 24 antibodies x 218,328 cells | antibodies x cells, first column unnamed | 3 cells missing; raw |
| `metadata/RNA_metadata.csv` | cell metadata (SCP format, `TYPE` row) | no | 12 MB | 218,331 x 5 | cells x columns, `NAME` | complete |
| `documentation/all_sgRNA_assignments.txt` | per-cell guide list | no | 5.4 MB | 218,331 x 1 (`Cell,sgRNAs`) | cells x columns | complete; derived (upstream guide calling) |
| `cluster/5fd0e449771a5b0db7207711/RNA_UMAP_cluster.csv` | authors' RNA UMAP (SCP cluster file) | no | 7.0 MB | 218,331 x 2 (`X`,`Y`) | cells x columns, `NAME` | complete; derived |
| `documentation/untreated_regulatory_matrix.csv` | perturbation x gene effect matrix (Control) | no | 9.6 MB | 4,481 x ~1,000 | genes x perturbations | derived result; not used |
| `documentation/treated_regulatory_matrix.csv` | same, IFNγ | no | 10.6 MB | ~4,480 x ~1,000 | genes x perturbations | derived result; not used |
| `documentation/cocx_regulatory_matrix.csv` | same, co-culture | no | 10.4 MB | 4,481 x ~1,000 | genes x perturbations | derived result; not used |
| `documentation/Fig4C_pvals_heatmap_values.csv` | figure 4C values | no | 5.1 MB | 2,195 x 176 | table | derived result; not used |
| `file_supplemental_info.tsv` | SCP file manifest | no | 1 KB | 16 rows | table | metadata only (no units / raw flags filled in) |

Verified facts about the RNA files:

* The six `expression/` chunks partition the cells of `other/RNA_expression.csv.gz`
  exactly (same order, no overlap, no gap; values byte-identical for the rows
  checked). The pipeline reads the six chunks in parallel; the combined file is
  redundant.
* Gene axis: 23,712 unique symbols (GRCh38 v3.0.0 annotation; e.g. `A1BG`,
  `AL589666.1`), no duplicates; **13 `MT-` genes and 102 `RPS/RPL` genes are
  present**, so mitochondrial and ribosomal fractions can be computed.
* Cells are numbered in condition blocks: `CELL_1–57627` = Control,
  `CELL_57628–145217` = IFNγ, `CELL_145218–218331` = Co-culture.

## 2. Modalities present locally

| modality | present | file | notes |
|---|---|---|---|
| RNA raw counts | **no** (but exactly reconstructable, see §5) | – | |
| RNA normalized expression | yes | `other/RNA_expression.csv.gz` / `expression/*` | ln(TPM+1), i.e. `log1p(1e6 * count / UMI_count)` |
| CITE-seq/ADT raw counts | yes | `other/raw_CITE_expression.csv.gz` | 24 antibodies incl. 4 isotype controls |
| CITE-seq/ADT normalized | yes | `documentation/Protein_expression.csv.gz` | 20 antibodies, isotype-ratio transform (reproduced exactly) |
| guide count matrix | **no** | – | guide calling done upstream by the authors |
| guide assignments | yes | `documentation/all_sgRNA_assignments.txt` | all guides per cell |
| guide -> target mapping | implicit | guide IDs `GENE_<n>` | 3 guides per target (controls: many) |
| perturbation labels | implicit | derived from single-guide cells | |
| cell metadata | yes | `metadata/RNA_metadata.csv` | 5 columns |
| sample / donor / patient metadata | **no** | – | single patient-derived model; no sample column |
| experimental condition | yes | `condition` | Control / IFNγ / Co-culture |
| batch / library / lane | **no** | – | `library_preparation_protocol` is constant; the raw-ADT column order shows 26 contiguous cell-ID blocks that probably reflect the original 10x libraries, but no label is provided |
| HTO / hashing | no | – | |
| QC metrics | partial | `UMI_count` in metadata | genes detected, mito fraction not provided (computed by the pipeline) |
| precomputed embedding | yes | `RNA_UMAP_cluster.csv` | authors' UMAP; no cluster labels |
| other tables | yes | 3 regulatory matrices, Fig4C values | downstream results, not inputs |

## 3. Cell identifiers

| set | total | unique | duplicates |
|---|---|---|---|
| RNA (6 chunks = combined) | 218,331 | 218,331 | 0 |
| protein normalized | 218,331 | 218,331 | 0 |
| protein raw counts | 218,328 | 218,328 | 0 |
| metadata | 218,331 | 218,331 | 0 |
| guide assignments | 218,331 | 218,331 | 0 |
| authors' UMAP | 218,331 | 218,331 | 0 |

* RNA ∩ protein-normalized ∩ metadata ∩ guides ∩ UMAP = **218,331 (100 %)**.
* RNA ∩ protein-raw = **218,328 (99.999 %)**. Missing from the raw ADT file:
  `CELL_27337`, `CELL_34511`, `CELL_56031` (all Control block). Those three
  cells still have normalized protein values, so they are kept with
  `has_protein_counts = False` and NaN raw counts.
* No suffix/prefix differences: identifiers match exactly across files; no
  transformation is applied. Order differs only in the raw ADT file
  (26 blocks, e.g. `7700–14735, 50005–56030, …`); the pipeline aligns by ID, not
  position.

## 4. Experimental design (from `metadata/RNA_metadata.csv`)

| column | type | values |
|---|---|---|
| `NAME` | cell ID | `CELL_1` … `CELL_218331` |
| `library_preparation_protocol` | constant | `10X 3' v3 sequencing` (218,331) |
| `condition` | group | `Control` 57,627 · `IFNγ` 87,590 (note: non-ASCII γ, UTF-8) · `Co-culture` 73,114 |
| `MOI` | numeric | number of assigned guides, 0–19 (0: 23,028 · 1: 126,966 · 2: 45,135 · ≥3: 23,202) |
| `sgRNA` | group | guide for **MOI = 1 cells only** (126,966 filled, 91,365 empty); 818 unique guides |
| `UMI_count` | numeric | per-cell RNA UMIs: min 500, median 10,988, mean 12,700, max 59,486 |

Design summary: one patient-derived melanoma model, three conditions, no
donor/replicate/batch/lane columns. Conditions are preserved as `obs['condition']`
and never pooled by the pipeline.

## 5. RNA state

* Values are floats, min non-zero 2.90, max 12.17 in the sample; every cell's
  `sum(expm1(x))` equals **1,000,000** (CV 6e-8), so the matrix is
  **ln(TPM+1)** with `TPM = 1e6 * count / UMI_count`.
* The smallest non-zero value per cell corresponds to exactly one UMI, and
  `round(expm1(x) * UMI_count / 1e6)` is integer to within 7e-4 for every
  entry, with row sums equal to `UMI_count`. **Integer counts are therefore
  reconstructed exactly** into `layers['reconstructed_counts']` (never
  `layers['counts']`, which is reserved for directly observed UMI counts;
  `uns['petrubseq_protein']['rna']['counts_source']` records the formula).
  The reconstruction is unique (each integer count maps to a distinct value)
  but it recovers only the 23,712 genes retained by the authors, so it is not
  the original UMI matrix.
* QC already applied by the authors (paper methods): cells with < 200 genes or
  > 18 % mitochondrial reads removed, genes in < 200 cells removed. Sampled
  cells: 448–7,311 genes detected (median 3,771); zero fraction 0.836.
* The pipeline **does not normalize the RNA again**; `X` is the provided ln(TPM+1).

## 6. Protein / ADT

### Panel (24 antibodies in the raw file; 20 targets in the normalized file)

| antibody | role | isotype control used by the authors | raw median / mean | zero frac (raw) | normalized zero frac | background-dominated¹ |
|---|---|---|---|---|---|---|
| CD117 | target | Mouse_IgG1 | 1 / 8.97 | 0.50 | 0.58 | yes |
| CD119 (IFNGR1) | target | Mouse_IgG1 | 5 / 8.36 | 0.04 | 0.07 | no |
| CD140a | target | Mouse_IgG1 | 0 / 3.16 | 0.71 | 0.80 | yes |
| CD140b | target | Mouse_IgG1 | 0 / 1.48 | 0.57 | 0.66 | yes |
| CD172a | target | Mouse_IgG2a | 1 / 3.97 | 0.27 | 0.35 | no |
| CD184 (CXCR4) | target | Mouse_IgG2a | 1 / 31.4 | 0.39 | 0.46 | no |
| CD202b | target | Mouse_IgG1 | 0 / 0.58 | 0.67 | 0.77 | yes |
| CD274 (PD-L1) | target | Mouse_IgG2b | 11 / 21.2 | 0.04 | 0.05 | no |
| CD29 | target | Mouse_IgG1 | 120 / 173 | 0.00 | 0.00 | no |
| CD309 | target | Mouse_IgG1 | 0 / 1.15 | 0.63 | 0.72 | yes |
| CD44 | target | Mouse_IgG1 | 523 / 674 | 0.00 | 0.00 | no |
| CD47 | target | Mouse_IgG1 | 32 / 40.6 | 0.00 | 0.00 | no |
| CD49f | target | Rat_IgG2a | 15 / 28.2 | 0.01 | 0.04 | no |
| CD58 | target | Mouse_IgG1 | 5 / 7.67 | 0.04 | 0.07 | no |
| CD59 | target | Mouse_IgG2a | 7 / 11.7 | 0.03 | 0.03 | no |
| CD61 | target | Mouse_IgG1 | 1 / 1.25 | 0.47 | 0.57 | yes |
| HLA_A | target | Mouse_IgG2a | 25 / 34.3 | 0.00 | 0.00 | no |
| HLA_E | target | Mouse_IgG1 | 0 / 7.29 | 0.53 | 0.63 | yes |
| CD9 | target | Mouse_IgG1 | 38 / 55.2 | 0.00 | 0.00 | no |
| CD279 (PD-1) | target | Mouse_IgG1 | 0 / 1.88 | 0.51 | 0.61 | yes |
| Rat_IgG2a | isotype control | – | 0 / 1.13 | 0.52 | – | – |
| Mouse_IgG1 | isotype control | – | 0 / 0.64 | 0.73 | – | – |
| Mouse_IgG2a | isotype control | – | 0 / 0.45 | 0.72 | – | – |
| Mouse_IgG2b | isotype control | – | 0 / 0.37 | 0.77 | – | – |

¹ "background-dominated" = fewer than 50 % of cells have more counts than their
matched isotype control (equivalently, normalized value is 0 in > 50 % of
cells). Eight antibodies (CD117, CD140a, CD140b, CD202b, CD309, CD61, HLA_E,
CD279) fall in this class and should be interpreted with care downstream.

Names: the normalized file uses `"<name> protein"` (e.g. `CD117 protein`,
`HLA_A protein`); the raw file uses bare names. The pipeline strips the
`" protein"` suffix so both matrices share canonical names.

### Values

* Raw file: integers, min 1 (non-zero), max 19,718 in the sample and up to
  878,564 total ADT counts in extreme cells; per-cell totals median 829
  (5–95 %: 292–2,570); proteins detected per cell median 16/24. A handful of
  cells carry six-figure counts of a single antibody (e.g. CD184, HLA_E) —
  likely antibody aggregates; they are flagged by the MAD outlier columns, not
  removed.
* Normalized file: floats in [0, 12.2]; row sums are not constant (not
  size-normalized, not CLR). Checking all 20 antibodies on 20,000 cells, the
  values equal **`max(0, ln((count + 1) / (isotype_count + 1)))`** with the
  isotype map in the table above, with **maximum absolute difference 0.0**.
  This is exactly `protein_expression.py` of the original repository. The
  provided values are therefore kept unchanged as `obsm['protein']`, and the
  pipeline additionally stores a CLR transform of the raw counts as
  `obsm['protein_clr']` (computed from raw counts, never from the provided
  values, so nothing is normalized twice).

## 7. Perturbations

Guide calling was done upstream; **no guide-count matrix is included**, so the
pipeline does not (and must not) call guides for SCP1064. It only harmonizes
the provided assignments. All numbers from `all_sgRNA_assignments.txt`, which
agrees with the metadata `sgRNA` column on every MOI = 1 cell.

| quantity | value |
|---|---|
| cells | 218,331 |
| unassigned (0 guides) | 23,028 (10.5 %) |
| single-guide cells | 126,966 (58.2 %) |
| multi-guide cells (2–19 guides) | 68,337 (31.3 %) |
| unique guides | 818 |
| unique targets (guide name without `_<n>`) | 250 = 248 genes + 2 control classes |
| guides per gene target | 3 (all genes) |
| control guides | `NO_SITE_<n>` (non-targeting; 37 guides) and `ONE_NON-GENE_SITE_<n>` (intergenic cutting controls; ~30 guides) |
| single-guide cells by class | targeting 111,345 · non-targeting 8,955 · intergenic 6,666 |
| cells per target (single-guide, genes) | median ≈ 509, min 12, max 829 (IFNGR2); < 10: 0, < 20: 2, < 30: 7, < 50: 13 |
| cells per guide (single-guide) | median 170, min 1, max 358; < 10: 34, < 20: 61, < 30: 77, < 50: 102 |
| min cells per target across conditions | median 121 (min 2) |

Single-guide class x condition:

| class | Control | IFNγ | Co-culture |
|---|---|---|---|
| targeting | 26,716 | 43,909 | 40,720 |
| non-targeting (NO_SITE) | 2,189 | 3,473 | 3,293 |
| intergenic (ONE_NON-GENE_SITE) | 1,581 | 2,671 | 2,414 |
| multi/unassigned | 27,141 | 37,537 | 26,687 |

Pipeline representation (`obs`): `guides` (all, `;`-joined), `n_guides`,
`guide`, `target`, `perturbation` (target / `multi` / `unassigned`),
`control_class` (`non_targeting` / `intergenic` / ''), `perturbation_class`
(`single_targeting` / `single_control` / `multi_targeting` / `multi_control` /
`mixed_control_targeting` / `unassigned`), `targets`, `n_targets`,
`is_control`, `is_targeting`, `is_single_guide`. Multi-guide cells are kept and
labelled, not dropped, so later analyses can decide how to use them (the
original paper fits all guides jointly).

## 8. Data-quality issues and decisions

1. **No raw RNA counts shipped** → reconstructed exactly (validated per run;
   the run report records the maximum deviation).
2. **Three cells lack raw ADT counts** → kept, flagged (`has_protein_counts`).
3. **No sample / donor / batch / lane labels** → `sample`, `donor`, `batch`
   are `null` in the config; the raw-ADT block structure is documented above
   but not used.
4. **Non-ASCII condition label** (`IFNγ`) → preserved verbatim (UTF-8) in
   `obs` and YAML.
5. **Eight background-dominated antibodies** and rare extreme ADT totals →
   flagged in `protein_qc_summary.tsv` / outlier columns, not filtered.
6. **31 % multi-guide cells, 10.5 % unassigned** → kept and labelled.
7. **Guides with very few cells** (34 guides < 10 cells) → reported in
   `guide_counts.tsv`; no minimum is enforced in v0.1.
8. Authors' QC (≥ 200 genes, ≤ 18 % mito) is already applied; the pipeline's
   defaults use the same thresholds for flagging and `qc.filter: false`.

## 9. Resources: estimate and measured full run

| item | estimate (audit) | measured (SLURM job 20838373, 2026-09-24) |
|---|---|---|
| non-zero RNA values | ≈ 0.85–0.9 x 10⁹ | ≈ 0.75 x 10⁹ (median 3,310 genes/cell) |
| peak RAM | ≈ 35–45 GB | **23.9 GB** (`/usr/bin/time`; sacct MaxRSS 22.0 GB) |
| CPUs | 8 (6 parser processes) | 8 allocated, 184 % average utilisation (52 min CPU) |
| wall time | 1.5–3 h | **28 min** (load 6.0 min, count reconstruction + QC 1.8 min, representations 15.2 min, h5ad write 4.2 min) |
| output | ≈ 4–6 GB | **3.6 GB** (h5ad 3.36 GB gzip; tables 23 MB; figures 11 MB) |

`slurm/run_scp1064.slurm` requests 8 CPUs, 48 GB and 2 h (≈ 2x the measured
peak and 4x the measured time). The login node (4 CPUs, 15 GB RAM) is not
suitable for the run or for loading the processed object; the audit itself
(headers + protein files + 2,000 sampled RNA cells) needs 1 GB and ~2 min.
