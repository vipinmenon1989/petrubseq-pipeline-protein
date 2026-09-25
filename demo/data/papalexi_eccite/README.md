# Bundled demo input: Papalexi ECCITE-seq subset (1,800 cells)

A frozen, deterministic subset of a public Perturb-CITE-seq (ECCITE-seq)
experiment, bundled here (~22 MB) so that a fresh clone can run

```bash
petrubseq-protein run --config config/demo_papalexi.yaml
```

without downloading anything. It exists solely as a reproducible software
demonstration / test dataset for this pipeline; it is not a reanalysis of
the publication.

## Source

* Publication: Papalexi E., Mimitou E.P., Butler A.W., Foster S., Bracken B.,
  Mauck W.M. III, Wessels H.-H., Hao Y., Yeung B.Z., Smibert P., Satija R.
  *Characterizing the molecular regulation of inhibitory immune checkpoints
  with multimodal single-cell screens.* Nature Genetics 53, 322–331 (2021).
* Public accession: GEO **GSE153056**
  (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153056).
* Samples used (the **pooled ECCITE-seq screen** in IFN-γ stimulated THP-1
  cells): GSM4633614 (cDNA / RNA counts), GSM4633615 (ADT counts + antibody
  barcodes), GSM4633616 and GSM4633617 (HTO counts and the two HTO barcode
  sets, used only to call hashing singlets), GSM4633618 (guide / GDO counts
  + guide barcodes). Sizes and sha256 of every source file are in
  `demo_manifest.json` (source paths there are written relative to
  `${DATA_ROOT}`; the bundle was prepared on the maintainers' HPC and copied
  verbatim, every other manifest field is the original).
* Original experiment: 20,729 cells (8 lanes, `l1_`–`l8_` barcode prefixes),
  18,649 genes, 4 antibodies, 111 sgRNAs.

## What is bundled (all raw integer UMI counts, values unchanged)

| file | content |
|---|---|
| `filtered_feature_bc_matrix/matrix.mtx.gz` | 18,764 features × 1,800 cells, Matrix Market, integer |
| `filtered_feature_bc_matrix/features.tsv.gz` | id, name, type: 18,649 `Gene Expression` (gene symbols), 4 `Antibody Capture` (CD86, PDL1, PDL2, CD366), 111 `CRISPR Guide Capture` (`<TARGET>g<n>`) |
| `filtered_feature_bc_matrix/barcodes.tsv.gz` | the 1,800 cell IDs (original `l<lane>_<barcode>` identifiers) |
| `cells.csv` | per-cell metadata: lane, lane group, HTO singlet call (top HTO ≥ 10 UMIs and > 3× the runner-up), `hto_sample` (label as deposited), `stimulation` (`IFNg`, constant), sampling stratum |
| `antibodies.csv` | antibody annotation for `protein.feature_table` (antigen and gene symbol; clone not deposited, left empty) |
| `guide_targets.csv` | guide → target (25 gene targets, `NT` non-targeting, `eGFP`) |
| `selected_cells.txt` | the selected cell IDs |
| `demo_manifest.json` | sources (path, size, sha256), parameters, strata, composition, derived dimensions |

No isotype-control antibody exists in the panel. Every hashed sample with
signal in the pooled screen is an IFN-γ (`tx`) sample, so `stimulation` is
constant; the deposited HTO count-matrix labels and the two HTO barcode
files disagree on replicate names, so `hto_sample` keeps the count-matrix
label verbatim (details: `docs/datasets/PAPALEXI_ECCITE_AUDIT.md`).

## Selection (deterministic)

`demo/prepare_papalexi_demo.py`, seed **0**, `--n-cells 1800`,
`--min-per-stratum 8`, prepared 2026-09-25 with the script in this
repository (version 0.2.0 of the pipeline). HTO singlets (17,473 of 20,729)
were stratified by (HTO sample × preliminary guide target from the dominant
rule: top guide ≥ 3 UMIs and > 2× the runner-up, else `ambiguous`), 81
strata, proportional quotas with a floor of 8 (capped by stratum size), the
largest strata trimmed to 1,800; cells drawn without replacement with
`numpy.random.default_rng(0)` from the sorted stratum members. Result: 684 /
571 / 545 cells from the three hashed samples, 199 non-targeting cells, 85
ambiguous cells, 25 gene targets with 21–102 cells each.

## Regenerating the bundle (maintainers / verification)

```bash
python demo/fetch_papalexi_data.py --dest ${DATA_ROOT}/ECCITE-seq          # accession-based GEO FTP, skips complete files
python demo/prepare_papalexi_demo.py --source ${DATA_ROOT}/ECCITE-seq --output /tmp/papalexi_eccite
python demo/check_demo_reproducibility.py /tmp/papalexi_eccite demo/data/papalexi_eccite   # byte-identical files, same manifest
```

Three independent preparations on the source data were byte-identical
(gzip members are written with a fixed timestamp). Normal users do not need
these steps.
