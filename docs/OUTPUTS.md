# Run outputs

Every `petrubseq-protein run` writes one run directory (`output.dir`, usually
`${RESULTS_ROOT}/<dataset>`):

```
results/<run>/
├── report.html                      # self-contained (report.embed_figures: true); primary deliverable
├── report.md                        # same content, figures by relative path (report.write_markdown)
├── processed/
│   ├── <h5ad_name>                  # the processed object (docs/PROCESSED_OBJECT.md)
│   └── <stem>_prefilter.h5ad        # only with output.write_prefilter_h5ad: true
├── figures/
│   ├── qc_rna/{before_filtering,after_filtering}/
│   ├── qc_protein/{before_filtering,after_filtering}/
│   ├── qc_perturbation/
│   ├── representations_rna/  representations_protein/  multimodal/
├── tables/                          # CSV (cell-level tables gzipped)
├── logs/  run.log  resolved_config.yaml  run_manifest.json
└── <run>_results.tar.gz             # only with output.archive: true (matrices excluded)
```

## Regression references

`results/baselines/v0.1/SCP1064` and `results/baselines/v0.1/SCP1064_smoke`
hold the frozen v0.1 objects (read-only); `results/SCP1064` and
`results/SCP1064_smoke` hold the latest validated current outputs and are the
standard destinations of `config/scp1064*.yaml`. `docs/REGRESSION.md` gives
the protocol; a config whose `output.dir` lies under a `baselines` directory
is rejected.

## Run lifecycle (stage banners in `logs/run.log`, timings in the report and manifest)

1 validate configuration · 2 create output tree · 3 logging / resolved config ·
4 input audit · 5 load data (format adapters -> canonical input, validation) · 6 align cells · 7 harmonize perturbations (guide calling when guide counts exist) ·
8 RNA input state / count layer (+ protein representations) · 9 permissive
prefilter · 10 pre-filter QC metrics · 11 BEFORE-filter figures · 12 strict
filtering · 13 filtering audit table · 14 post-filter QC summaries and
perturbation QC · 15 AFTER-filter figures · 16 representations (normalization
of raw-count input happens here, after filtering) · 17 cross-modality
diagnostics · 18 processed h5ad · 19 tables · 20 figure manifest · 21 report ·
22 run manifest · 23 optional archive. With `analysis.perturbation_effects.enabled` a
"perturbation effects" stage runs after the cross-modality diagnostics (24 stages) and the
report gains section 6 *Perturbation effects* (PS, lochNESS, gene programs and modules,
protein effects, RNA–protein concordance, integrated summary). With
`analysis.clustering.enabled` a "cell states and perturbation enrichment" stage follows
(one more stage) and the report gains the section *Cell states and perturbation
enrichment*; section numbers shift so that *Outputs and provenance* stays last.

## Tables

| file | content |
|---|---|
| `qc_filtering_steps.csv` | step, category (input/prefilter/rna/protein/perturbation/final), threshold, cells and genes before/after/removed |
| `qc_summary.csv` | median QC metrics per condition/lane group and overall, after filtering |
| `cell_qc_prefilter.csv.gz` | every cell that entered strict filtering: metrics, flags, perturbation class, `qc_retained`, `removed_by` |
| `cell_qc.csv.gz` | retained cells: design, perturbation and QC columns (+ guide-count diagnostics, `lane_id` / `barcode_original` when applicable) |
| `protein_qc.csv` | per antibody: role, isotype control, count/normalized statistics, `frac_cells_above_isotype`, `background_dominated` |
| `protein_features.csv`, `protein_pca_loadings.csv` | antibody metadata (role, isotype control, presence, `feature_id`, `antibody_name`, `protein_name`, `gene_symbol`, `clone`, `feature_type`, `isotype`, `annotation_source`; empty = not available); protein PC loadings |
| `guide_features.csv` | only with guide counts: per guide `target`, `control_class`, `total_umis`, `n_cells_detected`, `n_cells_dominant` (+ 10x feature columns) |
| `perturbation_qc.csv` | scalar perturbation QC (class counts, coverage, fractions, assignment source, guide-calling rule and provided-vs-counts agreement when applicable) |
| `perturbation_counts.csv` | per perturbation label: class, cells, per condition |
| `guide_coverage.csv`, `target_coverage.csv` | per guide / target: class, single-guide and any-guide cells, per condition, `low_coverage` flag |
| `condition_guide_coverage.csv`, `condition_target_coverage.csv` | long (guide or target, condition) tables with `low_coverage` |
| `moi_summary.csv`, `condition_summary.csv` (+ `sample_summary.csv` ...) | guides-per-cell by condition; per-group cell counts, class breakdown, median QC |
| `rna_gene_qc.csv` | per gene: mt/ribo flags, detection, HVG columns |
| `figure_manifest.csv` | section, stage, name, title, caption, in_report, path for every figure |
| `run_summary.csv` | headline numbers of the run |
| `cell_states/*.csv` | Stage F only: `cluster_summary`, `cluster_composition_by_class`, `cluster_composition_by_target`, `cluster_composition_by_<sample/lane/batch/...>`, `perturbation_cluster_enrichment`, `cluster_enrichment_vs_lochness`, `cluster_enrichment_skipped` (columns: docs/CELL_STATES.md) |
| `perturbation_effects/*.csv` | Stage E only: `ps_targets`, `ps_skipped`, `ps_signatures`, `lochness_targets`, `lochness_skipped`, `lochness_by_sample`, `perturbation_effect_matrix`, `perturbation_de_mask`, `gene_programs`, `perturbation_modules`, `module_program_strength`, `perturbation_program_effects`, `protein_effects`, `protein_effect_matrix`, `ps_protein_association`, `lochness_protein_association`, `lochness_protein_summary`, `program_protein_association_cells`, `program_protein_association_targets`, `perturbation_summary` (columns: docs/PERTURBATION_EFFECTS.md) |

## Figures

Before/after pairs (same drawing code, `stage` differs; thresholds drawn dashed):

| section | figure | when |
|---|---|---|
| qc_rna | `rna_qc_distributions` (genes, log counts, % mito, % ribo), `rna_counts_vs_genes`, `rna_cells_per_<group>`, `rna_<metric>_by_<group>` | always (group figures when a 2-30-level group exists) |
| qc_protein | `protein_qc_distributions` (total ADT, antibodies detected, % isotype), `protein_targeting_vs_isotype` (extreme cells marked), `protein_depth_vs_rna_depth`, `protein_antibody_distributions`, `protein_antibody_summary` | raw ADT counts available |
| qc_protein | `protein_normalized_distributions` | any normalized protein matrix (the only protein figure for normalized-only inputs) |
| qc_perturbation | `perturbation_class_composition`, `perturbation_class_by_condition`, `guides_per_cell`, `cells_per_perturbation`, `coverage_histograms`, `condition_target_coverage`; `guide_count_diagnostics` only when guide counts exist | post-filter |
| representations_rna | `rna_pca_variance`, `rna_pc_vs_depth`, `rna_umap_<color>`, `provided_embedding_<color>` | representation |
| representations_protein | `protein_pca_variance`, `protein_pca_loadings`, `protein_pc_vs_depth`, `protein_umap_<color>` | representation |
| multimodal | `knn_overlap`, `multimodal_umap_<color>` (when enabled) | representation |
| cell_states | `umap_leiden`, `cluster_composition` (stage `clusters`); `perturbation_cluster_enrichment` (`enrichment`) | only with `analysis.clustering.enabled` |
| perturbation_effects | `ps_target_median`, `ps_distributions` (stage `ps`); `lochness_targets`, `lochness_umap` (`lochness`); `perturbation_program_heatmap`, `effect_matrix_programs`, `module_program_strength` (`gene_programs`); `protein_effect_heatmap` (`protein_effects`); `ps_vs_protein`, `program_activity_vs_protein`, `rna_vs_protein_effects`, `perturbation_overview` (`concordance`) | only with `analysis.perturbation_effects.enabled` |

## report.html

Header/navigation (run, time, version) · Run summary cards · Warnings ·
1 Inputs and modality audit · 2 Cell and RNA QC (filtering steps, QC summary,
before, after) · 3 Protein QC (representations, features, isotype map,
background-dominated antibodies, extreme cells, before, after) · 4 Perturbation
QC · 5 Representations (RNA, protein, cross-modality) · 6 Outputs and
provenance (deliverables, object schema, timings, table and figure manifests,
resolved config, versions, provenance).
