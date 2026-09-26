"""Configuration schema for petrubseq-pipeline-protein.

One YAML file fully describes a run. Dataset-specific facts (paths, column
names, control naming) live in the YAML, never in Python. Every key has a
documented default below; unknown keys raise an error so typos are caught
before any data is read.

The layout follows the config-driven philosophy of weili-lab/perturbseq-pipeline
(dataclass per section, ``Config.from_yaml`` / ``Config.to_dict``), but the
schema is deliberately much smaller.
"""

from __future__ import annotations

import dataclasses
import os
from dataclasses import dataclass, field, fields, is_dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml


class ConfigError(ValueError):
    """Raised for malformed or inconsistent configuration."""


# ---------------------------------------------------------------------------
# Sections
# ---------------------------------------------------------------------------


@dataclass
class DatasetConfig:
    name: str = "dataset"
    #: Directory that relative input paths are resolved against.
    input_dir: str = "."
    description: str = ""


@dataclass
class MatrixInput:
    """A dense feature-by-cell (or cell-by-feature) text matrix.

    ``files`` may list several files that partition the cells (the SCP1064
    RNA matrix ships as six chunks); they are concatenated along cells and
    must share the same feature axis.
    """

    #: One file (or, for ``mtx``, one 10x directory); ``files`` lists several
    #: dense text files that partition the cells.
    file: Optional[str] = None
    files: List[str] = field(default_factory=list)
    #: Several lanes of the same modality, ``[{id: lane1, path: ...}, ...]``;
    #: cell IDs become ``<barcode>__<lane id>`` (see docs/DEFAULTS.md "Inputs").
    lanes: List[Dict[str, str]] = field(default_factory=list)
    #: ``dense_csv``: text table with a header row; first column = row names.
    #: ``mtx``: 10x Matrix Market directory (barcodes/features/matrix, gz or not).
    #: ``10x_h5``: 10x feature-barcode HDF5 file. ``h5ad``: AnnData file
    #: (``slot``/``key`` address the matrix inside it).
    format: str = "dense_csv"
    #: ``features_x_cells`` (SCP convention) or ``cells_x_features``; dense_csv only.
    orientation: str = "features_x_cells"
    #: Field separator; ``auto`` picks by extension (.csv -> ',', else tab).
    sep: str = "auto"
    #: Value state: ``auto`` | ``raw_counts`` | ``normalized``.
    state: str = "auto"
    #: Whether the modality must be present for a cell to be kept.
    required: bool = True
    #: h5ad only: where the matrix lives: ``X`` | ``layers`` | ``obsm`` (+ ``key``).
    slot: str = "X"
    key: Optional[str] = None
    #: mtx / 10x_h5 / h5ad: keep only features whose type is listed (empty = all).
    #: The type column is ``feature_type_column`` (10x: ``feature_types``).
    feature_types: List[str] = field(default_factory=list)
    feature_type_column: str = "feature_types"
    #: mtx / 10x_h5: which feature column becomes the feature ID (``id`` | ``name``).
    var_names: str = "id"

    def paths(self) -> List[str]:
        given = sum(bool(x) for x in (self.file, self.files, self.lanes))
        if given > 1:
            raise ConfigError("Give only one of 'file', 'files' or 'lanes' for a matrix input.")
        if self.file:
            return [self.file]
        if self.lanes:
            return [str(l.get("path", "")) for l in self.lanes]
        return list(self.files)

    def is_set(self) -> bool:
        return bool(self.file or self.files or self.lanes)


@dataclass
class TableInput:
    file: Optional[str] = None
    #: ``scp_metadata`` (second line is a TYPE row), ``csv``, ``tsv``, ``auto``.
    format: str = "auto"
    sep: str = "auto"
    required: bool = True


@dataclass
class GuideAssignmentInput:
    """Per-cell list of assigned guides (upstream guide calling already done)."""

    file: Optional[str] = None
    format: str = "auto"
    sep: str = "auto"
    cell_column: str = "Cell"
    guides_column: str = "sgRNAs"
    #: Separator between guides inside the guides column.
    list_separator: str = ","
    required: bool = False


@dataclass
class EmbeddingInput:
    """Optional precomputed embedding (e.g. the authors' UMAP) kept for diagnostics."""

    file: Optional[str] = None
    format: str = "auto"
    sep: str = "auto"
    key: str = "X_umap_provided"
    columns: List[str] = field(default_factory=lambda: ["X", "Y"])
    required: bool = False


@dataclass
class SlotAddress:
    """Where a modality lives inside a multiplexed h5ad: ``slot`` is ``X`` |
    ``layers`` | ``obsm`` | ``obs`` (guide assignments only), ``key`` names the
    layer / obsm matrix / obs column. ``feature_types`` (with the multiplexed
    ``feature_type_column``) selects features when the slot shares ``var``."""

    slot: Optional[str] = None
    key: Optional[str] = None
    #: obs-column guide lists only: separator between guides in one cell.
    list_separator: str = ";"


@dataclass
class MultiplexedSlots:
    rna: SlotAddress = field(default_factory=SlotAddress)
    protein_counts: SlotAddress = field(default_factory=SlotAddress)
    protein: SlotAddress = field(default_factory=SlotAddress)
    guide_counts: SlotAddress = field(default_factory=SlotAddress)
    guide_assignments: SlotAddress = field(default_factory=SlotAddress)

    def declared(self) -> List[str]:
        return [n for n in ("rna", "protein_counts", "protein", "guide_counts", "guide_assignments") if getattr(self, n).slot]


@dataclass
class MultiplexedFeatureTypes:
    """Feature-type labels that split one combined matrix into modalities.
    An empty list means "this modality is not in the multiplexed input"."""

    rna: List[str] = field(default_factory=lambda: ["Gene Expression"])
    protein: List[str] = field(default_factory=lambda: ["Antibody Capture"])
    guide: List[str] = field(default_factory=lambda: ["CRISPR Guide Capture"])


@dataclass
class MultiplexedInput:
    """One combined feature-barcode matrix (10x MTX directory, 10x H5 file or
    h5ad) that carries several modalities. It is read once and split by
    feature type (``feature_types``) or, for h5ad, by slot (``slots``).
    Protein features from a combined matrix are always raw ADT counts."""

    #: ``mtx`` | ``10x_h5`` | ``h5ad``.
    format: Optional[str] = None
    path: Optional[str] = None
    #: Several lanes: ``[{id: lane1, path: ...}, ...]`` (same feature set each).
    lanes: List[Dict[str, str]] = field(default_factory=list)
    feature_types: MultiplexedFeatureTypes = field(default_factory=MultiplexedFeatureTypes)
    feature_type_column: str = "feature_types"
    var_names: str = "id"
    #: h5ad only: explicit slot addressing (overrides feature-type splitting).
    slots: MultiplexedSlots = field(default_factory=MultiplexedSlots)
    #: h5ad only: value state of the RNA slot (``auto`` | ``raw_counts`` | ``normalized``).
    rna_state: str = "auto"

    def is_set(self) -> bool:
        return bool(self.path or self.lanes)

    def paths(self) -> List[str]:
        if self.path and self.lanes:
            raise ConfigError("inputs.multiplexed: give either 'path' or 'lanes', not both.")
        if self.path:
            return [self.path]
        return [str(l.get("path", "")) for l in self.lanes]

    def modalities(self) -> List[str]:
        """Modalities this block provides (``rna``, ``protein_counts``, ``protein``,
        ``guide_counts``, ``guide_assignments``)."""
        if self.format == "h5ad" and self.slots.declared():
            return self.slots.declared()
        out = []
        if self.feature_types.rna:
            out.append("rna")
        if self.feature_types.protein:
            out.append("protein_counts")
        if self.feature_types.guide:
            out.append("guide_counts")
        return out


@dataclass
class LaneMetadataInput:
    """Per-lane sample sheet joined onto cells by lane ID (multi-lane inputs)."""

    file: Optional[str] = None
    format: str = "auto"
    sep: str = "auto"
    lane_column: str = "lane_id"


@dataclass
class InputsConfig:
    #: A combined matrix carrying several modalities (see ``MultiplexedInput``).
    multiplexed: MultiplexedInput = field(default_factory=MultiplexedInput)
    rna: MatrixInput = field(default_factory=MatrixInput)
    #: Normalized protein values shipped by the authors (dense_csv, h5ad ...).
    protein: MatrixInput = field(default_factory=lambda: MatrixInput(required=True))
    #: Raw ADT UMI counts when they ship separately from the normalized matrix.
    protein_counts: MatrixInput = field(default_factory=lambda: MatrixInput(required=False))
    #: Guide (sgRNA) UMI counts per cell; enables guide calling and diagnostics.
    guide_counts: MatrixInput = field(default_factory=lambda: MatrixInput(required=False, state="raw_counts"))
    metadata: TableInput = field(default_factory=TableInput)
    lane_metadata: LaneMetadataInput = field(default_factory=LaneMetadataInput)
    guide_assignments: GuideAssignmentInput = field(default_factory=GuideAssignmentInput)
    embedding: EmbeddingInput = field(default_factory=EmbeddingInput)


@dataclass
class ColumnsConfig:
    """Column names inside the metadata table. ``null`` = not available."""

    cell_id: str = "NAME"
    perturbation: Optional[str] = None
    guide: Optional[str] = None
    target: Optional[str] = None
    condition: Optional[str] = None
    sample: Optional[str] = None
    donor: Optional[str] = None
    batch: Optional[str] = None
    replicate: Optional[str] = None
    lane: Optional[str] = None
    moi: Optional[str] = None
    #: Per-cell library size used to reconstruct counts from log-normalized data.
    rna_total_counts: Optional[str] = None
    #: Extra metadata columns to carry into ``obs`` unchanged.
    keep: List[str] = field(default_factory=list)


@dataclass
class PerturbationConfig:
    """Guide -> target parsing, control classes and perturbation-coverage QC."""

    #: Regex with a ``target`` group that extracts the target from a guide ID.
    guide_target_regex: str = r"^(?P<target>.+?)_(?P<index>\d+)$"
    #: Mapping control class -> list of regexes matched against the *target*.
    control_classes: Dict[str, List[str]] = field(
        default_factory=lambda: {
            "non_targeting": [r"^NO_SITE$", r"^NTC$", r"^non[-_ ]?targeting", r"scramble", r"^safe[-_ ]?harbor"],
        }
    )
    #: Optional explicit guide -> target table (csv/tsv with columns guide,target).
    guide_target_table: Optional[str] = None
    #: Coverage flags (nothing is removed): a guide/target with fewer cells is
    #: flagged ``low_coverage``. Wei Li pipeline: perturbation.min_cells_per_target = 10.
    min_cells_per_guide: int = 10
    min_cells_per_target: int = 10
    #: Keep multi-guide / unassigned cells in the processed object (recommended).
    preserve_multiguide: bool = True
    preserve_unassigned: bool = True
    #: Labels written to ``obs['perturbation']`` for non single-guide cells.
    multi_guide_label: str = "multi"
    unassigned_label: str = "unassigned"
    #: How the per-cell guide assignment used downstream is chosen / computed.
    assignment: "AssignmentConfig" = field(default_factory=lambda: AssignmentConfig())


@dataclass
class AssignmentConfig:
    """Guide-assignment source and the guide-calling rule for guide counts.

    ``source``: ``provided`` uses the supplied guide list / metadata column
    (guide counts, when present, only feed diagnostics and a comparison);
    ``guide_counts`` calls guides from the count matrix; ``auto`` picks the only
    available source and fails when both exist (choose explicitly).
    """

    source: str = "auto"
    #: ``dominant``: one guide per cell when the top guide has >= ``min_umi``
    #: UMIs and > ``dominance_ratio`` x the runner-up (weili-lab/perturbseq-pipeline
    #: rule); otherwise ``ambiguous`` (counts present) or ``unassigned`` (none).
    #: ``threshold``: every guide with >= ``detection_min_umi`` UMIs is assigned
    #: (multi-guide lists, as in high-MOI screens).
    method: str = "dominant"
    min_umi: int = 3
    dominance_ratio: float = 2.0
    #: Multiplet gate for ``dominant``: cells whose second guide exceeds this
    #: many UMIs are ``ambiguous``; -1 disables (reference default).
    max_second_umi: int = -1
    #: A guide counts as *detected* in a cell at >= this many UMIs (diagnostics
    #: ``n_guides_detected`` / ``guides_detected`` and the ``threshold`` rule).
    detection_min_umi: int = 3
    ambiguous_label: str = "ambiguous"


@dataclass
class AlignmentConfig:
    #: ``intersection``: keep cells present in every required modality.
    #: ``metadata``: keep the metadata cell set (required modalities must cover it).
    policy: str = "intersection"
    #: Modalities that must contain a cell for it to be kept.
    required: List[str] = field(default_factory=lambda: ["rna", "protein", "metadata"])
    #: Abort if any modality contains duplicate cell IDs.
    fail_on_duplicates: bool = True
    #: Abort if the retained set is smaller than this fraction of the RNA cells.
    min_overlap_fraction: float = 0.9


@dataclass
class HVGConfig:
    #: ``auto``: select HVGs when the matrix has more than ``n_top_genes`` genes.
    enabled: str = "auto"
    #: Wei Li pipeline default (cluster.n_top_genes: 3000).
    n_top_genes: int = 3000
    #: Scanpy flavor operating on log-normalized data (``seurat`` = Wei Li / Scanpy default).
    flavor: str = "seurat"
    batch_key: Optional[str] = None


@dataclass
class RNAConfig:
    #: ``auto``: detect raw vs log-normalized from the values; or force ``raw_counts`` / ``log_normalized``.
    input_state: str = "auto"
    #: ``auto``: normalize only raw counts; ``never``: keep values as-is;
    #: ``always``: normalize_total+log1p (refused on already-normalized input).
    normalize: str = "auto"
    #: Library-size target for raw-count input (Scanpy convention; Wei Li uses the median library size = null).
    target_sum: Optional[float] = 1.0e4
    log1p: bool = True
    #: Reconstruct integer counts from log-normalized values when
    #: ``columns.rna_total_counts`` is available: ``auto`` | ``never``.
    #: The result is stored as ``layers['reconstructed_counts']``, never as ``counts``.
    reconstruct_counts: str = "auto"
    #: Maximum tolerated |x - round(x)| for reconstructed counts to be accepted.
    reconstruct_tolerance: float = 0.05
    mito_prefix: str = "MT-"
    ribo_prefixes: List[str] = field(default_factory=lambda: ["RPS", "RPL"])
    hvg: HVGConfig = field(default_factory=HVGConfig)
    #: Scale the HVG matrix (zero-centre, unit variance, clipped) before PCA, as
    #: in the Wei Li pipeline (cluster.scale_max_value: 10). ``X`` is never scaled.
    scale: bool = True
    scale_max_value: Optional[float] = 10.0
    #: Wei Li / Scanpy default.
    n_pcs: int = 50


@dataclass
class ProteinConfig:
    #: Primary normalized matrix ``obsm['protein']``:
    #: ``auto``: keep provided normalized values when shipped, else CLR from counts;
    #: ``provided`` | ``clr`` | ``log1p`` | ``isotype_ratio`` | ``none``.
    normalization: str = "auto"
    #: Additional raw-derived matrices to store (``obsm['protein_<name>']``).
    extra_representations: List[str] = field(default_factory=lambda: ["clr"])
    #: Which matrix feeds protein PCA/UMAP: ``auto`` (clr if raw counts exist, else provided)
    #: | ``provided`` | ``clr`` | ``log1p``.
    embedding_representation: str = "auto"
    #: CLR geometric mean across ``cells`` (per protein; Seurat margin=2, recommended
    #: for small panels, Hao et al. 2021) or ``features`` (per cell; Stoeckius 2017).
    clr_axis: str = "cells"
    #: Substrings/regexes identifying isotype-control antibodies.
    isotype_patterns: List[str] = field(default_factory=lambda: ["IgG"])
    #: Explicit antibody -> isotype control map (optional; inferred when a provided
    #: isotype-ratio matrix can be reproduced from the counts).
    isotype_map: Dict[str, str] = field(default_factory=dict)
    #: Isotype controls never enter normalized matrices or embeddings ...
    exclude_isotypes_from_embedding: bool = True
    #: ... but stay in ``obsm['protein_counts']`` and drive the protein QC columns.
    use_isotypes_for_qc: bool = True
    #: Suffix stripped from feature names in the normalized matrix (e.g. " protein").
    strip_feature_suffix: Optional[str] = None
    #: Explicit rename map for protein names (normalized-matrix name -> canonical).
    rename: Dict[str, str] = field(default_factory=dict)
    #: Optional antibody annotation table (csv/tsv; relative to dataset.input_dir)
    #: with a ``feature_id`` column matching the canonical protein names and any
    #: of ``antibody_name``, ``protein_name``, ``gene_symbol``, ``clone``,
    #: ``isotype``, ``isotype_control``. Merged into ``uns['protein_features']``;
    #: values never given stay NA.
    feature_table: Optional[str] = None
    #: Scale (zero-centre, unit variance, clipped) before PCA.
    scale: bool = True
    scale_max_value: Optional[float] = 10.0
    #: Protein PCs: capped at n_features - 1; 10 PCs retain ~80 % of the variance of a ~20-plex panel (SCP1064).
    n_pcs: int = 10


@dataclass
class NeighborsConfig:
    #: Wei Li / Scanpy default.
    n_neighbors: int = 15
    #: RNA PCs used for the graph (null = all ``rna.n_pcs``).
    n_pcs: Optional[int] = None
    metric: str = "euclidean"


@dataclass
class UMAPConfig:
    enabled: bool = True
    #: Wei Li default (cluster.umap_min_dist: 0.5).
    min_dist: float = 0.5
    spread: float = 1.0
    #: UMAP is skipped for a modality with fewer features than this.
    min_features: int = 5
    #: obs columns to colour embeddings by (missing ones are skipped).
    color_by: List[str] = field(
        default_factory=lambda: ["condition", "perturbation_class", "sample", "donor", "batch"]
    )


@dataclass
class MultimodalConfig:
    """Optional joint RNA+protein representation.

    ``concat_pcs``: RNA PCs and protein PCs, each block divided by the square
    root of its total variance, protein block multiplied by ``protein_weight``,
    concatenated into ``obsm['X_multimodal']`` -> neighbors -> UMAP.
    Off by default: on SCP1064 the protein graph is dominated by ADT depth and
    shares ~2 % of neighbours with the RNA graph, so the joint space mostly
    dilutes RNA structure (see docs/DEFAULTS.md).
    """

    enabled: bool = False
    method: str = "concat_pcs"
    protein_weight: float = 1.0
    umap: bool = True
    #: Cells sampled for the RNA-vs-protein neighbourhood-agreement diagnostic (always computed).
    diagnostic_cells: int = 20000


@dataclass
class PrefilterConfig:
    """Permissive first pass: drop obvious empty droplets and never-detected genes
    *before* QC metrics are computed (Wei Li: min_genes_per_cell 200, min_cells_per_gene 3).
    Not a quality filter; recorded in the filtering audit as ``prefilter_*`` steps."""

    enabled: bool = True
    min_genes_per_cell: int = 200
    min_cells_per_gene: int = 3


@dataclass
class RNAFilterConfig:
    #: Strict cell thresholds (public defaults; see docs/DEFAULTS.md). ``null`` disables a step.
    min_genes: Optional[int] = 500
    min_counts: Optional[int] = None
    max_pct_mt: Optional[float] = 20.0


@dataclass
class ProteinFilterConfig:
    min_total_counts: Optional[int] = None
    min_proteins_detected: Optional[int] = None
    max_pct_isotype: Optional[float] = None
    #: Remove cells flagged ``protein_extreme_counts`` (total ADT > flags.extreme_fold x q99).
    remove_extreme_counts: bool = False


@dataclass
class PerturbationFilterConfig:
    #: ``all`` keeps every cell; ``assigned`` keeps cells with >= 1 guide;
    #: ``single_guide`` keeps single-guide cells only. Coverage is never a filter here.
    cells: str = "all"


@dataclass
class FilterConfig:
    """Strict, configured cell filtering. Every active step is written to
    ``tables/qc_filtering_steps.csv``. The same thresholds define the ``*_qc_fail``
    flags, so with ``enabled: false`` nothing is removed but the flags still exist."""

    enabled: bool = True
    rna: RNAFilterConfig = field(default_factory=RNAFilterConfig)
    protein: ProteinFilterConfig = field(default_factory=ProteinFilterConfig)
    perturbation: PerturbationFilterConfig = field(default_factory=PerturbationFilterConfig)


@dataclass
class FlagsConfig:
    """Diagnostic flags that never remove cells."""

    #: MAD outlier flags on log1p(total_counts) / log1p(n_genes) (Wei Li filters at 3; 5 flags only clear outliers).
    rna_n_mads: float = 5.0
    protein_n_mads: float = 5.0
    #: ``protein_extreme_counts``: total ADT > extreme_fold x the 99th percentile (antibody aggregates).
    extreme_fold: float = 10.0
    #: Antibody is ``background_dominated`` when fewer than this fraction of cells exceed its isotype control.
    background_min_fraction_above_isotype: float = 0.5


@dataclass
class QCConfig:
    prefilter: PrefilterConfig = field(default_factory=PrefilterConfig)
    filter: FilterConfig = field(default_factory=FilterConfig)
    flags: FlagsConfig = field(default_factory=FlagsConfig)


@dataclass
class ReportConfig:
    #: Report title (default: "<run name> - Perturb-CITE-seq preprocessing report").
    title: Optional[str] = None
    #: Embed figures as base64 data URIs so report.html is self-contained. Figures are always also written to figures/.
    embed_figures: bool = True
    write_markdown: bool = True
    figure_format: str = "png"
    figure_dpi: int = 120
    max_table_rows: int = 100


@dataclass
class RunConfig:
    #: Run name used in the report and archive name (default: dataset.name).
    name: Optional[str] = None


@dataclass
class ComputeConfig:
    #: Parallel workers for reading multi-file matrices (processes).
    n_jobs: int = 1
    #: Rows per chunk when streaming dense text matrices.
    chunk_rows: int = 500
    #: Single random seed for PCA, neighbors and UMAP (Wei Li run.seed = 0).
    seed: int = 0
    #: Optional: restrict to the first N cells of every modality (smoke tests only).
    max_cells: Optional[int] = None


@dataclass
class PSAnalysisConfig:
    """Perturbation-response score (Song et al. 2025, Nat Cell Biol; PS_python /
    pertps definition, re-implemented in ``analysis/ps_score.py``)."""

    enabled: bool = True
    #: response genes per target: top-N genes ranked by t-statistic target vs control
    top_n_genes: int = 100
    #: scores are clipped to [0, scale_factor] before max-normalization to [0, 1]
    scale_factor: float = 3.0
    #: quadrant cut on the normalized score
    ps_threshold: float = 0.5
    #: horizontal quadrant cut on the target's own expression: mean | median | quantile of controls
    expression_cut: str = "mean"
    expression_cut_quantile: float = 0.75
    min_cells_per_target: int = 10
    min_control_cells: int = 10
    #: skip targets expressed in fewer than this % of control cells (knockdown unmeasurable)
    min_pct_expressing_control: float = 1.0


@dataclass
class LochnessAnalysisConfig:
    """lochNESS local neighbourhood enrichment (Huang et al. 2023, Nature; pertTF
    port as in weili-lab/perturbseq-pipeline), ``analysis/lochness.py``."""

    enabled: bool = True
    #: k nearest neighbours in PCA space (reference: 300)
    n_neighbors: int = 300
    #: cap k at this fraction of the cell count (small objects); the cap is recorded when it applies
    max_k_fraction: float = 0.1
    n_pcs: int = 20
    #: obsm key for the neighbour space; null -> obsm['X_pca'] (never UMAP)
    use_rep: Optional[str] = None
    min_cells_per_target: int = 10
    #: score above which a cell counts as 'enriched' (pct_cells_enriched)
    enrichment_cut: float = 0.5
    #: label-permutation null per target for the own-cell mean (0 disables)
    n_permutations: int = 200


@dataclass
class ModulesAnalysisConfig:
    """Perturbation x response-gene effect matrix, gene programs (gene axis) and
    perturbation modules (perturbation axis), after Zhou et al. 2023 (Nature)
    as implemented in weili-lab/perturbseq-pipeline; ``analysis/modules.py``."""

    enabled: bool = True
    min_cells_per_perturbation: int = 20
    min_perturbations: int = 5
    min_genes: int = 10
    #: gene panel: ``response`` = union of each perturbation's top response genes
    #: (Welch t-test vs control, BH within perturbation); ``hvg`` = highly variable genes
    gene_selection: str = "response"
    n_response_genes_per_perturbation: int = 100
    response_fdr: float = 0.05
    #: correlation used to cluster genes (programs) and perturbations (modules)
    program_correlation: str = "pearson"
    module_correlation: str = "spearman"
    linkage_method: str = "average"
    #: number of clusters, or null for a dendrogram cut at cluster_distance_threshold
    n_programs: Optional[int] = 4
    n_modules: Optional[int] = 9
    cluster_distance_threshold: float = 0.7
    #: pseudocount added to the de-logged group means before log2 (normalized-count
    #: scale). 1.0 = Seurat FoldChange convention; the reference used 1e-9, which makes
    #: log2FC unbounded for genes undetected in a small group
    log2fc_pseudocount: float = 1.0
    #: panel genes must be detected (> 0) in at least this % of perturbed + control cells
    min_pct_cells_expressing: float = 5.0
    #: DE gate for 'responding' counts: |log2FC| > threshold and BH-FDR < alpha
    de_lfc_threshold: float = 0.5
    de_fdr_alpha: float = 0.05
    #: per-cell program activity (mean z-scored expression of the program's genes)
    score_programs: bool = True


@dataclass
class ProteinEffectsConfig:
    """Perturbation effects on each measured protein (``analysis/protein_effects.py``)."""

    enabled: bool = True
    #: obsm key with the normalized protein values (CLR by default; never counts)
    representation: str = "protein"
    min_cells_per_target: int = 10
    min_control_cells: int = 10
    fdr_alpha: float = 0.05


@dataclass
class ConcordanceConfig:
    """RNA-protein concordance: PS <-> protein, lochNESS <-> protein, program <-> protein
    (``analysis/concordance.py``)."""

    enabled: bool = True
    #: minimum cells for a within-target Spearman correlation
    min_cells: int = 20
    fdr_alpha: float = 0.05


@dataclass
class PerturbationEffectsConfig:
    """Stage E: downstream perturbation-effect analyses on the processed object.
    Off by default; the sub-analyses run when this block is enabled (each can be
    switched off individually). Controls are the ``single_control`` cells whose
    control class is listed in ``control_classes``; perturbed cells are
    ``single_targeting`` cells of one target. Ambiguous / multi-guide cells never
    enter either group."""

    enabled: bool = False
    control_classes: List[str] = field(default_factory=lambda: ["non_targeting"])
    #: target label -> gene symbol in var when they differ (e.g. PDL1 -> CD274);
    #: used for the PS expressed-in-controls guard and the knockdown quadrants
    target_gene_map: Dict[str, str] = field(default_factory=dict)
    ps: PSAnalysisConfig = field(default_factory=PSAnalysisConfig)
    lochness: LochnessAnalysisConfig = field(default_factory=LochnessAnalysisConfig)
    modules: ModulesAnalysisConfig = field(default_factory=ModulesAnalysisConfig)
    protein: ProteinEffectsConfig = field(default_factory=ProteinEffectsConfig)
    concordance: ConcordanceConfig = field(default_factory=ConcordanceConfig)
    #: targets shown in the report figures (tables are always complete)
    top_n_report: int = 12


@dataclass
class ClusterEnrichmentConfig:
    """Perturbation x cluster enrichment (``analysis/cluster_enrichment.py``)."""

    enabled: bool = True
    #: ``non_targeting``: single-guide cells of ``control_classes`` (default);
    #: ``other``: single-guide targeting cells of every other target (explicit
    #: alternative, the default of weili-lab/perturbseq-pipeline)
    control: str = "non_targeting"
    control_classes: List[str] = field(default_factory=lambda: ["non_targeting"])
    min_cells_per_target: int = 10
    #: clusters smaller than this are not tested
    min_cells_per_cluster: int = 20
    min_control_cells: int = 10
    #: pairs with fewer control cells in the cluster are flagged ``low_power``
    min_control_cells_in_cluster: int = 10
    #: guides with at least this many cells enter the guide-support count
    min_cells_per_guide: int = 5
    #: Haldane-Anscombe pseudocount for the finite log2 odds ratio (display/ranking only)
    odds_pseudocount: float = 0.5
    fdr_alpha: float = 0.05
    #: optional obs column (e.g. sample, lane): adds a Cochran-Mantel-Haenszel test across its levels
    stratify_by: Optional[str] = None
    #: label permutations for the omnibus target x cluster chi-square (0 disables)
    n_permutations: int = 1000
    top_n_report: int = 12


@dataclass
class ClusteringConfig:
    """Stage F: Leiden clustering of the RNA neighbour graph plus perturbation x
    cluster enrichment. Off by default; independent of ``perturbation_effects``."""

    enabled: bool = False
    #: obs column for the cluster labels
    key: str = "leiden"
    #: Leiden resolution (higher -> more clusters); reference default 1.0
    resolution: float = 1.0
    #: Leiden iterations (reference 2; -1 = until convergence)
    n_iterations: int = 2
    #: the RNA neighbour graph built in the representation stage
    neighbors_key: str = "rna"
    #: overwrite an existing obs column named ``key`` (otherwise an error)
    overwrite: bool = False
    enrichment: ClusterEnrichmentConfig = field(default_factory=ClusterEnrichmentConfig)


@dataclass
class AnalysisConfig:
    perturbation_effects: PerturbationEffectsConfig = field(default_factory=PerturbationEffectsConfig)
    clustering: ClusteringConfig = field(default_factory=ClusteringConfig)


@dataclass
class OutputConfig:
    dir: str = "results/run"
    h5ad_name: str = "processed.h5ad"
    compression: Optional[str] = "gzip"
    write_h5ad: bool = True
    #: Also write the pre-strict-filter object (processed/<stem>_prefilter.h5ad). Off: the
    #: before-filter figures and tables/cell_qc_prefilter.csv.gz cover this without doubling storage.
    write_prefilter_h5ad: bool = False
    #: Bundle report, figures, tables and logs into <run>_results.tar.gz (matrices excluded).
    archive: bool = False
    archive_name: Optional[str] = None
    archive_exclude: List[str] = field(default_factory=lambda: ["*.h5ad", "*.h5", "*.tar.gz", "*.csv.gz"])


@dataclass
class Config:
    run: RunConfig = field(default_factory=RunConfig)
    dataset: DatasetConfig = field(default_factory=DatasetConfig)
    inputs: InputsConfig = field(default_factory=InputsConfig)
    columns: ColumnsConfig = field(default_factory=ColumnsConfig)
    perturbation: PerturbationConfig = field(default_factory=PerturbationConfig)
    alignment: AlignmentConfig = field(default_factory=AlignmentConfig)
    rna: RNAConfig = field(default_factory=RNAConfig)
    protein: ProteinConfig = field(default_factory=ProteinConfig)
    neighbors: NeighborsConfig = field(default_factory=NeighborsConfig)
    umap: UMAPConfig = field(default_factory=UMAPConfig)
    multimodal: MultimodalConfig = field(default_factory=MultimodalConfig)
    qc: QCConfig = field(default_factory=QCConfig)
    report: ReportConfig = field(default_factory=ReportConfig)
    compute: ComputeConfig = field(default_factory=ComputeConfig)
    output: OutputConfig = field(default_factory=OutputConfig)
    analysis: AnalysisConfig = field(default_factory=AnalysisConfig)

    #: Path the config was loaded from (not part of the schema).
    source_path: Optional[str] = field(default=None, compare=False)

    # -- construction -----------------------------------------------------

    @classmethod
    def from_dict(cls, data: Dict[str, Any], source_path: Optional[str] = None) -> "Config":
        cfg = _build(cls, data or {}, path="")
        cfg.source_path = source_path
        cfg.validate()
        return cfg

    @classmethod
    def from_yaml(cls, path: str | Path) -> "Config":
        path = Path(path)
        if not path.exists():
            raise ConfigError(f"Config file not found: {path}")
        with open(path) as fh:
            data = yaml.safe_load(fh) or {}
        if not isinstance(data, dict):
            raise ConfigError(f"Config root must be a mapping: {path}")
        data = expand_roots(data, path_roots(path))
        return cls.from_dict(data, source_path=str(path.resolve()))

    def to_dict(self) -> Dict[str, Any]:
        d = dataclasses.asdict(self)
        d.pop("source_path", None)
        return d

    def to_yaml(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as fh:
            yaml.safe_dump(self.to_dict(), fh, sort_keys=False, allow_unicode=True)
        return path

    # -- helpers ----------------------------------------------------------

    def resolve(self, rel: Optional[str]) -> Optional[Path]:
        """Resolve a path from the config against ``dataset.input_dir``."""
        if rel is None:
            return None
        p = Path(rel).expanduser()
        if p.is_absolute():
            return p
        return (Path(self.dataset.input_dir).expanduser() / p)

    def output_dir(self) -> Path:
        return Path(self.output.dir).expanduser()

    def run_name(self) -> str:
        return self.run.name or self.dataset.name

    @classmethod
    def defaults_dict(cls) -> Dict[str, Any]:
        """Documented defaults as a plain dict (for ``init-config``); not validated."""
        return cls().to_dict()

    def validate(self) -> None:
        errs: List[str] = []
        mp = self.inputs.multiplexed
        mp_mods: List[str] = []
        try:
            mp_paths = mp.paths()
        except ConfigError as exc:
            errs.append(str(exc))
            mp_paths = []
        if mp.is_set():
            _choice(errs, "inputs.multiplexed.format", mp.format, ("mtx", "10x_h5", "h5ad"))
            _choice(errs, "inputs.multiplexed.rna_state", mp.rna_state, ("auto", "raw_counts", "normalized"))
            if mp.lanes:
                errs.extend(_check_lanes(mp.lanes, "inputs.multiplexed.lanes"))
            for n in mp.slots.declared():
                sa = getattr(mp.slots, n)
                allowed = ("obs",) if n == "guide_assignments" else ("X", "layers", "obsm")
                if sa.slot not in allowed:
                    errs.append(f"inputs.multiplexed.slots.{n}.slot must be one of {allowed}")
                if sa.slot in ("layers", "obsm", "obs") and not sa.key:
                    errs.append(f"inputs.multiplexed.slots.{n}.key is required for slot '{sa.slot}'")
            if mp.slots.declared() and mp.format != "h5ad":
                errs.append("inputs.multiplexed.slots is only meaningful for format 'h5ad'")
            mp_mods = mp.modalities() if not errs else []
            if not mp_mods:
                errs.append("inputs.multiplexed declares no modality (feature_types / slots are all empty)")
        elif mp.format or mp.slots.declared():
            errs.append("inputs.multiplexed: 'path' or 'lanes' is required when the block is used")
        for name in ("rna", "protein", "protein_counts", "guide_counts"):
            m: MatrixInput = getattr(self.inputs, name)
            try:
                m.paths()
            except ConfigError as exc:
                errs.append(f"inputs.{name}: {exc}")
                continue
            _choice(errs, f"inputs.{name}.format", m.format, ("dense_csv", "mtx", "10x_h5", "h5ad"))
            if m.orientation not in ("features_x_cells", "cells_x_features"):
                errs.append(f"inputs.{name}.orientation must be features_x_cells or cells_x_features")
            if m.state not in ("auto", "raw_counts", "normalized"):
                errs.append(f"inputs.{name}.state must be auto, raw_counts or normalized")
            _choice(errs, f"inputs.{name}.slot", m.slot, ("X", "layers", "obsm"))
            if m.slot in ("layers", "obsm") and not m.key:
                errs.append(f"inputs.{name}.key is required for slot '{m.slot}'")
            _choice(errs, f"inputs.{name}.var_names", m.var_names, ("id", "name"))
            if m.lanes:
                errs.extend(_check_lanes(m.lanes, f"inputs.{name}.lanes"))
            if m.format == "dense_csv" and (m.lanes or m.feature_types):
                errs.append(f"inputs.{name}: 'lanes' and 'feature_types' need format mtx, 10x_h5 or h5ad")
            if m.is_set() and name in mp_mods:
                errs.append(f"inputs.{name} is declared both in inputs.multiplexed and as a separate input; remove one (no override exists)")
        if self.inputs.guide_assignments.file and "guide_assignments" in mp_mods:
            errs.append("inputs.guide_assignments is declared both in inputs.multiplexed.slots and as a separate file; remove one")
        if not self.inputs.rna.is_set() and "rna" not in mp_mods:
            errs.append("inputs.rna: at least one file is required (or inputs.multiplexed must provide RNA)")
        if self.inputs.lane_metadata.file and not (mp.lanes or any(getattr(self.inputs, n).lanes for n in ("rna", "protein", "protein_counts", "guide_counts"))):
            errs.append("inputs.lane_metadata needs multi-lane inputs ('lanes')")
        if self.alignment.policy not in ("intersection", "metadata"):
            errs.append("alignment.policy must be 'intersection' or 'metadata'")
        for mod in self.alignment.required:
            if mod not in ("rna", "protein", "protein_counts", "guide_counts", "metadata", "guide_assignments"):
                errs.append(f"alignment.required: unknown modality '{mod}'")
        asg = self.perturbation.assignment
        _choice(errs, "perturbation.assignment.source", asg.source, ("auto", "provided", "guide_counts"))
        _choice(errs, "perturbation.assignment.method", asg.method, ("dominant", "threshold"))
        for k, v in (("perturbation.assignment.min_umi", asg.min_umi), ("perturbation.assignment.detection_min_umi", asg.detection_min_umi)):
            if v < 1:
                errs.append(f"{k} must be >= 1")
        if asg.dominance_ratio < 1:
            errs.append("perturbation.assignment.dominance_ratio must be >= 1")
        _choice(errs, "rna.input_state", self.rna.input_state, ("auto", "raw_counts", "log_normalized"))
        _choice(errs, "rna.normalize", self.rna.normalize, ("auto", "never", "always"))
        _choice(errs, "rna.reconstruct_counts", self.rna.reconstruct_counts, ("auto", "never"))
        _choice(errs, "rna.hvg.enabled", self.rna.hvg.enabled, ("auto", "true", "false"))
        _choice(errs, "rna.hvg.flavor", self.rna.hvg.flavor, ("seurat", "cell_ranger", "seurat_v3"))
        _choice(errs, "protein.normalization", self.protein.normalization, ("auto", "provided", "clr", "log1p", "isotype_ratio", "none"))
        for r in self.protein.extra_representations:
            _choice(errs, "protein.extra_representations", r, ("clr", "log1p", "isotype_ratio"))
        _choice(errs, "protein.embedding_representation", self.protein.embedding_representation, ("auto", "provided", "clr", "log1p"))
        _choice(errs, "protein.clr_axis", self.protein.clr_axis, ("features", "cells"))
        _choice(errs, "multimodal.method", self.multimodal.method, ("concat_pcs",))
        _choice(errs, "qc.filter.perturbation.cells", self.qc.filter.perturbation.cells, ("all", "assigned", "single_guide"))
        _choice(errs, "report.figure_format", self.report.figure_format, ("png", "pdf", "svg"))
        for k, v in (("qc.prefilter.min_genes_per_cell", self.qc.prefilter.min_genes_per_cell), ("qc.prefilter.min_cells_per_gene", self.qc.prefilter.min_cells_per_gene)):
            if v < 0:
                errs.append(f"{k} must be >= 0")
        for k, v in (("rna.n_pcs", self.rna.n_pcs), ("protein.n_pcs", self.protein.n_pcs), ("neighbors.n_neighbors", self.neighbors.n_neighbors), ("rna.hvg.n_top_genes", self.rna.hvg.n_top_genes), ("perturbation.min_cells_per_guide", self.perturbation.min_cells_per_guide), ("perturbation.min_cells_per_target", self.perturbation.min_cells_per_target)):
            if v < 1:
                errs.append(f"{k} must be >= 1")
        if self.multimodal.protein_weight < 0:
            errs.append("multimodal.protein_weight must be >= 0")
        pe = self.analysis.perturbation_effects
        _choice(errs, "analysis.perturbation_effects.ps.expression_cut", pe.ps.expression_cut, ("mean", "median", "quantile"))
        _choice(errs, "analysis.perturbation_effects.modules.gene_selection", pe.modules.gene_selection, ("response", "hvg"))
        for k, v in (("modules.program_correlation", pe.modules.program_correlation), ("modules.module_correlation", pe.modules.module_correlation)):
            _choice(errs, f"analysis.perturbation_effects.{k}", v, ("pearson", "spearman"))
        _choice(errs, "analysis.perturbation_effects.modules.linkage_method", pe.modules.linkage_method, ("average", "complete", "single", "ward"))
        for k, v in (("ps.top_n_genes", pe.ps.top_n_genes), ("ps.min_cells_per_target", pe.ps.min_cells_per_target), ("lochness.n_neighbors", pe.lochness.n_neighbors), ("lochness.min_cells_per_target", pe.lochness.min_cells_per_target), ("modules.min_cells_per_perturbation", pe.modules.min_cells_per_perturbation), ("protein.min_cells_per_target", pe.protein.min_cells_per_target), ("concordance.min_cells", pe.concordance.min_cells)):
            if v < 1:
                errs.append(f"analysis.perturbation_effects.{k} must be >= 1")
        if pe.ps.scale_factor <= 0:
            errs.append("analysis.perturbation_effects.ps.scale_factor must be > 0")
        if not 0 < pe.lochness.max_k_fraction <= 1:
            errs.append("analysis.perturbation_effects.lochness.max_k_fraction must be in (0, 1]")
        if not pe.control_classes:
            errs.append("analysis.perturbation_effects.control_classes must name at least one control class")
        cl = self.analysis.clustering
        if cl.resolution <= 0:
            errs.append("analysis.clustering.resolution must be > 0")
        if cl.n_iterations == 0 or cl.n_iterations < -1:
            errs.append("analysis.clustering.n_iterations must be >= 1 or -1 (until convergence)")
        if not str(cl.key).strip():
            errs.append("analysis.clustering.key must be a non-empty obs column name")
        ce = cl.enrichment
        _choice(errs, "analysis.clustering.enrichment.control", ce.control, ("non_targeting", "other"))
        if ce.control == "non_targeting" and not ce.control_classes:
            errs.append("analysis.clustering.enrichment.control_classes must name at least one control class")
        for k, v in (("min_cells_per_target", ce.min_cells_per_target), ("min_cells_per_cluster", ce.min_cells_per_cluster), ("min_control_cells", ce.min_control_cells), ("min_cells_per_guide", ce.min_cells_per_guide)):
            if v < 1:
                errs.append(f"analysis.clustering.enrichment.{k} must be >= 1")
        if ce.odds_pseudocount < 0 or not 0 < ce.fdr_alpha < 1 or ce.n_permutations < 0:
            errs.append("analysis.clustering.enrichment: odds_pseudocount >= 0, 0 < fdr_alpha < 1 and n_permutations >= 0 are required")
        if "baselines" in Path(self.output.dir).expanduser().parts:
            errs.append(f"output.dir {self.output.dir!r} lies under a 'baselines' directory, which holds immutable regression references (docs/REGRESSION.md); choose another output directory")
        if errs:
            raise ConfigError("Invalid configuration:\n  - " + "\n  - ".join(errs))


def _check_lanes(lanes: List[Dict[str, str]], key: str) -> List[str]:
    """Validate a ``[{id, path}, ...]`` lane list: mappings, both keys, unique IDs."""
    errs: List[str] = []
    ids: List[str] = []
    for i, l in enumerate(lanes):
        if not isinstance(l, dict) or not l.get("id") or not l.get("path"):
            errs.append(f"{key}[{i}] must be a mapping with 'id' and 'path'")
            continue
        ids.append(str(l["id"]))
    dup = sorted({i for i in ids if ids.count(i) > 1})
    if dup:
        errs.append(f"{key}: duplicate lane id(s) {dup}")
    return errs


def _choice(errs: List[str], key: str, value: Any, allowed: tuple) -> None:
    if value not in allowed:
        errs.append(f"{key} must be one of {list(allowed)} (got {value!r})")


# ---------------------------------------------------------------------------
# Path roots: ${DATA_ROOT} / ${RESULTS_ROOT} / ${REPO_ROOT} placeholders
# ---------------------------------------------------------------------------
#
# Dataset YAMLs should not hard-code one HPC's absolute paths. Three
# placeholders are expanded in every string value of a config file:
#
#   ${REPO_ROOT}     the repository (the parent of the ``config/`` directory
#                    holding the YAML; or the YAML's own directory otherwise)
#   ${DATA_ROOT}     $PETRUBSEQ_DATA_ROOT, else <REPO_ROOT>/../data
#   ${RESULTS_ROOT}  $PETRUBSEQ_RESULTS_ROOT, else <REPO_ROOT>/../results
#
# This is the whole path-management layer; absolute paths still work.

ROOT_ENV = {"DATA_ROOT": "PETRUBSEQ_DATA_ROOT", "RESULTS_ROOT": "PETRUBSEQ_RESULTS_ROOT"}


def path_roots(config_path: str | Path) -> Dict[str, str]:
    """Resolve the placeholder roots for a config file location."""
    cfg_dir = Path(config_path).resolve().parent
    repo = cfg_dir.parent if cfg_dir.name == "config" else cfg_dir
    roots = {"REPO_ROOT": str(repo)}
    for key, env in ROOT_ENV.items():
        default = (repo.parent / key.split("_")[0].lower()).resolve()  # ../data, ../results
        roots[key] = os.environ.get(env) or str(default)
    return roots


def expand_roots(data: Any, roots: Dict[str, str]) -> Any:
    """Recursively substitute ``${NAME}`` placeholders in string values."""
    if isinstance(data, dict):
        return {k: expand_roots(v, roots) for k, v in data.items()}
    if isinstance(data, list):
        return [expand_roots(v, roots) for v in data]
    if isinstance(data, str) and "${" in data:
        out = data
        for key, val in roots.items():
            out = out.replace("${" + key + "}", val)
        if "${" in out:
            raise ConfigError(f"Unknown path placeholder in {data!r}; known: {sorted(roots)}")
        return out
    return data


# ---------------------------------------------------------------------------
# Generic dict -> dataclass builder with unknown-key detection
# ---------------------------------------------------------------------------


def _build(cls, data: Dict[str, Any], path: str):
    if not isinstance(data, dict):
        raise ConfigError(f"Section '{path or '<root>'}' must be a mapping, got {type(data).__name__}")
    known = {f.name: f for f in fields(cls) if f.name != "source_path"}
    unknown = sorted(set(data) - set(known))
    if unknown:
        raise ConfigError(f"Unknown key(s) in '{path or '<root>'}': {', '.join(unknown)}")
    kwargs: Dict[str, Any] = {}
    for name, f in known.items():
        if name not in data:
            continue
        value = data[name]
        sub_path = f"{path}.{name}" if path else name
        ftype = f.type if not isinstance(f.type, str) else None
        default = f.default if f.default is not dataclasses.MISSING else (
            f.default_factory() if f.default_factory is not dataclasses.MISSING else None  # type: ignore[misc]
        )
        if is_dataclass(default) and value is not None:
            kwargs[name] = _build(type(default), value, sub_path)
        else:
            kwargs[name] = _coerce(value, default, sub_path)
    return cls(**kwargs)


def _coerce(value: Any, default: Any, path: str) -> Any:
    """Light type coercion so YAML scalars land in the expected Python type."""
    if value is None:
        return None
    if isinstance(default, bool):
        if isinstance(value, bool):
            return value
        if isinstance(value, str) and value.lower() in ("true", "false"):
            return value.lower() == "true"
        raise ConfigError(f"'{path}' must be a boolean")
    if isinstance(default, str):
        # Allow YAML booleans for keys like embeddings.rna.hvg ("auto"|"true"|"false").
        if isinstance(value, bool):
            return str(value).lower()
        return str(value)
    if isinstance(default, int) and not isinstance(default, bool):
        if isinstance(value, (int, float)) and float(value).is_integer():
            return int(value)
        raise ConfigError(f"'{path}' must be an integer")
    if isinstance(default, float):
        if isinstance(value, (int, float)):
            return float(value)
        if isinstance(value, str):
            try:
                return float(value)
            except ValueError:
                pass
        raise ConfigError(f"'{path}' must be a number")
    if isinstance(default, list):
        if isinstance(value, (list, tuple)):
            return list(value)
        return [value]
    if isinstance(default, dict):
        if not isinstance(value, dict):
            raise ConfigError(f"'{path}' must be a mapping")
        return dict(value)
    return value
