"""Downstream analyses on the processed multimodal object, in the reference order.

Reference RNA analyses (weili-lab/perturbseq-pipeline behaviour, see
docs/reference/REFERENCE_PIPELINE_COMPLETE_AUDIT.md), run in this order once the
representations exist:

    Leiden clustering  ->  perturbation strength  ->  perturbation x cluster enrichment
    ->  co-functional modules / gene programs  ->  PS  ->  lochNESS

followed by the protein extension (protein effects, RNA-protein concordance).
``run_clustering`` / ``run_enrichment`` (``analysis.clustering``) and
``run_perturbation_effects`` (``analysis.perturbation_effects``) are the entry
points; ``attach`` writes the results into the AnnData and ``tables`` returns every
result table keyed by its path under ``tables/``.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, Dict, Optional

import anndata as ad
import numpy as np
import pandas as pd

from ..config import Config
from .concordance import ConcordanceResults, compute_concordance
from .lochness import LochnessResults, compute_lochness
from .modules import ModulesResults, compute_modules
from .perturbation_strength import StrengthResults, compute_perturbation_strength
from .protein_effects import ProteinEffectsResults, compute_protein_effects
from .ps_score import PSResults, compare_with_strength, compute_ps

logger = logging.getLogger(__name__)


@dataclass
class PerturbationEffects:
    strength: Optional[StrengthResults] = None
    ps: Optional[PSResults] = None
    lochness: Optional[LochnessResults] = None
    modules: Optional[ModulesResults] = None
    protein: Optional[ProteinEffectsResults] = None
    concordance: Optional[ConcordanceResults] = None
    ps_vs_strength: Optional[pd.DataFrame] = None
    notes: Optional[list] = None

    def info(self) -> Dict[str, Any]:
        out: Dict[str, Any] = {"enabled": True}
        for name in ("strength", "ps", "lochness", "modules", "protein", "concordance"):
            r = getattr(self, name)
            out[name] = dict(r.info) if r is not None else {"status": "disabled"}
        return out


def run_perturbation_strength(adata: ad.AnnData, cfg: Config) -> Optional[StrengthResults]:
    pe = cfg.analysis.perturbation_effects
    if not pe.strength.enabled:
        return None
    return compute_perturbation_strength(adata, cfg)


def run_perturbation_effects(adata: ad.AnnData, cfg: Config, strength: Optional[StrengthResults] = None) -> PerturbationEffects:
    """Modules -> PS -> lochNESS -> protein effects -> concordance (the reference order after enrichment)."""
    pe = cfg.analysis.perturbation_effects
    res = PerturbationEffects(strength=strength, notes=[])
    if pe.modules.enabled:
        res.modules = compute_modules(adata, cfg)
    if pe.ps.enabled:
        res.ps = compute_ps(adata, cfg)
        res.ps_vs_strength = compare_with_strength(res.ps, strength)
    if pe.lochness.enabled:
        res.lochness = compute_lochness(adata, cfg)
    if pe.protein.enabled:
        res.protein = compute_protein_effects(adata, cfg)
    if pe.concordance.enabled:
        res.concordance = compute_concordance(adata, cfg, res.ps, res.lochness, res.modules, res.protein)
    return res


def _h5(df: pd.DataFrame) -> pd.DataFrame:
    """h5ad-safe copy: object columns as strings (NaN -> '')."""
    df = df.copy()
    for c in df.columns:
        if df[c].dtype == object:
            df[c] = df[c].where(df[c].notna(), "").astype(str)
    df.columns = [str(c) for c in df.columns]
    df.index = df.index.astype(str)
    return df


def attach(adata: ad.AnnData, res: PerturbationEffects) -> None:
    """Write per-cell results into obs/obsm and the tables into uns."""
    if res.strength is not None and not res.strength.empty:
        adata.uns["perturbation_strength"] = _h5(res.strength.table)
    if res.ps is not None and not res.ps.scores.empty:
        adata.obsm["ps_scores"] = res.ps.scores.astype("float32")
        adata.obsm["ps_scores_raw"] = res.ps.raw.astype("float32")
        adata.obs["ps_score"] = res.ps.own.to_numpy(dtype="float32")
        adata.obs["ps_quadrant"] = pd.Categorical(res.ps.own_quadrant.astype(str))
        adata.uns["ps_signatures"] = _h5(res.ps.signatures)
        if res.ps.lda_umap is not None:
            adata.obsm["X_lda_umap"] = np.asarray(res.ps.lda_umap, dtype=np.float32)
            if res.ps.lda_label is not None:
                adata.obs["lda_label"] = pd.Categorical(res.ps.lda_label.reindex(adata.obs_names).astype(str))
    if res.lochness is not None and not res.lochness.scores.empty:
        adata.obsm["lochness"] = res.lochness.scores.astype("float32")
        adata.obs["lochness_self"] = res.lochness.self_score.to_numpy(dtype="float32")
    if res.modules is not None and not res.modules.empty:
        if res.modules.program_activity is not None:
            adata.obsm["program_activity"] = res.modules.program_activity
        adata.uns["perturbation_effect_matrix"] = _h5(res.modules.effect.astype("float32"))
        adata.uns["gene_programs"] = _h5(res.modules.gene_programs)
        adata.uns["perturbation_modules"] = _h5(res.modules.perturbation_modules)
    if res.protein is not None and not res.protein.empty:
        adata.uns["protein_effects"] = _h5(res.protein.table)
    if res.concordance is not None and not res.concordance.summary.empty:
        adata.uns["perturbation_summary"] = _h5(res.concordance.summary)


def tables(res: PerturbationEffects) -> Dict[str, pd.DataFrame]:
    """Every result table, keyed by its path under tables/ (without extension)."""
    out: Dict[str, pd.DataFrame] = {}
    d = "perturbation_effects"
    if res.strength is not None:
        out["perturbation_strength/perturbation_full"] = res.strength.table
        out["perturbation_strength/perturbation"] = res.strength.display
        out["perturbation_strength/skipped"] = res.strength.skipped
    if res.ps is not None:
        out[f"{d}/ps_targets"] = res.ps.summary
        out[f"{d}/ps_skipped"] = res.ps.skipped
        out[f"{d}/ps_signatures"] = res.ps.signatures
        if res.ps_vs_strength is not None and not res.ps_vs_strength.empty:
            out[f"{d}/ps_vs_perturbation"] = res.ps_vs_strength
    if res.lochness is not None:
        out[f"{d}/lochness_targets"] = res.lochness.summary
        out[f"{d}/lochness_skipped"] = res.lochness.skipped
        if not res.lochness.by_sample.empty:
            out[f"{d}/lochness_by_sample"] = res.lochness.by_sample
        if not res.lochness.by_cluster.empty:
            out[f"{d}/lochness_by_cluster"] = res.lochness.by_cluster
    if res.modules is not None and not res.modules.empty:
        m = res.modules
        out[f"{d}/perturbation_effect_matrix"] = m.effect
        out[f"{d}/perturbation_de_mask"] = m.de_mask
        out[f"{d}/gene_programs"] = m.gene_programs
        out[f"{d}/perturbation_modules"] = m.perturbation_modules
        out[f"{d}/module_program_strength"] = m.module_program
        out[f"{d}/perturbation_program_effects"] = m.perturbation_program_effects
        if not m.program_activity_by_cluster.empty:
            out[f"{d}/program_activity_by_cluster"] = m.program_activity_by_cluster
        if not m.hubs.empty:
            out[f"{d}/tf_hubs"] = m.hubs
        if not m.tf_edges.empty:
            out[f"{d}/tf_edges"] = m.tf_edges
        if not m.module_connectivity.empty:
            out[f"{d}/module_connectivity"] = m.module_connectivity
    if res.protein is not None and not res.protein.empty:
        out[f"{d}/protein_effects"] = res.protein.table
        out[f"{d}/protein_effect_matrix"] = res.protein.matrix
    if res.concordance is not None:
        c = res.concordance
        for name, df in (("ps_protein_association", c.ps_protein), ("lochness_protein_association", c.lochness_protein), ("lochness_protein_summary", c.lochness_protein_summary),
                         ("program_protein_association_cells", c.program_protein_cells), ("program_protein_association_targets", c.program_protein_targets), ("perturbation_summary", c.summary)):
            if df is not None and not df.empty:
                out[f"{d}/{name}"] = df
    return {k: v for k, v in out.items() if v is not None}


# ---------------------------------------------------------------------------
# Cell states (Leiden) and perturbation x cluster enrichment
# ---------------------------------------------------------------------------


@dataclass
class CellStates:
    clustering: Any = None
    enrichment: Any = None

    def info(self) -> Dict[str, Any]:
        return {"clustering": dict(self.clustering.info) if self.clustering is not None else {"status": "disabled"},
                "enrichment": dict(self.enrichment.info) if self.enrichment is not None else {"status": "disabled"}}


def run_clustering(adata: ad.AnnData, cfg: Config) -> CellStates:
    """Leiden clustering (writes ``obs[analysis.clustering.key]``)."""
    from .clustering import compute_clustering

    return CellStates(clustering=compute_clustering(adata, cfg))


def run_enrichment(adata: ad.AnnData, cfg: Config, res: CellStates) -> CellStates:
    from .cluster_enrichment import compute_cluster_enrichment

    if cfg.analysis.clustering.enrichment.enabled:
        res.enrichment = compute_cluster_enrichment(adata, cfg, cfg.analysis.clustering.key)
    return res


def run_cell_states(adata: ad.AnnData, cfg: Config, lochness_summary: Optional[pd.DataFrame] = None) -> CellStates:
    """Clustering plus enrichment in one call (tests / notebooks); the pipeline runs them as separate stages."""
    from .cluster_enrichment import compare_with_lochness

    res = run_enrichment(adata, cfg, run_clustering(adata, cfg))
    if res.enrichment is not None and not res.enrichment.empty:
        res.enrichment.lochness_comparison = compare_with_lochness(res.enrichment, lochness_summary)
    return res


def cell_state_tables(res: CellStates) -> Dict[str, pd.DataFrame]:
    d = "cell_states"
    out: Dict[str, pd.DataFrame] = {}
    c = res.clustering
    if c is not None:
        out[f"{d}/cluster_summary"] = c.summary
        out[f"{d}/cluster_composition_by_class"] = c.by_class
        out[f"{d}/cluster_composition_by_target"] = c.by_target
        for col, tab in c.by_design.items():
            out[f"{d}/cluster_composition_by_{col}"] = tab
    e = res.enrichment
    if e is not None and not e.empty:
        from .cluster_enrichment import format_enrichment_table

        out[f"{d}/perturbation_cluster_enrichment"] = e.table
        out[f"{d}/enrichment"] = format_enrichment_table(e)
        out[f"{d}/enrichment_composition"] = e.composition.reset_index()
        out[f"{d}/enrichment_reference_composition"] = pd.DataFrame(e.reference_composition).T.rename_axis("control").reset_index()
        out[f"{d}/enrichment_effect_magnitude"] = e.effect_magnitude
        if not e.lochness_comparison.empty:
            out[f"{d}/cluster_enrichment_vs_lochness"] = e.lochness_comparison
    if e is not None and e.skipped is not None and len(e.skipped):
        out[f"{d}/cluster_enrichment_skipped"] = e.skipped
    return out
