"""Stage E: perturbation-effect analyses on the processed multimodal object.

``run_perturbation_effects`` runs, according to ``analysis.perturbation_effects``:
PS (``ps_score``), lochNESS (``lochness``), gene programs / perturbation
modules (``modules``), protein effects (``protein_effects``) and the
RNA-protein concordance + integrated summary (``concordance``). ``attach``
writes the results into the AnnData (obs / obsm / uns) and ``tables`` returns
every result table under ``perturbation_effects/``.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, Dict, Optional

import anndata as ad
import pandas as pd

from ..config import Config
from .concordance import ConcordanceResults, compute_concordance
from .lochness import LochnessResults, compute_lochness
from .modules import ModulesResults, compute_modules
from .protein_effects import ProteinEffectsResults, compute_protein_effects
from .ps_score import PSResults, compute_ps

logger = logging.getLogger(__name__)


@dataclass
class PerturbationEffects:
    ps: Optional[PSResults] = None
    lochness: Optional[LochnessResults] = None
    modules: Optional[ModulesResults] = None
    protein: Optional[ProteinEffectsResults] = None
    concordance: Optional[ConcordanceResults] = None
    notes: Optional[list] = None

    def info(self) -> Dict[str, Any]:
        out: Dict[str, Any] = {"enabled": True}
        for name in ("ps", "lochness", "modules", "protein", "concordance"):
            r = getattr(self, name)
            out[name] = dict(r.info) if r is not None else {"status": "disabled"}
        return out


def run_perturbation_effects(adata: ad.AnnData, cfg: Config) -> PerturbationEffects:
    pe = cfg.analysis.perturbation_effects
    res = PerturbationEffects(notes=[])
    if pe.ps.enabled:
        res.ps = compute_ps(adata, cfg)
    if pe.lochness.enabled:
        res.lochness = compute_lochness(adata, cfg)
    if pe.modules.enabled:
        res.modules = compute_modules(adata, cfg)
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
    if res.ps is not None and not res.ps.scores.empty:
        adata.obsm["ps_scores"] = res.ps.scores.astype("float32")
        adata.obsm["ps_scores_raw"] = res.ps.raw.astype("float32")
        adata.obs["ps_score"] = res.ps.own.to_numpy(dtype="float32")
        adata.obs["ps_quadrant"] = pd.Categorical(res.ps.own_quadrant.astype(str))
        adata.uns["ps_signatures"] = _h5(res.ps.signatures)
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
    if res.ps is not None:
        out[f"{d}/ps_targets"] = res.ps.summary
        out[f"{d}/ps_skipped"] = res.ps.skipped
        out[f"{d}/ps_signatures"] = res.ps.signatures
    if res.lochness is not None:
        out[f"{d}/lochness_targets"] = res.lochness.summary
        out[f"{d}/lochness_skipped"] = res.lochness.skipped
        if not res.lochness.by_sample.empty:
            out[f"{d}/lochness_by_sample"] = res.lochness.by_sample
    if res.modules is not None and not res.modules.empty:
        out[f"{d}/perturbation_effect_matrix"] = res.modules.effect
        out[f"{d}/perturbation_de_mask"] = res.modules.de_mask
        out[f"{d}/gene_programs"] = res.modules.gene_programs
        out[f"{d}/perturbation_modules"] = res.modules.perturbation_modules
        out[f"{d}/module_program_strength"] = res.modules.module_program
        out[f"{d}/perturbation_program_effects"] = res.modules.perturbation_program_effects
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
# Stage F: cell states (Leiden) and perturbation x cluster enrichment
# ---------------------------------------------------------------------------


@dataclass
class CellStates:
    clustering: Any = None
    enrichment: Any = None

    def info(self) -> Dict[str, Any]:
        return {"clustering": dict(self.clustering.info) if self.clustering is not None else {"status": "disabled"},
                "enrichment": dict(self.enrichment.info) if self.enrichment is not None else {"status": "disabled"}}


def run_cell_states(adata: ad.AnnData, cfg: Config, lochness_summary: Optional[pd.DataFrame] = None) -> CellStates:
    """Leiden clustering (writes ``obs[analysis.clustering.key]``) and, if enabled, enrichment."""
    from .cluster_enrichment import compute_cluster_enrichment
    from .clustering import compute_clustering

    cl = cfg.analysis.clustering
    res = CellStates(clustering=compute_clustering(adata, cfg))
    if cl.enrichment.enabled:
        res.enrichment = compute_cluster_enrichment(adata, cfg, cl.key, lochness_summary)
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
        out[f"{d}/perturbation_cluster_enrichment"] = e.table
        out[f"{d}/cluster_enrichment_vs_lochness"] = e.lochness_comparison
    if e is not None and e.skipped is not None and len(e.skipped):
        out[f"{d}/cluster_enrichment_skipped"] = e.skipped
    return out
