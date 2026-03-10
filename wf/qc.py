from os import PathLike

from latch.resources.tasks import small_task, medium_task
from latch.types.file import LatchFile
from latch.types.directory import LatchDir, LatchOutputDir
from pathlib import Path

from sup import loading_utils


## PARAM ideas:
# 1k RNA => min 20 ;; 6k RNA => min 50

@medium_task
def cosmx_qc_task(
    sample_h5ad: LatchFile,
    sample_name: str,
    min_counts: int,
    max_counts: int,
    min_genes: int,
    max_neg_probe_ratio: float,
    min_cells_per_gene: int,
    output_dir: LatchDir,
) -> LatchOutputDir:
    import scanpy as sc
    import pandas as pd
    import json

    local_out = Path("/root/cosmx_qc_output")
    local_out.mkdir(exist_ok=True)
    path_out = Path(local_out / "stats")
    path_out.mkdir(exist_ok=True)

    safe_name = sample_name # manipulate top level, not here!

    p = sample_h5ad.local_path
    adata = sc.read_h5ad(p)

    if "qcFlagsFOV" in adata.obs.columns:
        adata = adata[adata.obs["qcFlagsFOV"] == "Pass", :].copy()

    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_cells(adata, max_counts=max_counts)
    sc.pp.filter_cells(adata, min_genes=min_genes)

    neg_prob_filter = adata.obs["pct_counts_negprobes"] > (max_neg_probe_ratio * 100)
    adata = adata[~neg_prob_filter, :].copy()

    sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)

    _save_cell_stats(adata, path_out, safe_name, "postfilter")
    _save_fov_stats(adata, path_out, safe_name, "postfilter")
    _save_protein_stats(adata, path_out, safe_name, "postfilter")
    _save_summary(adata, path_out, safe_name, "postfilter")

    adata.write_h5ad(local_out / f"{safe_name}.h5ad")

    return LatchDir(str(local_out), output_dir)


def _save_cell_stats(
    adata, out_dir: Path, sample_name: str, stage: str
) -> None:
    import pandas as pd
    import numpy as np

    cols = ["total_counts", "n_genes_by_counts", "pct_counts_negprobes"]
    cols = [c for c in cols if c in adata.obs.columns]

    stats = adata.obs[cols].describe().T
    stats["median"] = adata.obs[cols].median()
    stats.to_csv(out_dir / f"{sample_name}_cell_stats_{stage}.csv")


def _save_fov_stats(
    adata, out_dir: Path, sample_name: str, stage: str
) -> None:
    import pandas as pd
    import numpy as np

    if "fov" not in adata.obs.columns:
        return

    agg_dict = {
        "n_cells": ("total_counts", "count"),
        "median_counts": ("total_counts", "median"),
        "mean_counts": ("total_counts", "mean"),
        "median_genes": ("n_genes_by_counts", "median"),
        "total_counts": ("total_counts", "sum"),
    }
    if "pct_counts_negprobes" in adata.obs.columns:
        agg_dict["median_pct_negprobes"] = ("pct_counts_negprobes", "median")

    fov_stats = adata.obs.groupby("fov").agg(**agg_dict)
    fov_stats.to_csv(out_dir / f"{sample_name}_fov_stats_{stage}.csv")
    

def _save_protein_stats(
    adata, out_dir: Path, sample_name: str, stage: str
) -> None:
    import pandas as pd

    protein_cols = [
        c for c in adata.obs.columns
        if c.startswith("Mean.") or c.startswith("Max.")
    ]
    if not protein_cols:
        return

    stats = adata.obs[protein_cols].describe().T
    stats["median"] = adata.obs[protein_cols].median()
    stats.to_csv(out_dir / f"{sample_name}_protein_stats_{stage}.csv")


def _save_summary(
    adata, out_dir: Path, sample_name: str, stage: str
) -> None:
    import pandas as pd
    import numpy as np

    summary = {
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "median_total_counts": float(np.median(adata.obs["total_counts"])),
        "mean_total_counts": float(np.mean(adata.obs["total_counts"])),
        "median_genes_per_cell": float(np.median(adata.obs["n_genes_by_counts"])),
    }

    if "pct_counts_negprobes" in adata.obs.columns:
        summary["median_pct_negprobes"] = float(
            np.median(adata.obs["pct_counts_negprobes"])
        )

    if "fov" in adata.obs.columns:
        summary["n_fovs"] = int(adata.obs["fov"].nunique())

    pd.DataFrame([summary]).to_csv(
        out_dir / f"{sample_name}_summary_{stage}.csv", index=False
    )