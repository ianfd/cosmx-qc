from typing import Dict, List, Optional

from latch.resources.workflow import workflow
from latch.resources.tasks import small_task
from latch.resources.map_tasks import map_task
from latch.resources.launch_plan import LaunchPlan
from latch.types.file import LatchFile
from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.metadata import LatchAuthor, LatchMetadata, LatchParameter

from .qc import QCInput, cosmx_qc_task


metadata = LatchMetadata(
    display_name="CosMx QC Filtering (Multi-Sample)",
    author=LatchAuthor(name="Ian"),
    parameters={
        "sample_h5ads": LatchParameter(
            display_name="Sample H5AD files",
            description="List of H5AD files from the conversion workflow. Sample name is read from the 'sample' column.",
            batch_table_column=True,
        ),
        "min_counts": LatchParameter(
            display_name="Min Counts per Cell (default)",
            description="Default minimum total transcript counts. Suggested: 20 (1K panel), 50 (6K panel).",
        ),
        "max_counts": LatchParameter(
            display_name="Max Counts per Cell (default)",
            description="Default maximum total transcript counts.",
        ),
        "min_genes": LatchParameter(
            display_name="Min Genes per Cell (default)",
            description="Default minimum genes detected per cell.",
        ),
        "max_neg_probe_ratio": LatchParameter(
            display_name="Max Negative Probe Ratio (default)",
            description="Default max fraction (0-1) of negative probe counts. E.g. 0.02 = 2%.",
        ),
        "min_cells_per_gene": LatchParameter(
            display_name="Min Cells per Gene (default)",
            description="Default minimum cells expressing a gene to keep it.",
        ),
        "per_sample_overrides_csv": LatchParameter(
            display_name="Per-Sample QC Overrides (CSV)",
            description=(
                "Optional CSV with sample name as the index column and QC param "
                "columns (min_counts, max_counts, min_genes, max_neg_probe_ratio, "
                "min_cells_per_gene). Only include rows/columns you want to override; "
                "missing values fall back to the top-level defaults."
            ),
        ),
        "output_dir": LatchParameter(
            display_name="Output Directory",
            description="Base Latch path for filtered H5ADs and post-QC statistics.",
            batch_table_column=True,
        ),
    },
)


VALID_OVERRIDE_KEYS = {
    "min_counts", "max_counts", "min_genes",
    "max_neg_probe_ratio", "min_cells_per_gene",
}


@small_task
def prep_qc_args(
    sample_h5ads: List[LatchFile],
    min_counts: int,
    max_counts: int,
    min_genes: int,
    max_neg_probe_ratio: float,
    min_cells_per_gene: int,
    per_sample_overrides_csv: Optional[LatchFile],
    output_dir: LatchDir,
) -> List[QCInput]:
    import pandas as pd
    from pathlib import Path

    defaults = {
        "min_counts": min_counts,
        "max_counts": max_counts,
        "min_genes": min_genes,
        "max_neg_probe_ratio": max_neg_probe_ratio,
        "min_cells_per_gene": min_cells_per_gene,
    }

    overrides: Dict[str, Dict[str, float]] = {}
    if per_sample_overrides_csv is not None:
        df = pd.read_csv(per_sample_overrides_csv.local_path, index_col=0)
        bad_cols = set(df.columns) - VALID_OVERRIDE_KEYS
        if bad_cols:
            raise ValueError(
                f"Unknown columns in overrides CSV: {bad_cols}. "
                f"Valid columns: {VALID_OVERRIDE_KEYS}"
            )
        for sample_name, row in df.iterrows():
            overrides[str(sample_name)] = {
                k: v for k, v in row.items() if pd.notna(v)
            }

    inputs: List[QCInput] = []
    for h5ad_file in sample_h5ads:
        sample_name = Path(h5ad_file.remote_path).stem

        params = {**defaults}
        if sample_name in overrides:
            for k, v in overrides[sample_name].items():
                params[k] = type(defaults[k])(v)

        inputs.append(
            QCInput(
                sample_name=sample_name,
                sample_h5ad=h5ad_file,
                min_counts=int(params["min_counts"]),
                max_counts=int(params["max_counts"]),
                min_genes=int(params["min_genes"]),
                max_neg_probe_ratio=float(params["max_neg_probe_ratio"]),
                min_cells_per_gene=int(params["min_cells_per_gene"]),
                output_dir=output_dir,
            )
        )

    return inputs


@workflow(metadata)
def cosmx_qc(
    sample_h5ads: List[LatchFile],
    min_counts: int = 20,
    max_counts: int = 5000,
    min_genes: int = 5,
    max_neg_probe_ratio: float = 0.02,
    min_cells_per_gene: int = 10,
    per_sample_overrides_csv: Optional[LatchFile] = None,
    output_dir: LatchDir = LatchDir("latch://40726.account/cosmx-test/out-dir/qc/"),
) -> List[LatchOutputDir]:
    qc_inputs = prep_qc_args(
        sample_h5ads=sample_h5ads,
        min_counts=min_counts,
        max_counts=max_counts,
        min_genes=min_genes,
        max_neg_probe_ratio=max_neg_probe_ratio,
        min_cells_per_gene=min_cells_per_gene,
        per_sample_overrides_csv=per_sample_overrides_csv,
        output_dir=output_dir,
    )

    return map_task(cosmx_qc_task)(input=qc_inputs)


LaunchPlan(
    cosmx_qc,
    "Test Data",
    {
        "sample_h5ads": [
            LatchFile("latch://40726.account/cosmx-test/out-dir/conv/GSE282193_Slide1/GSE282193_Slide1.h5ad"),
        ],
        "min_counts": 20,
        "max_counts": 5000,
        "min_genes": 5,
        "max_neg_probe_ratio": 0.02,
        "min_cells_per_gene": 10,
        "per_sample_overrides_csv": None,
        "output_dir": LatchDir("latch://40726.account/cosmx-test/out-dir/qc/"),
    },
)