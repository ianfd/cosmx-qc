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
        "per_sample_overrides": LatchParameter(
            display_name="Per-Sample QC Overrides",
            description=(
                "Optional dict mapping sample name -> dict of QC param overrides. "
                'E.g. {"Slide1": {"min_counts": 50, "max_counts": 800}}. '
                "Only the params you specify are overridden; the rest use defaults."
            ),
        ),
        "output_dir": LatchParameter(
            display_name="Output Directory",
            description="Base Latch path for filtered H5ADs and post-QC statistics.",
            batch_table_column=True,
        ),
    },
)


@small_task
def prep_qc_args(
    sample_h5ads: List[LatchFile],
    min_counts: int,
    max_counts: int,
    min_genes: int,
    max_neg_probe_ratio: float,
    min_cells_per_gene: int,
    per_sample_overrides: Dict[str, Dict[str, float]],
    output_dir: LatchDir,
) -> List[QCInput]:
    import scanpy as sc

    defaults = {
        "min_counts": min_counts,
        "max_counts": max_counts,
        "min_genes": min_genes,
        "max_neg_probe_ratio": max_neg_probe_ratio,
        "min_cells_per_gene": min_cells_per_gene,
    }

    inputs: List[QCInput] = []
    for h5ad_file in sample_h5ads:
        adata = sc.read_h5ad(h5ad_file.local_path, backed="r")
        sample_name = str(adata.obs["sample"].iloc[0])
        adata.file.close()

        params = {**defaults}
        if sample_name in per_sample_overrides:
            for k, v in per_sample_overrides[sample_name].items():
                if k not in defaults:
                    raise ValueError(
                        f"Unknown override key '{k}' for sample '{sample_name}'. "
                        f"Valid keys: {list(defaults.keys())}"
                    )
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
    per_sample_overrides: Dict[str, Dict[str, float]] = {},
    output_dir: LatchDir = LatchDir("latch://40726.account/cosmx-test/out-dir/qc/"),
) -> List[LatchOutputDir]:
    """
    ## CosMx QC Filtering (Multi-Sample)

    Takes a list of H5ADs from the conversion workflow, applies cell/gene
    filters (with optional per-sample overrides), and outputs filtered
    H5ADs with post-QC statistics.

    ### Filtering steps
    1. Remove cells from failed FOVs (`qcFlagsFOV != "Pass"`)
    2. Filter cells by min/max total counts
    3. Filter cells by min genes detected
    4. Filter cells by negative probe ratio
    5. Filter genes by min cells expressing

    ### Per-sample overrides
    Supply a dict like `{"Slide1": {"min_counts": 50, "max_counts": 800}}`
    to override specific thresholds for individual samples. Unspecified
    params fall back to the top-level defaults.

    ### Outputs (per sample)
    - `{sample_name}.h5ad` — filtered AnnData with raw counts
    - `stats/` — post-filter statistics CSVs
    """
    qc_inputs = prep_qc_args(
        sample_h5ads=sample_h5ads,
        min_counts=min_counts,
        max_counts=max_counts,
        min_genes=min_genes,
        max_neg_probe_ratio=max_neg_probe_ratio,
        min_cells_per_gene=min_cells_per_gene,
        per_sample_overrides=per_sample_overrides,
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
        "per_sample_overrides": {
            "GSE282193_Slide1": {"min_counts": 50, "max_counts": 400},
        },
        "output_dir": LatchDir("latch://40726.account/cosmx-test/out-dir/qc/"),
    },
)