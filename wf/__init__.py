from latch.resources.workflow import workflow
from latch.types.file import LatchFile
from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.metadata import (
    LatchAuthor,
    LatchMetadata,
    LatchParameter,
)
from latch.resources.launch_plan import LaunchPlan

from .qc import cosmx_qc_task


metadata = LatchMetadata(
    display_name="CosMx QC Filtering",
    author=LatchAuthor(name="Ian"),
    parameters={
        "sample_h5ad": LatchParameter(
            display_name="Sample H5AD",
            description="H5AD file from the conversion workflow (raw counts, QC metrics precomputed).",
            batch_table_column=True,
        ),
        "sample_name": LatchParameter(
            display_name="Sample Name",
            description="Name used for output file naming.",
            batch_table_column=True,
        ),
        "min_counts": LatchParameter(
            display_name="Min Counts per Cell",
            description="Minimum total transcript counts to keep a cell. Suggested: 20 for 1K panel, 50 for 6K panel.",
        ),
        "max_counts": LatchParameter(
            display_name="Max Counts per Cell",
            description="Maximum total transcript counts (removes doublets/artifacts).",
        ),
        "min_genes": LatchParameter(
            display_name="Min Genes per Cell",
            description="Minimum number of genes detected to keep a cell.",
        ),
        "max_neg_probe_ratio": LatchParameter(
            display_name="Max Negative Probe Ratio",
            description="Maximum fraction (0-1) of negative probe counts to total counts. E.g. 0.2 = 20%.",
        ),
        "min_cells_per_gene": LatchParameter(
            display_name="Min Cells per Gene",
            description="Minimum number of cells expressing a gene to keep it.",
        ),
        "output_dir": LatchParameter(
            display_name="Output Directory",
            description="Latch path for filtered H5AD and post-QC statistics.",
            batch_table_column=True,
        ),
    },
)


@workflow(metadata)
def cosmx_qc(
    sample_h5ad: LatchFile,
    sample_name: str,
    min_counts: int = 20,
    max_counts: int = 50000,
    min_genes: int = 5,
    max_neg_probe_ratio: float = 0.2,
    min_cells_per_gene: int = 10,
    output_dir: LatchOutputDir = LatchDir("latch://40726.account/cosmx-test/out-dir"),
) -> LatchOutputDir:
    """
    ## CosMx QC Filtering

    Takes an H5AD from the conversion workflow, applies cell/gene
    filters, and outputs a filtered H5AD with post-QC statistics.

    ### Filtering steps
    1. Remove cells from failed FOVs (`qcFlagsFOV != "Pass"`)
    2. Filter cells by min/max total counts
    3. Filter cells by min genes detected
    4. Filter cells by negative probe ratio
    5. Filter genes by min cells expressing

    ### Outputs
    - `{sample_name}.h5ad` — filtered AnnData with raw counts
    - `stats/` — post-filter statistics CSVs (cell stats, FOV stats, protein stats, summary)
    """
    return cosmx_qc_task(
        sample_h5ad=sample_h5ad,
        sample_name=sample_name,
        min_counts=min_counts,
        max_counts=max_counts,
        min_genes=min_genes,
        max_neg_probe_ratio=max_neg_probe_ratio,
        min_cells_per_gene=min_cells_per_gene,
        output_dir=output_dir,
    )


"""
Add test data with a LaunchPlan. Provide default values in a dictionary with
the parameter names as the keys. These default values will be available under
the 'Test Data' dropdown at console.latch.bio.
"""
LaunchPlan(
    cosmx_qc,
    "Test Data",
    {
        "expr_mat": LatchFile("latch://40726.account/cosmx-test/GSE282193_Slide1.tar.gz"),
        "sample_name": "GSE282193_Slide1",
        "min_counts": 50,
        "max_counts": 400,
        "min_genes": 10,
        "max_neg_probe_ratio": 0.01,
        "output_directory": LatchDir("latch://40726.account/cosmx-test/out-dir/GSE282193_Slide1"),
    }
)

