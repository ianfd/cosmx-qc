from __future__ import annotations

import os
import re
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

from pathlib import Path
import tarfile
import gzip
import shutil


import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from anndata import AnnData
from scipy.sparse import csr_matrix
from skimage.transform import estimate_transform
from spatialdata import SpatialData
from spatialdata.models import PointsModel, TableModel
from spatialdata.transformations.transformations import Affine, Identity
import squidpy as sq

if TYPE_CHECKING:
    from dask.dataframe import DataFrame as DaskDataFrame


class CosmxKeys:
    COUNTS_SUFFIX = "exprMat_file.csv"
    TRANSCRIPTS_SUFFIX = "tx_file.csv"
    METADATA_SUFFIX = "metadata_file.csv"
    FOV_SUFFIX = "fov_positions_file.csv"
    FOV = "fov"
    CELL_ID = "cell_ID"
    INSTANCE_KEY = "cell"
    REGION_KEY = "fov_labels"
    X_GLOBAL_CELL = "CenterX_global_px"
    Y_GLOBAL_CELL = "CenterY_global_px"
    X_LOCAL_CELL = "CenterX_local_px"
    Y_LOCAL_CELL = "CenterY_local_px"
    X_LOCAL_TRANSCRIPT = "x_local_px"
    Y_LOCAL_TRANSCRIPT = "y_local_px"
    TARGET_OF_TRANSCRIPT = "target"


def cosmx(
    path: str | Path,
    dataset_id: str | None = None,
    transcripts: bool = True,
) -> SpatialData:
    """Read *Cosmx Nanostring* data (no image/label loading).
        From the old spatialdata-io function
    This function reads the following files:

        - ``<dataset_id>_exprMat_file.csv``: Counts matrix.
        - ``<dataset_id>_metadata_file.csv``: Metadata file.
        - ``<dataset_id>_fov_positions_file.csv``: Field of view file.

    Parameters
    ----------
    path
        Path to the root directory containing *Nanostring* files.
    dataset_id
        Name of the dataset. Inferred from counts filename if not provided.
    transcripts
        Whether to also read in transcripts information.

    Returns
    -------
    :class:`spatialdata.SpatialData`
    """
    path = Path(path)

    if dataset_id is None:
        counts_files = [f for f in os.listdir(path) if f.endswith(CosmxKeys.COUNTS_SUFFIX)]
        if len(counts_files) == 1:
            found = re.match(rf"(.*)_{CosmxKeys.COUNTS_SUFFIX}", counts_files[0])
            if found:
                dataset_id = found.group(1)
    if dataset_id is None:
        raise ValueError("Could not infer `dataset_id` from the name of the counts file. Please specify it manually.")

    counts_file = path / f"{dataset_id}_{CosmxKeys.COUNTS_SUFFIX}"
    meta_file = path / f"{dataset_id}_{CosmxKeys.METADATA_SUFFIX}"
    fov_file = path / f"{dataset_id}_{CosmxKeys.FOV_SUFFIX}"
    transcripts_file = path / f"{dataset_id}_{CosmxKeys.TRANSCRIPTS_SUFFIX}" if transcripts else None

    for f, label in [(counts_file, "Counts"), (meta_file, "Metadata"), (fov_file, "FOV positions")]:
        if not f.exists():
            raise FileNotFoundError(f"{label} file not found: {f}")
    if transcripts_file is not None and not transcripts_file.exists():
        raise FileNotFoundError(f"Transcripts file not found: {transcripts_file}")

    counts_df = pd.read_csv(counts_file, header=0)
    obs_df = pd.read_csv(meta_file, header=0)

    cell_map = obs_df[[CosmxKeys.CELL_ID, CosmxKeys.FOV, CosmxKeys.INSTANCE_KEY]].drop_duplicates()
    counts_df = counts_df.merge(cell_map, on=[CosmxKeys.CELL_ID, CosmxKeys.FOV], how="left")
    counts_df = counts_df.set_index(CosmxKeys.INSTANCE_KEY).drop(columns=[CosmxKeys.CELL_ID, CosmxKeys.FOV])

    obs_df.drop(columns=[c for c in obs_df.columns if c.lower() == "cell_id"], inplace=True)
    obs_df[CosmxKeys.FOV] = obs_df[CosmxKeys.FOV].astype(str)
    obs_df = obs_df.set_index(CosmxKeys.INSTANCE_KEY)

    common_index = obs_df.index.intersection(counts_df.index)
    obs = obs_df.loc[common_index].copy()
    obs.index = obs.index.astype(str)
    obs[CosmxKeys.INSTANCE_KEY] = obs.index.values
    obs[CosmxKeys.FOV] = pd.Categorical(obs[CosmxKeys.FOV])

    adata = AnnData(
        csr_matrix(counts_df.loc[common_index].values),
        dtype=counts_df.values.dtype,
        obs=obs,
    )
    adata.var_names = counts_df.columns

    table = TableModel.parse(
        adata,
        instance_key=CosmxKeys.INSTANCE_KEY,
    )

    fovs = sorted(obs[CosmxKeys.FOV].cat.categories.tolist())
    affine_transforms_to_global: dict[str, Affine | Identity] = {}

    for fov in fovs:
        idx = table.obs[CosmxKeys.FOV] == fov
        loc = table[idx].obs[[CosmxKeys.X_LOCAL_CELL, CosmxKeys.Y_LOCAL_CELL]].values
        glob = table[idx].obs[[CosmxKeys.X_GLOBAL_CELL, CosmxKeys.Y_GLOBAL_CELL]].values
        tform = estimate_transform("affine", src=loc, dst=glob)
        if not tform:
            tform = estimate_transform("similarity", src=loc, dst=glob)
        if not tform:
            print(f"Warning: could not estimate transform for FOV {fov}, using identity.")
            affine_transforms_to_global[fov] = Identity()
        else:
            affine_transforms_to_global[fov] = Affine(
                tform.params, input_axes=("x", "y"), output_axes=("x", "y"),
            )

    table.obsm["global"] = table.obs[[CosmxKeys.X_GLOBAL_CELL, CosmxKeys.Y_GLOBAL_CELL]].to_numpy()
    table.obsm["spatial"] = table.obs[[CosmxKeys.X_LOCAL_CELL, CosmxKeys.Y_LOCAL_CELL]].to_numpy()
    table.obs.drop(
        columns=[CosmxKeys.X_LOCAL_CELL, CosmxKeys.Y_LOCAL_CELL, CosmxKeys.X_GLOBAL_CELL, CosmxKeys.Y_GLOBAL_CELL],
        inplace=True,
    )

    points: dict[str, DaskDataFrame] = {}
    if transcripts_file is not None:
        with tempfile.TemporaryDirectory() as tmpdir:
            print("converting transcripts .csv → .parquet ... ", end="", flush=True)
            pq_path = Path(tmpdir) / "transcripts.parquet"
            pd.read_csv(transcripts_file, header=0).to_parquet(pq_path)
            print("done")

            ptable = pq.read_table(pq_path)
            for fov in fovs:
                aff = affine_transforms_to_global[fov]
                sub = ptable.filter(pa.compute.equal(ptable.column(CosmxKeys.FOV), int(fov))).to_pandas()
                if len(sub) == 0:
                    continue
                sub[CosmxKeys.INSTANCE_KEY] = sub[CosmxKeys.INSTANCE_KEY].astype("category")
                sub.rename(columns={"z": "z_raw"}, inplace=True)
                sub.drop(columns=[c for c in sub.columns if c.lower() == "cell_id"], inplace=True, errors="ignore")
                points[f"{fov}_points"] = PointsModel.parse(
                    sub,
                    coordinates={"x": CosmxKeys.X_LOCAL_TRANSCRIPT, "y": CosmxKeys.Y_LOCAL_TRANSCRIPT},
                    feature_key=CosmxKeys.TARGET_OF_TRANSCRIPT,
                    instance_key=CosmxKeys.INSTANCE_KEY,
                    transformations={
                        fov: Identity(),
                        "global": aff,
                        "global_only_labels": aff,
                    },
                )

    return SpatialData(points=points, tables={"table": table})


def cosmx_simple(
    path: str | Path,
    sample_name: str,
    dataset_id: str | None = None,
) -> AnnData:

    path = Path(path)

    if dataset_id is None:
        counts_files = [f for f in os.listdir(path) if f.endswith(CosmxKeys.COUNTS_SUFFIX)]
        if len(counts_files) == 1:
            found = re.match(rf"(.*)_{CosmxKeys.COUNTS_SUFFIX}", counts_files[0])
            if found:
                dataset_id = found.group(1)
    if dataset_id is None:
        raise ValueError("Could not infer `dataset_id` from the name of the counts file. Please specify it manually.")

    adata = sq.read.nanostring(
        path=path,
        counts_file=f"{dataset_id}_{CosmxKeys.COUNTS_SUFFIX}",
        meta_file=f"{dataset_id}_{CosmxKeys.METADATA_SUFFIX}",
        fov_file=f"{dataset_id}_{CosmxKeys.FOV_SUFFIX}"
    )

    adata.obs.index = adata.obs['cell'].astype(str)
    adata.obs['sample'] = sample_name

    return adata


def extract_files(input_file: str, work_dir: str | Path) -> str:
    sample_path = Path(input_file)   
    work_path = Path(work_dir)

    temp_dir = work_path / f"{sample_path.name.replace('.tar.gz', '')}_tmp"
    tar_extract_dir = temp_dir / "tar_contents"
    tar_extract_dir.mkdir(parents=True, exist_ok=True)

    with tarfile.open(sample_path, "r:gz") as tar:
        tar.extractall(path=tar_extract_dir)

    csv_gz_files = list(tar_extract_dir.rglob("*.csv.gz"))
    for gz_file in csv_gz_files:
        out_file = (temp_dir / gz_file.relative_to(tar_extract_dir)).with_suffix("")
        out_file.parent.mkdir(parents=True, exist_ok=True)

        with gzip.open(gz_file, "rb") as f_in, open(out_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    return str(temp_dir)

def prep_fov_file(sample_folder: str):
    target_file = ""

    for f in os.listdir(sample_folder):
        print(f)
        if f.endswith("_fov_positions_file.csv"):
            target_file = f
    
    if target_file == "":
        return
    
    smpl_full_path = str(sample_folder) + "/" + target_file

    d = pd.read_csv(smpl_full_path)

    d = d.rename({"FOV": "fov"}, axis="columns")
    
    d.to_csv(smpl_full_path)

def read_full_sample(sample_name: str, data: str, work_dir: str | Path):
    sample_folder = extract_files(input_file=data, work_dir=work_dir)
    prep_fov_file(sample_folder)

    spatial_data = cosmx_simple(
        path=sample_folder,
        sample_name=sample_name
    )

    
    
    return spatial_data    




