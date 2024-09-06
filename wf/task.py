import logging
import os
import pandas as pd
import pychromvar as pc
import snapatac2 as snap

from typing import List
from pyjaspar import jaspardb

from latch.resources.tasks import large_task
from latch.types import LatchDir

import wf.preprocessing as pp
import wf.spatial as sp
import wf.features as ft

from wf.utils import Genome, Run


logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s",
    level=logging.INFO
)


@large_task
def snap_task(
    runs: List[Run],
    genome: Genome,
    min_tss: float,
    min_frags: int,
    tile_size: int,
    project_name: str
) -> LatchDir:

    samples = [run.run_id for run in runs]
    groups = ["cluster", "sample", "condition"]
    genome = genome.value  # Convert to str

    out_dir = f"/root/{project_name}"
    os.makedirs(out_dir, exist_ok=True)

    adatas = pp.make_anndatas(runs, genome, min_frags=min_frags)
    adatas = pp.filter_adatas(adatas, min_tss=min_tss)
    adatas = pp.add_tilematrix(adatas, tile_size=tile_size, n_features=25000)

    # Combine seperate runs into single analysis object
    adata = pp.combine_anndata(adatas, samples, filename="combined")
    snap.pp.select_features(adata, n_features=25000)
    adata = pp.add_clusters(adata)
    adata = sp.add_spatial(adata)

    adata_gene = ft.make_geneadata(adata, genome)
    


    adata.write(f"{out_dir}/combined.h5ad")

    return LatchDir(
        out_dir, f"latch:///snap_outs/{project_name}"
    )
