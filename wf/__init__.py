from typing import List

from latch.resources.workflow import workflow
from latch.types.metadata import (
    LatchAuthor, LatchMetadata, LatchParameter, LatchRule
)

# from wf.task import registry_task, snap_task
from wf.task import (
    genes_task,
    motifs_task,
    make_adata,
    registry_task,
)
from wf.utils import Genome, Run

metadata = LatchMetadata(
    display_name="atx_snap",
    author=LatchAuthor(
        name="James McGann",
        email="jamesm@atlasxomics.com",
        github="github.com/atlasxomics",
    ),
    repository="https://github.com/atlasxomics/ATX_snap",
    license="MIT",
    parameters={
        "runs": LatchParameter(
            display_name="runs",
            description="List of runs to be analyzed; each run must contain a \
                run_id and fragments.tsv file; optional: condition. Note that \
                multiple Conditions must be separted by '_'\
                      (i.e., Female-control).",
            batch_table_column=True,
            samplesheet=True,
        ),
        "genome": LatchParameter(
            display_name="genome",
            description="Reference genome for runs.",
            batch_table_column=True,
        ),
        "resolution": LatchParameter(
            display_name="clustering resolution",
            description="Clustering resolution for Leiden algorithm; higher \
                values result in more clusters.",
            batch_table_column=True,
        ),
        "leiden_iters": LatchParameter(
            display_name="leiden iterations",
            description="Number of iterations for the leiden algorithm to \
                perform when assigning labels to clusters. 'Positive values \
                above 2 define the total number of iterations to perform, -1 \
                has the algorithm run until it reaches its optimal \
                clustering.' - SnapATAC2 docs.",
            batch_table_column=True,
            hidden=True,
        ),
        "min_cluster_size": LatchParameter(
            display_name="minimum cells per cluster",
            description="Minimum number of cells in a cluster.",
            batch_table_column=True,
            hidden=True,
        ),
        "min_tss": LatchParameter(
            display_name="minimum TSS",
            description="The minimum numeric transcription start site (TSS) \
                enrichment score required for a cell to pass filtering.",
            batch_table_column=True,
            hidden=True,
        ),
        "min_frags": LatchParameter(
            display_name="minimum fragments",
            description="The minimum number of mapped fragments required  \
                per cell to pass filtering.",
            batch_table_column=True,
            hidden=True,
        ),
        "tile_size": LatchParameter(
            display_name="tile size",
            description="The size of the tiles used for binning counts in the \
                tile matrix.",
            batch_table_column=True,
            hidden=True,
        ),
        "clustering_iters": LatchParameter(
            display_name="clustering iterations",
            description="Iterations performed when selecting variable \
                features for tile matrix. 'If greater than 1, this function \
                will perform iterative clustering and feature selection based \
                on variable features found using previous clustering results. \
                This is similar to the procedure implemented in ArchR... \
                Default value is 1, which means no iterative clustering is \
                 performed.'- SnapATAC2 docs",
            batch_table_column=True,
        ),
        "n_features": LatchParameter(
            display_name="number of features",
            description="Number of features to be selected as 'most \
                accessible' in tile matrix.",
            batch_table_column=True,
        ),
        "n_comps": LatchParameter(
            display_name="number of components",
            description="Number of components/dimensions to keep during \
                dimensionality reduction with `snap.tl.spectral`.",
            batch_table_column=True,
        ),
        "project_name": LatchParameter(
            display_name="project name",
            description="Name of output directory in snap_outs/",
            batch_table_column=True,
            rules=[
                LatchRule(
                    regex="^[^/].*", message="project name cannot start with a '/'"
                )
            ],
        ),
    },
)


@workflow(metadata)
def snap_workflow(
    runs: List[Run],
    genome: Genome,
    project_name: str,
    resolution: float = 1.0,
    leiden_iters: int = -1,
    n_comps: int = 30,
    min_cluster_size: int = 20,
    min_tss: float = 2.0,
    min_frags: int = 10,
    tile_size: int = 5000,
    n_features: int = 25000,
    clustering_iters: int = 1,
) -> None:
    """
    SnapATAC2 analysis for DBiT-seq experiments

    **ATX snap** is a [LatchBio](https://latch.bio/) Workflow for generating Python objects and data for analysis of epigenomic [DBiT-seq](https://www.nature.com/articles/s41586-022-05094-1) experiments.  Provided a fragment file from a single-cell ATAC-seq preprocessing and alignment tool (i.e., [Chromap](https://github.com/haowenz/chromap)) and spatial information, **ATX snap** performs routine processing with [SnapATAC2](https://kzhang.org/SnapATAC2/) and [scanpy](https://scanpy.readthedocs.io/en/stable/), returning files that can be easily input into custom scripts for more neuanced analysis without the need to perform intensive computation.

    The Workflow utilizes SnapATAC2 for AnnData object generation, quality control, dimensionality reduction/clustering, peak calling, and gene experssion analysis.  Scanpy is used to identify differential features.  [pychromVAR](https://github.com/pinellolab/pychromVAR) is used to create a motif deviation matrix.  The Workflow can take data from either a single tissue-sample analyzed via DBiT-seq or multiple tissue-samples; in ATX parlance, tissue-samples analyzed via DBIT-seq are termed 'Runs'.  All Runs input to **ATX snap** are merged into a single AnnData object for analysis.


    ## Inputs
    All input files for **ATX snap** must be on the LatchBio [file system](https://wiki.latch.bio/wiki/data/overview).  Each run in the workflow takes the following parameters.

    * [fragments.tsv.gz file](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments): A BED-like, tab-delimited file in which each row contains an ATAC-seq fragment.

    * [Spatial folder](https://docs.atlasxomics.com/projects/AtlasXbrowser/en/latest/SpatialFolder.html): A directory containing tissue images and experiment metadata.

    * Run ID: An identifier for the run.

    * Condition (_optional_):  An experimental Condition descriptor (i.e., 'control', 'diseased').

    Individual runs are batched in a combined AnnData object with the following global parameters.

    * project name: A name for the output folder.

    * genome: A reference genome to be used for alignment.

    * clustering resolution: 'A parameter value controlling the coarseness of the clustering. Higher values lead to more clusters'; maps to the resolution parameter of [`snapatac2.tl.leiden`](https://kzhang.org/SnapATAC2/api/_autosummary/snapatac2.tl.leiden.html).

    * number of components: Number of dimensions to keep after dimensionality reduction via SnapATAC2's spectral embedding algorithm; maps to the n_comps parameter of [`snapatac2.tl.spectral`](https://kzhang.org/SnapATAC2/api/_autosummary/snapatac2.tl.spectral.html).

    * number of features: n most accessible tiles to be used for differential analysis; maps to the n_features parameter of [`https://kzhang.org/SnapATAC2/api/_autosummary/snapatac2.pp.select_features.html`](https://kzhang.org/SnapATAC2/api/_autosummary/snapatac2.pp.select_features.html),

    * clustering iterations: Iterations performed when selecting variable features for tile matrix. 'If greater than 1, this function will perform iterative clustering and feature selection based on variable features found using previous clustering results. This is similar to the procedure implemented in ArchR... Default value is 1, which means no iterative clustering is performed.'  Maps to the n_iterations parameter of [snapatac2.tl.leiden](https://kzhang.org/SnapATAC2/api/_autosummary/snapatac2.tl.umap.html).

    The Workflow also takes a series of additional parameters that can be found under the 'Hidden Parameters' dropdown.

    * leiden iterations: Number of iterations for the leiden algorithm to perform when assigning labels to clusters. 'Positive values above 2 define the total number of iterations to perform, -1 has the algorithm run until it reaches its optimal clustering.' Maps to the n_iterations parameter of [`snapatac2.tl.leiden`](https://kzhang.org/SnapATAC2/api/_autosummary/snapatac2.tl.leiden.html).

    * minimum cells per cluster: Filter defining the minimum number of cell a cluster can contain for it to be kept; maps to the min_cluster_size parameter of [`snapatac2.tl.leiden`](https://kzhang.org/SnapATAC2/api/_autosummary/snapatac2.tl.leiden.html).

    * minimum TSS: Filter defining minimum TSS for a cell to be kept; maps to the min_tsse parameter of [`snapatac2.pp.filter_cells`](https://kzhang.org/SnapATAC2/api/_autosummary/snapatac2.pp.filter_cells.html).

    * minimum fragments: Filter defining minimum number of fragments for a cell to be kept; maps to the min_counts parameter of [`snapatac2.pp.filter_cells`](https://kzhang.org/SnapATAC2/api/_autosummary/snapatac2.pp.filter_cells.html).

    * tile size: 'The size of consecutive genomic regions used to record the counts.'; maps to the bin_size parameter of [`snapatac2.pp.add_tile_matrix`](https://kzhang.org/SnapATAC2/api/_autosummary/snapatac2.pp.add_tile_matrix.html).


    ## Running the workflow

    The **ATX snap** workflow can be found in the [Workflows](https://wiki.latch.bio/workflows/overview) module in your Latch workspace. For access to an ATX-collaborator workspace, please contact your AtlasXomics Support Scientist or email support@atlasxomics.com.  See [here](https://wiki.latch.bio/workflows/overview) for general instructions for running Workflows in Latch.

    1. Navigate to the **ATX snap** workflow in the Workflows module in your LatchBio workspace.  Ensure you are on the 'Parameters' tab of the workflow.

    2. To add Runs to the Project, select the '+ runs' icon.  Add values for the Run parameters described above; repeat for each Run in the Project.

    3. Scroll to the bottom of the page and input values for global project parameters.

    4. Click the 'Hidden Parameters' button and change the global parameters as needed.

    5. Click the 'Launch Workflow' button on the bottom-right of the parameters page.  This will automatically navigate you to the Executions tab of the workflow.

    6. From the Executions tab, you can view the status of the launched Workflow.  Once the Workflow has completed running, the status will change to 'Succeeded'; if the workflow has the status 'Failed', please contact an AtlasXomics Support Scientist.  You can click on the workflow execution to view a more granular workflow status and see output logs.

    7. Workflow outputs are loaded into the LatchBio [Data module](https://wiki.latch.bio/wiki/data/overview) in the `snap_outs` directory.


    ## Outputs

    Outputs from **ATX snap** are loaded into LatchBio [Data module](https://wiki.latch.bio/wiki/data/overview) in the `snap_outs` directory.
    * combined.h5ad: combined (all runs) AnnData object with .X as a tile matrix.
    * combined_ge.h5ad: combined (all runs) AnnData object with .X as a gene accessibility matrix; created with [`snapatac2.pp.make_gene_matrix`](https://kzhang.org/SnapATAC2/api/_autosummary/snapatac2.pp.make_gene_matrix.html).
    * combined_motifs.h5ad: combined (all runs) AnnData object with .X as a motif deviation matrix; created with [`pychromvar.compute_deviations`](https://pychromvar.readthedocs.io/en/latest/generated/pychromvar.compute_deviations.html#pychromvar.compute_deviations);
    * [group]_coverage/: folder containing bedgraph files for coverages grouped by cluster, sample, or condition. Browser tracks can be natively viewed in Latch Data.

    Figures generated by the Workflow can be found in the 'figures/' directory.
    * matrixplot_genes.pdf: Heatmap of top 5 genes per group (cluster, sample, or condition), ranked by the scanpy 'scores' value from `scanpy.tl.rank_genes_groups`; see [`scanpy.pl.rank_genes_groups_matrixplot`](https://scanpy.readthedocs.io/en/latest/api/generated/scanpy.pl.rank_genes_groups_matrixplot.html).
    * matrixplot_motifs.pdf: Heatmap of top 5 motifs per group (cluster, sample, or condition), ranked by the scanpy 'scores' value from `scanpy.tl.rank_genes_groups`; see [`scanpy.pl.rank_genes_groups_matrixplot`](https://scanpy.readthedocs.io/en/latest/api/generated/scanpy.pl.rank_genes_groups_matrixplot.html).
    * neighborhood_enrichemnt.pdf: See [`squidpy.pl.nhood_enrichment`](https://squidpy.readthedocs.io/en/stable/api/squidpy.pl.nhood_enrichment.html).
    * ripleys_L.pdf: See [`squidpy.pl.ripley`](https://squidpy.readthedocs.io/en/stable/api/squidpy.pl.ripley.html).
    * spatial_dim.pdf: Spatial scatter plot for each Run colored by cluster; see [`squidpy.pl.spatial_scatter`](https://squidpy.readthedocs.io/en/stable/api/squidpy.pl.spatial_scatter.html#squidpy.pl.spatial_scatter).
    * spatial_qc.pdf: Spatial scatter plot for each Run colored by number of fragments, log10(number of fragments), and TSS Enrichment score; see [`squidpy.pl.spatial_scatter`](https://squidpy.readthedocs.io/en/stable/api/squidpy.pl.spatial_scatter.html#squidpy.pl.spatial_scatter).
    * umap.pdf: UMAP embedding from SnapATAC2 colored by cluster; see [`snapatac2.pl.umap`](https://kzhang.org/SnapATAC2/api/_autosummary/snapatac2.pl.umap.html).

    Tables generated by the Workflow can be found in the 'tables/' directory.
    * differential_genes_per_[group]_counts.csv: Count of differential (pval_adj<0.05, log2fc_min>0.1) genes per group (cluster, sample, or condition).
    * differential_motifs_per_[group]_counts.csv: Count of differential (pval_adj<0.05, log2fc_min>0.1) motifs per group (cluster, sample, or condition).
    * gene_metadata.csv: Cell metadata containing gene metrics generated via [`scanpy.pp.calculate_qc_metrics`](https://scanpy.readthedocs.io/en/latest/generated/scanpy.pp.calculate_qc_metrics.html).
    * marker_peaks_per_[group].csv: Defining peaks per group (cluster, sample, or condition) with pval<0.05, log2fc_min>0.1; see [`scanpy.tl.filter_rank_genes_groups`](https://scanpy.readthedocs.io/en/latest/generated/scanpy.tl.filter_rank_genes_groups.html).
    * medians.csv: Median n_fragment, TSS Enrichment, and fraction-of-reads-in-peaks (frip) scores for each Run.
    * metadata.csv: Workflow input parameters.
    * ranked_genes_per_[group].csv: Gene output of [`scanpy.get.rank_genes_groups_df`](https://scanpy.readthedocs.io/en/latest/generated/scanpy.get.rank_genes_groups_df.html) with no filters.
    * ranked_motifs_per_[group].csv: Motif output of [`scanpy.get.rank_genes_groups_df`](https://scanpy.readthedocs.io/en/latest/generated/scanpy.get.rank_genes_groups_df.html) with no filters.


    ## Next Steps

    Analysis can be performed locally or in a LatchBio [Pod](https://wiki.latch.bio/wiki/pods/overview).  Output from **ATX snap** can also analyzed in a ATX [Latch Plot Template](https://wiki.latch.bio/plots/overview). For access to ATX-specific Pods or Plots, please contact your AtlasXomics Support Scientist.


    ## Support
    Questions? Comments?  Contact support@atlasxomics.com or post in AtlasXomics [Discord](https://discord.com/channels/1004748539827597413/1005222888384770108).

    """

    outdir, groups = make_adata(
        runs=runs,
        genome=genome,
        project_name=project_name,
        resolution=resolution,
        leiden_iters=leiden_iters,
        n_comps=n_comps,
        min_cluster_size=min_cluster_size,
        min_tss=min_tss,
        min_frags=min_frags,
        tile_size=tile_size,
        n_features=n_features,
        clustering_iters=clustering_iters,
    )

    outdir_ge = genes_task(
        runs=runs,
        outdir=outdir,
        project_name=project_name,
        genome=genome,
        groups=groups
    )

    outdir_motifs = motifs_task(
        runs=runs,
        outdir=outdir_ge,
        project_name=project_name,
        groups=groups,
        genome=genome,
    )

    uploaded_results = registry_task(runs=runs, results=outdir_motifs)

    return uploaded_results


if __name__ == "__main__":

    from latch.types import LatchDir, LatchFile

    # genes_task(
    #     runs=[Run(
    #         run_id="demo",
    #         fragments_file=LatchFile("latch://13502.account/atac_outs/ds_D01033_NG01681/outs/ds_D01033_NG01681_fragments.tsv.gz"),
    #         spatial_dir=LatchDir("latch://atx-illumina.mount/Images_spatial/D1033/spatial/"),
    #         condition="control",
    #     )],
    #     outdir="latch://13502.account/snap_outs/demo_001716",
    #     genome=Genome("hg38"),
    #     groups=["cluster"],
    #     project_name="develop_archrGenes",
    # )
    motifs_task(
        runs=[Run(
            run_id="demo",
            fragments_file=LatchFile("latch://13502.account/atx-archive-latch-1762285824.0858006/atac_outs/ds_D01033_NG01681/outs/ds_D01033_NG01681_fragments.tsv.gz"),
            spatial_dir=LatchDir("latch://atx-illumina.mount/Images_spatial/D1033/spatial/"),
            condition="demo",
        )],
        outdir=LatchDir("latch://13502.account/snap_outs/demo_002212"),
        project_name="demo_002212",
        groups=["cluster"],
        genome=Genome("hg38"),
    )
