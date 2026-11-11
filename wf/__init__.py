from typing import List

from atx_common import Genome
from latch.resources.workflow import workflow
from latch.types import LatchDir
from latch.types.metadata import (
    LatchAuthor, LatchMetadata, LatchParameter, LatchRule
)

from wf.task import (
    cleanup_checkpoints_task,
    combine_gene_h5ads_task,
    complete_results_task,
    gene_project_task,
    gene_spatial_task,
    gene_stats_task,
    genes_task,
    make_adata,
    make_anndata_dataset_task,
    motif_coverages_task,
    motif_peaks_task,
    motifs_task,
    registry_task,
)
from wf.utils import Run


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
            description="List of runs to be analyzed; each run must contain a "
            "run_id and fragments.tsv file; optional: condition, alternative "
            "sample name. Spaces in condition labels are normalized to '_' for "
            "downstream ArchR grouping.",
            batch_table_column=True,
            samplesheet=True,
        ),
        "genome": LatchParameter(
            display_name="genome",
            description="Reference genome for runs.",
            batch_table_column=True,
        ),
        "project_name": LatchParameter(
            display_name="project name",
            description="Name of subfolder in output directory.",
            batch_table_column=True,
            rules=[
                LatchRule(
                    regex="^[^/].*", message="project name cannot start with a '/'"
                )
            ],
        ),
        "tile_size": LatchParameter(
            display_name="tile size",
            description="The size of the tiles used for binning counts in the "
            "tile matrix.",
            batch_table_column=True,
        ),
        "n_features": LatchParameter(
            display_name="number of features",
            description="Number of features to be selected as most accessible "
            "in the tile matrix.",
            batch_table_column=True,
        ),
        "n_comps": LatchParameter(
            display_name="number of components",
            description="Number of components/dimensions to keep during "
            "dimensionality reduction with `snap.tl.spectral`.",
            batch_table_column=True,
        ),
        "resolution": LatchParameter(
            display_name="clustering resolution",
            description="Clustering resolution for Leiden algorithm; higher "
            "values result in more clusters.",
            batch_table_column=True,
        ),
        "clustering_iters": LatchParameter(
            display_name="clustering iterations",
            description="Iterations performed when selecting variable features "
            "for tile matrix.",
            batch_table_column=True,
        ),
        "leiden_iters": LatchParameter(
            display_name="leiden iterations",
            description="Number of iterations for the leiden algorithm.",
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
            description="Minimum transcription start site enrichment score "
            "required for a cell to pass filtering.",
            batch_table_column=True,
            hidden=True,
        ),
        "min_frags": LatchParameter(
            display_name="minimum fragments",
            description="Minimum number of mapped fragments required per cell "
            "to pass filtering.",
            batch_table_column=True,
            hidden=True,
        ),
        "include_y_chromosome": LatchParameter(
            display_name="include y chromosome",
            description="Include Y chromosome regions in SnapATAC2 and ArchR "
            "matrices/peak calling. Defaults to excluding Y chromosome regions.",
            batch_table_column=True,
            hidden=True,
        ),
        "disable_harmony": LatchParameter(
            display_name="disable Harmony",
            description="Skip Harmony batch correction and use the uncorrected "
            "spectral embedding. Defaults to False.",
            batch_table_column=True,
            hidden=True,
        ),
        "svg_point_size": LatchParameter(
            display_name="SVG point size",
            description="Point size for spatially variable gene and motif plots. "
            "Defaults to 12.5; passed to the spatial scatter size argument.",
            batch_table_column=True,
            hidden=True,
        ),
        "output_dir": LatchParameter(
            display_name="output directory",
            description="Folder in Latch Data to save outputs; defaults to "
            "`epi_analysis_snap`. Outputs are saved in a subfolder named "
            "with the project name.",
            batch_table_column=True,
            hidden=True,
        ),
    },
)


@workflow(metadata)
def snap_workflow(
    runs: List[Run],
    genome: Genome,
    project_name: str,
    tile_size: int = 5000,
    n_features: int = 25000,
    n_comps: int = 30,
    resolution: float = 1.0,
    clustering_iters: int = 1,
    leiden_iters: int = -1,
    min_cluster_size: int = 20,
    min_tss: float = 2.0,
    min_frags: int = 10,
    include_y_chromosome: bool = False,
    output_dir: LatchDir = LatchDir("latch:///epi_analysis_snap/"),
    svg_point_size: float = 2.0,
    disable_harmony: bool = False,
) -> LatchDir:
    """Run the complete ATX Snap spatial ATAC analysis.

    SnapATAC2 generates tile matrices, clusters, embeddings, and coverage tracks.
    ArchR provides gene accessibility, peaks, motif deviations, and differential
    statistics. Squidpy provides spatial neighborhood and autocorrelation analyses.

    ## Outputs

    Results are saved in Latch Data under `[output_dir]/[project_name]/`, defaulting to `epi_analysis_snap/[project_name]/`. Paths below are relative to that project folder. Bracketed names are placeholders for run IDs, groups, or condition labels; `[i]` identifies a condition grouping and `[page]` is a three-digit page number such as `001`.

    ### Analysis objects and browser tracks

    | Output | Contents |
    | --- | --- |
    | `anndata/combined.h5ad` | Combined SnapATAC2 tile matrix, cell metadata, embeddings, and spatial coordinates. |
    | `anndata/combined_ge.h5ad` | Combined gene accessibility AnnData object derived from ArchR gene scores. |
    | `anndata/combined_motifs.h5ad` | Combined motif AnnData object derived from ArchR chromVAR deviation scores. |
    | `anndata/[run_id]_g_converted.h5ad`, `anndata/[run_id]_m_converted.h5ad` | Per-run gene accessibility and motif AnnData objects. |
    | `anndata/combined_sm_ge.h5ad`, `anndata/combined_sm_motifs.h5ad` | Reduced objects for plotting, with analysis tables stored in `.uns`; the reduced gene object receives the final gene differential statistics. |
    | `seurat_objects/[run_id]_SeuratObj.rds`, `seurat_objects/[run_id]_SeuratObjMotif.rds` | Per-run gene and motif Seurat objects. |
    | `seurat_objects/seqlogo.rds` | Motif probability matrices for sequence logos from the genome-specific ArchR motif annotation. |
    | `[project_name]_ArchRProject/` | Saved ArchR project containing the final peak and motif analysis. |
    | `[group]_coverages/` | BigWig (`.bw`) coverage tracks grouped by cluster, sample, or condition. Sample tracks may also have copies labeled with supplied sample names. |
    | `[group]_peak_beds/` | Peak BED (`.peaks.bed`) files exported from the ArchR peak set for available groupings. |
    | `Launch_Plots/artifact.json` | Configuration for opening the results in the ATX Latch Plot Template. |

    ### Figures

    | Output | Contents |
    | --- | --- |
    | `figures/umap.png` | UMAP embeddings colored by available analysis groups. |
    | `figures/spatial_dim_[page].png` | Spatial scatter plots colored by cluster. |
    | `figures/spatial_qc_[page].png` | Spatial scatter plots of fragment counts, log10 fragment counts, and TSS enrichment. |
    | `umap.html`, `spatial_dim.html`, `spatial_qc.html` | HTML galleries for the corresponding PNG pages, saved at the project root. |
    | `figures/all_neighborhoods.pdf` | Spatial neighborhood enrichment across clusters. |
    | `figures/heatmap_genes.pdf` | ArchR cluster marker-gene heatmap (`Pval <= 0.05`, `Log2FC >= 0.1`). |
    | `figures/heatmaps_peaks-motifs.pdf` | Cluster peak and motif enrichment heatmaps. |
    | `figures/volcano_plots_[condition].pdf` | Gene differential-analysis volcano plots for condition comparisons. |
    | `figures/volcano_plots_motifs_[i]_[condition].pdf` | Motif differential-analysis volcano plots for condition comparisons. |
    | `figures/svg_spatial_genes_[page].png`, `figures/svg_spatial_motifs_[page].png` | Spatial plots of top features ranked by Moran's I, with corresponding `svg_spatial_genes.html` and `svg_spatial_motifs.html` galleries in `figures/`. |

    ArchR peak-calling PDF reports are also preserved in `figures/` before checkpoint cleanup.

    ### Tables

    All tables below are saved in `tables/`.

    | Output | Contents |
    | --- | --- |
    | `input_parameters.csv` | Recorded project, genome, tile, feature-selection, clustering, and QC parameters. |
    | `obs.csv` | Cell metadata from the combined tile object, including sample/group labels, clusters, QC metrics, and spatial information. |
    | `X_umap.csv`, `spectral.csv`, `spatial.csv` | Cell-indexed UMAP embedding, spectral embedding used for clustering, and spatial coordinates. |
    | `medians.csv` | Per-run median TSS enrichment, fragment count, and fraction of reads in peaks (FRIP), plus condition. |
    | `ranked_genes_per_[group].csv` | ArchR gene-score marker statistics, including P values, FDR, and log2 fold changes; exported with `FDR <= 1` and no log2 fold-change threshold. Groups include `cluster`, `sample`, `sample_name`, and `condition_[i]` when applicable. |
    | `genes_per_[group]_hm.csv` | Gene heatmap matrices filtered at `Pval <= 0.05` and `Log2FC >= 0.1`. |
    | `marker_peaks_per_[group].csv` | ArchR marker peaks filtered at `Pval <= 0.05` and `Log2FC >= 0.1`. Groups are `cluster`, `sample`, or `condition-[i]` (with a hyphen). |
    | `complete_peak_list_[group].csv` | Marker peaks joined to peak-set coordinates and annotations, using the same group names as the marker-peak tables. |
    | `enrichedMotifs_[group].csv` | Motif enrichment in marker peaks for `cluster`, `sample`, or `condition_[i]`. |
    | `volcanoMarkers_genes_[i]_[condition].csv` | Gene differential statistics for a condition versus the other conditions, across all cells and eligible clusters, using log2 fold change. |
    | `volcanoMarkers_motifs_[i]_[condition].csv` | Motif differential statistics for a condition versus the other conditions, across all cells and eligible clusters, using mean deviation differences. |
    | `svg_genes.csv`, `svg_motifs.csv` | Spatial autocorrelation results ranked by Moran's I. |
    | `Clusters_Peak-Call-Summary.csv` and related ArchR reports | Peak-calling summaries preserved from the motif checkpoint. |

    Sample-level gene and motif enrichment comparisons require multiple runs; sample-name gene tables also require supplied sample names. Condition comparisons require multiple conditions. Spatial autocorrelation tables and plots are produced when that analysis succeeds. Marker tables or heatmaps can be empty when no usable results are available. Motif heatmap matrices (`motif_per_*_hm.csv`) are loaded into motif AnnData `.uns` and are not exported as standalone tables.

    Intermediate files remain under `checkpoints/` after a failed execution for recovery. The workflow deletes that folder after all result-producing tasks succeed, retaining the final objects, tables, figures, and ArchR project.
    """

    anndata_dataset = make_anndata_dataset_task(
        runs=runs,
        genome=genome,
        project_name=project_name,
        min_tss=min_tss,
        min_frags=min_frags,
        include_y_chromosome=include_y_chromosome,
        tile_size=tile_size,
        output_dir=output_dir,
    )

    results, _groups = make_adata(
        runs=runs,
        anndata_dataset=anndata_dataset,
        genome=genome,
        project_name=project_name,
        resolution=resolution,
        leiden_iters=leiden_iters,
        n_comps=n_comps,
        min_cluster_size=min_cluster_size,
        min_tss=min_tss,
        min_frags=min_frags,
        include_y_chromosome=include_y_chromosome,
        tile_size=tile_size,
        n_features=n_features,
        clustering_iters=clustering_iters,
        output_dir=output_dir,
        disable_harmony=disable_harmony,
    )

    gene_project = gene_project_task(
        runs=runs,
        results_dir=results,
        project_name=project_name,
        genome=genome,
        include_y_chromosome=include_y_chromosome,
    )

    gene_results = genes_task(
        runs=runs,
        results_dir=results,
        gene_project_dir=gene_project,
        project_name=project_name,
        genome=genome,
        include_y_chromosome=include_y_chromosome,
    )

    combined_gene = combine_gene_h5ads_task(
        runs=runs,
        results_dir=results,
        gene_results_dir=gene_results,
        project_name=project_name,
    )

    results_ge = gene_spatial_task(
        runs=runs,
        results_dir=results,
        gene_results_dir=gene_results,
        gene_combined_dir=combined_gene,
        project_name=project_name,
        svg_point_size=svg_point_size,
    )

    motif_coverages = motif_coverages_task(
        gene_results_dir=gene_results,
        project_name=project_name,
    )
    motif_peaks = motif_peaks_task(
        motif_coverages_dir=motif_coverages,
        project_name=project_name,
        genome=genome,
        include_y_chromosome=include_y_chromosome,
    )
    results_motifs = motifs_task(
        runs=runs,
        results_dir=results,
        motif_peaks_dir=motif_peaks,
        project_name=project_name,
        genome=genome,
        include_y_chromosome=include_y_chromosome,
        svg_point_size=svg_point_size,
    )

    results_with_gene_stats = gene_stats_task(
        runs=runs,
        gene_results_dir=gene_results,
        gene_expression_results_dir=results_ge,
        results_root=results,
        project_name=project_name,
    )

    final_results = complete_results_task(
        base_results_dir=results,
        gene_results_dir=gene_results,
        gene_expression_results_dir=results_ge,
        gene_stats_results_dir=results_with_gene_stats,
        motif_results_dir=results_motifs,
    )

    cleaned_results = cleanup_checkpoints_task(results=final_results)
    return registry_task(runs=runs, results=cleaned_results)
