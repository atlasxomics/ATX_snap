# ATX snap

**ATX snap** is a [LatchBio](https://latch.bio/) Workflow for generating Python objects and data for analysis of epigenomic [DBiT-seq](https://www.nature.com/articles/s41586-022-05094-1) experiments.  Provided a fragment file from a single-cell ATAC-seq preprocessing and alignment tool (i.e., [Chromap](https://github.com/haowenz/chromap)) and spatial information, **ATX snap** performs routine processing with [SnapATAC2](https://kzhang.org/SnapATAC2/) and [Scanpy](https://scanpy.readthedocs.io/en/stable/), returning files that can be easily input into custom scripts for more nuanced analysis without the need to perform intensive computation.

The Workflow uses SnapATAC2 for tile-matrix generation, quality control, dimensionality reduction, clustering, and coverage tracks. ArchR generates gene accessibility scores, calls peaks, and performs gene, peak, and motif analyses, including chromVAR motif deviations. Squidpy provides spatial neighborhood enrichment and spatial autocorrelation analysis. The Workflow accepts a single tissue sample or multiple tissue samples; in ATX parlance, tissue samples analyzed via DBiT-seq are termed 'Runs'. Runs are combined for analysis, with per-run and combined objects exported.


## Inputs
All input files for **ATX snap** must be on the LatchBio [file system](https://wiki.latch.bio/wiki/data/overview).  Each Run in the Workflow takes the following parameters.

* [fragments.tsv.gz file](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments): A BED-like, tab-delimited file in which each row contains an ATAC-seq fragment.

* [tissue_positions_list.csv](https://docs.atlasxomics.com/projects/AtlasXbrowser/en/latest/SpatialFolder.html): A comma-separated file in which each row contains a unique barcode, an indicator for whether the tixel is 'on-tissue' (1, 0), and a row/column index.

* [Spatial folder](https://docs.atlasxomics.com/projects/AtlasXbrowser/en/latest/SpatialFolder.html): A directory containing tissue images and experiment metadata.

* Run ID: An identifier for the Run.

* Condition (_optional_):  An experimental Condition descriptor (i.e., 'control', 'diseased').

Individual runs are batched in a combined AnnData object with the following global parameters.

* project name: A name for the output folder.

* genome: A reference genome to be used for alignment.

* clustering resolution: 'A parameter value controlling the coarseness of the clustering. Higher values lead to more clusters'; maps to the resolution parameter of [`snapatac2.tl.leiden`](https://kzhang.org/SnapATAC2/api/_autosummary/snapatac2.tl.leiden.html).
  
* number of components: Number of dimensions to keep after dimensionality reduction via SnapATAC2's spectral embedding algorithm; maps to the n_comps parameter of [`snapatac2.tl.spectral`](https://kzhang.org/SnapATAC2/api/_autosummary/snapatac2.tl.spectral.html).

* number of features: n most accessible tiles to be used for differential analysis; maps to the n_features parameter of [`https://kzhang.org/SnapATAC2/api/_autosummary/snapatac2.pp.select_features.html`](https://kzhang.org/SnapATAC2/api/_autosummary/snapatac2.pp.select_features.html).

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

6. From the Executions tab, you can view the status of the launched Workflow.  Once the Workflow has completed running, the status will change to 'Succeeded'; if the workflow has the status 'Failed', please contact an AtlasXomics Support Scientist.  You can click on the Workflow execution to view a more granular Workflow status and see output logs.

7. Workflow outputs are loaded into the LatchBio [Data module](https://wiki.latch.bio/wiki/data/overview) under `[output_dir]/[project_name]/`, defaulting to `epi_analysis_snap/[project_name]/`.


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

## Next Steps

Analysis can be performed locally or in a LatchBio [Pod](https://wiki.latch.bio/wiki/pods/overview).  Output from **ATX snap** can also analyzed in a ATX [Latch Plot Template](https://wiki.latch.bio/plots/overview). For access to ATX-specific Pods or Plots, please contact your AtlasXomics Support Scientist.
 

## Support
Questions? Comments?  Contact support@atlasxomics.com or post in AtlasXomics [Discord](https://discord.com/channels/1004748539827597413/1005222888384770108).
