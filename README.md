# ArchR motif recovery

This branch repairs missing or incomplete motif data in an existing ArchRProject
that already contains peaks and a PeakMatrix. It rebuilds motif annotations,
background peaks, and MotifMatrix (both `z` deviation scores and raw `deviations`)
once, using the existing active peak set.

## Inputs

- `archr_project`: the complete project folder containing
  `Save-ArchR-Project.rds` and `ArrowFiles/`. Select the project with restored peaks.
- `genome`: the same genome used by the original analysis.
- `project_name`: the name for the repaired output.
- `output_dir`: defaults to `latch:///archr_motif_recovery/`.

Existing genes, peak coordinates, PeakMatrix, cell metadata, reduced dimensions,
embeddings and available imputation weights are retained. There is no peak calling,
clustering, Harmony, gene-score calculation, enrichment analysis or spatial analysis.
The workflow fails clearly if the supplied project lacks peaks or readable gene/peak
matrices. The earlier peak-group and Y-chromosome parameters are no longer needed:
the supplied peak set determines the regions used for motif recovery.

## Outputs

Results are written to `<output_dir>/<project_name>/`:

- `<project_name>_ArchRProject/`: saved project with rebuilt motif annotations,
  background peaks and MotifMatrix in its Arrow files.
- `motif_recovery_summary.csv`: genome and peak, motif and cell counts.

The task works on a copy and does not overwrite the input. Before upload, it reloads
the output and checks motif feature metadata, completion flags, agreement across
Arrow files and coverage of all project cells. Actual matrix data are read for one
cell per sample, and gene/peak values for those cells are compared before and after
repair. Cell metadata, embeddings and peak coordinates are also checked.

Scores are recomputed from the existing active peaks. They need not equal historical
scores computed from a different peak set (for example, cluster peaks). This workflow
does not update separately exported motif H5AD or Seurat objects.
