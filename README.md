# ArchR peak recovery

This branch repairs an ArchRProject produced by the workflow that omitted peak
results from its final saved project. It runs only group coverages, reproducible
peak calling, PeakMatrix construction and peak motif annotations, using the
original peak-calling settings. It does not rerun preprocessing, Harmony,
clustering, gene analysis, motif deviations or spatial analysis.

## Inputs

- `archr_project`: the complete existing project folder containing
  `Save-ArchR-Project.rds` and `ArrowFiles/`. The RDS alone is insufficient.
- `genome`: the same reference genome as the original analysis.
- `project_name`: the name for the recovered output.
- `output_dir`: defaults to `latch:///archr_peak_recovery/`.
- `group_by`: leave empty to select the final grouping from the original
  workflow. This is the last `condition_*` column in metadata order when there
  are multiple conditions; otherwise `Sample` for multiple samples, otherwise
  `Clusters`. Only this one grouping is called.
- `include_y_chromosome`: hidden; defaults to false. Match the original setting.

Automatic grouping assumes the original sample/condition calls succeeded. If a
later grouping failed in the original run, use `group_by` to select its last
successful grouping explicitly. Recovery fails on peak-calling errors rather
than silently returning a different active set.

## Outputs

Results are written to `<output_dir>/<project_name>/`:

- `<project_name>_ArchRProject/`: complete saved project with restored peak set,
  PeakMatrix, motif annotations, coverage files and peak-calling reports.
- `<group_by>_peak_beds/`: exported peak BED files.
- `peak_recovery_summary.csv`: selected grouping, genome and peak count.

The task works on a copy of the supplied project. Its metadata and embeddings
are reused, and the input project is not overwritten. The saved output is
reloaded and checked for peaks and PeakMatrix before upload. This reconstructs
results from the Arrow files; it does not guarantee byte-identical results to
an earlier run.
