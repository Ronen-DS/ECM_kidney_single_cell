# ccRCC ECM multi-cohort analysis

## Overview

Three ccRCC single-cell/nucleus RNA cohorts (A: GSE227898 snRNA, B: GSE178481
scRNA, C: GSE159115 scRNA derivative) are QC'd per-cohort, jointly integrated
(Seurat v5 RPCA), and used to test whether a 10-gene ECM panel shifts between
tumor and adjacent tissue. Cohort C has no author cell labels, so it is
labeled by resolution=1 cluster majority from the harmonized A/B population
scheme. Endothelial cells are the only population with consistent identity
and paired tumor/adjacent coverage in all three cohorts; the panel's
Tumor-vs-Adjacent effect in Endothelial cells (summed-panel, SPARC, FN1) is
tested by pseudobulk limma-voom, checked for robustness to removing
mesenchymal-marker-positive cells, and compared against cohort A's linked
bulk RNA and proteomics. `Report.Ronen.pdf` (this repository's root) is the
final report; this repository is the code and evidence behind it - see
Verification below for exactly what was independently reproduced from it.

## Input data

Raw input data is **not committed** to this repository. Set `DATA_ROOT` (the
applicant-package folder) once per machine via a local, gitignored
`.Renviron` file at the repo root:

```
CCRCC_DATA_ROOT=/path/to/applicant-package
```

Expected layout under `DATA_ROOT` (see `config.R` for the exact paths used):

```
applicant-package/
  data/rna/cohort_A_GSE227898/   snRNA, 10x raw_feature_bc_matrix archives per specimen
  data/rna/cohort_B_GSE178481/   scRNA, count.csv.gz per specimen (GSE178481_RAW.tar)
  data/rna/cohort_C_GSE159115/   scRNA, one .h5 per specimen
  data/bulk/                     cohort A's linked bulk RNA (RSEM, log2 UQ-normalized)
                                  and proteomics (log2 reference-intensity-normalized), Tumor/Normal
  metadata/                      sample_manifest.tsv, label_mapping.tsv, multimodal_crosswalk.tsv,
                                  ecm_program.tsv, bulk_samples.tsv, bulk_column_manifest.tsv,
                                  original/ (author annotations, cohort_C_samples.tsv)
```

## Setup

No `renv.lock` is used in this project; packages are installed into a
per-user library (`~/R/library-x64` on Windows, added to `.libPaths()` by
`config.R` when present - no admin rights required). See `sessionInfo.txt`
(committed) for the exact R version and package versions this was run
against: **R 4.6.1**, Seurat 5.5.1, edgeR 4.10.5, limma 3.68.5, ggplot2
4.0.3, Matrix 1.7-6, plus yaml/xml2/knitr/rmarkdown for the notebook itself.
To reproduce the environment, install the package list in `sessionInfo.txt`
with `install.packages()`/`BiocManager::install()` into any writable library
on `.libPaths()`. `rmarkdown::render()` additionally needs pandoc >= 2.8
installed separately to produce the final HTML/PDF from `analysis.Rmd`
(not required to run the pipeline itself, only to knit the notebook output).

## Execution order

`analysis.Rmd` is the single consolidated notebook, in execution order: 01
load+QC -> 02 integration -> 03 C labeling -> 04 ECM overview + population
choice -> 05 pseudobulk models + robustness -> 06 linked modalities. Each
stage is checkpointed - `run_from_scratch <- TRUE` at the top forces every
stage to recompute from source data; `FALSE` (default) reuses a stage's
checkpoint when present and only (re)computes a stage whose checkpoint is
missing. Section 02 runs the full cluster-resolution sweep (0.3, 0.5, 1.0)
together for simplicity; only resolution=1 is used downstream (Section 03)
and reported in the text.

Runtime and peak memory per stage (clean-session rerun, see Verification):
`results/verification/runtime_memory_per_stage.tsv`.

## Parameters & seeds

Every threshold and setting used (QC per cohort, integration settings, C
labeling rule, pseudobulk threshold, mesenchymal markers, random seeds) is
in `params.yaml`, read by the notebook rather than hard-coded in any chunk.

## Report-to-output map

| Report item | Output file |
|---|---|
| Table 1 (cohort contents, QC) | `results/qc/qc_combined_ABC_counts.tsv`, `results/tables/qc_cohort_{a,b,c}_per_specimen.tsv`, `results/tables/post_qc/filter_cohort_{a,b,c}_per_specimen.tsv` |
| Fig. 1 (QC outcome per cohort) | `results/qc/qc_combined_ABC.png` |
| A's 2%/5% MT sensitivity (Adjacent Lymphoid retention) | `results/qc_population_retention/retention_by_population_condition.tsv`, `results/qc_population_retention_A5pct/retention_by_population_condition.tsv` |
| Cohort C %MT / non-MT-UMI threshold derivation | `results/tables/post_qc/cohort_c_mt_binning.tsv`, `cohort_c_mt_binning_by_condition.tsv` |
| Fig. 2a (unintegrated A+B+C UMAP by cohort) | `results/figures/integration/ABC_umap_unintegrated_by_cohort.png` |
| Fig. 2b (RPCA-integrated UMAP by cohort) | `results/figures/integration/ABC_umap_rpca_by_cohort.png` |
| Fig. 2c (final labels) | `results/figures/integration/ABC_umap_final_annotation.png` |
| Cluster purity (33/37 >=80% pure; 32/37 assigned) | `results/tables/integration/ABC_cluster_majority_assignment.tsv` |
| Marker verification (11/13 markers, NDUFA4L2/ACTA2 exceptions) | `results/tables/integration/marker_verification_by_label.tsv`, `marker_verification_by_cluster.tsv` |
| Fig. 3a (ECM panel by population, per cohort) | `results/figures/ecm_panel_by_population_per_cohort.png` |
| Fig. 3b (Endothelial Adjacent->Tumor, paired donors) | `results/ecm_endothelial/plots/paired_lines_SPARC_FN1_panelsum.png` |
| Population choice (paired-donor coverage) | `results/ecm_population_choice/paired_coverage_ge10.tsv`, `paired_coverage_ge20.tsv` |
| Table 2 (Endothelial models, before/after cleaning) | `results/ecm_endothelial/clean/before_after_comparison.tsv` (also: `panel_summary.tsv`, `model_*_panel_genes.tsv`, `clean/clean_models_panel_sum_SPARC_FN1.tsv`) |
| Interaction F-test (SPARC/FN1/COL3A1) | `results/ecm_endothelial/model_e_interaction_panel_genes.tsv` |
| Donor-level directionality (B 8/8, C 3/3) | `results/ecm_endothelial/donor_level_diff.tsv` |
| Mesenchymal-marker flagged fractions | `results/ecm_endothelial/mesenchymal_flag_summary.tsv` |
| Table 3 (linked modalities, within-donor) | `results/multimodal/ecm_direction_paired_summary.tsv` (full long table: `ecm_direction_bulk_vs_protein.tsv`) |
| Matched/excluded snRNA-bulk-protein specimens | `results/multimodal/matched_subset.tsv`, `excluded_specimens.tsv`, `donor_condition_coverage.tsv` |

## Supporting outputs

Not individually called out in the report text but produced and available
in the repository:

- Per-gene donor-level plot (all 10 panel genes, paired donors): `results/ecm_endothelial/plots/donor_diff_strip_all_genes.png`
- Pooled ECM-panel dot plot (all populations, all cohorts combined): `results/figures/ecm_panel_by_population.png`
- Per-condition ECM-panel dot plot (Tumor vs Adjacent): `results/figures/ecm_panel_by_population_per_condition.png`
- Integration elbow plots: `results/figures/integration/ABC_elbow_plot.png` (unintegrated), `ABC_elbow_plot_rpca.png` (RPCA) - both flatten at 10 PCs
- Per-cohort HVG plots (`FindVariableFeatures`, top 10 labeled): `results/figures/integration/{A,B,C}/variable_features.png`
- Duplicate-symbol collapse summary (A 24, C 34 symbols summed): `results/tables/duplicate_symbol_summary.tsv`, `duplicate_symbol_summary_by_cohort.tsv`
- Cohort B recomputed-vs-supplied QC consistency check (r=1.000): `results/tables/qc_cohort_b_consistency.tsv`
- Resolution-sweep cluster UMAPs: `results/figures/integration/ABC_umap_rpca_res0.3.png`, `res0.5.png`, `res1.png`; population overlay `ABC_umap_rpca_by_population.png`
- Marker verification dot plots: `results/figures/integration/marker_dotplot_by_label.png`, `marker_dotplot_by_cluster.png`

## Inclusion/exclusion records

- **QC**: `results/tables/qc_cohort_{a,b,c}_per_specimen.tsv` (raw per-specimen recovery/median QC), `results/tables/post_qc/filter_cohort_{a,b,c}_per_specimen.tsv` (post-filter), `results/qc/qc_combined_ABC_counts.tsv` (per-cohort QC-outcome counts), `results/qc_population_retention/` (per-population x condition retention, flagged groups, paired-donor coverage before/after QC), `results/qc_population_retention_A5pct/` (A's MT cap resimulated at 5%), `results/tables/post_qc/cohort_c_mt_binning*.tsv` (cohort C's %MT/non-MT-UMI threshold derivation).
- **C labeling**: `results/tables/integration/ABC_cluster_majority_assignment.tsv` - one row per resolution=1 cluster: majority A/B population, purity %, n A/B cells, assigned/unresolved and why (purity <80% or <50 A/B cells).
- **Population choice**: `results/ecm_population_choice/paired_coverage_ge10.tsv`, `paired_coverage_ge20.tsv` - paired-donor coverage per population x cohort.
- **Models**: `results/ecm_endothelial/pseudobulk_sample_table.tsv`, `pseudobulk_clean_sample_table.tsv`, `model_*_panel_genes.tsv`, `panel_summary.tsv`, `models_genes_tested_summary.tsv`, `donor_level_diff.tsv`.
- **Robustness**: `results/ecm_endothelial/mesenchymal_flag_summary.tsv`, `results/ecm_endothelial/clean/before_after_comparison.tsv`, `clean/clean_models_panel_sum_SPARC_FN1.tsv`.
- **Multimodal**: `results/multimodal/matched_subset.tsv`, `excluded_specimens.tsv`, `donor_condition_coverage.tsv`, `ecm_direction_paired_summary.tsv`, `ecm_direction_bulk_vs_protein.tsv`.

`exploratory/` holds analyses not in the report, excluded from the main
notebook: the earlier A+B-only integration/baseline feasibility pass
(`prepare_integration.R`, `baseline_scale_pca/clusters/umap.R`,
`prepare_rpca_pca.R`, `find_rpca_anchors.R`, `integrate_rpca*.R`,
`rpca_umap_screen.R`, `rpca_clusters_umap.R`, `replot_elbow_ab.R`), superseded
once A+B+C joint integration proved feasible, plus an early single-population
coverage check superseded by `ecm_population_choice_coverage.R`
(`stromal_donor_coverage.R`, `plot_by_population.R`).

## Verification

`results/report_check.tsv` cross-checks every number in Tables 1-3, the
Figure 1 captions, and the body text against the actual output files,
produced by independently rerunning the pipeline (not read back from values
already quoted in the report). **105/105 items matched exactly; 0 mismatches,
0 results changed to force agreement.**

**5 of 6 stages were reproduced from scratch; stage 02 was attempted but not
completed.** Stage 01 was rerun from the raw input files (not from any
checkpoint) and reproduced exactly (Table 1 retention: A 87.4%, B 76.5%, C
75.2%). Stage 02 (integration) was also attempted from raw inputs, but
`IntegrateData` crashed with an out-of-memory error (`std::bad_alloc`) after
~5 hours in one long R session, on a standard laptop configuration; this was
not retried, and stage 02's existing checkpoint (from an earlier successful
run, independently confirmed to reproduce the same numbers used throughout
this repository) was used unchanged. Stages 03-06 were rerun from that
checkpoint, each in a clean R session. See
`results/verification/runtime_memory_per_stage.tsv` for timing/memory per
stage, including the stage-02 crash.

## AI assistance

Claude Code wrote and ran most of the R code, from my specifications, and
assembled this repository. All analytical decisions were mine, and AI
suggestions were accepted only after I checked them against the data.
Perplexity was used to check for support in the literature for biological
aspects (for example endothelial involvement in kidney cancer) and for
citation verification. Claude (chat) was used to assist in writing the
report and Gemini helped find small textual errors. Consequential outputs
were checked by comparing recomputed QC metrics with cohort B's
author-supplied fields, cross-checking QC counts against Table 1, and
verifying reported numbers against the outputs (`results/report_check.tsv`).

## Time spent

Approximately 6 hours of hands-on work, excluding prolonged computational
runtime on a standard laptop configuration.
