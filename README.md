# Plasma eccDNA analysis

Analysis scripts for *Large-scale plasma extrachromosomal circular DNA profiling enables machine learning classification of urological disease*. The scripts generate tables and model outputs for the main and supplementary analyses.

## Scope and environment

The workflow starts from filtered eccDNA BED files, annotated BED files, and per-megabase window counts. Use the study's original reference annotations and sample split to reproduce the reported results. Read alignment, eccDNA calling, and BED annotation are upstream inputs to this release.

The analysis was developed in R 4.4.1. Install the CRAN packages required for the steps you intend to run:

```r
install.packages(c("data.table", "dplyr", "tidyr", "caret", "glmnet", "pROC",
                   "randomForestSRC", "gbm", "mboost", "plsRglm", "xgboost",
                   "e1071", "MASS", "R.utils", "digest", "kernelshap"))
```

For GO enrichment, install `clusterProfiler`, `org.Hs.eg.db` and `AnnotationDbi` through Bioconductor in a compatible R environment. The manuscript reports `clusterProfiler` 4.10.1. Package and annotation database versions can affect the terms returned. Model completion and predictions may also depend on package versions and the original per-fit timeout. For scripts run on Windows, use paths readable in the active R locale.

## Inputs

Use the following seven directory names under each BED root: `BCa`, `BCa_Inflammation`, `Normal`, `PCa`, `Prostatic_hyperplasia`, `RCC` and `Renal_cyst`. Here `Normal` denotes healthy controls. Put study inputs outside the public repository.

| Input | Format |
| --- | --- |
| Filtered calls | Per-sample BED files with chromosome, start and end in columns 1–3. Coordinates are 0-based and half-open. |
| Gene annotations | Sample files ending in `_coding.bed`. Columns 1–4 give the fragment coordinates and fragment ID. Gene start, end and symbol are in columns 6–8 for the eight-column layout or 8–10 for the ten-column layout. |
| Repeat annotations | Sample files ending in `_repeat.bed`. Repeat class is column 8 or 10 in the corresponding layout. |
| Genomic-element annotations | Sample files ending in `_elements.bed`. The final column gives the genomic element category. |
| Per-megabase counts | Sample files ending in `_perMB.bed`. Columns 1–4 are chromosome, window start, window end and count. |
| Reference annotations | hg38 chromosome sizes (chromosome and length), RepeatMasker BED (class in column 4) and genomic-element BED (category in column 5). |
| Sample split | Private CSV with `sample`, `Group` and `Dataset`. `sample` matches the gene-annotation BED basename without `.bed`. `Group` is 1 for disease and 0 for healthy controls. `Dataset` is `Train_70%` or `Test_30%`. |
| Clinical metadata | Private CSV with `sample_id`, `Major_group`, `Sub_group`, `Age` and `gender`. `Major_group` uses the study labels `Cancer`, `Inflammation` (non-malignant disease) and `Normal`. This input is used for cohort summaries. |

Input data, trained models, and sample-level outputs can contain participant information. Obtain inputs through the study's Data availability statement and keep them separate from the release. The model strategy list is included in `analysis/model_strategies.txt`.

## Run the analysis

Run from the repository root. Replace the capitalized arguments with the actual paths. The commands write CSV outputs to the named output directory.

```sh
Rscript analysis/01_fragment_lengths.R FILTERED_BED_DIR results/lengths
Rscript analysis/02_fragment_counts.R results/lengths/BED_reading_summary.csv results/counts
Rscript analysis/03_genomic_annotations.R ANNOTATION_DIR CHROM_SIZES REPEAT_REFERENCE ELEMENT_REFERENCE results/annotations
Rscript analysis/04_chromosome_counts.R WINDOW_BED_DIR CHROM_SIZES results/chromosomes
Rscript analysis/05_gene_associations.R ANNOTATION_DIR results/genes
Rscript analysis/06_go_enrichment.R results/genes results/go
Rscript analysis/07_candidate_detection.R results/genes/candidate_sample_detection.csv results/detection
Rscript analysis/08_wnt9b_coverage.R FILTERED_BED_DIR results/wnt9b
Rscript analysis/09_train_models.R ANNOTATION_DIR PRIVATE_SPLIT_CSV results/models
Rscript analysis/10_model_explanations.R results/models/model_results.rds results/models/bed_gene_ea_cache.rds results/explanations
Rscript analysis/11_held_out_metrics.R results/models/model_results.rds results/evaluation
Rscript analysis/12_cohort_summary.R PRIVATE_CLINICAL_CSV results/cohort
```

- **Results 1–2:** group and organ-level counts, empirical cumulative distributions, sample-level length densities and permutation tests, pooled length bins, mapping ratios, genomic window counts and chromosome correlations.
- **Results 3:** relative abundance and detection, Fisher tests and differential gene sets, GO biological processes, selected-gene detection and WNT9B coverage. The six-gene testing family contains 36 comparisons, comprising the three pooled comparisons and the three organ-specific cancer versus non-malignant comparisons for each gene.
- **Results 4:** the fixed cohort split, fold-specific feature selection, five-fold validation, held-out classification metrics and ROC coordinates, Kernel SHAP, and per-gene fragment counts in the training set. Figure 4 rankings are derived from the mean validation AUC in the training set.
- **Supplementary analyses:** aggregate cohort demographics, organ-stratified mapping ratios, diagnostic-subgroup detection rates, and the complete cross-validation performance table.
