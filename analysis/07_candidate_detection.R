args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript analysis/07_candidate_detection.R CANDIDATE_SAMPLE_DETECTION.csv OUTPUT_DIR",
       call. = FALSE)
}
suppressPackageStartupMessages(library(data.table))
out_dir <- args[2L]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
input <- fread(args[1L])
genes <- c("WNT9B", "DLL1", "IFNGR2", "CRHR1", "NKD2", "SPTBN5")
stopifnot(all(c("Gene", "sample", "subgroup", "detected") %in% names(input)),
          all(input$Gene %in% genes), all(input$detected %in% c(TRUE, FALSE)),
          !anyDuplicated(input[, .(Gene, sample)]))
groups <- list(Cancer = c("BCa", "PCa", "RCC"),
               Nonmalignant = c("BCa_Inflammation", "Prostatic_hyperplasia", "Renal_cyst"),
               Healthy = "Normal", PCa = "PCa", PH = "Prostatic_hyperplasia",
               RCC = "RCC", Cyst = "Renal_cyst", BCa = "BCa", BI = "BCa_Inflammation")
stopifnot(all(input$subgroup %in% unlist(groups)), all(genes %in% input$Gene))
comparison_labels <- c("Cancer vs Nonmalignant", "Cancer vs Healthy", "Nonmalignant vs Healthy",
                       "PCa vs PH", "RCC vs Cyst", "BCa vs BI")
summary <- input[, .(positive = sum(detected), total = .N, frequency = mean(detected)),
                 by = .(Gene, subgroup)]
fwrite(summary, file.path(out_dir, "candidate_detection_by_subgroup.csv"))

comparisons <- rbindlist(lapply(genes, function(gene) {
  x <- input[Gene == gene]
  rbindlist(lapply(comparison_labels, function(comparison) {
    pair <- strsplit(comparison, " vs ", fixed = TRUE)[[1L]]
    first <- x[subgroup %in% groups[[pair[1L]]], detected]
    second <- x[subgroup %in% groups[[pair[2L]]], detected]
    if (!length(first) || !length(second)) stop("Missing group for ", gene, " ", comparison)
    table <- matrix(c(sum(first), sum(!first), sum(second), sum(!second)), nrow = 2L)
    fit <- fisher.test(table, alternative = "two.sided")
    data.table(Gene = gene, comparison = comparison,
               detected1 = sum(first), total1 = length(first),
               detected2 = sum(second), total2 = length(second),
               difference_pp = 100 * (mean(first) - mean(second)),
               OR = unname(fit$estimate), OR_CI_low = fit$conf.int[1L],
               OR_CI_high = fit$conf.int[2L], P = fit$p.value)
  }))
}))
comparisons[, BH_P_36 := p.adjust(P, method = "BH")]
fwrite(comparisons, file.path(out_dir, "candidate_detection_tests.csv"))
