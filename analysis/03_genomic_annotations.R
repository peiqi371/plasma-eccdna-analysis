args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(paste("Usage: Rscript analysis/03_genomic_annotations.R",
             "ANNOTATION_BED_ROOT CHROM_SIZES REPEAT_REFERENCE ELEMENT_REFERENCE OUTPUT_DIR"),
       call. = FALSE)
}
suppressPackageStartupMessages(library(data.table))
annotation_dir <- normalizePath(args[1L], mustWork = TRUE)
chrom_sizes <- fread(args[2L], header = FALSE)
repeat_reference <- fread(args[3L], header = FALSE)
element_reference <- fread(args[4L], header = FALSE, quote = "")
out_dir <- args[5L]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

chromosomes <- paste0("chr", c(1:22, "X", "Y"))
repeat_classes <- c("LINE", "SINE", "LTR", "Simple_repeat", "Satellite", "Retroposon",
                    "RC", "scRNA", "snRNA", "srpRNA", "tRNA", "rRNA", "Low_complexity")
element_classes <- c("Gene2kU", "5'UTR", "Intron", "Exon", "3'UTR", "Gene2kD", "CpG")
repeat_background <- repeat_reference[V1 %in% chromosomes & !grepl("?", V4, fixed = TRUE),
  .(length_fraction = sum(V3 - V2) / sum(chrom_sizes[V1 %in% chromosomes, V2])),
  by = .(feature = V4)]
element_background <- element_reference[V5 %in% element_classes,
  .(length = sum(V3 - V2)), by = .(feature = V5)]
element_background[, length_fraction := length / sum(length)]

count_sample <- function(path, subgroup, kind, classes, background) {
  if (file.info(path)$size == 0) return(NULL)
  bed <- fread(path, header = FALSE, fill = TRUE, quote = "")
  if (ncol(bed) < 8L || nrow(bed) == 0L) return(NULL)
  if (kind == "Repeat") {
    bed <- bed[V1 %in% chromosomes]
    class_column <- if (ncol(bed) >= 10L) 10L else 8L
    categories <- trimws(bed[[class_column]])
    denominator <- nrow(bed)
  } else {
    categories <- bed[[ncol(bed)]]
    denominator <- sum(categories %in% classes)
  }
  if (!denominator) return(NULL)
  counts <- data.table(feature = categories[categories %in% classes])[, .(count = .N), by = feature]
  counts <- merge(counts, background[, .(feature, length_fraction)], by = "feature")
  if (!nrow(counts)) return(NULL)
  counts[, `:=`(kind = kind, subgroup = subgroup,
                sample_id = sub("_(repeat|elements)\\.bed$", "", basename(path)),
                ratio = count / denominator / length_fraction)]
  counts
}

disease <- list(Prostate = c("PCa", "Prostatic_hyperplasia", "Normal"),
                Renal = c("RCC", "Renal_cyst", "Normal"),
                Bladder = c("BCa", "BCa_Inflammation", "Normal"))
subgroups <- unique(unlist(disease, use.names = FALSE))
pooled <- list(Cancer = c("PCa", "RCC", "BCa"),
               `Non-malignant` = c("Prostatic_hyperplasia", "Renal_cyst", "BCa_Inflammation"),
               `Healthy controls` = "Normal")

read_kind <- function(kind, suffix, classes, background) {
  rbindlist(lapply(subgroups, function(subgroup) {
    files <- list.files(file.path(annotation_dir, subgroup), pattern = suffix,
                        recursive = TRUE, full.names = TRUE)
    if (!length(files)) stop("No matching annotation files for ", subgroup, " (", kind, ")")
    rbindlist(lapply(files, count_sample, subgroup = subgroup, kind = kind,
                     classes = classes, background = background))
  }))
}

sample_values <- rbindlist(list(
  read_kind("Repeat", "_repeat\\.bed$", repeat_classes, repeat_background),
  read_kind("Element", "_elements\\.bed$", element_classes, element_background)
))
fwrite(sample_values, file.path(out_dir, "normalized_mapping_ratios.csv"))

test_groups <- function(values, groups, stratum) {
  pairs <- combn(names(groups), 2L, simplify = FALSE)
  rbindlist(lapply(pairs, function(pair) {
    first <- values[subgroup %in% groups[[pair[1L]]], ratio]
    second <- values[subgroup %in% groups[[pair[2L]]], ratio]
    p <- if (length(first) && length(second)) {
      suppressWarnings(wilcox.test(log10(first), log10(second), exact = FALSE)$p.value)
    } else NA_real_
    data.table(stratum = stratum, group1 = pair[1L], group2 = pair[2L],
               n1 = length(first), n2 = length(second), P = p)
  }))
}

results <- sample_values[, {
  rbindlist(c(list(test_groups(.SD, pooled, "All")),
              lapply(names(disease), function(stratum) {
                strata <- as.list(setNames(disease[[stratum]], names(pooled)))
                test_groups(.SD, strata, stratum)
              })))
}, by = .(kind, feature)]
fwrite(results, file.path(out_dir, "pairwise_mapping_ratio_tests.csv"))
