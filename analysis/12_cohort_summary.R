args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript analysis/12_cohort_summary.R CLINICAL_METADATA.csv OUTPUT_DIR",
       call. = FALSE)
}
metadata <- read.csv(args[1L], fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
output_dir <- args[2L]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
required <- c("sample_id", "Major_group", "Sub_group", "Age", "gender")
stopifnot(all(required %in% names(metadata)), !anyDuplicated(metadata$sample_id))

groups <- c(Cancer = "Cancer", Inflammation = "Non-malignant", Normal = "Healthy controls")
metadata$group <- unname(groups[metadata$Major_group])
metadata$subgroup <- metadata$Sub_group
metadata$subgroup[metadata$Major_group == "Cancer" & grepl("^BCa_", metadata$subgroup)] <- "BCa"
metadata$subgroup[metadata$Major_group == "Normal"] <- "Healthy controls"
stopifnot(!anyNA(metadata$group), !anyNA(metadata$subgroup))
metadata$Age <- suppressWarnings(as.numeric(metadata$Age))
sex <- tolower(trimws(metadata$gender))
sex[!sex %in% c("male", "female")] <- "missing"
metadata$sex <- factor(sex, levels = c("male", "female", "missing"))

summarize <- function(group_field) {
  r <- split(metadata, metadata[[group_field]])
  do.call(rbind, lapply(names(r), function(group) {
    x <- r[[group]]
    counts <- table(x$sex)
    q <- quantile(x$Age, c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
    data.frame(group = group, total = nrow(x), age_available = sum(!is.na(x$Age)),
               median_age = q[2L], Q1 = q[1L], Q3 = q[3L],
               male = as.integer(counts["male"]), female = as.integer(counts["female"]),
               sex_missing = as.integer(counts["missing"]))
  }))
}

write.csv(summarize("group"), file.path(output_dir, "cohort_by_group.csv"), row.names = FALSE)
write.csv(summarize("subgroup"), file.path(output_dir, "cohort_by_subgroup.csv"), row.names = FALSE)
