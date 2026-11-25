args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript analysis/02_fragment_counts.R BED_READING_SUMMARY.csv OUTPUT_DIR", call. = FALSE)
input_file <- normalizePath(args[1L], mustWork = TRUE)
output_dir <- args[2L]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read.csv(input_file, stringsAsFactors = FALSE)
stopifnot(all(c("Disease", "Group", "N_records", "N_invalid") %in% names(dat)))
dat$Sub_group <- dat$Disease
dat$Major_group <- ifelse(dat$Disease == "Normal", "Normal", dat$Group)
dat$eccDNA_Total_Count <- dat$N_records - dat$N_invalid
stopifnot(all(dat$eccDNA_Total_Count > 0))

overall <- kruskal.test(eccDNA_Total_Count ~ Group, data = dat)
counts <- split(dat$eccDNA_Total_Count, dat$Group)
overall_summary <- do.call(rbind, lapply(names(counts), function(group) {
  x <- counts[[group]]
  data.frame(group = group, n = length(x), median = median(x),
             q1 = unname(quantile(x, 0.25)), q3 = unname(quantile(x, 0.75)))
}))
ecdf_table <- do.call(rbind, lapply(names(counts), function(group) {
  x <- counts[[group]]
  values <- sort(unique(x))
  data.frame(group = group, count = values, log10_count = log10(values),
             cumulative_proportion = ecdf(x)(values))
}))
write.csv(overall_summary, file.path(output_dir, "overall_fragment_count_summary.csv"), row.names = FALSE)
write.csv(ecdf_table, file.path(output_dir, "overall_ecdf.csv"), row.names = FALSE)
overall_pairwise <- do.call(rbind, lapply(combn(names(counts), 2, simplify = FALSE), function(pair) {
  result <- wilcox.test(counts[[pair[1L]]], counts[[pair[2L]]], exact = FALSE)
  data.frame(group1 = pair[1L], group2 = pair[2L], P = result$p.value)
}))
overall_pairwise$P_BH <- p.adjust(overall_pairwise$P, method = "BH")
write.csv(data.frame(statistic = unname(overall$statistic), P = overall$p.value),
          file.path(output_dir, "overall_kruskal_wallis.csv"), row.names = FALSE)
write.csv(overall_pairwise, file.path(output_dir, "overall_pairwise_wilcoxon.csv"), row.names = FALSE)
organ_definitions <- list(
  Prostate = list(
    cancer = "PCa",
    nonmalignant = "Prostatic_hyperplasia",
    labels = c("PCa", "Prostatic\nhyperplasia", "Healthy\ncontrols")
  ),
  Renal = list(
    cancer = "RCC",
    nonmalignant = "Renal_cyst",
    labels = c("RCC", "Renal cyst", "Healthy\ncontrols")
  ),
  Bladder = list(
    cancer = c("BCa", "BCa_NMIBC", "BCa_MIBC"),
    nonmalignant = "BCa_Inflammation",
    labels = c("BCa", "Bladder\ninflammation", "Healthy\ncontrols")
  )
)

group_levels <- c("Cancer", "Non-malignant", "Healthy controls")
organ_levels <- names(organ_definitions)

# 每个器官均使用相同的健康对照
plot_data <- do.call(
  rbind,
  lapply(organ_levels, function(organ) {
    
    definition <- organ_definitions[[organ]]
    
    selected <- dat[
      dat$Sub_group %in% c(
        definition$cancer,
        definition$nonmalignant
      ) |
        dat$Major_group == "Normal",
    ]
    
    selected$Group <- ifelse(
      selected$Major_group == "Normal",
      "Healthy controls",
      ifelse(
        selected$Sub_group %in% definition$cancer,
        "Cancer",
        "Non-malignant"
      )
    )
    
    selected$Organ <- organ
    selected
  })
)

rownames(plot_data) <- NULL

plot_data$Group <- factor(
  plot_data$Group,
  levels = group_levels
)

plot_data$Organ <- factor(
  plot_data$Organ,
  levels = organ_levels
)

plot_data$log10_count <- log10(plot_data$eccDNA_Total_Count)
plot_data$x <- as.integer(plot_data$Group)

# 4. 描述性统计
summary_table <- do.call(
  rbind,
  lapply(organ_levels, function(organ) {
    do.call(
      rbind,
      lapply(group_levels, function(group) {
        
        values <- plot_data$eccDNA_Total_Count[
          plot_data$Organ == organ &
            plot_data$Group == group
        ]
        
        stopifnot(length(values) > 0)
        
        data.frame(
          Organ = organ,
          Group = group,
          N = length(values),
          Median = median(values),
          Q1 = unname(quantile(values, 0.25)),
          Q3 = unname(quantile(values, 0.75))
        )
      })
    )
  })
)

print(summary_table)

# 5. 各器官的总体检验
kw_results <- do.call(
  rbind,
  lapply(organ_levels, function(organ) {
    
    subdat <- plot_data[plot_data$Organ == organ, ]
    
    result <- kruskal.test(
      eccDNA_Total_Count ~ Group,
      data = subdat
    )
    
    data.frame(
      Organ = organ,
      H = unname(result$statistic),
      P = result$p.value
    )
  })
)

# 三个总体检验之间进行BH校正
kw_results$P_BH <- p.adjust(kw_results$P, method = "BH")

# 6. 各器官内两两比较
comparisons <- list(
  c("Cancer", "Non-malignant"),
  c("Non-malignant", "Healthy controls"),
  c("Cancer", "Healthy controls")
)

pairwise_results <- do.call(
  rbind,
  lapply(organ_levels, function(organ) {
    
    subdat <- plot_data[plot_data$Organ == organ, ]
    
    do.call(
      rbind,
      lapply(seq_along(comparisons), function(i) {
        
        groups <- comparisons[[i]]
        
        x <- subdat$eccDNA_Total_Count[
          subdat$Group == groups[1]
        ]
        
        y <- subdat$eccDNA_Total_Count[
          subdat$Group == groups[2]
        ]
        
        result <- wilcox.test(
          x,
          y,
          alternative = "two.sided",
          exact = FALSE,
          correct = TRUE
        )
        
        data.frame(
          Organ = organ,
          Group1 = groups[1],
          Group2 = groups[2],
          Comparison = i,
          W = unname(result$statistic),
          P = result$p.value
        )
      })
    )
  })
)

# 图C共9次两两比较，统一BH校正
pairwise_results$P_BH <- p.adjust(
  pairwise_results$P,
  method = "BH"
)

print(kw_results)
print(pairwise_results)

# 7. 统计标注

write.csv(summary_table, file.path(output_dir, "organ_summary.csv"), row.names = FALSE)
write.csv(kw_results, file.path(output_dir, "organ_kruskal_wallis.csv"), row.names = FALSE)
write.csv(pairwise_results, file.path(output_dir, "organ_pairwise_wilcoxon.csv"), row.names = FALSE)
