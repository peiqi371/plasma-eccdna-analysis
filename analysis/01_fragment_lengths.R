args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript analysis/01_fragment_lengths.R BED_ROOT OUTPUT_DIR", call. = FALSE)
base_dir <- normalizePath(args[1L], mustWork = TRUE)
output_dir <- args[2L]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(data.table))

length_limits <- c(13L, 300L, 500L, 1000L, 10000L)
bins <- numeric(length(length_limits) - 1L)
min_length <- 13
max_length <- 1000

# 密度曲线的平滑带宽，单位bp
bandwidth_bp <- 5

# 置换次数：先运行999次，正式分析可改为9999
n_perm <- 9999L
set.seed(20261029)

folder_map <- c(
  RCC = "Cancer",
  PCa = "Cancer",
  BCa = "Cancer",
  Renal_cyst = "Non-malignant",
  Prostatic_hyperplasia = "Non-malignant",
  BCa_Inflammation = "Non-malignant",
  Normal = "Healthy controls"
)

groups <- c(
  "Cancer",
  "Non-malignant",
  "Healthy controls"
)

stopifnot(
  min_length > 0,
  is.finite(max_length),
  max_length > min_length,
  bandwidth_bp > 0,
  n_perm >= 1
)

# ============================================================
# 2. 查找BED文件并分组
# ============================================================

metadata <- rbindlist(
  lapply(names(folder_map), function(folder) {
    
    folder_path <- file.path(base_dir, folder)
    
    if (!dir.exists(folder_path)) {
      stop("文件夹不存在：", folder_path)
    }
    
    files <- sort(list.files(
      folder_path,
      pattern = "\\.bed$",
      full.names = TRUE,
      recursive = TRUE,
      ignore.case = TRUE
    ))
    
    if (length(files) == 0) {
      stop("文件夹中没有找到BED文件：", folder_path)
    }
    
    # Accept either sample.bed or sample_filtered.bed.
    sample_ids <- sub(
      "\\.bed$", "", basename(files),
      ignore.case = TRUE
    )
    
    sample_ids <- sub(
      "_filtered$", "", sample_ids,
      ignore.case = TRUE
    )
    
    data.table(
      sample_id = sample_ids,
      Disease = folder,
      Group = unname(folder_map[folder]),
      bed_file = files
    )
  })
)

metadata <- as.data.frame(metadata)

# 检查重复样本
duplicate_index <- duplicated(metadata$sample_id) |
  duplicated(metadata$sample_id, fromLast = TRUE)

if (any(duplicate_index)) {
  
  write.csv(
    metadata[duplicate_index, ],
    file.path(output_dir, "duplicate_samples.csv"),
    row.names = FALSE
  )
  
  stop(
    "存在重复样本编号，详见duplicate_samples.csv。",
    "请检查是否同时保存了原始和filtered版本。"
  )
}

cat("\n各疾病文件夹的样本数：\n")
print(table(metadata$Disease))

cat("\n三大组的样本数：\n")
print(table(metadata$Group))

# ============================================================
# 3. 读取BED并计算eccDNA长度
# ============================================================

length_list <- vector("list", nrow(metadata))

reading_summary <- data.frame(
  sample_id = metadata$sample_id,
  Disease = metadata$Disease,
  Group = metadata$Group,
  N_records = integer(nrow(metadata)),
  N_invalid = integer(nrow(metadata)),
  N_in_range = integer(nrow(metadata))
)

for (i in seq_len(nrow(metadata))) {
  
  path <- metadata$bed_file[i]
  
  if (file.info(path)$size == 0) {
    
    length_list[[i]] <- numeric(0)
    
  } else {
    
    bed <- tryCatch(
      fread(
        path,
        header = FALSE,
        select = 1:3,
        showProgress = FALSE
      ),
      error = function(e) {
        stop(
          "读取文件失败：", path,
          "\n", conditionMessage(e)
        )
      }
    )
    
    start <- suppressWarnings(as.numeric(bed[[2]]))
    end <- suppressWarnings(as.numeric(bed[[3]]))
    
    # 标准BED区间长度：end - start
    len <- end - start
    
    valid <- is.finite(start) &
      is.finite(end) &
      start >= 0 &
      len > 0
    
    keep <- valid &
      len >= min_length &
      len <= max_length
    
    length_list[[i]] <- len[keep]
    pooled <- len[valid & len >= 13L & len <= 10000L]
    bins <- bins + tabulate(findInterval(pooled - 1L, c(12L, length_limits[-1L])),
                            nbins = length(bins))
    
    reading_summary$N_records[i] <- length(len)
    reading_summary$N_invalid[i] <- sum(!valid)
    reading_summary$N_in_range[i] <- sum(keep)
  }
  
  if (i %% 50 == 0 || i == nrow(metadata)) {
    cat("已读取", i, "/", nrow(metadata), "个BED文件\n")
  }
}

write.csv(
  reading_summary,
  file.path(output_dir, "BED_reading_summary.csv"),
  row.names = FALSE
)

# 排除分析范围内没有有效eccDNA的样本
valid_samples <- lengths(length_list) > 0

write.csv(
  metadata[!valid_samples, ],
  file.path(output_dir, "excluded_samples.csv"),
  row.names = FALSE
)

metadata <- metadata[valid_samples, , drop = FALSE]
length_list <- length_list[valid_samples]

rownames(metadata) <- NULL
names(length_list) <- metadata$sample_id

metadata$N_eccDNA <- lengths(length_list)

sample_n <- vapply(
  groups,
  function(g) sum(metadata$Group == g),
  integer(1)
)

if (any(sample_n < 2)) {
  stop("每组至少需要两个具有有效长度数据的样本。")
}

cat("\n最终纳入分析的样本数：\n")
print(sample_n)

write.csv(
  metadata,
  file.path(output_dir, "included_samples.csv"),
  row.names = FALSE
)

# 样本层面长度统计
sample_summary <- metadata[
  , c("sample_id", "Disease", "Group", "N_eccDNA")
]

sample_summary$Median_length_bp <- vapply(
  length_list, median, numeric(1)
)

sample_summary$Mean_length_bp <- vapply(
  length_list, mean, numeric(1)
)

write.csv(
  sample_summary,
  file.path(output_dir, "sample_length_summary.csv"),
  row.names = FALSE
)

# ============================================================
# 4. 每个样本建立经验累积分布
# ============================================================

cat("\n正在构建各样本的CDF……\n")

# 在所有实际出现的长度处计算，不进行人为分箱
length_grid <- sort(unique(
  unlist(length_list, use.names = FALSE)
))

# 行为样本，列为长度
cdf_matrix <- t(vapply(
  length_list,
  function(x) ecdf(x)(length_grid),
  numeric(length(length_grid))
))

rownames(cdf_matrix) <- metadata$sample_id

# 保存组内平均样本CDF
mean_cdf <- data.frame(Length_bp = length_grid)

for (g in groups) {
  mean_cdf[[g]] <- colMeans(
    cdf_matrix[metadata$Group == g, , drop = FALSE]
  )
}

write.csv(
  mean_cdf,
  file.path(output_dir, "group_mean_sample_CDF.csv"),
  row.names = FALSE
)

# ============================================================
# 5. 样本层面的KS距离置换检验
#
# D：两组平均样本CDF的最大绝对差值
# 置换单位：样本，而不是单条eccDNA
# ============================================================

comparisons <- combn(groups, 2, simplify = FALSE)

sample_level_test <- function(group1, group2) {
  
  cat(
    "\n开始置换检验：",
    group1, "vs", group2, "\n"
  )
  
  index <- metadata$Group %in% c(group1, group2)
  
  mat <- cdf_matrix[index, , drop = FALSE]
  labels <- metadata$Group[index]
  
  n1 <- sum(labels == group1)
  n2 <- sum(labels == group2)
  n_total <- nrow(mat)
  
  observed_D <- max(abs(
    colMeans(mat[labels == group1, , drop = FALSE]) -
      colMeans(mat[labels == group2, , drop = FALSE])
  ))
  
  total_cdf <- colSums(mat)
  
  # 计算时抽取较小组，减少计算量
  small_n <- min(n1, n2)
  large_n <- max(n1, n2)
  
  permuted_D <- numeric(n_perm)
  progress_step <- max(1L, floor(n_perm / 10))
  
  for (b in seq_len(n_perm)) {
    
    selected <- sample.int(
      n_total,
      size = small_n,
      replace = FALSE
    )
    
    sum_selected <- colSums(
      mat[selected, , drop = FALSE]
    )
    
    mean_small <- sum_selected / small_n
    mean_large <- (total_cdf - sum_selected) / large_n
    
    permuted_D[b] <- max(
      abs(mean_small - mean_large)
    )
    
    if (b %% progress_step == 0 || b == n_perm) {
      cat("  已完成", b, "/", n_perm, "次置换\n")
    }
  }
  
  p_value <- (
    1 + sum(permuted_D >= observed_D - 1e-12)
  ) / (n_perm + 1)
  
  data.frame(
    Group1 = group1,
    Group2 = group2,
    N_samples_1 = n1,
    N_samples_2 = n2,
    D = observed_D,
    N_permutations = n_perm,
    P = p_value
  )
}

sample_results <- do.call(
  rbind,
  lapply(comparisons, function(z) {
    sample_level_test(z[1], z[2])
  })
)

# 三次两两比较统一BH校正
sample_results$P_BH <- p.adjust(
  sample_results$P,
  method = "BH"
)

write.csv(
  sample_results,
  file.path(output_dir, "sample_level_KS_permutation.csv"),
  row.names = FALSE
)

cat("\n样本层面KS距离置换检验结果：\n")
print(sample_results)

# ============================================================
# 7. 各样本计算密度，再计算组内平均密度
# ============================================================

cat("\n正在计算各组平均样本密度……\n")

# 所有样本使用相同网格和带宽
density_n <- max(
  2048L,
  as.integer(ceiling(max_length) + 1)
)

density_data <- do.call(
  rbind,
  lapply(groups, function(g) {
    
    indices <- which(metadata$Group == g)
    density_sum <- NULL
    density_grid <- NULL
    
    for (i in indices) {
      
      x <- length_list[[i]]
      
      if (length(x) == 1L) {
        x <- rep(x, 2L)
      }
      
      den <- density(
        x,
        bw = bandwidth_bp,
        kernel = "gaussian",
        from = 0,
        to = max_length,
        n = density_n,
        cut = 0
      )
      
      if (is.null(density_sum)) {
        density_sum <- numeric(length(den$y))
        density_grid <- den$x
      }
      
      density_sum <- density_sum + den$y
    }
    
    data.frame(
      Length_bp = density_grid,
      Density = density_sum / length(indices),
      Group = g
    )
  })
)

density_data$Group <- factor(
  density_data$Group,
  levels = groups
)

density_data$Length_kb <- density_data$Length_bp / 1000

write.csv(
  density_data,
  file.path(output_dir, "Figure1E_mean_density.csv"),
  row.names = FALSE
)

# ============================================================
if (!sum(bins)) stop("No valid fragment lengths in the pooled analysis range.")
length_bins <- data.frame(start_bp = c(13L, 301L, 501L, 1001L),
                          end_bp = length_limits[-1L],
                          count = bins, proportion = bins / sum(bins))
write.csv(length_bins, file.path(output_dir, "pooled_length_bins.csv"), row.names = FALSE)
