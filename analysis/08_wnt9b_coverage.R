args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript analysis/08_wnt9b_coverage.R BED_ROOT OUTPUT_DIR", call. = FALSE)
base_dir <- normalizePath(args[1L], mustWork = TRUE)
out_dir <- args[2L]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(data.table))

gene_chr <- "chr17"
gene_start <- 46833201L
gene_end <- 46886730L

# 图F候选区间范围，1-based坐标，用于绘图
candidate_start <- 46877244L
candidate_end <- 46877550L

group_order <- c(
  "Cancer", "Non-malignant", "Healthy controls"
)

group_map <- data.table(
  subgroup = c(
    "RCC", "PCa", "BCa",
    "Renal_cyst", "Prostatic_hyperplasia",
    "BCa_Inflammation", "Normal"
  ),
  group = c(
    rep("Cancer", 3),
    rep("Non-malignant", 3),
    "Healthy controls"
  )
)

# 1. 建立样本清单
# Accept either sample.bed or sample_filtered.bed.
registry <- rbindlist(lapply(seq_len(nrow(group_map)), function(i) {
  folder <- file.path(base_dir, group_map$subgroup[i])
  
  if (!dir.exists(folder)) {
    stop("文件夹不存在：", folder)
  }
  
  files <- list.files(
    folder,
    pattern = "\\.bed$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (!length(files)) {
    stop("没有找到 BED 文件：", folder)
  }
  
  sample_ids <- sub(
    "_filtered$", "",
    tools::file_path_sans_ext(basename(files))
  )
  
  data.table(
    sample = sample_ids,
    subgroup = group_map$subgroup[i],
    group = group_map$group[i],
    file = files
  )
}))

if (anyDuplicated(registry$sample)) {
  stop("发现重复样本 ID，请检查样本清单。")
}

totals <- registry[, .(total_samples = .N), by = group]
print(totals)

expected_totals <- c(
  Cancer = 489L,
  `Non-malignant` = 105L,
  `Healthy controls` = 307L
)

actual_totals <- totals$total_samples[
  match(names(expected_totals), totals$group)
]

if (anyNA(actual_totals) ||
    !all(actual_totals == unname(expected_totals))) {
  stop(
    "样本数与489、105、307不一致。",
    "请核查文件是否缺失或存在额外BED文件。"
  )
}

fwrite(
  registry,
  file.path(out_dir, "00_sample_registry.csv"),
  bom = TRUE
)

# 2. 读取一个样本的片段
read_intervals <- function(path) {
  size <- file.info(path)$size
  
  if (is.na(size)) stop("无法访问文件：", path)
  if (size == 0) return(NULL)
  
  # 兼容BED中可能存在的track/browser说明行
  lines <- readLines(path, warn = FALSE)
  lines <- lines[
    nzchar(trimws(lines)) &
      !grepl("^\\s*(#|track\\b|browser\\b)", lines)
  ]
  
  if (!length(lines)) return(NULL)
  
  bed <- fread(
    text = paste(lines, collapse = "\n"),
    header = FALSE,
    showProgress = FALSE
  )
  
  if (ncol(bed) < 3L) {
    stop("BED少于3列：", path)
  }
  
  coords <- data.table(
    chr = as.character(bed[[1L]]),
    start = suppressWarnings(as.numeric(bed[[2L]])),
    end = suppressWarnings(as.numeric(bed[[3L]]))
  )
  
  if (anyNA(coords)) {
    stop("发现缺失或无效坐标：", path)
  }
  
  if (any(coords$start < 0 | coords$end <= coords$start)) {
    stop("发现无效BED区间：", path)
  }
  
  dt <- coords[chr == gene_chr]
  if (!nrow(dt)) return(NULL)
  
  dt[, overlap := pmax(
    0,
    pmin(end, gene_end) - pmax(start, gene_start)
  )]
  
  # 沿用此前规则：
  # 与WNT9B重叠的长度，占eccDNA长度至少50%
  dt <- dt[overlap / (end - start) >= 0.5]
  
  if (!nrow(dt)) return(NULL)
  
  # 保留原始片段坐标，另外生成基因范围内的裁剪坐标
  dt[, `:=`(
    clipped_start = pmax(start, gene_start),
    clipped_end = pmin(end, gene_end)
  )]
  
  unique(dt[, .(
    chr, start, end, clipped_start, clipped_end
  )])
}

# 3. 逐样本累计覆盖
# 同一样本在同一位置最多贡献1次
gene_length <- gene_end - gene_start

coverage_counts <- setNames(
  lapply(group_order, function(x) integer(gene_length)),
  group_order
)

sample_QC <- copy(registry)
sample_QC[, `:=`(
  qualifying_intervals = 0L,
  covered_bases = 0L
)]

interval_list <- vector("list", nrow(registry))

for (i in seq_len(nrow(registry))) {
  if (i %% 100L == 0L) {
    message("已处理 ", i, "/", nrow(registry), " 个样本")
  }
  
  intervals <- read_intervals(registry$file[i])
  if (is.null(intervals)) next
  
  # BED区间转换为覆盖数组索引
  starts <- as.integer(
    intervals$clipped_start - gene_start + 1L
  )
  
  stops <- as.integer(
    intervals$clipped_end - gene_start + 1L
  )
  
  delta <- tabulate(starts, nbins = gene_length + 1L) -
    tabulate(stops, nbins = gene_length + 1L)
  
  sample_covered <- cumsum(delta)[seq_len(gene_length)] > 0
  
  group_name <- registry$group[i]
  
  coverage_counts[[group_name]] <-
    coverage_counts[[group_name]] +
    as.integer(sample_covered)
  
  sample_QC[i, `:=`(
    qualifying_intervals = nrow(intervals),
    covered_bases = sum(sample_covered)
  )]
  
  intervals[, `:=`(
    sample = registry$sample[i],
    group = group_name,
    subgroup = registry$subgroup[i]
  )]
  
  interval_list[[i]] <- intervals
}

# 4. 生成逐碱基覆盖比例表
coverage <- rbindlist(lapply(group_order, function(g) {
  denominator <- totals[group == g, total_samples]
  
  data.table(
    chr = gene_chr,
    position = seq.int(gene_start + 1L, gene_end),
    group = g,
    covered_samples = coverage_counts[[g]],
    total_samples = denominator,
    coverage_percent =
      100 * coverage_counts[[g]] / denominator
  )
}))

detection_check <- sample_QC[, .(
  total_samples = .N,
  detected_samples = sum(qualifying_intervals > 0)
), by = group]

message("从Filtered重新计算的WNT9B检出样本数：")
print(detection_check)

# 仅提示差异，不强行将健康组结果设为0
previous_counts <- c(
  Cancer = 124L,
  `Non-malignant` = 8L,
  `Healthy controls` = 0L
)

observed_counts <- detection_check$detected_samples[
  match(names(previous_counts), detection_check$group)
]

if (!all(observed_counts == unname(previous_counts))) {
  warning(
    "Filtered重新计算的检出数与此前coding.bed结果不同。",
    "请核对原始片段与注释文件对应关系，暂缓解释组间差异。"
  )
}

# 5. 保存数据
fwrite(
  coverage,
  file.path(out_dir, "01_per_base_sample_coverage.csv"),
  bom = TRUE
)

fwrite(
  sample_QC,
  file.path(out_dir, "02_sample_QC.csv"),
  bom = TRUE
)

fwrite(
  detection_check,
  file.path(out_dir, "03_detection_check.csv"),
  bom = TRUE
)

all_intervals <- rbindlist(interval_list)

if (nrow(all_intervals)) {
  fwrite(
    all_intervals,
    file.path(out_dir, "04_qualifying_intervals.csv"),
    bom = TRUE
  )
}


if (nrow(all_intervals)) {
  recurrent <- unique(all_intervals[end == 46877550L & start %in% 46877243:46877245])
  fwrite(recurrent, file.path(out_dir, "recurrent_candidate_intervals.csv"))
}
