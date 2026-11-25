args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript analysis/04_chromosome_counts.R WINDOW_BED_ROOT CHROM_SIZES OUTPUT_DIR", call. = FALSE)
base_dir <- normalizePath(args[1L], mustWork = TRUE)
chrom_size_file <- normalizePath(args[2L], mustWork = TRUE)
output_dir <- args[3L]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))

chrom_sizes <- read.table(chrom_size_file, header = FALSE, col.names = c("chr", "length")) %>%
  filter(chr %in% paste0("chr", c(1:22, "X", "Y"))) %>%
  mutate(length_10e8 = length / 1e8,
         chr_name = gsub("chr", "", chr))

# 2. \u8bfb\u53d6 7 \u4e2a\u6587\u4ef6\u5939\u4e2d\u7684 _perMB.bed \u6587\u4ef6 (Y\u8f74\u6570\u636e)
target_folders <- c("Normal", "RCC", "BCa", "PCa", "Prostatic_hyperplasia", "Renal_cyst", "BCa_Inflammation")

all_sample_data <- list()
window_summaries <- list()

for (folder in target_folders) {
  folder_path <- file.path(base_dir, folder)
  files <- list.files(folder_path, pattern = "_perMB\\.bed$", full.names = TRUE, recursive = TRUE)
  if (!length(files)) stop("No per-Mb window files found for ", folder)
  windows <- vector("list", length(files))
  
  for (i in seq_along(files)) {
    temp_df <- fread(files[i], header = FALSE, select = 1:4)
    windows[[i]] <- temp_df
    chr_counts <- temp_df %>%
      group_by(V1) %>%
      summarise(eccDNA_count = sum(V4), .groups = "drop") %>%
      rename(chr = V1) %>%
      mutate(Folder = folder)
    
    all_sample_data[[length(all_sample_data) + 1]] <- chr_counts
  }
  window_summaries[[folder]] <- rbindlist(windows)[,
    .(mean_count_per_sample = sum(V4) / length(files)),
    by = .(chromosome = V1, start = V2, end = V3)]
}

windows <- rbindlist(window_summaries, idcol = "subgroup")
fwrite(windows, file.path(output_dir, "window_counts_by_subgroup.csv"))
# 3. \u6570\u636e\u6574\u5408\u4e0e\u5206\u7ec4
all_data_df <- bind_rows(all_sample_data)

plot_data <- all_data_df %>%
  mutate(
    Status = case_when(
      Folder %in% c("RCC", "PCa", "BCa") ~ "Cancer",
      Folder %in% c("Renal_cyst", "Prostatic_hyperplasia", "BCa_Inflammation") ~ "Non-malignant",
      Folder == "Normal" ~ "Healthy controls"
    )
  ) %>%
  filter(chr %in% paste0("chr", c(1:22, "X", "Y"))) %>%
  group_by(Status, chr) %>%
  summarise(mean_count = mean(eccDNA_count), .groups = "drop") %>%
  mutate(mean_count_k = mean_count / 1000) %>%
  left_join(chrom_sizes, by = "chr")

# \u56fa\u5b9a\u9762\u677f\u7684\u663e\u793a\u987a\u5e8f
plot_data$Status <- factor(plot_data$Status, levels = c("Cancer", "Non-malignant", "Healthy controls"))

write.csv(plot_data, file.path(output_dir, "chromosome_counts.csv"), row.names = FALSE)
cor_table <- do.call(rbind,lapply(split(plot_data,plot_data$Status),function(z){
  ct<-cor.test(z$length_10e8,z$mean_count_k,method="pearson")
  data.frame(group=as.character(z$Status[1]),n_chromosomes=nrow(z),R=unname(ct$estimate),p=ct$p.value)
}))
write.csv(cor_table,file.path(output_dir,"FigureE_correlations.csv"),row.names=FALSE)
message("Saved to: ",output_dir)
