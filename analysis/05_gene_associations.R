args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript analysis/05_gene_associations.R ANNOTATION_BED_ROOT OUTPUT_DIR", call. = FALSE)
base_dir <- normalizePath(args[1L], mustWork = TRUE)
output_dir <- args[2L]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

cancer_folders <- c("RCC", "PCa", "BCa")
inflam_folders <- c("Renal_cyst", "Prostatic_hyperplasia", "BCa_Inflammation")
normal_folders <- c("Normal")

get_files <- function(folders) {
  unlist(lapply(folders, function(f) {
    list.files(file.path(base_dir, f), pattern = "_coding\\.bed$", recursive = TRUE, full.names = TRUE)
  }))
}

files_cancer <- get_files(cancer_folders)
files_inflam <- get_files(inflam_folders)
files_normal <- get_files(normal_folders)

# ==============================================================================
# 3. \u6838\u5fc3\u6570\u636e\u8bfb\u53d6 (50% \u6620\u5c04\u8fc7\u6ee4 + \u8ba1\u7b97 EA \u4e0e 0/1 \u68c0\u51fa)
# ==============================================================================
read_data_from_bed <- function(bed_path) {
  sample_name <- tools::file_path_sans_ext(basename(bed_path))
  bed <- tryCatch(fread(bed_path, sep = "\t", header = FALSE, fill = TRUE), error = function(e) NULL)
  if (is.null(bed) || nrow(bed) < 1) return(NULL)
  
  # \u517c\u5bb9 8 \u5217\u548c 10 \u5217 bed \u683c\u5f0f
  if (ncol(bed) >= 10) {
    dt <- data.table(eccStart = as.numeric(bed$V2), eccEnd = as.numeric(bed$V3), 
                     geneName = as.character(bed$V10), geneStart= as.numeric(bed$V8), geneEnd = as.numeric(bed$V9))
  } else if (ncol(bed) >= 8) {
    dt <- data.table(eccStart = as.numeric(bed$V2), eccEnd = as.numeric(bed$V3), 
                     geneName = as.character(bed$V8), geneStart= as.numeric(bed$V6), geneEnd = as.numeric(bed$V7))
  } else { return(NULL) }
  
  dt <- dt[!is.na(geneName) & geneName != "" & !grepl("^NA$|^\\.$|^unassigned$", geneName, ignore.case = TRUE)]
  if (nrow(dt) == 0) return(NULL)
  
  # \u7b2c\u4e00\u7ea7\u8fc7\u6ee4\uff1a\u5361 50% \u6620\u5c04\u9608\u503c
  dt[, eccLen := abs(eccEnd - eccStart)]
  dt[, overlap := pmax(0, pmin(geneEnd, eccEnd) - pmax(geneStart, eccStart))]
  dt <- dt[eccLen > 0 & (overlap / eccLen) >= 0.5]
  if (nrow(dt) == 0) return(NULL)
  
  # \u8ba1\u7b97\u76f8\u5bf9\u4e30\u5ea6 EA
  dt_ea <- dt[, .(X = .N, L = max(abs(geneEnd - geneStart))), by = geneName]
  dt_ea <- dt_ea[L > 0]
  if (nrow(dt_ea) == 0) return(NULL)
  dt_ea[, X_div_L := X / L]
  dt_ea[, EA := (X_div_L / sum(X_div_L)) * 1e6]
  
  # \u68c0\u51fa\u6807\u8bb0\u4e3a 1
  dt_ea[, present := 1]
  dt_ea[, sample := sample_name]
  
  return(dt_ea[, .(geneName, EA, present, sample)])
}

message("\u23f3 \u6b63\u5728\u8bfb\u53d6\u5e76\u6574\u5408\u6240\u6709\u6837\u672c\u6570\u636e...")
all_files <- c(files_cancer, files_inflam, files_normal)
all_list <- lapply(all_files, read_data_from_bed)
all_list <- all_list[!sapply(all_list, is.null)]
all_data <- rbindlist(all_list, use.names = TRUE, fill = TRUE)

# ==============================================================================
# 4. \u6784\u5efa\u4e30\u5ea6\u77e9\u9635 (EA) \u4e0e \u9891\u7387\u77e9\u9635 (0/1)
# ==============================================================================
ea_mat_df <- all_data %>% tidyr::pivot_wider(id_cols = geneName, names_from = sample, values_from = EA, values_fill = 0)
occ_mat_df <- all_data %>% tidyr::pivot_wider(id_cols = geneName, names_from = sample, values_from = present, values_fill = 0)

genes <- ea_mat_df$geneName
ea_mat <- as.matrix(ea_mat_df[, -1]); rownames(ea_mat) <- genes
occ_mat <- as.matrix(occ_mat_df[, -1]); rownames(occ_mat) <- genes

samples_cancer <- intersect(tools::file_path_sans_ext(basename(files_cancer)), colnames(ea_mat))
samples_inflam <- intersect(tools::file_path_sans_ext(basename(files_inflam)), colnames(ea_mat))
samples_normal <- intersect(tools::file_path_sans_ext(basename(files_normal)), colnames(ea_mat))

# ==============================================================================
# 5. Fisher \u68c0\u9a8c\u4e0e log2FC \u8ba1\u7b97\u6838\u5fc3\u51fd\u6570
# ==============================================================================
run_diff <- function(t_samps, c_samps) {
  total_t <- length(t_samps)
  total_c <- length(c_samps)
  
  mat_occ_t <- occ_mat[, t_samps, drop=FALSE]; mat_occ_c <- occ_mat[, c_samps, drop=FALSE]
  mat_ea_t <- ea_mat[, t_samps, drop=FALSE];   mat_ea_c <- ea_mat[, c_samps, drop=FALSE]
  
  # \u7b2c\u4e8c\u7ea7\u8fc7\u6ee4\uff1a\u764c\u75c7\u7ec4\u68c0\u51fa\u7387 >= 5%
  freq_t <- rowSums(mat_occ_t) / total_t
  keep <- freq_t >= 0.05
  
  mat_occ_t <- mat_occ_t[keep, , drop=FALSE]; mat_occ_c <- mat_occ_c[keep, , drop=FALSE]
  mat_ea_t <- mat_ea_t[keep, , drop=FALSE];   mat_ea_c <- mat_ea_c[keep, , drop=FALSE]
  curr_genes <- rownames(mat_occ_t)
  
  # \u8ba1\u7b97\u500d\u6570\u5dee\u5f02 log2FC
  log2FC <- log2((rowMeans(mat_ea_t) + 1) / (rowMeans(mat_ea_c) + 1))
  
  # Fisher \u7cbe\u786e\u68c0\u9a8c\u7b97 P \u503c
  p_vals <- vapply(seq_len(nrow(mat_occ_t)), function(i) {
    t_pos <- sum(mat_occ_t[i, ]); t_neg <- total_t - t_pos
    c_pos <- sum(mat_occ_c[i, ]); c_neg <- total_c - c_pos
    tryCatch(fisher.test(matrix(c(t_pos, t_neg, c_pos, c_neg), nrow=2))$p.value, error=function(e) 1)
  }, numeric(1))
  
  df <- data.frame(Gene = curr_genes, log2FC = log2FC, PValue = p_vals, stringsAsFactors = FALSE)
  df$padj <- p.adjust(df$PValue, method = "BH")
  return(df)
}

message("\u23f3 \u6b63\u5728\u8ba1\u7b97 Cancer vs Normal \u5dee\u5f02...")
res_N <- run_diff(samples_cancer, samples_normal)
message("\u23f3 \u6b63\u5728\u8ba1\u7b97 Cancer vs Inflammation \u5dee\u5f02...")
res_I <- run_diff(samples_cancer, samples_inflam)

# ==============================================================================
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
registry<-data.table(sample=tools::file_path_sans_ext(basename(all_files)),file=all_files)
registry[,group:=ifelse(file %in% files_cancer,"Cancer",ifelse(file %in% files_inflam,"Non-malignant","Healthy controls"))]
if(anyDuplicated(registry$sample)) stop("Duplicate sample filenames")
registry[,included:=sample %in% colnames(ea_mat)]
registry[,status:=ifelse(included,"Included","Excluded: empty/unreadable or no qualifying records; inspect source")]
fwrite(registry,file.path(output_dir,"01_sample_inclusion.csv"),bom=TRUE)
groups<-list(Cancer=samples_cancer,`Non-malignant`=samples_inflam,`Healthy controls`=samples_normal)
gene_summary<-rbindlist(lapply(names(groups),function(g){
  ss<-groups[[g]]
  data.table(Gene=rownames(ea_mat),group=g,n_samples=length(ss),n_detected=rowSums(occ_mat[,ss,drop=FALSE]),detection_rate=rowMeans(occ_mat[,ss,drop=FALSE]),mean_EA=rowMeans(ea_mat[,ss,drop=FALSE]),median_EA=apply(ea_mat[,ss,drop=FALSE],1,median))
}))
fwrite(gene_summary,file.path(output_dir,"02_all_gene_group_summary.csv"),bom=TRUE)
# Match the two gene matrices explicitly rather than assuming pivot column ordering.
stopifnot(identical(rownames(ea_mat),rownames(occ_mat)),identical(colnames(ea_mat),colnames(occ_mat)))
a<-as.data.table(res_N);b<-as.data.table(res_I)
setnames(a,setdiff(names(a),"Gene"),paste0(setdiff(names(a),"Gene"),"_Cancer_vs_Healthy"))
setnames(b,setdiff(names(b),"Gene"),paste0(setdiff(names(b),"Gene"),"_Cancer_vs_Nonmalignant"))
combined<-merge(a,b,by="Gene",all=TRUE)
for(g in names(groups)){
  tag<-gsub("[ -]","_",g)
  z<-gene_summary[group==g, .(Gene,n_samples,n_detected,detection_rate,mean_EA,median_EA)]
  setnames(z,setdiff(names(z),"Gene"),paste0(tag,"_",setdiff(names(z),"Gene")))
  combined<-merge(combined,z,by="Gene",all.x=TRUE)
}
combined[,pass_Healthy:=!is.na(padj_Cancer_vs_Healthy)&padj_Cancer_vs_Healthy<.05&abs(log2FC_Cancer_vs_Healthy)>1]
combined[,pass_Nonmalignant:=!is.na(padj_Cancer_vs_Nonmalignant)&padj_Cancer_vs_Nonmalignant<.05&abs(log2FC_Cancer_vs_Nonmalignant)>1]
combined[,same_direction:=sign(log2FC_Cancer_vs_Healthy)==sign(log2FC_Cancer_vs_Nonmalignant)]
combined[,candidate_direction:=ifelse(same_direction & log2FC_Cancer_vs_Healthy>0,"Higher EA in cancer",ifelse(same_direction & log2FC_Cancer_vs_Healthy<0,"Lower EA in cancer","Opposite directions"))]
fwrite(combined,file.path(output_dir,"03_two_comparisons_full.csv"),bom=TRUE)
intersection<-combined[pass_Healthy & pass_Nonmalignant]
intersection[,min_padj:=pmin(padj_Cancer_vs_Healthy,padj_Cancer_vs_Nonmalignant)]
setorder(intersection,min_padj)
# Reproduce original Top30 BEFORE applying direction consistency.
fwrite(head(intersection,30),file.path(output_dir,"04_original_heatmap_top30.csv"),bom=TRUE)
fwrite(intersection,file.path(output_dir,"05_all_intersection_genes.csv"),bom=TRUE)
consistent<-intersection[same_direction==TRUE]
fwrite(consistent,file.path(output_dir,"06_direction_consistent_candidates.csv"),bom=TRUE)
fwrite(consistent[log2FC_Cancer_vs_Healthy>0],file.path(output_dir,"07_higher_EA_candidates.csv"),bom=TRUE)
fwrite(consistent[log2FC_Cancer_vs_Healthy<0],file.path(output_dir,"08_lower_EA_candidates.csv"),bom=TRUE)
# Background lists for enrichment: genes actually tested, not the whole genome.
fwrite(data.table(Gene=res_N$Gene),file.path(output_dir,"09_tested_genes_Healthy.csv"),bom=TRUE)
fwrite(data.table(Gene=res_I$Gene),file.path(output_dir,"10_tested_genes_Nonmalignant.csv"),bom=TRUE)
print(registry[,.(available=.N,included=sum(included),excluded=sum(!included)),by=group])
message("Strict intersection: ",nrow(intersection),"; consistent direction: ",nrow(consistent))
message("Saved to: ",output_dir)

# Run table-only analysis first in the same session; ComplexHeatmap NOT needed.

comparisons <- list(Cancer_vs_Healthy = res_N, Cancer_vs_Nonmalignant = res_I,
                    Nonmalignant_vs_Healthy = run_diff(samples_inflam, samples_normal))
selected <- lapply(comparisons, function(x) x$Gene[x$padj < 0.05 & abs(x$log2FC) > 1])
for (tag in names(comparisons)) {
  fwrite(comparisons[[tag]], file.path(output_dir, paste0(tag, ".csv")))
}
sets <- data.table(gene = sort(unique(unlist(selected))))
for (tag in names(selected)) sets[, (tag) := gene %in% selected[[tag]]]
fwrite(sets, file.path(output_dir, "differential_gene_memberships.csv"))
top <- head(intersection$Gene, 30L)
z_scores <- t(apply(ea_mat[top, , drop = FALSE], 1L, function(x) {
  scale_factor <- sd(x)
  if (scale_factor == 0) return(rep(0, length(x)))
  (x - mean(x)) / scale_factor
}))
colnames(z_scores) <- colnames(ea_mat)
write.csv(z_scores, file.path(output_dir, "selected_gene_ea_z_scores.csv"))

# Candidate detection across six clinical subgroups and healthy controls.
candidates <- c("WNT9B", "DLL1", "IFNGR2", "CRHR1", "NKD2", "SPTBN5")
bed_root <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)
relative_paths <- substring(normalizePath(registry$file, winslash = "/", mustWork = TRUE),
                            nchar(bed_root) + 2L)
registry[, subgroup := sub("/.*$", "", relative_paths)]
stopifnot(all(candidates %in% rownames(occ_mat)))
candidate_samples <- rbindlist(lapply(candidates, function(gene) {
  data.table(Gene = gene, sample = colnames(occ_mat),
             detected = as.logical(occ_mat[gene, ]))
}))
candidate_samples <- merge(candidate_samples, registry[, .(sample, subgroup)], by = "sample")
fwrite(candidate_samples, file.path(output_dir, "candidate_sample_detection.csv"))
