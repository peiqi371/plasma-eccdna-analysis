args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript analysis/06_go_enrichment.R GENE_RESULTS_DIR OUTPUT_DIR", call. = FALSE)
input_dir <- normalizePath(args[1L], mustWork = TRUE)
out <- args[2L]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(data.table))

input_file<-"04_original_heatmap_top30.csv"
read_genes<-function(path){
 z<-fread(path,encoding="UTF-8")
 if(!"Gene" %in% names(z)) stop("Gene column missing: ",path)
 unique(trimws(as.character(z$Gene[!is.na(z$Gene)&nzchar(z$Gene)])))
}
selected_symbols<-read_genes(file.path(input_dir,input_file))
universe_symbols<-intersect(read_genes(file.path(input_dir,"09_tested_genes_Healthy.csv")),read_genes(file.path(input_dir,"10_tested_genes_Nonmalignant.csv")))
if(!length(selected_symbols)) stop("Empty selected gene set")
if(!all(selected_symbols %in% universe_symbols)) stop("Some selected genes are not in the tested background")
# Exact SYMBOL matching only; no automatic alias replacement or one-to-many expansion.
# Ambiguous symbols are excluded from BOTH selected list and background and audited.
valid_keys<-AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db,keytype="SYMBOL")
recognized<-intersect(universe_symbols,valid_keys)
if(!length(recognized)) stop("No recognized human gene symbols")
mp<-as.data.table(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,keys=recognized,keytype="SYMBOL",columns="ENTREZID"))
mp<-unique(mp[!is.na(ENTREZID)&nzchar(ENTREZID),.(SYMBOL,ENTREZID=as.character(ENTREZID))])
mc<-mp[,.(n_entrez=uniqueN(ENTREZID)),by=SYMBOL]
audit<-merge(data.table(SYMBOL=universe_symbols),mc,by="SYMBOL",all.x=TRUE)
audit[is.na(n_entrez),n_entrez:=0L]
audit[,status:=ifelse(n_entrez==0,"Unmapped",ifelse(n_entrez>1,"Ambiguous_excluded","Mapped"))]
audit[,selected:=SYMBOL %in% selected_symbols]
fwrite(audit,file.path(out,"01_symbol_mapping_audit.csv"),bom=TRUE)
fwrite(mp,file.path(out,"02_all_symbol_entrez_mappings.csv"),bom=TRUE)
clean<-mp[SYMBOL %in% audit[status=="Mapped",SYMBOL]]
universe_ids<-unique(clean$ENTREZID)
selected_ids<-unique(clean[SYMBOL %in% selected_symbols,ENTREZID])
fwrite(clean,file.path(out,"03_used_background_mappings.csv"),bom=TRUE)
fwrite(clean[SYMBOL %in% selected_symbols],file.path(out,"04_used_candidate_mappings.csv"),bom=TRUE)
if(length(selected_ids)<2) stop("Fewer than two mapped candidates; enrichment not meaningful")
min_size <- 10L
max_size <- 500L

bp <- clusterProfiler::enrichGO(
  gene = selected_ids, universe = universe_ids, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
  keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
  pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = min_size, maxGSSize = max_size
)
all_terms <- as.data.table(as.data.frame(bp))
if (nrow(all_terms)) {
  ratio_value <- function(ratios) {
    vapply(strsplit(ratios, "/", fixed = TRUE), function(parts) {
      as.numeric(parts[1L]) / as.numeric(parts[2L])
    }, numeric(1L))
  }
  symbols_by_id <- split(clean$SYMBOL, clean$ENTREZID)
  all_terms[, gene_symbols := vapply(strsplit(geneID, "/", fixed = TRUE), function(ids) {
    paste(unique(unlist(symbols_by_id[ids])), collapse = "/")
  }, character(1L))]
  all_terms[, fold_enrichment := ratio_value(GeneRatio) / ratio_value(BgRatio)]
  setorder(all_terms, p.adjust, pvalue)
}
fwrite(all_terms, file.path(out, "go_bp_all_terms.csv"))
exploratory <- if (nrow(all_terms)) all_terms[pvalue < 0.05] else all_terms
fwrite(exploratory, file.path(out, "go_bp_exploratory_terms.csv"))
capture.output(sessionInfo(), file = file.path(out, "go_session_info.txt"))
