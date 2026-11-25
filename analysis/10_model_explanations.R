args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript analysis/10_model_explanations.R MODEL_RDS EA_CACHE_RDS OUTPUT_DIR", call. = FALSE)
model_rds <- normalizePath(args[1L], mustWork = TRUE)
cache_file <- normalizePath(args[2L], mustWork = TRUE)
out_dir <- args[3L]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(kernelshap))
suppressPackageStartupMessages(library(randomForestSRC))


script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
code_path <- dirname(normalizePath(script_file, mustWork = TRUE))
best_model <- "Lasso+RF"
saved <- readRDS(model_rds)
bundle <- saved$final_models[[best_model]]
stopifnot(inherits(bundle$fit,"rfsrc"))
pred_key<-paste(6,best_model,sep="::")
if(is.null(saved$predictions[[pred_key]])) stop("RDS中找不到Lasso+RF的第6阶段Test30预测")
pred_saved<-as.data.frame(saved$predictions[[pred_key]])
required_pred<-c("Sample","Actual","Probability","Predicted")
if(!all(required_pred%in%names(pred_saved))) stop("Test30预测表字段不完整")
clean_id<-function(x) sub("_coding$","",trimws(as.character(x)))
split<-as.data.frame(saved$split)
if(!all(c("sample","Group","Dataset")%in%names(split))) stop("样本划分表字段不完整")
split$Dataset<-trimws(as.character(split$Dataset))
train_raw<-as.character(split$sample[split$Dataset%in%c("Train_70%","Train70")])
test_raw<-as.character(split$sample[split$Dataset%in%c("Test_30%","Test30")])
if(!length(train_raw)||!length(test_raw)) stop("Train70/Test30样本为空")
if(anyDuplicated(clean_id(split$sample))) stop("样本划分表存在重复样本ID")
if(length(train_raw)!=631L||length(test_raw)!=270L) stop("样本数量异常：Train70=",length(train_raw)," Test30=",length(test_raw))
labels<-setNames(as.integer(split$Group),clean_id(split$sample))
cat("Train70：",length(train_raw)," | Test30：",length(test_raw),"\n")
# 7. Lasso+RF Test30 ROC
pred_saved$Actual<-as.integer(as.character(pred_saved$Actual))
pred_saved$Probability<-as.numeric(pred_saved$Probability)
pred_saved$Predicted<-as.integer(as.character(pred_saved$Predicted))
idx_pred<-match(clean_id(test_raw),clean_id(pred_saved$Sample))
if(anyNA(idx_pred)) stop("Lasso+RF预测表缺少Test30样本")
pred_saved<-pred_saved[idx_pred,,drop=FALSE]
if(anyNA(pred_saved$Actual)||any(!is.finite(pred_saved$Probability))) stop("Test30预测结果存在缺失或非有限值")
stopifnot(identical(pred_saved$Actual,unname(labels[clean_id(test_raw)])))
ML_FILE <- file.path(code_path, "model_methods.R")
saved$cache_file <- cache_file
if(!file.exists(ML_FILE)) stop("找不到ML.R")
source(ML_FILE,encoding="UTF-8")
suppressPackageStartupMessages(library(randomForestSRC))
if(is.null(saved$cache_file)||!file.exists(saved$cache_file)) stop("找不到RDS对应的EA缓存文件")
cache<-readRDS(saved$cache_file)
if(is.null(cache$X)) stop("EA缓存缺少X矩阵")
X<-cache$X
if(is.null(rownames(X))||is.null(colnames(X))) stop("EA矩阵缺少样本名或基因名")
x_ids<-clean_id(rownames(X))
if(anyDuplicated(x_ids)) stop("EA矩阵样本ID规范化后重复")
idx_train<-match(clean_id(train_raw),x_ids)
idx_test<-match(clean_id(test_raw),x_ids)
if(anyNA(idx_train)||anyNA(idx_test)) stop("EA矩阵缺少Train70/Test30样本")
train_ids<-rownames(X)[idx_train]
test_ids<-rownames(X)[idx_test]
cat("EA样本匹配成功：Train70=",length(train_ids)," Test30=",length(test_ids),"\n")
# 10. 取得Lasso+RF实际输入的基因和特征名称
input_genes<-as.character(bundle$input_genes)
input_features<-as.character(bundle$input_model_features)
if(!length(input_genes)||!length(input_features)) stop("Lasso+RF缺少输入基因或模型特征名称")
if(length(input_genes)!=length(input_features)) stop("基因与模型特征数量不一致")
if(anyDuplicated(input_features)) stop("模型输入特征名称重复")
if(anyNA(match(input_genes,colnames(X)))) stop("EA矩阵缺少Lasso+RF输入基因")
gene_map<-as.data.frame(saved$gene_map)
if(!all(c("Gene","Model_Feature")%in%names(gene_map))) stop("基因映射表字段不完整")
expected_features<-gene_map$Model_Feature[match(input_genes,gene_map$Gene)]
if(!identical(as.character(expected_features),input_features)) stop("Lasso+RF模型特征顺序与gene_map不一致")
cat("\nLasso+RF输入基因：\n")
print(input_genes)
make_input<-function(ids){
  xx<-as.matrix(X[ids,input_genes,drop=FALSE])
  storage.mode(xx)<-"double"
  colnames(xx)<-input_features
  if(any(!is.finite(xx))) stop("EA矩阵存在非有限值")
  as.data.frame(xx,check.names=FALSE)
}
train_x<-make_input(train_ids)
test_x<-make_input(test_ids)
# 11. 使用原ML.R的预测函数，保持和原模型预测逻辑一致
if(!exists("CalPredictScore",mode="function")) stop("ML.R中找不到CalPredictScore函数")
predict_original<-function(object,newdata){
  xx<-as.matrix(newdata[,input_features,drop=FALSE])
  storage.mode(xx)<-"double"
  prob<-as.numeric(CalPredictScore(object,xx))
  if(length(prob)!=nrow(xx)||any(!is.finite(prob))||any(prob<0|prob>1)) stop("Lasso+RF预测概率异常")
  unname(prob)
}
# 12. 检查模型重新预测与RDS原Test30概率一致
cat("\n检查Lasso+RF预测一致性……\n")
test_prob<-predict_original(bundle$fit,test_x)
max_diff<-max(abs(test_prob-pred_saved$Probability))
cat("最大概率差值：",format(max_diff,digits=10),"\n")
if(!is.finite(max_diff)||max_diff>1e-4) stop("Lasso+RF重新预测与原RDS预测不一致，停止SHAP")
cat("预测一致性检查通过\n")
# 13. Kernel SHAP：Train70背景，Test30解释
set.seed(2026)
background_n<-min(50L,nrow(train_x))
background_x<-train_x[sample(seq_len(nrow(train_x)),background_n),,drop=FALSE]
explain_x<-test_x
cat("\n开始Lasso+RF SHAP计算\n")
cat("背景样本：",nrow(background_x),"\n")
cat("解释样本：",nrow(explain_x),"\n")
cat("输入基因：",ncol(explain_x),"\n")
shap_result<-kernelshap::kernelshap(
  object=bundle$fit,
  X=explain_x,
  bg_X=background_x,
  pred_fun=function(object,X) predict_original(object,X),
  exact=FALSE,
  tol=0.01,
  max_iter=30L,
  verbose=TRUE,
  seed=2026
)
shap_mat<-as.matrix(shap_result$S)
if(!identical(dim(shap_mat),dim(as.matrix(explain_x)))) stop("SHAP矩阵维度与模型输入不一致")
if(any(!is.finite(shap_mat))) stop("SHAP矩阵存在NA或Inf")
colnames(shap_mat)<-input_features
cat("SHAP计算完成：",nrow(shap_mat),"×",ncol(shap_mat),"\n")
# 14. SHAP重要性排序
gene_labels<-input_genes
names(gene_labels)<-input_features
imp<-data.frame(
  Feature=input_features,
  Gene=unname(gene_labels[input_features]),
  MeanAbsSHAP=colMeans(abs(shap_mat))
)
imp<-imp[order(-imp$MeanAbsSHAP),,drop=FALSE]
shap_long <- do.call(rbind, lapply(seq_along(input_features), function(j) {
  feature <- input_features[j]
  data.frame(test_sample_index = seq_len(nrow(explain_x)),
             Gene = unname(gene_labels[feature]),
             SHAP = as.numeric(shap_mat[, j]), EA = as.numeric(explain_x[[feature]]))
}))

write.csv(imp, file.path(out_dir, "shap_importance.csv"), row.names = FALSE)
write.csv(shap_long, file.path(out_dir, "shap_values.csv"), row.names = FALSE)
fragment_counts <- as.data.frame(cache$ann)
fragment_counts <- unique(fragment_counts[fragment_counts$sample %in% train_ids &
  fragment_counts$geneName %in% input_genes, c("sample", "eccDNA_ID", "geneName")])
counts <- data.frame(gene = input_genes,
  fragments = as.integer(table(factor(fragment_counts$geneName, levels = input_genes))))
write.csv(counts, file.path(out_dir, "training_gene_fragment_counts.csv"), row.names = FALSE)
