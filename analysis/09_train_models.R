args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript analysis/09_train_models.R ANNOTATION_BED_ROOT PRIVATE_SPLIT_CSV OUTPUT_DIR", call. = FALSE)
base_dir <- normalizePath(args[1L], mustWork = TRUE)
split_file <- normalizePath(args[2L], mustWork = TRUE)
output_root <- args[3L]
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(randomForestSRC))
suppressPackageStartupMessages(library(gbm))
suppressPackageStartupMessages(library(mboost))
suppressPackageStartupMessages(library(plsRglm))
suppressPackageStartupMessages(library(xgboost))
suppressPackageStartupMessages(library(e1071))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(digest))


script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
code_path <- dirname(normalizePath(script_file, mustWork = TRUE))
source(file.path(code_path, "model_methods.R"), encoding = "UTF-8")
seed <- 2026L
alpha <- 0.1
max_features <- 10L
time_limit <- 60
out_dir <- output_root
cache_file <- file.path(out_dir, "bed_gene_ea_cache.rds")
force_rebuild_cache <- FALSE

write_csv <- function(x, name) {
  data.table::fwrite(x, file.path(out_dir, name), bom = TRUE, na = "NA")
}
split <- read.csv(split_file, fileEncoding = "UTF-8-BOM", check.names = FALSE,
                  stringsAsFactors = FALSE)
stopifnot(all(c("sample", "Group", "Dataset") %in% names(split)))
if (anyDuplicated(split$sample)) stop("原划分存在重复样本名。")
if (!all(split$Group %in% 0:1)) stop("Group必须为0/1。")
if (!all(split$Dataset %in% c("Train_70%", "Test_30%"))) stop("Dataset标签异常。")
train_ids <- split$sample[split$Dataset == "Train_70%"]
test_ids <- split$sample[split$Dataset == "Test_30%"]
if (!length(train_ids) || !length(test_ids)) stop("划分为空。")
labels <- setNames(as.integer(split$Group), split$sample)
write_csv(split, "Task1_Sample_Split_Record.csv")

# 保留BED第4列eccDNA ID及已有基因注释，不重新进行基因组注释。
read_bed <- function(path) {
  sample_id <- tools::file_path_sans_ext(basename(path))
  bed <- fread(path, sep = "\t", header = FALSE, fill = TRUE)
  if (!nrow(bed) || ncol(bed) < 8L) stop("空BED或列数不足：", path)
  if (ncol(bed) >= 10L) {
    gene_col <- 10L; gs_col <- 8L; ge_col <- 9L
  } else {
    gene_col <- 8L; gs_col <- 6L; ge_col <- 7L
  }
  dt <- data.table(sample = sample_id, eccDNA_ID = as.character(bed[[4L]]),
                   eccStart = as.numeric(bed[[2L]]), eccEnd = as.numeric(bed[[3L]]),
                   geneName = as.character(bed[[gene_col]]),
                   geneStart = as.numeric(bed[[gs_col]]), geneEnd = as.numeric(bed[[ge_col]]))
  dt <- dt[!is.na(geneName) & geneName != "" &
             !grepl("^NA$|^\\.$|^unassigned$", geneName, ignore.case = TRUE)]
  dt[, eccLen := abs(eccEnd - eccStart)]
  dt[, overlap := pmax(0, pmin(geneEnd, eccEnd) - pmax(geneStart, eccStart))]
  dt <- dt[eccLen > 0 & overlap / eccLen >= 0.5]
  if (!nrow(dt)) stop("样本无通过筛选的记录：", sample_id)
  if (anyNA(dt$eccDNA_ID) || any(dt$eccDNA_ID %in% c("", "."))) {
    stop("BED第4列存在空eccDNA ID：", sample_id)
  }
  ea <- dt[, .(X = .N, L = max(abs(geneEnd - geneStart))), by = geneName]
  ea <- ea[L > 0]
  if (!nrow(ea)) stop("样本无有效基因长度：", sample_id)
  ea[, EA := (X / L) / sum(X / L) * 1e6]
  ea[, sample := sample_id]
  # 与EA中实际使用的有效基因保持一致。
  dt <- dt[geneName %in% ea$geneName]
  list(ea = ea[, .(sample, geneName, EA, X)],
       annotations = unique(dt[, .(sample, eccDNA_ID, geneName)]))
}
dirs <- c("BCa", "PCa", "RCC", "BCa_Inflammation", "Prostatic_hyperplasia", "Renal_cyst", "Normal")
files <- unlist(lapply(dirs, function(d) list.files(file.path(base_dir, d),
                                                    pattern = "_coding\\.bed$", recursive = TRUE, full.names = TRUE)), use.names = FALSE)
ids <- tools::file_path_sans_ext(basename(files))
files <- files[ids %in% split$sample]
ids <- tools::file_path_sans_ext(basename(files))
if (anyDuplicated(ids)) stop("BED文件存在重复样本名，需先核对。")
if (length(setdiff(split$sample, ids))) stop("部分划分样本找不到BED文件。")
files <- sort(files)
info <- file.info(files)
manifest <- data.frame(Path = normalizePath(files, winslash = "/", mustWork = TRUE),
                       Size = as.numeric(info$size), MTime = as.numeric(info$mtime))
cache_signature <- list(version = "BED_ID_EA_v1", manifest = manifest,
                        samples = split$sample, overlap_fraction = 0.5)
use_cache <- FALSE
if (file.exists(cache_file) && !force_rebuild_cache) {
  message("检查已保存的数据缓存……")
  cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
  if (!is.null(cached) && identical(cached$signature, cache_signature)) {
    ea <- cached$ea; ann <- cached$ann; X <- cached$X; gene_map <- cached$gene_map
    use_cache <- TRUE
    message("缓存匹配：直接读取EA矩阵和样本—eccDNA ID—基因表，不再解析BED。")
  }
  rm(cached)
}
if (!use_cache) {
  message("首次读取或数据已变化：读取BED、计算EA并保存共享缓存……")
  raw <- lapply(seq_along(files), function(i) {
    if (i %% 100L == 0L) message(i, "/", length(files))
    read_bed(files[i])
  })
  ea <- rbindlist(lapply(raw, `[[`, "ea"))
  ann <- rbindlist(lapply(raw, `[[`, "annotations"))
  rm(raw)
  genes <- sort(unique(ea$geneName))
  gene_map <- data.frame(Gene = genes, Model_Feature = make.names(genes, unique = TRUE))
  X <- matrix(0, nrow = nrow(split), ncol = length(genes),
              dimnames = list(split$sample, genes))
  X[cbind(match(ea$sample, split$sample), match(ea$geneName, genes))] <- ea$EA
  if (any(!is.finite(X))) stop("EA矩阵出现非有限值。")
  saveRDS(list(signature = cache_signature, ea = ea, ann = ann, X = X,
               gene_map = gene_map), cache_file, compress = FALSE)
  message("共享缓存已保存：", cache_file)
}

# 仅接收当前训练部分的样本ID，禁止验证标签参与筛选。
select_genes <- function(training_ids) {
  a <- X[training_ids[labels[training_ids] == 1L], , drop = FALSE]
  b <- X[training_ids[labels[training_ids] == 0L], , drop = FALSE]
  if (!nrow(a) || !nrow(b)) stop("当前训练部分缺少某一类别。")
  keep <- colMeans(a > 0) >= 0.05
  a <- a[, keep, drop = FALSE]; b <- b[, keep, drop = FALSE]
  if (!ncol(a)) stop("出现频率过滤后无基因。")
  pv <- vapply(seq_len(ncol(a)), function(j) {
    if (var(c(a[, j], b[, j])) == 0) return(1)
    wilcox.test(a[, j], b[, j], exact = FALSE)$p.value
  }, numeric(1))
  tab <- data.frame(Gene = colnames(a),
                    log2FC = log2((colMeans(a) + 1) / (colMeans(b) + 1)),
                    PValue = pv, Disease_Frequency = colMeans(a > 0))
  tab <- tab[is.finite(tab$PValue) & tab$PValue < 0.05 & abs(tab$log2FC) > 1, ]
  # P值并列时按基因名固定排序，确保结果可重复。
  tab <- tab[order(tab$PValue, tab$Gene), ]
  head(tab, max_features)
}
matrix_for <- function(ids, candidate_genes) {
  z <- X[ids, candidate_genes, drop = FALSE]
  colnames(z) <- gene_map$Model_Feature[match(candidate_genes, gene_map$Gene)]
  z
}

# ==================== 模型接口 ====================
# 下述兼容修正在本脚本内生效，不改动磁盘上的ML.R。
RunXGBoost <- function(Train_set, Train_label, mode, classVar) {
  y <- as.numeric(Train_label[[classVar]])
  folds <- caret::createFolds(factor(y), k = 5, list = TRUE)
  params <- list(objective = "binary:logistic", eval_metric = "logloss",
                 max_depth = 2L, eta = 1, nthread = 2L)
  best <- vapply(folds, function(v) {
    tr <- xgboost::xgb.DMatrix(Train_set[-v, , drop = FALSE], label = y[-v])
    va <- xgboost::xgb.DMatrix(Train_set[v, , drop = FALSE], label = y[v])
    b <- xgboost::xgb.train(params = params, data = tr, nrounds = 10L,
                            watchlist = list(test = va), verbose = 0)
    which.min(b$evaluation_log$test_logloss)
  }, integer(1))
  rounds <- as.integer(names(which.max(table(best))))
  fit <- xgboost::xgb.train(params = params,
                            data = xgboost::xgb.DMatrix(Train_set, label = y), nrounds = rounds, verbose = 0)
  # 特征名放在模型外层bundle，避免跨版本修改Booster对象。
  if (mode == "Model") return(fit)
  xgboost::xgb.importance(feature_names = colnames(Train_set), model = fit)$Feature
}
# 原plsRglm先运行的CV结果并未用于拟合，保留最终拟合参数，省去无效重复计算。
RunplsRglm <- function(Train_set, Train_label, mode, classVar) {
  fit <- plsRglm::plsRglm(Train_label[[classVar]], as.data.frame(Train_set),
                          modele = "pls-glm-logistic", verbose = FALSE, sparse = TRUE)
  fit$subFeature <- colnames(Train_set)
  if (mode == "Model") return(fit)
  extract_features(fit, colnames(Train_set))
}
# RF不生成预测所不需要的样本间proximity大矩阵；树数和nodesize不变。
RunRF <- function(Train_set, Train_label, mode, classVar) {
  d <- as.data.frame(Train_set); d[[classVar]] <- factor(Train_label[[classVar]], levels = 0:1)
  fit <- randomForestSRC::rfsrc(stats::reformulate(colnames(Train_set), response = classVar),
                                data = d, ntree = 1000, nodesize = 5, importance = TRUE, proximity = FALSE, forest = TRUE)
  fit$subFeature <- colnames(Train_set)
  if (mode == "Model") return(fit)
  # 保持原var.select方法；不擅自换成VIMP>0筛选。
  fun <- get0("var.select", envir = asNamespace("randomForestSRC"), inherits = FALSE)
  if (is.null(fun)) stop("当前randomForestSRC没有var.select；RF前置筛选不可用，未替换筛选规则。")
  fun(fit, verbose = FALSE)$topvars
}
extract_features <- function(fit, input) {
  if (inherits(fit, "glmnet")) {
    b <- as.matrix(coef(fit)); z <- rownames(b)[is.finite(b[, 1]) & b[, 1] != 0]
  } else if (inherits(fit, "glm")) {
    b <- coef(fit); z <- names(b)[is.finite(b) & b != 0]
  } else if (inherits(fit, "glmboost")) {
    b <- coef(fit); z <- names(b)[is.finite(b) & b != 0]
  } else if (inherits(fit, "plsRglmmodel")) {
    b <- fit$Coeffs; z <- rownames(b)[is.finite(b[, 1]) & b[, 1] != 0]
  } else if (inherits(fit, "gbm")) {
    imp <- summary(fit, plotit = FALSE); z <- as.character(imp$var[imp$rel.inf > 0])
  } else if (inherits(fit, "xgb.Booster")) {
    z <- xgboost::xgb.importance(feature_names = input, model = fit)$Feature
  } else {
    # RF/SVM/LDA/NB无非零系数筛选：统计真正提供给最终模型的全部特征。
    z <- input
  }
  intersect(input, setdiff(z, c("(Intercept)", "Intercept")))
}
predict_bundle <- function(bundle, ids) {
  z <- matrix_for(ids, bundle$input_genes); fit <- bundle$fit
  if (inherits(fit, "xgb.Booster")) {
    prob <- as.numeric(predict(fit, xgboost::xgb.DMatrix(z)))
    pred <- as.integer(prob > 0.5)
  } else if (inherits(fit, "svm")) {
    ans <- predict(fit, as.data.frame(z), probability = TRUE)
    prob <- as.numeric(attr(ans, "probabilities")[, "1"])
    pred <- as.integer(as.character(ans))
  } else if (inherits(fit, "train")) {
    prob <- as.numeric(predict(fit, as.data.frame(z), type = "prob")[["1"]])
    pred <- as.integer(as.character(predict(fit, as.data.frame(z), type = "raw")))
  } else {
    prob <- as.numeric(CalPredictScore(fit, z))
    pred <- as.integer(PredictClass(fit, z))
  }
  if (length(prob) != length(ids) || length(pred) != length(ids) ||
      any(!is.finite(prob)) || any(!pred %in% 0:1) || any(prob < 0 | prob > 1)) {
    stop("预测概率/类别异常，未据此计算指标。")
  }
  y <- unname(labels[ids])
  TP <- sum(pred == 1 & y == 1); TN <- sum(pred == 0 & y == 0)
  FP <- sum(pred == 1 & y == 0); FN <- sum(pred == 0 & y == 1)
  div <- function(a,b) if (b == 0) NA_real_ else a/b
  sens <- div(TP,TP+FN); spec <- div(TN,TN+FP)
  den <- sqrt(as.double(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  metrics <- data.frame(Accuracy = (TP+TN)/length(y), F1_Score = div(2*TP,2*TP+FP+FN),
                        Sens = sens, Spec = spec, MCC = div(TP*TN-FP*FN,den),
                        H_Score = if (sens+spec == 0) 0 else 2*sens*spec/(sens+spec),
                        AUC = as.numeric(pROC::auc(pROC::roc(y,prob,levels=c(0,1),direction="<",quiet=TRUE))),
                        TP=TP,TN=TN,FP=FP,FN=FN)
  list(metrics=metrics, predictions=data.frame(Sample=ids,Actual=y,Probability=prob,Predicted=pred))
}
count_ecc <- function(gs, ids) {
  dt <- ann[sample %in% ids & geneName %in% gs]
  list(Gene_Count=length(gs), eccDNA_Count=uniqueN(dt[, .(sample,eccDNA_ID)]),
       by_gene=dt[, .(Unique_eccDNA_Count=.N),by=geneName])
}
error_log <- data.frame(Model=character(),Fold=character(),Phase=character(),Reason=character())
record_error <- function(model,stage,phase,msg) {
  error_log[nrow(error_log)+1L,] <<- list(model,as.character(stage),phase,msg)
}
limited <- function(expr) R.utils::withTimeout(expr, timeout=time_limit, onTimeout="error")

methods <- unique(readLines(file.path(code_path, "model_strategies.txt"), warn = FALSE))
methods <- methods[!is.na(methods) & nzchar(methods)]
if (!length(methods)) stop("模型清单为空。")
parts <- strsplit(methods,"\\+")
if (any(lengths(parts)>2L)) stop("清单包含不支持的多级组合。")
pre_methods <- unique(vapply(parts,function(x) if(length(x)==2L)x[1] else "simple",character(1)))
set.seed(seed)
folds <- caret::createFolds(factor(labels[train_ids]),k=5,list=TRUE)
outer <- data.frame(Sample=train_ids,Fold=NA_integer_)
for(i in seq_along(folds)) outer$Fold[folds[[i]]] <- i
state <- list(models=methods,fold_models=list(),final_models=list(),fold_metrics=list(),
              test_metrics=list(),predictions=list(),selected=list(),candidate_tables=list(),
              outer_folds=outer,split=split,gene_map=gene_map,cache_file=cache_file,
              config=list(seed=seed,time_limit=time_limit,alpha=alpha,max_features=max_features,
                           code_path=code_path,base_dir=base_dir,split_file=split_file,output_root=output_root,
                           cache_file=cache_file,ranking="Mean_AUC",
                           feature_selection="Disease frequency >= 0.05; Wilcoxon raw P < 0.05; abs(log2FC) > 1; order by PValue then Gene; max 10",
                           source_ML=readLines(file.path(code_path,"model_methods.R"),warn=FALSE)),
              session_info=sessionInfo())
model_rds <- file.path(out_dir,"model_results.rds")
checkpoint <- function() {
  state$errors <<- error_log
  saveRDS(state,model_rds,compress=FALSE)
}
message("模型清单共",length(methods),"种组合。每阶段共享候选基因及前置筛选。")
successful <- methods
for(stage in 1:6) {
  if(stage<=5L) {
    fit_ids <- train_ids[-folds[[stage]]]; valid_ids <- train_ids[folds[[stage]]]
    active <- methods
  } else {
    fit_ids <- train_ids; valid_ids <- test_ids
    fm <- rbindlist(state$fold_metrics,fill=TRUE)
    successful <- if(nrow(fm)) fm[,.(N=uniqueN(Fold)),by=Model][N==5L,Model] else character()
    active <- successful
  }
  message("阶段",stage,"/6：",if(stage<=5L)"外层CV" else "最终Train70拟合及Test30评估")
  candidate <- select_genes(fit_ids)
  if(nrow(candidate)<2L) stop("当前阶段差异筛选后不足2个基因，停止；此前结果保存在RDS。")
  state$candidate_tables[[stage]] <- candidate
  z <- matrix_for(fit_ids,candidate$Gene)
  yl <- data.frame(Group=unname(labels[fit_ids]),row.names=fit_ids)
  pre <- list(simple=colnames(z))
  for(j in seq_along(pre_methods)) {
    pm <- pre_methods[j]; if(pm=="simple")next
    # 每阶段每种前置筛选只做一次；失败仅影响依赖它的组合。
    set.seed(seed+stage*10000L+j)
    pre[[pm]] <- tryCatch({
      v <- limited(RunML(pm,z,yl,mode="Variable",classVar="Group"))
      intersect(as.character(v),colnames(z))
    },error=function(e) {
      record_error(pm,stage,"前置筛选",conditionMessage(e)); character()
    })
  }
  stage_models <- list()
  # 每个完整模型流程独立拟合；缓存键包含完整模型名称，避免不同前置筛选流程误复用同一底层模型。
  fit_cache <- new.env(parent=emptyenv())
  for(m in active) {
    ix <- match(m,methods); pp <- parts[[ix]]
    pm <- if(length(pp)==2L)pp[1] else "simple"; md <- tail(pp,1)
    v <- pre[[pm]]
    message("[",stage,"/6] ",ix,"/",length(methods)," ",m)
    if(length(v)<2L) {record_error(m,stage,"前置筛选阻断","有效特征少于2个");next}
    phase <- "模型训练"
    ans <- tryCatch({
      key <- digest::digest(list(model=m,method=md,features=sort(v)),algo="xxhash64")
      if(exists(key,fit_cache,inherits=FALSE)) {
        obj <- get(key,fit_cache,inherits=FALSE)
        if(inherits(obj,"condition"))stop(obj)
      } else {
        set.seed(seed+stage*10000L+100L+ix)
        obj <- tryCatch(limited(RunML(md,z[,v,drop=FALSE],yl,mode="Model",classVar="Group")),error=identity)
        assign(key,obj,fit_cache)
        if(inherits(obj,"condition"))stop(obj)
      }
      used <- limited(extract_features(obj,v))
      if(!length(used))stop("最终模型无有效基因，未列入成功名单。")
      if(length(used)>10L || length(v)>10L)stop("基因数超过10，停止记录该模型。")
      bundle <- list(fit=obj,input_genes=gene_map$Gene[match(v,gene_map$Model_Feature)],
                     used_genes=gene_map$Gene[match(used,gene_map$Model_Feature)],method=m,
                     input_model_features=v,training_ids=fit_ids)
      phase <- "预测"
      ev <- limited(predict_bundle(bundle,valid_ids))
      list(bundle=bundle,ev=ev)
    },error=function(e){record_error(m,stage,phase,conditionMessage(e));NULL})
    if(is.null(ans))next
    stage_models[[m]] <- ans$bundle
    ev <- ans$ev
    if(stage<=5L) {
      state$fold_metrics[[length(state$fold_metrics)+1L]] <-
        cbind(Model=m,Fold=stage,ev$metrics,Final_Gene_Count=length(ans$bundle$used_genes))
    } else {
      co <- count_ecc(ans$bundle$used_genes,fit_ids)
      state$test_metrics[[m]] <- cbind(Model=m,ev$metrics,
                                       Final_Gene_Count=co$Gene_Count,Final_eccDNA_Count=co$eccDNA_Count)
      st <- data.frame(Model=m,Gene=ans$bundle$used_genes)
      st <- merge(st,as.data.frame(co$by_gene),by.x="Gene",by.y="geneName",all.x=TRUE,sort=FALSE)
      if(inherits(ans$bundle$fit,"glmnet") || inherits(ans$bundle$fit,"glm")) {
        cc <- coef(ans$bundle$fit)
        if(is.matrix(cc) || inherits(cc,"Matrix"))cc <- as.matrix(cc)[,1]
        st$Coefficient <- unname(cc[gene_map$Model_Feature[match(st$Gene,gene_map$Gene)]])
      }
      state$selected[[m]] <- st
    }
    state$predictions[[paste(stage,m,sep="::")]] <- cbind(Model=m,Stage=stage,ev$predictions)
  }
  if(stage<=5L)state$fold_models[[stage]] <- stage_models else state$final_models <- stage_models
  checkpoint()
  rm(fit_cache)
  invisible(gc())
}

# ==================== 原类型CSV + 一个Excel ====================
fm <- as.data.frame(rbindlist(state$fold_metrics,fill=TRUE))
tm <- as.data.frame(rbindlist(state$test_metrics,fill=TRUE))
err <- as.data.table(error_log)
failure <- if(nrow(err))err[,.(Fail_Count=.N,Reason=paste(unique(Reason),collapse=" | ")),by=.(Model,Phase)] else
  data.table(Model=character(),Phase=character(),Fail_Count=integer(),Reason=character())
write_csv(failure,"model_failures.csv")
if(!nrow(fm) || !nrow(tm)) {checkpoint();stop("没有完整可展示结果，请查看失败名单及模型RDS。")}
cv <- as.data.table(fm)[,.(Mean_H_Score=mean(H_Score),Mean_AUC=mean(AUC),Mean_Acc=mean(Accuracy),
                           Mean_Sens=mean(Sens),Mean_Spec=mean(Spec),Mean_F1=mean(F1_Score),Mean_MCC=mean(MCC),N_Folds=uniqueN(Fold)),by=Model]
cv <- cv[N_Folds==5L][order(-Mean_AUC, -Mean_Acc, -Mean_Sens)]
cv <- merge(cv,as.data.table(tm)[,.(Model,Final_Gene_Count,Final_eccDNA_Count)],by="Model",all.x=TRUE,sort=FALSE)
cv <- cv[order(-Mean_AUC, -Mean_Acc, -Mean_Sens)]
test_table <- as.data.table(tm)[,.(Model,Mean_H_Score=H_Score,AUC,Mean_Acc=Accuracy,
                                   Mean_Sens=Sens,Mean_Spec=Spec,Mean_F1=F1_Score,Mean_MCC=MCC,Final_Gene_Count,Final_eccDNA_Count)]
test_table <- test_table[order(match(Model,cv$Model))]
result <- merge(cv[,.(Model,CV_H_Score=Mean_H_Score)],
                as.data.table(tm)[,.(Model,Validation_H_Score=H_Score,AUC,Acc=Accuracy,Sensitivity=Sens,
                                     Specificity=Spec,F1=F1_Score,MCC,Final_Gene_Count,Final_eccDNA_Count)],by="Model",sort=FALSE)
result <- as.data.frame(result[order(match(Model, cv$Model))])
candidate <- state$candidate_tables[[6]]
te <- ea[sample %in% train_ids & geneName %in% candidate$Gene]
te[,Group:=labels[sample]]
gs <- te[,.(Total_Annotation_Rows=sum(X),
            Disease_Mean_eccDNA=sum(X[Group==1])/sum(labels[train_ids]==1),
            Normal_Mean_eccDNA=sum(X[Group==0])/sum(labels[train_ids]==0),
            Disease_Mean_EA=sum(EA[Group==1])/sum(labels[train_ids]==1),
            Normal_Mean_EA=sum(EA[Group==0])/sum(labels[train_ids]==0)),by=geneName]
feature_table <- merge(candidate,as.data.frame(gs),by.x="Gene",by.y="geneName",sort=FALSE)
write_csv(feature_table,"training_candidate_features.csv")
write_csv(cv,"cross_validation_metrics.csv")
write_csv(test_table,"held_out_test_metrics.csv")
write_csv(as.data.frame(rbindlist(state$predictions,fill=TRUE)), "predictions.csv")
write_csv(as.data.frame(rbindlist(state$selected,fill=TRUE)), "selected_model_features.csv")
write_csv(as.data.frame(rbindlist(state$candidate_tables,fill=TRUE,idcol="Stage_Index")), "candidate_features_by_split.csv")
state$result <- result; state$cv_summary <- cv; state$test_summary <- test_table
state$feature_table <- feature_table
checkpoint()
