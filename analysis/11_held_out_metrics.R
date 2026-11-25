args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript analysis/11_held_out_metrics.R MODEL_RDS OUTPUT_DIR", call. = FALSE)
}
output_dir <- args[2L]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
state <- readRDS(args[1L])
cv <- as.data.frame(state$cv_summary)
required <- c("Model", "Mean_AUC", "Mean_Acc", "Mean_Sens")
stopifnot(all(required %in% names(cv)), !anyDuplicated(cv$Model),
          all(is.finite(cv$Mean_AUC)))
cv <- cv[order(-cv$Mean_AUC, -cv$Mean_Acc, -cv$Mean_Sens), , drop = FALSE]
cv$auc_rank <- seq_len(nrow(cv))
write.csv(cv, file.path(output_dir, "cross_validation_auc_ranking.csv"), row.names = FALSE)

evaluate <- function(model) {
  pred <- as.data.frame(state$predictions[[paste(6L, model, sep = "::")]])
  stopifnot(nrow(pred) > 0L, all(c("Actual", "Probability", "Predicted") %in% names(pred)),
            all(pred$Actual %in% 0:1), all(pred$Predicted %in% 0:1),
            all(is.finite(pred$Probability)))
  labels <- pred$Actual
  scores <- pred$Probability
  n_pos <- sum(labels == 1L)
  n_neg <- sum(labels == 0L)
  stopifnot(n_pos > 0L, n_neg > 0L)

  # Ranking with average ties gives the same empirical AUC as the ROC integral.
  auc <- (sum(rank(scores)[labels == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
  sorted <- order(scores, decreasing = TRUE)
  end_of_tie <- which(!duplicated(scores[sorted], fromLast = TRUE))
  tp <- sum(labels == 1L & pred$Predicted == 1L)
  tn <- sum(labels == 0L & pred$Predicted == 0L)
  fp <- n_neg - tn
  fn <- n_pos - tp
  mcc_denominator <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  metrics <- data.frame(model = model, AUC = auc, TP = tp, TN = tn, FP = fp, FN = fn,
                        sensitivity = tp / n_pos, specificity = tn / n_neg,
                        accuracy = (tp + tn) / (n_pos + n_neg),
                        F1 = 2 * tp / (2 * tp + fp + fn),
                        MCC = if (mcc_denominator > 0) (tp * tn - fp * fn) / mcc_denominator else NA_real_)
  roc <- data.frame(model = model, score_cutoff = c(Inf, scores[sorted][end_of_tie]),
                    FPR = c(0, cumsum(labels[sorted] == 0L)[end_of_tie] / n_neg),
                    TPR = c(0, cumsum(labels[sorted] == 1L)[end_of_tie] / n_pos))
  list(metrics = metrics, roc = roc)
}

evaluated <- lapply(head(cv$Model, 10L), evaluate)
write.csv(do.call(rbind, lapply(evaluated, `[[`, "metrics")),
          file.path(output_dir, "top_10_test_metrics.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(evaluated, `[[`, "roc")),
          file.path(output_dir, "top_10_test_roc_coordinates.csv"), row.names = FALSE)
