#' run_survivalROC
#'
#' @param data 
#' @param time_point 
#' @param risk_col 
#'
#' @returns
#' @export
#'
#' @examples
run_survivalROC <- function(data, time_point, risk_col) {
  # Prepare inputs
  surv_time <- data$W
  surv_event <- data$delta
  marker <- data[[risk_col]]
  
  # Try-catch for convergence
  out <- tryCatch({
    roc <- survivalROC(Stime = surv_time,
                       status = surv_event,
                       # risk scores
                       marker = marker,
                       # follow-up time
                       predict.time = time_point,
                       # censoring-adjustment method
                       method = "KM")
    
    # Youden's J
    j_stats <- roc$TP - roc$FP
    max_idx <- which.max(j_stats)
    
    # Store follow-up, optimal cut-point, sens., spec., and AUC at optimal cut-point
    tibble(
      time = time_point,
      cut = roc$cut.values[max_idx],
      sens = roc$TP[max_idx],
      spec = 1 - roc$FP[max_idx],
      auc  = roc$AUC %||% NA_real_
    )
  }, error = function(e) {
    tibble(time = time_point, cut = NA, sens = NA, spec = NA, auc = NA_real_)
  })
  
  return(out)
}

#' sample_size_binary
#'
#' @param p1 
#' @param effect_size 
#' @param alpha 
#' @param power 
#'
#' @returns
#' @export
#'
#' @examples
sample_size_binary <- function(p1, effect_size = 0.3, alpha = 0.05, power = 0.8) {
  z_alpha <- qnorm(1 - alpha / 2)
  z_beta  <- qnorm(power)
  p2 <- p1 * (1 - effect_size)
  p_bar <- (p1 + p2) / 2
  n <- 2 * (z_alpha + z_beta)^2 * p_bar * (1 - p_bar) / (p1 - p2)^2
  return(ceiling(n))
}
