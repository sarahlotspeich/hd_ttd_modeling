#' compute_uno_c
#'
#' @param test_data 
#' @param G_hat 
#' @param t_star 
#'
#' @returns
#' @export
#'
#' @examples
compute_uno_c <- function(test_data, G_hat, t_star) {
  # Extract min times, event indicators, linear predictors, and sample size
  W  <- test_data$W
  d  <- test_data$delta
  lp <- test_data$lp
  n  <- length(W)
  
  # Create all pairwise combinations, dropping i == j
  pair_df <- expand.grid(i = 1:n, j = 1:n) %>%
    filter(i != j)
  
  # Extract values for each pair
  Ti <- W[pair_df$i]
  Tj <- W[pair_df$j]
  di <- d[pair_df$i]
  lpi <- lp[pair_df$i]
  lpj <- lp[pair_df$j]
  
  # Identify valid comparable pairs (Equation (6) in Uno et al. (2011))
  valid_pairs <- which(Ti < Tj & Ti <= t_star & di == 1)
  if (length(valid_pairs) == 0) return(NA_real_)
  
  # Calculate inverse probability of censoring weights
  G_Ti <- G_hat(Ti[valid_pairs])
  wi <- 1 / (G_Ti^2)
  
  # Concordance: a pair is concordant if the LP[i] > LP[j]
  concordant <- as.numeric(lpi[valid_pairs] > lpj[valid_pairs])
  
  # Compute Uno's C (Equation (6) in Uno et al. (2011))
  sum(wi * concordant) / sum(wi)
}
