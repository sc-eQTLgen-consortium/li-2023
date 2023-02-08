
#' Function for Rb analysis
#'
#' @param b1 Beta from first dataset.
#' @param se1 Standard error of beta from first dataset.
#' @param b2 Beta from second dataset.
#' @param se2 Standard error of beta from second dataset.
#' @param theta Variable representing sample overlap between two datasets. Should be set 0 if no sample overlap.
#'
#' @return Data frame with Rb, SE(Rb) and corresponding P-value.
#' @export
#'
#' @note This function is slightly adapted from the script shared by Ting Qi.
#'
#' @examples
calcu_cor_true <- function(b1, se1, b2, se2, theta) {
  idx <- which(is.infinite(b1) | is.infinite(b2) | is.infinite(se1) | is.infinite(se2))
  if (length(idx) > 0) {
    b1 <- b1[-idx]
    se1 <- se1[-idx]
    b2 <- b2[-idx]
    se2 <- se2[-idx]
    theta <- theta[-idx]
  }

  var_b1 <- var(b1, na.rm = T) - mean(se1^2, na.rm = T)
  var_b2 <- var(b2, na.rm = T) - mean(se2^2, na.rm = T)
  if (var_b1 < 0) {
    var_b1 <- var(b1, na.rm = T)
  }
  if (var_b2 < 0) {
    var_b2 <- var(b2, na.rm = T)
  }
  cov_b1_b2 <- cov(b1, b2, use = "complete.obs") - mean(theta, na.rm = T) * sqrt(mean(se1^2, na.rm = T) * mean(se2^2, na.rm = T))
  r <- cov_b1_b2 / sqrt(var_b1 * var_b2)

  r_jack <- c()
  n <- length(b1)
  for (k in 1:n) {
    b1_jack <- b1[-k]
    se1_jack <- se1[-k]
    var_b1_jack <- var(b1_jack, na.rm = T) - mean(se1_jack^2, na.rm = T)
    b2_jack <- b2[-k]
    se2_jack <- se2[-k]
    var_b2_jack <- var(b2_jack, na.rm = T) - mean(se2_jack^2, na.rm = T)
    if (var_b1_jack < 0) {
      var_b1_jack <- var(b1_jack, na.rm = T)
    }
    if (var_b2_jack < 0) {
      var_b2_jack <- var(b2_jack, na.rm = T)
    }
    theta_jack <- theta[-k]
    cov_e1_jack_e2_jack <- mean(theta_jack, na.rm = T) * sqrt(mean(se1_jack^2, na.rm = T) * mean(se2_jack^2, na.rm = T))
    cov_b1_b2_jack <- cov(b1_jack, b2_jack, use = "complete.obs") - cov_e1_jack_e2_jack
    r_tmp <- cov_b1_b2_jack / sqrt(var_b1_jack * var_b2_jack)
    r_jack <- c(r_jack, r_tmp)
  }
  r_mean <- mean(r_jack, na.rm = T)
  idx <- which(is.na(r_jack))
  if (length(idx) > 0) {
    se_r <- sqrt((n - 1) / n * sum((r_jack[-idx] - r_mean)^2))
  } else {
    se_r <- sqrt((n - 1) / n * sum((r_jack - r_mean)^2))
  }

  p <- pchisq((r / se_r)**2, df = 1, lower.tail = FALSE)

  res <- cbind(r, se_r, p)
  return(res)
}
