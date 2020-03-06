
neyman <- function(W, Y, alpha = 0.05) {

  tstar <- qnorm(1 - alpha / 2)
  
  Y0 <- Y[W == 0]
  Y1 <- Y[W == 1]

  Y0mean <- mean(Y0)
  Y1mean <- mean(Y1)

  N0 <- length(Y0)
  N1 <- length(Y1)

  s20 <- sum((Y0 - Y0mean)^2) / (N0 - 1)
  s21 <- sum((Y1 - Y1mean)^2) / (N1 - 1)

  est <- Y1mean - Y0mean
  
  vhat <- s20 / N0 + s21 / N1

  se <- sqrt(vhat)

  CI <- c(est - tstar * se, est + tstar * se)

  list(est = est,
       se = se,
       CI = CI,
       alpha = alpha)

}
