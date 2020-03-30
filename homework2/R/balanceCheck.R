
balanceCheck <- function(X, Z) {

  tsVec <- rep(NA, 8)
  pVec <- rep(NA, 8)
  ttVec <- rep(NA, 8)
  vtVec <- rep(NA, 8)

  for (i in 1:8) {
    if (i %in% c(1, 2, 3, 5)) {
      testI <- t.test(X[Z == 1, i],
                      X[Z == 0, i])

      vtVec[i] <- "continuous"

      ttVec[i] <- "t-test, difference in means"

    }

    if (i %in% c(4, 6, 7)) {
      testI <- prop.test(c(sum(X[Z == 1, i]),
                           sum(X[Z == 0, i])),
                         c(length(X[Z == 1, i]),
                           length(X[Z == 0, i])))

      vtVec[i] <- "binary"

      ttVec[i] <- "z-test, difference in proportion"

    }

    if (i == 8) {
      X[, i] <- factor(X[, i], levels = sort(unique(X[, i])))

      testtab <- table(X[, i], Z)

      testI <- chisq.test(testtab)

      vtVec[i] <- "categorical"

      ttVec[i] <- "chi-sq test of independence"

    }

    tsVec[i] <- sprintf("%3.3f", testI$statistic)
    pVec[i] <- ifelse(testI$p.value < 0.0001,
                      "< 0.0001",
                      sprintf("%1.4f", testI$p.value))

  }

  data.frame(
    variable = colnames(X),
    variable.type = vtVec,
    significance.test = ttVec,
    test.statistic = tsVec,
    p.value = pVec
  )

}

