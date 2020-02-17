
fisherExactTest <- function(W, Y, nperm = 10000, C = 0, verbose = TRUE,
                            xLab = "Test statistic") {

  require(ggplot2)

  Y[W == 1] <- Y[W == 1] - C

  rankStat <- normRank(W, Y, C = 0)

  nullrank <- rep(NA, nperm)

  for (j in 1:nperm) {
    nullrank[j] <- normRank(sample(W), Y, C = 0)
    if (verbose & (j %% 1000 == 0)) {
      cat(sprintf("permutation %i out of %i...\n", j, nperm))
    }
  }

  pvalue <- mean(nullrank > rankStat)

  permPlot <- ggplot() +
    geom_histogram(aes(nullrank, ..density..),
                   col = "snow", fill = "grey80", bins = 50) +
    geom_density(aes(nullrank), col = "black") +
    geom_vline(aes(xintercept = rankStat,
                   lty = sprintf("Obserrved test statistic (p = %1.3f)",
                                 pvalue))) +
    scale_linetype("") +
    theme(legend.position = "top") + 
    labs(x = xLab, y = "Null density") 

  list(pvalue = pvalue,
       rankStat = rankStat,
       nullrank = nullrank,
       plot = permPlot)
  
}

difTest <- function(W, Y, nperm = 10000, C = 0,
                    verbose = TRUE, xLab = "Test statistic") {

  Y[W == 1] <- Y[W == 1] - C

  difStat <- abs(mean(Y[W ==1]) - mean(Y[W == 0]))

  nulldif <- rep(NA, nperm)

  for (j in 1:nperm) {
    Wj <- sample(W)
    nulldif[j] <- abs(mean(Y[Wj ==1 ]) - mean(Y[Wj == 0]))
    if (verbose & (j %% 1000 == 0)) {
      cat(sprintf("permutation %i out of %i...\n", j, nperm))
    }
  }

  pvalue <- mean(nulldif > difStat)

  permPlot <- ggplot() +
    geom_histogram(aes(nulldif, ..density..),
                   col = "snow", fill = "grey80", bins = 50) +
    geom_density(aes(nulldif), col = "black") +
    geom_vline(aes(xintercept = difStat,
                   lty = sprintf("Obserrved test statistic (p = %1.3f)",
                                 pvalue))) + scale_linetype("") +
    theme(legend.position = "top") + 
    labs(x = xLab, y = "Null density")

  list(pvalue = pvalue,
       difStat = difStat,
       nulldif = nulldif,
       plot = permPlot)

}

