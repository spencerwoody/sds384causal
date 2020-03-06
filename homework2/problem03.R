
###############################################################################
                                        #      Spencer Woody 12 Feb 2020      #
###############################################################################

library(dplyr)
library(ggplot2)
library(cowplot)
library(latex2exp)
library(gridExtra)

source("R/normRank.R")
source("R/fisherExactTest.R")
source("R/neyman.R")

theme_set(theme_cowplot(font_size = 16))

data(lalonde, package = "MatchIt")

glimpse(lalonde)

W <- lalonde$treat
Y <- lalonde$re78
Ygain <- lalonde$re78 - lalonde$re74

Wnohs <- W[lalonde$nodegree == 1]
Ynohs <- Y[lalonde$nodegree == 1]
Ygainnohs <- Ygain[lalonde$nodegree == 1]

Whs <- W[lalonde$nodegree == 0]
Yhs <- Y[lalonde$nodegree == 0]
Ygainhs <- Ygain[lalonde$nodegree == 0]

N <- length(W)
Nt <- sum(W == 1)
Nc <- sum(W == 0)

Ybar <- mean(Y)
Ycbar <- mean(Y[W == 0])
Ytbar <- mean(Y[W == 1])

Ysum <- sum(Y)
Ycsum <- sum(Y[W == 0])
Ytsum <- sum(Y[W == 1])


###############################################################################
                                        #          Fisher exact test          #
###############################################################################

## All 
RankTest <- fisherExactTest(W, Y, nperm = 1e4,
                            xLab = TeX("$T^{rank}$ (all)"))


RankGainTest <- fisherExactTest(W, Ygain, nperm = 1e4,
                                xLab = TeX("$T^{rank-gain}$ (all)"))

diffTest <- difTest(W, Y, nperm = 1e4,
                    xLab = TeX("$T^{dif}$ (all)"))

## No HS
RankTestNohs <- fisherExactTest(Wnohs, Ynohs, nperm = 1e4,
                                xLab = TeX("$T^{rank}$ (No high school)"))

RankGainTestNohs <- fisherExactTest(Wnohs, Ygainnohs, nperm = 1e4,
                                    xLab = TeX("$T^{rank-gain}$ (No high school)"))

diffTestNohs <- difTest(Wnohs, Ynohs, nperm = 1e4,
                        xLab = TeX("$T^{dif}$ (No high school)"))

## HS
RankTestHS <- fisherExactTest(Whs, Yhs, nperm = 1e4,
                              xLab = TeX("$T^{rank}$ (High school degree)"))

RankGainTestHS <- fisherExactTest(Whs, Ygainhs, nperm = 1e4,
                                  xLab = TeX("$T^{rank-gain}$ (High school degree)"))

diffTestHS <- difTest(Whs, Yhs, nperm = 1e4,
                      xLab = TeX("$T^{dif}$ (High school degree)"))

## Histograms

myList <- list(RankTest$plot, RankGainTest$plot, diffTest$plot,
                      RankTestNohs$plot, RankGainTestNohs$plot, diffTestNohs$plot,
                      RankTestHS$plot, RankGainTestHS$plot, diffTestHS$plot)

mygrob <- arrangeGrob(grobs = myList)

grid.arrange(mygrob)

ggsave("figures/hist-all.pdf", mygrob,
       width = 14, height = 8, units = "in")

## Fisher interval

###############################################################################
                                        #           Fisher intervals          #
###############################################################################

## Rank
CvecRank1 <- seq(-1100, -750, length.out = 15)

CvecRank2 <- seq(-10, 10, length.out = 5)

CvecRank3 <- seq(-500, 0, length.out = 5)

CvecRank <- sort(c(CvecRank1, CvecRank2, CvecRank3, 0))

PvecRank <- rep(NA, length(CvecRank))

TestListRank <- vector("list", length(CvecRank))

for (j in 1:length(CvecRank)) {
  TestListRank[[j]] <- fisherExactTest(W, Y, nperm = 5000, C= CvecRank[j],
                                       xLab = TeX("$T^{rank}$ (all)"), verbose = FALSE)
  PvecRank[j] <- TestListRank[[j]]$pvalue
  cat(j)
  cat("\n")
}


CplotRank <- ggplot() +
  geom_hline(yintercept = 0.05, lty = "dotted") + 
  geom_point(aes(CvecRank, PvecRank))

CplotRank

data.frame(C = CvecRank, P = PvecRank)

ggplotly(CplotRank)




## Gain
## Gain 500 - 3500

CvecRankGain1 <- seq(500, 3500, length.out = 25)

CvecRankGain <- sort(c(CvecRankGain1))

PvecRankGain <- rep(NA, length(CvecRankGain))

TestListRankGain <- vector("list", length(CvecRankGain))

for (j in 1:length(CvecRankGain)) {
  TestListRankGain[[j]] <- fisherExactTest(W, Ygain, nperm = 5000, C= CvecRankGain[j],
                                       xLab = TeX("$T^{rank}$ (all)"), verbose = FALSE)
  PvecRankGain[j] <- TestListRankGain[[j]]$pvalue
  cat(j)
  cat("\n")
}


CplotRankGain <- ggplot() +
  geom_hline(yintercept = 0.05, lty = "dotted") + 
  geom_point(aes(CvecRankGain, PvecRankGain))

CplotRankGain

data.frame(C = CvecRankGain, P = PvecRankGain)

ggplotly(CplotRankGain)

## Dif


CvecDif2 <- seq(-2000, 1000, length.out = 50)

CvecDif <- sort(c(CvecDif2))

PvecDif <- rep(NA, length(CvecDif))

TestListDif <- vector("list", length(CvecDif))

for (j in 1:length(CvecDif)) {
  TestListDif[[j]] <- difTest(W, Y, nperm = 5000, C = CvecDif[j],
                                       xLab = TeX("$T^{rank}$ (all)"), verbose = FALSE)
  PvecDif[j] <- TestListDif[[j]]$pvalue
  cat(j)
  cat("\n")
}

CplotDif <- ggplot() +
  geom_hline(yintercept = 0.05, lty = "dotted") + 
  geom_point(aes(CvecDif, PvecDif))

CplotDif

data.frame(C = CvecDif, P = PvecDif)

save(list=ls(), file = "fisher-tests.Rdata")

###############################################################################
                                        #                Neyman               #
###############################################################################

neyman(W, Y)

neyman(Wnohs, Ynohs)

neyman(Whs, Yhs)

###############################################################################
                                        #         Regression analysis         #
###############################################################################

## Select covariates, center them, and make a dataframe
X <- lalonde %>%
  select(age, educ, black, hispan, married, nodegree)

Xs <- scale(X, scale = FALSE) 

WXs <- cbind(W, Xs) %>% as.data.frame()


## Compute linear model with interactions
mylm <- lm(Y ~ W * ., data = WXs)

summary(mylm)


## Impute on missing potential outcomes
WmisXs <- WXs %>%
  mutate(W = 1 - W)

predObs <- predict(mylm)
predMis <- predict(mylm, newdata = WmisXs)

tauhat <- mean(W * (predObs - predMis) + (1 - W) * (predMis - predObs))

tauhat


## Only using treatment indicator
mylmW <- lm(Y ~ W, data = WXs)

summary(mylmW)

summary(mylmW)$coefficients

###############################################################################
                                        #         (Bayesian) Model            #
###############################################################################

## No. of burn in and samples
nburn <- 5000
NMC <- 10000

totsim <- nburn + NMC + 1

rho <- 1

## Prior for Gamma
## a <- 0.005 / 2
## b <- 1 / 2

a <- 0
b <- 0

## a <- 1
## b <- 1

## Prior for mu
mc <- 0
nuc <- 100^2

mt <- 0
nut <- 100^2

## Create vectors & matrix for output
mucVec <- rep(NA, totsim)
mutVec <- rep(NA, totsim)

sigma2cVec <- rep(NA, totsim)
sigma2tVec <- rep(NA, totsim)

YmisMat <- matrix(nrow = N, ncol = totsim)

ITEMat <- matrix(nrow = N, ncol = totsim)

## Initialize
mucVec[1] <- -1
mutVec[1] <- 1

YmisMat[, 1] <- rnorm(N)

sigma2cVec[1] <- 1 / rgamma(1, 1, 1)
sigma2tVec[1] <- 1 / rgamma(1, 1, 1)


## Gibbs sampler
for (i in 2:totsim) {
  
  ## Sample Ymis

  YmisMeanc <- mucVec[i - 1] +
    rho * sqrt(sigma2cVec[i - 1] * sigma2tVec[i - 1]) / sigma2tVec[i - 1] *
    (W * Y - mutVec[i - 1])

  YmisMeant <- mutVec[i - 1] +
    rho * sqrt(sigma2cVec[i - 1] * sigma2tVec[i - 1]) / sigma2cVec[i - 1] *
    ((1 - W) * Y - mucVec[i - 1])

  YmisVarc <- rep((1 - rho^2) * sigma2cVec[i - 1], N)
  YmisVart <- rep((1 - rho^2) * sigma2tVec[i - 1], N)

  meanvecI <- W * YmisMeanc + (1 - W) * YmisMeant
  varvecI <- W * YmisVarc + (1 - W) * YmisVart

  YmisMat[, i] <- rnorm(N, meanvecI, sqrt(varvecI))

  ## Sample means
  mucVar <- 1 / (1 / nuc + Nc / sigma2cVec[i - 1])
  mucMean <- mucVar * (mc / nuc + Ycsum / sigma2cVec[i - 1])

  mutVar <- 1 / (1 / nut + Nt / sigma2tVec[i - 1])
  mutMean <- mutVar * (mt / nut + Ytsum / sigma2tVec[i - 1])

  mucVec[i] <- rnorm(1, mucMean, sqrt(mucVar))
  mutVec[i] <- rnorm(1, mutMean, sqrt(mutVar))

  ## Sample variance
  RSSc <- sum((Y[W == 0] - mucVec[i])^2)
  RSSt <- sum((Y[W == 1] - mutVec[i])^2)
  
  sigma2cVec[i] <- 1 / rgamma(1, a + Nc / 2, b + RSSc / 2)
  sigma2tVec[i] <- 1 / rgamma(1, a + Nt / 2, b + RSSt / 2)

  ## Impute treatment effect
  ITEMat[, i] <- W * (Y - YmisMat[, i]) + (1 - W) * (YmisMat[, i] - Y)

  if (i %% 100 == 0) cat(sprintf("Iteration %i out of %i...\n", i, totsim - 1))

}


## Remove burnin
mucVec <- mucVec[-(1:(nburn + 1))]
mutVec <- mutVec[-(1:(nburn + 1))]

sigma2cVec <- sigma2cVec[-(1:(nburn + 1))]
sigma2tVec <- sigma2tVec[-(1:(nburn + 1))]

YmisMat <- YmisMat[, -(1:(nburn + 1))]

ITEMat <- ITEMat[, -(1:(nburn + 1))]

## Check mixing
plot(sqrt(sigma2cVec), type = "l")


## Posterior for treatment effect
tauVec <- colMeans(ITEMat)

tauEst <- mean(tauVec)
tauCI <- quantile(tauVec, c(0.025, 0.975))

tauEst
tauCI

mean(Y[W==1]) - mean(Y[W==0])

tauPlot <- ggplot() +
  geom_histogram(aes(tauVec, ..density..),
                 bins = 30, col = "snow", fill = "grey50") +
  labs(x = TeX("$\\tau$")
       ## subtitle = sprintf("Posterior mean (95%% CI): %3.1f (%3.1f, %3.1f)",
       ##                    tauEst, tauCI[1], tauCI[2])
       )

tauPlot

ggsave("figures/tauplot.pdf", tauPlot,
       width = 6, height = 4)
  
