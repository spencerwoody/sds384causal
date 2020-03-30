
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(dbarts)
library(xtable)
library(gridExtra)

library(MatchIt)

theme_set(theme_half_open())

source("R/balanceCheck.R")

###############################################################################
                                        #          Read and prep data         #
###############################################################################

## Read and glimpse data
nox <- read.csv("annualEGUs.csv", header = TRUE, stringsAsFactors = FALSE)

glimpse(nox)

## Change region into a factor
nox <- nox %>%
  mutate(EPA.Region = factor(EPA.Region))

table(nox$Year)


## Subset data columns, and separate by year
nox_sub <- nox %>%
  select(Year, totOpTime, HeatInput, pctCapacity, Phase2, avgNOxControls,
         coal_no_scrubber, coal_with_scrubber, EPA.Region, Outcome, Tx) %>%
  mutate(Outcome = log(Outcome),
         Tx = as.integer(Tx * 1),
         coal_no_scrubber = as.integer(coal_no_scrubber * 1),
         coal_with_scrubber = as.integer(coal_with_scrubber * 1)) 

nox_sub02 <- nox %>%
  filter(Year == 2002) %>%
  select(Year, totOpTime, HeatInput, pctCapacity, Phase2, avgNOxControls,
         coal_no_scrubber, coal_with_scrubber, EPA.Region, Outcome, Tx) %>%
  mutate(Outcome = log(Outcome),
         Tx = as.integer(Tx * 1),
         coal_no_scrubber = as.integer(coal_no_scrubber * 1),
         coal_with_scrubber = as.integer(coal_with_scrubber * 1)) %>%
  select(-Year)

nox_sub14 <- nox %>%
  filter(Year == 2014) %>%
  select(Year, totOpTime, HeatInput, pctCapacity, Phase2, avgNOxControls,
         coal_no_scrubber, coal_with_scrubber, EPA.Region, Outcome, Tx) %>%
  mutate(Outcome = log(Outcome),
         Tx = as.integer(Tx * 1),
         coal_no_scrubber = as.integer(coal_no_scrubber * 1),
         coal_with_scrubber = as.integer(coal_with_scrubber * 1)) %>%
  select(-Year)

Z <- nox_sub$Tx
Y <- nox_sub$Outcome

Z02 <- nox_sub02$Tx
Z14 <- nox_sub14$Tx

Y02 <- nox_sub02$Outcome
Y14 <- nox_sub14$Outcome


## Covariate matrices
X <- nox_sub %>%
  select(-Year, -Outcome, -Tx)

X02 <- nox_sub02 %>%
  select(-Tx, -Outcome)

X14 <- nox_sub14 %>%
  select(-Tx, -Outcome)

###############################################################################
                                        #              Exercise 1             #
###############################################################################



## 
## 2002
##

## Linear model
crudelm02 <- lm(Y02 ~ Z02 + ., data = X02)

summary(crudelm02)

print(xtable(crudelm02,
             caption = "``Crude'' linear regression model for estimating of treatment effect using 2002 data",
             label = "tab:crude-lm-02", align = "lrrrr"), booktabs = TRUE)


## 2014
crudelm14 <- lm(Y14 ~ Z14 + ., data = X14)

summary(crudelm14)

crudelmtab <- xtable(crudelm14)

print(xtable(crudelm14,
             caption = "``Crude'' linear regression model for estimating of treatment effect using 2014 data",
             label = "tab:crude-lm-14", align = "lrrrr"), booktabs = TRUE)

## 2002

## Check covariate balance

glimpse(X02)

covbal02 <- balanceCheck(X02, Z02)

covbal14 <- balanceCheck(X14, Z14)

print(xtable(covbal02,
             caption = "Covariate balance check for crude analysis, year 2002",
             label = "tab-bal1-02", align = "llllrr"),
      include.rownames=FALSE,
      booktabs = TRUE)

print(xtable(covbal14,
             caption = "Covariate balance check for crude analysis, year 2014",
             label = "tab-bal1-14", align = "llllrr"),
      include.rownames=FALSE,
      booktabs = TRUE)


###############################################################################
                                        #              Exercise 2             #
###############################################################################

## 
## 
## Part A
##
##

logistic02 <- glm(Z02 ~ ., data = X02, family=binomial)
prop02 <- logistic02$fitted.values

logistic14 <- glm(Z14 ~ ., data = X14, family=binomial)
prop14 <- logistic14$fitted.values

propDf02 <- data.frame(Z = Z02, Y = Y02, prop = prop02, Year = 2002)
propDf14 <- data.frame(Z = Z14, Y = Y14, prop = prop14, Year = 2014)

propDf <- rbind(propDf02, propDf14) %>%
  mutate(Year = factor(Year),
         Z = factor(Z)) 

propPlot <- ggplot(propDf) +
  geom_histogram(aes(prop, fill = Z), col = "snow") +
  facet_grid(rows = vars(Z), cols = vars(Year), labeller = label_both) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5)) +
  scale_fill_brewer(palette = "Dark2") + 
  labs(x = "Propensity score, estimated using logistic regression")

propPlot

ggsave("figures/prop-count-logistic.pdf", propPlot,
       width = 6, height = 4, units = "in")



## ggsave("")

propDensPlot <- ggplot(propDf) +
  geom_histogram(aes(prop, ..density.., fill = Z), col = "snow") +
  facet_grid(rows = vars(Z), cols = vars(Year), labeller = label_both) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5)) +
  scale_fill_brewer(palette = "Dark2") + 
  labs(x = "Propensity score, estimated using logistic regression") 

propDensPlot

ggsave("figures/prop-density-logistic.pdf", propDensPlot,
       width = 6, height = 4, units = "in")

propGrob <- arrangeGrob(propPlot, propDensPlot, ncol = 2)

grid.arrange(propGrob)

ggsave("figures/propensity-logistic-1row.pdf", propGrob,
       width = 11.6, height = 5, units = "in")

##
## 
## Part B
##
##

set.seed(1)

matchB02 <- matchit(Z02 ~ totOpTime + HeatInput + pctCapacity +
                      Phase2 + avgNOxControls + coal_no_scrubber +
                      coal_with_scrubber + EPA.Region, distance = prop02,
                    data = X02, method = "nearest")

summary(matchB02)

difvecB02 <- (Y02[as.numeric(row.names(matchB02$match.matrix))] -
              Y02[as.numeric(matchB02$match.matrix)])

mean(difvecB02)

Y02B <- Y02[c(as.numeric(row.names(matchB02$match.matrix)),
              as.numeric(matchB02$match.matrix))]

Z02B <- Z02[c(as.numeric(row.names(matchB02$match.matrix)),
              as.numeric(matchB02$match.matrix))]

X02B <- X02[c(as.numeric(row.names(matchB02$match.matrix)),
              as.numeric(matchB02$match.matrix)), ]


length(Z02B)

lm02B <- lm(Y02B ~ Z02B + ., data = X02B)
summary(lm02B)

print(xtable(lm02B,
             caption = "Linear regression model for estimating of treatment effect using 2002 data, after using one-to-one matching (Problem 2(b))",
             label = "tab:lm-2b-02", align = "lrrrr"), booktabs = TRUE)

## Overlap
covbal02B <- balanceCheck(X02B, Z02B)

print(xtable(covbal02B,
             caption = "Covariate balance check for one-to-one propensity score matching (Problem 2(b)), year 2002",
             label = "tab-bal2b-02", align = "llllrr"),
      include.rownames=FALSE,
      booktabs = TRUE)




## 
## 2014
## 

matchB14 <- matchit(Z14 ~ totOpTime + HeatInput + pctCapacity +
                      Phase2 + avgNOxControls + coal_no_scrubber +
                      coal_with_scrubber + EPA.Region, distance = prop14,
                    data = X14, method = "nearest")

summary(matchB14)

difvecB14 <- (Y14[as.numeric(row.names(matchB14$match.matrix))] -
              Y14[as.numeric(matchB14$match.matrix)])

mean(difvecB14)

Y14B <- Y14[c(as.numeric(row.names(matchB14$match.matrix)),
              as.numeric(matchB14$match.matrix))]

Z14B <- Z14[c(as.numeric(row.names(matchB14$match.matrix)),
              as.numeric(matchB14$match.matrix))]

X14B <- X14[c(as.numeric(row.names(matchB14$match.matrix)),
              as.numeric(matchB14$match.matrix)), ]

lm14B <- lm(Y14B ~ Z14B + ., data = X14B)
summary(lm14B)


print(xtable(lm14B,
             caption = "Linear regression model for estimating of treatment effect using 2014 data, after using one-to-one matching (Problem 2(b))",
             label = "tab:lm-2b-14", align = "lrrrr"), booktabs = TRUE)

## Overlap
covbal14B <- balanceCheck(X14B, Z14B)

print(xtable(covbal14B,
             caption = "Covariate balance check for one-to-one propensity score matching (Problem 2(b)), year 2014",
             label = "tab-bal2b-14", align = "llllrr"),
      include.rownames=FALSE,
      booktabs = TRUE)

X[row.names(matchB14$match.matrix), 1]
X[matchB14$match.matrix, 1]

## 
## Part C
##

## 2002
set.seed(1)

matchC02 <- matchit(Z02 ~ totOpTime + HeatInput + pctCapacity +
                      Phase2 + avgNOxControls + coal_no_scrubber +
                      coal_with_scrubber + EPA.Region,
                    data = X02, method = "nearest", caliper = 0.1 * sd(prop02))

summary(matchC02)


matchC02$match.matrix <- matchC02$match.matrix %>% na.omit()

difvecC02 <- (Y02[as.numeric(row.names(matchC02$match.matrix))] -
              Y02[as.numeric(matchC02$match.matrix)])

mean(difvecC02)

Y02C <- Y02[c(as.numeric(row.names(matchC02$match.matrix)),
              as.numeric(matchC02$match.matrix))]

Z02C <- Z02[c(as.numeric(row.names(matchC02$match.matrix)),
              as.numeric(matchC02$match.matrix))]

X02C <- X02[c(as.numeric(row.names(matchC02$match.matrix)),
              as.numeric(matchC02$match.matrix)), ]


length(Z02C)

lm02C <- lm(Y02C ~ Z02C + ., data = X02C)
summary(lm02C)

print(xtable(lm02C,
             caption = "Linear regression model for estimating of treatment effect using 2002 data, after using one-to-one matching with a caliber (Problem 2(c))",
             label = "tab:lm-2c-02", align = "lrrrr"), booktabs = TRUE)

## Overlap
covbal02C <- balanceCheck(X02C, Z02C)

print(xtable(covbal02C,
             caption = "Covariate balance check for one-to-one propensity score matching with a caliber (Problem 2(c)), year 2002",
             label = "tab-bal2c-02", align = "llllrr"),
      include.rownames=FALSE,
      booktabs = TRUE)


## 2014
matchC14 <- matchit(Z14 ~ totOpTime + HeatInput + pctCapacity +
                      Phase2 + avgNOxControls + coal_no_scrubber +
                      coal_with_scrubber + EPA.Region,
                    data = X14, method = "nearest", caliper = 0.1 * sd(prop14))

summary(matchC14)


matchC14$match.matrix <- matchC14$match.matrix %>% na.omit()

difvecC14 <- (Y14[as.numeric(row.names(matchC14$match.matrix))] -
              Y14[as.numeric(matchC14$match.matrix)])

mean(difvecC14)

Y14C <- Y14[c(as.numeric(row.names(matchC14$match.matrix)),
              as.numeric(matchC14$match.matrix))]

Z14C <- Z14[c(as.numeric(row.names(matchC14$match.matrix)),
              as.numeric(matchC14$match.matrix))]

X14C <- X14[c(as.numeric(row.names(matchC14$match.matrix)),
              as.numeric(matchC14$match.matrix)), ]


length(Z14C)

lm14C <- lm(Y14C ~ Z14C + ., data = X14C)
summary(lm14C)

print(xtable(lm14C,
             caption = "Linear regression model for estimating of treatment effect using 2014 data, after using one-to-one matching with a caliber (Problem 2(c))",
             label = "tab:lm-2c-14", align = "lrrrr"), booktabs = TRUE)

## Overlap
covbal14C <- balanceCheck(X14C, Z14C)

print(xtable(covbal14C,
             caption = "Covariate balance check for one-to-one propensity score matching with a caliber (Problem 2(c)), year 2014",
             label = "tab-bal2c-14", align = "llllrr"),
      include.rownames=FALSE,
      booktabs = TRUE)

## 2014

set.seed(1)
matchC14 <- matchit(Z14 ~ totOpTime + HeatInput + pctCapacity +
                      Phase2 + avgNOxControls + coal_no_scrubber +
                      coal_with_scrubber + EPA.Region,
                    data = X14, method = "nearest", distance = prop14,
                    caliper = 0.1 * sd(prop14))

names(matchC02)

matchC02$match.matrix

summary(matchC02)


## 
## Part D
##

## 2002
matchD02 <- matchit(Z02 ~ totOpTime + HeatInput + pctCapacity +
                      Phase2 + avgNOxControls + coal_no_scrubber +
                      coal_with_scrubber + EPA.Region,
                    data = X02, method = "subclass", distance = prop02, subclass = 4)

summary(matchD02)

names(matchD02)

matchD02$subclass %>% table()

range(prop02)

range(prop02[matchD02$subclass == 1])
range(prop02[matchD02$subclass == 2])
range(prop02[matchD02$subclass == 3])
range(prop02[matchD02$subclass == 4])

numsubclass02 <- matchD02$subclass %>% unique() %>% length()

subdif02 <- rep(NA, numsubclass02)
subnum02 <- rep(NA, numsubclass02)

covbalList02 <- vector("list", numsubclass02)

TEvec02 <- rep(NA, numsubclass02)
SEvec02 <- rep(NA, numsubclass02)

for (i in 1:numsubclass02) {
  Y02sub <- Y02[matchD02$subclass == i]
  Z02sub <- Z02[matchD02$subclass == i]
  X02sub <- X02[matchD02$subclass == i, ]

  lm02sub <- lm(Y02sub ~ Z02sub + ., data = X02sub)

  TEvec02[i] <- summary(lm02sub)$coefficients[2, 1] 
  SEvec02[i] <- summary(lm02sub)$coefficients[2, 2] 

  covbalList02[[i]] <- balanceCheck(X02sub, Z02sub)

  subnum02[i] <- length(Y02sub)
  subdif02[i] <- mean(Y02sub[Z02sub == 1]) - mean(Y02sub[Z02sub == 0])
}


wj02 <- subnum02 / sum(subnum02)

sum(TEvec02 * wj02)
sqrt(sum(SEvec02^2 * wj02))

covbalList02[[i]]

covbalDf02 <- data.frame(variable = covbalList02[[i]][, "variable"])


for (i in 1:numsubclass02) {
  covbalDf02 <- cbind(covbalDf02, covbalList02[[i]][, 5])
}

colnames(covbalDf02)[-1] <- paste0("subgroup", 1:numsubclass02)

print(xtable(covbalDf02,
             caption = "Covariate balance check for subclassification using propensity score (4 subclasses), year 2002",
             label = "tab-bal2d-02", align = "llrrrr"),
      include.rownames=FALSE,
      booktabs = TRUE)



subdif02
sum(subdif02 * subnum02 / sum(subnum02))

## 2014
matchD14 <- matchit(Z14 ~ totOpTime + HeatInput + pctCapacity +
                      Phase2 + avgNOxControls + coal_no_scrubber +
                      coal_with_scrubber + EPA.Region,
                    data = X14, method = "subclass", distance = prop14, subclass = 4)

summary(matchD14)

names(matchD14)

matchD14$subclass %>% table()

range(prop14)

range(prop14[matchD14$subclass == 1])
range(prop14[matchD14$subclass == 2])
range(prop14[matchD14$subclass == 3])
range(prop14[matchD14$subclass == 4])

numsubclass14 <- matchD14$subclass %>% unique() %>% length()

subdif14 <- rep(NA, numsubclass14)
subnum14 <- rep(NA, numsubclass14)

covbalList14 <- vector("list", numsubclass14)

TEvec14 <- rep(NA, numsubclass14)
SEvec14 <- rep(NA, numsubclass14)

for (i in 1:numsubclass14) {
  Y14sub <- Y14[matchD14$subclass == i]
  Z14sub <- Z14[matchD14$subclass == i]
  X14sub <- X14[matchD14$subclass == i, ]

  lm14sub <- lm(Y14sub ~ Z14sub + ., data = X14sub)

  TEvec14[i] <- summary(lm14sub)$coefficients[2, 1] 
  SEvec14[i] <- summary(lm14sub)$coefficients[2, 2] 

  covbalList14[[i]] <- balanceCheck(X14sub, Z14sub)

  subnum14[i] <- length(Y14sub)
  subdif14[i] <- mean(Y14sub[Z14sub == 1]) - mean(Y14sub[Z14sub == 0])
}


wj14 <- subnum14 / sum(subnum14)

sum(TEvec14 * wj14)
sqrt(sum(SEvec14^2 * wj14))

covbalList14[[i]]

covbalDf14 <- data.frame(variable = covbalList14[[i]][, "variable"])


for (i in 1:numsubclass14) {
  covbalDf14 <- cbind(covbalDf14, covbalList14[[i]][, 5])
}

colnames(covbalDf14)[-1] <- paste0("subgroup", 1:numsubclass14)

print(xtable(covbalDf14,
             caption = "Covariate balance check for subclassification using propensity score (4 subclasses), year 2014",
             label = "tab-bal2d-14", align = "llrrrr"),
      include.rownames=FALSE,
      booktabs = TRUE)

matchD14 <- matchit(Z14 ~ totOpTime + HeatInput + pctCapacity +
                      Phase2 + avgNOxControls + coal_no_scrubber +
                      coal_with_scrubber + EPA.Region,
                    data = X14, method = "subclass", distance = prop14, subclass = 4)

summary(matchD14)

names(matchD14)

matchD14$subclass %>% table()

range(prop14)

range(prop14[matchD14$subclass == 1])
range(prop14[matchD14$subclass == 2])
range(prop14[matchD14$subclass == 3])
range(prop14[matchD14$subclass == 4])

numsubclass14 <- matchD14$subclass %>% unique() %>% length()

subdif14 <- rep(NA, numsubclass14)
subnum14 <- rep(NA, numsubclass14)

for (i in 1:numsubclass14) {
  Y14sub <- Y14[matchD14$subclass == i]
  Z14sub <- Z14[matchD14$subclass == i]
  subnum14[i] <- length(Y14sub)
  subdif14[i] <- mean(Y14sub[Z14sub == 1]) - mean(Y14sub[Z14sub == 0])
}

subdif14
sum(subdif14 * subnum14 / sum(subnum14))

## Part E



ipwE02 <- ipw::ipwpoint(Tx, family = "binomial", link = "logit",
                       denominator = ~totOpTime + HeatInput + pctCapacity +
                         Phase2 + avgNOxControls + coal_no_scrubber +
                         coal_with_scrubber + factor(EPA.Region),
                       data = nox_sub02)

mean((stabweights02 - myipwF02$ipw.weights)^2)

weights02 <- Z02 / (prop02) + (1 - Z02) / (1 - prop02)
weights14 <- Z14 / (prop14) + (1 - Z14) / (1 - prop14)


ipwE02$ipw.weights

ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(ipwE02$ipw.weights, weights02))


lmE02 <- lm(Y02 ~ Z02 + ., data = X02, weights = weights02)


print(xtable(lmE02,
             caption = "Linear regression model with inverse probability weights for estimating of treatment effect using 2002 data (Problem 2(e))",
             label = "tab:lm-2e-02", align = "lrrrr"), booktabs = TRUE)


lmE14 <- lm(Y14 ~ Z14 + ., data = X14, weights = weights14)


print(xtable(lmE14,
             caption = "Linear regression model with inverse probability weights for estimating of treatment effect using 2014 data (Problem 2(e))",
             label = "tab:lm-2e-14", align = "lrrrr"), booktabs = TRUE)

lmF02 <- lm(Y02 ~ Z02, weights = stabweights02)

summary(lmE02)
summary(lmF02)



library(survey)

myglmE02 <- survey::svyglm(Outcome ~ Tx ,
                           design = svydesign(~1, weights = ~ myipwE02$ipw.weights, data = nox_sub02))

summary(myglmE02)

plot(myipw$ipw.weights, weights02)

?fisher.test

plot(sort(weights02))

weightsDf <- rbind(
  data.frame(Weights = weights02, Z = Z02, Year = "2002"),
  data.frame(Weights = weights14, Z = Z14, Year = "2014")
)

ggplot(weightsDf) +
  geom_histogram(aes(Weights, col = Z)) +
  facet_grid(Z ~ Year)

## Part F

myipwF02   <- ipw::ipwpoint(Tx, family = "binomial", link = "logit",
                            numerator = ~ 1,
                       denominator = ~totOpTime + HeatInput + pctCapacity +
                         Phase2 + avgNOxControls + coal_no_scrubber +
                         coal_with_scrubber + factor(EPA.Region),
                       data = nox_sub02)

myipwE02$ipw.weights
myipwF02$ipw.weights



myipw$ipw.weights

myglmF02 <- survey::svyglm(Outcome ~ Tx ,
                           design = svydesign(~1, weights = ~ myipwF02$ipw.weights, data = nox_sub02))

library(mgcv)

summary(myglmE02)
summary(myglmF02)

Zmean02 <- mean(Z02)
Zmean14 <- mean(Z14)

stabweights02 <- Z02 * Zmean02 / (prop02) + (1 - Z02) * (1 - Zmean02) / (1 - prop02)
stabweights14 <- Z14 * Zmean14 / (prop14) + (1 - Z14) * (1 - Zmean14) / (1 - prop14)

lmF02 <- lm(Y02 ~ Z02 + ., data = X02, weights = stabweights02)


print(xtable(lmF02,
             caption = "Linear regression model with inverse probability weights for estimating of treatment effect using 2002 data (Problem 2(e))",
             label = "tab:lm-2e-02", align = "lrrrr"), booktabs = TRUE)


lmF14 <- lm(Y14 ~ Z14 + ., data = X14, weights = stabweights14)


print(xtable(lmF14,
             caption = "Linear regression model with inverse probability weights for estimating of treatment effect using 2014 data (Problem 2(e))",
             label = "tab:lm-2e-14", align = "lrrrr"), booktabs = TRUE)



var(weights02)
var(stabweights02)

plot(sort(stabweights02))
plot(sort(stabweights02))

weightsDf <- rbind(
  data.frame(Weights = stabweights02,
             Z = Z02, Year = "2002",
             type = "stabilized") %>% arrange(Weights) %>% mutate(Order = 1:n()),
  data.frame(Weights = stabweights14,
             Z = Z14, Year = "2014",
             type = "stabilized") %>% arrange(Weights) %>% mutate(Order = 1:n()),
    data.frame(Weights = weights02,
             Z = Z02, Year = "2002",
             type = "unstabilized") %>% arrange(Weights) %>% mutate(Order = 1:n()),
  data.frame(Weights = weights14,
             Z = Z14, Year = "2014",
             type = "unstabilized") %>% arrange(Weights) %>% mutate(Order = 1:n())
) %>%
  mutate(type = factor(type, levels = c("unstabilized", "stabilized")))



weightsPlot <- ggplot(weightsDf) +
  ## geom_point(aes(Order, Weights, col = factor(Z))) +
  geom_jitter(aes(Order, Weights, col = factor(Z)), width = 4, pch = 20) +
  facet_grid(rows = vars(Year), cols = vars(type), scales = "free_y") +
  scale_color_brewer("Z", palette = "Set1") + 
  theme(legend.position = "top")

weightsPlot

ggsave("figures/weights-all.pdf", weightsPlot, width = 10, height = 5)

plot(weights02, stabweights02)

qplot(1:length(weights02), sort(weights02)) +
  labs(x = "Order", y = "weight")

qplot(1:length(stabweights02), sort(stabweights02)) +
  labs(x = "Order", y = "weight")

plot(1:length(stabweights02), sort(stabweights02))

ggplot(stabweightsDf) +
  geom_histogram(aes(Weights, fill = factor(Z))) +
  facet_grid(Z ~ Year)

ggplot(stabweightsDf) +
  geom_point(aes(Order, Weights, col = Z)) +
  facet_grid(Z ~ Year)


###############################################################################
                                        #              Exercise 4             #
###############################################################################


Xmat02 <- makeModelMatrixFromDataFrame(X02)

Z02bart <- bart(x.train = Xmat02, y.train = Z02,
                nskip = 5000, ndpost = 1000)

Zhat02 <- pnorm(colMeans(Z02bart$yhat.train))


Xmat14 <- makeModelMatrixFromDataFrame(X14)

Z14bart <- bart(x.train = Xmat14, y.train = Z14,
                nskip = 5000, ndpost = 1000)

Zhat14 <- pnorm(colMeans(Z14bart$yhat.train))

propbart02 <- Zhat02
propbart14 <- Zhat14

weightsbart02 <- Z02 / (propbart02) + (1 - Z02) / (1 - propbart02)

weightsbart14 <- Z14 / (propbart14) + (1 - Z14) / (1 - propbart14)

lmG02 <- lm(Y02 ~ Z02 + ., data = X02, weights = weightsbart02)


print(xtable(lmG02,
             caption = "Linear regression model with inverse probability weights (estimated using BART) for estimating of treatment effect using 2002 data (Problem 4)",
             label = "tab:lm-2g-02", align = "lrrrr"), booktabs = TRUE)


lmG14 <- lm(Y14 ~ Z14 + ., data = X14, weights = weightsbart14)


print(xtable(lmG14,
             caption = "Linear regression model with inverse probability weights (estimated using BART) for estimating of treatment effect using 2014 data (Problem 4)",
             label = "tab:lm-2g-14", align = "lrrrr"), booktabs = TRUE)

plot(1:length(weights02), sort(weights02))
plot(1:length(weightsbart02), sort(weightsbart02))

plot(1:length(weights14), sort(weights14))
plot(1:length(weightsbart14), sort(weightsbart14))


var(weights02)
var(weightsbart02)

var(weights14)
var(weightsbart14)



propDfBART <- rbind(propDf02 %>%
                    mutate(prop_bart = Zhat02),
                    propDf14 %>%
                    mutate(prop_bart = Zhat14)) %>%
  rename(prop_logistic = prop) %>%
  mutate(Year = factor(Year),
         Z = factor(Z)) 

head(propDfBART)

propDfCat <- propDfBART %>%
  pivot_longer(-Z, -Y, -Year, cols = starts_with("prop"),
               names_to = "prop_model", values_to = "prop_score")

propDfCat <- propDfBART %>%
  pivot_longer(cols = starts_with("prop"), names_to = "prop_model", names_prefix = "prop_",
               values_to = "prop_score") %>%
  mutate(prop_model = factor(prop_model, levels = c("logistic", "bart")))

glimpse(propDfCat)

## propDfCat %>%
##   filter(Year == "2002") %>%
##   ggplot() + 
##   geom_histogram(aes(prop_score, fill = Z), col = "snow", bins = 40) +
##   facet_grid(cols = vars(Z), rows = vars(prop_model), labeller = label_both) +
##   theme(legend.position = "none",
##         axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5)) +
##   scale_fill_brewer(palette = "Dark2") + 
##   labs(title = "Year 2002", x = "Propensity score")

propDfCat %>%
  filter(Year == "2002") %>%
  ggplot() + 
  geom_histogram(aes(prop_score, ..density.., fill = Z), col = "snow", bins = 40) +
  facet_grid(cols = vars(Z), rows = vars(prop_model), labeller = label_both) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5)) +
  scale_fill_brewer(palette = "Dark2") + 
  labs(title = "Year 2002", x = "Propensity score")

## propComp14 <- propDfCat %>%
##   filter(Year == "2014") %>%
##   ggplot() + 
##   geom_histogram(aes(prop_score, fill = Z), col = "snow", bins = 40) +
##   facet_grid(cols = vars(Z), rows = vars(prop_model), labeller = label_both) +
##   theme(legend.position = "none",
##         axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5)) +
##   scale_fill_brewer(palette = "Dark2") + 
##   labs(title = "Year 2014", x = "Propensity score")



propDfCat %>%
  filter(Year == "2014") %>%
  ggplot() + 
  geom_histogram(aes(prop_score, ..density.., fill = Z), col = "snow", bins = 40) +
  facet_grid(cols = vars(Z), rows = vars(prop_model), labeller = label_both) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5)) +
  scale_fill_brewer(palette = "Dark2") + 
  labs(title = "Year 2014", x = "Propensity score")

propComp <- ggplot(propDfBART) +
  geom_abline(intercept = 0, slope = 1, col = "black") + 
  geom_point(aes(prop_logistic, prop_bart, col = prop_bart > prop_logistic)) +
  facet_grid(rows = vars(Z), cols = vars(Year), labeller = label_both) +
  scale_color_brewer(palette = "Set1") + 
  theme(legend.position = "top") +
  coord_equal()

propComp

ggsave("figures/prop-comp.pdf", propComp, width = 9, height = 6)



ggplot(propDfBART) +
  geom_vline(xintercept = 0) + 
  ## geom_abline(intercept = 0, slope = 1, col = "black") + 
  geom_histogram(aes(prop - prop_bart), bins = 50) +
  facet_grid(rows = vars(Z), cols = vars(Year), labeller = label_both)


propPlotBART <- ggplot(propDfBART) +
  geom_histogram(aes(prop_bart, fill = Z), col = "snow") +
  facet_grid(rows = vars(Z), cols = vars(Year), labeller = label_both) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5)) +
  scale_fill_brewer(palette = "Dark2") + 
  labs(x = "Propensity score, estimated using BART")

propPlotBART

propDensPlotBART <- ggplot(propDfBART) +
  geom_histogram(aes(prop_bart, ..density.., fill = Z), col = "snow") +
  facet_grid(rows = vars(Z), cols = vars(Year), labeller = label_both) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5)) +
  scale_fill_brewer(palette = "Dark2") + 
  labs(x = "Propensity score, estimated using BART") 

propDensPlotBART

propGrobBART <- arrangeGrob(propPlotBART, propDensPlotBART, ncol = 2)

grid.arrange(propGrobBART)

ggsave("figures/propensity-bart-1row.pdf", propGrobBART,
       width = 11.6, height = 5, units = "in")
