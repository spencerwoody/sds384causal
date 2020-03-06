
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

theme_set(theme_half_open())

###############################################################################
                                        #          Read and prep data         #
###############################################################################

nox <- read.csv("annualEGUs.csv", header = TRUE, stringsAsFactors = FALSE)

glimpse(nox)
