
## totOpTime
test1_02 <- t.test(X02$totOpTime[Z02 == 1],
                   X02$totOpTime[Z02 == 0])

as.numeric(test1_02$estimate[1] - test1_02$estimate[2])

test1_02$estimate

test1_02$statistic
test1_02$p.value

## HeatInput
test2_02 <- t.test(X02$HeatInput[Z02 == 1],
                   X02$HeatInput[Z02 == 0])

as.numeric(test2_02$estimate[1] - test2_02$estimate[2])

test2_02$statistic
test2_02$p.value

## pctCapacity
test3_02 <- t.test(X02$pctCapacity[Z02 == 1],
                   X02$pctCapacity[Z02 == 0])

test3_02$statistic
test3_02$p.value

## Phase2
test4_02 <- prop.test(c(sum(X02$Phase2[Z02 == 1]), sum(X02$Phase2[Z02 == 0])),
                      c(length(X02$Phase2[Z02 == 1]), length(X02$Phase2[Z02 == 0])))

test4_02$estimate

test4_02$statistic
test4_02$p.value

## avgNOxControls
test5_02 <- t.test(X02$avgNOxControls[Z02 == 1],
                   X02$avgNOxControls[Z02 == 0])

test5_02$statistic
test5_02$p.value


## coal_no_scrubber
test6_02 <- prop.test(c(sum(X02$coal_no_scrubber[Z02 == 1]),
                        sum(X02$coal_no_scrubber[Z02 == 0])),
                      c(length(X02$coal_no_scrubber[Z02 == 1]),
                        length(X02$coal_no_scrubber[Z02 == 0])))

test6_02$estimate

test6_02$statistic
test6_02$p.value

## coal_with_scrubber
test7_02 <- prop.test(c(sum(X02$coal_with_scrubber[Z02 == 1]),
                        sum(X02$coal_with_scrubber[Z02 == 0])),
                      c(length(X02$coal_with_scrubber[Z02 == 1]),
                        length(X02$coal_with_scrubber[Z02 == 0])))

test7_02$estimate

test7_02$statistic
test7_02$p.value


## EPA.Region
tab8 <- table(nox_sub02$EPA.Region, nox_sub02$Tx)

tab8

sweep(tab8, 2, colSums(tab8), "/")

test8_02 <- chisq.test(tab8)

test8_02$statistic

test8_02$statistic

