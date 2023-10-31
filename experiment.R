source("simulation.R")
###################################################################
# Please comment out the synthetic experiment in the following 
# Rfiles before the simulation test to avoid unnecessary outputs.
###################################################################
source("rldt.R")
source("ldt.R")


rate <- 0.5  # Training rate
nsample <- 50  # Sample size
niter <- 50 # Number of generated non-corrupted datasets
ntimes <- 3 # Times of repeating experiments (corruption) on a generated dataset
magnitude <- 10 # Corruption magnitude
len  <- 2
N <- 1
d <- 5
ratio <- 0.3
rldtacc_m <- c()
rldtmse_m <- c()
ldtacc_m <- c()
ldtmse_m <- c()

############################# Note #################################
# Due to the randomness of the experiments, the results can be be 
# changed when re-running. So here, to simplify reproducing the 
# results, we fix the seed to 1. 
####################################################################
set.seed(1)  # 
### Data generation
for (k in 1:niter){
  ds0 <- simdata(d, nsample, rate)
  rldtacc <- c()
  rldtmse <- c()
  ldtacc <- c()
  ldtmse <- c()
  for (i in 1:ntimes){
    ds <- corruptdata(ds0, d, ratio, magnitude, len)
    obsdata <- ds[[1]]
    truemodel <- ds[[2]]
    train_split <- ds[[3]]
    ntrain <- length(train_split)
    test_split <- ds[[4]]
    cor_split <- ds[[5]]
    alpha <- ds[[6]]
    bref <- 0*truemodel
    
    resRLDT <- expRLDT(obsdata, train_split, test_split, bref, truemodel, round(alpha*ntrain))
    resLDT <- expLDT(obsdata, train_split, test_split, bref, truemodel)
    rldtacc <- c(rldtacc, resRLDT[[2]])
    rldtmse <- c(rldtmse, resRLDT[[3]])
    ldtacc <- c(ldtacc, resLDT[[2]])
    ldtmse <- c(ldtmse, resLDT[[3]])
  }
  ldtacc_m <- c(ldtacc_m, mean(ldtacc))
  ldtmse_m <- c(ldtmse_m, mean(ldtmse))
  rldtacc_m <- c(rldtacc_m, mean(rldtacc))
  rldtmse_m <- c(rldtmse_m, mean(rldtmse))
}
print(paste("The test accuracy of R-LDT is: ", mean(rldtacc_m)))
print(paste("The mse of R-LDT is: ", mean(rldtmse_m)))
print(paste("The test accuracy of LDT is: ", mean(ldtacc_m)))
print(paste("The mse of LDT is: ", mean(ldtmse_m)))