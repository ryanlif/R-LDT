library(MASS)
library(gurobi)
library(Matrix)
library(glmnetUtils)
library(gridExtra)
library(xlsx)
source("ldt.R")

update_beta<-function(X, Pert, R, Y, lambda, bref){
  beta <- LDT(X, Pert, R, Y, lambda, bref)
  beta[is.na(beta)]<-0
  update_beta <- beta
}

update_pert <- function(beta, pert, X, R, Y, St, lambda2){
  t <- rep(0, nrow(X))
  for (j in 1:nrow(X)){
    if (Y[j] ==1){
      t[j] <- min(R[j]-sum(X[j,]*beta)-Y[j],0)
    }else{
      t[j] <- max(R[j]-sum(X[j,]*beta)-Y[j],0)
    }
  }
  update_pert <- t
}

HT <- function(ls, k){
  org <- 1:length(ls)
  reorg <- sample(org, length(org))
  ps <- ls[reorg]
  ord <- order(ps, decreasing=FALSE)
  pick <- ord[1:k]
  HT <- reorg[pick]
}

RLDT <- function(X, R, Y, kcor, lambda1, lambda2,bref){
  d <- ncol(X)
  nu <- nrow(X)
  beta <- rep(0, d)
  optval <- c()
  eps <- 0.1
  t <- 0
  lr0 <- 0.5
  St<- 1:nu
  lr <- lr0
  pert <- R - Y
  prepert <- rep(0, nu)
  ncor <- nu - kcor
  St <- 1:nu
  beta <- bref
  setc <- c()
  setc <- c(setc, list(beta))
  ### Alternating optimization
  while (sqrt(sum(pert[St]^2)/length(St))>eps && t<=100){
    al <- log(t+2)/2
    lr <-1/(al*sqrt(t+1)) * lr0
    print(paste("iteration " ,t,"...", sep=""))
    beta <- update_beta(X[St,], 0*prepert[St], R[St], Y[St], lambda1*lr, beta)
    pert <- update_pert(beta, pert, X, R, Y, St, lambda2)
    St <- HT(pert, ncor)
    setc <- c(setc, list(beta))
    prepert <- pert
    t <- t + 1
  }
  beta <- update_beta(X[St,], 0*prepert[St], R[St], Y[St], lambda1, bref)
  setc <- c(setc, list(beta))
  rpert <- rep(0, length(beta))
  for (dims in 1:length(beta)){
    rpert <- R - beta[dims]*X[,dims]
  }
  res <- list(beta, rpert, setc)
  return(res)
}

expRLDT <- function(obdata, train_split, test_split, bref, truemodel, kcor){
  perturb <- c()
  score <- 0
  train_acc <- c()
  test_acc <- c()
  abe_acc <- c()
  mse_acc <- c()
  data <- obdata
  d <- ncol(data) - 4
  nu <- nrow(data)
  X <- data[, 1:d]
  print(data)
  ones <- rep(1,nu)
  R <- data$r
  Y <- data$y
  Yt <- data$ytrue
  full <- 1:nu
  St<-train_split
  Se<-test_split
  ntrain <- length(St) 
  ntest <- length(Se)
  beta <- rep(0, d+1)
  eps <- 0.1
  lambda1 <- 1
  lambda2 <- 10
  res <- RLDT(X[St,], R[St], Y[St], kcor, lambda1,lambda2, bref)
  beta <- res[[1]]
  pert <- res[[2]]
  betac <- res[[3]]
  mse <- mean((beta-truemodel)^2) # MSE
  ytr <- 2*(R[St]> as.matrix(X[St,]) %*% beta)-1
  train_acc <- c(train_acc, sum(Yt[St]==ytr)/ntrain) # Training accuracy
  yte <- 2*(R[Se] > as.matrix(X[Se,]) %*% beta)-1
  test_acc <- c(test_acc, sum(Yt[Se]==yte)/ntest) # Testing accuracy
  expRLDT <- list(train_acc, test_acc, mse, beta)
}

################ Synthetic case ################
### True model: beta = (0.4, 0.6)
truemodel<-c(0.4,0.6)

## Toy 1
data1<-read.xlsx("synthetic/toy1.xlsx", sheetIndex=1, header=TRUE)
obsdata<-data1[,1:6]
train_split <- 1:6
test_split <- 7:8
bref <- c(0, 0)
res <- expRLDT(obsdata, train_split, test_split, bref, truemodel, 1)
print(paste("The estimated model of R-LDT is: beta=(", round(res[[4]][1],digits=2), ",",  round(res[[4]][2],digits=2), ")"))
sprintf("The test accuracy of R-LDT is: %.2f%%", res[[2]] * 100)
sprintf("The parameter mse of R-LDT is: %.2f", res[[3]])

## Toy 2
data2<-read.xlsx("synthetic/toy2.xlsx", sheetIndex=1, header=TRUE)
obsdata<-data2[,1:6]
train_split <- 1:7
test_split <- 8:9
bref <- c(0, 0)
res <- expRLDT(obsdata, train_split, test_split, bref, truemodel, 1)
print(paste("The estimated model R-LDT is: beta=(", round(res[[4]][1],digits=2), ",",  round(res[[4]][2],digits=2), ")"))
sprintf("The test accuracy of R-LDT is: %.2f%%", res[[2]] * 100)
sprintf("The parameter mse of R-LDT is: %.2f", res[[3]])

## Toy 3 
data3<-read.xlsx("synthetic/toy3.xlsx", sheetIndex=1, header=TRUE)
obsdata<-data3[,1:6]
train_split <- 1:7
test_split <- 8:9
bref <- c(0, 0)
res <- expRLDT(obsdata, train_split, test_split, bref, truemodel, 1)
print(paste("The estimated model of R-LDT is: beta=(", round(res[[4]][1],digits=2), ",",  round(res[[4]][2],digits=2), ")"))
sprintf("The test accuracy of R-LDT is: %.2f%%", res[[2]] * 100)
sprintf("The parameter mse of R-LDT is: %.2f", res[[3]])

