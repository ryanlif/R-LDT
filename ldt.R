LDT <- function(X, Pert, R, Y, lambda, bref){
  n1 <- ncol(X)
  n2 <- nrow(X)
  bi <- c()
  ai <- c()
  for (j in 1:n2){
    bi<-c(bi, 1-Y[j]*(R[j]-Pert[j]))
    for (i in 1:n1){
      ai <- c(ai, -Y[j]*X[j,i])
    }
    for (i in 1:n2){
      if (i == j){
        ai <- c(ai, 1)
      }
      else {
        ai <- c(ai, 0)
      }
    }
  }
  model <- list()
  model$modelsense <- 'min'
  model$A <- matrix(ai, nrow=n2, ncol=n1+n2, byrow=T)
  model$Q <- as.matrix(bdiag(0.5*diag(n1), matrix(0,n2,n2)))
  model$obj <- c(-bref, rep(lambda, n2))
  model$rhs<-bi
  model$sense <- rep('>=',n2)
  name <- c()
  for (i in 1:n1){
    name <- c(name, paste('a',i,sep=""))
  }
  for (i in 1:n2){
    name <- c(name, paste('e',i,sep=""))
  }
  model$varnames <- name
  model$lb <- c(rep(-Inf, n1), rep(0, n2))
  model$ub <- rep(Inf, n1+n2)
  res <- gurobi(model)
  para <- res$x
  return(para[1:n1])
}

expLDT <- function(obdata, train_split, test_split, bref, truemodel){
  train_acc <- c()
  test_acc <- c()
  mse_acc <- c()
  data <- obdata
  d <- ncol(data) - 4
  nu <- nrow(data)
  
  X <- data[, 1:d]
  ones <- rep(1,nu)
  R <- data$r
  Y <- data$y
  Yt <- data$ytrue
  St<-train_split
  Se<-test_split
  ntrain <- length(St) 
  ntest <- length(Se)
  beta <- rep(0, d)
  pert <- rep(0, nu)
  lambda1 <- 1
  beta <- LDT(X[St,], rep(0,length(St)), R[St], Y[St], lambda1, rep(0,d))
  mse <- mean((beta-truemodel)^2)
  ytr <- 2*(R[St]>as.matrix(X[St,]) %*% beta)-1
  train_acc <- sum(Yt[St]==ytr)/ntrain
  yte <- 2*(R[Se]>as.matrix(X[Se,]) %*% beta)-1
  test_acc <- sum(Yt[Se]==yte)/ntest
   expLDT <- list(train_acc, test_acc, mse, beta)
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
res <- expLDT(obsdata, train_split, test_split, bref, truemodel)
print(paste("The estimated model of LDT is: beta=(", round(res[[4]][1],digits=2), ",",  round(res[[4]][2],digits=2), ")"))
sprintf("The test accuracy of LDT is: %.2f%%", res[[2]] * 100)
sprintf("The parameter mse of LDT is: %.2f", res[[3]])

## Toy 2
data2<-read.xlsx("synthetic/toy2.xlsx", sheetIndex=1, header=TRUE)
obsdata<-data2[,1:6]
train_split <- 1:7
test_split <- 8:9
bref <- c(0, 0)
res <- expLDT(obsdata, train_split, test_split, bref, truemodel)
print(paste("The estimated model of LDT is: beta=(", round(res[[4]][1],digits=2), ",",  round(res[[4]][2],digits=2), ")"))
sprintf("The test accuracy of LDT is: %.2f%%", res[[2]] * 100)
sprintf("The parameter mse of LDT is: %.2f", res[[3]])

## Toy 3 
data3<-read.xlsx("synthetic/toy3.xlsx", sheetIndex=1, header=TRUE)
obsdata<-data3[,1:6]
train_split <- 1:7
test_split <- 8:9
bref <- c(0, 0)
res <- expLDT(obsdata, train_split, test_split, bref, truemodel)
print(paste("The estimated model of LDT is: beta=(", round(res[[4]][1],digits=2), ",",  round(res[[4]][2],digits=2), ")"))
sprintf("The test accuracy of LDT is: %.2f%%", res[[2]] * 100)
sprintf("The parameter mse of LDT is: %.2f", res[[3]])