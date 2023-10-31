require(MASS)
library(Matrix)

simdata <- function(d, nu, rate){
  ################ Non-corrupted data generation ################
  # d: Attribute dimension
  # nu:  Training sample size
  # rate: Training rate
  ###############################################################
  train_split <- c()
  test_split <- c()
  fulldata <- c()
  ### Generate the true coefficients
  coef <- mvrnorm(n=1, rep(0,d), diag(1,d))
  normcoef <- sqrt(sum(coef^2))
  coef <- coef / normcoef
  truemodel <- coef
  ### Generate the observations
  obs <- c()
  nuser <- round(nu/rate) # Total sample size
  sigma <- 1
  for (j in 1:d){
    col <- rnorm(nuser, 0, sigma)
    obs <- cbind(obs, col)
  }
  ### Generate the true model
  D <- obs %*% coef #Latent threshold generation
  ### Contribution reward
  #r <- D + runif(n=nuser, -1, 1)
  ## Random reward 
  r <- rnorm(n=nuser, 0, 2)
  ### Predictive reward
  #b0 <- rnorm(d, mean=0, sd=1)
  #r <- obs %*% b0 + runif(n=nuser, -1, 1)
  y <- c()
  for (k in 1:nuser){
    if (r[k]>=D[k]){
      y <- c(y, 1)
    }
    else {
      y <- c(y,-1)
    }
  }
  ### Data split
  ntrain <- nu
  ntest <- nuser - nu
  St<- 1:nu
  train_split <- St
  Se<- (nu+1):nuser
  test_split <- Se  
  
  userdata <- cbind(obs, r, D, y)
  colnames(userdata) <- c(paste("col",1:d,sep=""), "r","D", "y")
  simdata <- list(userdata, truemodel, train_split, test_split)
}

corruptdata <- function(simudata, d, ratio, magnitude, len){
  userdata <- data.frame(simudata[[1]])
  truemodel <- simudata[[2]]
  St <- simudata[[3]]
  Se <- simudata[[4]]
  traindat <- userdata[St, ]
  testdat <- userdata[Se, ]
  ni <- nrow(userdata)
  ntrain <- nrow(traindat)
  ntest <- nrow(testdat)
  userdata <- cbind(userdata, userdata$y)
  nameset <-c()
  for (t in 1:d){
    nameset <-c(nameset, paste("obs",t, sep=""))
  }
  nameset <- c(nameset, "r","D","y","ytrue")
  colnames(userdata) <- nameset
  pos <- which(userdata[St,]$ytrue==1)
  no <- round(length(pos) * ratio)
  if (no==0){
    Stc <- c()
    ds <- list(userdata, truemodel,St, Se, Stc, 0)
    return(ds)
  }
  Stc <- sample(St[pos], no) 
  print("Hi")
  print(Stc)
  D <- userdata$D
  r <- userdata$r
  y <- userdata$y
  b <- runif(ni, magnitude, len * magnitude)
  for (i in 1:no){
    k <- Stc[i]
    r[k] <- r[k] + b[k]
    # userdata[k,1:d] <- userdata[k, 1:d] + runif(d, 0, 2)
    D[k] = sum(userdata[k, 1:d] * truemodel)
    if (r[k] >= D[k]){
      y[k] <- -1
    }else{
      y[k] <- 1
    }
  }
  userdata$r <- r
  userdata$y <- y
  ds <- list(userdata, truemodel,St, Se, Stc, no / ntrain)
  return(ds)
}
