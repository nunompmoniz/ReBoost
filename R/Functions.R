#'
#' AdaBoost.R2
#'

AdaBoost.R2 <- function(form,train,test,t_final=100,power=2,...) {

  require(spatstat)

  models <- list()
  betas <- c()
  pred.mat <- c()

  ind.y <- which(colnames(train)==as.character(form[[2]]))

  n <- nrow(train) #size of train

  weights <- rep(1/n,n) #initialize weights

  err_t <- 0

  for (t in 1:t_final) {

    train.ind <- sample(1:n,n,replace=TRUE,prob=weights)
    m <- rpart(form,train[train.ind,],...)

    models[[t]] <- m

    f <- predict(m,train)

    ar <- abs(f-train[,ind.y])
    ar <- (ar/max(ar))^power

    err_t <- sum(weights*ar)

    if(err_t>=0.5) break

    beta_t <- err_t / (1-err_t)
    betas[[t]] <- beta_t

    weights <- weights*(beta_t^(1-err_t))
    weights <- weights/sum(weights)

  }

  if(t==t_final) t <- t+1

  for(i in 1:(t-1)) {

    pred.mat <- cbind(pred.mat,predict(models[[i]],test))

  }

  finalpreds <- c()
  for(i in 1:nrow(pred.mat)) {
    finalpreds <- c(finalpreds,weighted.median(pred.mat[i,],unlist(betas)))
  }

  names(finalpreds) <- rownames(test)

  finalpreds

}

#'
#' SMOTEd AdaBoost.R2
#'

SMOTEBoost.R2 <- function(form,train,test,t_final=100,power=2,perc.O=1.5,rel.thr=0.9,k=3,coef=1.5,...) {

  require(spatstat)

  models <- list()
  betas <- c()
  pred.mat <- c()

  ind.y <- which(colnames(train)==as.character(form[[2]]))

  n <- nrow(train) #size of train

  weights <- rep(1/n,n) #initialize weights

  err_t <- 0

  pc <- UBL::phi.control(y = train[,ind.y],method = "extremes",coef=coef)

  for (t in 1:t_final) {

    train.ind <- sample(1:n,n,replace=TRUE,prob=weights)
    new.train <- train[train.ind,]

    new.train <- adaSMOTE(form,new.train,perc.O,rel.thr,k,pc)

    m <- rpart(form,new.train,...)

    models[[t]] <- m

    f <- predict(m,train)

    ar <- abs(f-train[,ind.y])
    ar <- (ar/max(ar))^power

    err_t <- sum(weights*ar)

    if(err_t>=0.5) break

    beta_t <- err_t / (1-err_t)
    betas[[t]] <- beta_t

    weights <- weights*(beta_t^(1-err_t))
    weights <- weights/sum(weights)

  }

  if(t==t_final) t <- t+1

  for(i in 1:(t-1)) {

    pred.mat <- cbind(pred.mat,predict(models[[i]],test))

  }

  finalpreds <- c()
  for(i in 1:nrow(pred.mat)) {
    finalpreds <- c(finalpreds,weighted.median(pred.mat[i,],unlist(betas)))
  }

  names(finalpreds) <- rownames(test)

  finalpreds

}

#'
#' AdaBoost.RQ
#'

AdaBoost.RQ <- function(form,train,test,t_final=100,power=2,...) {

  models <- list()
  betas <- c()
  pred.mat <- c()

  ind.y <- which(colnames(train)==as.character(form[[2]]))

  n <- nrow(train) #size of train

  weights <- rep(1/n,n) #initialize weights

  err_t <- 0

  for (t in 1:t_final) {

    train.ind <- sample(1:n,n,replace=TRUE,prob=weights)
    m <- rpart(form,train[train.ind,],...)

    models[[t]] <- m

    f <- predict(m,train)

    abs.err <- abs(as.numeric((f-train[,ind.y])))

    large.err.ind <- which(abs.err>boxplot.stats(abs.err)$stats[3])

    err_t <- sum(weights[large.err.ind])

    beta_t <- err_t^power
    betas[[t]] <- beta_t

    weights[] <- beta_t; weights[large.err.ind] <- 1;
    weights <- weights/sum(weights)

    # weights <- softmax(weights)

  }

  num <- 0
  for(t in 1:t_final) {

    preds <- predict(models[[t]],test)
    pred.mat <- cbind(pred.mat,preds)
    num <- num + (log(1/betas[t]) * preds)


  }

  finalpreds <- num/sum(log(1/betas))

  names(finalpreds) <- rownames(test)

  finalpreds

}

#'
#' SMOTEd AdaBoost.RQ
#'

SMOTEBoost.RQ <- function(form,train,test,t_final=100,power=2,perc.O=1.5,rel.thr=0.9,k=3,coef=1.5,...) {

  models <- list()
  betas <- c()
  pred.mat <- c()

  ind.y <- which(colnames(train)==as.character(form[[2]]))

  n <- nrow(train) #size of train

  weights <- rep(1/n,n) #initialize weights

  err_t <- 0

  pc <- UBL::phi.control(y = train[,ind.y],method = "extremes",coef=coef)

  for (t in 1:t_final) {

    train.ind <- sample(1:n,n,replace=TRUE,prob=weights)
    new.train <- train[train.ind,]

    new.train <- adaSMOTE(form,new.train,perc.O,rel.thr,k=k,pc)

    m <- rpart(form,new.train,...)

    models[[t]] <- m

    f <- predict(m,train)

    abs.err <- abs(as.numeric((f-train[,ind.y])))

    large.err.ind <- which(abs.err>boxplot.stats(abs.err)$stats[3])

    err_t <- sum(weights[large.err.ind])
    beta_t <- err_t^power
    betas[[t]] <- beta_t

    weights[] <- beta_t; weights[large.err.ind] <- 1;
    weights <- weights/sum(weights)

  }

  num <- 0
  for(t in 1:t_final) {

    preds <- predict(models[[t]],test)
    pred.mat <- cbind(pred.mat,preds)
    num <- num + (log(1/betas[t]) * preds)


  }

  finalpreds <- num/sum(log(1/betas))

  names(finalpreds) <- rownames(test)

  finalpreds

}

#'
#' AdaBoost.RT
#'

AdaBoost.RT <- function(form,train,test,t_final=100,thr=0.1,power=2,...) {

  models <- list()
  betas <- c()
  pred.mat <- c()

  ind.y <- which(colnames(train)==as.character(form[[2]]))

  n <- nrow(train) #size of train

  weights <- rep(1/n,n) #initialize weights

  err_t <- 0

  for (t in 1:t_final) {

    train.ind <- sample(1:n,n,replace=TRUE,prob=weights)
    m <- rpart(form,train[train.ind,],...)

    models[[t]] <- m

    f <- predict(m,train)

    are <- abs((f-train[,ind.y])/train[,ind.y]); are[is.na(are)] <- 0

    err_t <- sum(weights[are>thr])

    beta_t <- err_t^power
    betas[[t]] <- beta_t

    weights[are <= thr] <- beta_t; weights[are > thr] <- 1
    weights <- weights/sum(weights)

  }

  num <- 0
  for(t in 1:t_final) {

    preds <- predict(models[[t]],test)
    pred.mat <- cbind(pred.mat,preds)
    num <- num + (log(1/betas[t]) * preds)

  }

  finalpreds <- num/sum(log(1/betas))

  names(finalpreds) <- rownames(test)

  finalpreds

}

#'
#' SMOTEd AdaBoost.RT
#'

SMOTEBoost.RT <- function(form,train,test,t_final=100,thr=0.1,power=2,perc.O=1.5,rel.thr=0.9,k=3,coef=1.5,...) {

  nrow(train)
  nrow(test)
  models <- list()
  betas <- c()
  pred.mat <- c()

  ind.y <- which(colnames(train)==as.character(form[[2]]))

  n <- nrow(train) #size of train

  weights <- rep(1/n,n) #initialize weights

  err_t <- 0

  pc <- UBL::phi.control(y = train[,ind.y],method = "extremes",coef=coef)

  for (t in 1:t_final) {

    train.ind <- sample(1:n,n,replace=TRUE,prob=weights)
    new.train <- train[train.ind,]

    new.train <- adaSMOTE(form,new.train,perc.O,rel.thr=rel.thr,k=k,pc=pc)

    m <- rpart(form,new.train,...)

    models[[t]] <- m

    f <- predict(m,train)

    are <- abs(as.numeric((f-train[,ind.y])/train[,ind.y])); are[is.na(are)] <- 0
    err_t <- sum(weights[are>thr])

    beta_t <- err_t^power
    betas[[t]] <- beta_t

    weights[are <= thr] <- beta_t; weights[are > thr] <- 1
    weights <- weights/sum(weights)

  }

  num <- 0
  for(t in 1:t_final) {

    preds <- predict(models[[t]],test)
    pred.mat <- cbind(pred.mat,preds)
    num <- num + (log(1/betas[t]) * preds)


  }

  finalpreds <- num/sum(log(1/betas))

  names(finalpreds) <- rownames(test)

  finalpreds

}

#'
#' AdaBoost.RT+
#'

AdaBoost.RTPlus <- function(form,train,test,t_final=100,thr=0.01,power=2,sigma=0.5,...) {

  require(MASS)

  models <- list()
  betas <- c()
  pred.mat <- c()
  train_pred.mat <- c()

  ind.y <- which(colnames(train)==as.character(form[[2]]))

  n <- nrow(train) #size of train

  weights <- rep(1/n,n) #initialize weights

  err_t <- 0

  for (t in 1:t_final) {

    train.ind <- sample(1:n,n,replace=TRUE,prob=weights)
    m <- rpart(form,train[train.ind,],...)

    models[[t]] <- m

    f <- predict(m,train)

    are <- abs((f-train[,ind.y])/train[,ind.y]); are[is.na(are)] <- 0

    err_t <- sum(weights[are>thr])

    beta_t <- err_t^power
    betas[[t]] <- beta_t

    weights[are <= thr] <- beta_t; weights[are > thr] <- 1
    weights <- weights/sum(weights)

  }

  for(t in 1:t_final) {
    preds <- predict(models[[t]],train)
    train_pred.mat <- cbind(train_pred.mat,preds)
    preds <- predict(models[[t]],test)
    pred.mat <- cbind(pred.mat,preds)
  }

  delta <- 0


  if(is.null(sigma)) { # no regularization
    delta <- t(ginv(t(train_pred.mat))) %*% train[,ind.y]
  } else {
    delta <- t(ginv(train_pred.mat %*% t(train_pred.mat) + sigma * diag(nrow(train))) %*% train_pred.mat) %*% train[,ind.y]
  }

  finalpreds <- pred.mat %*% delta

  names(finalpreds) <- rownames(test)

  finalpreds

}

#'
#' SMOTEd AdaBoost.RT+
#'

SMOTEBoost.RTPlus <- function(form,train,test,t_final=100,thr=0.01,power=2,sigma=0.5,perc.O=1.5,rel.thr=0.9,k=3,coef=1.5,...) {

  require(MASS)

  models <- list()
  betas <- c()
  pred.mat <- c()
  train_pred.mat <- c()

  ind.y <- which(colnames(train)==as.character(form[[2]]))

  n <- nrow(train) #size of train

  weights <- rep(1/n,n) #initialize weights

  err_t <- 0

  pc <- UBL::phi.control(y = train[,ind.y],method = "extremes",coef=coef)

  for (t in 1:t_final) {

    train.ind <- sample(1:n,n,replace=TRUE,prob=weights)
    new.train <- train[train.ind,]

    new.train <- adaSMOTE(form,new.train,perc.O,rel.thr,k=k,pc)

    m <- rpart(form,new.train,...)

    models[[t]] <- m

    f <- predict(m,train)

    are <- abs((f-train[,ind.y])/train[,ind.y]); are[is.na(are)] <- 0

    err_t <- sum(weights[are>thr])

    beta_t <- err_t^power
    betas[[t]] <- beta_t

    weights[are <= thr] <- beta_t; weights[are > thr] <- 1
    weights <- weights/sum(weights)

  }

  for(t in 1:t_final) {
    preds <- predict(models[[t]],train)
    train_pred.mat <- cbind(train_pred.mat,preds)
    preds <- predict(models[[t]],test)
    pred.mat <- cbind(pred.mat,preds)
  }

  delta <- 0


  if(is.null(sigma)) { # no regularization
    delta <- t(ginv(t(train_pred.mat))) %*% train[,ind.y]
  } else {
    delta <- t(ginv(train_pred.mat %*% t(train_pred.mat) + sigma * diag(nrow(train))) %*% train_pred.mat) %*% train[,ind.y]
  }

  finalpreds <- pred.mat %*% delta

  names(finalpreds) <- rownames(test)

  finalpreds

}

#'
#' BEMBoost
#'

BEMBoost <- function(form,train,test,t_final=100,BEM=0.5,...) {

  models <- list()
  betas <- c()
  pred.mat <- c()

  ind.y <- which(colnames(train)==as.character(form[[2]]))

  n <- nrow(train) #size of train

  weights <- rep(1/n,n) #initialize weights

  err_t <- 0

  for (t in 1:t_final) {

    train.ind <- sample(1:n,n,replace=TRUE,prob=weights)
    m <- rpart(form,train[train.ind,],...)

    models[[t]] <- m

    f <- predict(m,train)

    ae <- abs(f-train[,ind.y])

    grtBEM <- ae>BEM
    errCnt <- sum(grtBEM)

    #' Snippet
    #' For situations where the ad-hoc big error margin is too large,
    #' This snippet reduces it by half until the BEM error count is more than zero.
    if(t==1 & errCnt==0) {

      while(errCnt==0) {
        BEM <- BEM/2
        grtBEM <- ae>BEM
        errCnt <- sum(grtBEM)
      }

    }

    if(errCnt==0) break

    upfactor <- n/errCnt
    downfactor <- 1/upfactor

    lwrBEM <- !grtBEM

    weights[grtBEM] <- weights[grtBEM] * upfactor
    weights[lwrBEM] <- weights[lwrBEM] * downfactor

    weights <- weights/sum(weights)

  }

  if(t==t_final) t <- t+1

  for(i in 1:(t-1)) {

    preds <- predict(models[[i]],test)
    pred.mat <- cbind(pred.mat,preds)

  }

  finalpreds <- rowMeans(pred.mat)

  names(finalpreds) <- rownames(test)

  finalpreds

}

#'
#' SMOTEd BEMBoost
#'

SMOTEBoost.BEM <- function(form,train,test,t_final=100,BEM=0.5,perc.O=1.5,rel.thr=0.9,k=3,coef=1.5,...) {

  models <- list()
  betas <- c()
  pred.mat <- c()

  ind.y <- which(colnames(train)==as.character(form[[2]]))

  n <- nrow(train) #size of train

  weights <- rep(1/n,n) #initialize weights

  err_t <- 0

  pc <- UBL::phi.control(y = train[,ind.y],method = "extremes",coef=coef)

  for (t in 1:t_final) {

    train.ind <- sample(1:n,n,replace=TRUE,prob=weights)
    new.train <- train[train.ind,]

    new.train <- adaSMOTE(form,new.train,perc.O,rel.thr,k=k,pc)

    m <- rpart(form,new.train,...)

    models[[t]] <- m

    f <- predict(m,train)

    ae <- abs(f-train[,ind.y])

    grtBEM <- ae>BEM
    errCnt <- sum(grtBEM)

    #' Snippet
    #' For situations where the ad-hoc big error margin is too large,
    #' This snippet reduces it by half until the BEM error count is more than zero.
    if(t==1 & errCnt==0) {

      while(errCnt==0) {
        BEM <- BEM/2
        grtBEM <- ae>BEM
        errCnt <- sum(grtBEM)
      }

    }

    if(errCnt==0) break

    upfactor <- n/errCnt
    downfactor <- 1/upfactor

    lwrBEM <- !grtBEM

    weights[grtBEM] <- weights[grtBEM] * upfactor
    weights[lwrBEM] <- weights[lwrBEM] * downfactor

    weights <- weights/sum(weights)

  }

  if(t==t_final) t <- t+1

  for(i in 1:(t-1)) {

    preds <- predict(models[[i]],test)
    pred.mat <- cbind(pred.mat,preds)

  }

  finalpreds <- rowMeans(pred.mat)

  names(finalpreds) <- rownames(test)

  finalpreds

}
