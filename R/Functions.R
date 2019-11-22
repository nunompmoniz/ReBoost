#' AdaBoost.R2
#'
#' @references
#'
#' @param form The model formula.
#' @param train A data.frame with the training data.
#' @param test A data.frame with the test data.
#' @param t_final The number of maximum boosting iterations. Default is 100.
#' @param power Type of loss function, e.g. linear (1), squared (2). Default is 2.
#' @param ... Dots are passed to rpart
#'
#' @return Returns a vector with the predictions made by AdaBoost.R2.
#' @export
#'
#' @importFrom rpart rpart
#' @importFrom spatstat weighted.median
#'
#' @examples
#' data(Boston,package="MASS")
#'
#' idx <- sample(1:nrow(Boston),nrow(Boston)*0.75)
#' form <- medv ~ .
#'
#' train <- Boston[idx,]
#' test <- Boston[-idx,]
#'
#' preds <- AdaBoost.R2(form,train,test)
#'
AdaBoost.R2 <- function(form,train,test,t_final=100,power=2,...) {

  require(rpart)
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

#' SMOTEBoost variant of AdaBoost.R2
#'
#' @param form The model formula.
#' @param train A data.frame with the training data.
#' @param test A data.frame with the test data.
#' @param t_final The number of maximum boosting iterations. Default is 100.
#' @param power Type of loss function, e.g. linear (1), squared (2). Default is 2.
#' @param perc.O Percentage for Oversampling via SMOTE, i.e. percentage of extreme cases to be generated. Default is 1.5.
#' @param rel.thr Relevance threshold. Default is 0.9.
#' @param k Number of neighbours used in SMOTE. Defaults to 3.
#' @param coef Coefficient used in boxplot statistics, which is used to create the relevance function. Default is 1.5.
#' @param ... Dots are passed to rpart
#'
#' @return Returns a vector with the predictions made by SMOTEBoost.R2.
#' @export
#'
#' @importFrom rpart rpart
#' @importFrom spatstat weighted.median
#' @importFrom IRon phi.control
#'
#' @examples
#' data(Boston,package="MASS")
#'
#' idx <- sample(1:nrow(Boston),nrow(Boston)*0.75)
#' form <- medv ~ .
#'
#' train <- Boston[idx,]
#' test <- Boston[-idx,]
#'
#' preds <- SMOTEBoost.R2(form,train,test)
#'
SMOTEBoost.R2 <- function(form,train,test,t_final=100,power=2,perc.O=1.5,rel.thr=0.9,k=3,coef=1.5,...) {

  require(rpart)
  require(spatstat)
  require(IRon)

  models <- list()
  betas <- c()
  pred.mat <- c()

  ind.y <- which(colnames(train)==as.character(form[[2]]))

  n <- nrow(train) #size of train

  weights <- rep(1/n,n) #initialize weights

  err_t <- 0

  pc <- IRon::phi.control(y = train[,ind.y],method = "extremes",coef=coef)

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

#' RUSBoost variant of AdaBoost.R2
#'
#' @param form The model formula.
#' @param train A data.frame with the training data.
#' @param test A data.frame with the test data.
#' @param t_final The number of maximum boosting iterations. Default is 100.
#' @param power Type of loss function, e.g. linear (1), squared (2). Default is 2.
#' @param perc.U Percentage for Undersampling via Random Undersampling, i.e. percentage of cases with normal values to remain in the new dataset. Default is 0.9.
#' @param rel.thr Relevance threshold. Default is 0.9.
#' @param coef Coefficient used in boxplot statistics, which is used to create the relevance function. Default is 1.5.
#' @param ... Dots are passed to rpart
#'
#' @return Returns a vector with the predictions made by RUSBoost.R2.
#' @export
#'
#' @importFrom UBL RandUnderRegress
#' @importFrom IRon phi.control
#' @importFrom rpart rpart
#' @importFrom spatstat weighted.median
#'
#' @examples
#' data(Boston,package="MASS")
#'
#' idx <- sample(1:nrow(Boston),nrow(Boston)*0.75)
#' form <- medv ~ .
#'
#' train <- Boston[idx,]
#' test <- Boston[-idx,]
#'
#' preds <- RUSBoost.R2(form,train,test)
#'
RUSBoost.R2 <- function(form,train,test,t_final=100,power=2,perc.U=0.9,rel.thr=0.9,coef=1.5,...) {

  require(spatstat)
  require(rpart)
  require(UBL)
  require(IRon)

  models <- list()
  betas <- c()
  pred.mat <- c()

  ind.y <- which(colnames(train)==as.character(form[[2]]))

  n <- nrow(train) #size of train

  weights <- rep(1/n,n) #initialize weights

  err_t <- 0

  pc <- IRon::phi.control(y = train[,ind.y],method = "extremes",coef=coef)

  for (t in 1:t_final) {

    train.ind <- sample(1:n,n,replace=TRUE,prob=weights)
    new.train <- train[train.ind,]

    new.train <- UBL::RandUnderRegress(form,new.train,pc,rel.thr,perc.U)

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

#' ROSBoost variant of AdaBoost.R2
#'
#' @param form The model formula.
#' @param train A data.frame with the training data.
#' @param test A data.frame with the test data.
#' @param t_final The number of maximum boosting iterations. Default is 100.
#' @param power Type of loss function, e.g. linear (1), squared (2). Default is 2.
#' @param perc.O Percentage for Oversampling via Random Oversampling, i.e. percentage of extreme cases to be generated. Default is 1.5.
#' @param rel.thr Relevance threshold. Default is 0.9.
#' @param coef Coefficient used in boxplot statistics, which is used to create the relevance function. Default is 1.5.
#' @param ... Dots are passed to rpart
#'
#' @return Returns a vector with the predictions made by ROSBoost.R2.
#' @export
#'
#' @importFrom UBL RandOverRegress
#' @importFrom IRon phi.control
#' @importFrom rpart rpart
#' @importFrom spatstat weighted.median
#'
#' @examples
#' data(Boston,package="MASS")
#'
#' idx <- sample(1:nrow(Boston),nrow(Boston)*0.75)
#' form <- medv ~ .
#'
#' train <- Boston[idx,]
#' test <- Boston[-idx,]
#'
#' preds <- ROSBoost.R2(form,train,test)
#'
ROSBoost.R2 <- function(form,train,test,t_final=100,power=2,perc.O=1.5,rel.thr=0.9,coef=1.5,...) {

  require(spatstat)
  require(rpart)
  require(UBL)
  require(IRon)

  models <- list()
  betas <- c()
  pred.mat <- c()

  ind.y <- which(colnames(train)==as.character(form[[2]]))

  n <- nrow(train) #size of train

  weights <- rep(1/n,n) #initialize weights

  err_t <- 0

  pc <- IRon::phi.control(y = train[,ind.y],method = "extremes",coef=coef)

  for (t in 1:t_final) {

    train.ind <- sample(1:n,n,replace=TRUE,prob=weights)
    new.train <- train[train.ind,]

    new.train <- UBL::RandOverRegress(form,new.train,pc,rel.thr,perc.O)

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

#' AdaBoost.RQ
#'
#' @references
#'
#' @param form The model formula.
#' @param train A data.frame with the training data.
#' @param test A data.frame with the test data.
#' @param t_final The number of maximum boosting iterations. Default is 100.
#' @param power Type of loss function, e.g. linear (1), squared (2). Default is 2.
#' @param ... Dots are passed to rpart
#'
#' @return Returns a vector with the predictions made by AdaBoost.RQ.
#' @export
#'
#' @importFrom rpart rpart
#'
#' @examples
#' data(Boston,package="MASS")
#'
#' idx <- sample(1:nrow(Boston),nrow(Boston)*0.75)
#' form <- medv ~ .
#'
#' train <- Boston[idx,]
#' test <- Boston[-idx,]
#'
#' preds <- AdaBoost.RQ(form,train,test)
AdaBoost.RQ <- function(form,train,test,t_final=100,power=2,...) {

  require(rpart)

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

#' SMOTEBoost variant of AdaBoost.RQ
#'
#' @param form The model formula.
#' @param train A data.frame with the training data.
#' @param test A data.frame with the test data.
#' @param t_final The number of maximum boosting iterations. Default is 100.
#' @param power Type of loss function, e.g. linear (1), squared (2). Default is 2.
#' @param perc.O Percentage for Oversampling via SMOTE, i.e. percentage of extreme cases to be generated. Default is 1.5.
#' @param rel.thr Relevance threshold. Default is 0.9.
#' @param k Number of neighbours used in SMOTE. Defaults to 3.
#' @param coef Coefficient used in boxplot statistics, which is used to create the relevance function. Default is 1.5.
#' @param ... Dots are passed to rpart
#'
#' @return Returns a vector with the predictions made by SMOTEBoost.RQ.
#' @export
#'
#' @importFrom IRon phi.control
#' @importFrom rpart rpart
#' @importFrom grDevices boxplot.stats
#'
#' @examples
#' data(Boston,package="MASS")
#'
#' idx <- sample(1:nrow(Boston),nrow(Boston)*0.75)
#' form <- medv ~ .
#'
#' train <- Boston[idx,]
#' test <- Boston[-idx,]
#'
#' preds <- SMOTEBoost.RQ(form,train,test)
#'
SMOTEBoost.RQ <- function(form,train,test,t_final=100,power=2,perc.O=1.5,rel.thr=0.9,k=3,coef=1.5,...) {

  require(IRon)
  require(rpart)

  models <- list()
  betas <- c()
  pred.mat <- c()

  ind.y <- which(colnames(train)==as.character(form[[2]]))

  n <- nrow(train) #size of train

  weights <- rep(1/n,n) #initialize weights

  err_t <- 0

  pc <- IRon::phi.control(y = train[,ind.y],method = "extremes",coef=coef)

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

#' RUSBoost variant of AdaBoost.RQ
#'
#' @param form The model formula.
#' @param train A data.frame with the training data.
#' @param test A data.frame with the test data.
#' @param t_final The number of maximum boosting iterations. Default is 100.
#' @param power Type of loss function, e.g. linear (1), squared (2). Default is 2.
#' @param perc.U Percentage for Undersampling via Random Undersampling, i.e. percentage of cases with normal values to remain in the new dataset. Default is 0.9.
#' @param rel.thr Relevance threshold. Default is 0.9.
#' @param coef Coefficient used in boxplot statistics, which is used to create the relevance function. Default is 1.5.
#' @param ... Dots are passed to rpart
#'
#' @return Returns a vector with the predictions made by RUSBoost.RQ.
#' @export
#'
#' @importFrom UBL RandUnderRegress
#' @importFrom IRon phi.control
#' @importFrom rpart rpart
#' @importFrom grDevices boxplot.stats
#'
#' @examples
#' data(Boston,package="MASS")
#'
#' idx <- sample(1:nrow(Boston),nrow(Boston)*0.75)
#' form <- medv ~ .
#'
#' train <- Boston[idx,]
#' test <- Boston[-idx,]
#'
#' preds <- RUSBoost.RQ(form,train,test)
#'
RUSBoost.RQ <- function(form,train,test,t_final=100,power=2,perc.U=0.9,rel.thr=0.9,coef=1.5,...) {

  require(UBL)
  require(rpart)
  require(IRon)

  models <- list()
  betas <- c()
  pred.mat <- c()

  ind.y <- which(colnames(train)==as.character(form[[2]]))

  n <- nrow(train) #size of train

  weights <- rep(1/n,n) #initialize weights

  err_t <- 0

  pc <- IRon::phi.control(y = train[,ind.y],method = "extremes",coef=coef)

  for (t in 1:t_final) {

    train.ind <- sample(1:n,n,replace=TRUE,prob=weights)
    new.train <- train[train.ind,]

    new.train <- UBL::RandUnderRegress(form,new.train,pc,rel.thr,perc.U)

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

#' ROSBoost variant of AdaBoost.RQ
#'
#' @param form The model formula.
#' @param train A data.frame with the training data.
#' @param test A data.frame with the test data.
#' @param t_final The number of maximum boosting iterations. Default is 100.
#' @param power Type of loss function, e.g. linear (1), squared (2). Default is 2.
#' @param perc.O Percentage for Oversampling via Random Oversampling, i.e. percentage of extreme cases to be generated. Default is 1.5.
#' @param rel.thr Relevance threshold. Default is 0.9.
#' @param coef Coefficient used in boxplot statistics, which is used to create the relevance function. Default is 1.5.
#' @param ... Dots are passed to rpart
#'
#' @return Returns a vector with the predictions made by ROSBoost.RQ.
#' @export
#'
#' @importFrom UBL RandOverRegress
#' @importFrom IRon phi.control
#' @importFrom rpart rpart
#' @importFrom grDevices boxplot.stats
#'
#' @examples
#' data(Boston,package="MASS")
#'
#' idx <- sample(1:nrow(Boston),nrow(Boston)*0.75)
#' form <- medv ~ .
#'
#' train <- Boston[idx,]
#' test <- Boston[-idx,]
#'
#' preds <- ROSBoost.RQ(form,train,test)
#'
ROSBoost.RQ <- function(form,train,test,t_final=100,power=2,perc.O=1.5,rel.thr=0.9,coef=1.5,...) {

  require(UBL)
  require(rpart)
  require(IRon)

  models <- list()
  betas <- c()
  pred.mat <- c()

  ind.y <- which(colnames(train)==as.character(form[[2]]))

  n <- nrow(train) #size of train

  weights <- rep(1/n,n) #initialize weights

  err_t <- 0

  pc <- IRon::phi.control(y = train[,ind.y],method = "extremes",coef=coef)

  for (t in 1:t_final) {

    train.ind <- sample(1:n,n,replace=TRUE,prob=weights)
    new.train <- train[train.ind,]

    new.train <- UBL::RandOverRegress(form,new.train,pc,rel.thr,perc.O)

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

#' AdaBoost.RT
#'
#' @references
#'
#' @param form The model formula.
#' @param train A data.frame with the training data.
#' @param test A data.frame with the test data.
#' @param t_final The number of maximum boosting iterations. Default is 100.
#' @param thr The error threshold. Default is 0.1.
#' @param power Type of loss function, e.g. linear (1), squared (2). Default is 2.
#' @param ... Dots are passed to rpart
#'
#' @return Returns a vector with the predictions made by AdaBoost.RT.
#' @export
#'
#' @importFrom rpart rpart
#'
#' @examples
#' data(Boston,package="MASS")
#'
#' idx <- sample(1:nrow(Boston),nrow(Boston)*0.75)
#' form <- medv ~ .
#'
#' train <- Boston[idx,]
#' test <- Boston[-idx,]
#'
#' preds <- AdaBoost.RT(form,train,test)
AdaBoost.RT <- function(form,train,test,t_final=100,thr=0.1,power=2,...) {

  require(rpart)

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

#' SMOTEBoost variant of AdaBoost.RT
#'
#' @param form The model formula.
#' @param train A data.frame with the training data.
#' @param test A data.frame with the test data.
#' @param thr The error threshold. Default is 0.1.
#' @param t_final The number of maximum boosting iterations. Default is 100.
#' @param power Type of loss function, e.g. linear (1), squared (2). Default is 2.
#' @param perc.O Percentage for Oversampling via SMOTE, i.e. percentage of extreme cases to be generated. Default is 1.5.
#' @param rel.thr Relevance threshold. Default is 0.9.
#' @param k Number of neighbours used in SMOTE. Defaults to 3.
#' @param coef Coefficient used in boxplot statistics, which is used to create the relevance function. Default is 1.5.
#' @param ... Dots are passed to rpart
#'
#' @return Returns a vector with the predictions made by SMOTEBoost.RT.
#' @export
#'
#' @importFrom IRon phi.control
#' @importFrom rpart rpart
#'
#' @examples
#' data(Boston,package="MASS")
#'
#' idx <- sample(1:nrow(Boston),nrow(Boston)*0.75)
#' form <- medv ~ .
#'
#' train <- Boston[idx,]
#' test <- Boston[-idx,]
#'
#' preds <- SMOTEBoost.RT(form,train,test)
#'
SMOTEBoost.RT <- function(form,train,test,t_final=100,thr=0.1,power=2,perc.O=1.5,rel.thr=0.9,k=3,coef=1.5,...) {

  require(IRon)
  require(rpart)

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

#' RUSBoost variant of AdaBoost.RT
#'
#' @param form The model formula.
#' @param train A data.frame with the training data.
#' @param test A data.frame with the test data.
#' @param thr The error threshold. Default is 0.1.
#' @param t_final The number of maximum boosting iterations. Default is 100.
#' @param power Type of loss function, e.g. linear (1), squared (2). Default is 2.
#' @param perc.U Percentage for Undersampling via Random Undersampling, i.e. percentage of cases with normal values to remain in the new dataset. Default is 0.9.
#' @param rel.thr Relevance threshold. Default is 0.9.
#' @param coef Coefficient used in boxplot statistics, which is used to create the relevance function. Default is 1.5.
#' @param ... Dots are passed to rpart
#'
#' @return Returns a vector with the predictions made by RUSBoost.RT.
#' @export
#'
#' @importFrom UBL RandUnderRegress
#' @importFrom IRon phi.control
#' @importFrom rpart rpart
#'
#' @examples
#' data(Boston,package="MASS")
#'
#' idx <- sample(1:nrow(Boston),nrow(Boston)*0.75)
#' form <- medv ~ .
#'
#' train <- Boston[idx,]
#' test <- Boston[-idx,]
#'
#' preds <- RUSBoost.RT(form,train,test)
#'
RUSBoost.RT <- function(form,train,test,t_final=100,thr=0.1,power=2,perc.U=0.9,rel.thr=0.9,coef=1.5,...) {

  require(UBL)
  require(IRon)
  require(rpart)

  models <- list()
  betas <- c()
  pred.mat <- c()

  ind.y <- which(colnames(train)==as.character(form[[2]]))

  n <- nrow(train) #size of train

  weights <- rep(1/n,n) #initialize weights

  err_t <- 0

  pc <- IRon::phi.control(y = train[,ind.y],method = "extremes",coef=coef)

  for (t in 1:t_final) {

    train.ind <- sample(1:n,n,replace=TRUE,prob=weights)
    new.train <- train[train.ind,]

    new.train <- UBL::RandUnderRegress(form,new.train,pc,rel.thr,perc.U)

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

#' ROSBoost variant of AdaBoost.RT
#'
#' @param form The model formula.
#' @param train A data.frame with the training data.
#' @param test A data.frame with the test data.
#' @param thr The error threshold. Default is 0.1.
#' @param t_final The number of maximum boosting iterations. Default is 100.
#' @param power Type of loss function, e.g. linear (1), squared (2). Default is 2.
#' @param perc.O Percentage for Oversampling via Random Oversampling, i.e. percentage of extreme cases to be generated. Default is 1.5.
#' @param rel.thr Relevance threshold. Default is 0.9.
#' @param coef Coefficient used in boxplot statistics, which is used to create the relevance function. Default is 1.5.
#' @param ... Dots are passed to rpart
#'
#' @return Returns a vector with the predictions made by ROSBoost.RT.
#' @export
#'
#' @importFrom UBL RandOverRegress
#' @importFrom IRon phi.control
#' @importFrom rpart rpart
#'
#' @examples
#' data(Boston,package="MASS")
#'
#' idx <- sample(1:nrow(Boston),nrow(Boston)*0.75)
#' form <- medv ~ .
#'
#' train <- Boston[idx,]
#' test <- Boston[-idx,]
#'
#' preds <- ROSBoost.RT(form,train,test)
#'
ROSBoost.RT <- function(form,train,test,t_final=100,thr=0.1,power=2,perc.O=1.5,rel.thr=0.9,coef=1.5,...) {

  require(UBL)
  require(rpart)
  require(IRon)

  models <- list()
  betas <- c()
  pred.mat <- c()

  ind.y <- which(colnames(train)==as.character(form[[2]]))

  n <- nrow(train) #size of train

  weights <- rep(1/n,n) #initialize weights

  err_t <- 0

  pc <- IRon::phi.control(y = train[,ind.y],method = "extremes",coef=coef)

  for (t in 1:t_final) {

    train.ind <- sample(1:n,n,replace=TRUE,prob=weights)
    new.train <- train[train.ind,]

    new.train <- UBL::RandOverRegress(form,new.train,pc,rel.thr,perc.O)

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

#' AdaBoost.RTPlus
#'
#' @references
#'
#' @param form The model formula.
#' @param train A data.frame with the training data.
#' @param test A data.frame with the test data.
#' @param t_final The number of maximum boosting iterations. Default is 100.
#' @param thr The error threshold. Default is 0.01.
#' @param power Type of loss function, e.g. linear (1), squared (2). Default is 2.
#' @param sigma Regularization factor. Default is 0.5.
#' @param ... Dots are passed to rpart
#'
#' @return Returns a vector with the predictions made by AdaBoost.RTPlus.
#' @export
#'
#' @importFrom rpart rpart
#' @importFrom MASS ginv
#'
#' @examples
#' data(Boston,package="MASS")
#'
#' idx <- sample(1:nrow(Boston),nrow(Boston)*0.75)
#' form <- medv ~ .
#'
#' train <- Boston[idx,]
#' test <- Boston[-idx,]
#'
#' preds <- AdaBoost.RTPlus(form,train,test)
AdaBoost.RTPlus <- function(form,train,test,t_final=100,thr=0.01,power=2,sigma=0.5,...) {

  require(MASS)
  require(rpart)

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

#' SMOTEBoost variant of AdaBoost.RTPlus
#'
#' @param form The model formula.
#' @param train A data.frame with the training data.
#' @param test A data.frame with the test data.
#' @param thr The error threshold. Default is 0.01.
#' @param t_final The number of maximum boosting iterations. Default is 100.
#' @param power Type of loss function, e.g. linear (1), squared (2). Default is 2.
#' @param sigma Regularization factor. Default is 0.5.
#' @param perc.O Percentage for Oversampling via SMOTE, i.e. percentage of extreme cases to be generated. Default is 1.5.
#' @param rel.thr Relevance threshold. Default is 0.9.
#' @param k Number of neighbours used in SMOTE. Defaults to 3.
#' @param coef Coefficient used in boxplot statistics, which is used to create the relevance function. Default is 1.5.
#' @param ... Dots are passed to rpart
#'
#' @return Returns a vector with the predictions made by SMOTEBoost.RTPlus.
#' @export
#'
#' @importFrom MASS ginv
#' @importFrom rpart rpart
#' @importFrom IRon phi.control
#'
#' @examples
#' data(Boston,package="MASS")
#'
#' idx <- sample(1:nrow(Boston),nrow(Boston)*0.75)
#' form <- medv ~ .
#'
#' train <- Boston[idx,]
#' test <- Boston[-idx,]
#'
#' preds <- SMOTEBoost.RTPlus(form,train,test)
#'
SMOTEBoost.RTPlus <- function(form,train,test,t_final=100,thr=0.01,power=2,sigma=0.5,perc.O=1.5,rel.thr=0.9,k=3,coef=1.5,...) {

  require(MASS)
  require(IRon)
  require(rpart)

  models <- list()
  betas <- c()
  pred.mat <- c()
  train_pred.mat <- c()

  ind.y <- which(colnames(train)==as.character(form[[2]]))

  n <- nrow(train) #size of train

  weights <- rep(1/n,n) #initialize weights

  err_t <- 0

  pc <- IRon::phi.control(y = train[,ind.y],method = "extremes",coef=coef)

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

#' RUSBoost variant of AdaBoost.RTPlus
#'
#' @param form The model formula.
#' @param train A data.frame with the training data.
#' @param test A data.frame with the test data.
#' @param thr The error threshold. Default is 0.01.
#' @param t_final The number of maximum boosting iterations. Default is 100.
#' @param power Type of loss function, e.g. linear (1), squared (2). Default is 2.
#' @param sigma Regularization factor. Default is 0.5.
#' @param perc.U Percentage for Undersampling via Random Undersampling, i.e. percentage of cases with normal values to remain in the new dataset. Default is 0.9.
#' @param rel.thr Relevance threshold. Default is 0.9.
#' @param coef Coefficient used in boxplot statistics, which is used to create the relevance function. Default is 1.5.
#' @param ... Dots are passed to rpart
#'
#' @return Returns a vector with the predictions made by RUSBoost.RTPlus.
#' @export
#'
#' @importFrom MASS ginv
#' @importFrom rpart rpart
#' @importFrom IRon phi.control
#' @importFrom UBL RandUnderRegress
#'
#' @examples
#' data(Boston,package="MASS")
#'
#' idx <- sample(1:nrow(Boston),nrow(Boston)*0.75)
#' form <- medv ~ .
#'
#' train <- Boston[idx,]
#' test <- Boston[-idx,]
#'
#' preds <- RUSBoost.RTPlus(form,train,test)
#'
RUSBoost.RTPlus <- function(form,train,test,t_final=100,thr=0.01,power=2,sigma=0.5,perc.U=0.9,rel.thr=0.9,coef=1.5,...) {

  require(MASS)
  require(UBL)
  require(rpart)
  require(IRon)

  models <- list()
  betas <- c()
  pred.mat <- c()
  train_pred.mat <- c()

  ind.y <- which(colnames(train)==as.character(form[[2]]))

  n <- nrow(train) #size of train

  weights <- rep(1/n,n) #initialize weights

  err_t <- 0

  pc <- IRon::phi.control(y = train[,ind.y],method = "extremes",coef=coef)

  for (t in 1:t_final) {

    train.ind <- sample(1:n,n,replace=TRUE,prob=weights)
    new.train <- train[train.ind,]

    new.train <- UBL::RandUnderRegress(form,new.train,pc,rel.thr,perc.U)

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

#' ROSBoost variant of AdaBoost.RTPlus
#'
#' @param form The model formula.
#' @param train A data.frame with the training data.
#' @param test A data.frame with the test data.
#' @param thr The error threshold. Default is 0.01.
#' @param t_final The number of maximum boosting iterations. Default is 100.
#' @param power Type of loss function, e.g. linear (1), squared (2). Default is 2.
#' @param sigma Regularization factor. Default is 0.5.
#' @param perc.O Percentage for Oversampling via Random Oversampling, i.e. percentage of extreme cases to be generated. Default is 1.5.
#' @param rel.thr Relevance threshold. Default is 0.9.
#' @param coef Coefficient used in boxplot statistics, which is used to create the relevance function. Default is 1.5.
#' @param ... Dots are passed to rpart
#'
#' @return Returns a vector with the predictions made by ROSBoost.RTPlus.
#' @export
#'
#' @importFrom MASS ginv
#' @importFrom rpart rpart
#' @importFrom IRon phi.control
#' @importFrom UBL RandOverRegress
#'
#' @examples
#' data(Boston,package="MASS")
#'
#' idx <- sample(1:nrow(Boston),nrow(Boston)*0.75)
#' form <- medv ~ .
#'
#' train <- Boston[idx,]
#' test <- Boston[-idx,]
#'
#' preds <- ROSBoost.RTPlus(form,train,test)
#'
ROSBoost.RTPlus <- function(form,train,test,t_final=100,thr=0.01,power=2,sigma=0.5,perc.O=1.5,rel.thr=0.9,coef=1.5,...) {

  require(MASS)
  require(UBL)
  require(rpart)
  require(IRon)

  models <- list()
  betas <- c()
  pred.mat <- c()
  train_pred.mat <- c()

  ind.y <- which(colnames(train)==as.character(form[[2]]))

  n <- nrow(train) #size of train

  weights <- rep(1/n,n) #initialize weights

  err_t <- 0

  pc <- IRon::phi.control(y = train[,ind.y],method = "extremes",coef=coef)

  for (t in 1:t_final) {

    train.ind <- sample(1:n,n,replace=TRUE,prob=weights)
    new.train <- train[train.ind,]

    new.train <- UBL::RandOverRegress(form,new.train,pc,rel.thr,perc.O)

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

#' BEMBoost
#'
#' @references
#'
#' @param form The model formula.
#' @param train A data.frame with the training data.
#' @param test A data.frame with the test data.
#' @param t_final The number of maximum boosting iterations. Default is 100.
#' @param BEM Biggest error margin admissible. Defaults to 0.5.
#' @param ... Dots are passed to rpart
#'
#' @return Returns a vector with the predictions made by BEMBoost.
#' @export
#'
#' @importFrom rpart rpart
#'
#' @examples
#' data(Boston,package="MASS")
#'
#' idx <- sample(1:nrow(Boston),nrow(Boston)*0.75)
#' form <- medv ~ .
#'
#' train <- Boston[idx,]
#' test <- Boston[-idx,]
#'
#' preds <- BEMBoost(form,train,test)
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

#' SMOTEBoost variant of BEMBoost
#'
#' @param form The model formula.
#' @param train A data.frame with the training data.
#' @param test A data.frame with the test data.
#' @param t_final The number of maximum boosting iterations. Default is 100.
#' @param BEM Biggest error margin admissible. Defaults to 0.5.
#' @param perc.O Percentage for Oversampling via SMOTE, i.e. percentage of extreme cases to be generated. Default is 1.5.
#' @param rel.thr Relevance threshold. Default is 0.9.
#' @param k Number of neighbours used in SMOTE. Defaults to 3.
#' @param coef Coefficient used in boxplot statistics, which is used to create the relevance function. Default is 1.5.
#' @param ... Dots are passed to rpart
#'
#' @return Returns a vector with the predictions made by SMOTEBoost.BEM.
#' @export
#'
#' @importFrom IRon phi.control
#' @importFrom rpart rpart
#'
#' @examples
#' data(Boston,package="MASS")
#'
#' idx <- sample(1:nrow(Boston),nrow(Boston)*0.75)
#' form <- medv ~ .
#'
#' train <- Boston[idx,]
#' test <- Boston[-idx,]
#'
#' preds <- SMOTEBoost.BEM(form,train,test)
#'
SMOTEBoost.BEM <- function(form,train,test,t_final=100,BEM=0.5,perc.O=1.5,rel.thr=0.9,k=3,coef=1.5,...) {

  require(IRon)
  require(rpart)

  models <- list()
  betas <- c()
  pred.mat <- c()

  ind.y <- which(colnames(train)==as.character(form[[2]]))

  n <- nrow(train) #size of train

  weights <- rep(1/n,n) #initialize weights

  err_t <- 0

  pc <- IRon::phi.control(y = train[,ind.y],method = "extremes",coef=coef)

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

#' RUSBoost variant of BEMBoost
#'
#' @param form The model formula.
#' @param train A data.frame with the training data.
#' @param test A data.frame with the test data.
#' @param t_final The number of maximum boosting iterations. Default is 100.
#' @param BEM Biggest error margin admissible. Defaults to 0.5.
#' @param perc.U Percentage for Undersampling via Random Undersampling, i.e. percentage of cases with normal values to remain in the new dataset. Default is 0.9.
#' @param rel.thr Relevance threshold. Default is 0.9.
#' @param coef Coefficient used in boxplot statistics, which is used to create the relevance function. Default is 1.5.
#' @param ... Dots are passed to rpart
#'
#' @return Returns a vector with the predictions made by RUSBoost.BEM.
#' @export
#'
#' @importFrom IRon phi.control
#' @importFrom UBL RandUnderRegress
#' @importFrom rpart rpart
#'
#' @examples
#' data(Boston,package="MASS")
#'
#' idx <- sample(1:nrow(Boston),nrow(Boston)*0.75)
#' form <- medv ~ .
#'
#' train <- Boston[idx,]
#' test <- Boston[-idx,]
#'
#' preds <- RUSBoost.BEM(form,train,test)
#'
RUSBoost.BEM <- function(form,train,test,t_final=100,BEM=0.5,perc.U=0.9,rel.thr=0.9,coef=1.5,...) {

  require(UBL)
  require(rpart)
  require(IRon)

  models <- list()
  betas <- c()
  pred.mat <- c()

  ind.y <- which(colnames(train)==as.character(form[[2]]))

  n <- nrow(train) #size of train

  weights <- rep(1/n,n) #initialize weights

  err_t <- 0

  pc <- IRon::phi.control(y = train[,ind.y],method = "extremes",coef=coef)

  for (t in 1:t_final) {

    train.ind <- sample(1:n,n,replace=TRUE,prob=weights)
    new.train <- train[train.ind,]

    new.train <- UBL::RandUnderRegress(form,new.train,pc,rel.thr,perc.U)

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

#' ROSBoost variant of BEMBoost
#'
#' @param form The model formula.
#' @param train A data.frame with the training data.
#' @param test A data.frame with the test data.
#' @param t_final The number of maximum boosting iterations. Default is 100.
#' @param BEM Biggest error margin admissible. Defaults to 0.5.
#' @param perc.O Percentage for Oversampling via Random Oversampling, i.e. percentage of extreme cases to be generated. Default is 1.5.
#' @param rel.thr Relevance threshold. Default is 0.9.
#' @param coef Coefficient used in boxplot statistics, which is used to create the relevance function. Default is 1.5.
#' @param ... Dots are passed to rpart
#'
#' @return Returns a vector with the predictions made by ROSBoost.BEM.
#' @export
#'
#' @importFrom IRon phi.control
#' @importFrom UBL RandOverRegress
#' @importFrom rpart rpart
#'
#' @examples
#' data(Boston,package="MASS")
#'
#' idx <- sample(1:nrow(Boston),nrow(Boston)*0.75)
#' form <- medv ~ .
#'
#' train <- Boston[idx,]
#' test <- Boston[-idx,]
#'
#' preds <- ROSBoost.BEM(form,train,test)
#'
ROSBoost.BEM <- function(form,train,test,t_final=100,BEM=0.5,perc.O=1.5,rel.thr=0.9,coef=1.5,...) {

  require(UBL)
  require(rpart)
  require(IRon)

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

    new.train <- UBL::RandOverRegress(form,new.train,pc,rel.thr,perc.O)

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
