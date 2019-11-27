#' Synthetic Minority Oversampling Technique
#'
#' @param form The model formula.
#' @param dat A data.frame with the training data.
#' @param perc.o Percentage for Oversampling via SMOTE, i.e. percentage of extreme cases to be generated. Default is 1.5.
#' @param rel.thr Relevance threshold. Default is 0.9.
#' @param k Number of neighbours used in SMOTE. Defaults to 3.
#' @param pc Relevance function phi.
#'
#' @return Return a data frame with the data set results from the application of the SMOTE strategy
#' @export
#'
#' @importFrom IRon phi.control
#' @importFrom UBL SmoteRegress
#' @importFrom robustbase adjboxStats
#'
#' @examples
#'
#' require(robustbase)
#'
#' data(Boston,package="MASS")
#'
#' idx <- sample(1:nrow(Boston),nrow(Boston)*0.75)
#' form <- medv ~ .
#'
#' train <- Boston[idx,]
#'
#' new.train <- adaSMOTE(form,train)
#'
#' adjbox(train$medv)
#' adjbox(new.train$medv)
adaSMOTE <- function(form,dat,perc.o=1.5,rel.thr,k,pc=NULL) {

  require(UBL)

  y <- dat[,as.character(form[[2]])]

  if(length(pc)!=3) {
    pc <- IRon::phi.control(y = y,method = "extremes",coef=1.5)
  }

  pc.out <- pc.m <- matrix(pc$control.pts,nrow = (length(pc$control.pts)/3),ncol = 3,byrow = TRUE)

  new.dat <- c()

  if(any(pc.m[,2]==1)) {

    percs <- list()

    if(nrow(pc.out)==3) {

      if(pc.out[1,2]==1) {

        if(sum(y<=pc.out[1,1]) >= 1) {

          percs <- c(percs,perc.o)

        } else {

          pc.m <- pc.m[-1,]

        }
      }

      percs <- c(percs,1)

      if(pc.out[3,2]==1) {

        if(sum(y>=pc.out[3,1]) >= 1) {

          percs <- c(percs,perc.o)

        } else {

          pc.m <- pc.m[-3,]
        }
      }
    } else {

      if(pc.m[1,2]==1) {

        if(sum(y<=pc.m[1,1]) >= 1) {

          percs <- c(percs,perc.o)

        } else {

          percs <- c(percs,1)
        }

      } else {

        percs <- c(percs,1)
      }

      if(pc.m[2,2]==1) {

        if(sum(y>=pc.m[2,1]) >= 1) {

          percs <- c(percs,perc.o)

        } else {

          percs <- c(percs,1)

        }

      } else {

        percs <- c(percs,1)
      }

    }

    if(length(percs)>1) {

      #' Snippet: SMOTE requires distinct cases. If there's only one case and its repetitions,
      #' We randomly select a numerical column (except for the target) and add Gaussian noise (sd=0.001)

      if(nrow(pc.out)==3) {

        dat.high <- dat[which(y>=pc$control.pts[7]),]
        dat.low <- dat[which(y<=pc$control.pts[1]),]
        tgt <- which(colnames(dat)==as.character(form[[2]]))

        num.colname <- names(which(sapply(dat[,-tgt],is.numeric)))
        rnd.col <- as.numeric(which(colnames(dat)==sample(num.colname,1)))

        if(nrow(unique(dat.high))==1 & nrow(dat.high)>1) {

          dat[which(y>=pc$control.pts[7]),rnd.col] <- rnorm(n=nrow(dat.high),mean = dat.high[1,rnd.col],sd = 0.001)
          dat[which(y>=pc$control.pts[7]),tgt] <- rnorm(n=nrow(dat.high),mean = dat.high[1,tgt],sd = 0.001)
        }

        if(nrow(unique(dat.low))==1 & nrow(dat.low)>1) {

          dat[which(y<=pc$control.pts[1]),rnd.col] <- rnorm(n=nrow(dat.low),mean = dat.low[1,rnd.col],sd = 0.001)
          dat[which(y<=pc$control.pts[1]),tgt] <- rnorm(n=nrow(dat.low),mean = dat.low[1,tgt],sd = 0.001)
        }

        rm(dat.high,dat.low,rnd.col)

      } else {

      }

      #' End of snippet

      if(any(sapply(dat,is.numeric)==FALSE)){ #If there's any nominal predictor, use HEOM distance

        new.dat <- UBL::SmoteRegress(form,dat,rel=pc.m,thr.rel=rel.thr,C.perc=percs,k=k,dist="HEOM")

      } else { #If all predictors are numerical, use Euclidean distance

        new.dat <- UBL::SmoteRegress(form,dat,rel=pc.m,thr.rel=rel.thr,C.perc=percs,k=k,dist="Euclidean")

      }

    } else {

      warning("Did not found any extreme cases. Returning the original train set.")
      new.dat <- dat

    }


  } else {

    new.dat <- dat
    warning("Did not found any extreme cases. Returning the original train set.")
  }

  new.dat

}
