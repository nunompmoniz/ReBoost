#' Title
#'
#' @param form The model formula.
#' @param dat A data.frame with the training data.
#' @param perc.o Percentage for Oversampling via SMOTE, i.e. percentage of extreme cases to be generated. Default is 1.5.
#' @param rel.thr Relevance threshold. Default is 0.9.
#' @param k Number of neighbours used in SMOTE. Defaults to 3.
#' @param pc Relevance function phi.
#'
#' @return
#' @export
#'
#' @importFrom UBL phi.control SmoteRegress
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
    pc <- UBL::phi.control(y = y,method = "extremes",coef=1.5)
  }

  new.dat <- c()

  if(any(pc$control.pts[c(2,8)]==1)) {

    percs <- list()

    if(pc$control.pts[2]==1) {
      if(length(adjboxStats(y)$out<=pc$control.pts[1]) > 1) {
        percs <- c(percs,perc.o)
      } else {
        percs <- c(percs,1)
      }
    }

    percs <- c(percs,1)

    if(pc$control.pts[8]==1) {
      if(length(adjboxStats(y)$out>=pc$control.pts[7]) > 1) {
        percs <- c(percs,perc.o)
      } else {
        percs <- c(percs,1)
      }
    }

    if(any(sapply(dat,is.numeric)==FALSE)){ #If there's any nominal predictor, use HEOM distance

      new.dat <- UBL::SmoteRegress(form,dat,rel=pc,thr.rel=rel.thr,C.perc=percs,k=k,dist="HEOM")

    } else { #If all predictors are numerical, use Euclidean distance

      new.dat <- UBL::SmoteRegress(form,dat,rel=pc,thr.rel=rel.thr,C.perc=percs,k=k,dist="Euclidean")

    }

  } else {
    warning("Did not found any extreme cases. Returning the original train set.")
  }

  new.dat

}
