adaSMOTE <- function(form,dat,perc.o,rel.thr,k,pc=NULL) {

  require(UBL)

  y <- dat[,as.character(form[[2]])]

  if(length(pc)!=3) {
    pc <- UBL::phi.control(y = y,method = "extremes",coef=1.5)
  }

  new.dat <- c()

  if(any(pc$control.pts[c(2,8)]==1)) {

    percs <- list()

    if(pc$control.pts[2]==1) {
      if(length(boxplot.stats(y)$out<=pc$control.pts[1]) > 1) {
        percs <- c(percs,perc.o)
      } else {
        percs <- c(percs,1)
      }
    }

    percs <- c(percs,1)

    if(pc$control.pts[8]==1) {
      if(length(boxplot.stats(y)$out>=pc$control.pts[7]) > 1) {
        percs <- c(percs,perc.o)
      } else {
        percs <- c(percs,1)
      }
    }

    if(any(sapply(dat,is.numeric)==FALSE)){ #If there's any nominal predictor, use HEOM distance

      new.dat <- SmoteRegress(form,dat,rel=pc,thr.rel=rel.thr,C.perc=percs,k=k,dist="HEOM")

    } else { #If all predictors are numerical, use Euclidean distance

      new.dat <- SmoteRegress(form,dat,rel=pc,thr.rel=rel.thr,C.perc=percs,k=k,dist="Euclidean")

    }

  } else {
    warning("Did not found any extreme cases. Returning the original train set.")
  }

  new.dat

}
