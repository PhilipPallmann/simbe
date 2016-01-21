bootkern2D <- function(dat, alpha=0.1, nboot=10000, equi=1.25, plotrange=c(0.77, 1.3), axnames=NULL, main="Bootstrap"){
  
  if(ncol(dat)!=2){
    stop("Data must be bivariate.")
  }
  
  est <- matrix(colMeans(dat), 2)
  
  bivarmean <- function(x, d) {
    e <- x[d, ]
    return(c(mean(e[, 1]), mean(e[, 2])))
  } # werden 1. und 2. Spalte immer zusammen gesampelt?!
  
  b <- boot::boot(dat, bivarmean, R=nboot)
  bdat <- as.data.frame(b$t)
  
  kern <- KernSmooth::bkde2D(bdat, bandwidth=sapply(bdat, KernSmooth::dpik))
  
  alphasum <- sum(kern$fhat) * alpha
  
  while(sum(kern$fhat, na.rm=TRUE) > alphasum){
    kern$fhat[which.max(kern$fhat)] <- NA
  }
  
  KERN <- is.na(kern$fhat)
  K <- expand.grid(kern$x1, kern$x2)
  K$truefalse <- as.vector(KERN)
  K <- K[K$truefalse==1, ]
  
  cint <- rbind(range(K$Var1), range(K$Var2))
  
  hu <- chull(K[, -3])
  hull <- c(hu, hu[1])
  
  if(is.null(axnames)==TRUE){
    axisnames <- colnames(dat)
  }else{
    axisnames <- axnames
  }
  
  plot(0, xlim=log(plotrange), ylim=log(plotrange), las=1, xlab=axisnames[1], ylab=axisnames[2],
       cex.main=2.5, cex.axis=1.5, cex.lab=1.5, main=main, mar=c(5, 5, 4, 2))
  rect(log(1/equi), log(1/equi), log(equi), log(equi), col="gray95", border=NA)
  polygon(K[hull, ], col="black")
  points(est[1], est[2], pch=19, col="white")
  
}
