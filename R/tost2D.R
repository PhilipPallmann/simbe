tost2D <- function(dat, alpha=0.1, equi=1.25, plotrange=c(0.77, 1.3), axnames=NULL, main="TOST"){
  
  if(ncol(dat)!=2){
    stop("Data must be bivariate.")
  }
  
  n <- nrow(dat)
    
  est <- matrix(colMeans(dat), 2)
  sd <- c(sd(dat[, 1]), sd(dat[, 2]))
  
  ci <- c(est - sd * qt(1 - alpha, n - 1) / sqrt(n),
          est + sd * qt(1 - alpha, n - 1) / sqrt(n))
  
  if(is.null(axnames)==TRUE){
    axisnames <- colnames(dat)
  }else{
    axisnames <- axnames
  }
  
  par(mar=c(5, 5, 4, 2))
  plot(0, xlim=log(plotrange), ylim=log(plotrange), las=1, xlab=axisnames[1], ylab=axisnames[2],
       cex.main=2.5, cex.axis=1.5, cex.lab=1.7, main=main)
  rect(log(1/equi), log(1/equi), log(equi), log(equi), col="gray95", border=NA)
  segments(x0=ci[1], x1=ci[3], y0=est[2], y1=est[2])
  segments(y0=ci[2], y1=ci[4], x0=est[1], x1=est[1])
  par(mar=c(5, 4, 4, 2))
  
}