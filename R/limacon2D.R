casella2D <- function(dat, alpha=0.1, steps=400, equi=1.25, plotrange=c(0.77, 1.3), axnames=NULL){
  
  if(ncol(dat)!=2){
    stop("Data must be bivariate.")
  }
  
  n <- nrow(dat)
  df <- n - 1
  
  est <- matrix(colMeans(dat), 2)
  poolvar <- var(as.vector(as.matrix(dat)))
  cov <- cov(dat)
  
  togrid <- list()
  
  for(i in 1:2){
    togrid[[i]] <- seq(est[i] - 8 * poolvar, est[i] + 8 * poolvar, length.out=steps)
  }
  
  grid <- expand.grid(togrid)
  
  findcrLim <- apply(grid, 1, function(x){
    theta <- matrix(x, 2)
    ((t(theta) %*% est) / sqrt((t(theta) %*% cov %*% theta) / n) + qt(1 - alpha, df)) >
      ((t(theta) %*% theta) / sqrt((t(theta) %*% cov %*% theta) / n))
  })
  
  crLim <- cbind(grid, findcrLim)[findcrLim==1, ]
  
  if(is.null(axnames)==TRUE){
    axisnames <- colnames(dat)
  }else{
    axisnames <- axnames
  }
  
  plot(0, xlim=log(plotrange), ylim=log(plotrange), las=1, xlab=axisnames[1], ylab=axisnames[2],
       cex.main=2.5, cex.axis=1.5, cex.lab=1.5, main="Limacon")
  rect(log(1/equi), log(1/equi), log(equi), log(equi), col="gray95", border="gray95")
  points(crLim[, -3])
  points(est, pch=16, col="white")
  
}
