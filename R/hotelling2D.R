hotelling2D <- function(dat, alpha=0.1, steps=400, equi=1.25, plotrange=c(0.77, 1.3), axnames=NULL){
  
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
  
  findcrHot <- apply(grid, 1, function(x){
    theta <- matrix(x, 2)
    (n * t(est - theta) %*% solve(cov) %*% (est - theta)) < qf(p=alpha, df1=2, df2=df - 1) * 2 * df / (df - 1)
  })
  
  crHot <- cbind(grid, findcrHot)[findcrHot==1, ]
  
  if(is.null(axnames)==TRUE){
    axisnames <- colnames(dat)
  }else{
    axisnames <- axnames
  }
  
  plot(0, xlim=log(plotrange), ylim=log(plotrange), las=1, xlab=axisnames[1], ylab=axisnames[2],
       cex.main=2.5, cex.axis=1.5, cex.lab=1.5, main="Hotelling")
  rect(log(1/equi), log(1/equi), log(equi), log(equi), col="gray95", border=NA)
  points(crHot[, -3], pch=20)
  points(est, pch=19, col="white")
  
}
