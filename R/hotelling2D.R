hotelling2D <- function(dat, alpha=0.1, steps=400, searchwidth=8, equi=1.25, plotrange=c(0.77, 1.3), axnames=NULL, main="Hotelling"){
  
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
    togrid[[i]] <- seq(est[i] - searchwidth * poolvar, est[i] + searchwidth * poolvar, length.out=steps)
  }
  
  grid <- expand.grid(togrid)
  
  findcrHot <- apply(grid, 1, function(x){
    theta <- matrix(x, 2)
    (n * t(est - theta) %*% solve(cov) %*% (est - theta)) < qf(p=1 - alpha, df1=2, df2=df - 1) * 2 * df / (df - 1)
  })
  
  crHot <- cbind(grid, findcrHot)[findcrHot==1, ]
  
  if(min(crHot[, 1])==min(grid[, 1]) | max(crHot[, 1])==max(grid[, 1]) |
       min(crHot[, 2])==min(grid[, 2]) | max(crHot[, 2])==max(grid[, 2])){
    warning("The search grid is too narrow, please increase searchwidth.")
  }
  
  if(is.null(axnames)==TRUE){
    axisnames <- colnames(dat)
  }else{
    axisnames <- axnames
  }
  
  par(mar=c(5, 5, 4, 2))
  plot(0, xlim=log(plotrange), ylim=log(plotrange), las=1, xlab=axisnames[1], ylab=axisnames[2],
       cex.main=2.5, cex.axis=1.5, cex.lab=1.7, main=main)
  rect(log(1/equi), log(1/equi), log(equi), log(equi), col="gray95", border=NA)
  points(crHot[, -3], pch=20)
  points(est[1], est[2], pch=19, col="white")
  par(mar=c(5, 4, 4, 2))
  
}
