standardcorr2D <- function(dat, alpha=0.1, steps=400, searchwidth=8, equi=1.25, plotrange=c(0.77, 1.3), axnames=NULL,
                           main="Standard (corr.)"){
  
  if(ncol(dat)!=2){
    stop("Data must be bivariate.")
  }
  
  n <- nrow(dat)
  df <- n - 1
  
  est <- matrix(colMeans(dat), 2)
  LambdaInv <- solve(cov(dat))
  poolvar <- var(as.vector(as.matrix(dat)))
  #s2 <- poolvar / n
  
  togrid <- list()
  
  for(i in 1:2){
    togrid[[i]] <- seq(est[i] - searchwidth * poolvar, est[i] + searchwidth * poolvar, length.out=steps)
  }
  
  grid <- expand.grid(togrid)
  
  findcrUsu <- apply(grid, 1, function(x){
    theta <- x
    #sqrt(sum((est - theta)^2))^2 < (s2 * 2 * qf(p=1 - alpha, df1=2, df2=n - 1))
    t(est - theta) %*% LambdaInv %*% (est - theta) < (2 / (n - 1) * qf(p=1 - alpha, df1=2, df2=n - 1))
  })
  
  crUsu <- cbind(grid, findcrUsu)[findcrUsu==1, ]
  
  if(min(crUsu[, 1])==min(grid[, 1]) | max(crUsu[, 1])==max(grid[, 1]) |
       min(crUsu[, 2])==min(grid[, 2]) | max(crUsu[, 2])==max(grid[, 2])){
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
  points(crUsu[, -3], pch=20)
  points(est[1], est[2], pch=19, col="white")
  par(mar=c(5, 4, 4, 2))
  
}
