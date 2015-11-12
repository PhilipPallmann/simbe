tseng2D <- function(dat, alpha=0.1, steps=400, equi=1.25, plotrange=c(0.77, 1.3), axnames=NULL, main="Tseng"){
  
  if(ncol(dat)!=2){
    stop("Data must be bivariate.")
  }
  
  n <- nrow(dat)
  df <- n - 1
  
  est <- matrix(colMeans(dat), 2)
  poolvar <- var(as.vector(as.matrix(dat)))
  s2 <- poolvar / n
  
  togrid <- list()
  
  for(i in 1:2){
    togrid[[i]] <- seq(est[i] - 8 * poolvar, est[i] + 8 * poolvar, length.out=steps)
  }
  
  grid <- expand.grid(togrid)
  
  findcrTse <- apply(grid, 1, function(x){
    theta <- x
    (sqrt(sum(est^2))^2 / (2 * s2)) > qf(p=alpha, df1=2, df2=df, ncp=(sqrt(sum(theta^2))^2 / s2))
  })
  
  crTse <- cbind(grid, findcrTse)[findcrTse==1, ]
  
  if(is.null(axnames)==TRUE){
    axisnames <- colnames(dat)
  }else{
    axisnames <- axnames
  }
    
  plot(0, xlim=log(plotrange), ylim=log(plotrange), las=1, xlab=axisnames[1], ylab=axisnames[2],
       cex.main=2.5, cex.axis=1.5, cex.lab=1.5, main=main)
  rect(log(1/equi), log(1/equi), log(equi), log(equi), col="gray95", border=NA)
  points(crTse[, -3], pch=20)
  points(est[1], est[2], pch=19, col="white")
  
}
