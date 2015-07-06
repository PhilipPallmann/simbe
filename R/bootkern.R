bootkern <- function(dat, alpha=0.1, nboot=10000){
  
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
  
  Boo <- rbind(range(K$Var1), range(K$Var2))
  
  BooOut <- cbind(est, Boo)
  rownames(BooOut) <- colnames(dat)
  colnames(BooOut) <- c("estimate", "lower", "upper")
  
  return(BooOut)
  
}
