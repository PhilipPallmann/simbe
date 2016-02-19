he <- function(dat, alpha=0.1, denominator="p-1", steps=1e4, mini=-3, maxi=3){
  
  n <- nrow(dat)
  p <- ncol(dat)
  df <- n - 1
  
  est <- matrix(colMeans(dat), p)
  poolvar <- var(as.vector(as.matrix(dat)))
    
  deltaplus <- max((1 - (p - 2) / sum(est^2)), 0) * est # * poolvar
  
  denominator <- match.arg(denominator, choices=c("p-1", "p-0.5", "p"))
  
  if(denominator=="p-1"){
    denom <- p - 1
  }
  
  if(denominator=="p-0.5"){
    denom <- p - 0.5
  }
  
  if(denominator=="p"){
    denom <- p
  }
  
  if(sum(est^2) > (p - 1)){
    v2 <- (1 - (p - 2)/sum(est^2)) * (qchisq(1 - alpha, p) - log(1 - (p - 2)/sum(est^2)))
  }else{
    v2 <- (1 - (p - 2)/denom) * (qchisq(1 - alpha, p) - log(1 - (p - 2)/denom))
  }
  
  grid <- matrix(NA, ncol=p, nrow=steps)
  minmax <- matrix(NA, ncol=p, nrow=2)
  
  for(i in 1:p){
    for(s in 1:steps){
      x <- mini + (s - 1) * ((maxi - mini)/steps)
      grid[s, i] <- sqrt(sum((x - deltaplus[i])^2)) < sqrt(v2)
    }
    minmax[1, i] <- (min(which(grid[, i]==TRUE)) - 1) * (maxi - mini)/steps + mini
    minmax[2, i] <- (max(which(grid[, i]==TRUE)) - 1) * (maxi - mini)/steps + mini
  }
  
  if(sum(grid[1, ]) > 0 | sum(grid[steps, ]) > 0) warning("Choose a wider search window.")
  
  HeOut <- cbind(deltaplus, minmax[1, ], minmax[2, ])
  rownames(HeOut) <- colnames(dat)
  colnames(HeOut) <- c("estimate", "lower", "upper")
  
  return(HeOut)
  
}