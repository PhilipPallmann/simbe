casella <- function(dat, alpha=0.1, steps=100){
  
  n <- nrow(dat)
  p <- ncol(dat)
  df <- n - 1
  a <- (p - 2) * (df / (df + 2))
  
  est <- matrix(colMeans(dat), p)
  poolvar <- var(as.vector(as.matrix(dat)))
  s <- sqrt(poolvar / n)
  
  JSfactor <- (1 - (a * (poolvar/n)) / sum(est^2)) # or just poolvar???
  
  if(JSfactor < 0){
    JSplus <- est
  }else{
    JSplus <- JSfactor * est
  }
  
  cond <- (sum(est^2) / (poolvar/n)) < (p * qf(p=1 - alpha, df1=p, df2=df)) # or just poolvar???
  
  if(cond==TRUE){
    vE2 <- (1 - a/(p * qf(p=1 - alpha, df1=p, df2=df))) *
      (p * qf(p=1 - alpha, df1=p, df2=df) - p * log(1 - a/(p * qf(p=1 - alpha, df1=p, df2=df))))
  }else{
    vE2 <- (1 - a/(p * qf(p=1 - alpha, df1=p, df2=df))) *
      (p * qf(p=1 - alpha, df1=p, df2=df) - p * log(1 - a/(p * qf(p=1 - alpha, df1=p, df2=df))))
  }
  
  togrid <- list()
  
  for(i in 1:p){
    togrid[[i]] <- seq(est[i] - 2 * poolvar, est[i] + 2 * poolvar, length.out=steps)
  }
  
  grid <- expand.grid(togrid)
  
  findcrCas <- apply(grid, 1, function(x){
    theta <- matrix(x, p)
    sqrt(sum((theta - JSplus)^2)) < (s * sqrt(vE2))
  })
  
  crCas <- cbind(grid, findcrCas)[findcrCas==1, ]
  
  Cas0 <- t(apply(crCas[, -(p + 1)], 2, range, na.rm=TRUE))
  
  if(min(abs(Cas0[, 1] - est + 2 * poolvar)) < 0.001 | min(abs(Cas0[, 2] - est - 2 * poolvar)) < 0.001){
    
    togrid2 <- list()
    
    for(i in 1:p){
      togrid2[[i]] <- seq(est[i] - 8 * poolvar, est[i] + 8 * poolvar, length.out=4 * steps)
    }
    
    grid2 <- expand.grid(togrid2)
    
    findcrCas2 <- apply(grid2, 1, function(x){
      theta <- matrix(x, p)
      sqrt(sum((theta - JSplus)^2)) < (s * sqrt(vE2))
    })
    
    crCas2 <- cbind(grid2, findcrCas2)[findcrCas2==1, ]
    
    Cas0 <- t(apply(crCas2[, -(p + 1)], 2, range))
    
  }
  
  stepwidth <- 4 * poolvar / steps
  
  if(p==2){
    
    GridA <- expand.grid(seq(Cas0[1, 1] - stepwidth, Cas0[1, 1], length.out=steps),
                         seq(Cas0[2, 1] - stepwidth, Cas0[2, 2], length.out=steps))
    
    GridB <- expand.grid(seq(Cas0[1, 2], Cas0[1, 2] + stepwidth, length.out=steps),
                         seq(Cas0[2, 1], Cas0[2, 2] + stepwidth, length.out=steps))
    
    GridC <- expand.grid(seq(Cas0[1, 1] - stepwidth, Cas0[1, 2], length.out=steps),
                         seq(Cas0[2, 1] - stepwidth, Cas0[2, 1], length.out=steps))
    
    GridD <- expand.grid(seq(Cas0[1, 1], Cas0[1, 2] + stepwidth, length.out=steps),
                         seq(Cas0[2, 2], Cas0[2, 2] + stepwidth, length.out=steps))
    
    Grid <- rbind(GridA, GridB, GridC, GridD)
    
    FindcrCas <- apply(Grid, 1, function(x){
      theta <- matrix(x, 2)
      sqrt(sum((theta - JSplus)^2)) < (s * sqrt(vE2))
    })
    
    CrCas <- cbind(Grid, FindcrCas)[FindcrCas==1, ]
    
    Cas <- t(apply(CrCas[, -(p + 1)], 2, range))
    
  }else{
    
    Cas <- Cas0
    
  }
  
  CasOut <- cbind(JSplus, Cas)
  rownames(CasOut) <- colnames(dat)
  colnames(CasOut) <- c("estimate", "lower", "upper")
  
  return(CasOut)
  
}
