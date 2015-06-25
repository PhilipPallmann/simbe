tseng <- function(dat, alpha=0.1, steps=100){
  
  n <- nrow(dat)
  p <- ncol(dat)
  df <- n - 1
  
  est <- matrix(colMeans(dat), p)
  poolvar <- var(as.vector(as.matrix(dat)))
  s2 <- poolvar / n
  
  part1 <- sqrt(sum(est^2))^2 / (p * s2)
  
  togrid <- list()
  
  for(i in 1:p){
    togrid[[i]] <- seq(est[i] - 2 * poolvar, est[i] + 2 * poolvar, length.out=steps)
  }
  
  grid <- expand.grid(togrid)
  
  findcrTse <- apply(grid, 1, function(x){
    theta <- x
    part1 > qf(p=alpha, df1=p, df2=df, ncp=(sqrt(sum(theta^2))^2 / s2))
  })
  
  crTse <- cbind(grid, findcrTse)[findcrTse==1, ]
  
  if(nrow(crTse)==0){
    Tse0 <- matrix(rep(Inf, 2 * p), 2)
  }else{
    Tse0 <- t(apply(crTse[, -(p + 1)], 2, range))
  }
  
  if(Tse0[1, 1]==Inf | Tse0[1, 1]==-Inf){
    
    Tse <- matrix(rep(0, 2 * p), 2)
    
  }else{
    
    if(min(abs(Tse0[, 1] - est + 2 * poolvar)) < 0.001 | min(abs(Tse0[, 2] - est - 2 * poolvar)) < 0.001){
      
      togrid2 <- list()
      
      for(i in 1:p){
        togrid2[[i]] <- seq(est[i] - 8 * poolvar, est[i] + 8 * poolvar, length.out=4 * steps)
      }
      
      grid2 <- expand.grid(togrid2)
      
      findcrTse2 <- apply(grid2, 1, function(x){
        theta <- x
        part1 > qf(p=alpha, df1=p, df2=df, ncp=(sqrt(sum(theta^2))^2 / s2))
      })
      
      crTse2 <- cbind(grid2, findcrTse2)[findcrTse2==1, ]
      
      Tse0 <- t(apply(crTse2[, -(p + 1)], 2, range))
      
    }
    
    stepwidth <- 4 * poolvar / steps
    
    if(p==2){
      
      GridA <- expand.grid(seq(Tse0[1, 1] - stepwidth, Tse0[1, 1], length.out=steps),
                           seq(Tse0[2, 1] - stepwidth, Tse0[2, 2], length.out=steps))
      
      GridB <- expand.grid(seq(Tse0[1, 2], Tse0[1, 2] + stepwidth, length.out=steps),
                           seq(Tse0[2, 1], Tse0[2, 2] + stepwidth, length.out=steps))
      
      GridC <- expand.grid(seq(Tse0[1, 1] - stepwidth, Tse0[1, 2], length.out=steps),
                           seq(Tse0[2, 1] - stepwidth, Tse0[2, 1], length.out=steps))
      
      GridD <- expand.grid(seq(Tse0[1, 1], Tse0[1, 2] + stepwidth, length.out=steps),
                           seq(Tse0[2, 2], Tse0[2, 2] + stepwidth, length.out=steps))
      
      Grid <- rbind(GridA, GridB, GridC, GridD)
      
      FindcrTse <- apply(Grid, 1, function(x){
        theta <- x
        part1 > qf(p=alpha, df1=p, df2=df, ncp=(sqrt(sum(theta^2))^2 / s2))
      })
      
      CrTse <- cbind(Grid, FindcrTse)[FindcrTse==1, ]
      
      Tse <- t(apply(CrTse[, -(p + 1)], 2, range))
      
    }else{
      
      Tse <- Tse0
      
    }
        
  }
  
  TseOut <- cbind(est, Tse)
  rownames(TseOut) <- colnames(dat)
  colnames(TseOut) <- c("estimate", "lower", "upper")
  
  return(TseOut)
  
}
