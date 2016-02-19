limasym <- function(dat, alpha=0.1, steps=100){
  
  n <- nrow(dat)
  p <- ncol(dat)
  df <- n - 1
  
  est <- matrix(colMeans(dat), p)
  poolvar <- var(as.vector(as.matrix(dat)))
  cov <- cov(dat)
  
  togrid <- list()
  
  for(i in 1:p){
    togrid[[i]] <- seq(est[i] - 2 * poolvar, est[i] + 2 * poolvar, length.out=steps)
  }
  
  grid <- expand.grid(togrid)
  
  findcrLim <- apply(grid, 1, function(x){
    theta <- matrix(x, p)
    ((t(theta) %*% solve(cov) %*% est) / sqrt(t(theta) %*% solve(cov) %*% theta) * sqrt(n) + qnorm(1 - alpha)) >
      (sqrt(t(theta) %*% solve(cov) %*% theta) * sqrt(n))
  })
  
  crLim <- cbind(grid, findcrLim)[findcrLim==1, ]
  
  Lim0 <- t(apply(crLim[, -(p + 1)], 2, range, na.rm=TRUE))
  
  if(min(abs(Lim0[, 1] - est + 2 * poolvar)) < 0.001 | min(abs(Lim0[, 2] - est - 2 * poolvar)) < 0.001){
    
    togrid2 <- list()
    
    for(i in 1:p){
      togrid2[[i]] <- seq(est[i] - 8 * poolvar, est[i] + 8 * poolvar, length.out=4 * steps)
    }
    
    grid2 <- expand.grid(togrid2)
    
    findcrLim2 <- apply(grid2, 1, function(x){
      theta <- matrix(x, p)
      ((t(theta) %*% solve(cov) %*% est) / sqrt(t(theta) %*% solve(cov) %*% theta) * sqrt(n) + qnorm(1 - alpha)) >
        (sqrt(t(theta) %*% solve(cov) %*% theta) * sqrt(n))
    })
    
    crLim2 <- cbind(grid2, findcrLim2)[findcrLim2==1, ]
    
    if(min(crLim2[, 1])==min(grid2[, 1]) | max(crLim2[, 1])==max(grid2[, 1]) |
         min(crLim2[, 2])==min(grid2[, 2]) | max(crLim2[, 2])==max(grid2[, 2])){
      
      togrid3 <- list()
      
      for(i in 1:p){
        togrid3[[i]] <- seq(est[i] - 16 * poolvar, est[i] + 16 * poolvar, length.out=8 * steps)
      }
      
      grid3 <- expand.grid(togrid3)
      
      findcrLim3 <- apply(grid3, 1, function(x){
        theta <- matrix(x, p)
        ((t(theta) %*% solve(cov) %*% est) / sqrt(t(theta) %*% solve(cov) %*% theta) * sqrt(n) + qnorm(1 - alpha)) >
          (sqrt(t(theta) %*% solve(cov) %*% theta) * sqrt(n))
      })
      
      crLim3 <- cbind(grid3, findcrLim3)[findcrLim3==1, ]
      
      crLim2 <- crLim3
      
    }
    
    Lim0 <- t(apply(crLim2[, -(p + 1)], 2, range))
    
  }
  
  stepwidth <- 4 * poolvar / steps
  
  if(p==2){
    
    GridA <- expand.grid(seq(Lim0[1, 1] - stepwidth, Lim0[1, 1], length.out=steps),
                         seq(Lim0[2, 1] - stepwidth, Lim0[2, 2], length.out=steps))
    
    GridB <- expand.grid(seq(Lim0[1, 2], Lim0[1, 2] + stepwidth, length.out=steps),
                         seq(Lim0[2, 1], Lim0[2, 2] + stepwidth, length.out=steps))
    
    GridC <- expand.grid(seq(Lim0[1, 1] - stepwidth, Lim0[1, 2], length.out=steps),
                         seq(Lim0[2, 1] - stepwidth, Lim0[2, 1], length.out=steps))
    
    GridD <- expand.grid(seq(Lim0[1, 1], Lim0[1, 2] + stepwidth, length.out=steps),
                         seq(Lim0[2, 2], Lim0[2, 2] + stepwidth, length.out=steps))
    
    Grid <- rbind(GridA, GridB, GridC, GridD)
    
    FindcrLim <- apply(Grid, 1, function(x){
      theta <- matrix(x, 2)
      ((t(theta) %*% solve(cov) %*% est) / sqrt(t(theta) %*% solve(cov) %*% theta) * sqrt(n) + qnorm(1 - alpha)) >
        (sqrt(t(theta) %*% solve(cov) %*% theta) * sqrt(n))
    })
    
    CrLim <- cbind(Grid, FindcrLim)[FindcrLim==1, ]
    
    Lim <- t(apply(CrLim[, -(p + 1)], 2, range))
    
  }else{
    
    Lim <- Lim0
    
  }
  
  LimOut <- cbind(est, Lim)
  rownames(LimOut) <- colnames(dat)
  colnames(LimOut) <- c("estimate", "lower", "upper")
  
  return(LimOut)
  
}
