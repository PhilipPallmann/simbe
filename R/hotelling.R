hotelling <- function(dat, alpha=0.1, steps=100){
  
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
  
  findcrHot <- apply(grid, 1, function(x){
    theta <- matrix(x, p)
    (n * t(est - theta) %*% solve(cov) %*% (est - theta)) <
      qf(p=alpha, df1=2, df2=df - p + 1) * p * df / (df - p + 1)
  })
  
  crHot <- cbind(grid, findcrHot)[findcrHot==1, ]
  
  Hot0 <- t(apply(crHot[, -(p + 1)], 2, range, na.rm=TRUE))
  
  if(min(abs(Hot0[, 1] - est + 2 * poolvar)) < 0.001 | min(abs(Hot0[, 2] - est - 2 * poolvar)) < 0.001){
    
    togrid2 <- list()
    
    for(i in 1:p){
      togrid2[[i]] <- seq(est[i] - 8 * poolvar, est[i] + 8 * poolvar, length.out=4 * steps)
    }
    
    grid2 <- expand.grid(togrid2)
    
    findcrHot2 <- apply(grid2, 1, function(x){
      theta <- matrix(x, p)
      (n * t(est - theta) %*% solve(cov) %*% (est - theta)) <
        qf(p=alpha, df1=2, df2=df - p + 1) * p * df / (df - p + 1)
    })
    
    crHot2 <- cbind(grid2, findcrHot2)[findcrHot2==1, ]
    
    Hot0 <- t(apply(crHot2[, -(p + 1)], 2, range))
    
  }
  
  stepwidth <- 4 * poolvar / steps
  
  if(p==2){
    
    GridA <- expand.grid(seq(Hot0[1, 1] - stepwidth, Hot0[1, 1], length.out=steps),
                         seq(Hot0[2, 1] - stepwidth, Hot0[2, 2], length.out=steps))
    
    GridB <- expand.grid(seq(Hot0[1, 2], Hot0[1, 2] + stepwidth, length.out=steps),
                         seq(Hot0[2, 1], Hot0[2, 2] + stepwidth, length.out=steps))
    
    GridC <- expand.grid(seq(Hot0[1, 1] - stepwidth, Hot0[1, 2], length.out=steps),
                         seq(Hot0[2, 1] - stepwidth, Hot0[2, 1], length.out=steps))
    
    GridD <- expand.grid(seq(Hot0[1, 1], Hot0[1, 2] + stepwidth, length.out=steps),
                         seq(Hot0[2, 2], Hot0[2, 2] + stepwidth, length.out=steps))
    
    Grid <- rbind(GridA, GridB, GridC, GridD)
    
    FindcrHot <- apply(Grid, 1, function(x){
      theta <- matrix(x, 2)
      (n * t(est - theta) %*% solve(cov) %*% (est - theta)) <
        qf(p=alpha, df1=2, df2=df - p + 1) * p * df / (df - p + 1)
    })
    
    CrHot <- cbind(Grid, FindcrHot)[FindcrHot==1, ]
    
    Hot <- t(apply(CrHot[, -(p + 1)], 2, range))
    
  }else{
    
    Hot <- Hot0
    
  }
  
  HotOut <- cbind(est, Hot)
  rownames(HotOut) <- colnames(dat)
  colnames(HotOut) <- c("estimate", "lower", "upper")
  
  return(HotOut)
  
}
