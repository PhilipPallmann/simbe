standardcorr <- function(dat, alpha=0.1, steps=100){
  
  n <- nrow(dat)
  p <- ncol(dat)
  df <- n - 1
  
  est <- matrix(colMeans(dat), p)
  LambdaInv <- solve(cov(dat))
  poolvar <- var(as.vector(as.matrix(dat)))
  #s2 <- poolvar / n
  
  togrid <- list()
  
  for(i in 1:p){
    togrid[[i]] <- seq(est[i] - 2 * poolvar, est[i] + 2 * poolvar, length.out=steps)
  }
  
  grid <- expand.grid(togrid)
  
  findcrUsu <- apply(grid, 1, function(x){
    theta <- x
    #sqrt(sum((est - theta)^2))^2 < (s2 * p * qf(p=1 - alpha, df1=p, df2=n - 1))
    t(est - theta) %*% LambdaInv %*% (est - theta) < (p / (n - 1) * qf(p=1 - alpha, df1=p, df2=n - 1))
  })
  
  crUsu <- cbind(grid, findcrUsu)[findcrUsu==1, ]
  
  Usu0 <- t(apply(crUsu[, -(p + 1)], 2, range, na.rm=TRUE))
  
  if(min(abs(Usu0[, 1] - est + 2 * poolvar)) < 0.001 | min(abs(Usu0[, 2] - est - 2 * poolvar)) < 0.001){
    
    togrid2 <- list()
    
    for(i in 1:p){
      togrid2[[i]] <- seq(est[i] - 8 * poolvar, est[i] + 8 * poolvar, length.out=4 * steps)
    }
    
    grid2 <- expand.grid(togrid2)
    
    findcrUsu2 <- apply(grid2, 1, function(x){
      theta <- x
      #sqrt(sum((est - theta)^2))^2 < (s2 * p * qf(p=1 - alpha, df1=p, df2=n - 1))
      t(est - theta) %*% LambdaInv %*% (est - theta) < (p / (n - 1) * qf(p=1 - alpha, df1=p, df2=n - 1))
    })
    
    crUsu2 <- cbind(grid2, findcrUsu2)[findcrUsu2==1, ]
    
    Usu0 <- t(apply(crUsu2[, -(p + 1)], 2, range))
    
  }
  
  stepwidth <- 4 * poolvar / steps
  
  if(p==2){
    
    GridA <- expand.grid(seq(Usu0[1, 1] - stepwidth, Usu0[1, 1], length.out=steps),
                         seq(Usu0[2, 1] - stepwidth, Usu0[2, 2], length.out=steps))
    
    GridB <- expand.grid(seq(Usu0[1, 2], Usu0[1, 2] + stepwidth, length.out=steps),
                         seq(Usu0[2, 1], Usu0[2, 2] + stepwidth, length.out=steps))
    
    GridC <- expand.grid(seq(Usu0[1, 1] - stepwidth, Usu0[1, 2], length.out=steps),
                         seq(Usu0[2, 1] - stepwidth, Usu0[2, 1], length.out=steps))
    
    GridD <- expand.grid(seq(Usu0[1, 1], Usu0[1, 2] + stepwidth, length.out=steps),
                         seq(Usu0[2, 2], Usu0[2, 2] + stepwidth, length.out=steps))
    
    Grid <- rbind(GridA, GridB, GridC, GridD)
    
    FindcrUsu <- apply(Grid, 1, function(x){
      theta <- x
      #sqrt(sum((est - theta)^2))^2 < (s2 * p * qf(p=1 - alpha, df1=p, df2=n - 1))
      t(est - theta) %*% LambdaInv %*% (est - theta) < (p / (n - 1) * qf(p=1 - alpha, df1=p, df2=n - 1))
    })
    
    CrUsu <- cbind(Grid, FindcrUsu)[FindcrUsu==1, ]
    
    Usu <- t(apply(CrUsu[, -(p + 1)], 2, range))
    
  }else{
    
    Usu <- Usu0
    
  }
  
  UsuOut <- cbind(est, Usu)
  rownames(UsuOut) <- colnames(dat)
  colnames(UsuOut) <- c("estimate", "lower", "upper")
  
  return(UsuOut)
  
}
