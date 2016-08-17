plotMV2D <- function(dat, method, alpha=0.1, axnames=c("Mean", "Variance"),
                     main="Title", col="black", steps=400, searchwidth=1){
  
  if(is.vector(dat)!=TRUE){
    stop("dat must be a vector of numeric values.")
  }
  
  method <- match.arg(method, choices=c("Mood", "large", "plugin", "pluginF", "LRT"))
  
  if(method=="Mood"){
    
    n <- length(dat)
    df <- n - 1
    mea <- mean(dat)
    s <- sqrt(var(dat) * df / n)
    
    togrid <- list()
    togrid[[1]] <- seq(mea - qnorm(1 - alpha/16) * searchwidth * s / sqrt(n),
                       mea + qnorm(1 - alpha/16) * searchwidth * s / sqrt(n), length.out=steps)
    togrid[[2]] <- seq(s^2 * 1/searchwidth * n / qchisq(df=df, 1 - alpha/16),
                       s^2 * searchwidth * n / qchisq(df=df, alpha/16), length.out=steps)
    
    grid <- expand.grid(togrid)
    
    grid[, 3] <- mea - qnorm(1 - alpha/2) * sqrt(grid[, 2]) / sqrt(n) < grid[, 1] &
      grid[, 1] < mea + qnorm(1 - alpha/2) * sqrt(grid[, 2]) / sqrt(n) &
      s^2 * n / qchisq(df=df, 1 - alpha/2) < grid[, 2] &
      grid[, 2] < s^2 * n / qchisq(df=df, alpha/2)
    
    crFinal <- grid[grid[, 3]==TRUE, ]
    
  }
  
  if(method=="large"){
    
    n <- length(dat)
    df <- n - 1
    mea <- mean(dat)
    s <- sqrt(var(dat) * df / n)
    
    togrid <- list()
    togrid[[1]] <- seq(mea - qnorm(1 - alpha/16) * s / sqrt(n),
                       mea + qnorm(1 - alpha/16) * s / sqrt(n), length.out=steps)
    togrid[[2]] <- seq(s^2 * n / qchisq(df=df, 1 - alpha/16), s^2 * n / qchisq(df=df, alpha/16), length.out=steps)
    
    grid <- expand.grid(togrid)
    
    grid[, 3] <- n/grid[, 2] * (mea - grid[, 1])^2 + n/(2 * grid[, 2]^2) * (s^2 - grid[, 2])^2 < qchisq(1 - alpha, df=2)
    
    crFinal <- grid[grid[, 3]==TRUE, ]
    
  }
    
  if(method=="plugin"){
    
    n <- length(dat)
    df <- n - 1
    mea <- mean(dat)
    s <- sqrt(var(dat) * df / n)
    
    togrid <- list()
    togrid[[1]] <- seq(mea - qnorm(1 - alpha/16) * s / sqrt(n),
                       mea + qnorm(1 - alpha/16) * s / sqrt(n), length.out=steps)
    togrid[[2]] <- seq(s^2 * n / qchisq(df=df, 1 - alpha/16), s^2 * n / qchisq(df=df, alpha/16), length.out=steps)
    
    grid <- expand.grid(togrid)
    
    grid[, 3] <- n/s^2 * (mea - grid[, 1])^2 + n/(2 * s^4) * (s^2 - grid[, 2])^2 < qchisq(1 - alpha, df=2)
    
    crFinal <- grid[grid[, 3]==TRUE, ]
    
  }
  
  if(method=="pluginF"){
    
    n <- length(dat)
    df <- n - 1
    mea <- mean(dat)
    s <- sqrt(var(dat) * df / n)
    
    togrid <- list()
    togrid[[1]] <- seq(mea - qnorm(1 - alpha/16) * s / sqrt(n),
                       mea + qnorm(1 - alpha/16) * s / sqrt(n), length.out=steps)
    togrid[[2]] <- seq(s^2 * n / qchisq(df=df, 1 - alpha/16), s^2 * n / qchisq(df=df, alpha/16), length.out=steps)
    
    grid <- expand.grid(togrid)
    
    grid[, 3] <- n/s^2 * (mea - grid[, 1])^2 + n/(2 * s^4) * (s^2 - grid[, 2])^2 < qf(1 - alpha, df1=2, df2=n - 2)
    
    crFinal <- grid[grid[, 3]==TRUE, ]
    
  }
  
  if(method=="LRT"){
    
    n <- length(dat)
    df <- n - 1
    mea <- mean(dat)
    s <- sqrt(var(dat) * df / n)
    
    togrid <- list()
    togrid[[1]] <- seq(mea - qnorm(1 - alpha/16) * s / sqrt(n),
                       mea + qnorm(1 - alpha/16) * s / sqrt(n), length.out=steps)
    togrid[[2]] <- seq(s^2 * n / qchisq(df=df, 1 - alpha/16), s^2 * n / qchisq(df=df, alpha/16), length.out=steps)
    
    grid <- expand.grid(togrid)
    
    grid[, 3] <- n * log(grid[, 2] / s^2) + n * s^2 / grid[, 2] + n * (mea - grid[, 1])^2 / grid[, 2] - n < qchisq(1 - alpha, df=2)
    
    crFinal <- grid[grid[, 3]==TRUE, ]
    
  }

  if(min(crFinal[, 1])==min(grid[, 1]) | max(crFinal[, 1])==max(grid[, 1]) |
       min(crFinal[, 2])==min(grid[, 2]) | max(crFinal[, 2])==max(grid[, 2])){
    warning("The search grid is too narrow, please increase searchwidth.")
  }
  
  xlims <- range(grid[, 1])
  ylims <- range(grid[, 2])
  
  par(mar=c(5, 5, 4, 2))
  plot(0, xlim=xlims, ylim=ylims, las=1, xlab=axnames[1], ylab=axnames[2],
       cex.main=2.5, cex.axis=1.5, cex.lab=1.7, main=main)
  polygon(crFinal[chull(crFinal[, -3]), -3], col=col, border=col)
  points(mea, s^2, pch=19, col="white")
  par(mar=c(5, 4, 4, 2))
  
}
