fixseq2D <- function(dat, alpha=0.1, equi=1.25, plotrange=c(0.77, 1.3), axnames=NULL, main="Fixed Sequence"){
  
  if(ncol(dat)!=2){
    stop("Data must be bivariate.")
  }
  
  n <- nrow(dat)
  p <- ncol(dat)
  df <- n - 1
  
  est <- matrix(colMeans(dat), p)
  cov <- cov(dat)
  
  t_1st <- c(min(0, est[1] - sqrt(cov[1, 1]) * qt(1 - alpha, df) / sqrt(n)),
             max(0, est[1] + sqrt(cov[1, 1]) * qt(1 - alpha, df) / sqrt(n)))
  
  t_2nd <- c(min(0, est[2] - sqrt(cov[2, 2]) * qt(1 - alpha, df) / sqrt(n)),
             max(0, est[2] + sqrt(cov[2, 2]) * qt(1 - alpha, df) / sqrt(n)))
  
  if(t_1st[1] > log(0.8) & t_1st[2] < log(1.25)){
    
    if(t_2nd[1] > log(0.8) & t_2nd[2] < log(1.25)){
      
      # both bioequivalent
      
      ma <- max(c(max(abs(t_1st)), max(abs(t_2nd))))
      T_1st <- T_2nd <- c(-ma, ma)
      
    }else{
      
      # only 1st bioequivalent
      
      T_1st <- c(log(0.8), log(1.25))
      T_2nd <- c(min(log(0.8), t_2nd[1]), max(log(1.25), t_2nd[2]))
      
    }
    
  }else{
    
    # 1st not bioequivalent (thus 2nd not tested)
    
    T_1st <- c(min(log(0.8), t_1st[1]), max(log(1.25), t_1st[2]))
    T_2nd <- c(NA, NA)
    
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
  segments(x0=T_1st[1], x1=T_1st[2], y0=est[2], y1=est[2], lwd=2)
  segments(y0=T_2nd[1], y1=T_2nd[2], x0=est[1], x1=est[1], lwd=2)
  par(mar=c(5, 4, 4, 2))
  
}