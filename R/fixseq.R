fixseq <- function(dat, alpha=0.1){
  
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
  
  Fix <- rbind(T_1st, T_2nd)
  
  FixOut <- cbind(est, Fix)
  rownames(FixOut) <- colnames(dat)
  colnames(FixOut) <- c("estimate", "lower", "upper")
  
  return(FixOut)
  
}
