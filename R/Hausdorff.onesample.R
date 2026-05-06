H_test_1s_1d = function(x, CDF, CDFinverse = NULL, pdf = NULL, tol = 1e-10,  max.init = 1000){
  q = H_stat_1s_1d(x, CDF = CDF, pdf = pdf, tol = tol, max.init = max.init)
  n = length(x)
  PVAL = H_test_c_cdf(q, n, CDF=CDF, CDFinverse = CDFinverse, pdf = NULL, tol = tol, max.init = max.init)
  names(q) <- "H"
  result = list(p.value = PVAL, method = "One-sample Hausdorff test",
                statistic = q, alternative = "two-sided",
                data.name = paste(deparse(substitute(x))))
  class(result) = "htest"
  return(result)
}

H_test_c_cdf = function(q, n, CDF, CDFinverse = NULL, pdf = NULL,
                        tol = 1e-10, max.init = 1000) {
  
  if ((q >= 1) || (q <= 0)) return(ifelse(q >= 1, 1, 0))
  
  f_b = rep(0, n)
  f_a = rep(0, n)
  
  if (is.null(CDFinverse)) {
    f <- function(lambda) {
      function(x) CDF(x) - lambda
    }
    
    if (is.null(pdf)) {
      for (i in 1:n) {
        if (i/n - q < 0) {
          f_a[i] <- 0
        } else {
          fta    <- f(i/n - q)
          f_a[i] <- CDF(uniroot(fta, c(0, 1), extendInt = "yes",
                                tol = tol)$root - q)
        }
        
        if ((i-1)/n + q > 1) {
          f_b[i] <- 1
        } else {
          ftb    <- f((i-1)/n + q)
          f_b[i] <- CDF(uniroot(ftb, c(0, 1), extendInt = "yes",
                                tol = tol)$root + q)
        }
      }
    } else {
      for (i in 1:n) {
        num.init = 0
        if (i/n - q < 0) {
          f_a[i] <- 0
        } else {
          fta  <- f(i/n - q)
          temp <- 0
          h    <- fta(temp) / pdf(temp)
          while (abs(h) > tol) {
            h    <- fta(temp) / pdf(temp)
            temp <- temp - h
            num.init <- num.init + 1
            if (max.init < num.init) stop("Maximum number of iteration reached.")
          }
          f_a[i] <- CDF(temp - q)
        }
        
        num.init = 0
        if ((i-1)/n + q > 1) {
          f_b[i] <- 1
        } else {
          ftb  <- f((i-1)/n + q)
          temp <- 0
          h    <- ftb(temp) / pdf(temp)
          while (abs(h) > tol) {
            h    <- ftb(temp) / pdf(temp)
            temp <- temp - h
            num.init <- num.init + 1
            if (max.init < num.init) stop("Maximum number of iteration reached.")
          }
          f_b[i] <- CDF(temp + q)
        }
      }
    }
    
  } else {
    for (i in 1:n) {
      if (i/n - q < 0) {
        f_a[i] <- 0
      } else {
        f_a[i] <- CDF(CDFinverse(i/n - q) - q)
      }
      
      if ((i-1)/n + q > 1) {
        f_b[i] <- 1
      } else {
        f_b[i] <- CDF(CDFinverse((i-1)/n + q) + q)
      }
    }
  }
  
  f_a[f_a < 0] = 0
  f_b[f_b > 1] = 1
  df <- data.frame(rbind(f_b, f_a))
  write.table(df, "Boundary_Crossing_Time.txt",
              sep = ", ", row.names = FALSE, col.names = FALSE)
  PVAL <- KSgeneral::ks_c_cdf_Rcpp(n)
  file.remove("Boundary_Crossing_Time.txt")
  return(PVAL)
}

H_stat_1s_1d = function(x, CDF, pdf = NULL, tol = 1e-10, max.init = 1000){
  
  f <- function(lambda){
    function(u){
      u + CDF(u) - lambda
    }
  }
  x = sort(x)
  edfx = ecdf(x)
  y0 = edfx(x)
  yf = CDF(x)
  n = length(x)
  x = x[ceiling(1:(2*n)/2)]
  y0 = c(0,y0[ceiling(1:(2*n-1)/2)])
  yf = yf[ceiling(1:(2*n)/2)]
  Index = rep(c(T,F),n)
  Index = (Index == (yf>y0))
  x = x[Index]
  y0 = y0[Index]
  parameter = x + y0
  
  xs = rep(0,length(parameter))
  if(is.null(pdf)){
    for(i in 1:length(parameter)){
      ft = f(parameter[i])
      xs[i] = uniroot(ft, c(x[i]-1,x[i]+1), tol = tol)$root
    }
  } else {
    for(i in 1:length(parameter)){
      ft = f(parameter[i])
      h = ft(x[i])/(1+pdf(x[i]))
      xs[i] = x[i]
      num.init = 0
      while (abs(h)>tol) {
        h = ft(xs[i])/(1+pdf(xs[i]))
        xs[i] = xs[i] - h
        num.init = num.init + 1
        if(max.init < num.init) stop("Maximum number of iteration reached.")
      }
    }
  }
  H_vec = xs - x
  return(max(abs(H_vec)))
}