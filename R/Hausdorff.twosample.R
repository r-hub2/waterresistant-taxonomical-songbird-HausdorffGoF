
H_test_2s_1d = function(x1, x2, nboots = 2000, Exact = FALSE){
  test_stat = H_stat_2s_1d_tr(x1, x2) #Finds test stat
  comb = c(x1, x2)
  n = length(comb)
  na = min(length(x1),length(x2))
  if(choose(n,na)>2.147483646e9){
    Exact = F
    message("Sample sizes too large, switch to Monte Carlo")
  }
  if(!Exact){
    Method = "Monte Carlo"
    nboots = as.integer(nboots)					#Speeds up comparison below.
    reps = 0L
    bigger = 0L							  #Initializes Counter
    while (reps < nboots) {						#Loops over vector
      vec_labels = rep(F,n)
      vec_labels[sample.int(n,na,F)] = T #Samples indexes
      boot_t = H_stat_2s_1d_tr(comb[vec_labels],comb[!vec_labels]) #boot strap test stat
      bigger = bigger +(boot_t >= test_stat) #if new stat is bigger, increment
      reps = 1L+reps
    }
    out = bigger/nboots
    out[which(out==0)] = 1/(2*nboots)
    
  } else {
    Method = "Exact"
    bigger = 0L
    combn_exact = combn(n, na)
    for(i in 1:choose(n,na)){
      bigger = bigger + (H_stat_2s_1d_tr(comb[combn_exact[,i]], comb[-combn_exact[,i]]) >= test_stat)
    }
    out = bigger/choose(n,na)
    
  }
  names(test_stat) <- "H"
  result = list(p.value = out, method = paste("Two-sample Hausdorff Test (",Method, ")", sep = ""),
                statistic = test_stat, alternative = "two-sided",
                data.name = paste(deparse(substitute(x1)), "and", deparse(substitute(x2))))
  class(result) = "htest"
  return(result)
}


H_test_2s_2d <- function(x, y, nboots = 2000, Exact = FALSE,
                         invariant = FALSE, tol = 1e-6) {
  if (!is.matrix(x) || ncol(x) != 2)
    stop("`x` must be a two-column numeric matrix.")
  if (!is.matrix(y) || ncol(y) != 2)
    stop("`y` must be a two-column numeric matrix.")
  
  m    <- nrow(x)
  n    <- nrow(y)
  pool <- rbind(x, y)
  
  # Switch to Monte Carlo if exact enumeration would overflow
  if (Exact && choose(m + n, m) > 2.147483647e9) {
    Exact <- FALSE
    message("Sample sizes too large for exact enumeration; switching to Monte Carlo.")
  }
  
  # Statistic function: standard or order-invariant
  stat_fn <- if (!invariant) {
    function(a, b) H_stat_2s_2d(a, b, tol = tol)
  } else {
    function(a, b) {
      max(
        H_stat_2s_2d( a,                          b,                         tol = tol),
        H_stat_2s_2d(cbind( a[,1], -a[,2]), cbind( b[,1], -b[,2]),           tol = tol),
        H_stat_2s_2d(cbind(-a[,1],  a[,2]), cbind(-b[,1],  b[,2]),           tol = tol),
        H_stat_2s_2d(cbind(-a[,1], -a[,2]), cbind(-b[,1], -b[,2]),           tol = tol)
      )
    }
  }
  
  observed <- stat_fn(x, y)
  
  if (!Exact) {
    method_str <- "Two-sample bivariate Hausdorff test (Monte Carlo permutation)"
    perm_stats <- replicate(nboots, {
      idx <- sample.int(m + n, m, replace = FALSE)
      stat_fn(pool[idx, , drop = FALSE], pool[-idx, , drop = FALSE])
    })
    p_value <- mean(perm_stats >= observed)
    if (p_value == 0) p_value <- 1 / (2 * nboots)
  } else {
    method_str <- "Two-sample bivariate Hausdorff test (Exact)"
    combn_idx  <- combn(m + n, m)
    bigger     <- 0L
    for (i in seq_len(ncol(combn_idx))) {
      idx    <- combn_idx[, i]
      bigger <- bigger + (stat_fn(pool[idx, , drop = FALSE],
                                  pool[-idx, , drop = FALSE]) >= observed)
    }
    p_value <- bigger / choose(m + n, m)
  }
  
  if (invariant) method_str <- paste0(method_str, ", order-invariant")
  
  result <- list(
    statistic   = c(H = observed),
    p.value     = p_value,
    method      = method_str,
    alternative = "two-sided",
    data.name   = paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  )
  class(result) <- "htest"
  result
}

