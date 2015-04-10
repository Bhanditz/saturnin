edge.prob <-
function(W, log = TRUE, account.prior = FALSE, q0 = 0.5, nbcores = 1){
  p <- nrow(W)
  if (log){
    M <- exp(W - mean(W[upper.tri(W)]))
    diag(M) <- 0
  } else {
    M <- W
  }
  Delta <- -M+diag(apply(M,1,sum))
  
  if (requireNamespace("parallel",quietly = TRUE) &&
      (nbcores > 1)){
    RES <- parallel::mcmapply(function(i){
      Delta_i <- Delta[-i,-i]
      D2 <- (diag(chol2inv(chol(Delta_i))))
      append(D2,0,after=i-1)
    },
    1:p,
    mc.cores = nbcores)
  } else {
    RES <- mapply(function(i){
      Delta_i <- Delta[-i,-i]
      D2 <- (diag(chol2inv(chol(Delta_i))))
      append(D2,0,after=i-1)
    },1:p)
  }
  
  prob <- abs(RES*M)
  if (account.prior){return(account.for.prior(prob,q0))} else {return(prob)}
}
