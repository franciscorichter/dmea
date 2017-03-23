expectedLTT2 <- function(pars, ct=15, n_it= 10, median=FALSE){ # it should also include model
  # add message that checks that do not have extinct species
  a = vector("list", n_it)
  # diversity dependence model. parameters are par=c(lambda,mu,K)
  time=proc.time()
  for (i in 1:(n_it)){
    s = sim_phyl(seed=i, ct = ct, lambda0 = pars[1], mu0 = pars[2], K = pars[3])$newick.extant.p
    a[[i]] = s$wt
  }
  proc.time()-time
  n.obs <- sapply(a, length)
  #m.obs = round(mean(n.obs))
  seq.max <- seq_len(max(n.obs))
  mat <- t(sapply(a, "[", i = seq.max))
  #mat = mat[,1:m.obs]
  # Warning: the median represents better the 'average' of waiting time because the distribution is skewded
  expect = colMeans(mat,na.rm = TRUE)
  ce = cumsum(expect)
  #ce
  expect = list(bt=c(cumsum(expect[ce<15]),15), Ex = 1:length(cumsum(expect[ce<15])))
  return(expect)
}



