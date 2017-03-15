expectedLTT2 <- function(pars, ct=15, n_it= 10){ # it should also include model
  # add message that checks that do not have extinct species
  a = NULL # change this for a list with n_it length
  # diversity dependence model. parameters are par=c(lambda,mu,K)
  for (i in 1:(n_it)){
    s = sim_phyl(ct = ct, lambda0 = pars[1], mu0 = pars[2], K = pars[3])
    s = phylo2p(drop.fossil(s$newick))
    a[[i]] = s$wt
  }
  n.obs <- sapply(a, length)
  m.obs = round(mean(n.obs))
  seq.max <- seq_len(max(n.obs))
  mat <- t(sapply(a, "[", i = seq.max))
  mat = mat[,1:m.obs]
  # Warning: the median represents better the 'average' of waiting time because the distribution is skewded
  expect = colMeans(mat,na.rm = TRUE)
  expect = data.frame(t=expect, Ex = 1:length(expect))
  return(expect)
}



