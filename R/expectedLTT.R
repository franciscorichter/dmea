expectedLTT <- function(pars, dt = 0.1, ct=15, n_it= 100){ # it should also include model
  grid = seq(0,ct, by=dt)
  a = data.frame(t=grid)
  # diversity dependence model. parameters are par=c(lambda,mu,K)
  s = sim_phyl(ct = ct, lambda0 = pars[1], mu0 = pars[2], K = pars[3])
  a[[2]] = approx(cumsum(s$wt),s$n,xou=grid)$y
  n_it
  for (i in 1:(n_it)){
    s = sim_phyl(ct = ct, lambda0 = pars[1], mu0 = pars[2], K = pars[3])
    a[[i+2]] = approx(cumsum(s$wt), s$n,xou=grid)$y
  }
  a[is.na(a)]=0
  expect = data.frame(t=a[[1]], E = rowMeans(a[,2:100]))
  return(expect)
}
