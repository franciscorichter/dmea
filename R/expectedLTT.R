expectedLTT <- function(pars, dt = 0.1, ct=15, n_it= 100, value = 0, grid = 0, drop.extinct = FALSE){ # it should also include model
  if(length(grid) == 1 & grid[1] == 0) grid = seq(0,ct, by=dt)
  if(value>0){
    grid = value
  }
  a = data.frame(t=grid)
  # diversity dependence model. parameters are par=c(lambda,mu,K)
  s = sim_phyl(ct = ct, lambda0 = pars[1], mu0 = pars[2], K = pars[3])
  if(drop.extinct == TRUE) s = phylo2p(drop.fossil(s$newick))
  a[[2]] = approx(cumsum(s$wt),s$n,xou=grid,  rule = 2)$y
  for (i in 1:(n_it)){
    s = sim_phyl(ct = ct, lambda0 = pars[1], mu0 = pars[2], K = pars[3])
    if(drop.extinct == TRUE) s = phylo2p(drop.fossil(s$newick))
    a[[i+2]] = approx(cumsum(s$wt), s$n,xou=grid, rule = 2)$y
  }
  #a[is.na(a)]=0
  expect = data.frame(t=a[[1]], Ex = rowMeans(a[,2:100]))
#  expect_integers
  return(expect)
}
