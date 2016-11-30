sim_est <- function(n_trees, init_par=c(8,0.175,0.9)){ # simulate a tree, drop fossil, and estimate parameters back after bootstrap reconstruction
 # seed = round(runif(1,0,10000000))  # is this silly?
  st = dmea::sim_phyl(seed=seed)
  p <- subplex(par = init_par, fn = llik,n = st$n, E = st$E, t = st$t)$par
  sit = dmea::drop.fossil(st$newick)
  sit = dmea::phylo2p(sit)
  trees = sim_srt(tree=sit, pars=p, parallel = F, n_trees = n_trees)
  pars = subplex(par = init_par, fn = llik_st , setoftrees = trees, impsam = FALSE)$par
  return(data.frame(real=p,est=pars))
}


