mle_tree <- function(tree,init_par = c(2,1,60)){
  p = subplex(par = init_par, fn = llik, n = tree$n, E = tree$E, t = tree$wt)$par
  return(p)
}
