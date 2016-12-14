EM_phylo <- function(wt, init_par, n_trees=100, n_it=30, printpar=TRUE, mu=NULL, impsam=FALSE, dummy=0){
  pars = init_par
  for(j in 1:n_it){
    if(printpar) print(pars)
    trees <- sim_srt(wt=wt, pars=pars, parallel=F, n_trees=n_trees)
    if(length(init_par)==3){
      pars = subplex(par = c(8,0.175,0.9),fn = llik_st , setoftrees = trees, impsam = impsam)$par
    }
    else{
      pars = subplex(par = c(8,0.175),mu2par=mu,fn = llik_st , setoftrees = trees, impsam = impsam)$par
    }
    pars[3]=pars[3]+dummy
  }
  return(pars)
}

