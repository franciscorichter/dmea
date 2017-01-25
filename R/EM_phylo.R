EM_phylo <- function(wt, init_par, n_trees=100, n_it=30, printpar=TRUE, mu=NULL, impsam=FALSE, dummy=0){
  pars = init_par
  Pars = matrix(nrow=n_it,ncol=3)
  El = vector(mode='numeric',length = n_it)
  for(j in 1:n_it){
    #print(paste('iteration #',i))
    if(printpar) print(pars)
    trees <- sim_srt(wt=wt, pars=pars, parallel=F, n_trees=n_trees)
    if(length(init_par)==3){
      pars = subplex(par = pars, fn = llik_st, setoftrees = trees, impsam = impsam)$par
      p=subplex(par = pars, fn = llik_st, setoftrees = trees, impsam = impsam)$value
      if(printpar) print(paste('conditional E[log(f)]:',p))
      Pars[j,] = pars
      El[j] = p
    }
    else{
      pars = subplex(par = c(8,0.175),mu2par=mu,fn = llik_st , setoftrees = trees, impsam = impsam)$par
    }
    pars[3]=pars[3]+dummy
  }
  return(list(pars=pars,El=El,Pars=Pars))
}

