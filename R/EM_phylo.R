EM_phylo <- function(wt, init_par, n_trees=10, n_it=30, printpar=TRUE, mu=0.1, impsam=FALSE, tol=0.0001, parallel=F){
  pars = init_par
  Pars = matrix(nrow=n_it,ncol=3)
  El = vector(mode='numeric',length = n_it)
  x1=c(0,2,1)
  x2=c(20,2,1)
  i=1
  q = c(0,0,0)
  ltt = NULL
  i=1
  while(i <= n_it){ #dist(rbind(x1, x2)) > tol &
    print(paste('iteration #',i,':'))
    if(printpar) print(pars)
    #,'are parameters on iteration',i,'and took',q[3],'segs'))
    x1 = pars
    p = proc.time()
    trees <- sim_srt(wt=wt, pars=pars, n_trees=n_trees)
    pars = subplex(par = pars, fn = llik_st, setoftrees = trees, impsam = impsam)$par
    x2 = pars
    Pars[i,] = pars
    qt=proc.time()-p
    q[i] = qt[3]
    i=i+1
  }
  return(list(pars=Pars))
}





