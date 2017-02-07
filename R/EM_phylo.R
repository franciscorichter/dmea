EM_phylo <- function(wt, init_par, n_trees=100, n_it=30, printpar=TRUE, mu=0.1, impsam=FALSE, dummy=0,tol=0.001, parallel=F){
  pars = init_par
  Pars = matrix(nrow=n_it,ncol=3)
  El = vector(mode='numeric',length = n_it)
  x1=c(0,2,1)
  x2=c(20,2,1)
  i=1
  while(dist(rbind(x1, x2))>tol){
    #print(paste('iteration #',i))
    if(printpar) print(pars)
    x1 = pars
    p = proc.time()
    trees <- sim_srt(wt=wt, pars=pars, parallel=parallel, n_trees=n_trees)
    if(length(init_par)==3){
      pars = subplex(par = pars, fn = llik_st, setoftrees = trees, impsam = impsam)$par
      x2 = pars
      #Pars[j,] = pars
      #El[j] = p
      q=proc.time()-p
      print(paste('iteration',i,'took',q[3],'seg'))
      i=i+1
    }
    else{
      pars = c(subplex(par = c(8,0.175),fn = llik_st ,mu=mu , setoftrees = trees, impsam = impsam)$par,mu)
      x2 = pars
      #print('Lenght of parameters is not 3!')
    }
    pars[3]=pars[3]+dummy
  }
  return(list(pars=pars))
}





