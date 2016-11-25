EM_phylo <- function(bt,init_par,n_trees=100,n_it=10,printpar=TRUE,mu=NULL){
  par = init_par
  for(j in 1:n_it){
    if(printpar) print(par)
    trees = vector('list',length=n_trees)
    for (i in 1:n_trees){
      rec = rec_tree(w=bt,pars=par)
      trees[[i]] = rec
      #setofweith[i] = -rec$prob
    }
    if(length(init_par)==3){
      par = subplex(par = c(8,0.175,0.9),fn = llik_st , setoftrees = trees, impsam = FALSE)$par
    }
    else{
      par = subplex(par = c(8,0.175),mu2par=mu,fn = llik_st , setoftrees = trees, impsam = FALSE)$par
    }
  }
  return(par)
}

