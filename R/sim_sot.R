sim_sot <- function(pars,n_trees,ct,parallel=F,useallcores=F,fossil=F,...){
  lambda0=pars[1]
  mu0=pars[2]
  K=pars[3]
  if(parallel){
    no_cores <- detectCores()- 1
    if(useallcores) no_cores <- detectCores()
    cl <- makeCluster(no_cores)
    registerDoParallel(cl)
    if(fossil){
      trees <- foreach(i = 1:n_trees, combine = list, .packages='dmea') %dopar% dmea::sim_phyl(ct=ct, lambda0=lambda0, mu0=mu0, K=K, seed=i)$phylo
    }
    else{
      trees <- foreach(i = 1:n_trees, combine = list, .packages='dmea') %dopar% dmea::sim_phyl(ct=ct, lambda0=lambda0, mu0=mu0, K=K, seed=i)$phylo.extant
    }
    stopCluster(cl)
  }
  else{
    trees = vector('list',length=n_trees)
    for (i in 1:n_trees){
      sim = sim_phyl(ct=ct, lambda0=lambda0, mu0=mu0, K=K, seed=i)
      if(fossil){sim=sim$phylo.extant}
      else{sim=sim$phylo}
      trees[[i]] = sim
    }
  }
  return(trees)
}
