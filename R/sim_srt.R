sim_srt <- function(wt, pars, parallel=F, n_trees,rec_method=1){    # simulate set of reconstructed trees
  if(length(pars)==2){
    pars[3] = pars[2]
    pars[2] = mu
  }
  if(parallel){
    no_cores <- detectCores()- 1
    cl <- makeCluster(no_cores)
    registerDoParallel(cl)
    trees <- foreach(i = 1:n_trees, combine = list) %dopar% dmea::rec_tree(wt=wt, pars=pars,rec_method = rec_method)
    stopCluster(cl)
  }
  else{
    trees = vector('list',length=n_trees)
    for (i in 1:n_trees){
      rec = dmea::rec_tree(wt=wt, pars=pars,rec_method = rec_method)
      trees[[i]] = rec
    }
  }
  # w=0
  # for(i in 1:n_trees){
  #   w[i]=trees[[i]]$weight
  # }
  # for(i in 1:n_trees){
  #   trees[[i]]$weight = w[i]/max(w)
  # }
  return(trees)
}
