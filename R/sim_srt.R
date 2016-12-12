sim_srt <- function(wt, pars, parallel=F, n_trees){    # simulate set of reconstructed trees
  if(parallel){
    # todo: check if registerdoparallel is on. If not, warning or break.
    trees <- foreach(i = 1:n_trees, combine = list) %dopar% dmea::rec_tree(w=tree$t, pars=pars)
  }
  else{
    trees = vector('list',length=n_trees)
    for (i in 1:n_trees){
      rec = dmea::rec_tree(wt=wt, pars=pars)
      trees[[i]] = rec
    }
  }
  return(trees)
}
