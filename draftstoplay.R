st <- sim_phyl()
st2 = drop.fossil(st$newick)
st2 = phylo2p(st2)
for(i in 1:100){
  rec_tree(wt = st2$t)
}




p = proc.time()
no_cores <- detectCores()- 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)

n_sim = 1000
n_trees = 10
MP = matrix(nrow=n_sim,ncol=3)
RP = matrix(nrow=n_sim,ncol=3)

ests <- foreach(i = 1:n_sim, .combine=data.frame,.packages='dmea') %dopar% {
  set.seed(i)
  sim_est(n_trees=n_trees,impsam=FALSE)
}
for (i in 1:n_sim){
  RP[i,] = ests[,(2*i-1)]
  MP[i,] = ests[,2*i]
}
par_est_vis(P=MP,par=1,PR=RP)
par_est_vis(P=MP,par=2,PR=RP)
par_est_vis(P=MP,par=3,PR=RP)




p = proc.time()
no_cores <- detectCores()- 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)

n_sim = 1000
n_trees = 10
MP = matrix(nrow=n_sim,ncol=3)
RP = matrix(nrow=n_sim,ncol=3)

ests <- foreach(i = 1:n_sim, .combine=data.frame,.packages='dmea') %dopar% {
  set.seed(i)
  sim_est(n_trees=n_trees,impsam=TRUE)
}
for (i in 1:n_sim){
  RP[i,] = ests[,(2*i-1)]
  MP[i,] = ests[,2*i]
}
par_est_vis(P=MP,par=1,PR=RP)
par_est_vis(P=MP,par=2,PR=RP)
par_est_vis(P=MP,par=3,PR=RP)

stopCluster(cl)

Max=0
for (i in 1:10){
  ttt = trees[[i]]
  www[i] = ttt$weight
  Max[i] = max(ttt$n)
  print(paste('max of n',max(ttt$n)))
  print(paste('last n', tail(ttt$n,n=1)))
  print(paste('speciations',sum(E)))
  print(paste('extinctions',sum(1-E)))
}

