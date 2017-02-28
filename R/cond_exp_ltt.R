cond_exp_ltt <- function(obsPhylo,pars=0,n_it=100,est_tree=NULL){
  # get branching times of obsPhylo
  p = phylo2p(obsPhylo)
  ct = sum(p$wt)
  # get expected ltt with pars
  grid = cumsum(p$wt)
  if (length(pars)>1){
    n_exp = expectedLTT(pars=pars, grid=grid,n_it = n_it)}
  else{
 #   if(est_tree == NULL) stop('not estimated tree neither parameters provided')
    n_exp = est_tree
  }
  # aproax expected ltt to branchingtimes
#  n_exp = approx(exp$t, exp$Ex, xou=grid, rule = 2)$y
  #tt = approx(exp$Ex, exp$t, xou=1:45)$y
  n_exp = round(n_exp[,2]) # it may be better without round  (?)
  #Ex = expectedLTT(pars = pars, ct = ct)$E  # approximate is too expensive :( ...
  #added = 0
  wt = p$wt
  E = p$E
  obs_n = p$n
  n = p$n
  N = 1
  cbt=0
  pn = p
  for(i in 1:length(wt)){
    cwt = p$wt[i]
    if(i>1) cbt = grid[i-1]
    next_n = obs_n[i] + 1
    miss = n_exp[i] - next_n
    N = obs_n[i]
    if(miss>0){
      for(j in 1:miss){
        ## ADD ALL THOSE PRINT ON A LOG FILE
        #print(paste('iteration',j))
        #print(paste('current waiting time',cwt))
        lambda = pars[1]-((pars[1]-pars[2])/pars[3])*N
        suml = N*lambda
        t_spe = rtexp(1,suml,cwt)  # there is a more proper way to do this (?)
        t_ext = rtexp(1, pars[2], ct-(cbt + t_spe))
        #print(paste('new speciation',cbt + t_spe))
        #print(paste('new extinction',(cbt + t_spe)+t_ext))
        up = update_tree(pn, t_spe = (cbt + t_spe), t_ext = (cbt + t_spe)+t_ext)
        wt = up$wt
        n = up$n
        E = up$E
        pn = list(wt=wt,E=E,n=n)
        obs_n[grid > (cbt + t_spe) & grid < (cbt + t_ext)] = obs_n[grid > (cbt + t_spe) & grid < (cbt + t_ext)] + 1
        cwt = cwt - t_spe
        cbt = cbt + t_spe
        N = N+1 # this is not exact
      }
    }
  }
  return(list(wt=wt,n=n,E=E))
}

