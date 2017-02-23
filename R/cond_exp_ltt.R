# WORK IN PROGRESS

cond_exp_ltt <- function(obsPhylo,pars){
  # get branching times of obsPhylo
  p = phylo2p(obsPhylo)
  # get expected ltt with pars
  exp = expectedLTT(pars=pars)
  # aproax expected ltt to branchingtimes
  grid = cumsum(p$wt)
  nn = approx(exp$t,exp$E,xou=grid)$y
  #added = 0
  N = 1
  j = 1
  wt=E=n=0
  for(i in 1:length(nn)){
    lambda = pars[1]-((pars[1]-pars[2])/pars[3])*N
    suml = N*lambda
    exp_spec_time = 1/suml
    mu = pars[2]
   # summ = (max(0,N-p$n[i]))*mu
    exp_ext_time = 1/mu
 #   if(exp_ext_time==0) exp_ext_time=999999  # this means that there is not possible extinction
    cwt = p$wt[i]
    while(nn[i] > N & exp_spec_time < cwt){
      wt[j] = exp_spec_time
      n[j] = N
      E[j] = 1
      j = j + 1
      cwt = cwt - exp_spec_time
      N = N+1
      lambda = pars[1]-((pars[1]-pars[2])/pars[3])*N
      suml = N*lambda
      et[k] =
#      exp_spec_time = 1/suml
#      summ = (max(0,N-p$n[i]))*mu
#      exp_ext_time = 1/summ
#      if(exp_ext_time==0) exp_ext_time=999999  # this means that there is not possible extinction
    }
    while(nn[i] < N & exp_ext_time < cwt){
      wt[j] = exp_ext_time
      n[j] = N
      E[j] = 0
      j = j = 1
      cwt = cwt - exp_ext_time
      N = N-1
      summ = (max(0,N-p$n[i]))*mu
      exp_ext_time = 1/summ
      if(exp_ext_time==0) exp_ext_time=999999  # this means that there is not possible extinction
    }
    wt[j] = cwt
    E[j] = 1
    n[j] = N
    N = N + 1
    j = j+1
  }
  newp = list(wt = wt, n = n, E= E)
  return(newp)
}
