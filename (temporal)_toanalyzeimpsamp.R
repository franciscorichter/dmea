st <- sim_phyl(ct=15)
st2 <- drop.fossil(st$newick)
st3 <- phylo2p(st2)
W1 = NULL
W2 = NULL
for(k in 1:100){
wt = st3$t
pars = c(0.8,0.0175,0.1)

lambda0 = pars[1]
mu0 = pars[3]
K = (lambda0-mu0)/pars[2]
n = 1:length(wt)
i = 1
E = rep(1,(length(wt)-1))
fake = FALSE
ct = sum(wt)
prob = 0
while(i < length(wt)){
  N = n[i]
  if(model == "dd"){  # diversity-dependence model
    lambda = max(0,lambda0 - (lambda0-mu0)*N/K)
    mu = mu0
    lambda = rep(lambda,N)
    mu = rep(mu,N)
  }
  if(model == 'cr'){ # constant-rate model
    lambda = rep(lambda0,N)
    mu = rep(mu0,N)
  }
  s = sum(lambda)
  if(s==0){
    #print('s=0')
    break
  }
  if(fake){ # in case there was an speciation but not extinction previously
    cwt = cwt - t_spe
    cbt = cbt + t_spe
  }
  else{
    cwt = wt[i]
    cbt = cumsum(wt)[i] - cwt
  }
  t_spe = rexp(1,s)
  if (t_spe < cwt){
    t_ext = rexp(1,mu0) # this is not as general as trees with trait-dependance species yet,
    t_ext = cbt + t_spe + t_ext
    prob = prob + log(1-exp(-s*cwt))
    if (t_ext < ct){
      up = update_tree(wt=wt,t_spe = (cbt + t_spe), t_ext = t_ext, E = E, n = n)
      E = up$E
      n = up$n
      wt = up$wt
      fake = FALSE
      prob = prob + log(1-exp(-mu0*(ct-t_spe-cbt)))
    }else{
      prob = prob - mu0*(ct-t_spe-cbt)
      fake = TRUE
      i = i-1
    }
  }else{
    fake = FALSE
    prob = prob - s*cwt
  }
  i = i+1
}
#  newick = p2phylo(wt,E,)
L = create_L(wt,E)
n_prob = num_weigh(rec_tree=list(wt=wt,E=E,n=n,L=L), pars_rec = pars, ct=ct)
W1[k] = n_prob
W2[k] = prob
}


st4 <- rec_tree(wt = st3$t, pars = c(0.8,0.0175,0.1))
### implementar importance sampling...



n_sim = 1000
n_trees = 100
MP2 = matrix(nrow=1000,ncol=3)
#MPf = matrix(nrow=n_sim,ncol=3)
RP2 = matrix(nrow=1000,ncol=3)
setofweith = NULL
for(j in 1:1000){
  print(j)
  st = sim_phyl(seed=j)
  p <- subplex(par = c(8,0.175,0.9),fn = llik,n = st$n, E = st$E, t = st$t)$par
  RP2[j,] = p
  sit = drop.fossil(st$newick)
  sit = phylo2p(sit)
  #P = matrix(nrow=n_trees,ncol=3)
  trees = vector('list',length=n_trees)
  for (i in 1:n_trees){
    rec = rec_tree(w=sit$t, pars=p)
    setofweith[i] = -rec$prob
    #p2 <- subplex(par = c(8,0.175,0.9),fn = llik,n = rec$n, E = rec$E, t = rec$wt)$par
    trees[[i]] = rec
    #P[i,] = p2
  }
  #####
  #MPf[j,] = colMeans(P)
  pars = subplex(par = c(8,0.175,0.9),fn = llik_st , setoftrees = trees, impsam = FALSE)$par
  MP2[j,] = pars
}

