st = sim_phyl(seed=7)
plot(st$newick)
st2 = drop.fossil(st$newick)
plot(st2)
st2 = phylo2p(st2)
# Recontruction and estimations
rec = rec_tree(wt = st2$t)
p <- subplex(par = c(8,0.175,0.9),fn = llik,n = st$n, E = st$E, t = st$t)
c(p$par[1],(p$par[1]-p$par[3])/p$par[2],p$par[3])
p2 <- subplex(par = c(8,0.175,0.9),fn = llik,n = rec$n, E = rec$E, t = rec$wt)
c(p2$par[1],(p2$par[1]-p2$par[3])/p2$par[2],p2$par[3])
rec = rec_tree(wt = st2$t,pars = p$par)
p3 <- subplex(par = p$par,fn = llik,n = rec$n, E = rec$E, t = rec$wt)
c(p3$par[1],(p3$par[1]-p3$par[3])/p3$par[2],p3$par[3])
###

# Simulations
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


par_est_vis(P=MP3,par=1,PR=RP3) ## looks good


  ## EM
KK = matrix(nrow=1000,ncol=3)
RP = matrix(nrow=1000,ncol=3)
init = c(4,0.175,1)
for(j in 1:1000){
  print(j)
  st = sim_phyl(seed = round(runif(1,0,j*100)))
  p <- subplex(par = c(8,0.175,0.9),fn = llik,n = st$n, E = st$E, t = st$t)$par
  RP[j,] = p
  sit = drop.fossil(st$newick)
  sit = phylo2p(sit)
  parsis = EM_phylo(bt=sit$t,init_par = init,n_trees=500,n_it=2,printpar = F)
  KK[j,]=parsis
}



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

#EM 2 pars version
st = sim_phyl(seed = 5)
p <- subplex(par = c(8,0.175,0.9),fn = llik,n = st$n, E = st$E, t = st$t)$par
sit = drop.fossil(st$newick)
sit = phylo2p(sit)
init = c(4,0.175)
n_trees = 100
par = init
for(j in 1:30){
  print(par)

  trees = vector('list',length=n_trees)
  for (i in 1:n_trees){
    rec = rec_tree(w=sit$t,pars=par,mu2par=0.09458395)
    trees[[i]] = rec
    setofweith[i] = -rec$prob
  }
  par = subplex(par = c(8,0.175), fn = llik_st , mu2par=0.09458395, setoftrees = trees, impsam = TRUE, setofweith=setofweith)$par
}
###


n_sim = 1000
n_trees = 10
MP2 = matrix(nrow=n_sim,ncol=2)
#MPf = matrix(nrow=n_sim,ncol=3)
RP2 = matrix(nrow=n_sim,ncol=3)
setofweith = NULL
for(j in 1:n_sim){
  print(j)
  st = sim_phyl(seed=j)
  p <- subplex(par = c(8,0.175,0.9),fn = llik,n = st$n, E = st$E, t = st$t)$par
  RP2[j,] = p
  par = c(p[1],p[2])
  sit = drop.fossil(st$newick)
  sit = phylo2p(sit)
  trees = vector('list',length=n_trees)
  for (i in 1:n_trees){
    rec = rec_tree(w=sit$t, pars=par,mu2par=p[3])
    setofweith[i] = -rec$prob
    trees[[i]] = rec
  }
  pars = subplex(par = c(8,0.175), fn = llik_st, mu2par=p[3], setoftrees = trees, impsam = FALSE)$par
  MP2[j,] = pars
}


par_est_vis(P=MP2,par=2,PR=RP2) ## looks good


### What is the influence of the extinction rate in the reconstructed tree


#NOTE that this code gives non consistent ltt plots... check it!
n_it=1000
SMP = vector(mode='list',length=16)
for (h in 2:16){
  print(k)
  Vals = vector(mode='list',length = n_it)
  Vals_rec = vector(mode='list',length = n_it)
  for (k in 1:n_it){
    #print(k)
    s1 = sim_phyl(seed=runif(1,1,10000),mu0 = 0.05*(h))
    Vals[[k]]=data.frame(time=s1$br,n=s1$n)
    s2 = drop.fossil(s1$newick)
    if(length(s2$tip.label)<5){
      k = k-1
    }
    else{
      s2 = phylo2p(s2)
      Vals_rec[[k]]=data.frame(time=cumsum(s2$t),n=c(1,s2$n))
    }

  }

  AA = matrix(nrow=length(seq(0,15,by=0.1)),ncol=n_it)
  A = matrix(nrow=length(seq(0,15,by=0.1)),ncol=n_it)
  for (i in 1:n_it){
    AA[,i]=approx(Vals[[i]]$time,Vals[[i]]$n,xout = seq(0,15,by=0.1))$y
    A[,i]=approx(Vals_rec[[i]]$time,Vals_rec[[i]]$n,xout = seq(0,15,by=0.1))$y
  }
  tt=15
  smp1 = cbind(data.frame(time=(seq(0,15,by=0.1)),val=(rowMeans(A, na.rm = TRUE, dims = 1)),mu='incomplete phylogeny'))
  smp2 = cbind(data.frame(time=(seq(0,14.9,by=0.1)),val=diff(rowMeans(A, na.rm = TRUE, dims = 1)),mu='diff incomplete phylogeny'))
  smp3 = cbind(data.frame(time=(seq(0,14.8,by=0.1)),val=diff(diff(rowMeans(A, na.rm = TRUE, dims = 1))),mu='diff diff incomplete phylogeny'))

  smp = rbind(smp1,smp2,smp3)
  SMP[[h]] = smp
  print(ggplot(smp,aes(x=time, y=val,color=mu))+geom_smooth()+ylab('number of lineages (log)')+xlab('time (Myr)')+ scale_y_log10()+ggtitle(paste('Extinction rate',as.character(0.05*h))))
}




Msmp = data.frame(SMP[[h]],val_mu = as.character(0.5))
for (h in 2:6){
  smp = data.frame(SMP[[h]],val_mu = as.character(0.5*h))
  Msmp = rbind(Msmp,smp)
}

print(ggplot(Msmp[Msmp$mu=='diff incomplete phylogeny',],aes(x=time, y=val,color=val_mu))+geom_smooth()+ylab('number of lineages (log)')+xlab('time (Myr)')+ scale_y_log10()+ggtitle(paste('Extinction rate',as.character(0.05*h))))

