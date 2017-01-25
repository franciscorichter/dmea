sl = paste(letters[1],letters,":0",sep="")
for (i in 2:26){
  ll = paste(letters[i],letters,":0",sep="")
  sl = c(sl,ll)
}
sl=c(sl,sl)
compphyl <- function(newi,identf,ct){
  #set to extant species to the present time
  identf[,1] = as.character(identf[,1])
  identf[,2] = ct-identf[,2]
  for(i in 1:length(identf[,1])){
    ind = regexpr(identf[i,1],newi)[1] + 2
    newi = paste(substr(newi,1,ind),as.character(identf[i,2]),substring(newi,ind+2),sep="")
  }
  return(newi)
}

# TODO: add newick output to update_tree
# add warning message when spe or ext is beyond present time
update_tree <- function(wt, t_spe, t_ext, E, n){
  #adding speciation
  ct = sum(wt)
  k = length(wt[cumsum(wt)<t_spe])
   if(k<length(E)){
     lastbitE = E[(k+1):length(E)]
     lastbitN = c(n[k+1],n[(k+1):length(n)]+1)
     lastbitt = wt[(k+2):length(wt)]
   }else{
     lastbitE = NULL
     lastbitN = c(n[k+1],n[k+1]+1)#n[k]+1
     lastbitt = NULL
   }
  wt = c(wt[0:k], t_spe-sum(wt[0:k]), wt[k+1]-(t_spe-sum(wt[0:k])), lastbitt)
  E = c(E[0:k],1,lastbitE)
  n = c(n[0:k],lastbitN)
  #adding extinction
  k = length(wt[cumsum(wt)<t_ext])
  #print(k)
  if(k<length(E)){
    lastbitE = E[(k+1):length(E)]
    lastbitN = c(n[k+1],n[(k+1):length(n)]-1)
    lastbitt = wt[(k+2):length(wt)]
  }else{
    lastbitE = NULL
    lastbitN = c(n[k+1],n[(k+1):length(n)]-1)#n[k]-1
    lastbitt = NULL
  }
  E = c(E[0:k],0,lastbitE)
  n = c(n[0:k],lastbitN)
  wt = c(wt[0:k], t_ext-sum(wt[0:k]), wt[k+1]-(t_ext-sum(wt[0:k])), lastbitt)
  #print(n)
  return(list(wt=wt,E=E,n=n))
}



###


par_est_vis <- function(P,par,PR){
  #P is the recost estim values
  # par is the parameter you want to use
  #PR is the real estim
  if (par == 1){
    int = 0.8 #change for general case
    parname = 'lambda'
  }
  if (par==2){
    int= 0.0175
    parname = 'beta'
  }
  if (par == 3){
    int = 0.1
    parname = 'mu'
  }
  if (par ==4){
    int = 40
    parname = 'K'
    P = P[P[,4]<100,]
  }
  hist_top <- ggplot()+geom_histogram(aes(P[,par])) + geom_vline(xintercept=int) + xlab('MLE from incomplete tree')
  empty <- ggplot()+geom_point(aes(1,1), colour="white")+
    theme(axis.ticks=element_blank(),
          panel.background=element_blank(),
          axis.text.x=element_blank(), axis.text.y=element_blank(),
          axis.title.x=element_blank(), axis.title.y=element_blank())

  scatter <- ggplot()+geom_point(aes(P[,par], PR[,par]))+ geom_abline(intercept = 0, slope = 1) + ylab(TeX(paste('$\\hat{\\',parname,'}_{C}$',sep='')))+xlab(TeX(paste('$\\hat{\\',parname,'}_{I}$',sep='')))
  hist_right <- ggplot()+geom_histogram(aes(PR[,par]))+coord_flip()+ geom_vline(xintercept=int) +xlab('MLE of complete tree')

  grid.arrange(hist_top, empty, scatter, hist_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
}

drop.fossil <- function (phy, tol = 1e-08)
{
  n <- Ntip(phy)
  x <- dist.nodes(phy)[n + 1, ][1:n]
  drop.tip(phy,root.edge = T ,which(x < max(x) - tol))
}


num_weigh <- function(rec, pars_rec, ct){
  lambda0 = pars_rec[1]
  beta = pars_rec[2]
  mu = pars_rec[3]
  L = rec$L
  wt = rec$wt
  miss <- L[L[,3] != (-1),]
  miss_spe <- miss$spec
  time <- rbind(data.frame(brtimes = L[,2], E=1, spec = L[,1]),data.frame(brtimes = miss[,3], E=0, spec = miss[,1]))
  time <- time[order(time$brtimes),]
  time$n <- cumsum(time$E)-cumsum(1-time$E)-1
  wait_time = c(diff(time$brtimes))
  time <- time[2:dim(time)[1],]
  time$wt <- wait_time
  missing = time[which(is.element(time$spec,miss_spe)),]
  if(dim(missing)[1]<1){print('THERE IS NOT MISSING SPECIES')}
  m_spec = unique(missing$spec)
  P = matrix(nrow = length(m_spec), ncol=1)
  for(i in 1:length(m_spec)){
    sub_mat = missing[missing$spec == m_spec[i],]
    if(m_spec[i] == 'aa') sub_mat = rbind( data.frame(brtimes=0,E=1,spec='aa',n=1,wt=0),missing[missing$spec == m_spec[i],])
    lambda = max(0,lambda0-beta*sub_mat$n[1])
    t_ext = sub_mat$brtimes[2]-sub_mat$brtimes[1]
    P[i] = log(lambda) - sub_mat$wt[1]*sub_mat$n[1]*lambda +log(mu) - t_ext*mu-log(1-exp(-mu*(ct-sub_mat$brtimes[1])))
  }
  prob = sum(P)
  return(prob)
}


sim_est <- function(n_trees, init_par=c(8,0.175,0.9),impsam=FALSE,rec_method=1,seed=runif(1,1,100000)){ # simulate a tree, drop fossil, and estimate parameters back after bootstrap reconstruction
  set.seed(seed)
  st = dmea::sim_phyl()
  p <- subplex(par = init_par, fn = llik,n = st$n, E = st$E, t = st$t)$par
  sit = dmea::drop.fossil(st$newick)
  sit = dmea::phylo2p(sit)
  trees = sim_srt(wt=sit$t, pars=p, parallel = F, n_trees = n_trees,rec_method=rec_method)
  pars = subplex(par = init_par, fn = llik_st , setoftrees = trees, impsam = impsam)$par
  return(data.frame(real=p, est=pars))
}


number_missing <- function(st){
  total = length(st$t)
  drop = drop.fossil(st$newick)
  drop = phylo2p(drop)
  extant = length(drop$t)
  missing = (total - extant)/2
  return(missing)
}


create_L <- function(t,E){
  L = data.frame(spec='aa', spec_time=0, ext_time=-1)
  bt = cumsum(t)
  for(i in 1:(length(t)-1)){
    if(E[i] == 1){
      L = rbind(L,data.frame(spec=substr(sl[i+1],1,2), spec_time=bt[i], ext_time = -1))
    }
    if(E[i] == 0){
      set = L[L$ext_time == (-1),]
      ext_spec = sample(set$spec,1) # diversity-dependence
      #print(ext_spec)
      #print(L)
      L[L$spec == ext_spec , 3] = bt[i]
    }
  }
  return(L)
}
