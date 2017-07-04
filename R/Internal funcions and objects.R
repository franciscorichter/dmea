sl = paste(letters[1],letters,":0",sep="")
for (i in 2:26){
  ll = paste(letters[i],letters,":0",sep="")
  sl = c(sl,ll)
}
SL = paste(LETTERS[1],LETTERS,":0",sep="")
for (i in 2:26){
  ll = paste(LETTERS[i],LETTERS,":0",sep="")
  SL = c(SL,ll)
}
S1 = paste(LETTERS[1],0:9,":0",sep="")
for (i in 2:26){
  ll = paste(LETTERS[i],1:10,":0",sep="")
  SL = c(SL,ll)
}
s1 = paste(letters[1],0:9,":0",sep="")
for (i in 2:26){
  ll = paste(letters[i],1:10,":0",sep="")
  SL = c(SL,ll)
}
sl=c(sl,SL,S1,s1)
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



###


par_est_vis <- function(P,par,PR){
  # P is the recost estim values
  # par is the parameter you want to use
  # PR is the real estim
  if (par == 1){
    int = 0.8 #change for general case
    parname = 'lambda'
  }
  if (par==2){
    int= 0.1
    parname = 'mu'
  }
  if (par == 3){
    int = 40
    parname = 'K'
    n = dim(P)[1]
    P1 = P[P[,3]<100 & PR[,3]<100,]
    PR = PR[P[,3]<100 & PR[,3]<100,]
    P = P1
    n2 = dim(P)[1]
    print(paste(1-n2/n,'proportion of data was excluded for vizualization purposes'))
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
  p = phylo2p(phy)
  ct = sum(p$wt)
  n <- Ntip(phy)
  x <- dist.nodes(phy)[n + 1, ][1:n]
  dphy = drop.tip(phy, root.edge = T , which(x < max(x) - tol))
  p2 = phylo2p(dphy)
  if(sum(p2$wt) != ct){
    rem_time = ct - sum(p2$wt)
    p2$wt[1] = p2$wt[1] + rem_time
    dphy = p2phylo(p2)
  }
  return(dphy)
}

# obsolete?
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

sim_est <- function(n_trees, pars, init_par=c(1.8,0.13,60),impsam=FALSE,rec_method=1,seed=0,conditionToSurvival=FALSE){ # simulate a tree, drop fossil, and estimate parameters back after bootstrap reconstruction
  if (seed != 0) set.seed(seed)
  st = dmea::sim_phyl(lambda0 = pars[1], mu0 = pars[2], K = pars[3])
  p <- subplex(par = init_par, fn = llik, n = st$n, E = st$E, t = st$wt)$par
  sit = st$tree.extant
  trees = sim_srt(wt=sit$wt, pars=p, parallel = F, n_trees = n_trees, rec_method=rec_method)
  pars = subplex(par = init_par, fn = llik_st , setoftrees = trees, impsam = impsam)$par
  return(data.frame(real=p, est=pars))
}

number_missing <- function(st){
  total = length(st$tree$wt)
  extant = length(st$tree.extant$wt)
  missing = (total - extant)/2
  return(missing)
}

create_L <- function(t,E){
  ## THIS FUNCTION HAS A BUG WHEN THE LENGTH OF T IS SO HIGH
  L = data.frame(spec='aa', spec_time=0, ext_time=-1)
  bt = cumsum(t)
  for(i in 1:(length(t)-1)){
    if(E[i] == 1){
      L = rbind(L,data.frame(spec=substr(sl[i+1],1,2), spec_time=bt[i], ext_time = -1))
    }
    if(E[i] == 0){
      set = L[L$ext_time == (-1),]
      ext_spec = sample(set$spec,1) # diversity-dependence
      L[L$spec == ext_spec , 3] = bt[i]
    }
  }
  return(L)
}

# obsolete?
get_g <- function(wt,wt_obs){
  Sp = split(cumsum(wt), findInterval(cumsum(wt), cumsum(wt_obs), left.open = TRUE))
  for (i in 1:length(Sp)){

  }
}

convol <-function(wt,lambda,mu,remt){
  out = 1-exp(-lambda*wt)-lambda*exp(-mu*remt)/(mu-lambda)*(exp(wt*(mu-lambda))-1)
  return(out)
}



get_comb_ltt <- function(phylo1,phylo2,phylo=TRUE){
  # Trabajar aca para generalizarlo a medias de arboles...con phylo=F
  if(phylo){
    t1 = data.frame(ltt.plot.coords(phylo1))[1,1]
    t2 = data.frame(ltt.plot.coords(phylo2))[1,1]
    if(round(t1) != round(t2)){print('DIFFERENT LENGHT IN PHYLOGENIES!')}
    ltt1 = data.frame(ltt.plot.coords(phylo1)[-1,])
    ltt2 = data.frame(ltt.plot.coords(phylo2)[-1,])
  }
  i = j = 1
  N1 = dim(ltt1)[1]
  N2 = dim(ltt2)[1]
  comb_ltt = data.frame(time=NULL, n1=NULL, n2=NULL)
  while(i <= N1 & j <= N2 ){
    check = T
    if(ltt1$time[i] < ltt2$time[j]){
      comb_ltt = rbind(comb_ltt,data.frame(time=ltt1$time[i]-t1,n1=ltt1$N[i],n2=ltt2$N[j]))
      t1 = ltt1$time[i]
      i = i+1
      check = F
    }
    if(ltt1$time[i] > ltt2$time[j] & check){
      comb_ltt = rbind(comb_ltt,data.frame(time=ltt2$time[j]-t1,n1=ltt1$N[i],n2=ltt2$N[j]))
      t1 = ltt2$time[j]
      j=j+1
      check = F
    }
    if(ltt1$time[i] == ltt2$time[j] & check){
      comb_ltt = rbind(comb_ltt,data.frame(time=ltt2$time[j]-t1,n1=ltt1$N[i],n2=ltt2$N[j]))
      t1 = ltt2$time[j]
      j=j+1
      i=i+1
    }
  }
  return(comb_ltt)
}



get_comb_ltt2 <- function(sot,phyl=NULL){  #pd=phylo_data
  m = length(sot)
  time = NULL
  ind=NULL
  for(i in 1:m){
    s=sot[[i]]
    s=phylo2p(s)
    s=s$wt
    s = cumsum(s)
    s = s[-length(s)]
    time =c(time,s)
    ind = c(ind,rep(i,length(s)))
  }
  if(!is.null(phyl)){
    s = cumsum(phyl)
    s = s[-length(s)]
    time = c(time,s)
    ind = c(ind,rep((m+1),length(s)))
    m=m+1
  }
  df = data.frame(time=time,ind=ind)
  df = df[order(time),]
  df = cbind(df,as.data.frame(matrix(ncol=m,nrow=dim(df)[1])))
  d2 = as.data.frame(matrix(ncol=(m+2),nrow=1))
  names(d2) = names(df)
  vec=rep(1,m)
  d2[1,] =  c(0,0,vec)
  df = rbind(d2,df)
  for(i in 2:dim(df)[1]){
    # esto puede ser escrito de forma mas eficiente sin cambiar todos los valores cada vez
    j = df$ind[i]
    vec[j] = vec[j] + 1
    df[i,3:(m+2)] = vec
  }
  df['mean'] = rowMeans(df[,3:(m+2)])
  return(df)
}



ltt_mu <- function(mu,phylo,prior_pars,n_trees=10){
  po = phylo2p(phylo)
  trees = sim_srt(wt=po$wt,pars = prior_pars,n_trees = n_trees,mu = mu)
  lambda_K = subplex(par =c(2,60), fn = llik_st , setoftrees = trees, mu = mu, impsam = FALSE)$par
  pars = c(lambda_K[1],mu,lambda_K[2])
  expe = expectedLTT2(pars,n_it=n_trees)
  wt = c(expe$bt[1],diff(expe$bt))
  p = list(wt=wt,E=rep(1,length(expe$bt)),n=expe$Ex)
  #Ltt2$mu = as.factor(Ltt2$mu)
  phylo2 = p2phylo(p)
  ltt = ltt_stat(phylo,phylo2)
  return(list(ltt=ltt,pars=pars))
}



