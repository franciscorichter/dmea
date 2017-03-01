rec_tree1 <- function(wt, pars, model='dd', v2=F){
  lambda0 = pars[1]
  mu0 = pars[2]
  K = pars[3]
  n = 1:length(wt)
  i = 1
  tails = tail(n,n=1)
  E = rep(1,(length(wt)-1))
  fake = FALSE
  ct = sum(wt)
  prob = NaN #this has a better alternative, no?
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
      t_ext = rexp(1,mu0)
      t_ext = cbt + t_spe + t_ext
      if (t_ext < ct){
        p = list(wt=wt,E=E,n=n)
        up = update_tree(p,t_spe = (cbt + t_spe), t_ext = t_ext)
        E = up$E
        n = up$n
        wt = up$wt
        fake = FALSE
        prob[i] = dexp(t_ext,rate=mu0,log=TRUE) + dexp(t_spe,rate=s,log = TRUE)
      }else{
        prob[i] = pexp(q = ct,rate = mu0,lower.tail = F,log.p = TRUE) + dexp(t_spe,rate=s,log = TRUE)
        if(v2){ prob[i] = log(convol(wt = t_spe,lambda = s,mu = mu,remt = ct-cbt))}
        fake = TRUE
        i = i-1
      }
    }else{
      fake = FALSE
      prob[i] = pexp(q = cwt,rate = s,lower.tail = F,log.p = TRUE)
      if(v2){prob[i] = log(convol(wt = t_spe,lambda = s,mu = mu,remt = ct-cbt))}
    }
    i = i+1
  }
  L = create_L(wt,E)
  f_n = -llik(pars=pars,n=n,E=E,t=wt)
  logweight = f_n-sum(prob)
  return(list(wt=wt,E=E,n=n,weight=logweight,L=L,g=prob,f_n=f_n))
}
