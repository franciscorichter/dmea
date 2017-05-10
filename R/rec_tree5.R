rec_tree5 <- function(wt, pars, model='dd', v2=F){
  lambda0 = pars[1]
  mu0 = pars[2]
  K = pars[3]
  limit = floor((lambda0*K)/(lambda0-mu0))
  n = 1:length(wt)
  i = 1
  tails = tail(n,n=1)
  E = rep(1,(length(wt)-1))
  fake = FALSE
  ct = sum(wt)
  prob = NaN #this has a better alternative, no?
  m = 0 #number of missing species
  M = data.frame(spe=numeric(0),ext=numeric(0))
  acumulated_time = 0
  while(i < length(wt)){
    N = n[i]+m
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
    s = sum(lambda+mu)
    spext = TRUE
    if(s==0){
      #print('s=0')
      break
    }
    cwt = wt[i]+acumulated_time # current waiting time
    cbt = cumsum(wt)[i] - cwt  # current branching time
    rt = ct - cbt # remaining time
    nwt = rexp(1,s) # new waiting time
    if (nwt < cwt){
      if(m==0){ #speciation
        M2 = data.frame(spe=nwt+cbt,ext=NaN)
        M=rbind(M,M2)
        m=sum(is.nan(M[,2]))
        spex = FALSE
      }
      if(N==limit){
        M[is.nan(M[,2]),][1,2] = nwt+cbt
        m=sum(is.nan(M[,2]))
        spext = FALSE
      }
      if(spex){
        # draw event
        tem_mu = mu/pexp(rt,rate=mu)
        lambdat = lambda[1]/(lambda[1]+tem_mu[1]) # not general
        mut = tem_mu[1]/(lambda[1]+tem_mu[1])
        event = sample(c('s','e'),size=1,prob = c(lambdat,mut))
        if(event=='s'){
          M2 = data.frame(spe=nwt+cbt,ext=NaN)
          M=rbind(M,M2)
          m=sum(is.nan(M[,2]))
        }
        if(event=='e'){
          M[is.nan(M[,2]),][1,2] = nwt+cbt
          m=sum(is.nan(M[,2]))
        }
      }
      spex=TRUE
      acumulated_time = acumulated_time + nwt

    }
    else{
      i=i+1
      acumulated_time = 0
    }

  }
  #
  #     if (t_ext < ct){
  #       p = list(wt=wt,E=E,n=n)
  #       up = update_tree(p,t_spe = (cbt + t_spe), t_ext = t_ext)
  #       E = up$E
  #       n = up$n
  #       wt = up$wt
  #       fake = FALSE
  #       prob[i] = dexp(t_ext,rate=mu0,log=TRUE) + dexp(t_spe,rate=s,log = TRUE)
  #     }else{
  #       prob[i] = pexp(q = ct,rate = mu0,lower.tail = F,log.p = TRUE) + dexp(t_spe,rate=s,log = TRUE)
  #       if(v2){ prob[i] = log(convol(wt = t_spe,lambda = s,mu = mu,remt = ct-cbt))}
  #       fake = TRUE
  #       i = i-1
  #     }
  #   i = i+1
  # }
  # f_n = -llik(pars=pars,n=n,E=E,t=wt)
  # logweight = f_n-sum(prob)
  #return(list(wt=wt,E=E,n=n,weight=logweight,L=L,g=prob,f_n=f_n))
  return(M)
}

