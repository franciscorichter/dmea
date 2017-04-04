itexp <- function(u, m, t) { -log(1-u*(1-exp(-t*m)))/m }
rtexp <- function(n, m, t) { itexp(runif(n), m, t) }


rec_tree3 <- function(wt, model='dd',pars){
  lambda0 = pars[1]
  mu0 = pars[2]
  K = pars[3]
  n = 1:length(wt)
  i = 1
  E = rep(1,(length(wt)-1))
  ct = sum(wt)
  prob = 1
  p=list(wt=wt,E=E,n=n)
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
      break
    }
    cwt = wt[i]
    cbt = cumsum(wt)[i] - cwt
    t_spe = rexp(1,s)
    if (t_spe < cwt){
      prob = prob*dexp(x = t_spe,rate = s)
      t_ext = rtexp(1, mu0, ct-(cbt + t_spe))
      t_ext = cbt + t_spe + t_ext
      up = update_tree(p=p, t_spe = (cbt + t_spe), t_ext = t_ext)
      E = up$E
      n = up$n
      wt = up$wt
      p=list(wt=wt,E=E,n=n)
    }else{
      prob = prob*pexp(q = cwt,rate = s,lower.tail = F)
    }
    i = i+1
  }
  #L = create_L(wt,E)
  f_n = exp(-llik(pars=pars,n=n,E=E,t=wt))
  weight = f_n/prob
  return(list(wt=wt,E=E,n=n,f_n=f_n,weight=weight,prob=prob))
}
