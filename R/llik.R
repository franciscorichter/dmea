#TODO: create the general one
llik = function(pars,n,E,t,conditionToSurvival=FALSE){
  b = c(pars[1],(pars[1]-pars[2])/pars[3],pars[2])
  t=t[1:(length(t)-1)]
  n = n[1:(length(n)-1)]
  sigma = n*(b[1]-b[2]*n + b[3]) #n-dimentional
  rho = pmax(b[1]*E-b[2]*n*E+b[3]*(1-E),0)
  if(conditionToSurvival){
    pars1 =  c(pars[1],(pars[1]-pars[2])/pars[3],pars[2])
    lx = 100 # lx is the size of your ode system. Set it to 10 times thenumber of species, should be at least the number of species
    soc = 1 # soc = 1 when you start at stem age, and 2 when you start atcrown age
    probsn = rep(0,lx)
    probsn[1] = 1 # change if other species at stem or crown age
    k = soc
    t1 = t[1]  # stem or crown age
    t2 = sum(t) # usually the present
    y = ode(y = probsn, times = c(t1,t2), func = dd_loglik_rhs, parms =c(pars1,k,ddep = 1),rtol = 1E-16,atol = 1E-16,method = 'ode45'); # ifode45 does not work, use lsoda
    probsn = y[2,2:(lx+1)]
    if(soc == 1) { aux = 1:lx }
    if(soc == 2) { aux = (2:(lx+1)) * (3:(lx+2))/6 }
    probsc = probsn/aux
    logliknorm = log(sum(probsc))
  }
  else{
    logliknorm = 0
  }
  l = -sum(-sigma*t+log(rho)) - logliknorm
  if(min(b)<0){l = Inf}
  return(l)
}

