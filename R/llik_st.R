llik_st = function(pars, setoftrees, impsam = F, mu=0.1){
  if(length(pars)==2){
    pars[3] = mu
  }
  m = length(setoftrees)
  l = NULL # this should be changed
  for(i in 1:m){
    s = setoftrees[[i]]
    if(impsam){
      weight = exp(s$weight)
      if(weight == 0){
        l[i] = 0
      }
      else{
        l[i] = llik(b=pars,n=s$n,E=s$E,t=s$wt)*weight
      }
    }
    else{
      l[i] = llik(b=pars,n=s$n,E=s$E,t=s$wt)
    }
  }
  L = sum(l)
  return(L)
}
