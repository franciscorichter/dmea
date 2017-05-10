llik_st = function(pars, setoftrees, impsam = F, mu=NULL, conditionToSurvival=FALSE){
  if(length(pars)==2){
    pars[3] = pars[2]
    pars[2] = mu
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
        weight = 1
        l[i] = llik(pars=pars,n=s$n,E=s$E,t=s$wt,conditionToSurvival=conditionToSurvival)*weight
      }
    }
    else{
      l[i] = llik(pars=pars,n=s$n,E=s$E,t=s$wt,conditionToSurvival=conditionToSurvival)
    }
  }
  L = sum(l)
  return(L)
}
