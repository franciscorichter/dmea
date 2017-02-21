#TODO: create the general one
llik = function(b,n,E,t){
  t=t[1:(length(t)-1)]
  n = n[1:(length(n)-1)]
  sigma = n*(b[1]-b[2]*n + b[3]) #n-dimentional
  rho = pmax(b[1]*E-b[2]*n*E+b[3]*(1-E),0)
  l = -sum(-sigma*t+log(rho))#+0.04/(b[3]*b[3])
  if(min(b)<0){l = Inf}
  return(l)
}

