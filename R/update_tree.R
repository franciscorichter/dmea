update_tree <- function(p,t_spe,t_ext){
  wt = p$wt
  n = p$n
  E = p$E
  ct = sum(wt)
  if(t_ext > ct){
    stop('Extinction beyond present!')
  }
  if(t_ext < t_spe){
    stop('Speciation after extinction!')
  }
  # speciation
  K = length(wt)
  k = length(wt[cumsum(wt) < t_spe])
  n = c(n[1:(k+1)],(n[(k+1):K]+1))
  if((k+1)<K){
    lastbit = E[(k+1):(K-1)]
  }
  else{
    lastbit = NULL
  }
  E = c(E[0:k],1,lastbit)
  # (this can be witten in one line, is convinent?)
  if(k==0){
    wt = c(t_spe, wt[1]-t_spe, wt[2:K])
  }
  if(k > 0 & k < (K-1)){
    wt = c(wt[1:k], t_spe - sum(wt[1:k]), wt[k+1]-(t_spe - sum(wt[1:k])), wt[(k+2):K])
  }
  if(k == (K-1)){
    wt = c(wt[1:(K-1)],t_spe-sum(wt[1:(K-1)]),ct-t_spe)
  }
  #extinction
  K = length(wt)
  k = length(wt[cumsum(wt) < t_ext])
  n = c(n[1:(k+1)],(n[(k+1):K]-1))
  if((k+1)<K){
    lastbit = E[(k+1):(K-1)]
  }
  else{
    lastbit = NULL
  }
  E = c(E[0:k],0,lastbit)
  # (this can be witten in one line, is convinent?)
  if(k==0){
    wt = c(t_ext, wt[1]-t_ext, wt[2:K])
  }
  if(k > 0 & k < (K-1)){
    wt = c(wt[1:k], t_ext - sum(wt[1:k]), wt[k+1]-(t_ext - sum(wt[1:k])), wt[(k+2):K])
  }
  if(k == (K-1)){
    wt = c(wt[1:(K-1)],t_ext-sum(wt[1:(K-1)]),ct-t_ext)
  }
  return(list(wt=wt,E=E,n=n))
}




# add warning message when spe or ext is beyond present time
# update_tree <- function(wt, t_spe, t_ext, E, n, ct){
#   #adding speciation
#   #ct = sum(wt)
#   k = length(wt[cumsum(wt)<t_spe])
#   if(k<length(E)){
#   #  lastbitE = E[(k+1):length(E)]
#     lastbitN = c(n[k+1],n[(k+1):length(n)]+1)
#     lastbitt = wt[(k+2):length(wt)]
#   }else{
#     lastbitE = NULL
#     lastbitN = c(n[k+1],n[k+1]+1)#n[k]+1
#     lastbitt = NULL
#   }
#   wt = c(wt[0:k], t_spe-sum(wt[0:k]), wt[k+1]-(t_spe-sum(wt[0:k])), lastbitt)
#   E = c(E[0:k],1,lastbitE)
#   n = c(n[0:k],lastbitN)
#   #adding extinction
#   k = length(wt[cumsum(wt)<t_ext])
#   #print(k)
#   if(k<length(E)){
#     lastbitE = E[(k+1):length(E)]
#     lastbitN = c(n[k+1],n[(k+1):length(n)]-1)
#     lastbitt = wt[(k+2):length(wt)]
#   }else{
#     ## THIS IS NOT WORKING!!
#     lastbitE = NULL
#     lastbitN = c(n[k+1],n[(k+1):length(n)]-1)#n[k]-1
#     lastbitt = NULL
#   }
#   E = c(E[0:k],0,lastbitE)
#   n = c(n[0:k],lastbitN)
#   wt = c(wt[0:k], t_ext-sum(wt[0:k]), wt[k+1]-(t_ext-sum(wt[0:k])), lastbitt)
#   #print(n)
#   return(list(wt=wt,E=E,n=n))
# }
#
