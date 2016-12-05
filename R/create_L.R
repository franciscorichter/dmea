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
      L[L$spec == ext_spec , 3] = bt[i]
    }
  }
  return(L)
}
