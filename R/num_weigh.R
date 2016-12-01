num_weigh <- function(rec_tree,pars_rec){
  L = rec_tree$L
  miss <- L[L[,3]!=(-1),]
  miss_spe <- miss$spec
  time <- rbind(data.frame(brtimes = L[,2], E=1, spec = L[,1]),data.frame(brtimes = miss[,3], E=0, spec = miss[,1]))
  time <- time[order(time$brtimes),]
  time$n <- c(0,time$n[1:length(time$n)-1])
  missing = time[which(is.element(time$spec,miss_spe)),]
  if(dim(missing)[1]<1){print('THERE IS NOT MISSING SPECIES')}
  m_spec = unique(missing$spec)
  P = matrix(nrow = length(m_spec), ncol=1)
  for(j in 1:length(m_spec)){
    sub_mat = missing[missing$spec == m_spec[j],]
    # DESDE ACA HAY QUE seguir.. ALMOST DONE!
    P[i] = 
  }
  s = n*(lambda-beta*n)*(1/n)
  formula = s*(1/n)
}


