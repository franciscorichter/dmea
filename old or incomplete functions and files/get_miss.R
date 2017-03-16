get_miss <- function(tobs,trec){
  trec=trec[2:length(trec)]
  count = NULL
  for (i in 1:(length(tobs)-1)){
    count[i] = length( trec[cumsum(trec)>cumsum(tobs)[i] & cumsum(trec) <= cumsum(tobs)[i+1] ] )
  }
  return(count)
}
