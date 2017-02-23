ltt_stat <- function(phylo1,phylo2){
  comb = get_comb_ltt(phylo1,phylo2)
  ltt_stat = sum(abs((comb$n1-comb$n2)*comb$time))
  return(ltt_stat)
}
