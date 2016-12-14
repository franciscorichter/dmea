p2phylo <- function(t,n,E){
  newick = paste(sl[1],";",sep="")
  for (i in 1:length(t)){
    nt = read.tree(text=newick)
    # speciation
    if (E[i] == 1){
      spec = sample(1:n[i],1)
      species = nt$tip.label[spec]
      ind = regexpr(species,newick)[1]-1
      atm = ...
      newick = paste(substr(newick,1,ind),"(",substr(newick,ind+1,ind+4),",",sl[i+1],"):",as.character(atm),substring(newick,ind+5),sep="")
    }
    # extinction
    if (E[i]==0){
      species = identf[BD-N,1]
      ind = regexpr(species,newick)[1] + 2
      atm = sumt-identf[which(identf[,1]==species),2]
      identf = identf[-(BD-N),]
      L[L$spec == species,]$ext_time = sumt
      newick = paste(substr(newick,1,ind),as.character(atm),substring(newick,ind+2),sep="")
    }
  }

}
