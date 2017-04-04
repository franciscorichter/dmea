expectedLTT2 <- function(pars, ct=15, n_it= 10, sn_it = 100, median=FALSE, dim=NULL){ # it should also include model
  # add message that checks that do not have extinct species
  a = vector("list", n_it)
  # diversity dependence model. parameters are par=c(lambda,mu,K)
  time=proc.time()
  for (i in 1:(n_it)){
    s = sim_phyl(seed=i, ct = ct, lambda0 = pars[1], mu0 = pars[2], K = pars[3])$newick.extant.p
    a[[i]] = s$wt
  }
  proc.time()-time
  n.obs <- sapply(a, length)
  #m.obs = round(mean(n.obs))
  seq.max <- seq_len(max(n.obs))
  same_dim = 0
  if(!is.null(dim) & sum(n.obs ==  dim)>1 ){
    p_it = n_it #previous it
    n_it = sn_it
    a = a[1:n_it]
    for (i in (p_it+1):(n_it)){
      s = sim_phyl(seed=i, ct = ct, lambda0 = pars[1], mu0 = pars[2], K = pars[3])$newick.extant.p
      a[[i]] = s$wt
    }
    n.obs <- sapply(a, length)
    a = a[n.obs ==  dim]
    seq.max = seq_len(dim)
    same_dim = sum(n.obs ==  dim)
  }
  mat <- t(sapply(a, "[", i = seq.max))
  #mat = mat[,1:m.obs]
  expect = colMeans(mat,na.rm = TRUE)
  if(median){expect = colMedians(mat,na.rm = TRUE)}
  ce = cumsum(expect)
  bt=c(cumsum(expect[ce<15]),15)
  if(!is.null(dim) & sum(n.obs ==  dim)>1){
    bt = ce
  }
  expect = list(bt=bt, Ex = 1:length(cumsum(expect[ce<15])),same_dim=same_dim)
  return(expect)
}



