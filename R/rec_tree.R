rec_tree <- function(wt, pars=c(0.8,0.0175,0.1), model='dd',rec_method=1,v2=F){
  if(rec_method==1){
    rec = rec_tree1(wt=wt,pars=pars,model=model,v2=v2)
  }
  if(rec_method==2){
    rec = rec_tree2(wt=wt,pars=pars,model=model)
  }
  if(rec_method==3){
    rec = rec_tree3(wt=wt,pars=pars,model=model)
  }
  if(is.na(rec$weight)) rec$weight = 0
  return(rec)
}
