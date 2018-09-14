# Purpose: Objective function for binary regression models
# Updated: 180914

#' Objective Function for Binary Regression
#' 
#' @param y Outcome vector.
#' @param h Linear predictor. 
#' @param model Selected from among logistic, probit, and robit.
#' @param df Degrees of freedom if the robit model is selected. 
#' @import stats

Q.BinReg = function(y,h,model,df=NULL){
  # Partition
  key0 = (y==0);
  key1 = !(key0);
  # Parition
  h0 = h[key0];
  h1 = h[key1];
  # Calculate objective
  if(model=="logistic"){
    Out = sum(plogis(h1,log.p=T,lower.tail=T))+sum(plogis(h0,log.p=T,lower.tail=F));
  } else if(model=="probit"){
    Out = sum(pnorm(h1,log.p=T,lower.tail=T))+sum(pnorm(h0,log.p=T,lower.tail=F));
  } else if(model=="robit"){
    Out = sum(pt(q=h1,df=df,log.p=T,lower.tail=T))+sum(pt(q=h0,df=df,log.p=T,lower.tail=F));
  }
  # Output
  return(Out);
}