# Purpose: Activation for binary regression models
# Updated: 180914

#' Activation for Binary Regression
#' 
#' @param h Linear predictor. 
#' @param model Selected from among logistic, probit, and robit.
#' @param df Degrees of freedom if the robit model is selected. 
#' @importFrom stats plogis pnorm pt

act.BinReg = function(h,model,df=NULL){
  # Calculate weights
  if(model=="logistic"){
    Out = plogis(q=h);
  } else if(model=="probit"){
    Out = pnorm(q=h);
  } else if(model=="robit"){
    Out = pt(q=h,df=df)
  }
  # Output
  return(Out);
}