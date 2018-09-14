# Purpose: Weights for binary regression models
# Updated: 180914

#' Weights for Binary Regression
#' 
#' @param h Linear predictor. 
#' @param model Selected from among logistic, probit, and robit.
#' @param df Degrees of freedom if the robit model is selected. 
#' @importFrom stats dlogis dnorm pnorm dt pt

weight.BinReg = function(h,model,df=NULL){
  # Calculate weights
  if(model=="logistic"){
    Out = dlogis(x=h);
  } else if(model=="probit"){
    Out = dnorm(x=h)^2/(pnorm(q=h)*pnorm(q=-h));
  } else if(model=="robit"){
    Out = dt(x=h,df=df)^2/(pt(q=h,df=df)*pt(q=-h,df=df));
  }
  # Output
  return(Out);
}