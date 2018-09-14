# Purpose: Binary data generation
# Updated: 180913

#' Random Binary Outcome
#' 
#' Simulates outcomes from a binary regression model. Model choices include
#' logistic, probit, and robit.
#' 
#' @param X Design matrix.
#' @param b Regression coefficient.
#' @param model Selected from among logistic, probit, and robit.
#' @param df Degrees of freedom if the robit model is selected. 
#' 
#' @importFrom stats rnorm rt runif
#' @export
#' 
#' @return Numeric vector of binary outcomes, one for each row of the design
#'   matrix.

rBinReg = function(X,b,model="logistic",df=NULL){
  # Input check
  n = nrow(X);
  if(!is.matrix(X)){stop("A numeric matrix is expected for X.")};
  if(!is.vector(b)){stop("A numeric vector is expected for b.")};
  Choices = c("logistic","probit","robit");
  if(!(model%in%Choices)){stop("Select model from among: logistic, probit, robit.")};
  if((model=="robit")&&is.null(df)){stop("If using the robit model, specify degrees of freedom.")};
  
  # Linear predictor
  h = as.numeric(MMP(X,b));
  # Latent variable
  if(model=="logistic"){
    u = runif(n);
    z = h+log(u/(1-u));
  } else if(model=="probit"){
    z = h+rnorm(n);
  } else {
    z = rt(n=n,df=df,ncp=h);
  };
  # Dichotomize
  y = 1*(z>0);
  return(y);
}