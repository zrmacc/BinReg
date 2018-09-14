# Purpose: Master fitting function
# Updated: 180914

#' Fit Binary Regression Model
#' 
#' @param y Binary 0/1 outcome vector.
#' @param X Numeric model matrix. Include an intercept. 
#' @param model Selected from among logistic, probit, and robit.
#' @param df Degrees of freedom. 
#' @param sig Significance level, for CIs.
#' @param eps Tolerance for Newton-Raphson iterations.
#' @param maxit Maximum number of NR iterations.
#' @param report Report fitting progress? 
#' 
#' @export 

fit.BinReg = function(y,X=NULL,model="logistic",df=NULL,sig=0.05,eps=1e-8,report=T){
  # Intput check
  n = length(y);
  if(!is.numeric(y)){stop("A numeric vector is expected for y.")};
  if(is.null(X)){X = array(1,dim=c(n,1)); colnames(X)="int";};
  if(!is.matrix(X)){stop("A numeric matrix is expected for X.")};
  Choices = c("logistic","probit","robit");
  if(!(model%in%Choices)){stop("Select model from among: logistic, probit, robit.")};
  if((model=="robit")&&is.null(df)){stop("If using the robit model, specify degrees of freedom.")};
  
  # Fitting
  if(model=="logistic"){
    return(fit.logistic(y=y,X=X,sig=sig,eps=eps,report=report));
  } else if(model=="probit"){
    return(fit.probit(y=y,X=X,sig=sig,eps=eps,report=report));
  } else {
    return(fit.robit(y=y,X=X,df=df,sig=sig,eps=eps,report=report));
  }
};