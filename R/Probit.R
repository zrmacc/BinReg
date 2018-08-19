# Purpose: Fitting function for Probit regression
# Updated: 180818

#' @useDynLib Probit
#' @importFrom Rcpp sourceCpp
NULL

########################
# Probit IRLS
########################

#' Fit Probit Regression
#' 
#' @param y Binary 0/1 outcome vector.
#' @param Z Numeric model matrix.
#' @param alpha Significance level, for CIs.
#' @param eps Tolerance for Newton-Raphson iterations.
#' @param maxit Maximum number of NR iterations.
#' @importFrom stats dnorm pnorm qnorm
#' @importFrom methods new
#' @export 

fit.Probit = function(y,Z,alpha=0.05,eps=1e-6,maxit=10){
  # Objective function
  Q = function(b){
    # Linear predictor
    eta = fastMMp(Z,b);
    # Parition
    h0 = eta[y==0];
    h1 = eta[y==1];
    # Output
    Out = sum(pnorm(h1,log.p=T,lower.tail=T))+sum(pnorm(h0,log.p=T,lower.tail=F));
    return(Out);
  }
  # Observed information
  obsInfo = function(b){
    # Linear predictor
    eta = as.numeric(fastMMp(Z,b));
    # Current weights
    w = dnorm(eta)^2/(pnorm(eta)*pnorm(-eta));
    # Observed information
    J = diagQF(Z=Z,w=w);
    return(J);
  }
  # IRLS update
  Update = function(b){
    # Linear predictor
    eta = as.numeric(fastMMp(Z,b));
    # Current weights
    w = dnorm(eta)^2/(pnorm(eta)*pnorm(-eta));
    # Current response
    u = eta + (y-pnorm(eta))/dnorm(eta);
    # Proposed update
    Prop = WLS(Z=Z,w=w,y=u);
    # Output
    Out = list("b"=Prop);
    return(Out);
  }
  # Initialize
  theta0 = list();
  theta0$b = OLS(Z=Z,y=y);
  theta0$ll = Q(b=theta0$b);
  ## Newton-Raphson
  for(i in 1:maxit){
    # Propose update
    theta1 = Update(b=theta0$b);
    # Proposed objective
    theta1$ll = Q(b=theta1$b);
    # Accept first update,
    # Otherwise, check for improvement
    if(i==1){
      theta0 = theta1;
    } else {
      delta = theta1$ll-theta0$ll;
      if(delta>eps){
        theta0 = theta1;
      } else {
        break;
      }
    }
  }; # End NR
  ## Report
  if(i<maxit){
    cat(paste0(i-1," update(s) performed before tolerance limit."),"\n");
  } else {
    cat(paste0(i," update(s) performed without reaching tolerance limit."));
  }
  # Observed information
  J = obsInfo(b=theta0$b);
  Ji = fastInv(J);
  se = sqrt(diag(Ji));
  # Coefficient frame
  if(is.null(colnames(Z))){colnames(Z)=paste0("z",seq(1:ncol(Z)))};
  B = data.frame(colnames(Z),theta0$b,se);
  colnames(B) = c("Coeff","Est","SE");
  # CIs
  z = qnorm(1-alpha/2);
  B$L = B$Est-z*B$SE;
  B$U = B$Est+z*B$SE;
  B$p = 2*pnorm(abs(B$Est/B$SE),lower.tail=F);
  # Residuals
  E = y-pnorm(fastMMp(Z,B$Est));
  # Output
  Out = new(Class="fit",Model="Probit",Coefficients=B,Information=J,Residuals=E);
  return(Out);
}