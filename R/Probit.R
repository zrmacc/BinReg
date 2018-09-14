# Purpose: Fitting function for probit regression
# Updated: 180913

#' Fit Probit Regression
#' 
#' @param y Binary 0/1 outcome vector.
#' @param X Numeric model matrix.
#' @param sig Significance level, for CIs.
#' @param eps Tolerance for Newton-Raphson iterations.
#' @param maxit Maximum number of NR iterations.
#' @param report Report fitting progress? 
#' 
#' @importFrom stats dnorm pnorm qnorm
#' @importFrom methods new
#' 
#' @return An object of class \code{fit} containing the regression coefficients,
#' information, and residuals.
#' @examples 
#' \dontrun{
#' set.seed(101);
#' # Design matrix
#' n = 1e3;
#' X = matrix(rnorm(4*n),nrow=n);
#' # Coefficient
#' b = c(1,-0.5,-0.5,0);
#' # Probit outcomes
#' y = rBinary(X=X,b=b,model="probit");
#' # Fit probit model
#' M = fit.probit(y=y,X=X);
#' show(M);
#' } 

fit.probit = function(y,X,sig=0.05,eps=1e-8,maxit=10,report=T){
  # Intput check
  if(!is.numeric(y)){stop("A numeric vector is expected for y.")};
  if(!is.matrix(X)){stop("A numeric matrix is expected for X.")};
  # Partition
  key0 = (y==0);
  key1 = !(key0);
  
  # Objective function
  Q = function(b){
    # Linear predictor
    eta = MMP(X,b);
    # Parition
    h0 = eta[key0];
    h1 = eta[key1];
    # Output
    Out = sum(pnorm(h1,log.p=T,lower.tail=T))+sum(pnorm(h0,log.p=T,lower.tail=F));
    return(Out);
  }
  
  # Information
  Info = function(b){
    # Linear predictor
    eta = as.numeric(MMP(X,b));
    # Current weights
    w = dnorm(eta)^2/(pnorm(eta)*pnorm(-eta));
    # Calculate A=X'WX
    A = diagQF(X=X,w=w);
    return(A);
  }
  
  # IRLS update
  Update = function(theta){
    # Current beta
    b0 = theta$b;
    # Initial objective
    q0 = Q(b0);
    # Linear predictor
    eta = as.numeric(MMP(X,b0));
    # Current weights
    w = dnorm(eta)^2/(pnorm(eta)*pnorm(-eta));
    # Current working vector
    u = eta+(y-pnorm(eta))/dnorm(eta);
    # Update beta
    b1 = fitWLS(y=u,X=X,w=w)$Beta;
    # Final objective
    q1 = Q(b1);
    # Increment
    d = q1-q0;
    # Output
    Out = list("b"=b1,"d"=d);
    return(Out);
  }
  
  # Initialize
  theta0 = list();
  theta0$b = fitOLS(y=y,X=X)$Beta;
  
  ## Maximzation
  for(i in 1:maxit){
    # Update
    theta1 = Update(theta0);
    # Accept if increment is positive
    if(theta1$d>0){
      theta0 = theta1;
      if(report){cat("Objective increment: ",signif(theta1$d,digits=3),"\n")}
    }
    # Terminate if increment is below tolerance
    if(theta1$d<eps){
      rm(theta1);
      break;
    }
  };
  
  ## Fitting report
  if(report){
    if(i<maxit){
      cat(paste0(i-1," update(s) performed before tolerance limit.\n\n"));
    } else {
      cat(paste0(i," update(s) performed without reaching tolerance limit.\n\n"));
    }
  };
  
  # Information
  J = Info(b=theta0$b);
  Ji = matInv(J);
  se = sqrt(diag(Ji));
  
  # Coefficient frame
  if(is.null(colnames(X))){colnames(X)=paste0("x",seq(1:ncol(X)))};
  B = data.frame(colnames(X),theta0$b,se);
  colnames(B) = c("Coeff","Point","SE");
 
  # CIs
  z = qnorm(1-sig/2);
  B$L = B$Point-z*B$SE;
  B$U = B$Point+z*B$SE;
  B$p = 2*pnorm(abs(B$Point/B$SE),lower.tail=F);
  
  # Residuals
  e = y-pnorm(MMP(X,theta0$b));
  # Output
  Out = new(Class="fit",Model="Probit",Coefficients=B,Information=J,Residuals=e);
  return(Out);
}