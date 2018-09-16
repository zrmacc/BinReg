# Purpose: Master fitting function
# Updated: 180914

#' Fit Binary Regression Model
#' 
#' @param y Binary 0/1 outcome vector.
#' @param X Numeric model matrix. Include an intercept. 
#' @param model Selected from among logistic, probit, and robit.
#' @param offset Fixed component added to the linear predictor. 
#' @param df Degrees of freedom, if using the robit model. 
#' @param sig Significance level, for CIs.
#' @param eps Tolerance for Newton-Raphson iterations.
#' @param maxit Maximum number of NR iterations.
#' @param report Report fitting progress? 
#' 
#' @importFrom methods new
#' @export 
#' @examples
#' \dontrun{
#' set.seed(100);
#' # Design matrix
#' X = cbind(1,matrix(rnorm(3e3),nrow=1e3));
#' # Coefficient
#' b = c(1,-1,1,0);
#' 
#' # Logistic observations
#' y = rBinReg(X,b,model="logistic");
#' # Estimate logistic model
#' M = Fit.BinReg(y,X,model="logistic");
#' 
#' # Probit observations
#' y = rBinReg(X,b,model="probit");
#' # Estimate probit model
#' M = Fit.BinReg(y,X,model="probit");
#' 
#' # Robit observations
#' y = rBinReg(X,b,model="robit",df=5);
#' # Estimate robit model
#' M = Fit.BinReg(y,X,model="robit",df=5);
#' }

Fit.BinReg = function(y,X=NULL,model="logistic",offset=0,df=NULL,sig=0.05,eps=1e-8,maxit=10,report=T){
  # Intput check
  n = length(y);
  if(!is.numeric(y)){stop("A numeric vector is expected for y.")};
  if(is.null(X)){X = array(1,dim=c(n,1)); colnames(X)="int";};
  if(!is.matrix(X)){stop("A numeric matrix is expected for X.")};
  Choices = c("logistic","probit","robit");
  if(!(model%in%Choices)){stop("Select model from among: logistic, probit, robit.")};
  if((model=="robit")&&is.null(df)){stop("If using the robit model, specify degrees of freedom.")};
  
  # Objective function
  Q = function(h){Q.BinReg(y=y,h=h,model=model,df=df)};
  # Information function
  Info = function(h){
    # Current weights
    w = as.numeric(weight.BinReg(h=h,model=model,df=df));
    # Calculate A=X'WX
    A = matIP(X,w*X);
    return(A);
  }
  # Upate function
  Update = function(theta){
    # Current beta
    b0 = theta$b;
    # Current linear predictor
    h0 = MMP(X,b0)+offset;
    # Initial objective
    q0 = Q(h0);
    # Current weights
    w = weight.BinReg(h=h0,model=model,df=df);
    # Current means
    m = act.BinReg(h=h0,model=model,df=df);
    # Current working vector
    delta = delta.BinReg(h=h0,model=model,df=df);
    u = (h0-offset)+(y-m)/(delta);
    # Update beta
    b1 = fitWLS(y=u,X=X,w=w)$Beta;
    # Updated linear predictor
    h1 = MMP(X,b1)+offset;
    # Final objective
    q1 = Q(h1);
    # Increment
    d = q1-q0;
    # Output
    Out = list("b"=b1,"d"=d);
    return(Out);
  }
  
  # Initialize
  theta0 = list();
  theta0$b = fitOLS(y=(y-offset),X=X)$Beta;
  
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
  
  # Final coefficient
  b1 = theta0$b;
  # Final linear predictor
  h1 = MMP(X,b1)+offset;
  # Information
  J = Info(h=h1);
  Ji = matInv(J);
  se = sqrt(diag(Ji));
  
  # Coefficient frame
  if(is.null(colnames(X))){colnames(X)=paste0("x",seq(1:ncol(X)))};
  B = data.frame(colnames(X),b1,se);
  colnames(B) = c("Coeff","Point","SE");
  
  # CIs
  z = qnorm(1-sig/2);
  B$L = B$Point-z*B$SE;
  B$U = B$Point+z*B$SE;
  B$p = 2*pnorm(abs(B$Point/B$SE),lower.tail=F);
  
  # Residuals
  m = act.BinReg(h=h1,model=model,df=df);
  e = (y-m);
  # Output
  Out = new(Class="fit",Model=model,Df=list("df"=df),Coefficients=B,Eta=h1,Information=J,Residuals=e);
  return(Out);
};