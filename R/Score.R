# Purpose: Score Test for Binary Regression Models
# Updated: 180913

#' Score Test for Binary Regression
#' 
#' Tests the hypothesis that a subset of the regression coefficients are fixed
#' at a reference value. Specifically, let \eqn{\beta} denote the regression
#' coefficient. Partition \eqn{\beta=(\beta_{1},\beta_{2})}. Suppose that 
#' interest lies in testing that \eqn{\beta_{1}} is fixed at \eqn{\beta_{10}}. 
#' \code{Score.BinReg} performs a score test of
#' \eqn{H_{0}:\beta_{1}=\beta_{10}}. The test is specified using a logical vector
#' \code{L}, with as many entries as columns in the model matrix \code{X}. The
#' values of \code{L} set to \code{T} are constrained under the null, while
#' values of \code{L} set to \code{F} are estimated under the null.
#' 
#' @param y Numeric response vector.
#' @param X Numeric model matrix. 
#' @param L Logical vector, with as many entires as columns in the model matrix,
#'   indicating which columns have fixed coefficients under the null.
#' @param b10 Value of the regression coefficient for the selected columns under
#'   the null. Defaults to zero.
#' @param model Selected from among logistic, probit, and robit.
#' @param df Degrees of freedom, if using the robit model. 
#' @param sig Significance level, for CIs.
#' @param eps Tolerance for Newton-Raphson iterations.
#' @param maxit Maximum number of NR iterations.
#' @param report Report fitting progress? 
#' 
#' @return A numeric vector containing the score statistic, the degrees of 
#'   freedom, and a p-value. 
#'   
#' @export
#' 
#' @examples 
#' \dontrun{
#' set.seed(101);
#' # Design matrix
#' X = cbind(1,matrix(rnorm(n=4*1e3),nrow=1e3));
#' # Regression coefficient
#' b = c(1,-1,2,-1,0);
#' # Logistic outcome
#' y = rBinReg(X,b,model="logistic");
#' # Test b1=b2=b3=b4=0, which is false.
#' Score.BinReg(y=y,X=X,L=c(F,T,T,T,T),model="logistic",report=F);
#' # Test b4=0, which is true.
#' Score.BinReg(y=y,X=X,L=c(F,F,F,F,T),model="logistic",report=F);
#' # Test b2=0 and b4=2, which if false.
#' Score.BinReg(y=y,X=X,L=c(F,F,T,F,T),b10=c(0,2),model="logistic",report=F);
#' # Test b1=b3=-1, which is true. 
#' Score.BinReg(y=y,X=X,L=c(F,T,F,T,F),b10=c(-1,-1),model="logistic",report=F);
#' }

Score.BinReg = function(y,X,L,b10=NULL,model="logistic",df=NULL,
                        sig=0.05,eps=1e-8,maxit=10,report=T){
  ## Input checks
  if(!is.vector(y)){stop("A numeric vector is expected for y.")};
  if(!is.matrix(X)){stop("A numeric matrix is expected for X.")};
  if(!is.logical(L)){stop("A logical vector is expected for L.")};
  # Test specification
  q = ncol(X);
  k = sum(L);
  if(length(L)!=q){stop("L should have as many entries as columns in X.")};
  if(k==0){stop("At least 1 entry of L should be TRUE.")};
  if(k==q){stop("At least 1 entry of L should be FALSE.")};
  # Check for missingness
  Miss = sum(is.na(y))+sum(is.na(X));
  if(Miss>0){stop("Neither y nor Z, should contain missing values.")};
  # Null coefficient
  if(is.null(b10)){b10=rep(0,times=k)};
  if(length(b10)!=k){stop("b10 should contain a reference value for each TRUE element of L.")};
  # Model
  Choices = c("logistic","probit","robit");
  if(!(model%in%Choices)){stop("Select model from among: logistic, probit, robit.")};
  if((model=="robit")&&is.null(df)){stop("If using the robit model, specify degrees of freedom.")};
  
  # Partition target design
  Xa = X[,L,drop=F];
  Xb = X[,!L,drop=F];
  # Null model
  M0 = Fit.BinReg(y=y,X=Xb,model=model,offset=MMP(Xa,b10),df=df,sig=sig,eps=eps,maxit=maxit,report=report);
  # Extract residuals
  e0 = resid(M0,type="raw");
  # Fitted linear predictor
  h0 = predict(M0,type="eta");
  # Calculate weights
  w0 = weight.BinReg(h=h0,model=model,df=df);
  # Calculate delta
  d0 = delta.BinReg(h=h0,model=model,df=df);
  # Score vector
  u = matIP(Xa,w0*e0/d0);
  # Information
  Ibb = matIP(Xa,w0*Xa);
  Iaa = M0@Information;
  Iba = matIP(Xa,w0*Xb);
  # Schur complement
  V = SchurC(Ibb,Iaa,Iba);
  # Score statistic
  Ts = matQF(X=u,A=matInv(V));
  # P value
  p = pchisq(q=Ts,df=k,lower.tail=F);
  # Output
  Out = c(Ts,k,p);
  names(Out) = c("Score","df","p");
  return(Out);
}