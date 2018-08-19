// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

//' Weighted Least Squares
//'
//' Estimates the weighted least squares coefficient 
//' \deqn{\hat{\beta}=(Z'WZ)^{-1}Z'Wy}
//' 
//' The weight matrix is assumed diagonal, and is supplied
//' as a vector. 
//'
//' @param Z Design matrix.
//' @param w Weight vector.
//' @param y Response vector.
//' @return Numeric vector. 
//' @export 
//'
// [[Rcpp::export]]
SEXP WLS(const Eigen::Map<Eigen::MatrixXd> Z, const Eigen::Map<Eigen::VectorXd> w,
          const Eigen::Map<Eigen::VectorXd> y){
  // Dimensions
  const int n = Z.rows();
  const int p = Z.cols();
  // Calculate A=Z'WZ
  Eigen::MatrixXd A(p,p);
  A.setZero();
  for(int i=0; i<n; i++){
    A += Z.row(i).transpose()*w(i)*Z.row(i);
  };
  // Calculate b=Z'Wy
  Eigen::VectorXd b(p);
  b.setZero();
  for(int i=0; i<n; i++){
    b += Z.row(i).transpose()*w(i)*y(i);
  };
  // Solve
  const Eigen::VectorXd Out = A.ldlt().solve(b);
  return Rcpp::wrap(Out);
}

//' Ordiary Least Squares
//'
//' Estimates the ordinary least squares coefficient 
//' \deqn{\hat{\beta}=(Z'Z)^{-1}Z'y}
//'
//' @param Z Design matrix.
//' @param y Response vector.
//' @return Numeric vector. 
//' @export 
//'
// [[Rcpp::export]]
SEXP OLS(const Eigen::Map<Eigen::MatrixXd> Z, const Eigen::Map<Eigen::VectorXd> y){
  // Components
  const Eigen::MatrixXd A = (Z.transpose()*Z);
  const Eigen::VectorXd b = (Z.transpose()*y);
  const Eigen::VectorXd Out = A.llt().solve(b);
  return Rcpp::wrap(Out);
}