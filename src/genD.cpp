#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double findFunc(arma::vec funcdetail,
                arma::rowvec data1,
                arma::rowvec data2,
                arma::vec pars){
  //generate Euclidean distance
  double val = 0;
  arma::uword f1 = funcdetail(0);
  
  for(arma::uword j=0;j<4;j++){
    arma::sword idx = funcdetail(j+1);
    if(idx >= 0){
      val += pow(data1(idx) - data2(idx),2);
    }
  }
  val = sqrt(val);
  double out;
  //group
  if(f1==1){
    if(val==0){
      out = pars(funcdetail(5));
    } else {
      out = 0;
    }
  }
  
  //exponential 1
  if(f1==2){
    out = pars(funcdetail(5))*exp(-1*val*pars(funcdetail(6)));
  }
  
  //power exponential
  if(f1==3){
    out = pow(pars(funcdetail(5)),val);
  }
  
  return(out);
}

// [[Rcpp::export]]
arma::mat genSubD(arma::mat &func,
                  arma::mat &data,
                  arma::vec pars){
  arma::uword nF = func.n_cols;
  arma::uword N = data.n_rows;
  
  //initialise the matrix
  arma::mat D(N,N);
  D.fill(1);
  double val;
  //loop through the matrix and produce all the elements
  for(arma::uword j=0; j<N; j++){
    for(arma::uword k=j+1; k < N; k++){
      for(arma::uword l=0; l<nF; l++){
        val = findFunc(func.col(l),
                       data.row(j),
                       data.row(k),
                       pars);
        D(j,k) = D(j,k)*val;
        D(k,j) = D(k,j)*val;
      }
      
    }
    for(arma::uword l=0; l<nF; l++){
      val = findFunc(func.col(l),
                     data.row(j),
                     data.row(j),
                     pars);
      D(j,j) = D(j,j)*val;
    }
  }
  return(D);
}

// [[Rcpp::export]]
arma::mat blockMat(arma::mat mat1,
                   arma::mat mat2){
  arma::uword n1 = mat1.n_rows;
  arma::uword n2 = mat2.n_rows;
  arma::mat dmat(n1+n2,n1+n2);
  dmat.fill(0);
  dmat.submat(0,0,n1-1,n1-1) = mat1;
  dmat.submat(n1,n1,n1+n2-1,n1+n2-1) = mat2;
  return(dmat);
}

// [[Rcpp::export]]
arma::mat genD(Rcpp::List &func,
               Rcpp::List &data,
               arma::vec pars){
  //arma::uword n = func.n_slices;
  int n = func.length();
  arma::mat f1 = func[0];
  arma::mat d1 = data[0];
  arma::mat D1 = genSubD(f1,d1,pars);
  if(n > 1){
    for(int j=1;j<n;j++){
      arma::mat f1 = func[j];
      arma::mat d1 = data[j];
      arma::mat D2 = genSubD(f1,d1,pars);
      D1 = blockMat(D1,D2);
    }
  }
  return(D1);
}

