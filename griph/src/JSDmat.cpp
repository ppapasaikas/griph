#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' @title Jensen-Shannon divergence between rows of a matrix.
//'
//' @description
//' \code{JSDmat} returns a matrix of Jensen-Shannon divergences between the rows
//' of the matrix passed as argument.
//'
//' @details
//' This function calculates the Jensen-Shannon divergences between row of the matrix
//' given as argument. It is implemented in C++ for efficiency. For an argument
//' matrix with r rows and c columns, the return value is an r-by-r matrix with
//' all pairwise distances, where only the upper-triangle of the matrix is filled.
//'
//' @param A Numeric matrix.
//'
//' @return A square matrix of dimensions \code{nrow(A})-by-\code{nrow(A)}.
//'
//' @examples
//' x <- matrix(1:12, nrow=3, ncol=4)
//' JSDmat(t(x))
//'
// [[Rcpp::export]]
NumericMatrix JSDmat(NumericMatrix A){

   NumericMatrix A2 = A; // don't overwrite inputs

    int rows = A2.nrow();
    int cols = A2.ncol();
    NumericMatrix answer(rows,rows);

    // add 10^(-4) so that there are no zero entries
    for( int j = 0; j < rows; j++){
      for( int k = 0; k < cols; k++){
        A2(j, k) = A2(j, k) + 0.0001;
      }
    }

    // reweight so each row sums to one
    for( int j = 0; j < rows; j++){
      double row_tot = 0.0;

      for(int k = 0; k < cols; k++){
       row_tot = row_tot + A2(j, k);
      }

      for(int k = 0; k < cols; k++){
        A2(j , k) = A2(j , k) / row_tot;
      }
    }


    // Do the main calculations
    for(int j = 0; j < rows - 1; j++){

      for(int k = j + 1; k < rows; k++){

        NumericVector p(cols), q(cols);

        for(int i = 0; i < cols; i++){
          p[ i ] = A2(j , i);
          q[ i ] = A2(k , i);
        }

        NumericVector m = 0.5 * (p + q);

        // get the KL divergence for p||m and q||m
        double kl_pm = 0;
        double kl_qm = 0;

        for( int i = 0; i < cols; ++i ){
          kl_pm = kl_pm + std::log(p[ i ] / m[ i ]) * p[ i ];
          kl_qm = kl_qm + std::log(q[ i ] / m[ i ]) * q[ i ];
        }

        answer(j , k) = 0.5 * (kl_pm + kl_qm);

      }
    }

    return(answer);
}