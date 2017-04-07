#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;


//' @title Hellinger distance between columns of a matrix.
//'
//' @description
//' \code{HellingerMat} returns a matrix of Hellinger distances between
//' the columns of the matrix passed as argument.
//'
//' @details
//' This function calculates the Hellinger distances between columns of the matrix
//' given as argument. It is implemented in C++ for efficiency. For an argument
//' matrix with r rows and c columns, the return value is an c-by-c matrix with
//' all pairwise distances.
//'
//' @param A Numeric matrix.
//'
//' @return A square matrix of dimensions \code{ncol(A})-by-\code{ncol(A)}.
//'
//' @examples
//' x <- matrix(1:12, nrow=3, ncol=4)
//' HellingerMat(x)
//'
// [[Rcpp::export]]
NumericMatrix HellingerMat(NumericMatrix A) {
    NumericMatrix A2 = A; // don't overwrite inputs

    unsigned int rows = A2.nrow();
    unsigned int cols = A2.ncol();
    NumericMatrix answer(cols,cols);

    // reweight so each column sums to one, and take sqrt of each value
    for( unsigned int k = 0; k < cols; k++){
        double col_tot = 0.0;

        for( unsigned int j = 0; j < rows; j++){
            col_tot = col_tot + A2(j, k);
        }

        for( unsigned int j = 0; j < rows; j++){
            A2(j , k) = A2(j , k) / col_tot;
        }

        for( unsigned int j = 0; j < rows; j++){
            A2(j , k) = std::sqrt(double(A2(j , k)));
        }
    }


    // Do the main calculations
    const double sqrtHalf = std::sqrt(double(0.5));
    for( unsigned int j = 0; j < cols - 1; j++){

        for( unsigned int k = j + 1; k < cols; k++){

            double result = 0.0;

            for( unsigned int i = 0; i < rows; i++){
                result = result + (A2(i, j) - A2(i, k)) * (A2(i, j) - A2(i, k));
            }

            answer(j , k) = answer(k , j) = sqrtHalf * std::sqrt(result);

        }
    }

    return(answer);
}
