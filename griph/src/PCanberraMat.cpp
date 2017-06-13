#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;


//' @title Canberra distance between columns of two matrices.
//'
//' @description
//' \code{CanberraMat} returns a matrix of Canberra distances between
//' the columns of the matrix passed as argument.
//'
//' @details
//' This function calculates the Canberra distances between the columns of two matrices
//' given as argument. It is implemented in C++ for efficiency. For two argument
//' matrices with the same number of r rows and c1 , c2 columns, the return value is an c1-by-c2 matrix with
//' the pairwise c1i c2j distances.
//'
//' @param A  Numeric matrix (variables by features)
//' @param B  Numeric matrix (variables by features)
//'
//' @return A matrix of dimensions \code{ncol(A})-by-\code{ncol(B)}.
//'
//' @examples
//' x <- matrix(1:12, nrow=3, ncol=4)
//' y <- matrix(1:15, nrow=3, ncol=5)
//' PHellingerMat(x, y)
//'
// [[Rcpp::export]]
NumericMatrix PCanberraMat(NumericMatrix A, NumericMatrix B) {
    NumericMatrix A2 = clone(A); // don't overwrite inputs
    NumericMatrix B2 = clone(B); // don't overwrite inputs
    
    unsigned int rows = A2.nrow();
    unsigned int Acols = A2.ncol();
    unsigned int Bcols = B2.ncol();
    
    NumericMatrix answer(Acols,Bcols);
    
    // Do the main calculations
    double denominator=0;
    double result=0;
    for( unsigned int j = 0; j < Acols; j++){
        
        for( unsigned int k = 0; k < Bcols; k++){

            result = 0;
            
            for( unsigned int i = 0; i < rows; i++){
                denominator=  A2(i, j)  +  B2(i, k)   ; //Only for speedup in case of log2 transformed count matrices with pseudocount. Safer if line above
                //result += denominator > 0 ? (std::abs(A2(i, j) - B2(i, k)) / ( denominator )  ) : 0;
                denominator += denominator <= 0; //Same as line above but avoids branching (see: https://stackoverflow.com/questions/16777456/what-is-the-fastest-integer-division-supporting-division-by-zero-no-matter-what)
                result +=  (std::abs(A2(i, j) - B2(i, k)) /  denominator   ) ;
            }
            answer(j , k) =  result;
        }
    }
    
    return(answer);
}
