#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;


//' @title Hellinger distance between columns of two matrices.
//'
//' @description
//' \code{HellingerMat} returns a matrix of Hellinger distances between
//' the columns of the matrix passed as argument.
//'
//' @details
//' This function calculates the Hellinger distances between the columns of two matrices
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
//' y <- matrix(1:12, nrow=3, ncol=5)
//' PHellingerMat(x, y)
//'
// [[Rcpp::export]]
NumericMatrix PHellingerMat(NumericMatrix A,NumericMatrix B) {
    NumericMatrix A2=clone(A); // don't overwrite inputs
    NumericMatrix B2=clone(B); // don't overwrite inputs
    
    unsigned int rows = A2.nrow();
    unsigned int Acols = A2.ncol();
    unsigned int Bcols = B2.ncol();
    
    NumericMatrix answer(Acols,Bcols);
    
    // reweight matrix A so each column sums to one, and take sqrt of each value
    for( unsigned int k = 0; k < Acols; k++){
        double Acol_tot = 0.0;
        
        for( unsigned int j = 0; j < rows; j++){
            Acol_tot = Acol_tot + A2(j, k);
        }
        
        for( unsigned int j = 0; j < rows; j++){
            A2(j , k) = A2(j , k) / Acol_tot;
        }
        
        for( unsigned int j = 0; j < rows; j++){
            A2(j , k) = std::sqrt(double(A2(j , k)));
        }
    }
    
    // same for matrix B
    for( unsigned int k = 0; k < Bcols; k++){
        double Bcol_tot = 0.0;
        
        for( unsigned int j = 0; j < rows; j++){
            Bcol_tot = Bcol_tot + B2(j, k);
        }
        
        for( unsigned int j = 0; j < rows; j++){
            B2(j , k) = B2(j , k) / Bcol_tot;
        }
        
        for( unsigned int j = 0; j < rows; j++){
            B2(j , k) = std::sqrt(double(B2(j , k)));
        }
    }
    
    
    
    
    // Do the main calculations
    const double sqrtHalf = std::sqrt(double(0.5));
    for( unsigned int j = 0; j < Acols; j++){
        
        for( unsigned int k = 0; k < Bcols; k++){
            
            double result = 0.0;
            
            for( unsigned int i = 0; i < rows; i++){
                result = result + (A2(i, j) - B2(i, k)) * (A2(i, j) - B2(i, k));
            }
            
            answer(j , k) = sqrtHalf * std::sqrt(result);
            
        }
    }
    
    return(answer);
}
