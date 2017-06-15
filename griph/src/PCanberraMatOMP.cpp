#include <Rcpp.h>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

using namespace Rcpp;


//' @title Canberra distance between columns of two matrices (OpenMP version).
//'
//' @description
//' \code{PCanberraMatOMP} returns a matrix of Canberra distances between
//' the columns of the two matrices passed as arguments.
//'
//' @details
//' This function calculates the Canberra distances between the columns of two
//' matrices given as argument. It is implemented in C++ for efficiency. For two
//' argument matrices with the same number of r rows and c1 , c2 columns, the
//' return value is an c1-by-c2 matrix with the pairwise c1i c2j distances.
//'
//' @param A  Numeric matrix (variables by features)
//' @param B  Numeric matrix (variables by features)
//' @param nthreads Integer specifying the number of OpenMP threads to use in
//'     parallel parts (defaults to two, ignored on systems not supporting
//'     OpenMP).
//'
//' @return A matrix of dimensions \code{ncol(A})-by-\code{ncol(B)}.
//'
//' @examples
//' x <- matrix(1:12, nrow=3, ncol=4)
//' y <- matrix(1:15, nrow=3, ncol=5)
//' PHellingerMatOMP(x, y)
//'
// [[Rcpp::export]]
NumericMatrix PCanberraMatOMP(NumericMatrix A, NumericMatrix B, int nthreads = 2) {
    unsigned int rows = A.nrow();
    unsigned int Acols = A.ncol();
    unsigned int Bcols = B.ncol();
    
    NumericMatrix answer(Acols,Bcols);
    
#ifdef _OPENMP
    // set the number of threads
    int nthreads_old = omp_get_num_threads();
    omp_set_num_threads(nthreads);
#endif
    
    // Do the main calculations
    double denominator=0;
    double result=0;
#ifdef _OPENMP
#pragma omp parallel for private(denominator, result)
#endif
    for( unsigned int j = 0; j < Acols; j++){
        
        for( unsigned int k = 0; k < Bcols; k++){

            result = 0;
            
            for( unsigned int i = 0; i < rows; i++){
                denominator=  A(i, j)  +  B(i, k)   ; //Only for speedup in case of log2 transformed count matrices with pseudocount. Safer if line above
                //result += denominator > 0 ? (std::abs(A(i, j) - B(i, k)) / ( denominator )  ) : 0;
                denominator += denominator <= 0; //Same as line above but avoids branching (see: https://stackoverflow.com/questions/16777456/what-is-the-fastest-integer-division-supporting-division-by-zero-no-matter-what)
                result +=  (std::abs(A(i, j) - B(i, k)) /  denominator   ) ;
            }
            answer(j , k) =  result;
        }
    }
    
#ifdef _OPENMP
    // reset threads
    omp_set_num_threads(nthreads_old);
#endif
    
    return(answer);
}
