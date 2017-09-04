/* Pearson correlation coefficients between columns of a sparse matrix */

#include <R.h>
#define USE_RINTERNALS
#include <Rinternals.h>
#include <Rdefines.h>

//#include <Rcpp.h>
#include <cmath>

#include "Matrix.h"
#include "Matrix_stubs.c"

#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

// using namespace Rcpp;

// Based on code from R  (R/src/cov.c:corcov and R/src/library/stats/src/cov.c:cov_complete2),
// with modifications:
// - expect a matrix of class Matrix::dgCMatrix
// - only supporting single matrix version
// - only supporting Pearson correlation
// - not supporting NAs
// - assuming that there is long double
// - assuming that columns standard deviations are non-zero
// For details on how to work with a dgCMatrix in C, see t_gCMatrix_colSums.c in Matrix sources
//
//' @title Pearson's correlation coefficient between columns of a sparse matrix
//'     (OpenMP version).
//'
//' @description
//' \code{SPearsonMatOMP} returns a symmetric matrix of Pearson's correlation
//' coefficients between all pairs of columns of the sparse matrix passed as argument.
//'
//' @details
//' This function calculates the Pearson's correlation coefficients between all
//' pairs of columns of the sparse matrix given as argument. It is implemented
//' in C++ for efficiency. For a matrix with c columns, the return
//' value is an c-by-c matrix with the pairwise ci cj correlation coefficients.
//'
//' @param A  Numeric matrix (variables by features). Must be a \code{dgCMatrix}.
//' @param nthreads Integer specifying the number of OpenMP threads to use in
//'     parallel parts (defaults to two, ignored on systems not supporting
//'     OpenMP).
//'
//' @return A matrix of dimensions \code{ncol(A})-by-\code{ncol(B)}.
//'
//' @examples
//' x <- matrix(c(1,0,4,0,0,1,2,5,0,3,0,0,1,3,0), nrow=5)
//' cor(x)
//' xs <- as(x, "sparseMatrix")
//' SPearsonMatOMP(xs)
//'
// [[Rcpp::export]]
SEXP SPearsonMatOMP(SEXP A, int nthreads = 2) {
    Rf_warning("SPearsonMatOMP is not functional yet - use PPearsonMatOMP on dense matrices instead");
    
    CHM_SP cx = AS_CHM_SP__(A);
    // R_CheckStack();

    int i, i1, i2, j, k, nc = cx->ncol, nr = cx->nrow;
    int nr1 = nr - 1;
    int *xp = (int *)(cx -> p), *xi = (int *)(cx -> i);
    double *xx = (double *)(cx -> x), *ans, *xm, xm1, xm2;
    long double sum;
    SEXP ans_;

    // get data pointers and allocate memory
    PROTECT(ans_ = Rf_allocMatrix(REALSXP, nc, nc));
    ans = REAL(ans_);
    xm = (double *) malloc(nc * sizeof(double));

#ifdef _OPENMP
    // set the number of threads
    int nthreads_old = omp_get_num_threads();
    omp_set_num_threads(nthreads);
#endif

    // just count number of pairwise non-zero elements
    for (j = 0 ; j < nc ; j++) {
        for (k = 0 ; k <= j ; k++) {
            double cnt = 0.0;
            for (i1 = xp[j], i2 = xp[k];
                 i1 < xp[j+1] && i2 < xp[k+1]; ) {
                if (xi[i1] < xi[i2]) {
                    i1++;
                } else if (xi[i2] < xi[i1]) {
                    i2++;
                } else {
                    cnt++;
                    i1++;
                    i2++;
                }
            }
            ans[ j + k*nc ] = ans[ k + j*nc ] = cnt;
        }
    }
    /*    
    // calculate column means -> xm[]
#ifdef _OPENMP
#pragma omp parallel for private(sum, j, i)
#endif
    for (j = 0; j < nc; j++) {
        sum = 0.;
        for(i = xp[j]; i < xp[j + 1]; i++)
            sum += xx[i];
        xm[j] = (double)(sum / nr);
    }

    // calculate correlation coefficients
#ifdef _OPENMP
#pragma omp parallel for private(sum, xm1, xm2, j, k, i, i1, i2)
#endif
    for (j = 0 ; j < nc ; j++) {
        xm1 = xm[j];
        for (k = 0 ; k <= j ; k++) {
            xm2 = xm[k];
            sum = 0.;
            for (i1 = xp[j], i2 = xp[k], i = 0 ; i < nr && i1 < xp[j+1] && i2 < xp[k+1]; i++) {
                if (xi[i1] == i) {
                    if (xi[i2] == i) {
                        sum += (xx[i1] - xm1) * (xx[i2] - xm2);
                        i1++;
                        i2++;
                    } else {
                        sum += (xx[i1] - xm1) * (0. - xm2);
                        i1++;
                    }
                } else if (xi[i2] == i) {
                    sum += (0. - xm1) * (xx[i2] - xm2);
                    i2++;
                }
            }
            //printf("j:%d (%d-%d) k:%d (%d-%d) = %.3f\n", j+1,xp[j],xp[j+1]-1,k+1,xp[k],xp[k+1]-1,(double)(sum / nr1));
            ans[j + k * nc] = ans[k + j * nc] = (double)(sum / nr1);
        }
    }

    // clean up
#ifdef _OPENMP
    // reset threads
    omp_set_num_threads(nthreads_old);
#endif
    free(xm);
     */
    UNPROTECT(1);

    return ans_;
}
#undef CLAMP

/*

    // BELOW HERE IS OLD    
    int n = Rf_nrows(A), ncx = Rf_ncols(A);
    int i, j, k, n1 = n - 1;
    double *x, *y, *xm, *ym, *xx, *yy, *ans;
    long double sum, tmp, xxm, yym;
    SEXP ans_;
    
    
    // calculate standard deviations --> xm[] and ym[]
#ifdef _OPENMP
#pragma omp parallel for private(sum, xx, xxm, i, k)
#endif
    for (i = 0 ; i < ncx ; i++) { // sd(A[,i])
        xx = &x[i * n];
        sum = 0.;
        xxm = xm [i];
        for (k = 0 ; k < n ; k++)
            sum += (xx[k] - xxm) * (xx[k] - xxm);
        sum /= n1;
        xm [i] = (double)sqrtl(sum);
    }
    
#ifdef _OPENMP
#pragma omp parallel for private(sum, yy, yym, i, k)
#endif
    for (i = 0 ; i < ncy ; i++) { // sd(B[,i])
        yy = &y[i * n];
        sum = 0.;
        yym = ym [i];
        for (k = 0 ; k < n ; k++)
            sum += (yy[k] - yym) * (yy[k] - yym);
        sum /= n1;
        ym [i] = (double)sqrtl(sum);
    }

#define CLAMP(X)  (X >= 1. ? 1. : (X <= -1. ? -1. : X))
#ifdef _OPENMP
#pragma omp parallel for private(i, j)
#endif
    for (i = 0 ; i < ncx ; i++) {
        for (j = 0 ; j < ncy ; j++) {
            ANS(i,j) /= (xm[i] * ym[j]);
            ANS(i,j) = CLAMP(ANS(i,j));
        }
    }
     
    // set dimnames() - just copied from covcor; needs rewriting
    // x = getAttrib(x, R_DimNamesSymbol);
    // y = getAttrib(y, R_DimNamesSymbol);
    // if ((length(x) >= 2 && !isNull(VECTOR_ELT(x, 1))) ||
    //     (length(y) >= 2 && !isNull(VECTOR_ELT(y, 1)))) {
    //     PROTECT(ind = allocVector(VECSXP, 2));
    //     if (length(x) >= 2 && !isNull(VECTOR_ELT(x, 1)))
    //         SET_VECTOR_ELT(ind, 0, duplicate(VECTOR_ELT(x, 1)));
    //     if (length(y) >= 2 && !isNull(VECTOR_ELT(y, 1)))
    //         SET_VECTOR_ELT(ind, 1, duplicate(VECTOR_ELT(y, 1)));
    //     setAttrib(ans, R_DimNamesSymbol, ind);
    //     UNPROTECT(1);
    // }

    // clean up
#ifdef _OPENMP
    // reset threads
    omp_set_num_threads(nthreads_old);
#endif
    free(xm);
    free(ym);
    UNPROTECT(3);
    
    return ans_;
}
#undef ANS

 */


//' @title sparse-sparse verion of PCanberraMatOMP
//' 
//' @description
//' \code{ssPCanberraMatOMP} is an alternative version of \code{PCanberraMatOMP}
//' for sparse inputs.
//'
//' @param A  Sparse numeric matrix (variables by features). Must be a \code{dgCMatrix}.
//' @param B  Sparse numeric matrix (variables by features). Must be a \code{dgCMatrix}.
//' @param nthreads Integer specifying the number of OpenMP threads to use in
//'     parallel parts (defaults to two, ignored on systems not supporting
//'     OpenMP).
//'
//' @return A matrix of dimensions \code{ncol(A})-by-\code{ncol(B)}.
//'
//' @examples
//' x <- matrix(c(1,0,4,0,0,1,2,5,0,3,0,0,1,3,0), nrow=5)
//' xs <- as(x, "sparseMatrix")
//' 
//' PCanberraMatOMP(x, x)
//' ssPCanberraMatOMP(xs, xs)
//'
// [[Rcpp::export]]
SEXP ssPCanberraMatOMP(SEXP A, SEXP B, int nthreads = 2) {
    CHM_SP cA = AS_CHM_SP__(A);
    CHM_SP cB = AS_CHM_SP__(B);
    unsigned int rows = cA->nrow, Acols = cA->ncol, Bcols = cB->ncol;
    
    int *Ap = (int *)(cA -> p), *Ai = (int *)(cA -> i);
    int *Bp = (int *)(cB -> p), *Bi = (int *)(cB -> i);
    double *Ax = (double *)(cA -> x), *Bx = (double *)(cB -> x);
    
    SEXP answer;
    PROTECT(answer = Rf_allocMatrix(REALSXP, Acols, Bcols));
    double *ans = REAL(answer);
    
#ifdef _OPENMP
    // set the number of threads
    int nthreads_old = omp_get_num_threads();
    omp_set_num_threads(nthreads);
#endif
    
    // Do the main calculations
    // idea 1: on-the fly create dense vectors for each column (~20% slower than dense)
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (unsigned int j = 0; j < Acols; j++) {
        
        double denominator, result;
        double *a = (double *)malloc(rows * sizeof(double));
        double *b = (double *)malloc(rows * sizeof(double));
        
        memset(a, 0, rows * sizeof(double));
        for (unsigned int n = Ap[j]; n < Ap[j+1]; n++)
            a[Ai[n]] = Ax[n];  // a = A[,j]
        
        for (unsigned int k = 0; k < Bcols; k++) {
            memset(b, 0, rows * sizeof(double));
            for (unsigned int n = Bp[k]; n < Bp[k+1]; n++)
                b[Bi[n]] = Bx[n];  // b = B[,k]
            result = 0.0;
            
            for (unsigned int i = 0; i < rows; i++) {
                denominator = a[i] + b[i];
                // result += (denominator == 0 ? 0 : std::abs(a[i] - b[i]) / denominator);
                denominator += denominator <= 0; // avoids branching (see: https://stackoverflow.com/questions/16777456/what-is-the-fastest-integer-division-supporting-division-by-zero-no-matter-what)
                result += (std::abs(a[i] - b[i]) / denominator);
            }
            ans[j + k * Acols] = result;
        }
        free(a);
        free(b);
    }

    /* // idea 2: iterate only over non-zero elements, need if-else; seems ~4-fold slower than dense and currently not correct yet
     double result;
#ifdef _OPENMP
#pragma omp parallel for private(result)
#endif
     for( unsigned int j = 0; j < Acols; j++){
        
        for( unsigned int k = 0; k < Bcols; k++){
            
            result = 0.0;
            for (unsigned int i1 = Ap[j], i2 = Bp[k]; i1 < Ap[j+1] || i2 < Bp[k+1]; ) {
                if (Ai[i1] < Bi[i2] && i1 < Ap[j+1]) {
                    result += 1.0;
                    // printf("%d, %d|%d : 1.00\n",Ai[i1]+1,j+1,k+1);
                    i1++;
                } else if (Ai[i1] > Bi[i2] && i2 < Bp[k+1]) {
                    result += 1.0;
                    // printf("%d, %d|%d : 1.00\n",Bi[i2]+1,j+1,k+1);
                    i2++;
                } else if (i1 < Ap[j+1] && i2 < Bp[k+1]) { // without this if, result is a bit off; with it, it runs indefenetly...
                    result += (std::abs(Ax[i1] - Bx[i2]) / (Ax[i1] + Bx[i2]));
                    // printf("%d, %d|%d : %.2f\n",Ai[i1]+1,j+1,k+1,(std::abs(Ax[i1] - Bx[i2]) / (Ax[i1] + Bx[i2])));
                    i1++;
                    i2++;
                }
            }
            ans[j + k * Acols] =  result;
        }
    }
     */

#ifdef _OPENMP
    // reset threads
    omp_set_num_threads(nthreads_old);
#endif
    
    UNPROTECT(1);
    return(answer);
}

//' @title sparse-dense verion of PCanberraMatOMP
//' 
//' @description
//' \code{sdPCanberraMatOMP} is an alternative version of \code{PCanberraMatOMP}
//' for sparse and dense inputs.
//'
//' @param A  Sparse numeric matrix (variables by features). Must be a \code{dgCMatrix}.
//' @param B  Dense numeric matrix (variables by features). Must be a \code{matrix}.
//' @param nthreads Integer specifying the number of OpenMP threads to use in
//'     parallel parts (defaults to two, ignored on systems not supporting
//'     OpenMP).
//'
//' @return A matrix of dimensions \code{ncol(A})-by-\code{ncol(B)}.
//'
//' @examples
//' x <- matrix(c(1,0,4,0,0,1,2,5,0,3,0,0,1,3,0), nrow=5)
//' xs <- as(x, "sparseMatrix")
//' 
//' PCanberraMatOMP(x, x)
//' sdPCanberraMatOMP(xs, x)
//'
// [[Rcpp::export]]
SEXP sdPCanberraMatOMP(SEXP A, SEXP B, int nthreads = 2) {
    CHM_SP cA = AS_CHM_SP__(A);
    B = PROTECT(Rf_coerceVector(B, REALSXP));
    unsigned int rows = cA->nrow, Acols = cA->ncol, Bcols = Rf_ncols(B);
    
    int *Ap = (int *)(cA -> p), *Ai = (int *)(cA -> i);
    double *Ax = (double *)(cA -> x), *Bx = REAL(B);
    
    SEXP answer;
    PROTECT(answer = Rf_allocMatrix(REALSXP, Acols, Bcols));
    double *ans = REAL(answer);
    
#ifdef _OPENMP
    // set the number of threads
    int nthreads_old = omp_get_num_threads();
    omp_set_num_threads(nthreads);
#endif
    
    // Do the main calculations
    // on-the fly create dense vector for each sparse column (~as fast as dense)
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (unsigned int j = 0; j < Acols; j++) {
        
        double denominator, result;
        double *a = (double *)calloc(rows, sizeof(double)), *b;
        for (unsigned int n = Ap[j]; n < Ap[j+1]; n++)
            a[Ai[n]] = Ax[n];  // a = A[,j]
        
        unsigned int k;
        for (k = 0, b = Bx; k < Bcols; k++, b += rows) { // b = B[,k]
            result = 0.0;
            
            for (unsigned int i = 0; i < rows; i++) {
                denominator = a[i] + b[i];
                // result += (denominator == 0 ? 0 : std::abs(a[i] - b[i]) / denominator);
                denominator += denominator <= 0; // avoids branching (see: https://stackoverflow.com/questions/16777456/what-is-the-fastest-integer-division-supporting-division-by-zero-no-matter-what)
                result += (std::abs(a[i] - b[i]) / denominator);
            }
            ans[j + k * Acols] = result;
        }
        free(a);
    }
    
#ifdef _OPENMP
    // reset threads
    omp_set_num_threads(nthreads_old);
#endif
    
    UNPROTECT(2);
    return(answer);
}
