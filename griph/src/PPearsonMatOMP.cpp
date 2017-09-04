#include <Rcpp.h>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

using namespace Rcpp;

/* Based on code from R  (R/src/cov.c:corcov and R/src/library/stats/src/cov.c:cov_complete2),
 with modifications:
 - only supporting two matrix version
 - only supporting Pearson correlation
 - not supporting NAs
 - assuming that there is long double
 - assuming that columns standard deviations are non-zero
 */
//' @title Pearson's correlation coefficient between columns of two matrices
//'     (OpenMP version).
//'
//' @description
//' \code{PPearsonMatOMP} returns a matrix of Pearson's correlation coefficients
//' between the columns of the two matrices passed as arguments.
//'
//' @details
//' This function calculates the Pearson's correlation coefficients between the
//' columns of two matrices given as argument. It is implemented in C++ for
//' efficiency. For two argument matrices with the same number of r rows and c1 ,
//' c2 columns, the return value is an c1-by-c2 matrix with the pairwise c1i c2j
//' correlation coefficients.
//'
//' @param A  Numeric matrix (variables by features).
//' @param B  Numeric matrix (variables by features).
//' @param nthreads Integer specifying the number of OpenMP threads to use in
//'     parallel parts (defaults to two, ignored on systems not supporting
//'     OpenMP).
//'
//' @return A matrix of dimensions \code{ncol(A})-by-\code{ncol(B)}.
//'
//' @examples
//' x <- matrix(1:12, nrow=3, ncol=4)
//' y <- matrix(1:15, nrow=3, ncol=5)
//' PPearsonMatOMP(x, y)
//'
// [[Rcpp::export]]
SEXP PPearsonMatOMP(SEXP A, SEXP B, int nthreads = 2) {
    A = PROTECT(Rf_coerceVector(A, REALSXP));
    B = PROTECT(Rf_coerceVector(B, REALSXP));
    
    int n = Rf_nrows(A), ncx = Rf_ncols(A), ncy = Rf_ncols(B);
    int i, j, k, n1 = n - 1;
    double *x, *y, *xm, *ym, *xx, *yy, *ans;
    long double sum, tmp, xxm, yym;
    SEXP ans_;
    
    if (n != Rf_nrows(B))
        Rf_error("incompatible dimensions");

    // get data pointers and allocate memory
    x = REAL(A);
    y = REAL(B);
    PROTECT(ans_ = Rf_allocMatrix(REALSXP, ncx, ncy));
    ans = REAL(ans_);
    xm = (double *) malloc(ncx * sizeof(double));
    ym = (double *) malloc(ncy * sizeof(double));

#ifdef _OPENMP
    // set the number of threads
    int nthreads_old = omp_get_num_threads();
    omp_set_num_threads(nthreads);
#endif
    
    // calculate column means -> xm[] and ym[] (two passes for better accuarcy)
#ifdef _OPENMP
#pragma omp parallel for private(sum, tmp, xx, i, k)
#endif
    for (i = 0 ; i < ncx ; i++) {
        xx = &x[i * n];
        sum = 0.;
        for (k = 0 ; k < n ; k++)
            sum += xx[k];
        tmp = sum / n;
        sum = 0.;
        for (k = 0 ; k < n ; k++)
            sum += (xx[k] - tmp);
        tmp = tmp + sum / n;
        xm [i] = (double)tmp;
    }
#ifdef _OPENMP
#pragma omp parallel for private(sum, tmp, yy, i, k)
#endif
    for (i = 0 ; i < ncy ; i++) {
        yy = &y[i * n];
        sum = 0.;
        for (k = 0 ; k < n ; k++)
            sum += yy[k];
        tmp = sum / n;
        sum = 0.;
        for (k = 0 ; k < n ; k++)
            sum += (yy[k] - tmp);
        tmp = tmp + sum / n;
        ym [i] = (double)tmp;
    }

    // calculate correlation coefficients
#define ANS(I,J)  ans[I + J * ncx]
#ifdef _OPENMP
#pragma omp parallel for private(sum, xx, yy, xxm, yym, i, j, k)
#endif
    for (i = 0 ; i < ncx ; i++) {
        xx = &x[i * n];
        xxm = xm[i];
        for (j = 0 ; j < ncy ; j++) {
            yy = &y[j * n];
            yym = ym[j];
            sum = 0.;
            for (k = 0 ; k < n ; k++)
                sum += (xx[k] - xxm) * (yy[k] - yym);
            ANS(i,j) = (double)(sum / n1);
        }
    }
    
    // calculate standard deviations --> xm[] and ym[]
#ifdef _OPENMP
#pragma omp parallel for private(sum, xx, xxm, i, k)
#endif
    for (i = 0 ; i < ncx ; i++) { /* sd(A[,i]) */
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
    for (i = 0 ; i < ncy ; i++) { /* sd(B[,i]) */
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
     
    /*   
    // set dimnames() - just copied from covcor; needs rewriting
    x = getAttrib(x, R_DimNamesSymbol);
    y = getAttrib(y, R_DimNamesSymbol);
    if ((length(x) >= 2 && !isNull(VECTOR_ELT(x, 1))) ||
        (length(y) >= 2 && !isNull(VECTOR_ELT(y, 1)))) {
        PROTECT(ind = allocVector(VECSXP, 2));
        if (length(x) >= 2 && !isNull(VECTOR_ELT(x, 1)))
            SET_VECTOR_ELT(ind, 0, duplicate(VECTOR_ELT(x, 1)));
        if (length(y) >= 2 && !isNull(VECTOR_ELT(y, 1)))
            SET_VECTOR_ELT(ind, 1, duplicate(VECTOR_ELT(y, 1)));
        setAttrib(ans, R_DimNamesSymbol, ind);
        UNPROTECT(1);
    }
     */

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
#undef CLAMP
