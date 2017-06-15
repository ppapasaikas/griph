#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* 
 registering of native routines is now (R 3.4) required
 by R CMD check.
 
 the content of this file was auto-generated by calling:
 tools::package_native_routine_registration_skeleton(".")
 in the package root directory
 
 furthermore, the following was added to NAMESPACE:
 useDynLib(griph, .registration = TRUE)
*/

/* .Call calls */
extern SEXP griph_checkBits();
extern SEXP griph_checkOpenMP();
extern SEXP griph_fastCDistance(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP griph_fastDistance(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP griph_fastSDistance(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP griph_PCanberraMat(SEXP, SEXP);
extern SEXP griph_PCanberraMatOMP(SEXP, SEXP, SEXP);
extern SEXP griph_PHellingerMat(SEXP, SEXP);
extern SEXP griph_PHellingerMatOMP(SEXP, SEXP, SEXP);
extern SEXP griph_PPearsonMatOMP(SEXP, SEXP, SEXP);
extern SEXP griph_referenceWij(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP griph_searchTrees(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP griph_searchTreesCSparse(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP griph_searchTreesTSparse(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP griph_sgd(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"griph_checkBits",          (DL_FUNC) &griph_checkBits,           0},
    {"griph_checkOpenMP",        (DL_FUNC) &griph_checkOpenMP,         0},
    {"griph_fastCDistance",      (DL_FUNC) &griph_fastCDistance,       8},
    {"griph_fastDistance",       (DL_FUNC) &griph_fastDistance,        6},
    {"griph_fastSDistance",      (DL_FUNC) &griph_fastSDistance,       8},
    {"griph_PCanberraMat",       (DL_FUNC) &griph_PCanberraMat,        2},
    {"griph_PCanberraMatOMP",    (DL_FUNC) &griph_PCanberraMatOMP,     3},
    {"griph_PHellingerMat",      (DL_FUNC) &griph_PHellingerMat,       2},
    {"griph_PHellingerMatOMP",   (DL_FUNC) &griph_PHellingerMatOMP,    3},
    {"griph_PPearsonMatOMP",     (DL_FUNC) &griph_PPearsonMatOMP,      3},
    {"griph_referenceWij",       (DL_FUNC) &griph_referenceWij,        5},
    {"griph_searchTrees",        (DL_FUNC) &griph_searchTrees,         9},
    {"griph_searchTreesCSparse", (DL_FUNC) &griph_searchTreesCSparse, 11},
    {"griph_searchTreesTSparse", (DL_FUNC) &griph_searchTreesTSparse, 11},
    {"griph_sgd",                (DL_FUNC) &griph_sgd,                15},
    {NULL, NULL, 0}
};

void R_init_griph(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
