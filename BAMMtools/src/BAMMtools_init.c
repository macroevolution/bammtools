#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

// Created with tools::package_native_routine_registration_skeleton(".")

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void cohort_matrix(void *, void *, void *, void *);
extern void fetchmrca(void *, void *, void *, void *, void *, void *, void *, void *);
extern void jenksBrks(void *, void *, void *, void *);
extern void postorder_tree_traverse(void *, void *, void *, void *, void *);
extern void setphylotreecoords(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void setpolartreecoords(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void setrecursivesequence(void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP dtrates(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP seq_root2tip(SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"cohort_matrix",           (DL_FUNC) &cohort_matrix,            4},
    {"fetchmrca",               (DL_FUNC) &fetchmrca,                8},
    {"jenksBrks",               (DL_FUNC) &jenksBrks,                4},
    {"postorder_tree_traverse", (DL_FUNC) &postorder_tree_traverse,  5},
    {"setphylotreecoords",      (DL_FUNC) &setphylotreecoords,      11},
    {"setpolartreecoords",      (DL_FUNC) &setpolartreecoords,       9},
    {"setrecursivesequence",    (DL_FUNC) &setrecursivesequence,     6},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"dtrates",      (DL_FUNC) &dtrates,      5},
    {"seq_root2tip", (DL_FUNC) &seq_root2tip, 3},
    {NULL, NULL, 0}
};

void R_init_BAMMtools(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
