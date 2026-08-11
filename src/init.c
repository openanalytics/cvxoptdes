#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
void R_init_cvxoptdes(DllInfo *dll);
extern SEXP R_scs_solve(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dist_impl(SEXP, SEXP, SEXP);
extern SEXP delta_min_impl(SEXP, SEXP);
extern SEXP rank1_update(SEXP, SEXP, SEXP);
extern SEXP d_delta(SEXP, SEXP, SEXP);
extern SEXP a_delta(SEXP, SEXP, SEXP);
extern SEXP i_delta(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP g_delta(SEXP, SEXP, SEXP, SEXP);
extern SEXP alias_delta(SEXP, SEXP, SEXP, SEXP);
extern SEXP delta_eta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP crit_impl(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_highs_solve(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_highs_available(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"R_scs_solve", (DL_FUNC)&R_scs_solve, 14},
    {"dist_impl", (DL_FUNC)&dist_impl, 3},
    {"delta_min_impl", (DL_FUNC)&delta_min_impl, 2},
    {"rank1_update", (DL_FUNC)&rank1_update, 3},
    {"d_delta", (DL_FUNC)&d_delta, 3},
    {"a_delta", (DL_FUNC)&a_delta, 3},
    {"i_delta", (DL_FUNC)&i_delta, 5},
    {"g_delta", (DL_FUNC)&g_delta, 4},
    {"alias_delta", (DL_FUNC)&alias_delta, 4},
    {"delta_eta", (DL_FUNC)&delta_eta, 10},
    {"crit_impl", (DL_FUNC)&crit_impl, 5},
    {"R_highs_solve", (DL_FUNC)&R_highs_solve, 14},
    {"R_highs_available", (DL_FUNC)&R_highs_available, 1},
    {NULL, NULL, 0}};

void R_init_cvxoptdes(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
