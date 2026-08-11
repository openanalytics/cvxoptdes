#ifndef CVXOPTDES_H
#define CVXOPTDES_H

#ifndef HIGHS_AVAILABLE
#define HIGHS_AVAILABLE 0
#endif

#ifndef DEBUG_INFO
#define DEBUG_INFO 0
#endif

#ifndef USE_SPECTRAL_CONES
#define USE_SPECTRAL_CONES 0
#endif

#include <math.h>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#include "scs/include/glbopts.h"
#include "scs/include/scs.h"
#include "scs/include/util.h"
#include "scs/linsys/scs_matrix.h"

typedef struct
{
    SEXP J;        // design or Jacobian matrix
    SEXP f;        // cost vector
    SEXP ortho;    // (optional) orthogonality constraint matrix
    SEXP strata;   // (optional) strata + weights
    SEXP zeros;    // (optional) zero weight constraint vector
    SEXP Xzero;    // (optional) zero label binary matrix (same size as J)
    SEXP colmeans; // (optional) column mean constraint vector
    SEXP env;      // function environment
    SEXP w0;       // (optional) vector of initial design weights
    SEXP K;        // K matrix for A-, I-optimality
    SEXP v_sca;    // (optional) v matrix needed for mixed-model SCA
    SEXP start;    // (optional) list with x, y, s vectors used for warm-start
    SEXP ctrl_int; // integer control paramaters
    SEXP ctrl_dbl; // double control paramaters
    ScsCone *k;
    ScsData *d;
    ScsSettings *stgs;
    ScsSolution *sol;
    ScsInfo *info;
    ScsWork *w;
} pdata;

// scs_optimize.c
SEXP R_scs_solve(SEXP J, SEXP f, SEXP ortho, SEXP strata, SEXP zeros, SEXP Xzero, SEXP colmeans, SEXP env, SEXP w0, SEXP K, SEXP v_sca, SEXP start, SEXP ctrl_int, SEXP ctrl_dbl);

// criteria.c
void D_opt_populate_Abc(pdata *pars, R_len_t p, R_len_t n, R_len_t nnz, scs_int mtot, scs_int ntot, double alpha, double lambda, double gamma, double upr, scs_int nlvls, R_len_t z_ortho, R_len_t q_ortho, R_len_t mzeros, R_len_t mmeans);
void D_opt_populate_K(pdata *pars, R_len_t p, R_len_t n, double upr, scs_int nlvls, R_len_t z_ortho, R_len_t qsize, R_len_t mzeros, R_len_t mmeans, double lambda, R_len_t blocks);
void D_spec_opt_populate_Abc(pdata *pars, R_len_t p, R_len_t n, R_len_t nnz, scs_int mtot, scs_int ntot, double alpha, double lambda, double gamma, double upr, scs_int nlvls, R_len_t z_ortho, R_len_t q_ortho, R_len_t mzeros, R_len_t mmeans);
void D_spec_opt_populate_K(pdata *pars, R_len_t p, R_len_t n, double upr, scs_int nlvls, R_len_t z_ortho, R_len_t qsize, R_len_t mzeros, R_len_t mmeans);
void A_opt_populate_Abc(pdata *pars, R_len_t p, R_len_t n, R_len_t nnz, scs_int mtot, scs_int ntot, double alpha, double lambda, double gamma, double upr, scs_int nlvls, R_len_t z_ortho, R_len_t q_ortho, R_len_t mzeros, R_len_t mmeans);
void A_opt_populate_K(pdata *pars, R_len_t p, R_len_t n, double upr, scs_int nlvls, R_len_t z_ortho, R_len_t q_ortho, R_len_t qsize, R_len_t mzeros, R_len_t mmeans, double lambda, R_len_t blocks);
void G_opt_populate_Abc(pdata *pars, R_len_t p, R_len_t n, R_len_t nnz, scs_int mtot, scs_int ntot, double alpha, double lambda, double gamma, double upr, scs_int nlvls, R_len_t z_ortho, R_len_t q_ortho, R_len_t mzeros, R_len_t mmeans);
void G_opt_populate_K(pdata *pars, R_len_t p, R_len_t n, double upr, scs_int nlvls, R_len_t z_ortho, R_len_t qsize, R_len_t mzeros, R_len_t mmeans, double lambda, R_len_t blocks);
void alias_opt_populate_Abc(pdata *pars, R_len_t p, R_len_t n, R_len_t nnz, scs_int mtot, scs_int ntot, double alpha, double gamma, double upr, scs_int nlvls, R_len_t z_ortho, R_len_t q_ortho, R_len_t mzeros, R_len_t mmeans);
void alias_opt_populate_K(pdata *pars, R_len_t p, R_len_t n, double upr, scs_int nlvls, R_len_t z_ortho, R_len_t qsize, R_len_t mzeros, R_len_t mmeans);

// crit.c
double crit_internal(double *M, double *X, double *XtX, R_len_t p, R_len_t n, int crit, Rboolean return_trans);

// exchange.c
SEXP dist_impl(SEXP X1_in, SEXP X2_in, SEXP X1_idx);
SEXP delta_min_impl(SEXP X_in, SEXP tol);
SEXP rank1_update(SEXP Minv, SEXP u, SEXP v);
SEXP d_delta(SEXP wloc, SEXP X_in, SEXP Minv_in);
SEXP a_delta(SEXP wloc, SEXP X_in, SEXP Minv_in);
SEXP i_delta(SEXP wloc, SEXP X_in, SEXP XtX_in, SEXP M_in, SEXP Minv_in);
SEXP g_delta(SEXP wloc, SEXP X_in, SEXP M_in, SEXP Minv_in);
SEXP alias_delta(SEXP wloc, SEXP X_in, SEXP M_in, SEXP ortho_in);
SEXP delta_eta(SEXP wloc, SEXP X_in, SEXP A_in, SEXP iseq, SEXP crit, SEXP crit0, SEXP Zind, SEXP Vli, SEXP Xi, SEXP Ai);
SEXP crit_impl(SEXP M_in, SEXP X_in, SEXP XtX_in, SEXP crit_in, SEXP return_trans);

// highs.c
SEXP R_highs_solve(SEXP L, SEXP lwr, SEXP upr, SEXP A, SEXP lhs, SEXP rhs, SEXP types, SEXP maximum, SEXP offset, SEXP ctrl_dbl, SEXP ctrl_int, SEXP ctrl_bool, SEXP ctrl_str, SEXP write_mps);
SEXP R_highs_available(SEXP null);

#endif
