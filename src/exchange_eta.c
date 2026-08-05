#define R_NO_REMAP

#include "cvxoptdes.h"

/**
 * @brief Delta matrix (exchanging 1 design point) for mixed-effects models
 * with block-diagonal information matrix
 *
 * @param[in] wloc R integer vector, indices of design points
 * @param[in] X_in R matrix, fixed-effects design matrix
 * @param[in] XtX_in R
 * @param[in] crit R integer, optimality criterion
 * @param[in] crit0 R numeric, current value of optimality criterion
 * @param[in] Zind R vector, column indicators of Z-matrix
 * @param[in] iseq R vector,
 * @param[in] Vli R list
 * @param[in] Xi R list
 * @param[in] Ai R list
 * @param[in] A R matrix
 * @return delta matrix (length(wloc) x n)
 */
SEXP delta_eta(SEXP wloc, SEXP X_in, SEXP A_in, SEXP iseq, SEXP crit, SEXP crit0, SEXP Zind, SEXP Vli, SEXP Xi, SEXP Ai)
{
    int *iwloc = INTEGER(wloc);
    int *iiseq = INTEGER(iseq);
    int icrit = INTEGER_ELT(crit, 0);
    int m = Rf_length(wloc);

    int *dims = INTEGER(Rf_getAttrib(X_in, R_DimSymbol));
    int n = dims[0];
    int p = dims[1];
    int pp = p * p;
    const char transN = 'N';
    const char transT = 'T';
    const char sideL = 'L';
    const char uplo = 'U';
    const double one = 1.0, zero = 0.0, negone = -1.0;
    const int inc = 1;

    double *X = REAL(X_in);
    double *A = REAL(A_in);
    double *XtX = NULL;

    /* workspaces */
    SEXP delta_sxp = PROTECT(Rf_allocMatrix(REALSXP, m, n));
    double *delta = REAL(delta_sxp);
    double *XVX = (double *)R_alloc(p * p, sizeof(double));
    double *A_tmp = (double *)R_alloc(p * p, sizeof(double));

    /* X'X only needed for I-optimality */
    if(icrit == 1)
    {
        XtX = (double *)S_alloc(p * p, sizeof(double));
        F77_CALL(dsyrk)(&uplo, &transT, &p, &n, &one, X, &n, &zero, XtX, &p FCONE FCONE);
        for (int j = 0; j < p; j++)
            for (int i = j + 1; i < p; i++)
                XtX[i + j * p] = XtX[j + i * p];
    }

    for (int i = 0; i < n; i++)
    {
        R_len_t gi = INTEGER_ELT(Zind, i) - 1;
        for (int j = 0; j < m; j++)
        {
            if(i == (iwloc[j] - 1))
            {
                delta[j + i * m] = REAL_ELT(crit0, 0);
            }
            else 
            {
                R_len_t gj = INTEGER_ELT(Zind, iwloc[j] - 1) - 1;
                SEXP A1i_sxp = PROTECT(VECTOR_ELT(Ai, gi));
                SEXP X1i_sxp = PROTECT(Rf_duplicate(VECTOR_ELT(Xi, gi)));
                double *A1i = NULL;
                double *X1i = NULL;
                int nx = 0;
                if (!Rf_isNull(X1i_sxp))
                {
                    X1i = REAL(X1i_sxp);
                    nx = Rf_nrows(X1i_sxp);
                }
                if (!Rf_isNull(A1i_sxp))
                {
                    A1i = REAL(A1i_sxp);
                }
                if(gi == gj)
                {
                    SEXP V1_sxp = PROTECT(VECTOR_ELT(Vli, 3 * gi));
                    double *V1 = REAL(V1_sxp);
                    int nv = Rf_nrows(V1_sxp);
                    for (int k = 0; k < p; k++)
                        X1i[(iiseq[j] - 1) + k * nx] = X[i + k * n];

                    /* VX = V1 * X1 */
                    double *VX = (double *)S_alloc(nv * p, sizeof(double));
                    F77_CALL(dsymm)(&sideL, &uplo, &nv, &p, &one, V1, &nv, X1i, &nx, &zero, VX, &nv FCONE FCONE);

                    /* XVX = X1' * V1 * X1 */
                    memset(XVX, 0, p * p * sizeof(double));
                    F77_CALL(dgemm)(&transT, &transN, &p, &p, &nx, &one, X1i, &nx, VX, &nv, &zero, XVX, &p FCONE FCONE);

                    /* A := A - A[gi] + XVX */
                    memcpy(A_tmp, A, p * p * sizeof(double));
                    if(!Rf_isNull(A1i_sxp)) 
                    {
                        F77_CALL(daxpy)(&pp, &negone, A1i, &inc, A_tmp, &inc);
                    }
                    F77_CALL(daxpy)(&pp, &one, XVX, &inc, A_tmp, &inc);

                    UNPROTECT(1);
                }
                else
                {
                    SEXP V1_sxp = PROTECT(VECTOR_ELT(Vli, 3 * gi + 1));
                    SEXP A1j_sxp = PROTECT(VECTOR_ELT(Ai, gj));
                    double *V1 = REAL(V1_sxp);
                    int nv1 = Rf_nrows(V1_sxp);  // nv + 1
                    int nx1 = nx + 1;
                    double *X1i_new = (double *)R_alloc(nx1 * p, sizeof(double));
                    for (int k2 = 0; k2 < p; k2++)
                    {
                        if(nx > 0)
                        {
                            for (int k1 = 0; k1 < nx; k1++)
                                X1i_new[k1 + k2 * nx1] = X1i[k1 + k2 * nx];
                        }
                        X1i_new[nx + k2 * nx1] = X[i + k2 * n];
                    }

                    /* VX = V1 * X1 */
                    double *VX1 = (double *)S_alloc(nv1 * p, sizeof(double));
                    F77_CALL(dsymm)(&sideL, &uplo, &nv1, &p, &one, V1, &nv1, X1i_new, &nx1, &zero, VX1, &nv1 FCONE FCONE);

                    /* XVX = X1' * V1 * X1 */
                    memset(XVX, 0, p * p * sizeof(double));
                    F77_CALL(dgemm)(&transT, &transN, &p, &p, &nx1, &one, X1i_new, &nx1, VX1, &nv1, &zero, XVX, &p FCONE FCONE);

                    /* A := A - A[gi] - A[gj] + XVX */
                    memcpy(A_tmp, A, p * p * sizeof(double));
                    if(!Rf_isNull(A1i_sxp))
                    {
                        F77_CALL(daxpy)(&pp, &negone, A1i, &inc, A_tmp, &inc);
                    }
                    if (!Rf_isNull(A1j_sxp))
                    {
                        double *A1j = REAL(A1j_sxp);
                        F77_CALL(daxpy)(&pp, &negone, A1j, &inc, A_tmp, &inc);
                    }
                    F77_CALL(daxpy)(&pp, &one, XVX, &inc, A_tmp, &inc);

                    SEXP V2_sxp = PROTECT(VECTOR_ELT(Vli, 3 * gj + 2));
                    if(!Rf_isNull(V2_sxp))
                    {
                        SEXP X1j_sxp = PROTECT(Rf_duplicate(VECTOR_ELT(Xi, gj)));
                        double *X1j = REAL(X1j_sxp);
                        double *V2 = REAL(V2_sxp);
                        int nv2 = Rf_nrows(V2_sxp); // nv - 1
                        int nx2 = Rf_nrows(X1j_sxp) - 1;
                        double *X1j_new = (double *)R_alloc(nx2 * p, sizeof(double));
                        R_len_t irow = 0;
                        for (int k1 = 0; k1 < (nx2 + 1); k1++)
                        {
                            if(k1 != (iiseq[j] - 1)) 
                            {
                                for (int k2 = 0; k2 < p; k2++)
                                    X1j_new[irow + k2 * nx2] = X1j[k1 + k2 * (nx2 + 1)];
                                irow++;
                            }
                        }

                        /* VX = V2 * X2 */
                        double *VX2 = (double *)S_alloc(nv2 * p, sizeof(double));
                        F77_CALL(dsymm)(&sideL, &uplo, &nv2, &p, &one, V2, &nv2, X1j_new, &nx2, &zero, VX2, &nv2 FCONE FCONE);

                        /* XVX = X2' * V2 * X2 */
                        memset(XVX, 0, p * p * sizeof(double));
                        F77_CALL(dgemm)(&transT, &transN, &p, &p, &nx2, &one, X1j_new, &nx2, VX2, &nv2, &zero, XVX, &p FCONE FCONE);

                        /* A := A - A[gi] - A[gj] + X2' * V2 * X2 */
                        F77_CALL(daxpy)(&pp, &one, XVX, &inc, A_tmp, &inc);

                        UNPROTECT(1);
                    }
                    UNPROTECT(3);
                }
                // criterion
                delta[j + i * m] = crit_internal(A_tmp, X, XtX, p, n, icrit, FALSE);
                // rescale I-/D-criteria
                if(icrit == 1)
                    delta[j + i * m] *= n;
                else if(icrit == 2)
                    delta[j + i * m] *= p;
                
                UNPROTECT(2);
            }
        }
    }
    UNPROTECT(1);
    return delta_sxp;
}
