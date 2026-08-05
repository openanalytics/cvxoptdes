#define R_NO_REMAP

#include <float.h>
#include "cvxoptdes.h"

static double d_crit(double *M, R_len_t p, Rboolean return_exp) 
{
        /* M = L' * L */
        double *L = (double *)R_alloc(p * p, sizeof(double));
        memcpy(L, M, p * p * sizeof(double));
        int info = -1;
        const char uplo = 'L'; /* store lower triangle */
        F77_CALL(dpotrf)(&uplo, &p, L, &p, &info FCONE);
        if (info != 0)
        {
            return return_exp ? 0.0: R_NegInf;
        }
        double crit = 0.0;
        for (R_len_t i = 0; i < p; i++)
        {
            double Lii = L[i + i * p] < DBL_EPSILON ? DBL_EPSILON : L[i + i * p];
            crit += 2.0 * log(Lii) / p; // |M|^{1/p}
        }
        return return_exp ? exp(crit) : crit;
}

static double i_crit(double *M, double *XtX, R_len_t p, R_len_t n, Rboolean return_inv)
{
    /* M = L' * L */
    double *L1 = (double *)R_alloc(p * p, sizeof(double));  // placeholder matrix
    double *L = (double *)R_alloc(p * p, sizeof(double));
    memcpy(L1, XtX, p * p * sizeof(double));
    memcpy(L, M, p * p * sizeof(double));
    int info = -1;
    double one = 1.0;
    const char uplo = 'L'; /* store lower triangle */
    const char sideL = 'L';
    const char sideR = 'R';
    const char diag = 'N';
    const char transN = 'N';
    const char transT = 'T';

    F77_CALL(dpotrf)(&uplo, &p, L, &p, &info FCONE);
    if (info != 0)
    {
        return return_inv ? 0.0 : R_PosInf;
    }
    /* L1 := L^{-1} * XtX */
    F77_CALL(dtrsm)(&sideL, &uplo, &transN, &diag, &p, &p, &one, L, &p, L1, &p FCONE FCONE FCONE FCONE);
    /* L1 := L^{-1} * XtX * L^{-1}' */
    F77_CALL(dtrsm)(&sideR, &uplo, &transT, &diag, &p, &p, &one, L, &p, L1, &p FCONE FCONE FCONE FCONE);
    /* tr(X * M^{-1} * X') = tr(L^{-1} * X'X * L^{-1}') */
    double crit = 0.0;
    for (R_len_t i = 0; i < p; i++)
    {
        crit += L1[i + i * p] / n;
    }
    return return_inv ? 1.0 / crit : crit;
}

static double a_crit(double *M, R_len_t p, Rboolean return_inv)
{
    /* M = L' * L */
    double *L = (double *)R_alloc(p * p, sizeof(double));
    memcpy(L, M, p * p * sizeof(double));
    int info = -1;
    const char uplo = 'L'; /* store lower triangle */
    const char diag = 'N';
    F77_CALL(dpotrf)(&uplo, &p, L, &p, &info FCONE);
    if (info != 0)
    {
        return return_inv ? 0.0 : R_PosInf;
    }
    /* L := L^{-1} */
    F77_CALL(dtrtri)(&uplo, &diag, &p, L, &p, &info FCONE FCONE);
    /* tr(M^{-1}) = ||L^{-1}||^2 */
    double crit = 0.0;
    for (R_len_t j = 0; j < p; j++)
    {
        for (R_len_t i = j; i < p; i++)
            crit += L[i + j * p] * L[i + j * p];
    }
    return return_inv ? 1.0 / crit : crit;
}

static double g_crit(double *M, double *X, R_len_t p, R_len_t n)
{
    /* M = L' * L */
    double *L1 = (double *)R_alloc(n * p, sizeof(double)); // placeholder matrix
    double *L = (double *)R_alloc(p * p, sizeof(double));
    memcpy(L1, X, n * p * sizeof(double));
    memcpy(L, M, p * p * sizeof(double));
    int info = -1;
    double one = 1.0;
    const char uplo = 'L'; /* store lower triangle */
    const char sideR = 'R';
    const char diag = 'N';
    const char transT = 'T';

    F77_CALL(dpotrf)(&uplo, &p, L, &p, &info FCONE);
    if (info != 0)
    {
        return 0.0;
    }
    /* L1 := X * L^{-1}' */
    F77_CALL(dtrsm)(&sideR, &uplo, &transT, &diag, &n, &p, &one, L, &p, L1, &n FCONE FCONE FCONE FCONE);

    /* diag(X * M^{-1} * X') = ||(X * L{-1}')_i||^2 */
    double crit = 0.0;
    for (int i = 0; i < n; i++)
    {
        double s = 0.0;
        for (int k = 0; k < p; k++)
        {
            double v = L1[i + k * n];
            s += v * v;
        }
        if(s > crit)
            crit = s;
    }
    return 1.0 / crit;
}

static double alias_crit(double *M, double *ortho, R_len_t p, Rboolean return_one_min)
{
    double dj, di, r, r_rho;
    double sum = 0.0;
    R_len_t n = 0;

    for (R_len_t j = 0; j < p; j++)
    {
        if (M[j + j * p] < DBL_EPSILON)
            return R_PosInf;

        dj = sqrt(M[j + j * p]);
        for (R_len_t i = j + 1; i < p; i++)
        {
            if(!ortho)
            {
                di = sqrt(M[i + i * p]);
                r = M[i + j * p] / (di * dj);
                sum += fabs(r);
            }
            else if (ortho[i + j * p] != R_PosInf)
            {
                    di = sqrt(M[i + i * p]);
                    r = M[i + j * p] / (di * dj);
                    r_rho = fabs(r) - ortho[i + j * p];
                    if (r_rho > 0.0)
                        sum += r_rho;
            }
            n++;
        }
    }
    return return_one_min ? 1.0 - sum / n : sum / n;
}

double crit_internal(double *M, double *X, double *XtX, R_len_t p, R_len_t n, int crit, Rboolean return_trans)
{
    switch (crit)
    {
        case 0:
            return a_crit(M, p, return_trans);
        case 1:
            return i_crit(M, XtX, p, n, return_trans);
        case 2:
            return d_crit(M, p, return_trans);
        case 3:
            return g_crit(M, X, p, n);
        case 4:
            return alias_crit(M, XtX, p, return_trans);
        default:
            return R_NegInf;
    }
}

/*
* @brief Return D-, I-, A-, G-, alias-optimality criterion 
*
* @param[in] M_in R matrix, information matrix
* @param[in] X_in R matrix, design matrix. Only used for G-optimality
* @param[in] XtX_in R matrix, product X'X. Only used for I-optimality
* @param[in] crit_in R integer, 0 = A-, 1 = I-, 2 = D-, 3 = G-, 4 = alias-criterion
* @param[in] return_trans R boolean, if TRUE returns transformed (ordinary) criterion 
* for A-, I- and D-optimality. If FALSE, returns inverse criterion for A- and I-optimality, 
* log-criterion for D-optimality, and 1 minus criterion for alias-optimality.
* @return numeric optimality criterion
*/
SEXP crit_impl(SEXP M_in, SEXP X_in, SEXP XtX_in, SEXP crit_in, SEXP return_trans)
{
    int crit = INTEGER_ELT(crit_in, 0);
    R_len_t p = Rf_nrows(M_in);
    R_len_t n = 1;
    Rboolean rt = FALSE;
    double *M = REAL(M_in);
    double *X = NULL;
    double *XtX = NULL;

    if(crit == 0)
    {
        rt = LOGICAL_ELT(return_trans, 0);
    }
    else if(crit == 1)
    {
        XtX = REAL(XtX_in);
        n = Rf_nrows(X_in);
        rt = LOGICAL_ELT(return_trans, 0);
    }
    else if(crit == 2)
    {
        rt = LOGICAL_ELT(return_trans, 0);
    }
    else if(crit == 3)
    {
        X = REAL(X_in);
        n = Rf_nrows(X_in);
    } 
    else if(crit == 4)
    {
        // placeholder for ortho-constraint matrix
        if(!Rf_isNull(XtX_in))
            XtX = REAL(XtX_in);
        rt = LOGICAL_ELT(return_trans, 0);
    }

    double crit_val = crit_internal(M, X, XtX, p, n, crit, rt);
    return Rf_ScalarReal(crit_val);
}
