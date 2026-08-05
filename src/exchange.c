#define R_NO_REMAP

#include <float.h>
#include "cvxoptdes.h"

static inline double max_dbl(const double *x, R_len_t n)
{
    double maxv = x[0];
    for (R_len_t i = 1; i < n; ++i)
        if (x[i] > maxv)
            maxv = x[i];
    return maxv;
}

static double mean_abs_corr(const double *M, R_len_t p)
{
    double dj, di, r;
    double sum = 0.0;
    R_len_t n = 0;

    for (R_len_t j = 0; j < p; j++)
    {
        if(M[j + j * p] < DBL_EPSILON)
            return R_PosInf;    

        dj = sqrt(M[j + j * p]);
        for (R_len_t i = j + 1; i < p; i++)
        {
            di = sqrt(M[i + i * p]);
            r = M[i + j * p] / (di * dj);
            sum += fabs(r);
            n++;
        }
    }
    return sum / n;
}

static double mean_abs_corr_ortho(const double *M, const double *ortho, R_len_t p)
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
            if(ortho[i + j * p] != R_PosInf)
            {
                di = sqrt(M[i + i * p]);
                r = M[i + j * p] / (di * dj);
                r_rho = fabs(r) - ortho[i + j * p];
                if(r_rho > 0.0)
                    sum += r_rho;
            }
            n++;
        }
    }
    return sum / n;
}

// /**
//  * @brief Multiply design matrix X by inverse information matrix M^{-1}
//  *
//  * @param[in] M R matrix, information matrix
//  * @param[in] X R matrix, design matrix
//  * @param[in] transpose int, whether to transpose the output matrix
//  * @return Y = X * M^{-1} if transpose == 0, else Y = (X * M^{-1})'
//  */
// static SEXP X_Minv(SEXP M, SEXP X, int transpose)
// {
//     int *dims = INTEGER(Rf_getAttrib(X, R_DimSymbol));
//     int n = dims[0];
//     int p = dims[1];
//     /* M = L' * L */
//     SEXP L_sxp = PROTECT(Rf_duplicate(M));
//     double *L = REAL(L_sxp);
//     int info = -1;
//     const char uplo = 'U'; /* store upper triangle */
//     F77_CALL(dpotrf)(&uplo, &p, L, &p, &info FCONE);
//     if (info != 0)
//     {
//         UNPROTECT(1);
//         return R_NilValue;
//     }
//     /* solve M * Y' = X'  <=>  Y = X * M^{-1} */
//     SEXP Yt_sxp = PROTECT(Rf_allocMatrix(REALSXP, p, n));
//     double *Yt = REAL(Yt_sxp);
//     double *X_dbl = REAL(X);
//     for (R_len_t j = 0; j < p; j++)
//         for (R_len_t i = 0; i < n; i++)
//             Yt[j + p * i] = X_dbl[i + n * j];

//     F77_CALL(dpotrs)(&uplo, &p, &n, L, &p, Yt, &p, &info FCONE);
//     if (info != 0)
//     {
//         UNPROTECT(2);
//         return R_NilValue;
//     }
//     /* transpose Y' */
//     if (!transpose)
//     {
//         UNPROTECT(2);
//         return Yt_sxp;
//     }
//     else
//     {
//         SEXP Y_sxp = PROTECT(Rf_allocMatrix(REALSXP, n, p));
//         double *Y = REAL(Y_sxp);
//         for (R_len_t j = 0; j < p; j++)
//             for (R_len_t i = 0; i < n; i++)
//                 Y[i + n * j] = Yt[j + p * i];
//         UNPROTECT(3);
//         return Y_sxp;
//     }
// }

/**
 * @brief Compute Euclidean distance matrix between two sets of points,
 * where the first set of points is a subset of the second set.
 *
 * @param[in] X1_in R matrix, first set of points (n1 x p)
 * @param[in] X2_in R matrix, second set of points (n2 x p)
 * @param [in] X1_idx R integer vector, row locations of X1 in X2
 * @return distance matrix D (n1 x n2), D[i,j] = ||X1[i,] - X2[j,]||
 */
SEXP dist_impl(SEXP X1_in, SEXP X2_in, SEXP X1_idx)
{
    int *dims1 = INTEGER(Rf_getAttrib(X1_in, R_DimSymbol));
    int n1 = dims1[0];
    int p = dims1[1];

    int *dims2 = INTEGER(Rf_getAttrib(X2_in, R_DimSymbol));
    int n2 = dims2[0];

    double *X1 = REAL(X1_in);
    double *X2 = REAL(X2_in);

    SEXP dist_sxp = PROTECT(Rf_allocMatrix(REALSXP, n1, n2));
    double *dist = REAL(dist_sxp);
    double diff, sum_ij;
    R_len_t idx;

    for (R_len_t i = 0; i < n1; i++)
    {
        idx = INTEGER_ELT(X1_idx, i) - 1;
        for (R_len_t j = 0; j < n2; j++)
        {
            if(idx != j)
            {
                sum_ij = 0.0;
                for (R_len_t k = 0; k < p; k++)
                {
                    diff = X1[i + n1 * k] - X2[j + n2 * k];
                    sum_ij += diff * diff;
                }
                dist[i + n1 * j] = sqrt(sum_ij);
            }
            else 
            {
                dist[i + n1 * j] = 0.0;
            }
        }
    }

    UNPROTECT(1);
    return dist_sxp;
}

/**
 * @brief Compute column-wise minimum distances replacing minima
 * with the second minimum if the minimum occurs more than once.
 * @param[in] X_in R matrix, pre-computed distance matrix (n1 x n2)
 * @param[in] tol R numeric, tolerance added to minimum distances (to avoid early stopping)
 * @return R matrix of same dimension as X_in
 */
SEXP delta_min_impl(SEXP X_in, SEXP tol)
{
    int *dims = INTEGER(Rf_getAttrib(X_in, R_DimSymbol));
    int n1 = dims[0];
    int n2 = dims[1];
    double tol_dbl = REAL_ELT(tol, 0);

    SEXP delta_sxp = PROTECT(Rf_allocMatrix(REALSXP, n1, n2));
    SEXP is_min_sxp = PROTECT(Rf_allocVector(REALSXP, n1));

    double *X = REAL(X_in);
    double *delta = REAL(delta_sxp);
    double *is_min = REAL(is_min_sxp);

    double min1, min2;
    R_len_t count_min1;

    for (int j = 0; j < n2; j++)
    {
        /* first pass: find minimum */
        min1 = R_PosInf;
        for (R_len_t i = 0; i < n1; i++)
        {
            double v = X[i + j * n1];
            if (v != R_PosInf && v < min1)
                min1 = v;
        }
        /* second pass: number of occurrences of minimum + second minimum */
        count_min1 = 0;
        min2 = R_PosInf;
        for (R_len_t i = 0; i < n1; i++)
        {
            double v = X[i + j * n1];
            if (v <= min1 + 1e-8)
            {
                is_min[i] = 1.0;
                count_min1++;
            }
            else 
            {
                is_min[i] = 0.0;
                if (min2 == R_PosInf || v < min2)
                    min2 = v;
            }
        }
        /* third pass: populate delta */
        if (count_min1 > 1)
        {
            for (R_len_t i = 0; i < n1; i++)
                delta[i + j * n1] = min1 + is_min[i] * 1.1 * tol_dbl;
        }
        else
        {
            for (R_len_t i = 0; i < n1; i++)
                delta[i + j * n1] = is_min[i] ? min2 : min1;
        }
    }

    UNPROTECT(2);
    return delta_sxp;
}

/**
 * @brief Rank-1 update of inverse information matrix
 *
 * @param[in] Minv R matrix, initial inverse information matrix
 * @param[in] u R vector, design point to add/remove
 * @param[in] sgn numeric value, +1 to add design point, -1 to remove design point
 * @return updated inverse information matrix of same dimension as Minv
 */
SEXP rank1_update(SEXP Minv, SEXP u, SEXP sgn)
{
    int *dims = INTEGER(Rf_getAttrib(Minv, R_DimSymbol));
    int n = dims[0];
    double sign = Rf_asReal(sgn);

    double *Minv_dbl = REAL(Minv);
    const double *u_dbl = REAL(u);

    /* v = Minv * u */
    SEXP v = PROTECT(Rf_allocVector(REALSXP, n));
    double *v_dbl = REAL(v);
    const char trans = 'N';
    const double one = 1.0, zero = 0.0;
    const int ione = 1;
    F77_CALL(dgemv)(&trans, &n, &n, &one, Minv_dbl, &n,
                    u_dbl, &ione, &zero, v_dbl, &ione FCONE);

    /* uMu = u' * v */
    double uMu = F77_CALL(ddot)(&n, u_dbl, &ione, v_dbl, &ione);
    double denom = 1.0 + sign * uMu;
    if (fabs(denom) < DBL_EPSILON)
    {
        UNPROTECT(1);
        return R_NilValue;
    }

    /* MuuM = Minv * u * u' * Minv */
    /* Minv := Minv - sgn * MuuM / (1 + sgn * uMu) */
    SEXP Minv1 = PROTECT(Rf_allocMatrix(REALSXP, n, n));
    double *Minv1_dbl = REAL(Minv1);
    for (R_len_t i = 0; i < n; i++)
    {
        for (R_len_t j = 0; j < n; j++)
        {
            Minv1_dbl[i * n + j] = Minv_dbl[i * n + j];
        }
    }
    double alpha = -sign / denom;
    F77_CALL(dger)(&n, &n, &alpha, v_dbl, &ione, v_dbl, &ione, Minv1_dbl, &n);

    UNPROTECT(2);
    return Minv1;
}

/**
 * @brief D-optimality delta matrix (exchanging 1 design point)
 *
 * @param[in] wloc R integer vector, indices of design points
 * @param[in] X_in R matrix, design matrix
 * @param[in] Minv_in R matrix, inverse information matrix
 * @return delta matrix (length(wloc) x n)
 */
SEXP d_delta(SEXP wloc, SEXP X_in, SEXP Minv_in)
{
    int *iwloc = INTEGER(wloc);
    int m = Rf_length(wloc);

    int *dims = INTEGER(Rf_getAttrib(X_in, R_DimSymbol));
    int n = dims[0];
    int p = dims[1];

    double *X = REAL(X_in);
    double *Minv = REAL(Minv_in);

    /* XM = X * Minv */
    SEXP XM_sxp = PROTECT(Rf_allocMatrix(REALSXP, n, p));
    double *XM = REAL(XM_sxp);
    const char transN = 'N';
    const char transT = 'T';
    const double one = 1.0, zero = 0.0;

    F77_CALL(dgemm)(&transN, &transN, &n, &p, &p, &one, X, &n, Minv, &p, &zero, XM, &n FCONE FCONE);

    /* uMu = sum_j X[i,j] * XM[i,j] */
    SEXP uMu_sxp = PROTECT(Rf_allocVector(REALSXP, n));
    double *uMu = REAL(uMu_sxp);

    for (int i = 0; i < n; i++)
    {
        double sum_i = 0.0;
        for (int j = 0; j < p; j++)
            sum_i += X[i + n * j] * XM[i + n * j];
        uMu[i] = sum_i;
    }

    /* vMv = uMu[wloc] */
    SEXP vMv_sxp = PROTECT(Rf_allocVector(REALSXP, m));
    double *vMv = REAL(vMv_sxp);
    for (int i = 0; i < m; i++)
        vMv[i] = uMu[iwloc[i] - 1];

    /* extract X[wloc,] */
    SEXP Xw_sxp = PROTECT(Rf_allocMatrix(REALSXP, m, p));
    double *Xw = REAL(Xw_sxp);

    for (int i = 0; i < m; i++)
        for (int j = 0; j < p; j++)
            Xw[i + m * j] = X[(iwloc[i] - 1) + n * j];

    /* XwM = X[wloc,] * Minv */
    SEXP XwM_sxp = PROTECT(Rf_allocMatrix(REALSXP, m, p));
    double *XwM = REAL(XwM_sxp);

    F77_CALL(dgemm)(&transN, &transN, &m, &p, &p, &one, Xw, &m, Minv, &p, &zero, XwM, &m FCONE FCONE);

    /* vMu = X[wloc,] * Minv * t(X) -> m x n */
    SEXP vMu_sxp = PROTECT(Rf_allocMatrix(REALSXP, m, n));
    double *vMu = REAL(vMu_sxp);

    F77_CALL(dgemm)(&transN, &transT, &m, &n, &p, &one, XwM, &m, X, &n, &zero, vMu, &m FCONE FCONE);

    /* delta function : add xi, remove xj
       delta(xi, xj) = d(xi) - d(xj) - d(xi) d(xj) + d(xj, xi) ^ 2 */
    SEXP delta_sxp = PROTECT(Rf_allocMatrix(REALSXP, m, n));
    double *delta = REAL(delta_sxp);

    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            delta[i + m * j] = uMu[j] - (vMv[i] + vMv[i] * uMu[j]) + vMu[i + m * j] * vMu[i + m * j];

    UNPROTECT(7);
    return delta_sxp;
}

/**
 * @brief A-optimality delta matrix (exchanging 1 design point)
 *
 * @param[in] wloc R integer vector, indices of design points
 * @param[in] X_in R matrix, design matrix
 * @param[in] Minv_in R matrix, inverse information matrix
 * @return delta matrix (length(wloc) x n)
 */
SEXP a_delta(SEXP wloc, SEXP X_in, SEXP Minv_in)
{
    int *iwloc = INTEGER(wloc);
    int m = Rf_length(wloc);

    int *dims = INTEGER(Rf_getAttrib(X_in, R_DimSymbol));
    int n = dims[0];
    int p = dims[1];
    const char transN = 'N';
    const char transT = 'T';
    const char sideR = 'R';
    const char uplo = 'U';
    const double one = 1.0, zero = 0.0;

    double *X = REAL(X_in);
    double *Minv = REAL(Minv_in);

    /* Minv2 = Minv' * Minv */
    SEXP Minv2_sxp = PROTECT(Rf_allocMatrix(REALSXP, p, p));
    double *Minv2 = REAL(Minv2_sxp);

    F77_CALL(dsymm)(&sideR, &uplo, &p, &p, &one, Minv, &p, Minv, &p, &zero, Minv2, &p FCONE FCONE);

    /* XM = X * Minv */
    /* XM2 = X * Minv2 */
    SEXP XM_sxp = PROTECT(Rf_allocMatrix(REALSXP, n, p));
    SEXP XM2_sxp = PROTECT(Rf_allocMatrix(REALSXP, n, p));
    double *XM = REAL(XM_sxp);
    double *XM2 = REAL(XM2_sxp);

    F77_CALL(dsymm)(&sideR, &uplo, &n, &p, &one, Minv, &p, X, &n, &zero, XM, &n FCONE FCONE);
    F77_CALL(dsymm)(&sideR, &uplo, &n, &p, &one, Minv2, &p, X, &n, &zero, XM2, &n FCONE FCONE);

    /* uMu = sum_j X[i,j] * XM[i,j] */
    /* uM2u = sum_j X[i,j] * XM2[i,j] */
    SEXP uMu_sxp = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP uM2u_sxp = PROTECT(Rf_allocVector(REALSXP, n));
    double *uMu = REAL(uMu_sxp);
    double *uM2u = REAL(uM2u_sxp);

    for (int i = 0; i < n; i++)
    {
        double sum_i = 0.0, sum_i2 = 0.0;
        for (int j = 0; j < p; j++)
        {
            sum_i += X[i + n * j] * XM[i + n * j];
            sum_i2 += X[i + n * j] * XM2[i + n * j];
        }
        uMu[i] = sum_i;
        uM2u[i] = sum_i2;
    }

    /* vMv = uMu[wloc] */
    /* vM2v = uM2u[wloc] */
    SEXP vMv_sxp = PROTECT(Rf_allocVector(REALSXP, m));
    SEXP vM2v_sxp = PROTECT(Rf_allocVector(REALSXP, m));
    double *vMv = REAL(vMv_sxp);
    double *vM2v = REAL(vM2v_sxp);

    for (int i = 0; i < m; i++)
    {
        vMv[i] = uMu[iwloc[i] - 1];
        vM2v[i] = uM2u[iwloc[i] - 1];
        if(fabs(vMv[i] - 1.0) < 1e-8) {
            vMv[i] = 1.0 - 1e-8; // prevent division by zero later
        }
    }

    /* extract X[wloc,] */
    SEXP Xw_sxp = PROTECT(Rf_allocMatrix(REALSXP, m, p));
    double *Xw = REAL(Xw_sxp);

    for (int i = 0; i < m; i++)
        for (int j = 0; j < p; j++)
            Xw[i + m * j] = X[(iwloc[i] - 1) + n * j];

    /* XwM = X[wloc,] * Minv */
    /* XwM2 = X[wloc,] * Minv2 */
    SEXP XwM_sxp = PROTECT(Rf_allocMatrix(REALSXP, m, p));
    SEXP XwM2_sxp = PROTECT(Rf_allocMatrix(REALSXP, m, p));
    double *XwM = REAL(XwM_sxp);
    double *XwM2 = REAL(XwM2_sxp);

    F77_CALL(dsymm)(&sideR, &uplo, &m, &p, &one, Minv, &p, Xw, &m, &zero, XwM, &m FCONE FCONE);
    F77_CALL(dsymm)(&sideR, &uplo, &m, &p, &one, Minv2, &p, Xw, &m, &zero, XwM2, &m FCONE FCONE);

    /* vMu = X[wloc,] * Minv * t(X) -> m x n */
    /* vM2u = X[wloc,] * Minv2 * t(X) -> m x n */
    SEXP vMu_sxp = PROTECT(Rf_allocMatrix(REALSXP, m, n));
    SEXP vM2u_sxp = PROTECT(Rf_allocMatrix(REALSXP, m, n));
    double *vMu = REAL(vMu_sxp);
    double *vM2u = REAL(vM2u_sxp);

    F77_CALL(dgemm)(&transN, &transT, &m, &n, &p, &one, XwM, &m, X, &n, &zero, vMu, &m FCONE FCONE);
    F77_CALL(dgemm)(&transN, &transT, &m, &n, &p, &one, XwM2, &m, X, &n, &zero, vM2u, &m FCONE FCONE);

    /* delta function: add xi, remove xj
	    u2 = xj.M2.xj
	    v2 = xi.M2.xi
	    uv = xi.M2.xj
	    s = xi.M1.xi + (xj.M1.xi)^2 / (1 - xj.M1.xj)
	    w2 = v2 + 2 * xj.M1.xi / (1 - xj.M1.xj) * (uv) + (xj.M1.xi)^2 / (1 - xj.M1.xj)^2 * u2
	    delta(xi, xj) = u2 / (1 - xj.M1.xj) - w^2 / (1 + s) */
    SEXP delta_sxp = PROTECT(Rf_allocMatrix(REALSXP, m, n));
    double *delta = REAL(delta_sxp);
    int idx;
    double denom1, vMu2, s_ij, wnorm;

    for (int i = 0; i < m; i++)
    {
        denom1 = 1.0 - vMv[i];
        for (int j = 0; j < n; j++)
        {
            idx = i + m * j; /* column-major indexing */
            vMu2 = vMu[idx] * vMu[idx];
            s_ij = uMu[j] + vMu2 / denom1;
            wnorm = uM2u[j] + 2.0 * vMu[idx] / denom1 * vM2u[idx] + vMu2 / (denom1 * denom1) * vM2v[i];
            delta[idx] = -(vM2v[i] / denom1 - wnorm / (1.0 + s_ij));
        }
    }

    UNPROTECT(13);
    return delta_sxp;
}

/**
 * @brief I-optimality delta matrix (exchanging 1 design point)
 *
 * @param[in] wloc R integer vector, indices of design points
 * @param[in] X_in R matrix, design matrix
 * @param[in] XtX_in R matrix, X'X matrix
 * @param[in] M_in R matrix, information matrix
 * @param[in] Minv_in R matrix, inverse information matrix
 * @return delta matrix (length(wloc) x n)
 */
SEXP i_delta(SEXP wloc, SEXP X_in, SEXP XtX_in, SEXP M_in, SEXP Minv_in)
{
    int *iwloc = INTEGER(wloc);
    int m = Rf_length(wloc);

    int *dims = INTEGER(Rf_getAttrib(X_in, R_DimSymbol));
    int n = dims[0];
    int p = dims[1];
    const char transN = 'N';
    const char transT = 'T';
    const char sideL = 'L';
    const char sideR = 'R';
    const char uplo = 'U';
    const double one = 1.0, zero = 0.0;

    double *X = REAL(X_in);
    double *Minv = REAL(Minv_in);
    double *XtX = REAL(XtX_in);

    SEXP Xt_sxp = PROTECT(Rf_allocMatrix(REALSXP, p, n));
    double *Xt = REAL(Xt_sxp);
    for (R_len_t j = 0; j < p; j++)
        for (R_len_t i = 0; i < n; i++)
            Xt[j + p * i] = X[i + n * j];

    /* uMinv = Minv * X' */
    /* XuM = XtX * uMinv */
    SEXP uMinv_sxp = PROTECT(Rf_allocMatrix(REALSXP, p, n));
    SEXP XuM_sxp = PROTECT(Rf_allocMatrix(REALSXP, p, n));
    double *uMinv = REAL(uMinv_sxp);
    double *XuM = REAL(XuM_sxp);

    F77_CALL(dsymm)(&sideL, &uplo, &p, &n, &one, Minv, &p, Xt, &p, &zero, uMinv, &p FCONE FCONE);
    F77_CALL(dsymm)(&sideL, &uplo, &p, &n, &one, XtX, &p, uMinv, &p, &zero, XuM, &p FCONE FCONE);

    /* uMu = sum_j X[i,j] * uMinv[j,i] */
    /* uHu = sum_j uMinv[j,i] * XuM[j, i] */
    SEXP uMu_sxp = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP uHu_sxp = PROTECT(Rf_allocVector(REALSXP, n));
    double *uMu = REAL(uMu_sxp);
    double *uHu = REAL(uHu_sxp);

    for (int i = 0; i < n; i++)
    {
        double sum_i = 0.0, sum_i2 = 0.0;
        for (int j = 0; j < p; j++)
        {
            sum_i += X[i + n * j] * uMinv[j + p * i];
            sum_i2 += uMinv[j + p * i] * XuM[j + p * i];
        }
        uMu[i] = sum_i;
        uHu[i] = sum_i2;
    }

    /* vMv = uMu[wloc] */
    /* vHv = uHu[wloc] */
    SEXP vMv_sxp = PROTECT(Rf_allocVector(REALSXP, m));
    SEXP vHv_sxp = PROTECT(Rf_allocVector(REALSXP, m));
    double *vMv = REAL(vMv_sxp);
    double *vHv = REAL(vHv_sxp);

    for (int i = 0; i < m; i++)
    {
        vMv[i] = uMu[iwloc[i] - 1];
        vHv[i] = uHu[iwloc[i] - 1];
        if (fabs(vMv[i] - 1.0) < 1e-8)
        {
            vMv[i] = 1.0 - 1e-8; // prevent division by zero later
        }
    }

    /* extract X[wloc,], uMinv[wloc, ]' */
    SEXP Xw_sxp = PROTECT(Rf_allocMatrix(REALSXP, m, p));
    SEXP uMinvw_sxp = PROTECT(Rf_allocMatrix(REALSXP, m, p));
    double *Xw = REAL(Xw_sxp);
    double *uMinvw = REAL(uMinvw_sxp);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < p; j++)
        {
            Xw[i + m * j] = X[(iwloc[i] - 1) + n * j];
            uMinvw[i + m * j] = uMinv[j + p * (iwloc[i] - 1)];
        }
    }

    /* vM = X[wloc,] * Minv */
    /* vH = uMinv[wloc, ]' * XtX */
    SEXP vM_sxp = PROTECT(Rf_allocMatrix(REALSXP, m, p));
    SEXP vH_sxp = PROTECT(Rf_allocMatrix(REALSXP, m, p));
    double *vM = REAL(vM_sxp);
    double *vH = REAL(vH_sxp);

    F77_CALL(dsymm)(&sideR, &uplo, &m, &p, &one, Minv, &p, Xw, &m, &zero, vM, &m FCONE FCONE);
    F77_CALL(dsymm)(&sideR, &uplo, &m, &p, &one, XtX, &p, uMinvw, &m, &zero, vH, &m FCONE FCONE);

    /* vMu = X[wloc,] * Minv * X' -> m x n */
    /* vHu = uMinv[wloc, ]' * XtX * uMinv -> m x n */
    SEXP vMu_sxp = PROTECT(Rf_allocMatrix(REALSXP, m, n));
    SEXP vHu_sxp = PROTECT(Rf_allocMatrix(REALSXP, m, n));
    double *vMu = REAL(vMu_sxp);
    double *vHu = REAL(vHu_sxp);

    F77_CALL(dgemm)(&transN, &transT, &m, &n, &p, &one, vM, &m, X, &n, &zero, vMu, &m FCONE FCONE);
    F77_CALL(dgemm)(&transN, &transN, &m, &n, &p, &one, vH, &m, uMinv, &p, &zero, vHu, &m FCONE FCONE);

    /* delta function: add xi, remove xj
       uHu = xj.M2.H.M2.xj
       vHv = xi.M2.A.M2.xi
       uHv = xi.M2.A.M2.xj
       s = xi.M1.xi + (xj.M1.xi)^2 / (1 - xj.M1.xj)
       wHw = vHv + 2 * xj.M1.xi / (1 - xj.M1.xj) * uHv + (xj.M1.xi)^2 / (1 - xj.M1.xj)^2 * uHu
       delta(xi, xj) = u.H.u / (1 - xj.M1.xj) - w.H.w / (1 + s) */
    SEXP delta_sxp = PROTECT(Rf_allocMatrix(REALSXP, m, n));
    double *delta = REAL(delta_sxp);
    int idx;
    double denom1, vMu2, s_ij, wnorm;

    for (int i = 0; i < m; i++)
    {
        denom1 = 1.0 - vMv[i];
        for (int j = 0; j < n; j++)
        {
            idx = i + m * j; /* column-major indexing */
            vMu2 = vMu[idx] * vMu[idx];
            s_ij = uMu[j] + vMu2 / denom1;
            wnorm = uHu[j] + 2.0 * vMu[idx] / denom1 * vHu[idx] + vMu2 / (denom1 * denom1) * vHv[i];
            delta[idx] = -(vHv[i] / denom1 - wnorm / (1.0 + s_ij));
        }
    }

    UNPROTECT(14);
    return delta_sxp;
}

/**
 * @brief G-optimality delta matrix (exchanging 1 design point)
 *
 * @param[in] wloc R integer vector, indices of design points
 * @param[in] X_in R matrix, design matrix
 * @param[in] M_in R matrix, information matrix
 * @param[in] Minv_in R matrix, inverse information matrix
 * @return delta matrix (length(wloc) x n)
 */
SEXP g_delta(SEXP wloc, SEXP X_in, SEXP M_in, SEXP Minv_in)
{
    int *iwloc = INTEGER(wloc);
    int m = Rf_length(wloc);

    int *dims = INTEGER(Rf_getAttrib(X_in, R_DimSymbol));
    int n = dims[0];
    int p = dims[1];
    const char transN = 'N';
    const char transT = 'T';
    const char sideR = 'R';
    const char uplo = 'U';
    const double one = 1.0, zero = 0.0;

    double *X = REAL(X_in);
    double *Minv = REAL(Minv_in);

    /* uMinv = X * Minv */
    SEXP uMinv_sxp = PROTECT(Rf_allocMatrix(REALSXP, n, p));
    double *uMinv = REAL(uMinv_sxp);

    F77_CALL(dsymm)(&sideR, &uplo, &n, &p, &one, Minv, &p, X, &n, &zero, uMinv, &n FCONE FCONE);

    /* zMz = uMinv * X' */
    SEXP zMz_sxp = PROTECT(Rf_allocMatrix(REALSXP, n, n));
    double *zMz = REAL(zMz_sxp);

    F77_CALL(dgemm)(&transN, &transT, &n, &n, &p, &one, uMinv, &n, X, &n, &zero, zMz, &n FCONE FCONE);

    /* uMu = diag(zMz) */
    SEXP uMu_sxp = PROTECT(Rf_allocVector(REALSXP, n));
    double *uMu = REAL(uMu_sxp);
    for (int i = 0; i < n; i++)
        uMu[i] = zMz[i + n * i];

    /* vMv = uMu[wloc] */
    SEXP vMv_sxp = PROTECT(Rf_allocVector(REALSXP, m));
    double *vMv = REAL(vMv_sxp);
    for (int i = 0; i < m; i++)
    {
        vMv[i] = uMu[iwloc[i] - 1];
        if (fabs(vMv[i] - 1.0) < 1e-8)
            vMv[i] = 1.0 - 1e-8; // prevent division by zero later
    }

    /* vMu = zMz[wloc, ] */
    SEXP vMu_sxp = PROTECT(Rf_allocMatrix(REALSXP, m, n));
    double *vMu = REAL(vMu_sxp);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            vMu[i + m * j] = zMz[(iwloc[i] - 1) + n * j];

    /* delta function: add xi, remove xj
	   s = xi.M1.xi + (xj.M1.xi)^2 / (1 - xj.M1.xj)
	   z.w =  xi.M1.z + xj.M1.xi / (1 - xj.M1.xj) * xj.M1.z
	   delta(xi, xj) = z.M1.z + xj.M1.z / (1 - xj.M1.xj) - (z.w)^2 / (1 + s) */
    SEXP delta_sxp = PROTECT(Rf_allocMatrix(REALSXP, m, n));
    SEXP s_sxp = PROTECT(Rf_allocMatrix(REALSXP, m, n));
    SEXP zw2_sxp = PROTECT(Rf_allocMatrix(REALSXP, m, n));
    
    double *delta = REAL(delta_sxp);
    double *s = REAL(s_sxp);
    double *zw2 = REAL(zw2_sxp);

    int idx;
    double denom1, zw, wnorm, sw_max;
    double *s_wnorm = (double *)R_alloc(n, sizeof(double)); // temporary storage
    double max_uMu = 1.0 / max_dbl(uMu, n);

    for (int i = 0; i < m; i++)
    {
        denom1 = 1.0 - vMv[i];
        for (int j = 0; j < n; j++)
        {
            idx = i + m * j; 
            zw2[idx] = vMu[idx] / denom1;
            s[idx] = uMu[j] + vMu[idx] * zw2[idx];
        }
        for (int j = 0; j < n; j++)
        {
            idx = i + m * j;
            for (int k = 0; k < n; k++)
            {
                zw = zMz[k + n * j] + zw2[i + m * k] * vMu[idx];
                wnorm = (zw * zw) / (1.0 + s[idx]);
                s_wnorm[k] = s[i + m * k] - wnorm;
            }
            sw_max = max_dbl(s_wnorm, n);
            delta[idx] = 1.0 / sw_max - max_uMu;
        }
    }

    UNPROTECT(8);
    return delta_sxp;
}

/**
 * @brief Alias-optimality delta matrix (exchanging 1 design point)
 *
 * @param[in] wloc R integer vector, indices of design points
 * @param[in] X_in R matrix, design matrix
 * @param[in] M_in R matrix, information matrix
 * @param[in] ortho_in R matrix (optional), orthogonality constraint matrix
 * @return delta matrix (length(wloc) x n)
 */
SEXP alias_delta(SEXP wloc, SEXP X_in, SEXP M_in, SEXP ortho_in)
{
    int *iwloc = INTEGER(wloc);
    int m = Rf_length(wloc);

    int *dims = INTEGER(Rf_getAttrib(X_in, R_DimSymbol));
    int n = dims[0];
    int p = dims[1];
    const char uplo = 'L';
    const double one = 1.0;
    const int ione = 1; 

    double *X = REAL(X_in);
    double *M = REAL(M_in);
    double *ortho = NULL;
    
    /* initial correlation */
    double corr0, corr1;
    if(!Rf_isNull(ortho_in))
    {
        ortho = REAL(ortho_in);
        corr0 = mean_abs_corr_ortho(M, ortho, p);
    }
    else 
        corr0 = mean_abs_corr(M, p);

    SEXP delta_sxp = PROTECT(Rf_allocMatrix(REALSXP, m, n));
    double *delta = REAL(delta_sxp);

    /* workspaces */
    double *M_tmp = (double *)R_alloc(p * p, sizeof(double));
    double *M1 = (double *)R_alloc(p * p, sizeof(double));
    double *Xw = (double *)R_alloc(p * p * m, sizeof(double));
    double *Xi = (double *)R_alloc(p, sizeof(double));
    double *Xj = (double *)R_alloc(p, sizeof(double));

    /* pre-compute X[wi, ] * t(X[wi, ]) */
    for(R_len_t i = 0; i < m; i++) 
    {
        R_len_t wi = iwloc[i] - 1;
        for(R_len_t k = 0; k < p; k++)
            Xi[k] = X[wi + n * k];

        memset(M_tmp,  0, p * p * sizeof(double));
        F77_CALL(dsyr)(&uplo, &p, &one, Xi, &ione, M_tmp, &p FCONE);
        
        for (R_len_t k = 0; k < p * p; k++) 
            Xw[i * p * p + k] = M_tmp[k];
    }

    /* populate delta */
    for (R_len_t j = 0; j < n; j++) 
    {
        for (R_len_t k = 0; k < p; k++)
            Xj[k] = X[j + n * k];

        /* M1 = M + X[j, ] * t(X[j, ]) */
        memcpy(M_tmp, M, p * p * sizeof(double));
        F77_CALL(dsyr)(&uplo, &p, &one, Xj, &ione, M_tmp, &p FCONE);

        for(R_len_t i = 0; i < m; i++) 
        {
            /* M1 -= X[wi, ] * t(X[wi, ]) */
            memcpy(M1, M_tmp, p * p * sizeof(double));
            for(R_len_t k = 0; k < p * p; k++) 
                M1[k] -= Xw[i * p * p + k];

            if(ortho)
                corr1 = mean_abs_corr_ortho(M1, ortho, p);
            else 
                corr1 = mean_abs_corr(M1, p);

            delta[i + j * m] = corr0 - corr1;
        }
    }

    UNPROTECT(1);
    return delta_sxp;
}
