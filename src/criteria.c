#define R_NO_REMAP

#include "cvxoptdes.h"

static inline int mat_idx(int i, int j, int n)
{
    return i + n * j;
}

static inline int lower_tri_idx(int i, int j, int k)
{
    return i > j ? (k - 1) * j + (i - j - 1) - (j * (j - 1)) / 2 : -1;
}

static inline int lower_tri_idx_diag(int i, int j, int k)
{
    return i >= j ? k * j + i - (j * (j + 1)) / 2 : -1;
}

static inline int lower_tri_diag(int j, int k)
{
    return j * (2 * k - j + 1) / 2;
}

static inline int array_idx(int i, int j, int k, int p)
{
    return i + p * j + k * p * p;
}

void D_opt_populate_Abc(pdata *pars, R_len_t p, R_len_t n, R_len_t nnz, scs_int mtot, scs_int ntot, double alpha, double lambda,
                        double gamma, double upr, scs_int nlvls, R_len_t z_ortho, R_len_t q_ortho, R_len_t mzeros, R_len_t mmeans)
{
    R_len_t np = n * p;
    R_len_t pp = p * p;
    R_len_t pp2 = p * (p + 1) / 2;
    R_len_t brow = 1 + pp + z_ortho;
    R_len_t ntot0 = ntot - 2 * pp;
    /* mixed-model SCA */
    double *Vptr = NULL;
    R_len_t blocks = 0;
    if (!Rf_isNull(pars->v_sca))
    {
        Vptr = REAL(pars->v_sca);
        blocks = Rf_ncols(pars->v_sca);
    }
    R_len_t ntotb = ntot - 2 * blocks * p;

    /* populate b and c vectors */
    ((pars->d)->b)[0] = 1.0;
    for (R_len_t i = 0; i < p; i++)
    {
        ((pars->d)->c)[n + np + pp2 + np + i] = -1.0 / p;
        // ((pars->d)->b)[(mtot - 3 * p) + 3 * i + 1] = 1.0;   // still needed?
    }
    // optional: cost penalty
    if (!Rf_isNull(pars->f) && alpha > 0)
    {
        for (R_len_t s = 0; s < n; s++)
        {
            ((pars->d)->c)[s] = alpha * REAL_ELT(pars->f, s);
        }
    }
    // optional: strata constraints
    if (nlvls > 0)
    {
        double *prop = REAL(VECTOR_ELT(pars->strata, 1));
        for (R_len_t s = 0; s < nlvls; s++)
            ((pars->d)->b)[brow + s] = prop[s];
        brow += nlvls;
    }
    // optional: zero entry constraints
    if (mzeros > 0)
    {
        for (R_len_t i = 0; i < p; i++)
        {
            double wzero = REAL_ELT(pars->zeros, i);
            if (!R_IsNA(wzero))
            {
                ((pars->d)->b)[brow] = wzero;
                brow++;
            }
        }
    }
    // optional: mean constraints
    if (mmeans > 0)
    {
        for (R_len_t i = 0; i < p; i++)
        {
            double meanval = REAL_ELT(pars->colmeans, i);
            if (!R_IsNA(meanval))
            {
                ((pars->d)->b)[brow] = meanval;
                brow++;
            }
        }
    }
    // optional: augmented design
    if (!Rf_isNull(pars->w0) && gamma > 0)
    {
        for (R_len_t s = 0; s < n; s++)
        {
            ((pars->d)->b)[brow + s] = -(1.0 - gamma) * REAL_ELT(pars->w0, s);
        }
    }
    // optional: weight upper bound
    if (upr < 1)
    {
        for (R_len_t s = 0; s < n; s++)
            ((pars->d)->b)[brow + n + s] = upr;
    }

    /* populate A matrix */
    double *Jptr = REAL(pars->J);
    double *ortho = NULL;
    R_len_t ni = n;
    R_len_t irow = 0;
    R_len_t icol = 0;
    R_len_t *Ac = (R_len_t *)S_alloc(nnz, sizeof(R_len_t)); // COO column indices
    R_len_t *Ai = (R_len_t *)S_alloc(nnz, sizeof(R_len_t)); // COO row indices
    double *Ax = (double *)S_alloc(nnz, sizeof(double));    // COO nz values

    // COO format
    // sum(w) = 1
    for (R_len_t s = 0; s < n; s++)
    {
        Ax[s] = 1.0;
        Ai[s] = 0;
        Ac[s] = s;
    }
    irow += 1;
    // Ai * Aj = 0 (exact orthogonality)
    if (z_ortho > 0)
    {
        ortho = REAL(pars->ortho);
        R_len_t li = 0;
        for (R_len_t i = 1; i < p; i++)
        {
            for (R_len_t j = 0; j <= i; j++)
            {
                if (ortho[mat_idx(j, i, p)] != R_PosInf && ortho[mat_idx(j, i, p)] < 1e-4)
                {
                    for (R_len_t s = 0; s < n; s++)
                    {
                        Ax[ni] = Jptr[mat_idx(s, i, n)] * Jptr[mat_idx(s, j, n)]; // A[, i] * A[, j] = 0
                        Ai[ni] = irow + li;
                        Ac[ni] = s;
                        ni++;
                    }
                    li++;
                }
            }
        }
        irow += z_ortho;
    }
    // A * Z = J
    icol = n;
    for (R_len_t j = 0; j < p; j++)
    {
        for (R_len_t i = 0; i < p; i++)
        {
            for (R_len_t s = 0; s < n; s++)
            {
                Ax[ni] = Jptr[mat_idx(s, i, n)]; // A[i] * Z[j]
                Ai[ni] = irow + j * p + i;
                Ac[ni] = icol + j * n + s;
                ni++;
            }
            if (lambda > 0)
            {
                Ax[ni] = 1.0; // e_i * Z[j]
                Ai[ni] = irow + j * p + i;
                Ac[ni] = ntot0 + j * p + i;
                ni++;
            }
            if (blocks > 0)
            {
                for (R_len_t b = 0; b < blocks; b++)
                {
                    Ax[ni] = Vptr[mat_idx(i, b, p)];  // v[i] * Z[j]
                    Ai[ni] = irow + j * p + i;
                    Ac[ni] = ntotb + j * blocks + b;
                    ni++;
                }
            }
            if (i >= j)
            {
                Ax[ni] = -1.0; // J[i, j]
                Ai[ni] = irow + j * p + i;
                Ac[ni] = icol + np + lower_tri_idx_diag(i, j, p);
                ni++;
            }
        }
    }
    // strata
    irow += pp;
    if (nlvls > 0)
    {
        int *strata_id = INTEGER(VECTOR_ELT(pars->strata, 0));
        for (R_len_t s = 0; s < n; s++)
        {
            Ax[ni] = 1.0;
            Ai[ni] = irow + strata_id[s];
            Ac[ni] = s;
            ni++;
        }
        irow += nlvls;
    }
    // zero entry weights
    if (mzeros > 0)
    {
        int *Xzero_ptr = INTEGER(pars->Xzero);
        for (R_len_t i = 0; i < p; i++)
        {
            double wzero = REAL_ELT(pars->zeros, i);
            if (!R_IsNA(wzero))
            {
                for (R_len_t s = 0; s < n; s++)
                {
                    Ax[ni] = (double)(Xzero_ptr[mat_idx(s, i, n)]);
                    Ai[ni] = irow;
                    Ac[ni] = s;
                    ni++;
                }
                irow++;
            }
        }
    }
    // column means
    if (mmeans > 0)
    {
        for (R_len_t i = 0; i < p; i++)
        {
            double meanval = REAL_ELT(pars->colmeans, i);
            if (!R_IsNA(meanval))
            {
                for (R_len_t s = 0; s < n; s++)
                {
                    Ax[ni] = Jptr[mat_idx(s, i, n)];
                    Ai[ni] = irow;
                    Ac[ni] = s;
                    ni++;
                }
                irow++;
            }
        }
    }
    // w >= 0
    for (R_len_t s = 0; s < n; s++)
    {
        Ax[ni] = -1.0;
        Ai[ni] = irow + s;
        Ac[ni] = s;
        ni++;
    }
    // w <= upr
    if (upr < 1)
    {
        irow += n;
        for (R_len_t s = 0; s < n; s++)
        {
            Ax[ni] = 1.0;
            Ai[ni] = irow + s;
            Ac[ni] = s;
            ni++;
        }
    }
    // t >= 0
    irow += n;
    icol = n + np + pp2;
    for (R_len_t s = 0; s < np; s++)
    {
        Ax[ni] = -1.0;
        Ai[ni] = irow + s;
        Ac[ni] = icol + s;
        ni++;
    }
    irow += np;
    if (lambda > 0)
    {
        for (R_len_t s = 0; s < pp; s++)
        {
            Ax[ni] = -1.0;
            Ai[ni] = irow + s;
            Ac[ni] = ntot0 + pp + s;
            ni++;
        }
        irow += pp;
    }
    if(blocks > 0) 
    {
        for (R_len_t s = 0; s < blocks * p; s++)
        {
            Ax[ni] = -1.0;
            Ai[ni] = irow + s;
            Ac[ni] = ntotb + blocks * p + s;
            ni++;
        }
        irow += blocks * p;
    }
    // sum(t[j]) <= J[j]
    icol = n + np;
    for (R_len_t j = 0; j < p; j++)
    {
        for (R_len_t s = 0; s < n; s++)
        {
            Ax[ni] = 1.0; // sum(t[j])
            Ai[ni] = irow + j;
            Ac[ni] = icol + pp2 + j * n + s;
            ni++;
        }
        if (lambda > 0)
        {
            for (R_len_t k = 0; k < p; k++)
            {
                Ax[ni] = 1.0; // sum(t_lam[j])
                Ai[ni] = irow + j;
                Ac[ni] = ntot0 + pp + j * p + k;
                ni++;
            }
        }
        if(blocks > 0)
        {
            for (R_len_t b = 0; b < blocks; b++)
            {
                Ax[ni] = 1.0; // sum(t_b[j])
                Ai[ni] = irow + j;
                Ac[ni] = ntotb + blocks * p + j * blocks + b;
                ni++;
            }
        }
        Ax[ni] = -1.0; // J[j]
        Ai[ni] = irow + j;
        Ac[ni] = icol + lower_tri_diag(j, p);
        ni++;
    }
    // || (2 * Z[ij], t[ij] - w[i]) || <= t[ij] + w[i]
    irow += p;
    for (R_len_t j = 0; j < p; j++)
    {
        for (R_len_t s = 0; s < n; s++)
        {
            // t[ij] + w[i]
            Ax[ni] = -1.0;
            Ax[ni + 1] = -1.0;
            Ai[ni] = irow;
            Ai[ni + 1] = irow;
            Ac[ni] = s;
            Ac[ni + 1] = n + np + pp2 + j * n + s;
            // 2 * Z[ij]
            Ax[ni + 2] = -2.0;
            Ai[ni + 2] = irow + 1;
            Ac[ni + 2] = n + j * n + s;
            // t[ij] - w[i]
            Ax[ni + 3] = 1.0;
            Ax[ni + 4] = -1.0;
            Ai[ni + 3] = irow + 2;
            Ai[ni + 4] = irow + 2;
            Ac[ni + 3] = s;
            Ac[ni + 4] = n + np + pp2 + j * n + s;
            ni += 5;
            irow += 3;
        }
    }
    if (lambda > 0)
    {
        // || (2 * Z_[ij], t[ij] - lambda) || <= t[ij] + lambda
        for (R_len_t j = 0; j < p; j++)
        {
            for (R_len_t i = 0; i < p; i++)
            {
                // t[ij] + lambda
                Ax[ni] = -1.0;
                Ai[ni] = irow;
                Ac[ni] = ntot0 + pp + j * p + i;
                ((pars->d)->b)[irow] = lambda;
                // 2 * Z[ij]
                Ax[ni + 1] = -2.0;
                Ai[ni + 1] = irow + 1;
                Ac[ni + 1] = ntot0 + j * p + i;
                // t[ij] - lambda
                Ax[ni + 2] = -1.0;
                Ai[ni + 2] = irow + 2;
                Ac[ni + 2] = ntot0 + pp + j * p + i;
                ((pars->d)->b)[irow + 2] = -lambda;
                ni += 3;
                irow += 3;
            }
        }
    }
    if(blocks > 0)
    {
        // || (2 * Z_[ij], t[ij] - 1) || <= t[ij] + 1
        for (R_len_t j = 0; j < p; j++)
        {
            for (R_len_t b = 0; b < blocks; b++)
            {
                // t[ij] + 1
                Ax[ni] = -1.0;
                Ai[ni] = irow;
                Ac[ni] = ntotb + blocks * p + j * blocks + b;
                ((pars->d)->b)[irow] = 1.0;
                // 2 * Z[ij]
                Ax[ni + 1] = -2.0;
                Ai[ni + 1] = irow + 1;
                Ac[ni + 1] = ntotb + j * blocks + b;
                // t[ij] - 1
                Ax[ni + 2] = -1.0;
                Ai[ni + 2] = irow + 2;
                Ac[ni + 2] = ntotb + blocks * p + j * blocks + b;
                ((pars->d)->b)[irow + 2] = -1.0;
                ni += 3;
                irow += 3;
            }
        }
    }
    // |M[ij]| / sqrt(M[ii] * M[jj]) <= rho
    // RSOC: || (2 * M[ij], rho * M[ii] - rho * M[jj]) || <= rho * M[ii] + rho * M[jj]
    if (q_ortho > 0)
    {
        ortho = REAL(pars->ortho);
        for (R_len_t i = 1; i < p; i++)
        {
            for (R_len_t j = 0; j <= i; j++)
            {
                if (ortho[mat_idx(j, i, p)] != R_PosInf && ortho[mat_idx(j, i, p)] >= 1e-4)
                {
                    // rho * (M[ii] + M[jj])
                    double rho = ortho[mat_idx(j, i, p)];
                    for (R_len_t s = 0; s < n; s++)
                    {
                        Ax[ni] = -rho * (Jptr[mat_idx(s, i, n)] * Jptr[mat_idx(s, i, n)] + Jptr[mat_idx(s, j, n)] * Jptr[mat_idx(s, j, n)]);
                        Ai[ni] = irow;
                        Ac[ni] = s;
                        ni++;
                    }
                    // 2 * M[ij]
                    for (R_len_t s = 0; s < n; s++)
                    {
                        Ax[ni] = -2.0 * Jptr[mat_idx(s, i, n)] * Jptr[mat_idx(s, j, n)];
                        Ai[ni] = irow + 1;
                        Ac[ni] = s;
                        ni++;
                    }
                    // rho * (M[ii] - M[jj])
                    for (R_len_t s = 0; s < n; s++)
                    {
                        Ax[ni] = -rho * (Jptr[mat_idx(s, i, n)] * Jptr[mat_idx(s, i, n)] - Jptr[mat_idx(s, j, n)] * Jptr[mat_idx(s, j, n)]);
                        Ai[ni] = irow + 2;
                        Ac[ni] = s;
                        ni++;
                    }
                    irow += 3;
                }
            }
        }
    }

    // exp(v[j]) <= J[j]
    icol = n + np + pp2 + np;
    for (R_len_t j = 0; j < p; j++)
    {
        Ax[ni] = -1.0;     // v[j]
        Ax[ni + 1] = -1.0; // J[j]
        Ai[ni] = irow + 3 * j;
        Ai[ni + 1] = irow + 3 * j + 2;
        Ac[ni] = icol + j;
        Ac[ni + 1] = n + np + lower_tri_diag(j, p);
        ni += 2;
    }

    // CSC format
    // column pointers
    R_len_t *ccount = (R_len_t *)S_alloc(ntot, sizeof(R_len_t));
    for (R_len_t i = 0; i < nnz; i++)
        (((pars->d)->A)->p)[Ac[i] + 1] += 1;
    for (R_len_t s = 0; s < ntot; s++)
        (((pars->d)->A)->p)[s + 1] += (((pars->d)->A)->p)[s];
    // re-assign value/row indices
    // initialize column counts
    for (R_len_t s = 0; s < ntot; s++)
        ccount[s] = (((pars->d)->A)->p)[s];
    for (R_len_t i = 0; i < nnz; i++)
    {
        icol = Ac[i];
        ni = ccount[icol];
        (((pars->d)->A)->i)[ni] = Ai[i];
        (((pars->d)->A)->x)[ni] = Ax[i];
        (ccount[icol])++;
    }
}

void D_opt_populate_K(pdata *pars, R_len_t p, R_len_t n, double upr, scs_int nlvls, R_len_t z_ortho, R_len_t qsize, R_len_t mzeros, R_len_t mmeans, double lambda, R_len_t blocks)
{
    R_len_t np = n * p;
    R_len_t pp = p * p;
    (pars->k)->z = 1 + pp + nlvls + z_ortho + mzeros + mmeans;
    (pars->k)->l = (1 + (upr < 1)) * n + np + p + (lambda > 0) * pp + blocks * p;
    for (R_len_t s = 0; s < qsize; s++)
        ((pars->k)->q)[s] = 3;
    (pars->k)->qsize = qsize;
    (pars->k)->ep = p;
}

void D_spec_opt_populate_Abc(pdata *pars, R_len_t p, R_len_t n, R_len_t nnz, scs_int mtot, scs_int ntot, double alpha, double lambda,
                             double gamma, double upr, scs_int nlvls, R_len_t z_ortho, R_len_t q_ortho, R_len_t mzeros, R_len_t mmeans)
{
    // initialize
    R_len_t brow = 2 + z_ortho;
    double sqrt2 = sqrt(2.0);
    R_len_t pp2 = p * (p + 1) / 2;

    /* populate b and c vectors */
    ((pars->d)->b)[0] = 1.0;
    ((pars->d)->b)[1] = 1.0;     // v = 1
    ((pars->d)->c)[n] = 1.0 / p; // minimize t / p

    // optional: cost penalty
    if (!Rf_isNull(pars->f) && alpha > 0)
    {
        for (R_len_t s = 0; s < n; s++)
        {
            ((pars->d)->c)[s] = alpha * REAL_ELT(pars->f, s);
        }
    }
    // optional: strata constraints
    if (nlvls > 0)
    {
        double *prop = REAL(VECTOR_ELT(pars->strata, 1));
        for (R_len_t s = 0; s < nlvls; s++)
            ((pars->d)->b)[brow + s] = prop[s];
        brow += nlvls;
    }
    // optional: zero entry constraints
    if (mzeros > 0)
    {
        for (R_len_t i = 0; i < p; i++)
        {
            double wzero = REAL_ELT(pars->zeros, i);
            if (!R_IsNA(wzero))
            {
                ((pars->d)->b)[brow] = wzero;
                brow++;
            }
        }
    }
    // optional: mean constraints
    if (mmeans > 0)
    {
        for (R_len_t i = 0; i < p; i++)
        {
            double meanval = REAL_ELT(pars->colmeans, i);
            if (!R_IsNA(meanval))
            {
                ((pars->d)->b)[brow] = meanval;
                brow++;
            }
        }
    }
    // optional: augmented design
    if (!Rf_isNull(pars->w0) && gamma > 0)
    {
        for (R_len_t s = 0; s < n; s++)
        {
            ((pars->d)->b)[brow] = -(1.0 - gamma) * REAL_ELT(pars->w0, s);
            brow++;
        }
    }
    // optional: weight upper bound
    if (upr < 1)
    {
        brow += n;
        for (R_len_t s = 0; s < n; s++)
        {
            ((pars->d)->b)[brow] = upr;
            brow++;
        }
    }
    // optional: sca penalty (M + vv')
    if (!Rf_isNull(pars->v_sca))
    {
        double *Vptr = REAL(pars->v_sca);
        brow = 0;
        for (R_len_t i = 0; i < p; i++)
        {
            for (R_len_t j = i; j < p; j++)
            {
                for (R_len_t b = 0; b < Rf_ncols(pars->v_sca); b++)
                    ((pars->d)->b)[mtot - pp2 + brow] += Vptr[mat_idx(i, b, p)] * Vptr[mat_idx(j, b, p)];
                if (j > i)
                    ((pars->d)->b)[mtot - pp2 + brow] *= sqrt2;
                brow++;
            }
        }
    }
    // optional: ridge penalty
    if (lambda > 0)
    {
        brow = 0;
        for (R_len_t i = 0; i < p; i++)
        {
            // add lambda to diagonal terms
            ((pars->d)->b)[mtot - pp2 + brow] -= lambda;
            brow += p - i;
        }
    }

    /* populate A matrix */
    double *Jptr = REAL(pars->J);
    double *ortho = NULL;
    R_len_t ni = 0;
    R_len_t irow = 0;
    R_len_t icol = 0;
    R_len_t *Ac = (R_len_t *)S_alloc(nnz, sizeof(R_len_t)); // COO column indices
    R_len_t *Ai = (R_len_t *)S_alloc(nnz, sizeof(R_len_t)); // COO row indices
    double *Ax = (double *)S_alloc(nnz, sizeof(double));    // COO nz values

    // COO format
    // sum(w) = 1
    for (R_len_t s = 0; s < n; s++)
    {
        Ax[s] = 1.0;
        Ai[s] = 0;
        Ac[s] = s;
    }
    ni += n;
    // v = 1
    Ax[ni] = 1.0;
    Ai[ni] = 1;
    Ac[ni] = n + 1;
    ni++;
    // Ai * Aj = 0 (exact orthogonality)
    irow += 2;
    if (z_ortho > 0)
    {
        ortho = REAL(pars->ortho);
        R_len_t li = 0;
        for (R_len_t i = 1; i < p; i++)
        {
            for (R_len_t j = 0; j <= i; j++)
            {
                if (ortho[mat_idx(j, i, p)] != R_PosInf && ortho[mat_idx(j, i, p)] < 1e-4)
                {
                    for (R_len_t s = 0; s < n; s++)
                    {
                        Ax[ni] = Jptr[mat_idx(s, i, n)] * Jptr[mat_idx(s, j, n)]; // A[, i] * A[, j] = 0
                        Ai[ni] = irow + li;
                        Ac[ni] = s;
                        ni++;
                    }
                    li++;
                }
            }
        }
        irow += z_ortho;
    }
    // strata
    if (nlvls > 0)
    {
        int *strata_id = INTEGER(VECTOR_ELT(pars->strata, 0));
        for (R_len_t s = 0; s < n; s++)
        {
            Ax[ni] = 1.0;
            Ai[ni] = irow + strata_id[s];
            Ac[ni] = s;
            ni++;
        }
        irow += nlvls;
    }
    // zero entry weights
    if (mzeros > 0)
    {
        int *Xzero_ptr = INTEGER(pars->Xzero);
        for (R_len_t i = 0; i < p; i++)
        {
            double wzero = REAL_ELT(pars->zeros, i);
            if (!R_IsNA(wzero))
            {
                for (R_len_t s = 0; s < n; s++)
                {
                    Ax[ni] = (double)(Xzero_ptr[mat_idx(s, i, n)]);
                    Ai[ni] = irow;
                    Ac[ni] = s;
                    ni++;
                }
                irow++;
            }
        }
    }
    // column means
    if (mmeans > 0)
    {
        for (R_len_t i = 0; i < p; i++)
        {
            double meanval = REAL_ELT(pars->colmeans, i);
            if (!R_IsNA(meanval))
            {
                for (R_len_t s = 0; s < n; s++)
                {
                    Ax[ni] = Jptr[mat_idx(s, i, n)];
                    Ai[ni] = irow;
                    Ac[ni] = s;
                    ni++;
                }
                irow++;
            }
        }
    }
    // w >= 0
    for (R_len_t s = 0; s < n; s++)
    {
        Ax[ni] = -1.0;
        Ai[ni] = irow + s;
        Ac[ni] = s;
        ni++;
    }
    irow += n;
    // w <= upr
    if (upr < 1)
    {
        for (R_len_t s = 0; s < n; s++)
        {
            Ax[ni] = 1.0;
            Ai[ni] = irow + s;
            Ac[ni] = s;
            ni++;
        }
        irow += n;
    }
    // |M[ij]| / sqrt(M[ii] * M[jj]) <= rho
    // RSOC: || (2 * M[ij], rho * M[ii] - rho * M[jj]) || <= rho * M[ii] + rho * M[jj]
    if (q_ortho > 0)
    {
        ortho = REAL(pars->ortho);
        for (R_len_t i = 1; i < p; i++)
        {
            for (R_len_t j = 0; j <= i; j++)
            {
                if (ortho[mat_idx(j, i, p)] != R_PosInf && ortho[mat_idx(j, i, p)] >= 1e-4)
                {
                    // rho * (M[ii] + M[jj])
                    double rho = ortho[mat_idx(j, i, p)];
                    for (R_len_t s = 0; s < n; s++)
                    {
                        Ax[ni] = -rho * (Jptr[mat_idx(s, i, n)] * Jptr[mat_idx(s, i, n)] + Jptr[mat_idx(s, j, n)] * Jptr[mat_idx(s, j, n)]);
                        Ai[ni] = irow;
                        Ac[ni] = s;
                        ni++;
                    }
                    // 2 * M[ij]
                    for (R_len_t s = 0; s < n; s++)
                    {
                        Ax[ni] = -2.0 * Jptr[mat_idx(s, i, n)] * Jptr[mat_idx(s, j, n)];
                        Ai[ni] = irow + 1;
                        Ac[ni] = s;
                        ni++;
                    }
                    // rho * (M[ii] - M[jj])
                    for (R_len_t s = 0; s < n; s++)
                    {
                        Ax[ni] = -rho * (Jptr[mat_idx(s, i, n)] * Jptr[mat_idx(s, i, n)] - Jptr[mat_idx(s, j, n)] * Jptr[mat_idx(s, j, n)]);
                        Ai[ni] = irow + 2;
                        Ac[ni] = s;
                        ni++;
                    }
                    irow += 3;
                }
            }
        }
    }
    // -log det(M) <= t
    Ax[ni] = -1.0;     // t
    Ax[ni + 1] = -1.0; // v
    Ai[ni] = irow;
    Ai[ni + 1] = irow + 1;
    Ac[ni] = n;
    Ac[ni + 1] = n + 1;
    ni += 2;
    irow += 2;
    // X = vec(M) = (M[11], sqrt(2) * M[21], ..., M[pp])
    for (R_len_t i = 0; i < p; i++)
    {
        for (R_len_t j = i; j < p; j++)
        {
            for (R_len_t s = 0; s < n; s++)
            {
                Ax[ni] = -Jptr[mat_idx(s, j, n)] * Jptr[mat_idx(s, i, n)];
                if (j > i)
                    Ax[ni] *= sqrt2;
                Ai[ni] = irow;
                Ac[ni] = s;
                ni++;
            }
            irow++;
        }
    }

    // CSC format
    // column pointers
    R_len_t *ccount = (R_len_t *)S_alloc(ntot, sizeof(R_len_t));
    for (R_len_t i = 0; i < nnz; i++)
        (((pars->d)->A)->p)[Ac[i] + 1] += 1;
    for (R_len_t s = 0; s < ntot; s++)
        (((pars->d)->A)->p)[s + 1] += (((pars->d)->A)->p)[s];
    // re-assign value/row indices
    // initialize column counts
    for (R_len_t s = 0; s < ntot; s++)
        ccount[s] = (((pars->d)->A)->p)[s];
    for (R_len_t i = 0; i < nnz; i++)
    {
        icol = Ac[i];
        ni = ccount[icol];
        (((pars->d)->A)->i)[ni] = Ai[i];
        (((pars->d)->A)->x)[ni] = Ax[i];
        (ccount[icol])++;
    }
}

void D_spec_opt_populate_K(pdata *pars, R_len_t p, R_len_t n, double upr, scs_int nlvls, R_len_t z_ortho, R_len_t qsize, R_len_t mzeros, R_len_t mmeans)
{
    (pars->k)->z = 2 + nlvls + z_ortho + mzeros + mmeans;
    (pars->k)->l = (1 + (upr < 1)) * n;
    if (qsize)
    {
        for (R_len_t s = 0; s < qsize; s++)
            ((pars->k)->q)[s] = 3;
        (pars->k)->qsize = qsize;
    }
    (pars->k)->dsize = 1;
    (pars->k)->d[0] = p;
}

void A_opt_populate_Abc(pdata *pars, R_len_t p, R_len_t n, R_len_t nnz, scs_int mtot, scs_int ntot, double alpha, double lambda,
                        double gamma, double upr, scs_int nlvls, R_len_t z_ortho, R_len_t q_ortho, R_len_t mzeros, R_len_t mmeans)
{
    R_len_t np = n * p;
    R_len_t pp = p * p;
    R_len_t brow = 1 + pp + z_ortho;
    /* mixed-model SCA */
    double *Vptr = NULL;
    R_len_t blocks = 0;
    R_len_t ntotb = 0;
    if (!Rf_isNull(pars->v_sca))
    {
        Vptr = REAL(pars->v_sca);
        blocks = Rf_ncols(pars->v_sca);
        ntotb = ntot - blocks * (1 + p);
    }

    /* populate b and c vectors */
    ((pars->d)->b)[0] = 1.0;
    for (R_len_t s = 0; s < n; s++)
    {
        ((pars->d)->c)[n + np + s] = -1.0;
    }
    if (lambda > 0)
    {
        for (R_len_t k = 0; k < p; k++)
            ((pars->d)->c)[2 * n + np + pp + k] = -1.0;
    }
    if (blocks > 0)
    {
        for (R_len_t b = 0; b < blocks; b++)
            ((pars->d)->c)[ntotb + p * blocks + b] = -1.0;
    }
    // optional: cost penalty
    if (!Rf_isNull(pars->f) && alpha > 0)
    {
        for (R_len_t s = 0; s < n; s++)
        {
            ((pars->d)->c)[s] = alpha * REAL_ELT(pars->f, s);
        }
    }
    // optional: strata constraints
    if (nlvls > 0)
    {
        double *prop = REAL(VECTOR_ELT(pars->strata, 1));
        for (R_len_t s = 0; s < nlvls; s++)
            ((pars->d)->b)[brow + s] = prop[s];
        brow += nlvls;
    }
    // optional: zero entry constraints
    if (mzeros > 0)
    {
        for (R_len_t i = 0; i < p; i++)
        {
            double wzero = REAL_ELT(pars->zeros, i);
            if (!R_IsNA(wzero))
            {
                ((pars->d)->b)[brow] = wzero;
                brow++;
            }
        }
    }
    // optional: mean constraints
    if (mmeans > 0)
    {
        for (R_len_t i = 0; i < p; i++)
        {
            double meanval = REAL_ELT(pars->colmeans, i);
            if (!R_IsNA(meanval))
            {
                ((pars->d)->b)[brow] = meanval;
                brow++;
            }
        }
    }
    // optional: augmented design
    if (!Rf_isNull(pars->w0) && gamma > 0)
    {
        for (R_len_t s = 0; s < n; s++)
        {
            ((pars->d)->b)[brow + s] = -(1.0 - gamma) * REAL_ELT(pars->w0, s);
        }
    }
    // optional: weight upper bound
    if (upr < 1)
    {
        for (R_len_t s = 0; s < n; s++)
            ((pars->d)->b)[brow + n + s] = upr;
    }

    /* populate A matrix */
    double *Jptr = REAL(pars->J);
    double *Kptr = REAL(pars->K);
    double *ortho = NULL;
    R_len_t ni = n;
    R_len_t irow = 0;
    R_len_t icol = 0;
    R_len_t *Ac = (R_len_t *)S_alloc(nnz, sizeof(R_len_t)); // COO column indices
    R_len_t *Ai = (R_len_t *)S_alloc(nnz, sizeof(R_len_t)); // COO row indices
    double *Ax = (double *)S_alloc(nnz, sizeof(double));    // COO nz values

    // COO format
    // sum(w) = 1
    for (R_len_t s = 0; s < n; s++)
    {
        Ax[s] = 1.0;
        Ai[s] = 0;
        Ac[s] = s;
    }
    irow += 1;
    // Ai * Aj = 0 (exact orthogonality)
    if (z_ortho > 0)
    {
        ortho = REAL(pars->ortho);
        R_len_t li = 0;
        for (R_len_t i = 1; i < p; i++)
        {
            for (R_len_t j = 0; j <= i; j++)
            {
                if (ortho[mat_idx(j, i, p)] != R_PosInf && ortho[mat_idx(j, i, p)] < 1e-4)
                {
                    for (R_len_t s = 0; s < n; s++)
                    {
                        Ax[ni] = Jptr[mat_idx(s, i, n)] * Jptr[mat_idx(s, j, n)]; // A[, i] * A[, j] = 0
                        Ai[ni] = irow + li;
                        Ac[ni] = s;
                        ni++;
                    }
                    li++;
                }
            }
        }
        irow += z_ortho;
    }
    // A * Y = sum(mu) * K
    icol = n;
    for (R_len_t j = 0; j < p; j++)
    {
        for (R_len_t i = 0; i < p; i++)
        {
            for (R_len_t s = 0; s < n; s++)
            {
                Ax[ni] = Jptr[mat_idx(s, i, n)]; // A[i] * Y[j]
                Ai[ni] = irow + j * p + i;
                Ac[ni] = icol + j * n + s;
                Ax[ni + 1] = -Kptr[mat_idx(i, j, p)]; // K[ij]
                Ai[ni + 1] = irow + j * p + i;
                Ac[ni + 1] = icol + n * p + s;
                ni += 2;
            }
            if (lambda > 0)
            {
                Ax[ni] = 1.0; // e_i * Y[j]
                Ai[ni] = irow + j * p + i;
                Ac[ni] = 2 * n + np + j * p + i;
                ni++;
                for (R_len_t k = 0; k < p; k++)
                {
                    Ax[ni] = -Kptr[mat_idx(i, j, p)]; // K[ij]
                    Ai[ni] = irow + j * p + i;
                    Ac[ni] = 2 * n + np + pp + k;
                    ni++;
                }
            }
            if (blocks > 0)
            {
                for (R_len_t b = 0; b < blocks; b++)
                {
                    Ax[ni] = Vptr[mat_idx(i, b, p)]; // V[i] * Y[j]
                    Ai[ni] = irow + j * p + i;
                    Ac[ni] = ntotb + j * blocks + b;
                    Ax[ni + 1] = -Kptr[mat_idx(i, j, p)]; // mu[b] * K[ij]
                    Ai[ni + 1] = irow + j * p + i;
                    Ac[ni + 1] = ntotb + p * blocks + b;
                    ni += 2;
                }
            }
        }
    }
    irow += pp;
    // strata
    if (nlvls > 0)
    {
        int *strata_id = INTEGER(VECTOR_ELT(pars->strata, 0));
        for (R_len_t s = 0; s < n; s++)
        {
            Ax[ni] = 1.0;
            Ai[ni] = irow + strata_id[s];
            Ac[ni] = s;
            ni++;
        }
        irow += nlvls;
    }
    // zero entry weights
    if (mzeros > 0)
    {
        int *Xzero_ptr = INTEGER(pars->Xzero);
        for (R_len_t i = 0; i < p; i++)
        {
            double wzero = REAL_ELT(pars->zeros, i);
            if (!R_IsNA(wzero))
            {
                for (R_len_t s = 0; s < n; s++)
                {
                    Ax[ni] = (double)(Xzero_ptr[mat_idx(s, i, n)]);
                    Ai[ni] = irow;
                    Ac[ni] = s;
                    ni++;
                }
                irow++;
            }
        }
    }
    // column means
    if (mmeans > 0)
    {
        for (R_len_t i = 0; i < p; i++)
        {
            double meanval = REAL_ELT(pars->colmeans, i);
            if (!R_IsNA(meanval))
            {
                for (R_len_t s = 0; s < n; s++)
                {
                    Ax[ni] = Jptr[mat_idx(s, i, n)];
                    Ai[ni] = irow;
                    Ac[ni] = s;
                    ni++;
                }
                irow++;
            }
        }
    }
    // w >= 0
    for (R_len_t s = 0; s < n; s++)
    {
        Ax[ni] = -1.0;
        Ai[ni] = irow + s;
        Ac[ni] = s;
        ni++;
    }
    // w <= upr
    if (upr < 1)
    {
        irow += n;
        for (R_len_t s = 0; s < n; s++)
        {
            Ax[ni] = 1.0;
            Ai[ni] = irow + s;
            Ac[ni] = s;
            ni++;
        }
    }
    // mu >= 0
    irow += n;
    for (R_len_t s = 0; s < n; s++)
    {
        Ax[ni] = -1.0;
        Ai[ni] = irow + s;
        Ac[ni] = n + np + s;
        ni++;
    }
    irow += n;
    if (lambda > 0)
    {
        // mu_lam >= 0
        for (R_len_t k = 0; k < p; k++)
        {
            Ax[ni] = -1.0;
            Ai[ni] = irow + k;
            Ac[ni] = 2 * n + np + pp + k;
            ni++;
        }
        irow += p;
    }
    if (blocks > 0)
    {
        // mu_b >= 0
        for (R_len_t b = 0; b < blocks; b++)
        {
            Ax[ni] = -1.0;
            Ai[ni] = irow + b;
            Ac[ni] = ntotb + p * blocks + b;
            ni++;
        }
        irow += blocks;
    }

    // || (2 * Y[i1], ..., 2 * Y[ik], mu[i] - w[i]) || <= mu[i] + w[i]
    icol = 0;
    for (R_len_t s = 0; s < n; s++)
    {
        // mu[i] + w[i]
        Ax[ni] = -1.0;
        Ax[ni + 1] = -1.0;
        Ai[ni] = irow;
        Ai[ni + 1] = irow;
        Ac[ni] = s;
        Ac[ni + 1] = n + np + s;
        ni += 2;
        // 2 * Y[11], ... 2 * Y[ik]
        for (R_len_t j = 0; j < p; j++)
        {
            Ax[ni + j] = -2.0;
            Ai[ni + j] = irow + 1 + j;
            Ac[ni + j] = n + j * n + s;
        }
        ni += p;
        // mu[i] - w[i]
        Ax[ni] = 1.0;
        Ax[ni + 1] = -1.0;
        Ai[ni] = irow + p + 1;
        Ai[ni + 1] = irow + p + 1;
        Ac[ni] = s;
        Ac[ni + 1] = n + np + s;
        ni += 2;
        irow += p + 2;
    }
    if (lambda > 0)
    {
        // || (2 * Y_[i1], ..., 2 * Y[ik], mu_lam[i] - lambda) || <= mu_lam[i] + lambda
        icol = 2 * n + np;
        for (R_len_t k = 0; k < p; k++)
        {
            // mu[i] + lambda
            Ax[ni] = -1.0;
            Ai[ni] = irow;
            Ac[ni] = icol + pp + k;
            ((pars->d)->b)[irow] = lambda;
            // 2 * Y[11], ... 2 * Y[ik]
            for (R_len_t j = 0; j < p; j++)
            {
                Ax[ni + 1 + j] = -2.0;
                Ai[ni + 1 + j] = irow + 1 + j;
                Ac[ni + 1 + j] = icol + j * p + k;
            }
            // mu[i] - lambda
            Ax[ni + p + 1] = -1.0;
            Ai[ni + p + 1] = irow + p + 1;
            Ac[ni + p + 1] = icol + pp + k;
            ((pars->d)->b)[irow + p + 1] = -lambda;
            ni += p + 2;
            irow += p + 2;
        }
    }
    if (blocks > 0)
    {
        // || (2 * Y[i1], ..., 2 * Y[ik], mu[i] - 1) || <= mu[i] + 1
        for (R_len_t b = 0; b < blocks; b++)
        {
            // mu[i] + 1
            Ax[ni] = -1.0;
            Ai[ni] = irow;
            Ac[ni] = ntotb + p * blocks + b;
            ((pars->d)->b)[irow] = 1.0;
            // 2 * Y[11], ... 2 * Y[ik]
            for (R_len_t j = 0; j < p; j++)
            {
                Ax[ni + 1 + j] = -2.0;
                Ai[ni + 1 + j] = irow + 1 + j;
                Ac[ni + 1 + j] = ntotb + j * blocks + b;
            }
            // mu[i] - 1
            Ax[ni + p + 1] = -1.0;
            Ai[ni + p + 1] = irow + p + 1;
            Ac[ni + p + 1] = ntotb + p * blocks + b;
            ((pars->d)->b)[irow + p + 1] = -1.0;
            ni += p + 2;
            irow += p + 2;
        }
    }

    // |M[ij]| / sqrt(M[ii] * M[jj]) <= rho
    // RSOC: || (2 * M[ij], rho * M[ii] - rho * M[jj]) || <= rho * M[ii] + rho * M[jj]
    if (q_ortho > 0)
    {
        ortho = REAL(pars->ortho);
        for (R_len_t i = 1; i < p; i++)
        {
            for (R_len_t j = 0; j < i; j++)
            {
                if (ortho[mat_idx(j, i, p)] != R_PosInf && ortho[mat_idx(j, i, p)] >= 1e-4)
                {
                    // rho * (M[ii] + M[jj])
                    double rho = ortho[mat_idx(j, i, p)];
                    for (R_len_t s = 0; s < n; s++)
                    {
                        Ax[ni] = -rho * (Jptr[mat_idx(s, i, n)] * Jptr[mat_idx(s, i, n)] + Jptr[mat_idx(s, j, n)] * Jptr[mat_idx(s, j, n)]);
                        Ai[ni] = irow;
                        Ac[ni] = s;
                        ni++;
                    }
                    // 2 * M[ij]
                    for (R_len_t s = 0; s < n; s++)
                    {
                        Ax[ni] = -2.0 * Jptr[mat_idx(s, i, n)] * Jptr[mat_idx(s, j, n)];
                        Ai[ni] = irow + 1;
                        Ac[ni] = s;
                        ni++;
                    }
                    // rho * (M[ii] - M[jj])
                    for (R_len_t s = 0; s < n; s++)
                    {
                        Ax[ni] = -rho * (Jptr[mat_idx(s, i, n)] * Jptr[mat_idx(s, i, n)] - Jptr[mat_idx(s, j, n)] * Jptr[mat_idx(s, j, n)]);
                        Ai[ni] = irow + 2;
                        Ac[ni] = s;
                        ni++;
                    }
                    irow += 3;
                }
            }
        }
    }

    // CSC format
    // column pointers
    R_len_t *ccount = (R_len_t *)S_alloc(ntot, sizeof(R_len_t));
    for (R_len_t i = 0; i < nnz; i++)
        (((pars->d)->A)->p)[Ac[i] + 1] += 1;
    for (R_len_t s = 0; s < ntot; s++)
        (((pars->d)->A)->p)[s + 1] += (((pars->d)->A)->p)[s];
    // re-assign value/row indices
    // initialize column counts
    for (R_len_t s = 0; s < ntot; s++)
        ccount[s] = (((pars->d)->A)->p)[s];
    for (R_len_t i = 0; i < nnz; i++)
    {
        icol = Ac[i];
        ni = ccount[icol];
        (((pars->d)->A)->i)[ni] = Ai[i];
        (((pars->d)->A)->x)[ni] = Ax[i];
        (ccount[icol])++;
    }
}

void A_opt_populate_K(pdata *pars, R_len_t p, R_len_t n, double upr, scs_int nlvls, R_len_t z_ortho, R_len_t q_ortho, 
                      R_len_t qsize, R_len_t mzeros, R_len_t mmeans, double lambda, R_len_t blocks)
{
    (pars->k)->z = 1 + p * p + nlvls + z_ortho + mzeros + mmeans;
    (pars->k)->l = (2 + (upr < 1)) * n + (lambda > 0) * p + blocks;
    R_len_t q1 = n + (lambda > 0) * p;
    R_len_t q2 = q1 + blocks;
    for (R_len_t s = 0; s < qsize; s++)
    {
        if(s < q1)
            ((pars->k)->q)[s] = 2 + p;
        else if (s < q2)
            ((pars->k)->q)[s] = 2 + p;
        else
            ((pars->k)->q)[s] = 3;
    }
    (pars->k)->qsize = qsize;
}

void G_opt_populate_Abc(pdata *pars, R_len_t p, R_len_t n, R_len_t nnz, scs_int mtot, scs_int ntot, double alpha, double lambda,
                        double gamma, double upr, scs_int nlvls, R_len_t z_ortho, R_len_t q_ortho, R_len_t mzeros, R_len_t mmeans)
{
    /* mixed-model SCA */
    double *Vptr = NULL;
    R_len_t blocks = 0;
    if (!Rf_isNull(pars->v_sca))
    {
        Vptr = REAL(pars->v_sca);
        blocks = Rf_ncols(pars->v_sca);
    }
    R_len_t ucol = n + (lambda > 0) * p + blocks;
    R_len_t ucolb = ucol - blocks;
    R_len_t usize = ucol * ucol;
    R_len_t brow = 1 + n * ucol + z_ortho;

    /* populate b and c vectors */
    ((pars->d)->b)[0] = 1.0;
    ((pars->d)->c)[n + 2 * usize] = -1.0;
    // optional: cost penalty
    if (!Rf_isNull(pars->f) && alpha > 0)
    {
        for (R_len_t s = 0; s < n; s++)
        {
            ((pars->d)->c)[s] = alpha * REAL_ELT(pars->f, s);
        }
    }
    // optional: strata constraints
    if (nlvls > 0)
    {
        double *prop = REAL(VECTOR_ELT(pars->strata, 1));
        for (R_len_t s = 0; s < nlvls; s++)
            ((pars->d)->b)[brow + s] = prop[s];
        brow += nlvls;
    }
    // optional: zero entry constraints
    if (mzeros > 0)
    {
        for (R_len_t i = 0; i < p; i++)
        {
            double wzero = REAL_ELT(pars->zeros, i);
            if (!R_IsNA(wzero))
            {
                ((pars->d)->b)[brow] = wzero;
                brow++;
            }
        }
    }
    // optional: mean constraints
    if (mmeans > 0)
    {
        for (R_len_t i = 0; i < p; i++)
        {
            double meanval = REAL_ELT(pars->colmeans, i);
            if (!R_IsNA(meanval))
            {
                ((pars->d)->b)[brow] = meanval;
                brow++;
            }
        }
    }
    // optional: augmented design
    if (!Rf_isNull(pars->w0) && gamma > 0)
    {
        for (R_len_t s = 0; s < n; s++)
        {
            ((pars->d)->b)[brow + s] = -(1.0 - gamma) * REAL_ELT(pars->w0, s);
        }
    }
    // optional: weight upper bounds
    if (upr < 1)
    {
        for (R_len_t s = 0; s < n; s++)
            ((pars->d)->b)[brow + n + s] = upr;
    }

    /* populate A matrix */
    double *Jptr = REAL(pars->J);
    double *ortho = NULL;
    R_len_t ni = n;
    R_len_t irow = 0;
    R_len_t icol = 0;
    R_len_t *Ac = (R_len_t *)S_alloc(nnz, sizeof(R_len_t)); // COO column indices
    R_len_t *Ai = (R_len_t *)S_alloc(nnz, sizeof(R_len_t)); // COO row indices
    double *Ax = (double *)S_alloc(nnz, sizeof(double));    // COO nz values

    // COO format
    // sum(w) = 1
    for (R_len_t s = 0; s < n; s++)
    {
        Ax[s] = 1.0;
        Ai[s] = 0;
        Ac[s] = s;
    }
    // Ai * Aj = 0 (exact orthogonality)
    irow += 1;
    if (z_ortho > 0)
    {
        ortho = REAL(pars->ortho);
        R_len_t li = 0;
        for (R_len_t i = 1; i < p; i++)
        {
            for (R_len_t j = 0; j < i; j++)
            {
                if (ortho[mat_idx(j, i, p)] != R_PosInf && ortho[mat_idx(j, i, p)] < 1e-4)
                {
                    for (R_len_t s = 0; s < n; s++)
                    {
                        Ax[ni] = Jptr[mat_idx(s, i, n)] * Jptr[mat_idx(s, j, n)]; // A[, i] * A[, j] = 0
                        Ai[ni] = irow + li;
                        Ac[ni] = s;
                        ni++;
                    }
                    li++;
                }
            }
        }
        irow += z_ortho;
    }
    // A[j] * H = sum(u) * A[i]
    icol = n;
    for (R_len_t i = 0; i < ucol; i++)
    {
        for (R_len_t j = 0; j < p; j++)
        {
            for (R_len_t s = 0; s < n; s++)
            {
                Ax[ni] = Jptr[mat_idx(s, j, n)]; // A[ij] * H[ij]
                Ai[ni] = irow + i * p + j;
                Ac[ni] = icol + i * ucol + s;
                ni++;
            }
            if (lambda > 0)
            {
                Ax[ni] = 1.0; // e[j] * H[ij]
                Ai[ni] = irow + i * p + j;
                Ac[ni] = icol + i * ucol + (n + j);
                ni++;
            }
            if (blocks > 0)
            {
                for (R_len_t b = 0; b < blocks; b++)
                {
                    Ax[ni] = Vptr[mat_idx(j, b, p)]; // v[bj] * H[bj]
                    Ai[ni] = irow + i * p + j;
                    Ac[ni] = icol + i * ucol + (ucolb + j);
                    ni++;
                }
            }
            if (i < n || (i - n == j) || i > ucolb)
            {
                for (R_len_t s = 0; s < ucol; s++)
                {
                    if (i < n)
                        Ax[ni] = -Jptr[mat_idx(i, j, n)]; // A[ij] * u[si]
                    else if (i > ucolb)
                        Ax[ni] = -Vptr[mat_idx(j, i - ucolb, p)];
                    else if (i - n == j)
                        Ax[ni] = -1.0;
                    Ai[ni] = irow + i * p + j;
                    Ac[ni] = icol + usize + i * ucol + s;
                    ni++;
                }
            }
        }
    }
    irow += n * ucol;
    // strata
    if (nlvls > 0)
    {
        int *strata_id = INTEGER(VECTOR_ELT(pars->strata, 0));
        for (R_len_t s = 0; s < n; s++)
        {
            Ax[ni] = 1.0;
            Ai[ni] = irow + strata_id[s];
            Ac[ni] = s;
            ni++;
        }
        irow += nlvls;
    }
    // zero entry weights
    if (mzeros > 0)
    {
        int *Xzero_ptr = INTEGER(pars->Xzero);
        for (R_len_t i = 0; i < p; i++)
        {
            double wzero = REAL_ELT(pars->zeros, i);
            if (!R_IsNA(wzero))
            {
                for (R_len_t s = 0; s < n; s++)
                {
                    Ax[ni] = (double)(Xzero_ptr[mat_idx(s, i, n)]);
                    Ai[ni] = irow;
                    Ac[ni] = s;
                    ni++;
                }
                irow++;
            }
        }
    }
    // column means
    if (mmeans > 0)
    {
        for (R_len_t i = 0; i < p; i++)
        {
            double meanval = REAL_ELT(pars->colmeans, i);
            if (!R_IsNA(meanval))
            {
                for (R_len_t s = 0; s < n; s++)
                {
                    Ax[ni] = Jptr[mat_idx(s, i, n)];
                    Ai[ni] = irow;
                    Ac[ni] = s;
                    ni++;
                }
                irow++;
            }
        }
    }
    // w >= 0
    for (R_len_t s = 0; s < n; s++)
    {
        Ax[ni] = -1.0;
        Ai[ni] = irow + s;
        Ac[ni] = s;
        ni++;
    }
    // w <= upr
    if (upr < 1)
    {
        irow += n;
        for (R_len_t s = 0; s < n; s++)
        {
            Ax[ni] = 1.0;
            Ai[ni] = irow + s;
            Ac[ni] = s;
            ni++;
        }
    }
    // u >= 0
    irow += n;
    for (R_len_t s = 0; s < usize; s++)
    {
        Ax[ni] = -1.0;
        Ai[ni] = irow + s;
        Ac[ni] = n + usize + s;
        ni++;
    }
    // rho <= sum(u[j])
    irow += usize;
    icol = n + usize;
    for (R_len_t s = 0; s < ucol; s++)
    {
        for (R_len_t r = 0; r < ucol; r++)
        {
            Ax[ni] = -1.0;
            Ai[ni] = irow + s;
            Ac[ni] = icol + s * ucol + r;
            ni++;
        }
        Ax[ni] = 1.0;
        Ai[ni] = irow + s;
        Ac[ni] = icol + usize;
        ni++;
    }
    // || (2 * H[ij], u[ij] - w[j]) || <= u[ij] + w[j]
    irow += ucol;
    icol = n + ucol;
    for (R_len_t i = 0; i < ucol; i++)
    {
        for (R_len_t j = 0; j < ucol; j++)
        {
            // u[ij] + w[j] (or lambda)
            Ax[ni] = -1.0;
            Ai[ni] = irow;
            Ac[ni] = icol + i * ucol + j;
            ni++;
            if (j < n)
            {
                Ax[ni] = -1.0;
                Ai[ni] = irow;
                Ac[ni] = j;
                ni++;
            }
            else 
            {
                ((pars->d)->b)[irow] = (j > ucolb) ? 1.0 : lambda;
            }
            // 2 * H[ij]
            Ax[ni] = -2.0;
            Ai[ni] = irow + 1;
            Ac[ni] = n + i * ucol + j;
            // u[ij] - w[j] (or lambda)
            Ax[ni + 1] = -1.0;
            Ai[ni + 1] = irow + 2;
            Ac[ni + 1] = icol + i * ucol + j;
            ni += 2;
            if (j < n)
            {
                Ax[ni] = 1.0;
                Ai[ni] = irow + 2;
                Ac[ni] = j;
                ni++;
            }
            else 
            {
                ((pars->d)->b)[irow + 2] = (j > ucolb) ? -1.0 : -lambda;
            }
            irow += 3;
        }
    }

    // |M[ij]| / sqrt(M[ii] * M[jj]) <= rho
    // RSOC: || (2 * M[ij], rho * M[ii] - rho * M[jj]) || <= rho * M[ii] + rho * M[jj]
    if (q_ortho > 0)
    {
        ortho = REAL(pars->ortho);
        for (R_len_t i = 1; i < p; i++)
        {
            for (R_len_t j = 0; j < i; j++)
            {
                if (ortho[mat_idx(j, i, p)] != R_PosInf && ortho[mat_idx(j, i, p)] >= 1e-4)
                {
                    // rho * (M[ii] + M[jj])
                    double rho = ortho[mat_idx(j, i, p)];
                    for (R_len_t s = 0; s < n; s++)
                    {
                        Ax[ni] = -rho * (Jptr[mat_idx(s, i, n)] * Jptr[mat_idx(s, i, n)] + Jptr[mat_idx(s, j, n)] * Jptr[mat_idx(s, j, n)]);
                        Ai[ni] = irow;
                        Ac[ni] = s;
                        ni++;
                    }
                    // 2 * M[ij]
                    for (R_len_t s = 0; s < n; s++)
                    {
                        Ax[ni] = -2.0 * Jptr[mat_idx(s, i, n)] * Jptr[mat_idx(s, j, n)];
                        Ai[ni] = irow + 1;
                        Ac[ni] = s;
                        ni++;
                    }
                    // rho * (M[ii] - M[jj])
                    for (R_len_t s = 0; s < n; s++)
                    {
                        Ax[ni] = -rho * (Jptr[mat_idx(s, i, n)] * Jptr[mat_idx(s, i, n)] - Jptr[mat_idx(s, j, n)] * Jptr[mat_idx(s, j, n)]);
                        Ai[ni] = irow + 2;
                        Ac[ni] = s;
                        ni++;
                    }
                    irow += 3;
                }
            }
        }
    }

    // CSC format
    // column pointers
    R_len_t *ccount = (R_len_t *)S_alloc(ntot, sizeof(R_len_t));
    for (R_len_t i = 0; i < nnz; i++)
        (((pars->d)->A)->p)[Ac[i] + 1] += 1;
    for (R_len_t s = 0; s < ntot; s++)
        (((pars->d)->A)->p)[s + 1] += (((pars->d)->A)->p)[s];
    // re-assign value/row indices
    // initialize column counts
    for (R_len_t s = 0; s < ntot; s++)
        ccount[s] = (((pars->d)->A)->p)[s];
    for (R_len_t i = 0; i < nnz; i++)
    {
        icol = Ac[i];
        ni = ccount[icol];
        (((pars->d)->A)->i)[ni] = Ai[i];
        (((pars->d)->A)->x)[ni] = Ax[i];
        (ccount[icol])++;
    }
}

void G_opt_populate_K(pdata *pars, R_len_t p, R_len_t n, double upr, scs_int nlvls, R_len_t z_ortho, R_len_t qsize, R_len_t mzeros, R_len_t mmeans, double lambda, R_len_t blocks)
{
    R_len_t ucol = n + (lambda > 0) * p + blocks;
    R_len_t usize = ucol * ucol;
    (pars->k)->z = 1 + n * ucol + nlvls + z_ortho + mzeros + mmeans;
    (pars->k)->l = (1 + (upr < 1)) * n + usize + ucol;
    for (R_len_t s = 0; s < qsize; s++)
        ((pars->k)->q)[s] = 3;
    (pars->k)->qsize = qsize;
}

void alias_opt_populate_Abc(pdata *pars, R_len_t p, R_len_t n, R_len_t nnz, scs_int mtot, scs_int ntot, double alpha,
                            double gamma, double upr, scs_int nlvls, R_len_t z_ortho, R_len_t q_ortho, R_len_t mzeros, R_len_t mmeans)
{
    R_len_t brow = 1 + z_ortho;
    /* populate b and c vectors */
    ((pars->d)->b)[0] = 1.0;
    R_len_t pp2 = p * (p - 1) / 2;
    for (R_len_t s = n; s < ntot; s++)
    {
        ((pars->d)->c)[s] = -1.0 / pp2;
    }
    // optional: cost penalty
    if (!Rf_isNull(pars->f) && alpha > 0)
    {
        for (R_len_t s = 0; s < n; s++)
        {
            ((pars->d)->c)[s] = alpha * REAL_ELT(pars->f, s);
        }
    }
    // optional: strata constraints
    if (nlvls > 0)
    {
        double *prop = REAL(VECTOR_ELT(pars->strata, 1));
        for (R_len_t s = 0; s < nlvls; s++)
            ((pars->d)->b)[brow + s] = prop[s];
        brow += nlvls;
    }
    // optional: zero entry constraints
    if (mzeros > 0)
    {
        for (R_len_t i = 0; i < p; i++)
        {
            double wzero = REAL_ELT(pars->zeros, i);
            if (!R_IsNA(wzero))
            {
                ((pars->d)->b)[brow] = wzero;
                brow++;
            }
        }
    }
    // optional: mean constraints
    if (mmeans > 0)
    {
        for (R_len_t i = 0; i < p; i++)
        {
            double meanval = REAL_ELT(pars->colmeans, i);
            if (!R_IsNA(meanval))
            {
                ((pars->d)->b)[brow] = meanval;
                brow++;
            }
        }
    }
    // optional: augmented design
    if (!Rf_isNull(pars->w0) && gamma > 0)
    {
        for (R_len_t s = 0; s < n; s++)
        {
            ((pars->d)->b)[brow + s] = -(1.0 - gamma) * REAL_ELT(pars->w0, s);
        }
    }
    // optional: weight upper bound
    if (upr < 1)
    {
        for (R_len_t s = 0; s < n; s++)
            ((pars->d)->b)[brow + n + s] = upr;
    }

    /* populate A matrix */
    double *Jptr = REAL(pars->J);
    double *ortho = NULL;
    R_len_t ni = n;
    R_len_t irow = 0;
    R_len_t icol = 0;
    R_len_t *Ac = (R_len_t *)S_alloc(nnz, sizeof(R_len_t)); // COO column indices
    R_len_t *Ai = (R_len_t *)S_alloc(nnz, sizeof(R_len_t)); // COO row indices
    double *Ax = (double *)S_alloc(nnz, sizeof(double));    // COO nz values

    // COO format
    // sum(w) = 1
    for (R_len_t s = 0; s < n; s++)
    {
        Ax[s] = 1.0;
        Ai[s] = 0;
        Ac[s] = s;
    }
    // Ai * Aj = 0 (exact orthogonality)
    irow += 1;
    if (z_ortho > 0)
    {
        ortho = REAL(pars->ortho);
        R_len_t li = 0;
        for (R_len_t i = 1; i < p; i++)
        {
            for (R_len_t j = 0; j < i; j++)
            {
                if (ortho[mat_idx(j, i, p)] != R_PosInf && ortho[mat_idx(j, i, p)] < 1e-4)
                {
                    for (R_len_t s = 0; s < n; s++)
                    {
                        Ax[ni] = Jptr[mat_idx(s, i, n)] * Jptr[mat_idx(s, j, n)]; // A[, i] * A[, j] = 0
                        Ai[ni] = irow + li;
                        Ac[ni] = s;
                        ni++;
                    }
                    li++;
                }
            }
        }
        irow += z_ortho;
    }
    // strata
    if (nlvls > 0)
    {
        int *strata_id = INTEGER(VECTOR_ELT(pars->strata, 0));
        for (R_len_t s = 0; s < n; s++)
        {
            Ax[ni] = 1.0;
            Ai[ni] = irow + strata_id[s];
            Ac[ni] = s;
            ni++;
        }
        irow += nlvls;
    }
    // zero entry weights
    if (mzeros > 0)
    {
        int *Xzero_ptr = INTEGER(pars->Xzero);
        for (R_len_t i = 0; i < p; i++)
        {
            double wzero = REAL_ELT(pars->zeros, i);
            if (!R_IsNA(wzero))
            {
                for (R_len_t s = 0; s < n; s++)
                {
                    Ax[ni] = (double)(Xzero_ptr[mat_idx(s, i, n)]);
                    Ai[ni] = irow;
                    Ac[ni] = s;
                    ni++;
                }
                irow++;
            }
        }
    }
    // column means
    if (mmeans > 0)
    {
        for (R_len_t i = 0; i < p; i++)
        {
            double meanval = REAL_ELT(pars->colmeans, i);
            if (!R_IsNA(meanval))
            {
                for (R_len_t s = 0; s < n; s++)
                {
                    Ax[ni] = Jptr[mat_idx(s, i, n)];
                    Ai[ni] = irow;
                    Ac[ni] = s;
                    ni++;
                }
                irow++;
            }
        }
    }
    // w >= 0
    for (R_len_t s = 0; s < n; s++)
    {
        Ax[ni] = -1.0;
        Ai[ni] = irow + s;
        Ac[ni] = s;
        ni++;
    }
    // w <= upr
    if (upr < 1)
    {
        irow += n;
        for (R_len_t s = 0; s < n; s++)
        {
            Ax[ni] = 1.0;
            Ai[ni] = irow + s;
            Ac[ni] = s;
            ni++;
        }
    }
    // u[ij] >= 0
    irow += n;
    for (R_len_t i = 1; i < p; i++)
    {
        for (R_len_t j = 0; j < i; j++)
        {
            Ax[ni] = -1.0;
            Ai[ni] = irow + j;
            Ac[ni] = n + lower_tri_idx(i, j, p);
            ni++;
        }
        irow += i;
    }
    // (|M[ij]|^2 + |u_ij|^2) / (M[ii] * M[jj]) <= rho^2
    // RSOC: || (2 * M[ij], 2 * u[ij], rho * M[ii] - rho * M[jj]) || <= rho * M[ii] + rho * M[jj]
    double rho = 1.0;
    if (q_ortho > 0)
        ortho = REAL(pars->ortho);
    for (R_len_t i = 1; i < p; i++)
    {
        for (R_len_t j = 0; j < i; j++)
        {
            if (q_ortho > 0 && ortho[mat_idx(i, j, p)] != R_PosInf && ortho[mat_idx(i, j, p)] >= 1e-4)
                rho = ortho[mat_idx(i, j, p)];
            else
                rho = 1.0;

            for (R_len_t s = 0; s < n; s++)
            {
                // rho * (M[ii] + M[jj])
                Ax[ni] = -rho * (Jptr[mat_idx(s, i, n)] * Jptr[mat_idx(s, i, n)] + Jptr[mat_idx(s, j, n)] * Jptr[mat_idx(s, j, n)]);
                Ai[ni] = irow;
                Ac[ni] = s;
                // 2 * M[ij]
                Ax[ni + 1] = -2.0 * Jptr[mat_idx(s, i, n)] * Jptr[mat_idx(s, j, n)];
                Ai[ni + 1] = irow + 1;
                Ac[ni + 1] = s;
                // rho * (M[ii] - M[jj])
                Ax[ni + 2] = -rho * (Jptr[mat_idx(s, i, n)] * Jptr[mat_idx(s, i, n)] - Jptr[mat_idx(s, j, n)] * Jptr[mat_idx(s, j, n)]);
                Ai[ni + 2] = irow + 2;
                Ac[ni + 2] = s;
                ni += 3;
            }
            // 2 * u[ij]
            Ax[ni] = -2.0;
            Ai[ni] = irow + 3;
            Ac[ni] = n + lower_tri_idx(i, j, p);
            ni++;
            irow += 4;
        }
    }

    // CSC format
    // column pointers
    R_len_t *ccount = (R_len_t *)S_alloc(ntot, sizeof(R_len_t));
    for (R_len_t i = 0; i < nnz; i++)
        (((pars->d)->A)->p)[Ac[i] + 1] += 1;
    for (R_len_t s = 0; s < ntot; s++)
        (((pars->d)->A)->p)[s + 1] += (((pars->d)->A)->p)[s];
    // re-assign value/row indices
    // initialize column counts
    for (R_len_t s = 0; s < ntot; s++)
        ccount[s] = (((pars->d)->A)->p)[s];
    for (R_len_t i = 0; i < nnz; i++)
    {
        icol = Ac[i];
        ni = ccount[icol];
        (((pars->d)->A)->i)[ni] = Ai[i];
        (((pars->d)->A)->x)[ni] = Ax[i];
        (ccount[icol])++;
    }
}

void alias_opt_populate_K(pdata *pars, R_len_t p, R_len_t n, double upr, scs_int nlvls, R_len_t z_ortho, R_len_t qsize, R_len_t mzeros, R_len_t mmeans)
{
    (pars->k)->z = 1 + nlvls + z_ortho + mzeros + mmeans;
    (pars->k)->l = (1 + (upr < 1)) * n + p * (p - 1) / 2;
    for (R_len_t s = 0; s < qsize; s++)
        ((pars->k)->q)[s] = 4;
    (pars->k)->qsize = qsize;
}
