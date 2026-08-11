#define R_NO_REMAP

#include "cvxoptdes.h"

// #define USE_SPECTRAL_CONES 1

static inline int lower_tri_idx(int i, int j, int k)
{
    return i > j ? (k - 1) * j + (i - j - 1) - (j * (j - 1)) / 2 : -1;
}

static inline int mat_idx(int i, int j, int n)
{
    return i + n * j;
}

static double scs_g_crit(scs_float *xsol, SEXP J, R_len_t n, R_len_t p)
{
    const char uplo = 'L';
    const double one = 1.0, zero = 0.0;
    const char trans = 'T';

    double *X = REAL(J);
    double *M = (double *)R_alloc(p * p, sizeof(double));
    double *Xw = (double *)R_alloc(n * p, sizeof(double));
    double *w = (double *)S_alloc(n, sizeof(double));

    // M = X' * W * X
    for (R_len_t i = 0; i < n; i++)
    {
        if(xsol[i] > 0)
            w[i] = xsol[i];
        for (R_len_t k = 0; k < p; k++)
            Xw[i + n * k] = X[i + n * k] * sqrt(w[i]);
    }

    memset(M, 0., p * p * sizeof(double));
    F77_CALL(dsyrk)(&uplo, &trans, &p, &n, &one, Xw, &n, &zero, M, &p FCONE FCONE);

    // g-criterion
    double crit = crit_internal(M, X, NULL, p, n, 3, TRUE);
    return crit;
}

static SEXP C_scs_solve(void *data);

/**
 * @brief Update SCS settings
 *
 * Updates SCS settings using values passed from R
 *
 * @param[in] ctrl_int R vector containing integer settings
 * @param[in] ctrl_dbl R vector containing double settings
 * @param[in, out] stgs pointer to existing ScsSettings struct
 */
static void update_scs_settings(SEXP ctrl_int, SEXP ctrl_dbl, ScsSettings *stgs)
{
    /* New parameters for scs 3.0; see scs/src/include/glbopts.h for details*/
    scs_set_default_settings(stgs);

    stgs->max_iters = INTEGER_ELT(ctrl_int, 1);
    stgs->eps_rel = REAL_ELT(ctrl_dbl, 0);
    stgs->eps_abs = REAL_ELT(ctrl_dbl, 1);
    stgs->eps_infeas = REAL_ELT(ctrl_dbl, 2);
    stgs->alpha = REAL_ELT(ctrl_dbl, 3);
    stgs->rho_x = REAL_ELT(ctrl_dbl, 4);
    stgs->scale = REAL_ELT(ctrl_dbl, 5);
    stgs->verbose = INTEGER_ELT(ctrl_int, 2);
    stgs->normalize = INTEGER_ELT(ctrl_int, 3);
    stgs->acceleration_lookback = INTEGER_ELT(ctrl_int, 4);
    stgs->acceleration_interval = INTEGER_ELT(ctrl_int, 5);
    stgs->adaptive_scale = INTEGER_ELT(ctrl_int, 6);
    stgs->write_data_filename = (char *)NULL; // disabled
    stgs->log_csv_filename = (char *)NULL;    // disabled
    stgs->time_limit_secs = REAL_ELT(ctrl_dbl, 6);
}

/* cleanup memory */
static void R_scs_solve_cleanup(void *data, Rboolean jump)
{
    if (data)
    {
        pdata *pars = data;

        /* free memory */
        if (pars->d)
        {
            scs_free((pars->d)->b);
            scs_free((pars->d)->c);
            if ((pars->d)->A)
            {
                scs_free(((pars->d)->A)->x);
                scs_free(((pars->d)->A)->i);
                scs_free(((pars->d)->A)->p);
                scs_free((pars->d)->A);
            }
            scs_free(pars->d);
        }
        if (pars->k)
        {
            if ((pars->k)->q)
                scs_free((pars->k)->q);
            if ((pars->k)->d)
                scs_free((pars->k)->d);
            scs_free(pars->k);
        }
        if (pars->info)
            scs_free(pars->info);
        if (pars->stgs)
            scs_free(pars->stgs);
        if (pars->sol)
        {
            scs_free((pars->sol)->x);
            scs_free((pars->sol)->y);
            scs_free((pars->sol)->s);
            scs_free(pars->sol);
        }
        if (pars->w)
            scs_finish(pars->w);
    }
    if (jump)
    {
        Rprintf("SCS solver interrupted by user.\n");
    }
}

/**
 * @brief Solve second-order cone program (SOCP) using SCS
 *
 * @param[in] J numeric R matrix, design matrix for linear models or Jacobian matrix for nonlinear models
 * @param[in] f numeric R vector, design point costs
 * @param[in] ortho (optional) numeric R matrix, orthogonality constraints between design columns
 * @param[in] strata (optional) R list, contains strata definitions and strata weights
 * @param[in] zeros (optional) numeric R vector, proportion of zero entries in each design column
 * @param[in] Xzero (optional) integer R matrix, (practically) zero entries in design matrix, must be of the same size as J
 * @param[in] colmeans (optional) numeric R vector, constraints on mean values for each design column
 * @param[in] env R environment object, currently not used
 * @param[in] w0 (optional) numeric R vector, initial design weights
 * @param[in] K numeric R matrix, K matrix used to define A- and I-optimality criteria
 * @param[in] v_sca (optional) R matrix, v matrix needed for mixed-model SCA
 * @param[in] start (optional) list with x, y, s vectors used to warm-start SCS
 * @param[in] ctrl_dbl integer R vector, integer SCS control parameters
 * @param[in] ctrl_int numeric R vector, numeric SCS control parameters
 * @return R list
 */
SEXP R_scs_solve(SEXP J, SEXP f, SEXP ortho, SEXP strata, SEXP zeros, SEXP Xzero, SEXP colmeans, SEXP env, SEXP w0, SEXP K, SEXP v_sca, SEXP start, SEXP ctrl_int, SEXP ctrl_dbl)
{
    /* function arguments */
    pdata pars = {J, f, ortho, strata, zeros, Xzero, colmeans, env, w0, K, v_sca, start, ctrl_int, ctrl_dbl, NULL, NULL, NULL, NULL, NULL, NULL};

    /* safe function call */
    SEXP ans = R_UnwindProtect(C_scs_solve, &pars, R_scs_solve_cleanup, &pars, NULL);

    return ans;
}

static SEXP C_scs_solve(void *data)
{
    /* function arguments */
    pdata *pars = data;

    /* initialize parameters */
    int crit = INTEGER_ELT(pars->ctrl_int, 0);
    int crit0 = crit;
    int return_sol = INTEGER_ELT(pars->ctrl_int, 7);
    Rboolean warm_start = !Rf_isNull(pars->start);
    R_len_t p = Rf_ncols(pars->J);
    R_len_t n = Rf_nrows(pars->J);
    R_len_t np = n * p;
    R_len_t pp = p * p;
    R_len_t pp2 = p * (p + 1) / 2;
    scs_int mtot, ntot, nlvls;
    R_len_t z_ortho = 0;
    R_len_t q_ortho = 0;
    R_len_t nnz;
    R_len_t qsize = 0;
    R_len_t dsize = 0;
    R_len_t mzeros = 0;
    R_len_t mmeans = 0;
    double alpha = REAL_ELT(pars->ctrl_dbl, 7);
    double lambda = REAL_ELT(pars->ctrl_dbl, 8);
    double gamma = REAL_ELT(pars->ctrl_dbl, 9);
    double upr = REAL_ELT(pars->ctrl_dbl, 10);

    /* mixed-model SCA */
    R_len_t blocks = 0;
    if (!Rf_isNull(pars->v_sca))
    {
        blocks = Rf_ncols(pars->v_sca);
    }

    /* use D- instead of G-optimality */
    if (crit == 3 && !(Rf_isNull(pars->v_sca) && !Rf_isNull(pars->f) && alpha > 0))
    {
        crit = 2;
    }

    /* Data */
    pars->d = (ScsData *)scs_calloc(1, sizeof(ScsData));
    pars->k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));

    if (crit == 0 || crit == 1) // A-, I-optimality
    {
        mtot = (scs_int)(1 + pp + (2 + (upr < 1)) * n + (2 + p) * n); // constraints
        ntot = (scs_int)(2 * n + np);                                 // variables
        nnz = (3 + (upr < 1)) * n + 2 * n * pp + (4 + p) * n;         // non-zero values
        qsize = n;
        if (lambda > 0)
        {
            mtot += p + (2 + p) * p;
            ntot += pp + p;
            nnz += (1 + p) * pp + p + (2 + p) * p;
            qsize += p;
        }
        if (blocks > 0)
        {
            mtot += blocks + (2 + p) * blocks;
            ntot += (1 + p) * blocks;
            nnz += 2 * blocks * pp + blocks + (2 + p) * blocks;
            qsize += blocks;
        }
    }
    else if (crit == 2) // D-, unpenalized G-optimality
    {
#if USE_SPECTRAL_CONES
        mtot = (scs_int)(4 + (1 + (upr < 1)) * n + pp2);     // constraints
        ntot = (scs_int)(n + 2);                             // variables
        nnz = 3 + (2 + (upr < 1)) * n + n * pp2;             // non-zero values
        dsize = 1;
#else
        mtot = (scs_int)(1 + pp + (1 + (upr < 1)) * n + np + p + 3 * np + 3 * p); // constraints
        ntot = (scs_int)(n + np + pp2 + np + p);                                  // variables
        nnz = (2 + (upr < 1)) * n + 3 * p + (2 + p + 5) * np + pp2;               // non-zero values
        qsize = np;
        if (lambda > 0)
        {
            mtot += 4 * pp;
            ntot += 2 * pp;
            nnz += 6 * pp;
            qsize += pp;
        }
        if (blocks > 0)
        {
            mtot += 4 * blocks * p;
            ntot += 2 * blocks * p;
            nnz += blocks * pp + 5 * blocks * p;
            qsize += blocks * p;
        }
#endif
    }
    else if (crit == 3) // penalized G-optimality
    {
        R_len_t ucol = n + (lambda > 0) * p + blocks;
        R_len_t usize = ucol * ucol;
        mtot = (scs_int)(1 + n * ucol + (1 + (upr < 1)) * n + ucol + 4 * usize);  // constraints
        ntot = (scs_int)(1 + n + 2 * usize);                                      // variables
        nnz = (2 + (upr < 1)) * n + (1 + 2 * n + 2 * np) * ucol + 5 * usize;
        qsize = usize;
        if (lambda > 0)
            nnz += 2 * p * ucol;
        if (blocks > 0)
            nnz += 2 * ucol * p * blocks;
    }
    else if (crit == 4) // alias-optimality
    {
        mtot = (scs_int)(1 + (1 + (upr < 1)) * n + 5 * (pp2 - p)); // constraints
        ntot = (scs_int)(n + (pp2 - p));                           // variables
        nnz = (2 + (upr < 1)) * n + (3 * n + 2) * (pp2 - p);       // non-zero values
        qsize = pp2 - p;
    }
    else
    {
        // normally cannot be reached
        Rf_error("unknown optimality criterion");
    }

    if (!Rf_isNull(pars->strata))
    {
        nlvls = (scs_int)(Rf_length(VECTOR_ELT(pars->strata, 1)));
        mtot += nlvls;
        nnz += n;
    }
    else
    {
        nlvls = 0;
    }
    if (!Rf_isNull(pars->zeros))
    {
        for (R_len_t i = 0; i < p; i++)
            mzeros += !R_IsNA(REAL_ELT(pars->zeros, i));
        mtot += mzeros;
        nnz += n * mzeros;
    }
    if (!Rf_isNull(pars->ortho))
    {
        double *ortho = REAL(pars->ortho);
        for (R_len_t i = 1; i < p; i++)
            for (R_len_t j = 0; j < i; j++)
                if (ortho[j + p * i] != R_PosInf)
                {
                    if (ortho[j + p * i] < 1e-4)
                        z_ortho += 1;
                    else
                        q_ortho += 1;
                }
        mtot += (scs_int)(z_ortho);
        nnz += n * z_ortho;
        if (crit < 4)
        {
            mtot += (scs_int)(3 * q_ortho);
            nnz += 3 * n * q_ortho;
            qsize += q_ortho;
        }
    }
    if (!Rf_isNull(pars->colmeans))
    {
        for (R_len_t i = 0; i < p; i++)
            mmeans += !R_IsNA(REAL_ELT(pars->colmeans, i));
        mtot += mmeans;
        nnz += n * mmeans;
    }

    (pars->d)->m = mtot;
    (pars->d)->n = ntot;
    (pars->d)->P = NULL; // not used

    scs_float *c = (scs_float *)scs_calloc(ntot, sizeof(scs_float));
    scs_float *b = (scs_float *)scs_calloc(mtot, sizeof(scs_float));
    ScsMatrix *A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
    scs_float *Ax = (scs_float *)scs_calloc(nnz, sizeof(scs_float)); // CSC nz values
    scs_int *Ai = (scs_int *)scs_calloc(nnz, sizeof(scs_int));       // CSC row indices
    scs_int *Ap = (scs_int *)scs_calloc(ntot + 1, sizeof(scs_int));  // CSC column pointers

    (pars->d)->c = c;
    (pars->d)->b = b;
    A->m = mtot;
    A->n = ntot;
    A->x = Ax;
    A->i = Ai;
    A->p = Ap;
    (pars->d)->A = A;
    if(qsize)
    {
        scs_int *q = (scs_int *)scs_calloc(qsize, sizeof(scs_int));
        (pars->k)->q = q;
    }
#if USE_SPECTRAL_CONES
    if(dsize)
    {
        scs_int *d = (scs_int *)scs_calloc(dsize, sizeof(scs_int));
        (pars->k)->d = d;
    }
#endif

    if (crit == 0 || crit == 1) // A-, I-optimality
    {
        A_opt_populate_Abc(pars, p, n, nnz, mtot, ntot, alpha, lambda, gamma, upr, nlvls, z_ortho, q_ortho, mzeros, mmeans);
        A_opt_populate_K(pars, p, n, upr, nlvls, z_ortho, q_ortho, qsize, mzeros, mmeans, lambda, blocks);
    }
    else if (crit == 2) // D-, unpenalized G-optimality
    {
#if USE_SPECTRAL_CONES
        D_spec_opt_populate_Abc(pars, p, n, nnz, mtot, ntot, alpha, lambda, gamma, upr, nlvls, z_ortho, q_ortho, mzeros, mmeans);
        D_spec_opt_populate_K(pars, p, n, upr, nlvls, z_ortho, qsize, mzeros, mmeans);
#else
        D_opt_populate_Abc(pars, p, n, nnz, mtot, ntot, alpha, lambda, gamma, upr, nlvls, z_ortho, q_ortho, mzeros, mmeans);
        D_opt_populate_K(pars, p, n, upr, nlvls, z_ortho, qsize, mzeros, mmeans, lambda, blocks);
#endif
    }
    else if (crit == 3) // penalized G-optimality
    {
        G_opt_populate_Abc(pars, p, n, nnz, mtot, ntot, alpha, lambda, gamma, upr, nlvls, z_ortho, q_ortho, mzeros, mmeans);
        G_opt_populate_K(pars, p, n, upr, nlvls, z_ortho, qsize, mzeros, mmeans, lambda, blocks);
    }
    else if (crit == 4) // alias-optimality
    {
        alias_opt_populate_Abc(pars, p, n, nnz, mtot, ntot, alpha, gamma, upr, nlvls, z_ortho, q_ortho, mzeros, mmeans);
        alias_opt_populate_K(pars, p, n, upr, nlvls, z_ortho, qsize, mzeros, mmeans);
    }

    // -- return A matrix, b vector and c vector
    // SEXP ans = PROTECT(Rf_allocVector(VECSXP, 5));
    // SEXP ax = PROTECT(Rf_allocVector(REALSXP, nnz));
    // for (R_len_t s = 0; s < nnz; s++)
    //     SET_REAL_ELT(ax, s, (((pars->d)->A)->x)[s]);
    // SEXP ai = PROTECT(Rf_allocVector(INTSXP, nnz));
    // for (R_len_t s = 0; s < nnz; s++)
    //     SET_INTEGER_ELT(ai, s, (((pars->d)->A)->i)[s]);
    // SEXP ap = PROTECT(Rf_allocVector(INTSXP, ntot + 1));
    // for (R_len_t s = 0; s < ntot + 1; s++)
    //     SET_INTEGER_ELT(ap, s, (((pars->d)->A)->p)[s]);
    // SEXP bvec = PROTECT(Rf_allocVector(REALSXP, mtot));
    // for (R_len_t s = 0; s < mtot; s++)
    //     SET_REAL_ELT(bvec, s, ((pars->d)->b)[s]);
    // SEXP cvec = PROTECT(Rf_allocVector(REALSXP, ntot));
    // for (R_len_t s = 0; s < ntot; s++)
    //     SET_REAL_ELT(cvec, s, ((pars->d)->c)[s]);
    // SET_VECTOR_ELT(ans, 0, ax);
    // SET_VECTOR_ELT(ans, 1, ai);
    // SET_VECTOR_ELT(ans, 2, ap);
    // SET_VECTOR_ELT(ans, 3, bvec);
    // SET_VECTOR_ELT(ans, 4, cvec);
    // UNPROTECT(6);

        /* Solve */
        pars->sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
        pars->info = (ScsInfo *)scs_calloc(1, sizeof(ScsInfo));
        pars->stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));

        update_scs_settings(pars->ctrl_int, pars->ctrl_dbl, pars->stgs);

        if (warm_start)
        {
            SEXP startx = PROTECT(VECTOR_ELT(pars->start, 0));
            SEXP starty = PROTECT(VECTOR_ELT(pars->start, 1));
            SEXP starts = PROTECT(VECTOR_ELT(pars->start, 2));
            if(Rf_length(startx) == ntot)
            {
                scs_float *solx = (scs_float *)scs_calloc(ntot, sizeof(scs_float));
                for (R_len_t i = 0; i < ntot; i++)
                    solx[i] = REAL_ELT(startx, i);
                (pars->sol)->x = solx;
            }
            if(Rf_length(starty) == mtot)
            {
                scs_float *soly = (scs_float *)scs_calloc(mtot, sizeof(scs_float));
                for (R_len_t i = 0; i < mtot; i++)
                    soly[i] = REAL_ELT(starty, i);
                (pars->sol)->y = soly;
            }
            if(Rf_length(starts) == mtot)
            {
               scs_float *sols = (scs_float *)scs_calloc(mtot, sizeof(scs_float));
               for (R_len_t i = 0; i < mtot; i++)
                   sols[i] = REAL_ELT(starts, i);
               (pars->sol)->s = sols;
            }
            UNPROTECT(3);
        }

        pars->w = scs_init(pars->d, pars->k, pars->stgs);
        (void)scs_solve(pars->w, pars->sol, pars->info, warm_start); // discard output status

        /* Results */
        SEXP ans = PROTECT(Rf_allocVector(VECSXP, 5 + return_sol));
        SEXP ansnms = PROTECT(Rf_allocVector(STRSXP, 5 + return_sol));
        SET_STRING_ELT(ansnms, 0, Rf_mkChar("w"));
        SET_STRING_ELT(ansnms, 1, Rf_mkChar("crit_val"));
        SET_STRING_ELT(ansnms, 2, Rf_mkChar("info_int"));
        SET_STRING_ELT(ansnms, 3, Rf_mkChar("info_dbl"));
        SET_STRING_ELT(ansnms, 4, Rf_mkChar("info_chr"));
        if (return_sol)
        {
            SET_STRING_ELT(ansnms, 5, Rf_mkChar("scs_solution"));
        }
        Rf_setAttrib(ans, R_NamesSymbol, ansnms);

        // weights
        SEXP w = PROTECT(Rf_allocVector(REALSXP, n));
        for (R_len_t s = 0; s < n; s++)
        {
            SET_REAL_ELT(w, s, ((pars->sol)->x)[s]);
        }
        SET_VECTOR_ELT(ans, 0, w);

        // criterion value
        if (crit0 == 0 || crit0 == 1) // A-/I-optimality
        {
            double musum = 0.0;
            for (R_len_t s = 0; s < n; s++)
                musum += ((pars->sol)->x)[n + np + s];
            if (!Rf_isNull(pars->f) && alpha > 0)
            {
                for (R_len_t s = 0; s < n; s++)
                    musum -= alpha * REAL_ELT(pars->f, s) * ((pars->sol)->x)[s];
            }
            if(lambda > 0)
            {
                for (R_len_t k = 0; k < p; k++)
                    musum += ((pars->sol)->x)[2 * n + np + pp + k];
            }
            if(blocks > 0)
            {
                for (R_len_t b = 0; b < blocks; b++)
                    musum += ((pars->sol)->x)[ntot - blocks + b];
            }
            SET_VECTOR_ELT(ans, 1, Rf_ScalarReal(musum));
        }
        else if (crit0 == 2) // D-optimality
        {
    #if USE_SPECTRAL_CONES
            double vsum = -((pars->sol)->x)[n] / p;
    #else
            double vsum = 0.0;
            for (R_len_t j = 0; j < p; j++)
                vsum += ((pars->sol)->x)[n + np + pp2 + np + j] / p;
    #endif
            if (!Rf_isNull(pars->f) && alpha > 0)
            {
                for (R_len_t s = 0; s < n; s++)
                    vsum -= alpha * REAL_ELT(pars->f, s) * ((pars->sol)->x)[s];
            }
            SET_VECTOR_ELT(ans, 1, Rf_ScalarReal(exp(vsum)));
        }
        else if (crit0 == 3) // G-optimality
        {
            double rho = scs_g_crit((pars->sol)->x, pars->J, n, p);
            if (!Rf_isNull(pars->f) && alpha > 0)
            {
                for (R_len_t s = 0; s < n; s++)
                    rho -= alpha * REAL_ELT(pars->f, s) * ((pars->sol)->x)[s];
            }
            SET_VECTOR_ELT(ans, 1, Rf_ScalarReal(rho));
        }
        else if (crit0 == 4) // alias-optimality
        {
            double rhosum = 0;
            double rho, mij, mii, mjj;
            double *Jptr = REAL(pars->J);
            for (R_len_t i = 1; i < p; i++)
            {
                for (R_len_t j = 0; j < i; j++)
                {
                    mij = mii = mjj = 0;
                    for (R_len_t s = 0; s < n; s++)
                    {
                        mij += ((pars->sol)->x)[s] * Jptr[mat_idx(s, i, n)] * Jptr[mat_idx(s, j, n)];
                        mii += ((pars->sol)->x)[s] * Jptr[mat_idx(s, i, n)] * Jptr[mat_idx(s, i, n)];
                        mjj += ((pars->sol)->x)[s] * Jptr[mat_idx(s, j, n)] * Jptr[mat_idx(s, j, n)];
                    }
                    rho = (mij * mij) / (mii * mjj);
                    rhosum += (rho > 0 ? sqrt(rho) : 0);
                }
            }
            rhosum = 1.0 - rhosum / (double)(p * (p - 1.0) / 2.0);
            if (!Rf_isNull(pars->f) && alpha > 0)
            {
                for (R_len_t s = 0; s < n; s++)
                    rhosum -= alpha * REAL_ELT(pars->f, s) * ((pars->sol)->x)[s];
            }
            SET_VECTOR_ELT(ans, 1, Rf_ScalarReal(rhosum));
        }

        // info
        SEXP info_int = PROTECT(Rf_allocVector(INTSXP, 5));
        SET_INTEGER_ELT(info_int, 0, (pars->info)->iter);
        SET_INTEGER_ELT(info_int, 1, (pars->info)->status_val);
        SET_INTEGER_ELT(info_int, 2, (pars->info)->scale_updates);
        SET_INTEGER_ELT(info_int, 3, (pars->info)->rejected_accel_steps);
        SET_INTEGER_ELT(info_int, 4, (pars->info)->accepted_accel_steps);
        SET_VECTOR_ELT(ans, 2, info_int);

        SEXP info_dbl = PROTECT(Rf_allocVector(REALSXP, 15));
        SET_REAL_ELT(info_dbl, 0, (pars->info)->pobj);
        SET_REAL_ELT(info_dbl, 1, (pars->info)->dobj);
        SET_REAL_ELT(info_dbl, 2, (pars->info)->res_pri);
        SET_REAL_ELT(info_dbl, 3, (pars->info)->res_dual);
        SET_REAL_ELT(info_dbl, 4, (pars->info)->gap);
        SET_REAL_ELT(info_dbl, 5, (pars->info)->res_infeas);
        SET_REAL_ELT(info_dbl, 6, (pars->info)->res_unbdd_a);
        SET_REAL_ELT(info_dbl, 7, (pars->info)->res_unbdd_p);
        SET_REAL_ELT(info_dbl, 8, (pars->info)->setup_time);
        SET_REAL_ELT(info_dbl, 9, (pars->info)->solve_time);
        SET_REAL_ELT(info_dbl, 10, (pars->info)->scale);
        SET_REAL_ELT(info_dbl, 11, (pars->info)->comp_slack);
        SET_REAL_ELT(info_dbl, 12, (pars->info)->lin_sys_time);
        SET_REAL_ELT(info_dbl, 13, (pars->info)->cone_time);
        SET_REAL_ELT(info_dbl, 14, (pars->info)->accel_time);
        SET_VECTOR_ELT(ans, 3, info_dbl);

        SEXP info_chr = PROTECT(Rf_allocVector(STRSXP, 2));
        SET_STRING_ELT(info_chr, 0, Rf_mkChar((pars->info)->status));
        SET_STRING_ELT(info_chr, 0, Rf_mkChar((pars->info)->lin_sys_solver));
        SET_VECTOR_ELT(ans, 4, info_chr);

        // full scs solution (for warm-start)
        if (return_sol)
        {
            SEXP solvec = PROTECT(Rf_allocVector(VECSXP, 3));
            SEXP solnms = PROTECT(Rf_allocVector(STRSXP, 3));
            SET_STRING_ELT(solnms, 0, Rf_mkChar("x"));
            SET_STRING_ELT(solnms, 1, Rf_mkChar("y"));
            SET_STRING_ELT(solnms, 2, Rf_mkChar("s"));
            Rf_setAttrib(solvec, R_NamesSymbol, solnms);
            SEXP solx = PROTECT(Rf_allocVector(REALSXP, ntot));
            SEXP soly = PROTECT(Rf_allocVector(REALSXP, mtot));
            SEXP sols = PROTECT(Rf_allocVector(REALSXP, mtot));
            for (R_len_t i = 0; i < ntot; i++)
            {
                SET_REAL_ELT(solx, i, ((pars->sol)->x)[i]);
            }
            for (R_len_t i = 0; i < mtot; i++)
            {
                SET_REAL_ELT(soly, i, ((pars->sol)->y)[i]);
                SET_REAL_ELT(sols, i, ((pars->sol)->s)[i]);
            }
            SET_VECTOR_ELT(solvec, 0, solx);
            SET_VECTOR_ELT(solvec, 1, soly);
            SET_VECTOR_ELT(solvec, 2, sols);
            SET_VECTOR_ELT(ans, 5, solvec);
            UNPROTECT(5);
        }

        UNPROTECT(6);

    return ans;
}
