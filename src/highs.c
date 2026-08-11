#define R_NO_REMAP

#include "cvxoptdes.h"

// #define HIGHS_AVAILABLE 1

#if HIGHS_AVAILABLE

#include "highs_c_api.h"

/* R_highs_solve function arguments */
typedef struct
{
    SEXP L;
    SEXP lwr;
    SEXP upr;
    SEXP A;
    SEXP lhs;
    SEXP rhs;
    SEXP types;
    SEXP maximum;
    SEXP offset;
    SEXP ctrl_dbl;
    SEXP ctrl_int;
    SEXP ctrl_bool;
    SEXP ctrl_str;
    SEXP write_mps;
    void *highs;
} hdata;

static SEXP C_highs_solve(void *data);

/* cleanup memory */
static void R_highs_solve_cleanup(void *data, Rboolean jump)
{
    if(data)
    {
        hdata *pars = data;
        /* free memory */
        if (pars->highs != NULL)
        {
            Highs_destroy(pars->highs);
            pars->highs = NULL; // prevent double-free
        }
    }
    if(jump)
    {
        Rprintf("HiGHS solver interrupted by user.\n");
    }
}

// /* interrupt highs from R */
// static void R_interrupt_callback(const int callback_type, const char *message, const HighsCallbackDataOut* data_out, HighsCallbackDataIn* data_in, void *user_callback_data)
// {
//     R_CheckUserInterrupt();
// }

/**
 * @brief Update HiGHS control settings
 *
 * Updates HiGHS options using values passed from R
 *
 * @param[in] ctrl named R vector containing integer settings
 * @param[in] type integer, option type: 1 for double, 2 for integer, 3 for boolean, 4 for character
 * @param[in, out] highs pointer to existing HiGHS instance
 *
 **/
static void set_highs_options(void *highs, SEXP ctrl, int type)
{
    R_len_t n = Rf_length(ctrl);
    SEXP nms = PROTECT(Rf_getAttrib(ctrl, R_NamesSymbol));
    if (!Rf_isNull(nms))
    {
        for (R_len_t i = 0; i < n; i++)
        {
            const char *option_name = CHAR(STRING_ELT(nms, i));
            if (type == 1)
            {
                Highs_setDoubleOptionValue(highs, option_name, REAL_ELT(ctrl, i));
            }
            else if (type == 2)
            {
                Highs_setIntOptionValue(highs, option_name, INTEGER_ELT(ctrl, i));
            }
            else if (type == 3)
            {
                Highs_setBoolOptionValue(highs, option_name, LOGICAL_ELT(ctrl, i));
            }
            else if (type == 4)
            {
                Highs_setStringOptionValue(highs, option_name, CHAR(STRING_ELT(ctrl, i)));
            }
#if DEBUG_INFO
            if (type == 1)
            {
                double option_val;
                Highs_getDoubleOptionValue(highs, option_name, &option_val);
                Rprintf("double option \"%s\" set to: %g\n", option_name, option_val);
            }
            else if (type == 2)
            {
                HighsInt option_val;
                Highs_getIntOptionValue(highs, option_name, &option_val);
                Rprintf("integer option \"%s\" set to: %d\n", option_name, option_val);
            }
            else if (type == 3)
            {
                HighsInt option_val;
                Highs_getBoolOptionValue(highs, option_name, &option_val);
                Rprintf("boolean option \"%s\" set to: %d\n", option_name, option_val);
            }
            else if (type == 4)
            {
                char option_val[256];
                Highs_getStringOptionValue(highs, option_name, option_val);
                Rprintf("string option \"%s\" set to: \"%s\"\n", option_name, option_val);
            }
#endif
        }
    }
    UNPROTECT(1);
}

SEXP R_highs_available(SEXP null)
{
    return Rf_ScalarLogical(TRUE);
}

/**
 * @brief Solve mixed-integer linear program using HiGHS
 *
 * @param[in] L numeric R vector, linear part of objective function
 * @param[in] lwr numeric R vector, lower bounds of variables, must be of the same length as L
 * @param[in] upr numeric R vector, upper bounds of variables, must be of the same length as L
 * @param[in] A dgCMatrix, linear part of the constraints. Rows are constraints, columns correspond to variables.
 * The number of columns must be identical to the length of L
 * @param[in] lhs numeric R vector, left hand-side of the linear constraints, with length equal to the number of rows in A
 * @param[in] rhs numeric R vector, right hand-side of the linear constraints, with length equal to the number of rows in A
 * @param[in] types integer R vector, variable types: 1 for continuous, 2 for integer, 3 for semi-continuous,
 * 4 for semi-integer, 5 for implicit integer. Must be of the same length as L
 * @param[in] maximum logical, if TRUE solve for the maximum, otherwise solve for the minimum
 * @param[in] offset (optional) numeric, offset included in the objective, defaults to 0
 * @param[in] ctrl_dbl (optional) named integer R vector, integer HiGHS control parameters
 * @param[in] ctrl_int (optional) named numeric R vector, numeric HiGHS control parameters
 * @param[in] ctrl_bool (optional) named logical R vector, boolean HiGHS control parameters
 * @param[in] ctrl_str (optional) named character R vector, string HiGHS control parameters
 * @param[in] write_mps (optional) .mps file location to export HiGHS problem  
 * @return R list
 */
SEXP R_highs_solve(SEXP L, SEXP lwr, SEXP upr, SEXP A, SEXP lhs, SEXP rhs, SEXP types, SEXP maximum, 
    SEXP offset, SEXP ctrl_dbl, SEXP ctrl_int, SEXP ctrl_bool, SEXP ctrl_str, SEXP write_mps)
{
    /* function arguments */
    hdata pars = {L, lwr, upr, A, lhs, rhs, types, maximum, offset, ctrl_dbl, ctrl_int, ctrl_bool, ctrl_str, write_mps, NULL};

    pars.highs = Highs_create();

    /* safe function call */
    SEXP ans = R_UnwindProtect(C_highs_solve, &pars, R_highs_solve_cleanup, &pars, NULL);

    return ans;
}

static SEXP C_highs_solve(void *data)
{
    /* unpack function arguments */
    hdata *pars = data;

    if (!(pars->highs))
    {
        Rf_error("failed to create HiGHS instance, verify installation at $HIGHS_HOME");
    }

    // Highs_setCallback(pars->highs, R_interrupt_callback, NULL);

    /* set highs options */
    if (!Rf_isNull(pars->ctrl_dbl))
    {
        set_highs_options(pars->highs, pars->ctrl_dbl, 1);
    }
    if (!Rf_isNull(pars->ctrl_int))
    {
        set_highs_options(pars->highs, pars->ctrl_int, 2);
    }
    if (!Rf_isNull(pars->ctrl_bool))
    {
        set_highs_options(pars->highs, pars->ctrl_bool, 3);
    }
    if (!Rf_isNull(pars->ctrl_str))
    {
        set_highs_options(pars->highs, pars->ctrl_str, 4);
    }

    /* initialize MIP */
    SEXP A_i = PROTECT(R_do_slot(pars->A, Rf_install("i")));
    SEXP A_p = PROTECT(R_do_slot(pars->A, Rf_install("p")));
    SEXP A_x = PROTECT(R_do_slot(pars->A, Rf_install("x")));

    const HighsInt ntot = Rf_length(pars->L);   // variables
    const HighsInt mtot = Rf_length(pars->lhs); // constraints
    const HighsInt nnz = Rf_length(A_x);  // non-zero entries
    const double offset = REAL_ELT(pars->offset, 0);
    const HighsInt sense = LOGICAL_ELT(pars->maximum, 0) > 0 ? kHighsObjSenseMaximize : kHighsObjSenseMinimize;
    const HighsInt a_format = kHighsMatrixFormatColwise;
 
    double *col_cost = (double *)R_alloc(ntot, sizeof(double));
    double *col_lwr = (double *)R_alloc(ntot, sizeof(double));
    double *col_upr = (double *)R_alloc(ntot, sizeof(double));
    HighsInt *integrality = (HighsInt *)R_alloc(ntot, sizeof(HighsInt));

    for (R_len_t i = 0; i < ntot; i++)
    {
        col_cost[i] = REAL_ELT(pars->L, i);
        col_lwr[i] = REAL_ELT(pars->lwr, i);
        col_upr[i] = REAL_ELT(pars->upr, i);
        integrality[i] = INTEGER_ELT(pars->types, i) - 1;
    }

    double *row_lwr = (double *)R_alloc(mtot, sizeof(double));
    double *row_upr = (double *)R_alloc(mtot, sizeof(double));

    for (R_len_t j = 0; j < mtot; j++)
    {
        row_lwr[j] = REAL_ELT(pars->lhs, j);
        row_upr[j] = REAL_ELT(pars->rhs, j);
    }


    R_len_t psize = Rf_length(A_p) - 1; // last index not needed
    HighsInt *a_start = (HighsInt *)R_alloc(psize, sizeof(HighsInt));
    HighsInt *a_index = (HighsInt *)R_alloc(nnz, sizeof(HighsInt));
    double *a_value = (double *)R_alloc(nnz, sizeof(double));
    for (R_len_t k = 0; k < psize; k++)
    {
        a_start[k] = INTEGER_ELT(A_p, k);
    }
    for (R_len_t l = 0; l < nnz; l++)
    {
        a_index[l] = INTEGER_ELT(A_i, l);
        a_value[l] = REAL_ELT(A_x, l);
    }
    UNPROTECT(3);

    int highs_status = 0;
    highs_status = Highs_passMip(
        pars->highs,
        ntot,
        mtot,
        nnz,
        a_format,
        sense,
        offset,
        col_cost,
        col_lwr,
        col_upr,
        row_lwr,
        row_upr,
        a_start,
        a_index,
        a_value,
        integrality);

#if DEBUG_INFO
    Rprintf("MIP initialization status code: %d\n", highs_status);
#endif

    if (highs_status != 0)
    {
        Rf_error("failed to initialize MIP: status %d\n", highs_status);
    }
    highs_status = Highs_run(pars->highs);

#if DEBUG_INFO
    Rprintf("MIP solver status code: %d\n", highs_status);
#endif

    /* mps export */
    if (!Rf_isNull(pars->write_mps)) 
    {
        const char *mps_filename = CHAR(STRING_ELT(pars->write_mps, 0));
        highs_status = Highs_writeModel(pars->highs, mps_filename);
#if DEBUG_INFO
        Rprintf("MPS export status code: %d\n", highs_status);
#endif
    }

    /* collect results */
    SEXP ans = PROTECT(Rf_allocVector(VECSXP, 3));
    SEXP ansnms = PROTECT(Rf_allocVector(STRSXP, 3));
    SET_STRING_ELT(ansnms, 0, Rf_mkChar("sol"));
    SET_STRING_ELT(ansnms, 1, Rf_mkChar("hi_model_status"));
    SET_STRING_ELT(ansnms, 2, Rf_mkChar("hi_obj_value"));
    Rf_setAttrib(ans, R_NamesSymbol, ansnms);

    HighsInt primal_solution_status;
    double objective_value;
    Highs_getIntInfoValue(pars->highs, "primal_solution_status", &primal_solution_status);
    Highs_getDoubleInfoValue(pars->highs, "objective_function_value", &objective_value);

    SET_VECTOR_ELT(ans, 1, Rf_ScalarInteger(Highs_getModelStatus(pars->highs)));
    SET_VECTOR_ELT(ans, 2, Rf_ScalarReal(objective_value));

    // primal solution vector
    SEXP sol = PROTECT(Rf_allocVector(REALSXP, ntot));
    double *sol_dbl = REAL(sol);
    if (primal_solution_status == 2)
    {
        Highs_getSolution(pars->highs, sol_dbl, NULL, NULL, NULL);
    }
    else
    {
        for (R_len_t i = 0; i < ntot; i++)
            sol_dbl[i] = NA_REAL;
    }
    SET_VECTOR_ELT(ans, 0, sol);
    UNPROTECT(3);

    return ans;

}

#else

SEXP R_highs_available(SEXP null)
{
    return Rf_ScalarLogical(FALSE);
}

SEXP R_highs_solve(SEXP L, SEXP lwr, SEXP upr, SEXP A, SEXP lhs, SEXP rhs, SEXP types, SEXP maximum, SEXP offset, SEXP ctrl_dbl, SEXP ctrl_int, SEXP ctrl_bool, SEXP ctrl_str, SEXP write_mps)
{
    return R_NilValue;
}

#endif  /* HIGHS_AVAILABLE */
