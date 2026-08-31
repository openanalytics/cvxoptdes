utils::globalVariables(
		c("alpha", "cost", "lambda", "orthogonal", "eta", "Zt")
)

#' Tunable SCS solver parameters
#'
#' Allow the user to tune the parameters of the SCS solver procedure.
#'
#' @param maxiter positive integer, maximum number of iterations to run.
#' @param eps_rel positive numeric, relative feasibility tolerance.
#' @param eps_abs positive numeric, absolute feasibility tolerance.
#' @param eps_infeas positive numeric, infeasibility tolerance (primal and dual).
#' @param alpha_relax numeric value in (0, 2), douglas-Rachford relaxation parameter.
#' @param rho_x positive numeric, primal scale factor.
#' @param scale positive numeric, initial dual scale factor, updated if \code{adaptive_scale} is \code{TRUE}.
#' @param normalize logical, whether to perform heuristic data rescaling.
#' @param acceleration_lookback positive integer, relative amount of memory used for Anderson acceleration.
#' More memory requires more time to compute but can give more reliable steps. If \code{acceleration_lookback < 1}, Anderson acceleration is disabled.
#' @param acceleration_interval positive integer, run Anderson acceleration every \code{acceleration_interval} iteration(s).
#' @param adaptive_scale logical, whether to heuristically adapt dual scale through the solve.
#' @param adaptive_diag_scale logical, whether to refine the metric R per row, using the primal residual profile. Requires \code{adaptive_scale} (silently disabled without it).
#' @param time_limit_secs positive numeric, time limit for solve run in seconds (can be fractional). The value 0 is interpreted as no limit. The default is 60.
#' @seealso \url{https://www.cvxgrp.org/scs/api/settings}
#' @examples
#' ## default tuning parameters
#' scs_control_dflt()
#' @return A \code{list} with components:
#' \itemize{
#' \item maxiter
#' \item eps_rel
#' \item eps_abs
#' \item eps_infeas
#' \item alpha_relax
#' \item rho_x
#' \item scale
#' \item normalize
#' \item acceleration_lookback
#' \item acceleration_interval
#' \item adaptive_scale
#' \item adaptive_diag_scale
#' \item time_limit_secs
#' }
#' with meanings as explained under 'Arguments'.
#' @export
scs_control_dflt <- function(maxiter = 100000, eps_rel = 1e-5, eps_abs = 1e-5, eps_infeas = 1e-7,
                             alpha_relax = 1.5, rho_x = 1e-6, scale = 0.1, normalize = TRUE,
                             acceleration_lookback = 10, acceleration_interval = 5,
                             adaptive_scale = TRUE, adaptive_diag_scale = TRUE, time_limit_secs = 60) {

  stopifnot(
    is.numeric(maxiter), length(maxiter) == 1, maxiter >= 1,
    is.numeric(eps_rel), length(eps_rel) == 1, eps_rel > 0,
    is.numeric(eps_abs), length(eps_abs) == 1, eps_abs > 0,
    is.numeric(eps_infeas), length(eps_infeas) == 1, eps_infeas > 0,
    is.numeric(alpha_relax), length(alpha_relax) == 1, alpha_relax > 0, alpha_relax < 2,
    is.numeric(rho_x), length(rho_x) == 1, rho_x > 0,
    is.numeric(scale), length(scale) == 1, scale > 0,
    is.logical(normalize),
    is.numeric(acceleration_lookback), length(acceleration_lookback) == 1,
    is.numeric(acceleration_interval), length(acceleration_interval) == 1, acceleration_interval >= 1,
    is.logical(adaptive_scale), is.logical(adaptive_diag_scale),
    is.numeric(time_limit_secs), length(time_limit_secs) == 1, time_limit_secs >= 0
  )

  list(maxiter = as.integer(maxiter), eps_rel = eps_rel, eps_abs = eps_abs,
       eps_infeas = eps_infeas, alpha_relax = alpha_relax, rho_x = rho_x, scale = scale,
       normalize = normalize, acceleration_lookback = max(as.integer(acceleration_lookback), 0),
       acceleration_interval = as.integer(acceleration_interval),
       adaptive_scale = adaptive_scale, adaptive_diag_scale = adaptive_scale && adaptive_diag_scale, 
       time_limit_secs = time_limit_secs)
}

parse_ctrl <- function(control = NULL, alpha = 1, lambda = 0, gamma = 1, upper = 1, criterion = "D", verbose = FALSE) {
  criterion <- match.arg(criterion, c("A", "I", "D", "G", "alias"))
  ctrl <- do.call(scs_control_dflt, args = if(!is.null(control)) as.list(control) else list())
  stopifnot(
    is.numeric(ctrl$maxiter), length(ctrl$maxiter) == 1, ctrl$maxiter >= 1,
    is.numeric(ctrl$eps_rel), length(ctrl$eps_rel) == 1, ctrl$eps_rel > 0,
    is.numeric(ctrl$eps_abs), length(ctrl$eps_abs) == 1, ctrl$eps_abs > 0,
    is.numeric(ctrl$eps_infeas), length(ctrl$eps_infeas) == 1, ctrl$eps_infeas > 0,
    is.numeric(ctrl$alpha_relax), length(ctrl$alpha_relax) == 1, ctrl$alpha_relax > 0, ctrl$alpha_relax < 2,
    is.numeric(ctrl$rho_x), length(ctrl$rho_x) == 1, ctrl$rho_x > 0,
    is.numeric(ctrl$scale), length(ctrl$scale) == 1, ctrl$scale > 0,
    is.logical(ctrl$normalize),
    is.numeric(ctrl$acceleration_lookback), length(ctrl$acceleration_lookback) == 1,
    is.numeric(ctrl$acceleration_interval), length(ctrl$acceleration_interval) == 1, ctrl$acceleration_interval >= 1,
    is.logical(ctrl$adaptive_scale), is.logical(ctrl$adaptive_diag_scale),
    is.numeric(ctrl$time_limit_secs), length(ctrl$time_limit_secs) == 1, ctrl$time_limit_secs >= 0,
    is.numeric(alpha), length(alpha) == 1,
    is.numeric(gamma), length(gamma) == 1, gamma >= 0, gamma <= 1,
    is.numeric(lambda), length(lambda) == 1, lambda >= 0,
    is.numeric(upper), length(upper) == 1, upper > 0,
    is.logical(verbose), length(verbose) == 1
  )
  ctrl_int <- c(
    match(criterion, c("A", "I", "D", "G", "alias")) - 1L,
    as.integer(ctrl$maxiter),
    as.integer(verbose),
    as.integer(ctrl$normalize),
    max(as.integer(ctrl$acceleration_lookback), 0L),
    as.integer(ctrl$acceleration_interval),
    as.integer(ctrl$adaptive_scale),
    as.integer(ctrl$adaptive_scale && ctrl$adaptive_diag_scale)
  )
  ctrl_dbl <- c(
    unlist(ctrl[c("eps_rel", "eps_abs", "eps_infeas", "alpha_relax", "rho_x", "scale", "time_limit_secs")]),
    alpha,
    lambda,
    gamma,
    upper
  )
  list(ctrl_int = ctrl_int, ctrl_dbl = ctrl_dbl)
}
