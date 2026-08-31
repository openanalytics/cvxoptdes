#' Penalized nonlinear locally optimal design
#'
#' @description
#' An R6 class for managing penalized locally optimal experiment designs for nonlinear models.
#' This class encapsulates the nonlinear model at a set of nominal parameter values, the discrete design grid,
#' optimization of an approximate continuous design, augmentation of an existing design, and optimal rounding
#' to exact integer designs.
#'
#' @section Optimal design problem:
#' The following penalized optimization problem is solved by \code{$optimize()}:
#' \deqn{w^* \ = \ \arg \max_{w} \Phi(M(w | X, v, \theta) + \lambda I) - \alpha \cdot c^T w}
#' where,
#' \itemize{
#' \item \eqn{w = \{w_1, \ldots, w_n \}} is the vector of design weights.
#' \item \eqn{c = \{c_1, \ldots, c_n \}} is the cost assigned to each design weight with scaling/penalty parameter \eqn{\alpha}.
#' \item \eqn{v = \{v_1, \ldots, v_n \}} is a vector of fixed user weights.
#' \item \eqn{M(w | X, v, \theta)} is the information matrix depending on the weights \eqn{w} and \eqn{v},
#' the design matrix \eqn{X} and the nominal parameter values \eqn{\theta}.
#' \item \eqn{\lambda I} is a diagonal matrix to regularize the information matrix with ridge penalty parameter \eqn{\lambda}.
#' \item \eqn{\Phi} is the chosen optimality criterion, generally a function of the information matrix.
#' }
#' Given the nonlinear model \code{f} and nominal parameter values \code{theta}, the information matrix (conditional on the error variance) is given by:
#' \deqn{M(w | X, v, \theta) \ \propto \ \sum_{i = 1}^n w_i v_i F(x_i, \theta) F(x_i, \theta)^T}
#' where \eqn{F(x_i, \theta)} is the gradient vector at \eqn{f(x_i, \theta)}, i.e. the \eqn{i}-th row of the Jacobian matrix.
#'
#' @section Augmented optimal design:
#' Given a vector of initial design weights \eqn{w_0} and associated contribution weight \eqn{\gamma},
#' the following augmented penalized optimization problem is solved by \code{$augment()}:
#' \deqn{w^* \ = \ \arg \max_{w} \Phi(M(w | \gamma, w_0, X, v, \theta) + \lambda I) - \alpha \cdot c^T w}
#' using the augmented information matrix,
#' \deqn{M(w | \gamma, w_0, X, v, \theta) \ = (1 - \gamma) M(w_0 | X, v, \theta) + \gamma M(w | X, v, \theta)}
#'
#' @section Stratified optimal design:
#' The \code{$strata()} method stratifies the design by requiring that individual strata \eqn{X_1, \ldots, X_S}
#' are represented proportionally with given weights \eqn{p_1, \ldots, p_S}, (\eqn{\sum_i p_i = 1}). If stratification is applied,
#' the following optimization problem is solved by \code{$optimize()}:
#' \deqn{
#'  \begin{array}{rl}
#'    w^* \ = \ & \arg \max_{w} \Phi(M(w | X, v, \theta) + \lambda I) - \alpha \cdot c^T w, \\
#'    & \mathrm{subject\ to} \ \sum_{i \in X_s} w_i = p_s, \quad \quad \mathrm{for} \ s = 1,\ldots, S
#'  \end{array}
#' }
#' using the same notation as in the \sQuote{Optimal design problem} section above. This means that the total weight
#' assigned to the design points in stratum \eqn{X_i} is equal to \eqn{p_i} for each \eqn{i = 1,\ldots,S}.
#' Similarly, if stratification is applied, \code{$augment()} solves the augmented optimization problem \emph{with} stratification
#' constraints.
#'
#' @section Sparse optimal design:
#' The \code{$sparsify()} method solves the penalized optimization problem with an extra sparseness penalty \eqn{P(w)}
#' in the objective to promote sparsity in the approximate design \eqn{w^*}. See \link{lm_design} for the available penalty functions.
#'
#' @section Optimality criteria:
#' The optimality criteria \eqn{\Phi(\cdot)} supported by \code{$optimize()}, \code{$augment()} and \code{$sparsify()} are: \emph{D}-,
#' \emph{A}-, \emph{I}-, \emph{G}-, and \emph{alias}-optimality, see \link{lm_design} for the definitions.
#'
#' @field formula a nonlinear model defined as a one-sided \link{formula}.
#' @field theta a named vector or list of nominal parameter values at which the information matrix \eqn{M(w, \theta)} is evaluated.
#' @field data a data frame or matrix containing all design points (rows) at which \code{formula} can be evaluated.
#' @field X the Jacobian matrix obtained from all design points in \code{data} and the model \code{formula}.
#' @field cost an optional numeric vector containing costs associated to each design point (row) in \code{data}.
#' @field weights an optional numeric vector containing fixed user weights associated to each design point (row) in \code{data}.
#' @field design_weights a numeric vector with optimal design weights populated after optimization with \code{$optimize()}, \code{$augment()}, or \code{$sparsify()}.
#'  The attribute \code{criterion} lists the used design criterion and the achieved optimal criterion value.
#' @field strata a list with elements \code{strata} and \code{proportions} imposing a partition of the design points into individual (disjoint) strata
#' with given weight proportions.
#' @field alpha the cost scale/penalty parameter in the penalized optimization problem, only used if \code{$cost} is available.
#' @field lambda the ridge penalty added to the diagonal of the information matrix in the optimization problem. The ridge penalty is not used
#' for the \code{alias} optimality criterion.
#' @field upper the design weight upper limit. The default is 1, corresponding to no upper limit constraint on the design weights.
#' @param data a data frame or matrix containing all design points (rows) at which \code{formula} can be evaluated.
#' @param cost an optional numeric vector containing costs associated to each design point (row) in \code{data}.
#' Also accepts a one-sided \link{formula} made up of column names present in \code{data}. If \code{NULL},
#' the standard (non-penalized) optimal design problem is solved.
#' @param weights an optional numeric vector containing fixed user weights associated to each design point (row) in \code{data}.
#' Also accepts a one-sided \link{formula} made up of column names present in \code{data}. If \code{NULL},
#' the unweighted optimal design problem is solved.
#' @param criterion character string specifying the optimality criterion to use, see also \sQuote{Details}. The following choices are supported:
#' \tabular{l}{
#' \code{"D"} D-optimality (default).\cr
#' \code{"A"} A-optimality.\cr
#' \code{"I"} I-optimality.\cr
#' \code{"G"} G-optimality.\cr
#' \code{"alias"} alias-optimality.
#' }
#' @param scs_control an optional list of control parameters to tune the SCS optimization. See \code{\link{scs_control_dflt}}
#' for the available control parameters and their default values.
#' @param verbose display verbose messages printed by the scs solver. Defaults to \code{FALSE}.
#' @param return_value return the optimal design weights. If set to \code{FALSE}, the design weights should be retrieved from the \code{design_weights} field.
#' @param alpha numeric value, cost scale/penalty parameter in the penalized optimization problem, only used if \code{$cost} is available. The default is 1.
#' @param lambda numeric value, ridge penalty added to the diagonal of the information matrix in the penalized optimization problem. The default is 0.
#' @param upper numeric value in (0, 1], design weight upper limit. The default is 1, corresponding to no constraint on the design weights.
#' @importFrom stats deriv
#' @examples
#'
#' ## Lubricant example from [Bates & Watts (1998), Appendix 1, A1.8]
#'
#' ## design grid
#' data <- expand.grid(
#'   x1 = seq(from = 0, to = 100, by = 10), ## temp. (C)
#'   x2 = seq(from = 0, to = 7, by = 0.5)   ## pressure (atm)
#' )
#' ## initialize model, nominal parameters and design grid
#' design <- nlm_design$new(
#'   formula = ~b1/(b2 + x1) + b3 * x2 + b4 * x2^2 + b5 * x2^3 +
#'               (b6 * x2 + b7 * x2^3) * exp(-x1/(b8 + b9 * x2^2)),
#'   theta = c(b1 = 1054.541, b2 = 206.546, b3 = 1.46, b4 = -0.26,
#'             b5 = 0.023, b6 = 0.401, b7 = 0.035, b8 = 57.405, b9 = -0.477),
#'   data = data
#' )
#' ## approximate I-optimal design
#' design$optimize(criterion = "I")
#' ## integer design w/ 11 design points
#' (exact_design <- design$round(m = 11, method = "optimal", seed = 1))
#' design$subset(replicates = exact_design)
#'
#' ## Penalized optimal design
#'
#' design$update(
#'   cost = ~abs(x1 - 25) / 100 + abs(x2 - 1) / 7,                  ## condition weights
#'   alpha = 0.1,
#'   upper = 0.1                                                    ## max weight per point
#' )
#' ## re-optimize
#' design$optimize(criterion = "I")
#' ## optimal rounding w/ penalization
#' design$round(m = 14, method = "optimal", seed = 1)
#'
#' ## Stratified optimal design
#'
#' design$update(
#'   cost = NULL,   ## drop cost penalties
#'   data = subset(data, x1 %in% c(0, 20, 40, 100))
#' )
#' ## partition design space by temperature conditions
#' ## by default, equal weight proportions are assigned to each condition
#' design$stratify(
#'   strata = ~x1
#' )
#' ## stratified optimal design
#' design$optimize(criterion = "I")
#' ## verify stratification
#' tapply(design$design_weights, design$data$x1, sum)
#' ## optimal rounding w/ stratification
#' design$round(m = 12, method = "optimal", seed = 1)
#'
#' @export
nlm_design <- R6Class(
  "nlm_design",
  inherit = lm_design,
  public = list(
    formula = NULL,
    theta = NULL,
    data = NULL,
    weights = NULL,
    cost = NULL,
    design_weights = NULL,
    strata = NULL,
    alpha = 1,
    lambda = 0,
    upper = 1,

    #' @description
    #' Constructor. Initializes the experimental design to optimize with a nonlinear model formula and discrete design grid.
    #' @param formula a nonlinear model defined as a one-sided \link{formula}.
    #' @param theta a named vector or list of nominal parameter values at which the information matrix \eqn{M(w, \theta)} is evaluated.
    initialize = function(formula, data, theta, weights = NULL, cost = NULL, alpha = 1, lambda = 0, upper = 1) {

      self$formula <- as.formula(formula)
      self$theta <- theta
      self$alpha <- 1
      self$lambda <- 0
      self$upper <- 1
      if(is.data.frame(data) || is.matrix(data)) {
        self$data <- as.data.frame(data)
      } else if(!is.data.frame(data)) {
        stop("'data' must be a data.frame or a matrix")
      }
      private$.X <- jac_impl(
        frm = self$formula,
        X = self$data,
        theta = self$theta,
        env = parent.frame()
      )
      if(!is.null(weights)) {
        self$weights <- weights_impl(weights, self$data, name = "weights")
      }
      if(!is.null(cost)) {
        self$cost <- weights_impl(cost, self$data, name = "cost")
      }
      if(!identical(alpha, 1)) {
        stopifnot("'alpha' must be a nonnegative numeric value" = is.numeric(alpha) && length(alpha) == 1 && alpha >= 0)
        self$alpha <- alpha
      }
      if(!identical(lambda, 0)) {
        stopifnot("'lambda' must be a nonnegative numeric value" = is.numeric(lambda) && length(lambda) == 1 && lambda >= 0)
        self$lambda <- lambda
      }
      if(!identical(upper, 1)) {
        stopifnot("'upper' must be a positive numeric value" = is.numeric(upper) && length(upper) == 1 && upper > 0)
        self$upper <- min(upper, 1)
      }

    },

    #' @description
    #' Updates the fields of an existing object.
    #' Inputs should be named \code{formula}, \code{data}, \code{theta}, \code{weights}, \code{cost}, \code{alpha}, \code{lambda} and/or \code{upper}.
    #' @param ... named arguments \code{formula}, \code{data}, \code{theta}, \code{weights}, \code{cost}, \code{alpha}, \code{lambda} and/or \code{upper}.
    update = function(...) {

      dots <- list(...)

      if(is.element("formula", names(dots)) && !is.null(dots$formula)) {
        self$formula <- as.formula(dots$formula)
        ## cascade updates
        if(!is.element("data", names(dots))) {
          dots$data <- self$data
        }
      }
      if(is.element("theta", names(dots)) && !is.null(dots$theta)) {
        self$theta <- dots$theta
        ## cascade updates
        if(!is.element("data", names(dots))) {
          dots$data <- self$data
        }
      }
      if(is.element("data", names(dots)) && !is.null(dots$data)) {
        if(is.data.frame(dots$data) || is.matrix(dots$data)) {
          self$data <- as.data.frame(dots$data)
        } else if(!is.data.frame(dots$data)) {
          stop("'data' must be a data.frame or a matrix")
        }
        private$.X <- jac_impl(
          frm = self$formula,
          X = self$data,
          theta = self$theta
        )
        ## cascade updates
        if(!is.element("cost", names(dots)) && inherits(attr(self$cost, "formula"), "formula")) {
          dots$cost <- attr(self$cost, "formula")
        }
        if(!is.element("weights", names(dots)) && inherits(attr(self$weights, "formula"), "formula")) {
          dots$weights <- attr(self$weights, "formula")
        }
      }
      if(is.element("weights", names(dots))) {
        if(is.null(dots$weights)) {
          self$weights <- NULL
        } else {
          self$weights <- weights_impl(dots$weights, self$data, name = "weights")
        }
      }
      if(is.element("cost", names(dots))) {
        if(is.null(dots$cost)) {
          self$cost <- NULL
        } else {
          self$cost <- weights_impl(dots$cost, self$data, name = "cost")
        }
      }
      if(is.element("alpha", names(dots))) {
        stopifnot(
          "'alpha' must be a nonnegative numeric value" = is.numeric(dots$alpha) && length(dots$alpha) == 1 && dots$alpha >= 0
        )
        self$alpha <- dots$alpha
      }
      if(is.element("lambda", names(dots))) {
        stopifnot(
          "'lambda' must be a nonnegative numeric value" = is.numeric(dots$lambda) && length(dots$lambda) == 1 && dots$lambda >= 0
        )
        self$lambda <- dots$lambda
      }
      if(is.element("upper", names(dots))) {
        stopifnot(
          "'upper' must be a positive numeric value" = is.numeric(dots$upper) && length(dots$upper) == 1 && dots$upper > 0
        )
        self$upper <- min(dots$upper, 1)
      }
      return(invisible(self))
    },

    #' @description
    #' Stratify the design by requiring that individual strata are represented proportionally. Call with \code{strata = NULL} to
    #' unset previously defined \code{strata}.
    #' @param strata a \link{factor} variable of the same length as the number of rows in \code{data}. The \link{factor} levels determine
    #' the partition of the design points into individual (disjoint) strata. Also accepts a one-sided \link{formula} made up of
    #' column names present in \code{data}, in which case the strata are defined by the \link{interaction} between terms in the
    #' \link{model.frame} obtained from the \link{formula}.
    #' @param proportions an optional numeric vector containing the weight proportions of the individual strata. Must be of the same length
    #' as the number of factor levels in \code{strata}. If \code{NULL}, each stratum is assigned the same weight proportion.
    stratify = function(strata, proportions = NULL) {

      super$stratify(
        strata = strata,
        proportions = proportions
      )

    },

    #' @description
    #' Optimizes the augmented penalized approximate design given a vector of initial design weights
    #' @param design_weights a numeric vector with initial design weights, must be the same length as the number of design points in \code{data}.
    #' Also accepts a one-sided \link{formula} made up of parameters in \code{theta} and column names in \code{data}.
    #' @param gamma numeric contribution weight of the initial design in the augmented design problem, see also \sQuote{Augmented optimal design}
    #' above. Defaults to 0.5.
    augment = function(design_weights, gamma = 0.5, criterion = c("D", "A", "I", "G", "alias"), scs_control = scs_control_dflt(), verbose = FALSE, return_value = TRUE) {

      if(inherits(design_weights, "formula")) {
        design_weights <- tryCatch(eval(design_weights[[2L]], envir = c(as.list(self$theta), as.list(self$data))))
        if(inherits(design_weights, "error"))
          stop(sprintf("failed to evaluate 'design_weights' formula: %s", design_weights$message))
      }

      super$augment(
        design_weights = design_weights,
        gamma = gamma,
        criterion = criterion,
        scs_control = scs_control,
        verbose = verbose,
        return_value = return_value
      )

    },

    #' @description
    #' Optimizes the penalized approximate design based on the selected criterion.
    optimize = function(criterion = c("D", "A", "I", "G", "alias"), scs_control = scs_control_dflt(), verbose = FALSE, return_value = TRUE) {

      super$optimize(
        criterion = criterion,
        scs_control = scs_control,
        verbose = verbose,
        return_value = return_value
      )

    },

    #' @description
    #' Sparsify the penalized approximate design by iterative reweighting of the optimization problem.
    #' @details
    #' The sparseness cost penalties are scaled by the same \code{alpha} tuning parameter as for a user-specific cost vector, see \code{\link{scs_control_dflt}}.
    #' If a cost vector \code{c_1,\ldots,c_n} associated to the design weights is already specified. The sparseness penalty costs are multiplied
    #' by the existing cost vector.
    #' @param max_iter positive integer, maximum number of reweighted optimizations/iterations to run.
    #' @param reweight_fun character reweighting function, one of \code{"entropy", "log-sum", "power"}, see \link{lm_design} for the definitions.
    #' @param eps positive numeric value to ensure finite weights in reweighting functions.
    #' @param p postive numeric power with \eqn{0 < p < 1}. Only used if \code{reweight_fun = "power"}.
    #' @param tol numeric tolerance used as a threshold for mean absolute relative difference in design weights between iterations to stop early.
    #' @param design_weights (optional) vector of design weights used as initial approximate design.
    #' @param show_progress logical, display a progress bar during the iterative reweighting procedure. Defaults to \code{TRUE}.
    sparsify = function(max_iter, criterion = c("D", "A", "I", "G", "alias"), reweight_fun = c("entropy", "log-sum", "power"), eps = .Machine$double.eps,
                        tol = .Machine$double.eps^(1/4), p = 0.5, design_weights = NULL, scs_control = scs_control_dflt(),
                        return_value = TRUE, show_progress = TRUE) {

      super$sparsify(
        max_iter = max_iter,
        criterion = criterion,
        reweight_fun = reweight_fun,
        eps = eps,
        tol = tol,
        p = p,
        design_weights = design_weights,
        recenter = FALSE,
        scs_control = scs_control,
        show_progress = show_progress,
        return_value = return_value
      )

    },

    #' @description
    #' Round penalized approximate design to an exact (integer) design.
    #' By default, the approximate design is initialized from the \code{design_weights} field, which is the optimal (approximate) design calculated by \code{$optimize()}.
    #' For method \code{"optimal"}, ridge and cost penalties and upper limit constraints present in \code{$lambda}, \code{$cost} and \code{$upper}
    #' (if specified) are included in the criterion optimized by the point-exchange algorithm.
    #' To override the initial approximate (optimal) design, a vector of design weights can be passed using the \code{design_weights} argument.
    #' The vector must be the same length as the number of design points in \code{data}. If \code{design_weights} is set to \code{FALSE},
    #' the initial design is generated by random sampling.
    #'
    #' @param m integer number of design points included in the exact design.
    #' @param method the method used to calculate the exact design. The following choices are supported:
    #' \tabular{l}{
    #'  \code{"optimal"} optimality criterion point-exchange algorithm (default).\cr
    #'  \code{"efficient"} efficient rounding method of Pukelsheim & Reider (1992).\cr
    #'  \code{"dist-sum"} maximize total pairwise (Euclidean) distance between design points.\cr
    #'  \code{"dist-min"} maximize minimum pairwise (Euclidean) distance between design points.
    #' }
    #' @param max_iter integer number of maximum iterations in point-exchange algorithm. Only used when
    #' \code{method = "optimal"}, \code{method = "dist-sum"} or \code{method = "dist-min"}.
    #' @param n_repeats integer number of (random) repeats of the rounding/exchange algorithm. The exact design
    #' with the highest criterion value is returned. Defaults to 10.
    #' @param replicates logical, allow replicate design points in the exact design, defaults to \code{TRUE}.
    #' Can also be an integer number, limiting the maximum number of allowed replicates per design point.
    #' This argument overrides the \code{$upper} value constraint if existing.
    #' @param tol numeric tolerance used as lower threshold to continue iterations in point-exchange algorithm.
    #' @param seed (optional) random seed used to generate jitter to break ties in rounding the initial exact design.
    #' @param design_weights (optional) vector of design weights used as initial approximate design, with length equal
    #' to the number of available design points. If \code{TRUE} (default), the approximate design is initialized from the
    #' \code{design_weights} field. If \code{FALSE}, the initial approximate design is generated by random sampling.
    #' @param augment_design (optional) integer design to augment, with length equal to the number of available design points.
    #' @param show_progress logical, display a progress bar during the rounding procedure (only if
    #' \code{n_repeats > 1}). Defaults to \code{TRUE}.
    #' @param ... additional arguments to tune algorithm settings:
    #' \itemize{
    #' \item \itemize{
    #' \item \code{criterion}, character string specifying the optimality criterion for the rounding procedure in case
    #' \code{design_weights} has no \code{"criterion"} attribute, or to override the criterion specified in the \code{"criterion"} attribute.
    #' If specified, must be one of \code{"A"}, \code{"I"}, \code{"D"}, \code{"G"} or \code{"alias"}.}
    #' }
    round = function(m, method = c("optimal", "efficient", "dist-sum", "dist-min"), max_iter = 100, n_repeats = 10, replicates = TRUE,
                     tol = .Machine$double.eps^(1/3), seed = NULL, design_weights = TRUE, augment_design =  NULL, show_progress = TRUE, ...) {

      method <- match.arg(method, c("optimal", "efficient", "dist-sum", "dist-min"))

      super$round(
        m = m,
        method = method,
        max_iter = max_iter,
        n_repeats = n_repeats,
        replicates = replicates,
        tol = tol,
        seed = seed,
        design_weights = design_weights,
        augment_design = augment_design,
        show_progress = show_progress,
        ...
      )

    },

    #' @description
    #' Evaluate optimality criterion of a penalized approximate or rounded (exact) design.
    #' The approximate design is initialized from the \code{design_weights} field, which is the optimal (approximate) design calculated by \code{$optimize()}.
    #' To override the default approximate (optimal) design, a vector of design weights can be passed using the \code{design_weights} argument.
    #' @param design_weights (optional) vector of design weights for which to evaluate optimality criterion.
    #' @param lambda logical, whether to include the ridge penalty in the calculated optimality criterion (see \sQuote{Optimality criteria}). Only relevant if \code{$lambda} is positive.
    #' @param alpha logical, whether to include cost penalties in the calculated optimality criterion (see \sQuote{Optimality criteria}). Only relevant if cost penalties are present
    #' in the \code{$cost} field.
    #' @return numeric criterion value
    crit = function(design_weights = NULL, criterion = c("D", "A", "I", "G", "alias"), lambda = TRUE, alpha = TRUE) {

      super$crit(
        design_weights = design_weights,
        criterion = criterion,
        lambda = lambda,
        alpha = alpha
      )

    },

    #' @description
    #' Calculates the variance-covariance matrix of the model parameters based on an approximate or rounded (exact) design \eqn{(w_1, \ldots, w_n)} according to:
    #' \deqn{\Sigma \ = \ (F(X, \theta)^T \Omega F(X, \theta))^{-1}}
    #' where \eqn{\Omega} is \eqn{\mathrm{diag}(\omega_1, \ldots, \omega_n)} with \eqn{\omega_i = w_i v_i}, and \eqn{F(X, \theta)} is the Jacobian matrix at the nominal parameter values.
    #' The design weights are initialized from the \code{design_weights} field, which is the optimal (approximate) design calculated by \code{$optimize()}.
    #' Use the \code{m} argument to evaluate the covariance matrix for a specific number of design points, in which case the sum of the design weights is scaled to \code{m}.
    #' To override the default approximate (optimal) design weights, a vector of design weights can be passed using the \code{design_weights} argument.
    #' The vector must be the same length as the number of design points in \code{data}. If the vector of design weights is already an integer design, leave
    #' \code{m} unspecified, in which case the design is not rescaled.
    #' @param m (optional) integer number of design points at which to evaluate the estimated coviarance matrix. If missing, the
    #' design weights remain unscaled, which for the default approximate (optimal) design is equivalent to \code{m = 1}.
    #' @param design_weights (optional) vector of design weights for which to evaluate the estimated covariance matrix.
    #' @param sigma2 numeric error variance used to scale the covariance matrix. Defaults to 1, returning the unscaled covariance matrix.
    vcov = function(m, design_weights = NULL, sigma2 = 1) {

      super$vcov(
        m = m,
        design_weights = design_weights,
        sigma2 = sigma2
      )

    },

    #' @description
    #' Calculates the standard error of prediction or the square root prediction variance across the design grid based on an approximate or rounded (exact)
    #' design \eqn{(w_1, \ldots, w_n)} according to:
    #' \deqn{\mathrm{SEP}(x_i) \ = \ \sigma \sqrt{F(x_i, \theta)^T \Sigma F(x_i, \theta)}, \quad \mathrm{for } i = 1,\ldots,n}
    #' where \eqn{\Sigma} is the covariance matrix as returned by \code{$vcov()}.
    #' Use the \code{newdata} argument to evaluate the standard error of prediction at design points not in the candidate design grid defined by \code{$data}.
    #' The design weights are initialized from the \code{design_weights} field, which is the optimal (approximate) design calculated by \code{$optimize()}.
    #' To override the default approximate (optimal) design weights, a vector of design weights can be passed using the \code{design_weights} argument.
    #' The vector must be the same length as the number of design points in \code{data}.
    #' @param m (optional) integer number of design points at which to evaluate the estimated coviarance matrix. If missing, the
    #' design weights remain unscaled, which for the default approximate (optimal) design is equivalent to \code{m = 1}.
    #' @param newdata (optional) data frame or matrix with design points (rows) at which to return the standard error of prediction and which
    #' can be evaluated by \code{$formula}.
    #' @param design_weights (optional) vector of design weights for which to evaluate the estimated covariance matrix.
    #' @param sigma2 numeric error variance used to scale the standard errors. Defaults to 1, returning the standardized standard error of prediction.
    sep = function(m, newdata, design_weights = NULL, sigma2 = 1) {
      if(!missing(newdata)) {
        if(is.data.frame(newdata) || is.matrix(newdata)) {
          newdata <- as.data.frame(newdata)
        } else if(!is.data.frame(newdata)) {
          stop("'newdata' must be a data.frame or a matrix")
        }
        xpred <- jac_impl(
          frm = self$formula,
          X = newdata,
          theta = self$theta,
          env = parent.frame()
        )
      } else {
        xpred <- private$.X
      }
      Sigma <- self$vcov(m, design_weights, sigma2)
      ses <- sqrt(rowSums(xpred * (xpred %*% Sigma)))
      return(ses)
    },

    #' @description
    #' Calculates the correlation matrix of the columns of the Jacobian matrix based on an approximate or rounded (exact) design \eqn{(w_1, \ldots, w_n)} according to:
    #' \deqn{R \ = \ D^{-1/2}(J_c^T \Omega J_c)D^{-1/2}}
    #' where \eqn{\Omega} is \eqn{\mathrm{diag}(\omega_1, \ldots, \omega_n)} with \eqn{\omega_i = w_i v_i}, and \eqn{J_c} is the Jacobian matrix \eqn{J = F(X, \theta)} centered by the (weighted)
    #' column means:
    #' \deqn{\frac{\sum_{i=1}^n \omega_i J_{ij}}{\sum_{i=1}^n \omega_i}, \quad j = 1,\ldots, p}
    #' Furthermore, \eqn{D^{-1/2}} is \eqn{\mathrm{diag}(1/\sqrt{D_1},\ldots, 1/\sqrt{D_p})}, with \eqn{D_1,\ldots,D_p} the diagonal elements of \eqn{J_c^T \Omega J_c}.
    #' If \code{center = FALSE}, the uncentered correlation matrix is returned, which corresponds to the scale-invariant (unit diagonal) version of the information matrix \eqn{M(w | X, v, \theta)}.
    #' The design weights are initialized from the \code{design_weights} field, which is the optimal (approximate) design calculated by \code{$optimize()}.
    #' To override the default approximate (optimal) design weights, a vector of design weights can be passed using the \code{design_weights} argument.
    #' The vector must be the same length as the number of design points in \code{data}.
    #' @param design_weights (optional) vector of design weights for which to evaluate the information matrix.
    #' @param center logical, whether to center the Jacobian columns when evaluating the correlation matrix, defaults to \code{TRUE}.
    corr = function(design_weights = NULL, center = TRUE) {

      super$corr(
        design_weights = design_weights,
        center = center
      )

    }

  ),

  active = list(
    X = function(value) {
      if (missing(value)) {
        private$.X
      } else {
        stop("To modify the Jacobian matrix, update the data or model formula with `$update()`.", call. = FALSE)
      }
    }
  )
)

jac_impl <- function(frm, X, theta, env) {

  ## function call
  fn <- function(par, .data) eval(frm[[2L]], envir = c(as.list(par), as.list(.data)))
  fcall <- tryCatch(fn(theta, X[1L, , drop = FALSE]), error = function(err) err)
  if(inherits(fcall, "error"))
    stop(sprintf("failed to evaluate 'formula' at nominal parameter values: %s", fcall$message))
  if(!is.numeric(fcall) || !identical(length(fcall), 1L))
    stop("'formula' failed to return a numeric value at nominal parameter values")
  if(any(is.na(fcall)))
    stop("missing values returned by 'formula' at nominal parameter values")

  ## jac call
  if(!is.null(attr(fcall, "gradient"))) {
    jac <- function(par, .data) attr(fn(par, .data), "gradient")
  } else {
    jac <- NULL
    exprjac <- tryCatch(stats::deriv(frm[[2L]], namevec = names(theta)), error = function(err) err)
    if(inherits(exprjac, "error")) {
      warning(sprintf("failed to symbolically derive 'formula': %s", exprjac$message))
    } else if(is.expression(exprjac)){
      jac <- function(par, .data) {
        grad <- eval(exprjac, envir = c(as.list(par), .data))
        attr(grad, "gradient")
      }
    }
  }

  if(is.function(jac)) {
    J <- tryCatch(jac(theta, X), error = function(err) err)
    if(inherits(J, "error"))
      stop(sprintf("failed to evaluate jacobian at nominal parameter values: %s", J$message))
    if(!is.numeric(J) || !is.matrix(J) || !identical(dim(J), c(nrow(X), length(theta))))
      stop("Jacobian failed to return a numeric matrix of the expected dimensions at nominal parameter values")
    if(any(is.na(J)))
      stop("missing values detected in jacobian at nominal parameter values")
  } else {
    stop("no jacobian matrix available")
  }

  return(J)

}
