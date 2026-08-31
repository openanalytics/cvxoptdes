#' Penalized generalized linear model locally optimal design
#'
#' @description
#' An R6 class for managing penalized locally optimal experiment designs for generalized linear models.
#' This class encapsulates the generalized linear model at a set of nominal parameter values, the discrete design grid, 
#' optimization of an approximate continuous design, augmentation of an existing design, and optimal rounding 
#' to exact integer designs.
#' 
#' @section Optimal design problem:
#' The following penalized optimization problem is solved by \code{$optimize()}:
#' \deqn{w^* \ = \ \arg \max_{w} \Phi(M(w | X, v, \theta, F, g) + \lambda I) - \alpha \cdot c^T w}
#' where,
#' \itemize{
#' \item \eqn{w = \{w_1, \ldots, w_n \}} is the vector of design weights.
#' \item \eqn{c = \{c_1, \ldots, c_n \}} is the cost assigned to each design weight with scaling/penalty parameter \eqn{\alpha}.
#' \item \eqn{v = \{v_1, \ldots, v_n \}} is a vector of fixed user weights.
#' \item \eqn{M(w | X, v, \theta, F, g)} is the information matrix depending on the weights \eqn{w} and \eqn{v},
#' the design matrix \eqn{X}, the nominal parameter values \eqn{\theta}, the chosen (exponential) distribution family \eqn{F} and the link function
#' \eqn{g}.
#' \item \eqn{\Phi} is the chosen optimality criterion, generally a function of the information matrix.
#' }
#' In standard canonical form, the probability density (or mass) function of the exponential family \eqn{F} is:
#' \deqn{
#' f(y\ |\ \zeta, \phi) = \exp\left( \frac{y\zeta - b(\zeta)}{a(\phi)} + c(y, \phi) \right)
#' }
#' The cumulant function \eqn{b(\zeta)} admits \eqn{\mu = E[Y] = b'(\zeta)} and \eqn{V(\mu) = b''(\zeta)}, the variance function expressed as a
#' function of the mean. Given the link function \eqn{g}, the systematic component is \eqn{g(\mu) = \eta = x^T \theta}.
#' The information matrix is given by:
#' \deqn{
#' M(w | X, v, \theta, F, g) \ = \ \sum_{i = 1}^n w_i v_i u(x_i, \theta) x_i x_i^T
#' }
#' with weights \eqn{u(x_i, \theta) = \frac{1}{V(\mu_i)[g'(\mu_i)]^2}} causing the information matrix to depend on the parameter values \eqn{\theta},
#' unlike the information matrix \eqn{M(w | X, v)} for a standard linear model.
#'
#' @section Augmented optimal design:
#' Given a vector of initial design weights \eqn{w_0} and associated contribution weight \eqn{\gamma},
#' the following augmented penalized optimization problem is solved by \code{$augment()}:
#' \deqn{w^* \ = \ \arg \max_{w} \Phi(M(w | \gamma, w_0, X, v, \theta, F, g) + \lambda I) - \alpha \cdot c^T w}
#' using the augmented information matrix,
#' \deqn{M(w | \gamma, w_0, X, v, \theta, F, g) \ = (1 - \gamma) M(w_0 | X, v, \theta, F, g) + \gamma M(w | X, v, \theta, F, g)}
#'
#' @section Stratified optimal design:
#' The \code{$strata()} method stratifies the design by requiring that individual strata \eqn{X_1, \ldots, X_S}
#' are represented proportionally with given weights \eqn{p_1, \ldots, p_S}, (\eqn{\sum_i p_i = 1}). If stratification is applied,
#' the following optimization problem is solved by \code{$optimize()}:
#' \deqn{
#'  \begin{array}{rl}
#'    w^* \ = \ & \arg \max_{w} \Phi(M(w | X, v, \theta, F, g) + \lambda I) - \alpha \cdot c^T w, \\
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
#' @field formula a linear model defined as a one-sided \link{formula}.
#' @field theta a named vector or list of nominal parameter values at which the information matrix \eqn{M(w | X, \theta, F, g)} is evaluated.
#' @field family a description of the model distribution family and link function. See \link{family} for details of family functions.
#' @field data a data frame or matrix containing all design points (rows) at which \code{formula} can be evaluated.
#' @field X the design matrix obtained from all design points in \code{data} and the model \code{formula}.
#' @field weights an optional numeric vector containing fixed user weights associated to each design point (row) in \code{data}.
#' @field cost an optional numeric vector containing costs associated to each design point (row) in \code{data}.
#' @field design_weights a numeric vector with optimal design weights populated after optimization with \code{$optimize()}, \code{$augment()}, or \code{$sparsify()}.
#'  The attribute \code{criterion} lists the used design criterion and the achieved optimal criterion value.
#' @field strata a list with elements \code{strata} and \code{proportions} imposing a partition of the design points into individual (disjoint) strata
#' with given weight proportions.
#' @field alpha the cost scale/penalty parameter in the penalized optimization problem, only used if \code{$cost} is available.
#' @field lambda the ridge penalty added to the diagonal of the information matrix in the optimization problem. The ridge penalty is not used
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
#' @importFrom stats family gaussian
#' @seealso \code{\link[stats]{family}}
#' @references Lewis, S. L., Montgomery D. C., and Myers R. H. 2001. \emph{Examples of Designed Experiments with Nonnormal Responses}. Journal of Quality Technology 3 (33): 265–78.
#' @examples
#'
#' ## Non-conforming tiles experiment
#' ## Example 3 from [Lewis et al. (2001)]
#'
#' ## design grid
#' data <- expand.grid(
#'   x1 = c(-1, 1),
#'   x2 = c(-1, 1),
#'   x3 = c(-1, 1),
#'   x4 = c(-1, 1),
#'   x5 = c(-1, 1),
#'   x6 = c(-1, 1),
#'   x7 = c(-1, 1)
#' )
#' ## initialize model, design grid, family and nominal parameters
#' design <- glm_design$new(
#'   formula = ~x1 + x2 + x3 + x4 + x5 + x6 + x7,
#'   data = data,
#'   family = binomial(link = "logit"),
#'   theta = c(-1.419, -0.578, 0.109, 0.266, -0.262, 0.434, -0.633, 0.425)
#' )
#' ## approximate D-optimal design
#' design$optimize(criterion = "D")
#' ## integer design w/ 8 design points
#' (exact_design <- design$round(m = 8, method = "optimal", n_repeats = 25, seed = 1))
#' design$subset(replicates = exact_design)
#'
#' ## Defect semiconductor wafers
#' ## Example 5 from [Lewis et al. (2001)]
#'
#' ## design grid
#' data <- rbind(
#'   expand.grid(x1 = -1:1, x2 = -1:1, x3 = -1:1),
#'   data.frame(
#'     x1 = c(-1.732, 1.732, 0, 0, 0, 0),
#'     x2 = c(0, 0, -1.732, 1.732, 0, 0),
#'     x3 = c(0, 0, 0, 0, -1.732, 1.732)
#'   )
#' )
#' ## initialize model, design grid, family and nominal parameters
#' design <- glm_design$new(
#'   formula = ~(x1 + x2 + x3)^2 + I(x1^2) + I(x2^2) + I(x3^2),
#'   data = data,
#'   family = poisson(link = "log"),
#'   theta = c(1.9293, 0.6103, 0.7028, -0.3531, 0, 0.2471, 0, -0.5209, 0, 0)
#' )
#' ## approximate D-optimal design
#' design$optimize(criterion = "D")
#' ## integer design w/ 18 design points
#' (exact_design <- design$round(m = 18, method = "optimal", seed = 1))
#' design$subset(replicates = exact_design)
#'
#' @export
glm_design <- R6Class(
  "glm_design",
  inherit = lm_design,
  public = list(
    formula = NULL,
    theta = NULL,
    family = NULL,
    data = NULL,
    weights = NULL,
    cost = NULL,
    design_weights = NULL,
    strata = NULL,
    alpha = 1,
    lambda = 0,
    upper = 1,

    #' @description
    #' Constructor. Initializes the experimental design to optimize with a linear model formula, discrete design grid, nominal parameter values
    #' and distribution family with link function.
    #' @param formula a linear model defined as a one-sided \link{formula}.
    #' @param theta a named or unnamed vector or list of nominal parameter values at which the information matrix \eqn{M(w | X, \theta, F, g)} is evaluated.
    #' @param family a \link{family} function describing the error distribution and link function in the model. Defaults to the \link{gaussian}
    #' family with \code{"identity"} link function.
    #' @param contrasts optional contrasts used when evaluating the linear model formula.
    #' @param ... any additional arguments passed to \link{model.matrix}.
    initialize = function(formula, data, theta, family = gaussian(), weights = NULL, cost = NULL, contrasts = NULL, alpha = 1, lambda = 0, upper = 1, ...) {

      stopifnot("'family' must be an object of class \"family\", (e.g. gaussian(link = \"identity\")" = inherits(family, "family"))
      self$formula <- as.formula(formula)
      self$theta <- theta
      self$family <- family
      self$alpha <- 1
      self$lambda <- 0
      self$upper <- 1
      private$dots <- list(...)
      if(is.data.frame(data) || is.matrix(data)) {
        self$data <- as.data.frame(data)
      } else if(!is.data.frame(data)) {
        stop("'data' must be a data.frame or a matrix")
      }
      if(!is.null(contrasts) && is.list(contrasts)) {
        for(nm in names(contrasts)) {
          contrasts(self$data[[nm]]) <- contrasts[[nm]]
        }
        private$contrasts <- contrasts
      }
      private$.Xinit <- do.call(stats::model.matrix.default, args = c(list(object = self$formula, data = self$data), private$dots))
      private$.X <- Xv_impl(private$.Xinit, theta = self$theta, family = self$family)
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
    #' Inputs should be named \code{formula}, \code{data}, \code{theta}, \code{family}, \code{weights}, \code{cost}, 
    #' \code{alpha}, \code{lambda} and/or \code{upper}.
    #' @param ... named arguments \code{formula}, \code{data}, \code{theta}, \code{family}, \code{weights}, \code{cost}, 
    #' \code{alpha}, \code{lambda} and/or \code{upper}.
    update = function(...) {

      dots <- list(...)

      if(is.element("family", names(dots)) && !is.null(dots$family)) {
		  stopifnot("'family' must be an object of class \"family\", (e.g. gaussian(link = \"identity\")" = inherits(family, "family"))
		  self$family <- dots$family
        ## cascade updates
        if(!is.element("data", names(dots))) {
          dots$data <- self$data
        }
      }
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
        if(!is.null(private$contrasts)) {
          for(nm in names(private$contrasts)) {
            contrasts(self$data[[nm]]) <- private$contrasts[[nm]]
          }
        }
        mf <- do.call(model.frame, args = c(list(formula = self$formula, data = self$data), private$dots))
        private$.Xinit <- stats::model.matrix.default(attr(mf, "terms"), mf)
        private$.X <- Xv_impl(private$.X, theta = self$theta, family = self$family)
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
    #' If specified, must be one of \code{"D"}, \code{"A"}, \code{"I"}, \code{"G"} or \code{"alias"}.}
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
    #' \deqn{\Sigma \ = \ (X^T U^{1/2} \Omega U^{1/2} X)^{-1}}
    #' where \eqn{\Omega} is \eqn{\mathrm{diag}(\omega_1, \ldots, \omega_n)} with \eqn{\omega_i = w_i v_i}, and \eqn{U^{1/2}} is
    #' \eqn{\mathrm{diag}(\sqrt{u(x_1, \theta)}, \ldots, \sqrt{u(x_n, \theta)})}, with GLM weights \eqn{u(x_i, \theta) = \frac{1}{V(\mu_i)[g'(\mu_i)]^2}} 
    #' evaluated at the nominal parameter values. 
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
    #' \deqn{\mathrm{SEP}(x_i) \ = \ \sigma \sqrt{x_i^T \Sigma x_i}, \quad \mathrm{for } i = 1,\ldots,n}
    #' where \eqn{\Sigma} is the covariance matrix as returned by \code{$vcov()}.
    #' If \code{type = "response"}, the standard errors of prediction are transformed to the response scale by an application of the Delta method:
    #' \deqn{\mathrm{SEP}_{y}(x_i) \ = \ \mathrm{SEP}(x_i) \Big| \frac{\partial \mu}{\partial \eta }\Big|}
    #' Use the \code{newdata} argument to evaluate the standard error of prediction at design points not in the candidate design grid defined by \code{$data}.
    #' The design weights are initialized from the \code{design_weights} field, which is the optimal (approximate) design calculated by \code{$optimize()}.
    #' To override the default approximate (optimal) design weights, a vector of design weights can be passed using the \code{design_weights} argument.
    #' The vector must be the same length as the number of design points in \code{data}.
    #' @param m (optional) integer number of design points at which to evaluate the estimated coviarance matrix. If missing, the
    #' design weights remain unscaled, which for the default approximate (optimal) design is equivalent to \code{m = 1}.
    #' @param newdata (optional) data frame or matrix with design points (rows) at which to return the standard error of prediction and which
    #' can be evaluated by \code{$formula}.
    #' @param design_weights (optional) vector of design weights for which to evaluate the estimated covariance matrix.
    #' @param type the type of prediction required. The default choice \code{"link"} returns the standard error of prediction on the scale of the linear predictors, 
    #' whereas \code{"response"} returns the standard error of prediction on the scale of the response variable.
    #' @param sigma2 numeric error variance used to scale the standard errors. Defaults to 1, returning the standardized standard error of prediction.
    sep = function(m, newdata, design_weights = NULL, type = c("link", "response"), sigma2 = 1) {
      type <- match.arg(type, c("link", "response"))
      if(!missing(newdata)) {
        if(is.data.frame(newdata) || is.matrix(newdata)) {
          newdata <- as.data.frame(newdata)
        } else if(!is.data.frame(newdata)) {
          stop("'newdata' must be a data.frame or a matrix")
        }
        if(!is.null(private$contrasts)) {
          for(nm in names(private$contrasts)) {
            contrasts(newdata[[nm]]) <- private$contrasts[[nm]]
          }
        }
        mf <- do.call(model.frame, args = c(list(formula = self$formula, data = newdata), private$dots))
        xpred <- stats::model.matrix.default(attr(mf, "terms"), mf)
      } else {
        xpred <- private$.Xinit
      }
      Sigma <- self$vcov(m, design_weights, sigma2)
      ses <- sqrt(rowSums(xpred * (xpred %*% Sigma)))
      if(type == "response") {
        eta <- c(xpred %*% self$theta)   
        ses <- abs(self$family$mu.eta(eta)) * ses
      }
      return(ses)
    },

    #' @description
    #' Calculates the correlation matrix of the columns of the design matrix based on an approximate or rounded (exact) design \eqn{(w_1, \ldots, w_n)} according to:
    #' \deqn{R \ = \ D^{-1/2}(X_c^T U^{1/2} \Omega U^{1/2} X_c)D^{-1/2}}
    #' where \eqn{\Omega} is \eqn{\mathrm{diag}(\omega_1, \ldots, \omega_n)} with \eqn{\omega_i = w_i v_i}, and \eqn{X_c} is the design matrix \eqn{X} centered by the (weighted)
    #' column means:
    #' \deqn{\frac{\sum_{i=1}^n \omega_i X_{ij}}{\sum_{i=1}^n \omega_i}, \quad j = 1,\ldots, p}
    #' \eqn{D^{-1/2}} is \eqn{\mathrm{diag}(1/\sqrt{D_1},\ldots, 1/\sqrt{D_p})}, with \eqn{D_1,\ldots,D_p} the diagonal elements of \eqn{X_c^T U^{1/2} \Omega U^{1/2} X_c}.
    #' If \code{center = FALSE}, the uncentered correlation matrix is returned, which corresponds to the scale-invariant (unit diagonal) version of the information matrix 
    #' \eqn{M(w | X, v, \theta, F, g)}. 
    #' The design weights are initialized from the \code{design_weights} field, which is the optimal (approximate) design calculated by \code{$optimize()}.
    #' To override the default approximate (optimal) design weights, a vector of design weights can be passed using the \code{design_weights} argument.
    #' The vector must be the same length as the number of design points in \code{data}.
    #' @param design_weights (optional) vector of design weights for which to evaluate the information matrix.
    #' @param center logical, whether to center the design columns when evaluating the correlation matrix, defaults to \code{TRUE}.
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
        private$.Xinit
      } else {
        stop("To modify the design matrix, update the data or model formula with `$update()`.", call. = FALSE)
      }
    }
  ),

  private = list(
    dots = list(),
    .X = NULL,
    .Xinit = NULL,
    contrasts = NULL
  )

)

Xv_impl <- function(X, theta, family) {
  if(is.list(theta)) {
	  theta <- unlist(theta, recursive = FALSE)
  }
  if(!is.null(names(theta))) {
	  idx <- match(colnames(X), names(theta))
	  if(any(is.na(idx))) {
		  stop("names of 'theta' do not match with the column names of '$X'")
	  }
	  theta <- theta[idx]
  }	else {
	  stopifnot("'theta' must be the same length as the number of columns in '$X'" = identical(length(theta), ncol(X)))
  }
  eta <- c(X %*% theta)    ## eta = X * theta
  mu <- family$linkinv(eta)   ## mu = g^{-1}(eta)
  dg_mu <- 1 / family$mu.eta(eta)    ## g'(mu) = d.eta/d.mu = 1 / (d.mu/d.eta)
  v <- 1 / (family$variance(mu) * dg_mu^2)   ##v = 1 / (V(mu) * g'(mu)^2)
  Xv <- sqrt(v) * X
  return(Xv)
}
