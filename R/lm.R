#' Penalized linear optimal design
#'
#' @description
#' An R6 class for managing penalized optimal experiment designs for linear models.
#' This class encapsulates the linear model formula, the discrete design grid,
#' optimization of an approximate continuous design, augmentation of an existing design,
#' and optimal rounding to exact integer designs.
#'
#' @section Optimal design problem:
#' The following penalized optimization problem is solved by \code{$optimize()}:
#' \deqn{w^* \ = \ \arg \max_{w} \Phi(M(w | X, v) + \lambda I) - \alpha \cdot c^T w}
#' where,
#' \itemize{
#' \item \eqn{w = \{w_1, \ldots, w_n \}} is the vector of design weights.
#' \item \eqn{c = \{c_1, \ldots, c_n \}} is the cost assigned to each design weight with scaling/penalty parameter \eqn{\alpha}.
#' \item \eqn{v = \{v_1, \ldots, v_n \}} is a vector of fixed user weights.
#' \item \eqn{M(w | X, v)} is the information matrix depending on the weights \eqn{w},
#' the design matrix \eqn{X} created by \link{model.matrix}, and in the case of weighted regression
#' the heterogeneous weights \eqn{v}.
#' \item \eqn{\lambda I} is a diagonal matrix to regularize the information matrix with ridge penalty parameter \eqn{\lambda}.
#' \item \eqn{\Phi} is the chosen optimality criterion, generally a function of the information matrix.
#' }
#' Given the linear model \code{formula}, the information matrix (conditional on the error variance) is:
#' \deqn{M(w | X) \ \propto \ \sum_{i = 1}^n w_i v_i x_i x_i^T}
#'
#' @section Augmented optimal design:
#' Given a vector of initial design weights \eqn{w_0} and associated contribution weight \eqn{\gamma \in [0, 1]},
#' the following augmented penalized optimization problem is solved by \code{$augment()}:
#' \deqn{w^* \ = \ \arg \max_{w} \Phi(M(w | \gamma, w_0, X, v) + \lambda I) - \alpha \cdot c^T w}
#' using the augmented information matrix,
#' \deqn{M(w | \gamma, w_0, X, v) \ = (1 - \gamma) M(w_0 | X, v) + \gamma M(w | X, v)}
#'
#' @section Stratified optimal design:
#' The \code{$strata()} method stratifies the design by requiring that individual strata \eqn{X_1, \ldots, X_S}
#' are represented proportionally with given weights \eqn{p_1, \ldots, p_S}, (\eqn{\sum_i p_i = 1}). If stratification is applied,
#' the following optimization problem is solved by \code{$optimize()}:
#' \deqn{
#'  \begin{array}{rl}
#'    w^* \ = \ & \arg \max_{w} \Phi(M(w | X, v) + \lambda I) - \alpha \cdot c^T w, \\
#'    & \mathrm{subject\ to} \ \sum_{i \in X_s} w_i = p_s, \quad \quad \mathrm{for} \ s = 1,\ldots, S
#'  \end{array}
#' }
#' using the same notation as in the \sQuote{Optimal design problem} section above. This means that the total weight
#' assigned to the design points in stratum \eqn{X_i} is equal to \eqn{p_i} for each \eqn{i = 1,\ldots,S}.
#' Similarly, if stratification is applied, \code{$augment()} solves the augmented optimization problem with stratification
#' constraints.
#'
#' @section Sparse optimal design:
#' The \code{$sparsify()} method solves the penalized optimization problem using a custom penalty function \eqn{P(w)}
#' to promote sparsity in the approximate design \eqn{w^*}. Because the sparse opimization problem is
#' not convex, a sequential convex approximation is used substituting the (concave) sparseness penalty by linearized (convex) costs
#' on the design weights \eqn{w} similar to the cost penalties above. The following options are supported:
#' \subsection{\code{"entropy"}-penalty}{
#' Shannon entropy penalty, \eqn{P(w) = -\sum_i^n w_i \log(w_i)}, with linearized costs (iteration \eqn{k}):
#' \deqn{c_i^{(k)} = -(1 + \log(w_i^{(k)} + \epsilon)) }
#' }
#' \subsection{\code{"log-sum"}-penalty}{
#' Log-sum penalty, \eqn{P(w) = \sum_i^n \log(w_i + \epsilon)}, with linearized costs (iteration \eqn{k}):
#' \deqn{c_i^{(k)} = \frac{1}{w_i^{(k)} + \epsilon}}
#' }
#' \subsection{\code{"power"}-penalty}{
#' Power (bridge) penalty, \eqn{P(w) = \sum_i^n (w_i + \epsilon)^p}, with linearized costs (iteration \eqn{k}):
#' \deqn{c_i^{(k)} = p (w_i^{(k)} + \epsilon)^{p - 1}, \quad 0 < p < 1}
#' }
#'
#' @section Optimality criteria:
#' The following optimality criteria \eqn{\Phi(\cdot)} are supported in \code{$optimize()}, \code{$augment()} and \code{$sparsify()}:
#' \subsection{\code{"D"}-optimality}{
#' Maximize the volume of the information matrix \eqn{M}, or equivalently minimize the volume of
#' the confidence ellipsoid of the parameter estimates:
#' \deqn{\Phi(M) = \det(M)^{1/p}}
#' where \eqn{p} is the number of columns in the design matrix \eqn{X}.
#' With non-zero ridge/cost penalties, the (penalized) optimality criterion is calculated as:
#' \deqn{\Phi_c(w, M) = \exp\left(\log(\Phi(M + \lambda I)) - \alpha \cdot c^T w \right)}
#' }
#' \subsection{ \code{"A"}-optimality}{
#' Maximize the inverse trace of the inverse information \eqn{M^{-1}}, or equivalently minimize the
#' average variance of the parameter estimates:
#' \deqn{\Phi(M) = \left(\mathrm{Tr}(M^{-1})\right)^{-1}}
#' With non-zero ridge/cost penalties, the (penalized) optimality criterion is calculated as:
#' \deqn{\Phi_c(w, M) = \Phi(M + \lambda I) - \alpha \cdot c^T w}
#' }
#' \subsection{\code{"I"}-optimality}{
#' Maximize the inverse of the average prediction variance across the design space, (equivalently
#' minimize the average prediction variance):
#' \deqn{\Phi(M) = \left(\sum_{i=1}^n x_i^T M^{-1} x_i\right)^{-1}}
#' With non-zero ridge/cost penalties, the (penalized) optimality criterion is calculated as:
#' \deqn{\Phi_c(w, M) = \Phi(M + \lambda I) - \alpha \cdot c^T w}
#' }
#' \subsection{\code{"G"}-optimality}{
#' Maximize the inverse of the maximum prediction variance across the design space, (equivalently
#' minimize the maximum prediction variance):
#' \deqn{\Phi(M) = \left(\max_{x_i} x_i^T M^{-1} x_i\right)^{-1}}
#' With non-zero ridge/cost penalties, the (penalized) optimality criterion is calculated as:
#' \deqn{\Phi_c(w, M) = \Phi(M + \lambda I) - \alpha \cdot c^T w}
#' }
#' \subsection{\code{"alias"}-optimality}{
#' Maximize one minus the average absolute cosine similarity (also uncentered correlation or normalized dot product)
#' between columns of the design matrix \eqn{X}:
#' \deqn{\Phi(w, v, X) = 1 - \frac{2}{p(p - 1)} \sum_{i < j} \left| \frac{\langle \omega x_i, x_j \rangle}{\Vert \sqrt{\omega} x_i \Vert \Vert \sqrt{\omega} x_j \Vert} \right|}
#' where \eqn{\omega x_i} is the elementwise-product of \eqn{\omega = (w_1 v_1, \ldots, w_n v_n)} and the \eqn{i}-th design column \eqn{x_i = (x_{1,i},\ldots,x_{n,i})}.
#' If all columns in the design matrix \eqn{X} are centered, this is equivalent to minimizing the average correlation between columns in \eqn{X}.
#' With non-zero cost penalties, the (penalized) optimality criterion is calculated as:
#' \deqn{\Phi_c(w, v, X) = \Phi(w, v, X) - \alpha \cdot c^T w}
#' }
#' @field formula a linear model defined as a one-sided \link{formula}.
#' @field data a data frame or matrix containing all design points (rows) at which \code{formula} can be evaluated.
#' @field X the design matrix obtained from all design points in \code{data} and the model \code{formula}.
#' @field orthogonal the dot product (orthogonality) constraint matrix. This matrix contains all dot product (upper limit)
#' constraints between design columns parsed from the \code{orthogonal} formula(s) in \code{$new(orthogonal = ...)} or \code{$update(orthogonal = ...)}.
#' For design columns with no dot product constraint, the matrix entry is \code{Inf}.
#' @field zeros the proportion of zeros vector. This vector contains the required proportion of zeros for each column in the design matrix
#' parsed from the \code{zeros} formula(s) in \code{$new(zeros = ...)} or \code{$update(zeros = ...)}.
#' For design columns with no constraint, the value is \code{NA}.
#' @field col_means the imposed column mean vector. This vector contains the imposed column mean constraints for each column in the design matrix
#' parsed from the \code{col_means} formula in \code{$new(col_means = ...)} or \code{$update(col_means = ...)}.
#' For design columns with no mean constraint, the value is \code{NA}.
#' @field cost an optional numeric vector containing costs associated to each design point (row) in \code{data}.
#' @field design_weights a numeric vector with optimal design weights populated after optimization with \code{$optimize()}, \code{$augment()}, or \code{$sparsify()}.
#'  The attribute \code{criterion} lists the used design criterion and the achieved optimal criterion value.
#' @field weights an optional numeric vector containing fixed user weights associated to each design point (row) in \code{data}.
#' @field strata a list with elements \code{strata} and \code{proportions} imposing a partition of the design points into individual (disjoint) strata
#' with given weight proportions.
#' @field alpha the cost scale/penalty parameter in the penalized optimization problem, only used if \code{$cost} is available.
#' @field lambda the ridge penalty added to the diagonal of the information matrix in the optimization problem. The ridge penalty is not used
#' for the \code{alias} optimality criterion.
#' @field upper the design weight upper limit. The default is 1, corresponding to no upper limit constraint on the design weights.
#' @param data a data frame or matrix containing all design points (rows) at which \code{formula} can be evaluated.
#' @param cost an optional numeric vector containing costs associated to each design point (row) in \code{data}.
#' Also accepts a one-sided \link{formula} using the column names present in \code{data}. If \code{NULL},
#' the standard (non-penalized) optimal design problem is solved.
#' @param weights an optional numeric vector containing fixed user weights associated to each design point (row) in \code{data}.
#' Also accepts a one-sided \link{formula} using the column names present in \code{data}. If \code{NULL},
#' the unweighted optimal design problem is solved.
#' @param orthogonal an optional two- or three-part formula or list of formulas imposing dot product (orthogonality)
#' constraints between the columns of the design matrix. In the three-part formula \code{x1 + x2 ~ y1 + y2 ~ 0},
#' dot products are constrained between all combinations of columns in the first part of the formula (\code{x1}, \code{x2})
#' and the second part of the formula (\code{y1}, \code{y2}). The third part (\code{0}) specifies the
#' upper limit on the dot product(s). If the formula only has two parts, e.g. \code{x1 + x2 ~ y1 + y2}, zero
#' dot product constraints are imposed. Also accepts a list of formulas, in which case constraints are added in the
#' order in which they appear in the list.
#' @param zeros an optional two-part formula or list of formulas imposing constraints on the proportion of zero values
#' in the columns of the design matrix. In the two-part formula \code{x1 + x2 ~ 1/3}, the proportion of zeros in the
#' columns specified by the first part of the formula (\code{x1}, \code{x2}) is constrained to be equal to the
#' (nonnegative) numeric value in the second part of the formula (1/3). Note that imposing constraints on the proportion
#' of zeros in the optimized design is only relevant if the design grid actually contain zero entries. For a list of formulas,
#' zero constraints are added in the order in which they appear in the list, analogous to the \code{orthogonal} argument.
#' @param col_means an optional two-part formula imposing constraints on the column means of the design matrix.
#' In the two-part formula \code{x1 + x2 ~ 0}, the means of the columns specified by the first part of the formula (\code{x1}, \code{x2})
#' are constrained to be equal to the numeric value in the second part of the formula (0).
#' @param criterion character string specifying the optimality criterion to use, see also \sQuote{Details}. The following choices are supported:
#' \tabular{l}{
#' \code{"D"} D-optimality (default). \cr
#' \code{"A"} A-optimality. \cr
#' \code{"I"} I-optimality. \cr
#' \code{"G"} G-optimality. \cr
#' \code{"alias"} alias-optimality.
#' }
#' @param scs_control an optional list of control parameters to tune the SCS optimization. See \code{\link{scs_control_dflt}}
#' for the available control parameters and their default values.
#' @param verbose display verbose messages printed by the scs solver. Defaults to \code{FALSE}.
#' @param return_value return the optimal design weights. If set to \code{FALSE}, the design weights can be obtained from the \code{design_weights} field.
#' @param alpha numeric value, cost scale parameter in the penalized optimization problem, only used if \code{$cost} is available. The default is 1.
#' @param lambda numeric value, ridge penalty added to the diagonal of the information matrix in the penalized optimization problem. The default is 0.
#' @param upper numeric value in (0, 1], design weight upper limit. The default is 1, corresponding to no constraint on the design weights.
#' @import Matrix
#' @importFrom stats lm model.frame model.matrix as.formula dist runif na.omit terms quantile contrasts
#' @importFrom R6 R6Class
#' @importFrom utils installed.packages combn modifyList globalVariables
#' @importFrom progress progress_bar
#' @importFrom methods as
#' @importFrom mirai daemons_set mirai race_mirai everywhere
#' @references Pukelsheim, F., and Rieder S. (1992) \emph{"Efficient rounding of approximate designs"}, Biometrika.
#' @seealso \url{https://www.cvxgrp.org/scs/}
#' @useDynLib cvxoptdes, .registration = TRUE
#' @examples
#'
#' ## Getting started
#'
#' ## 3 factors w/ 5, 5 and 3 levels
#' data <- expand.grid(
#'     x1 = -2:2,
#'     x2 = -2:2,
#'     x3 = -1:1
#' )
#' design <- lm_design$new(
#'     formula = ~(x1 + x2 + x3)^2 + I(x1^2) + I(x2^2) + I(x3^2),   ## full quadratic model
#'     data = data    ## design points
#' )
#' ## approximate D-optimal design
#' design$optimize(criterion = "D")
#' ## integer D-optimal design w/ 13 design points
#' (exact_design <- design$round(m = 13, method = "optimal", seed = 1))
#' design$subset(replicates = exact_design)
#' ## integer design covariance and correlation matrices
#' design$vcov(design_weights = exact_design, sigma2 = 1)
#' design$corr(design_weights = exact_design)
#'
#' if(FALSE) {
#'  ## parallel evaluation
#'  mirai::daemons(n = parallel::detectCores() / 2, seed = 1)
#'  design$round(m = 100, method = "optimal", n_repeats = 1000)
#'  mirai::daemons(0)
#' }
#'
#' ## Penalized optimal design
#'
#' ## assign cost c = 1 if x3 = 1, and c = 0 otherwise
#' design$update(
#'   cost = ~pmax(x3, 0),
#'   alpha = 1
#' )
#' ## penalized approximate D-optimal design
#' design$optimize(criterion = "D")
#'
#' ## assign v = 1 if x3 = 1, and v = 0.1 otherwise
#' design$update(
#'  weights = ~pmax(x3, 0.1),
#'  cost = NULL
#' )
#' ## penalized approximate D-optimal design
#' design$optimize(criterion = "D")
#' ## exact design w/ 15 design points
#' design$round(m = 15, method = "optimal", seed = 1)
#'
#' ## drop weights term
#' design$update(weights = NULL)
#'
#' ## Stratified optimal design
#'
#' ## add stratification constraints
#' design$stratify(
#'   strata = ~x3,
#'   proportions = c(1/4, 1/2, 1/4)
#' )
#' ## approximate D-optimal design w/ stratification
#' design$optimize(criterion = "D")
#' ## verify stratification
#' tapply(design$design_weights, design$data$x3, sum)
#' ## exact design w/ stratification
#' design$round(m = 12, method = "optimal", seed = 1)
#' ## drop stratification
#' design$stratify(strata = NULL)
#'
#' ## Orthogonality constraints
#'
#' ## orthogonal main effects wrt all other columns
#' design$update(
#'   orthogonal = x1 + x2 + x3 ~ . ~ 0,
#'   alpha = 0.075  ## sparseness penalty
#' )
#' design$orthogonal
#' ## approximate D-optimal design w/ orthogonality constraints
#' design$optimize(criterion = "D")
#' ## sparsify support by iterative entropy reweighting
#' design$sparsify(max_iter = 5, criterion = "D", reweight_fun = "entropy")
#'
#' if(FALSE) {
#'  ## integer design w/ orthogonality constraints
#'  (ortho_design <- design$round(m = 8, method = "orthogonal", seed = 1))
#'  ## verify correlations
#'  design$corr(ortho_design)
#' }
#'
#' ## Optimal design augmentation
#'
#' ## Example from [Goos & Jones (2011), Chapter 3]
#' ## 6 factors w/ 2 levels each
#' data <- expand.grid(
#'   x1 = c(-1, 1),
#'   x2 = c(-1, 1),
#'   x3 = c(-1, 1),
#'   x4 = c(-1, 1),
#'   x5 = c(-1, 1),
#'   x6 = c(-1, 1)
#' )
#'
#' ## linear model w/ main effects and selected interactions
#' design <- lm_design$new(
#'   formula = ~1 + x1 + x2 + x3 + x4 + x5 + x6 + x2:x3 + x1:x6 + x2:x6 + x3:x6 + x4:x6 + x1:x4,
#'   data = data
#' )
#' ## initial design of size n=12
#' init_design <- c(0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
#'  0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
#'  0, 1, 0, 0, 0, 0, 1)
#' ## approximate optimal augmented design
#' design$augment(design_weights = init_design, gamma = 0.4, criterion = "D")
#' ## exact augmented design
#' (augment_design <- design$round(m = 8, method = "optimal", seed = 1, augment_design = init_design))
#' ## proposed design points
#' design$subset(replicates = augment_design - init_design)
#'
#' ## augmented design w/ orthogonality constraints
#' design_ortho <- design$clone()
#' design_ortho$update(
#'   orthogonal = list(
#'     x1 + x2 + x3 + x4 + x5 + x6 ~ x1 + x2 + x3 + x4 + x5 + x6 ~ 0  ## uncorrelated main effects
#'   )
#' )
#' ## approximate augmented design
#' design_ortho$augment(
#'   design_weights = init_design, gamma = 0.4, criterion = "D", return_value = FALSE
#' )
#'
#' if(FALSE) {
#'   ## integer design w/ orthogonality constraints
#'   ortho_augment_design <- design_ortho$round(
#'     m = 8, method = "orthogonal", augment_design = init_design
#'   )
#'   ## verify correlations
#'   design_ortho$corr(ortho_augment_design)
#' }
#'
#' @export
lm_design <- R6Class(
  "lm_design",
  public = list(
    formula = NULL,
    data = NULL,
    weights = NULL,
    cost = NULL,
    design_weights = NULL,
    strata = NULL,
    alpha = 1,
    lambda = 0,
    upper = 1,

    #' @description
    #' Constructor. Initializes the experimental design to optimize with a linear model formula and discrete design grid.
    #' @param formula a linear model defined as a one-sided \link{formula}.
    #' @param contrasts optional contrasts used when evaluating the linear model formula.
    #' @param ... any additional arguments passed to \link{model.matrix}.
    initialize = function(formula, data, orthogonal = NULL, weights = NULL, cost = NULL, zeros = NULL, col_means = NULL,
                          contrasts = NULL, alpha = 1, lambda = 0, upper = 1, ...) {

      self$formula <- as.formula(formula)
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
      private$.X <- do.call(stats::model.matrix.default, args = c(list(object = self$formula, data = self$data), private$dots))
      if(!is.null(weights)) {
        self$weights <- weights_impl(weights, self$data, name = "weights")
      }
      if(!is.null(orthogonal)) {
        private$.orthogonal <- parse_ortho_constraints(orthogonal, self$data, self$formula, private$dots)
      }
      if(!is.null(zeros)) {
        private$.zeros <- parse_zero_constraints(zeros, private$.X, self$data, self$formula, private$dots)
        private$.Xzero <- array(0L, dim = dim(private$.X))
        private$.Xzero[abs(private$.X) < 1e-6] <- 1L
      } else if(all(apply(private$.X, 2L, function(x) any(abs(private$.X) < 1e-6)))) {
        private$.Xzero <- array(0L, dim = dim(private$.X))
        private$.Xzero[abs(private$.X) < 1e-6] <- 1L
      }
      if(!is.null(col_means)) {
        stopifnot("'col_means' must be a two-part formula" = inherits(col_means, "formula") && length(col_means) == 3L)
        private$.colmeans <- parse_mean_constraints(col_means, private$.X, self$data, self$formula, private$dots)
      }
      if(!is.null(cost)) {
        self$cost <- weights_impl(cost, self$data, name = "cost")
      }
      if(!identical(alpha, 0)) {
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
    #' Inputs should be named \code{formula}, \code{data}, \code{weights}, \code{orthogonal}, \code{cost}, \code{zeros},
    #' \code{col_means}, \code{alpha}, \code{lambda} and/or \code{upper}.
    #' @param ... named arguments \code{formula}, \code{data}, \code{weigts}, \code{orthogonal}, \code{cost}, \code{zeros},
    #' \code{col_means}, \code{alpha}, \code{lambda} and/or \code{upper}.
    update = function(...) {

      dots <- list(...)

      if(is.element("formula", names(dots)) && !is.null(dots$formula)) {
        self$formula <- as.formula(dots$formula)
        ## cascade updates
        if(!is.element("data", names(dots))) {
          dots$data <- self$data
        }
        if(!is.element("orthogonal", names(dots)) && !is.null(private$.orthogonal)) {
          warning("removing existing orthogonality constraints, call `$update(orthogonal = ...)` to set new constraints.")
          private$.orthogonal <- NULL
        }
        if(!is.element("zeros", names(dots)) && !is.null(private$.zeros)) {
          warning("removing existing constraints on zero entries, call `$update(zeros = ...)` to set new constraints.")
          private$.zeros <- NULL
          private$.Xzero <- NULL
        }
        if(!is.element("col_means", names(dots)) && !is.null(private$.colmeans) ) {
          private$.colmeans <- NULL
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
        private$.Xzero <- NULL
        private$.colmeans <- NULL
        mf <- do.call(model.frame, args = c(list(formula = self$formula, data = self$data), private$dots))
        private$.X <- stats::model.matrix.default(attr(mf, "terms"), mf)
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
      if(is.element("orthogonal", names(dots))) {
        if(is.null(dots$orthogonal)) {
          private$.orthogonal <- NULL
        } else {
          private$.orthogonal <- parse_ortho_constraints(dots$orthogonal, self$data, self$formula, private$dots)
        }
      }
      if(is.element("zeros", names(dots))) {
        if(is.null(dots$zeros)) {
          private$.zeros <- NULL
          private$.Xzero <- NULL
        } else {
          private$.zeros <- parse_zero_constraints(dots$zeros, private$.X, self$data, self$formula, private$dots)
          private$.Xzero <- array(0L, dim = dim(private$.X))
          private$.Xzero[abs(private$.X) < 1e-6] <- 1L
        }
      }
      if(is.null(private$.Xzero) && all(apply(private$.X, 2L, function(x) any(abs(private$.X) < 1e-6)))) {
        private$.Xzero <- array(0L, dim = dim(private$.X))
        private$.Xzero[abs(private$.X) < 1e-6] <- 1L
      }
      if(is.element("col_means", names(dots))) {
        if(is.null(dots$col_means)) {
          private$.colmeans <- NULL
        } else {
          stopifnot(
            "'col_means' must be a two-part formula" = inherits(dots$col_means, "formula") && length(dots$col_means) == 3L
          )
          private$.colmeans <- parse_mean_constraints(dots$col_means, private$.X, self$data, self$formula, private$dots)
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
    #' the partition of the design points into individual (disjoint) strata. Also accepts a one-sided \link{formula} using
    #' column names present in \code{data}, in which case the strata are defined by the \link{interaction} between terms in the
    #' \link{model.frame} obtained from the \link{formula}.
    #' @param proportions an optional numeric vector containing the weight proportions of the individual strata. Must be of the same length
    #' as the number of factor levels in \code{strata}. If \code{NULL}, each stratum is assigned the same weight proportion.
    stratify = function(strata, proportions = NULL) {

      if(is.null(strata)) {
        self$strata <- NULL
      } else {
        if(!is.null(private$.zeros)) {
          warning("stratification constraints override existing constraints on zero entries")
          private$.zeros <- NULL
          private$.Xzero <- NULL
        }
        if(inherits(strata, "formula")) {
          mf <- do.call(model.frame, args = c(list(formula = strata, data = self$data), private$dots))
          strata <- do.call(interaction, args = c(as.list(mf), drop = TRUE))
        }
        stopifnot(
          is.factor(strata), nlevels(strata) > 0, length(strata) == nrow(self$data)
        )
        if(is.null(proportions)) {
          proportions <- rep(1 / nlevels(strata), nlevels(strata))
        } else {
          stopifnot(
            is.numeric(proportions), length(proportions) == nlevels(strata),
            all(proportions >= 0), !any(is.infinite(proportions)), any(proportions > 0)
          )
          proportions <- proportions / sum(proportions)
        }
        if(!is.null(names(proportions))) {
          nms_order <- match(names(proportions), levels(strata))
          if(any(is.na(nms_order))) {
            warning("'proportions' is a named vector, but its names do not match with the levels of 'strata'")
          } else {
            proportions <- proportions[nms_order]
          }
        }
        stopifnot(
          is.numeric(proportions), length(proportions) == nlevels(strata)
        )
        self$strata <- list(
          strata = strata,
          proportions = proportions
        )
      }
      return(invisible(self))
    },

    #' @description
    #' Optimizes the augmented penalized approximate design given a vector of initial design weights
    #' @param design_weights a numeric vector with initial design weights, must be the same length as the number of design points in \code{data}.
    #' Also accepts a one-sided \link{formula} made up of column names present in \code{data}.
    #' @param gamma numeric contribution weight of the initial design in the augmented design problem, see also \sQuote{Augmented optimal design}
    #' above. Defaults to 0.5.
    augment = function(design_weights, gamma = 0.5, criterion = c("D", "A", "I", "G", "alias"), scs_control = scs_control_dflt(), verbose = FALSE, return_value = TRUE) {

      if(inherits(design_weights, "formula")) {
        design_weights <- tryCatch(eval(design_weights[[2L]], envir = as.list(self$data)))
        if(inherits(design_weights, "error"))
          stop(sprintf("failed to evaluate 'design_weights' formula: %s", design_weights$message))
        if(!is.numeric(design_weights) || !identical(length(design_weights), nrow(private$.X)))
          stop("'design_weights' formula failed to return a numeric vector of the same length as number of rows in '$data'")
        if(any(is.na(design_weights) | design_weights < 0))
          stop("missing or negative values returned by 'design_weights' formula")
      } else if(!is.numeric(design_weights) || any(is.na(design_weights) | design_weights < 0) || !identical(length(design_weights), nrow(private$.X))) {
        stop("'design_weights' must be a nonnegative vector of the same length as number of rows in '$data'")
      }
      stopifnot(
        is.numeric(gamma), length(gamma) == 1, gamma > 0, gamma < 1,
        any(design_weights > 0)
      )
      gamma <- as.numeric(gamma)
      w0 <- design_weights / sum(design_weights)

      ## optimize augmented design
      design <- optimize_impl(
        X = private$.X,
        criterion = criterion,
        cost = self$cost,
        weights = self$weights,
        orthogonal = private$.orthogonal,
        strata = self$strata,
        zeros = private$.zeros,
        Xzero = private$.Xzero,
        colmeans = private$.colmeans,
        control = scs_control,
        alpha = self$alpha,
        lambda = self$lambda,
        upper = self$upper,
        verbose = verbose,
        gamma = gamma,
        w0 = w0
      )

      if(check_status(design)) {
        self$design_weights <- structure(
          design$w,
          criterion = list(
            crit.value = design$crit_val,
            crit = criterion
          )
        )
        if(isTRUE(return_value)) {
          return(self$design_weights)
        }
      }

    },

    #' @description
    #' Optimizes the penalized approximate design based on the selected criterion.
    optimize = function(criterion = c("D", "A", "I", "G", "alias"), scs_control = scs_control_dflt(), verbose = FALSE, return_value = TRUE) {

      stopifnot(is.logical(verbose), is.logical(return_value))

      design <- optimize_impl(
        X = private$.X,
        criterion = criterion,
        cost = self$cost,
        weights = self$weights,
        orthogonal = private$.orthogonal,
        strata = self$strata,
        zeros = private$.zeros,
        Xzero = private$.Xzero,
        colmeans = private$.colmeans,
        control = scs_control,
        alpha = self$alpha,
        lambda = self$lambda,
        upper = self$upper,
        verbose = verbose
      )

      if(check_status(design)) {
        self$design_weights <- structure(
          design$w,
          criterion = list(
            crit.value = design$crit_val,
            crit = criterion
          )
        )
        if(isTRUE(return_value)) {
          return(self$design_weights)
        }
      }

    },

    #' @description
    #' Sparsify the penalized approximate design by iterative reweighting of the optimization problem.
    #'
    #' @details
    #' The sparseness cost penalties are scaled by the same \code{alpha} tuning parameter as for a user-specific cost vector, see \code{\link{scs_control_dflt}}.
    #' If a cost vector \code{c_1,\ldots,c_n} associated to the design weights is already specified. The sparseness penalty costs are multiplied
    #' by the existing cost vector.
    #'
    #' @param max_iter positive integer, maximum number of reweighted optimizations/iterations to run.
    #' @param reweight_fun character reweighting function, one of \code{"entropy", "log-sum", "power"}, see \sQuote{Details} for the definitions.
    #' @param eps positive numeric value to ensure finite weights in reweighting functions, see \sQuote{Details}.
    #' @param p postive numeric power with \eqn{0 < p < 1}. Only used if \code{reweight_fun = "power"}.
    #' @param tol numeric tolerance used as a threshold for mean absolute relative difference in design weights between iterations to stop early.
    #' @param recenter logical, recenter design columns based on the current design weights after each iteration. Defaults to \code{FALSE}.
    #' @param design_weights (optional) vector of design weights used as initial approximate design.
    #' @param show_progress logical, display a progress bar during the iterative reweighting procedure. Defaults to \code{TRUE}.
    sparsify = function(max_iter, criterion = c("D", "A", "I", "G", "alias"), reweight_fun = c("entropy", "log-sum", "power"), eps = .Machine$double.eps,
                        tol = .Machine$double.eps^(1/4), p = 0.5, design_weights = NULL, recenter = FALSE, scs_control = scs_control_dflt(),
                        return_value = TRUE, show_progress = TRUE) {

      if(is.null(design_weights)) {
        design_weights <- rep(1 / nrow(private$.X), nrow(private$.X))
      } else {
        design_weights <- design_weights / sum(design_weights)
      }
      reweight_fun <- match.arg(reweight_fun, c("entropy", "power", "log-sum"))
      stopifnot(
        is.numeric(design_weights), length(design_weights) == nrow(private$.X),
        is.numeric(max_iter), length(max_iter) == 1, max_iter >= 1,
        is.numeric(eps), length(eps) == 1, eps > 0,
        is.numeric(tol), length(tol) == 1, tol > 0,
        reweight_fun != "power" || (length(p == 1) && p > 0 && p < 1)
      )
      if(reweight_fun == "entropy") {
        cost_fun <- function(w) -(1 + (log(pmax(w, 0) + eps)))
        penalty_fun <- function(w) -sum((pmax(w, 0) + eps) * log(pmax(w, 0) + eps))
      } else if(reweight_fun == "log-sum") {
        cost_fun <- function(w) 1 / (pmax(w, 0) + sqrt(eps))
        penalty_fun <- function(w) sum(log(pmax(w, 0) + sqrt(eps)))
      } else if(reweight_fun == "power") {
        cost_fun <- function(w) p * (pmax(w, 0) + eps)^(p - 1)
        penalty_fun <- function(w) sum((pmax(w, 0) + eps)^p)
      }
      if(is.null(attr(design_weights, "criterion"))) {
        design_weights <- update_crit(w0 = NULL, w1 = design_weights, X = private$.X, m = length(design_weights), crit = criterion, lambda = self$lambda)
      }
      penalty0 <- penalty_fun(design_weights)
      ## iterative reweighting scheme
      if(isTRUE(show_progress)) {
        pb <- progress_bar$new(
          format = "Iter: :iter (P(w)/P(w0): :rel_cost) [:bar] :percent (:elapsed)",
          total = max_iter + 1,
          clear = FALSE,
          show_after = 0
        )
        pb$tick(
          tokens = list(
            iter = 0,
            rel_cost = "1"
          )
        )
      }
      warm_start <- NULL   ## warm-start solution
      if(!is.null(self$weights)) {
        Xw <- sqrt(self$weights) * private$.X
      } else {
        Xw <- private$.X
      }
      for(iter in seq_len(max_iter)) {
        if(isTRUE(recenter)) {
          xmean <- colSums(design_weights / sum(design_weights) * private$.X)
          X <- sweep(private$.X, 2, xmean, "-")
        } else {
          X <- private$.X
        }
        ## recalculate cost penalties
        cost <- cost_fun(design_weights)
        ## avoid exploding costs
        if(reweight_fun != "entropy") {
          cost <- pmin(cost, 5.0e4 / self$alpha)
        }
        if(!is.null(self$cost)) {
          cost <- cost * self$cost
        }
        design <- optimize_impl(
          X = X,
          criterion = criterion,
          cost = cost,
          weights = self$weights,
          orthogonal = private$.orthogonal,
          strata = self$strata,
          zeros = private$.zeros,
          Xzero = private$.Xzero,
          colmeans = private$.colmeans,
          control = scs_control,
          alpha = self$alpha,
          lambda = self$lambda,
          upper = self$upper,
          start = warm_start,
          verbose = FALSE,
          return_sol = TRUE
        )
        if(isTRUE(show_progress)) {
          pb$tick(
            tokens = list(
              iter = iter,
              rel_cost = sprintf("%.3g", penalty_fun(design$w) / penalty0)
            )
          )
        }
        if(suppressWarnings(check_status(design))) {
          mae <- mean(abs((design_weights - design$w) / design_weights))
          if(mae < tol) {
            break
          }
          design_weights <- design$w
          ## update warm-start solution
          warm_start <- design$scs_solution
        } else {
          break
        }
      }
      design_weights <- update_crit(w1 = design_weights, X = Xw, crit = criterion, lambda = self$lambda)
      attr(design_weights, "criterion")[["crit.rel.eff"]] <- NULL
      self$design_weights <- design_weights
      if(isTRUE(return_value)) {
        return(self$design_weights)
      }

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
    #' @details
    #' \subsection{Parallel evaluation}{
    #' The repeated rounding/exchange algorithms can be run in parallel using the \code{mirai} asynchronous evaluation framework.
    #' The parallel background R processes are configured by the user outside of \code{$round()} with \code{mirai::daemons}.
    #' If daemon connections are available, \code{$round()} automatically runs the repeated exchange algorithms in parallel
    #' using the configured daemons. If the underlying BLAS/LAPACK libraries are already multi-threaded, this may
    #' inhibit the use of additional threads by the \code{mirai} asynchronous evaluation. In this case, it is recommended to restrict
    #' BLAS/LAPACK to a single thread with the \code{OPENBLAS_NUM_THREADS}, \code{MKL_NUM_THREADS} or \code{OMP_NUM_THREADS}
    #' environment variables, depending on the BLAS library that is used. In addition, note that there is a small overhead in populating the
    #' environments of the background R processes, which may result in increased runtime compared to sequential evaluation in the main R process,
    #' especially when \code{n_repeats} is small or the rounding/exchange algorithm does not take much time to complete.
    #' }
    #' \subsection{HiGHS solver}{
    #' The \code{"orthogonal"} rounding procedure uses HiGHS to find an exact design subject to orthogonality constraints.
    #' This requires either a local installation of the HiGHS solver or the \code{highs} package to be installed.
    #' To point to a local HiGHS installation, the environment variable \code{HIGHS_HOME} must be set to the \code{highs}
    #' root directory when installing \code{cvxoptdes}. If installation fails, verify that \code{highs_c_api.h} exists, typically
    #' present in \code{$HIGHS_HOME/include/highs/interfaces}, which is required to compile the source code.
    #' }
    #'
    #' @param m integer number of design points included in the exact design.
    #' @param method the method used to calculate the exact design. The following choices are supported:
    #' \tabular{l}{
    #'  \code{"optimal"} optimality criterion point-exchange algorithm (default). \cr
    #'  \code{"efficient"} efficient rounding method of Pukelsheim & Reider (1992). \cr
    #'  \code{"dist-sum"} maximize total pairwise (Euclidean) distance between design points. \cr
    #'  \code{"dist-min"} maximize minimum pairwise (Euclidean) distance between design points. \cr
    #'  \code{"orthogonal"} minimizes aliasing subject to orthogonality constraints, see \sQuote{Details}.
    #' }
    #' @param max_iter integer number of maximum iterations in point-exchange algorithm. Only used when
    #' \code{method = "optimal"}, \code{method = "dist-sum"} or \code{method = "dist-min"}.
    #' @param n_repeats integer number of (random) repeats of the rounding/exchange algorithm. The exact design
    #' with the highest criterion value is returned. Not used for \code{method = "orthogonal"}. Defaults to 10.
    #' @param replicates logical, allow replicate design points in the exact design, defaults to \code{TRUE}.
    #' Can also be an integer number, limiting the maximum number of allowed replicates per design point.
    #' This argument overrides the \code{$upper} value constraint if existing.
    #' @param tol numeric tolerance used as lower threshold to continue iterations in point-exchange algorithm.
    #' Also used as threshold to ensure positive diagonal entries of \eqn{M(w | X)} when \code{method = "orthogonal"}.
    #' @param seed (optional) random seed used to generate jitter to break ties in rounding the initial exact design.
    #' Not used in the case of (\code{mirai}) asynchronous evaluation, which requires the random seed to be set via the
    #' \code{seed} argument in \code{mirai::daemons}.
    #' @param design_weights (optional) vector of design weights used as initial approximate design, with length equal
    #' to the number of available design points. If \code{TRUE} (default), the approximate design is initialized from the
    #' \code{design_weights} field. If \code{FALSE}, the initial approximate design is generated by random sampling.
    #' @param augment_design (optional) integer design to augment, with length equal to the number of available design points.
    #' @param show_progress logical, display a progress bar monitoring the repeats of the rounding/exchange algorithm.
    #' Defaults to \code{TRUE}.
    #' @param ... additional arguments to tune algorithm/solver settings:
    #' \itemize{
    #' \item \itemize{
    #' \item \code{criterion}, character string specifying the optimality criterion for the rounding procedure in case \code{design_weights}
    #' has no \code{"criterion"} attribute, or to override the criterion specified in the \code{"criterion"} attribute.
    #' If specified, must be one of \code{"D"}, \code{"A"}, \code{"I"}, \code{"G"} or \code{"alias"}.}
    #' \item \itemize{
    #' \item \code{highs_control}, a named list with control arguments supported by the HiGHS solver, only used if \code{method = "orthogonal"}.
    #' Defaults to \code{highs_control = list(time_limit = 180, log_to_console = TRUE)}.}
    #' \item \itemize{
    #' \item \code{highs_write_mps}, \code{.mps} file path. If provided, the HiGHS MILP is exported in MPS (Mathematical Programming System)
    #' format to the provided file path. Can also be \code{TRUE}, in which case a new \code{.mps} file is created in the current working directory.
    #' Only used if \code{method = "orthogonal"}.}
    #' }
    #' @seealso https://ergo-code.github.io/HiGHS/dev/installation/
    #' @seealso https://mirai.r-lib.org/
    round = function(m, method = c("optimal", "efficient", "dist-sum", "dist-min", "orthogonal"), max_iter = 100, n_repeats = 10, replicates = TRUE,
                     tol = .Machine$double.eps^(1/3), seed = NULL, design_weights = TRUE, augment_design =  NULL, show_progress = TRUE, ...) {
      ## check arguments
      method <- match.arg(method, c("optimal", "efficient", "dist-sum", "dist-min", "orthogonal"))
      dots <- list(...)
      if(isTRUE(design_weights)) {
        design_weights <- self$design_weights
        if(is.null(design_weights) && method != "orthogonal") {
          stop("no initial design available to round. Call `$optimize()` first to initialize an approximate optimal design.")
        }
      } else if(isFALSE(design_weights)) {
        design_weights <- NULL
      } else {
        design_weights <- design_weights / sum(design_weights)
      }
      if(!is.null(augment_design)) {
        stopifnot(
          is.numeric(augment_design), length(augment_design) == nrow(private$.X)
        )
        augment_design <- round(augment_design)
      }
      stopifnot(
        is.null(design_weights) || is.numeric(design_weights), is.null(design_weights) || length(design_weights) == nrow(private$.X),
        is.null(design_weights) || !any(is.na(design_weights)),
        is.numeric(m), length(m) == 1, m > 0,
        is.numeric(tol), length(tol) == 1, tol > 0,
        is.numeric(max_iter), length(max_iter) == 1, max_iter > 0
      )
      if(method == "orthogonal") {
        stopifnot("method 'orthogonal' can only be used if orthogonality constraints are present." = !is.null(private$.orthogonal))
        if(!is.null(dots$criterion)) {
          crit <- match.arg(dots$criterion, c("I", "A", "D", "G", "alias"))
        } else if(!is.null(attr(design_weights, "criterion"))) {
          crit <- match.arg(attr(design_weights, "criterion")$crit, c("I", "A", "D", "G", "alias"))
        } else {
          crit <- "alias"
        }
        if(!is.null(dots$highs_write_mps)) {
          stopifnot(isTRUE(dots$highs_write_mps) || (is.character(dots$highs_write_mps) && length(dots$highs_write_mps) == 1))
          if(isTRUE(dots$highs_write_mps)) {
            highs_write_mps <- tempfile(pattern = "highs_export", tmpdir = getwd(), fileext = ".mps")
          } else {
            highs_write_mps <- normalizePath(dots$highs_write_mps, mustWork = FALSE)
            stopifnot("'highs_write_mps' must be an .mps file." = grepl("\\.mps$", basename(highs_write_mps)))
          }
        } else {
          highs_write_mps <- NULL
        }
      } else if(method == "optimal" || method == "efficient") {
        if(method == "efficient" && is.null(design_weights)) {
          stop("method 'efficient' can only be used if approximate design weights are available.")
        }
        if(!is.null(dots$criterion)) {
          crit <- match.arg(dots$criterion, c("I", "A", "D", "G", "alias"))
        } else if(!is.null(attr(design_weights, "criterion"))) {
          crit0 <- attr(design_weights, "criterion")
          crit <- match.arg(crit0$crit, c("I", "A", "D", "G", "alias"))
          attr(crit, "crit0") <- crit0
        } else {
          stop("no 'criterion' attribute found in 'design_weights'. Specify the 'criterion' argument to set an optimality criterion for the rounded design.")
        }
        if(!is.null(private$.orthogonal) || !is.null(private$.zeros)) {
          warning("orthogonality and/or zero constraints are not taken into account when using method 'optimal' or 'efficient'.")
        }
        if(!is.null(private$.colmeans) && any(!is.na(private$.colmeans))) {
          warning("column mean constraints are not taken into account when using method 'optimal' or 'efficient'.")
        }
        if(!is.null(dots$highs_write_mps)) {
          warning("'highs_write_mps' is skipped when using method 'optimal' or 'efficient'.")
        }
      } else if(method == "dist-min" || method == "dist-sum") {
        crit <- method
        if(!is.null(private$.orthogonal) || !is.null(private$.zeros)) {
          warning("orthogonality and/or zero constraints are not taken into account when using method 'dist-min' or 'dist-sum'.")
        }
        if(!is.null(private$.colmeans) && any(!is.na(private$.colmeans))) {
          warning("column mean constraints are not taken into account when using method 'dist-min' or 'dist-sum'.")
        }
        if(!is.null(dots$highs_write_mps)) {
          warning("'highs_write_mps' is skipped when using method 'dist-min' or 'dist-sum'.")
        }
      }
      m <- as.integer(m)
      max_iter <- as.integer(max_iter)
      n_repeats <- as.integer(n_repeats)
      if(!is.null(self$weights)) {
        X <- sqrt(self$weights) * private$.X
      } else {
        X <- private$.X
      }
      ## determine upper limit
      if(!isTRUE(replicates)) {
        if(isFALSE(replicates)) {
          replicates <- 1
        } else {
          stopifnot(is.numeric(replicates), length(replicates) == 1)
          replicates <- max(min(replicates, nrow(X)), 1)
        }
        upper <- replicates / m
        if(upper < 1 && nrow(X) * replicates < m + sum(augment_design)) {
          stop("maximum number of replicates cannot be satisfied, increase the allowed 'replicates' number or reduce the number of design points 'm'.")
        }
      } else {
        upper <- self$upper
        if(upper < 1 && nrow(X) * max(floor(m * upper), 1) < m + sum(augment_design)) {
          stop("upper limit constraint in `$upper` cannot be satisfied, reduce the number of design points 'm'.")
        }
      }
      ## round strata proportions
      if(!is.null(self$strata)) {
        if(method == "orthogonal") {
          warning("stratification constraints are not optimized when using method 'orthogonal'.")
          strata <- NULL
        } else {
          if(identical(method, "efficient")) {
            warning("stratification constraints are not optimized when using method 'efficient'.")
          }
          strata <- self$strata
          if(m < sum(strata$proportions > 0)) {
            warning("'m' is smaller than the number of individual strata, increase 'm' to (minimally) cover all strata.")
          }
          strata$proportions <- round_puk(strata$proportions, m)
        }
      } else {
        strata <- NULL
      }
      w_opt <- NULL
      crit_opt <- -Inf
      strata_opt <- Inf
      if(isTRUE(show_progress) && n_repeats > 1) {
        pb <- progress_bar$new(
          format = sprintf("Iter: :n (%s-crit: :crit) [:bar] :percent (:elapsed)", crit),
          total = n_repeats,
          clear = FALSE,
          show_after = 0
        )
      }
      if(method != "orthogonal") {
        if(!is.null(design_weights)) {
          supp_tol <- ifelse(crit == "alias" && method != "efficient", 0, (.Machine$double.eps)^(1/4))
          supp <- design_weights > supp_tol
          if(!is.null(augment_design)) {
            supp <- supp & (augment_design < .Machine$double.eps)
          }
          if(sum(supp) < m && method == "optimal") {
            supp_tol <- -Inf
            supp <- rep(TRUE, nrow(X))
          }
        } else {
          supp_tol <- NA
          supp <- rep(TRUE, nrow(X))
        }
        n_supp <- sum(supp)
        if(method == "efficient" && m < n_supp) {
          stop(sprintf("method 'efficient' requires m >= %d the size of the approximate design support.", n_supp))
        }
        if(!mirai::daemons_set() || n_repeats < 2) {
          if(!is.null(seed)) {
            if(exists(".Random.seed", envir = globalenv(), mode = "integer", inherits = FALSE)) {
              seed0 <- get(".Random.seed", envir = globalenv(), mode = "integer", inherits = FALSE)
              on.exit(set.seed(seed0))
            } else {
              on.exit({
                set.seed(seed = NULL)
                rm(".Random.seed", envir = globalenv())
              })
            }
            set.seed(seed)
          }
          for(n in 1:n_repeats) {
            if(isTRUE(show_progress) && n_repeats > 1) {
              pb$tick(
                tokens = list(
                  n = n,
                  crit = sprintf("%.3g", crit_opt)
                )
              )
            }
            w <- round_impl(
              design_weights = design_weights,
              X = X,
              m = m,
              n_supp = n_supp,
              supp_tol = supp_tol,
              supp = supp,
              tol = tol,
              max_iter = max_iter,
              crit = crit,
              method = method,
              augment_design = augment_design,
              upper = upper,
              strata = strata,
              cost = self$cost,
              lambda = self$lambda,
              alpha = self$alpha,
              orthogonal = private$.orthogonal
            )
            crit_val <- attr(w, "criterion")$crit.value
            if(!is.null(strata)) {
              strata_diff <- sum(abs(attr(w, "strata") - strata$proportions))
              if(strata_diff > strata_opt) {
                crit_val <- -Inf
              } else {
                strata_opt <- strata_diff
              }
            }
            if(crit_val > crit_opt) {
              w_opt <- w
              crit_opt <- crit_val
            }
          }
        } else {
          if(!is.null(seed)) {
            warning("'seed' argument is skipped for asynchronous evaluation, set the random seed with `mirai::daemons(seed = ...)`.")
          }
          mirai::everywhere(
            .expr = {
              requireNamespace("cvxoptdes", quietly = TRUE)
            },
            round_impl = round_impl,
            design_weights = design_weights,
            X = X,
            m = m,
            n_supp = n_supp,
            supp_tol = supp_tol,
            supp = supp,
            tol = tol,
            max_iter = max_iter,
            crit = crit,
            method = method,
            augment_design = augment_design,
            upper = upper,
            strata = strata,
            cost = as.vector(self$cost),
            lambda = self$lambda,
            alpha = self$alpha,
            orthogonal = private$.orthogonal
          )
          m_list <- lapply(1:n_repeats, function(n) {
            mirai::mirai(
              .expr = {
                round_impl(
                  design_weights, X, m, n_supp, supp_tol, supp, tol, max_iter, crit, method,
                  augment_design, upper, strata, cost, lambda, alpha, orthogonal
                )
              }
            )
          })
          while (length(m_list) > 0) {
            m_id <- mirai::race_mirai(m_list)
            w <- m_list[[m_id]]$data
            m_list[[m_id]] <- NULL
            crit_val <- attr(w, "criterion")$crit.value
            if(!is.null(strata)) {
              strata_diff <- sum(abs(attr(w, "strata") - strata$proportions))
              if(strata_diff > strata_opt) {
                crit_val <- -Inf
              } else {
                strata_opt <- strata_diff
              }
            }
            if(crit_val > crit_opt) {
              w_opt <- w
              crit_opt <- crit_val
            }
            if(isTRUE(show_progress) && n_repeats > 1) {
              pb$tick(
                tokens = list(
                  n = n_repeats - length(m_list),
                  crit = sprintf("%.3g", crit_opt)
                )
              )
            }
          }
        }
      } else if(method == "orthogonal") {
        ## highs optimization
        w_opt <- round_highs(
          w = design_weights,
          m = m,
          X = private$.X,
          tol = tol,
          orthogonal = private$.orthogonal,
          zeros = private$.zeros,
          Xzero = private$.Xzero,
          colmeans = private$.colmeans,
          augment = augment_design,
          control = dots$highs_control,
          seed = seed,
          write_mps = highs_write_mps
        )
        if(!is.null(w_opt)) {
          if(is.null(augment_design)) {
            w_opt <- update_crit(design_weights, w_opt, X, m = m, crit = crit, lambda = self$lambda)
          } else {
            w_opt <- update_crit(design_weights, w_opt + augment_design, X, m = m + sum(augment_design), crit = crit, lambda = self$lambda)
          }
        }
      }
      return(w_opt)
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
      ## check arguments
      if(is.null(design_weights)) {
        w <- self$design_weights
        if(is.null(w)) {
          stop("No design weights available to evaluate. Call `$optimize()` first to initialize an approximate optimal design.")
        }
      } else {
        w <- design_weights / sum(design_weights)
      }
      stopifnot(
        is.numeric(w), length(w) == nrow(private$.X), any(w > 0),
        is.logical(lambda), length(lambda) == 1L,
        is.logical(alpha), length(alpha) == 1L
      )
      criterion <- match.arg(criterion, c("A", "I", "D", "G", "alias"))
      if(!isFALSE(lambda)) {
        lambda <- self$lambda
      } else {
        lambda <- 0
      }
      if(!isFALSE(alpha)) {
        alpha <- self$alpha
        cost <- self$cost
      } else {
        alpha <- 0
        cost <- NULL
      }
      if(!is.null(self$weights)) {
        X <- sqrt(self$weights) * private$.X
      } else {
        X <- private$.X
      }
      crit <- switch(criterion,
                     D = d_crit(w, X, cost = cost, lambda = lambda, alpha = alpha),
                     A = a_crit(w, X, cost = cost, lambda = lambda, alpha = alpha),
                     I = i_crit(w, X, cost = cost, lambda = lambda, alpha = alpha),
                     G = g_crit(w, X, cost = cost, lambda = lambda, alpha = alpha),
                     alias = alias_crit(w, X, cost = cost, alpha = alpha)
      )
      return(crit)
    },

    #' @description
    #' Calculates the variance-covariance matrix of the model coefficients based on an approximate or rounded (exact) design \eqn{(w_1, \ldots, w_n)} according to:
    #' \deqn{\Sigma \ = \ (X^T \Omega X)^{-1}}
    #' where \eqn{\Omega} is \eqn{\mathrm{diag}(\omega_1, \ldots, \omega_n)} with \eqn{\omega_i = w_i v_i}.
    #' The design weights are initialized from the \code{design_weights} field, which is the optimal (approximate) design calculated by \code{$optimize()}.
    #' Use the \code{m} argument to evaluate the covariance matrix for a specific number of design points, in which case the sum of the design weights is scaled to \code{m}.
    #' To override the default approximate (optimal) design weights, a vector of design weights can be passed using the \code{design_weights} argument.
    #' The vector must be the same length as the number of design points in \code{data}. If the vector of design weights is already an integer design, leave
    #' \code{m} unspecified.
    #' @param m (optional) integer number of design points at which to evaluate the estimated coviarance matrix. If missing, the
    #' design weights remain unscaled, which for the default approximate (optimal) design is equivalent to \code{m = 1}.
    #' @param design_weights (optional) vector of design weights for which to evaluate the estimated covariance matrix.
    #' @param sigma2 numeric error variance used to scale the covariance matrix. Defaults to 1, returning the unscaled covariance matrix.
    vcov = function(m, design_weights = NULL, sigma2 = 1) {
      ## check arguments
      stopifnot(
        is.numeric(sigma2), length(sigma2) == 1, sigma2 > 0
      )
      if(is.null(design_weights)) {
        w <- self$design_weights
        if(is.null(w)) {
          stop("No design weights available to evaluate. Call `$optimize()` first to initialize an approximate optimal design.")
        }
        w <- pmax(w, 0)
        if(missing(m)) {
          m <- 1
        } else {
          m <- as.numeric(m)
        }
      } else {
        w0 <- pmax(design_weights, 0)
        w <- w0 / sum(w0)
        if(missing(m)) {
          m <- sum(w0)
        } else {
          m <- as.numeric(m)
        }
      }
      stopifnot(
        is.numeric(w), length(w) == nrow(private$.X), any(w > 0),
        is.numeric(m), length(m) == 1, m > 0
      )
      if(!is.null(self$weights)) {
        X <- sqrt(self$weights) * private$.X
      } else {
        X <- private$.X
      }
      M <- crossprod(sqrt(w * m) * X)
      Minv <- tryCatch(chol2inv(chol(M)), error = function(e) e)
      if(inherits(Minv, "error")) {
        Minv <- array(Inf, dim = dim(M))
        warning("Information matrix is not positive definite, unable to calculate its inverse.")
      }
      Sigma <- sigma2 * Minv
      rownames(Sigma) <- colnames(Sigma) <- colnames(private$.X)
      return(Sigma)
    },

    #' @description
    #' Calculates the standard error of prediction or the square root prediction variance across the design grid based on an approximate or rounded (exact)
    #' design \eqn{(w_1, \ldots, w_n)} according to:
    #' \deqn{\mathrm{SEP}(x_i) \ = \ \sigma \sqrt{x_i^T \Sigma x_i}, \quad \mathrm{for } i = 1,\ldots,n}
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
        if(!is.null(private$contrasts)) {
          for(nm in names(private$contrasts)) {
            contrasts(newdata[[nm]]) <- private$contrasts[[nm]]
          }
        }
        mf <- do.call(model.frame, args = c(list(formula = self$formula, data = newdata), private$dots))
        xpred <- stats::model.matrix.default(attr(mf, "terms"), mf)
      } else {
        xpred <- private$.X
      }
      Sigma <- self$vcov(m, design_weights, sigma2)
      ses <- sqrt(rowSums(xpred * (xpred %*% Sigma)))
      return(ses)
    },

    #' @description
    #' Calculates the correlation matrix of the columns of the design matrix based on an approximate or rounded (exact) design \eqn{(w_1, \ldots, w_n)} according to:
    #' \deqn{R \ = \ D^{-1/2}(X_c^T \Omega X_c)D^{-1/2}}
    #' where \eqn{\Omega} is \eqn{\mathrm{diag}(\omega_1, \ldots, \omega_n)} with \eqn{\omega_i = w_i v_i}, and \eqn{X_c} is the design matrix \eqn{X} centered by the (weighted)
    #' column means:
    #' \deqn{\frac{\sum_{i=1}^n \omega_i X_{ij}}{\sum_{i=1}^n \omega_i}, \quad j = 1,\ldots, p}
    #' \eqn{D^{-1/2}} is \eqn{\mathrm{diag}(1/\sqrt{D_1},\ldots, 1/\sqrt{D_p})}, with \eqn{D_1,\ldots,D_p} the diagonal elements of \eqn{X_c^T \Omega X_c}.
    #' If \code{center = FALSE}, the uncentered correlation matrix is returned, which corresponds to the scale-invariant (unit diagonal) version of the information matrix \eqn{M(w | X, v)}.
    #' The design weights are initialized from the \code{design_weights} field, which is the optimal (approximate) design calculated by \code{$optimize()}.
    #' To override the default approximate (optimal) design weights, a vector of design weights can be passed using the \code{design_weights} argument.
    #' The vector must be the same length as the number of design points in \code{data}.
    #' @param design_weights (optional) vector of design weights for which to evaluate the information matrix.
    #' @param center logical, whether to center the design columns when evaluating the correlation matrix, defaults to \code{TRUE}.
    corr = function(design_weights = NULL, center = TRUE) {
      ## check arguments
      if(is.null(design_weights)) {
        w <- self$design_weights
        if(is.null(w)) {
          stop("No design weights available to evaluate. Call `$optimize()` first to initialize an approximate optimal design.")
        }
        w <- pmax(w, 0)
      } else {
        w0 <- pmax(design_weights, 0)
        w <- w0 / sum(w0)
      }
      stopifnot(
        is.numeric(w), length(w) == nrow(private$.X), any(w > 0), is.logical(center)
      )
      if(!is.null(self$weights)) {
        omega <- w * self$weights
      } else {
        omega <- w
      }
      xcols <- apply(private$.X, 2, function(x) length(unique(x)) > 1)
      X <- private$.X[, xcols, drop = FALSE]
      if(isTRUE(center)) {
        wmean <- colSums(omega * X) / sum(omega)
        Xc <- sweep(X, 2, wmean, "-")
        M <- crossprod(sqrt(omega) * Xc)
      } else {
        M <- crossprod(sqrt(omega) * X)
      }
      D <- sqrt(pmax(diag(M), .Machine$double.eps))
      R <- M / outer(D, D)
      rownames(R) <- colnames(R) <- colnames(X)
      return(R)
    },

    #' @description
    #' Subsets the design points associated to an exact design vector with integer replicates as returned by \code{$round()}.
    #' The vector must be the same length as the number of design points in \code{data}.
    #' @param replicates integer vector of design point replicates of the same length as the candidate grid in \code{data}
    subset = function(replicates) {
        n <- nrow(self$data)
        stopifnot(
          is.numeric(replicates),
          length(replicates) == n,
          "'replicates' must be an integer design vector" = all(as.integer(replicates) == replicates),
          all(replicates >= 0),
          !any(is.na(replicates))
        )
        return(self$data[rep(1:n, times = as.integer(replicates)), ])
    }

  ),

  active = list(
    X = function(value) {
      if (missing(value)) {
        private$.X
      } else {
        stop("To modify the design matrix, update the data or model formula with `$update()`.", call. = FALSE)
      }
    },
    orthogonal = function(value) {
      if(missing(value)) {
        private$.orthogonal
      } else {
        stop("To modify the orthogonality constraint matrix, call `$update(orthogonal = ...)`.", call. = FALSE)
      }
    },
    zeros = function(value) {
      if(missing(value)) {
        private$.zeros
      } else {
        stop("To modify the weight constraints on zero entries, call `$update(zeros = ...)`.", call. = FALSE)
      }
    },
    col_means = function(value) {
      if(missing(value)) {
        private$.colmeans
      } else {
        stop("To modify the column mean constraints, call `$update(col_means = ...)`.", call. = FALSE)
      }
    }
  ),

  private = list(
    dots = list(),
    .X = NULL,
    .Xzero = NULL,
    .orthogonal = NULL,
    .zeros = NULL,
    .colmeans = NULL,
    contrasts = NULL
  )

)

optimize_impl <- function(X, criterion, cost, weights = NULL, orthogonal = NULL, strata, zeros = NULL, Xzero = NULL, colmeans = NULL, control,
                          alpha = 1, lambda = 0, upper = 1, verbose = FALSE, return_sol = FALSE, gamma = 0, w0 = NULL, start = NULL, v_sca = NULL) {

  ## prepare configuration
  ctrl <- parse_ctrl(control, alpha, lambda, gamma, upper, criterion, verbose)
  ctrl_dbl <- ctrl$ctrl_dbl
  ctrl_int <- c(ctrl$ctrl_int, isTRUE(return_sol))

  if(!is.null(weights)) {
    X <- sqrt(weights) * X
  }
  ## K matrix, (only used for A-, I-optimality)
  K <- diag(ncol(X))
  if(identical(ctrl_int[1], 1L)) {
    ## I-optimality
    svdX <- svd(X)
    K <- 1 / sqrt(nrow(X)) * ((svdX$v %*% diag(svdX$d)) %*% t(svdX$v))
  }
  if(!is.null(strata)) {
    strata <- list(
      as.integer(strata$strata) - 1L,
      strata$proportions
    )
  }
  ## FIXME
  if(is.factor(zeros)) {
    zeros <- NULL
    Xzero <- NULL
  }

  ## optimize design
  design <- .Call(
    "R_scs_solve",
    X, cost, orthogonal, strata, zeros, Xzero, colmeans, environment(), w0, K, v_sca, start, ctrl_int, ctrl_dbl,
    PACKAGE = "cvxoptdes"
  )
  ## shift augmented design weights
  if(!is.null(w0) && gamma > 0) {
    design$w <- (design$w - (1 - gamma) * w0) / gamma
  }

  return(design)

}

weights_impl <- function(weights, data, name = "weights") {
  weights_formula <- NULL
  if(inherits(weights, "formula")) {
    weights_formula <- weights
    weights <- tryCatch(eval(weights[[2L]], envir = as.list(data)))
    if(inherits(weights, "error"))
      stop(sprintf("failed to evaluate '%s' formula: %s", name, weights$message))
    if(!is.numeric(weights) || !identical(length(weights), nrow(data)))
      stop(sprintf("'%s' formula failed to return a numeric vector of the same length as number of rows in '$data'", name))
    if(any(is.na(weights)))
      stop(sprintf("missing values returned by '%s' formula", name))
  } else if(!is.numeric(weights) || any(is.na(weights)) || !identical(length(weights), nrow(data))) {
    stop(sprintf("'%s' must be a numeric vector of the same length as number of rows in '$data'", name))
  }
  stopifnot("all elements of 'weights' must be positive" = name != "weights" || all(weights > 0))
  structure(as.numeric(weights), formula = weights_formula)
}

parse_ortho_constraints <- function(orthogonal, data, formula, dots) {
  if(is.matrix(orthogonal)) {
    stopifnot(
      is.numeric(orthogonal), !any(is.na(orthogonal))
    )
    if(is.integer(orthogonal)) {
      orthogonal[] <- as.numeric(orthogonal)
    }
    if(any(orthogonal < 0)) {
      stop("'orthogonal' matrix must contain only nonnegative values")
    }
    diag(orthogonal) <- Inf
    colnms <- colnames(do.call(stats::model.matrix.default, args = c(list(object = formula, data = data), dots)))
    dimnames(orthogonal) <- list(colnms, colnms)
    return(orthogonal)
  }
  if(inherits(orthogonal, "formula")) {
    orthogonal <- list(orthogonal)
  }
  if(is.list(orthogonal)) {
    stopifnot(
      length(orthogonal) > 0, all(sapply(orthogonal, inherits, "formula")), all(lengths(orthogonal) == 3L)
    )
    colnms <- colnames(do.call(stats::model.matrix.default, args = c(list(object = formula, data = data), dots)))
    ortho <- matrix(Inf, nrow = length(colnms), ncol = length(colnms), dimnames = list(colnms, colnms))
    for(frm in orthogonal) {
      sub_frm <- frm[[2]]
      if(!is.name(sub_frm) && identical(sub_frm[[1]], quote(`~`))) {
        ## three-part formula
        lhs_vars <- sub_frm[[2]]
        rhs_vars <- sub_frm[[3]]
        if(is.language(frm[[3]])) {
          upr <- as.numeric(eval(frm[[3]], envir = parent.frame()))
        } else {
          upr <- as.numeric(frm[[3]])  ## fails if not numeric
        }
      } else {
        ## two-part formula
        lhs_vars <- frm[[2]]
        rhs_vars <- frm[[3]]
        upr <- 0
      }
      lhs_frm <- update(formula, call("~", call("+", lhs_vars, 0)))
      lhs <- colnames(do.call(stats::model.matrix.default, args = c(list(object = lhs_frm, data = data), dots)))
      lhs <- intersect(lhs, colnms)
      rhs_frm <- update(formula, call("~", rhs_vars))
      rhs <- colnames(do.call(stats::model.matrix.default, args = c(list(object = rhs_frm, data = data), dots)))
      rhs <- intersect(rhs, colnms)
      ortho[lhs, rhs] <- upr
      ortho[rhs, lhs] <- upr
    }
  }
  diag(ortho) <- Inf
  return(ortho)
}

parse_zero_constraints <- function(zeros, X, data, formula, dots) {
  if(inherits(zeros, "formula")) {
    zeros <- list(zeros)
  }
  if(is.list(zeros)) {
    stopifnot(
      length(zeros) > 0, all(sapply(zeros, inherits, "formula")), all(lengths(zeros) == 2L) || all(lengths(zeros) == 3L)
    )
    colnms <- colnames(do.call(stats::model.matrix.default, args = c(list(object = formula, data = data), dots)))

    if(all(lengths(zeros) == 3L)) {
      zerovec <- rep(NA_real_, times = length(colnms))
    } else {
      zerovec <- factor(rep(NA_integer_, times = length(colnms)), levels = paste0("n", seq_along(zeros)))
      n <- 1L
    }
    names(zerovec) <- colnms
    for(frm in zeros) {
      lhs_vars <- frm[[2]]
      lhs_frm <- update(formula, call("~", call("+", lhs_vars, 0)))
      lhs <- colnames(do.call(stats::model.matrix.default, args = c(list(object = lhs_frm, data = data), dots)))
      lhs <- intersect(lhs, colnms)
      if(length(frm) == 3L) {
        ## two-part formula
        if(is.language(frm[[3]])) {
          zerovec[lhs] <- as.numeric(eval(frm[[3]], envir = parent.frame()))
        } else {
          zerovec[lhs] <- as.numeric(frm[[3]])  ## fails if not numeric
        }
      } else {
        ## one-part formula
        zerovec[lhs] <- paste0("n", n)
        n <- n + 1
      }
    }
    nozero <- apply(X, 2, function(x) !any(abs(x) < 1e-6))
    if(any(!is.na(zerovec[nozero]))) {
      zerovec[nozero] <- NA
      warning("Skipping zero entry weight constraints for design columns without zero entries", call. = FALSE)
    }
  }
  return(zerovec)
}

parse_mean_constraints <- function(frm, X, data, formula, dots) {
  colnms <- colnames(do.call(stats::model.matrix.default, args = c(list(object = formula, data = data), dots)))
  meanvec <- rep(NA_real_, times = length(colnms))
  names(meanvec) <- colnms
  lhs_vars <- frm[[2]]
  lhs_frm <- update(formula, call("~", call("+", lhs_vars, 0)))
  lhs <- colnames(do.call(stats::model.matrix.default, args = c(list(object = lhs_frm, data = data), dots)))
  lhs <- intersect(lhs, colnms)
  if(is.language(frm[[3]])) {
    meanvec[lhs] <- as.numeric(eval(frm[[3]], envir = parent.frame()))
  } else {
    meanvec[lhs] <- as.numeric(frm[[3]])  ## fails if not numeric
  }
  return(meanvec)
}

check_status <- function(design) {
  scs_status <- design$info_int[2]
  scs_msg <- switch(as.character(scs_status),
                    `2` = "SCS solver did not reach full convergence",
                    `0` = "SCS solver did not finish",
                    `-1` = "SCS optimization problem is unbounded (dual infeasible)",
                    `-2` = "SCS optimization problem has no feasible solution (dual unbounded)",
                    `-3` = "SCS optimization problem is indeterminate (numerical errors)",
                    `-4` = "SCS solver general failure",
                    `-5` = "SCS solver interrupted by SIGINT",
                    `-6` = "SCS solver did not fully converge, but returned an (approximate) unbounded certificate",
                    `-7` = "SCS solver did not fully converge, but returned an (approximate) infeasible certificate",
                    "SCS solver general failure"
  )
  if(scs_status < 1)
    stop(scs_msg)
  if(scs_status > 1)
    warning(scs_msg)
  return(TRUE)
}

round_impl <- function(design_weights, X, m, n_supp, supp_tol, supp, tol, max_iter, crit, method, augment_design = NULL,
                       upper = 1, strata = NULL, cost = NULL, lambda = 0, alpha = 1, orthogonal = NULL) {

  if(is.null(design_weights) || m < n_supp || crit == "alias") {
    ## missing initial design/incomplete support
    w <- round_random(design_weights, nrow(X), m, tol = supp_tol, prob_min = ifelse(crit == "alias", 1 / nrow(X), 0),
                      augment = augment_design, upper = upper)
  } else {
    ## complete support
    w <- round_puk(design_weights, m, augment = augment_design, upper = upper)
  }
  ## exchange optimization
  if(method == "optimal") {
    ## initialize G- by D-exchange
    w <- exchange_impl(w, X, ifelse(crit == "G", "D", crit), supp = supp, strata = strata, max_iter = max_iter,
                       tol = tol, augment = augment_design, cost = cost, lambda = lambda, alpha = alpha,
                       upper = upper, ortho = orthogonal, add_crit = crit != "G")
    if(crit == "G") {
      if(!is.null(augment_design)) {
        w <- w - augment_design
      }
      if(!is.null(attr(crit, "crit0"))) {
        w <- structure(w, criterion = attr(crit, "crit0"))
        attr(crit, "crit0") <- NULL
      }
      w <- exchange_impl(w, X, crit, strata = strata, max_iter = max_iter, tol = tol, augment = augment_design,
                         cost = cost, lambda = lambda, alpha = alpha, upper = upper, ortho = orthogonal)
    }
  } else if(method == "dist-min" || method == "dist-sum") {
    ## no weights/costs
    w <- dist_exchange(w, X, strata = strata, max_iter = max_iter, tol = tol, type = method,
                       augment = augment_design, upper = upper)
  } else if(method == "efficient") {
    if(!is.null(augment_design)) {
      w <- w + augment_design
      w <- update_crit(design_weights, w, X, m = sum(w), crit = crit, cost = cost, lambda = lambda, alpha = alpha, strata = strata)
      attr(w, "criterion")$crit.rel.eff <- NA
    } else {
      w <- update_crit(design_weights, w, X, m = sum(w), crit = crit, cost = cost, lambda = lambda, alpha = alpha, strata = strata)
    }
  }

  return(w)
}
