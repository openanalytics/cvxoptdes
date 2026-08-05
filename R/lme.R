#' Penalized linear mixed model optimal design
#'
#' @description
#' An R6 class for managing penalized optimal experiment designs for random intercept linear mixed models.
#' This class encapsulates a random intercept linear mixed model, the discrete design grid,
#' optimization of an approximate continuous design, augmentation of an existing design,
#' and optimal rounding to exact integer designs including split-plot designs.
#'
#' @section Optimal design problem:
#' The following penalized optimization problem is solved by \code{$optimize()}:
#' \deqn{w^* \ = \ \arg \max_{w} \Phi(M(w | X, Z, \eta, v) + \lambda I) - \alpha \cdot c^T w }
#' where,
#' \itemize{
#' \item \eqn{w = \{w_1, \ldots, w_n \}} is the vector of design weights.
#' \item \eqn{c = \{c_1, \ldots, c_n \}} is the cost assigned to each design weight with scaling/penalty parameter \eqn{\alpha}.
#' \item \eqn{v = \{v_1, \ldots, v_n \}} is a vector of fixed user weights.
#' \item \eqn{M(w | X, Z, \eta, v)} is the information matrix depending on the weights \eqn{w},
#' the fixed-effects design matrix \eqn{X}, the random-effects design matrix \eqn{Z}, the variance component ratio
#' \eqn{\eta = \sigma_u^2 / \sigma_e^2}, and in the case of weighted regression the heterogeneous weights \eqn{v}.
#' \item \eqn{\Phi} is the chosen optimality criterion, generally a function of the information matrix.
#' }
#' Given a random intercept linear mixed model (\eqn{X\beta + Zu + \epsilon}) defined by \code{formula}, the information matrix
#' of the fixed-effects model parameters, conditional on the random effects covariance matrix and the error variance, is given by:
#' \deqn{M(w | X, Z, \eta, v) \ = \ X^T\Omega^{1/2}V^{-1}\Omega^{1/2}X }
#' where \eqn{\Omega^{1/2}} is the diagonal matrix with \eqn{(\sqrt{w_1v_1}, \ldots, \sqrt{w_nv_n})} on the diagonal, and
#' \eqn{V^{-1}} is the inverse of the block-diagonal covariance matrix \eqn{V}. Partitioned into random effects levels \eqn{1, \ldots, S},
#' where group \eqn{s} has \eqn{n_s} observations, \eqn{V} consists of the level-specific blocks:
#' \deqn{V_s = \sigma_e^2 \left( I_{n_s} + \eta \boldsymbol{1}_{n_s} \boldsymbol{1}_{n_s}^T \right) }
#'
#' @section Augmented optimal design:
#' Given a vector of initial design weights \eqn{w_0} and associated contribution weight \eqn{\gamma},
#' the following augmented penalized optimization problem is solved by \code{$augment()}:
#' \deqn{w^* \ = \ \arg \max_{w} \Phi(M(w | \gamma, w_0, X, Z, \eta, v) + \lambda I) - \alpha \cdot c^T w}
#' using the augmented information matrix,
#' \deqn{M(w | \gamma, w_0, X, Z, \eta, v) \ = (1 - \gamma) M(w_0 | X, Z, \eta, v) + \gamma M(w | X, Z, \eta, v)}
#'
#' @section Stratified optimal design:
#' The \code{$strata()} method stratifies the design by requiring that individual strata \eqn{X_1, \ldots, X_S}
#' are represented proportionally with given weights \eqn{p_1, \ldots, p_S}, (\eqn{\sum_i p_i = 1}). If stratification is applied,
#' the following optimization problem is solved by \code{$optimize()}:
#' \deqn{
#'  \begin{array}{rl}
#'    w^* \ = \ & \arg \max_{w} \Phi(M(w | X, Z, \eta, v) + \lambda I) - \alpha \cdot c^T w, \\
#'    & \mathrm{subject\ to} \ \sum_{i \in X_s} w_i = p_s, \quad \quad \mathrm{for} \ s = 1,\ldots, S
#'  \end{array}
#' }
#' This means that the total weight assigned to the design points in stratum \eqn{X_i} is equal to \eqn{p_i} for each \eqn{i = 1,\ldots,S}.
#' Similarly, if stratification is applied, \code{$augment()} solves the augmented optimization problem \emph{with} stratification
#' constraints.
#'
#' @section Exact split-plot design:
#' The \code{$round_split_plot()} method rounds an approximate design to an exact (integer) split-plot design.
#' The rounding procedure enforces fixed levels of hard-to-change variables within each invidual \emph{whole plot},
#' where the whole plots are the blocks in the random effects factors (i.e. whole plot variables).
#' Mathematically, for an exact design \eqn{m^* = \{m_1, \ldots, m_n \}}, the split-plot rounding procedure solves:
#' \deqn{
#'  \begin{array}{rl}
#'    m^* \ = \ & \arg \max_{m} \Phi(M(m | X, Z, \eta, v) + \lambda I) - \alpha \cdot c^T (m / \sum m_i), \\
#'    & \mathrm{subject\ to} \ H(m) = Z(m) H_b \mathrm{\ \ and\ \ } \boldsymbol{1}_n^T Z(m) = M \boldsymbol{1}_b^T
#'  \end{array}
#' }
#' \eqn{H(m)} and \eqn{Z(m)} are the design matrices of the hard-to-change and random effects factors subject to \eqn{m}.
#' \eqn{H_b} is a matrix with \eqn{b} rows selected from the levels of the hard-to-change variables, where \eqn{b} is the
#' number of whole plots. The second constraint ensures that the number of design points at each whole plot level is
#' exactly \eqn{M}, the whole plot size. The split-plot designs are only available in the form of \emph{exact} designs, and not approximate designs
#' (as returned by e.g. \code{$optimize()}).
#'
#' @section Optimality criteria:
#' The optimality criteria \eqn{\Phi(\cdot)} supported by \code{$optimize()} and \code{$augment()} are: \emph{D}-,
#' \emph{A}-, \emph{I}-, \emph{G}-, and \emph{alias}-optimality, see \link{lm_design} for the definitions.
#'
#' @field formula a two-part \link{formula} of the form \code{~ fix.eff ~ ran.eff} representing a random intercept linear mixed model,
#' where \code{fix.eff} is a linear model formula encoding the fixed effects design matrix and \code{ran.eff} is an expression
#' encoding the levels of the random effects factor.
#' @field data a data frame or matrix containing all design points (rows) at which \code{formula} can be evaluated.
#' @field eta numeric variance component ratio \eqn{\sigma_u^2 / \sigma_e^2}, where \eqn{\sigma_u^2} is the
#' (single) variance component of the random effects and \eqn{\sigma_e^2} is the error variance.
#' @field X the fixed effects design matrix obtained from all design points in \code{data} and the model \code{formula}.
#' @field Z the random effects design matrix obtained from all design points in \code{data} and the model \code{formula}.
#' @field cost an optional numeric vector containing costs associated to each design point (row) in \code{data}.
#' @field design_weights a numeric vector with optimal design weights populated after optimization with \code{$optimize()} or \code{$augment()}.
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
#' Also accepts a one-sided \link{formula} made up of column names present in \code{data}. If \code{NULL},
#' the standard (non-penalized) optimal design problem is solved.
#' @param weights an optional numeric vector containing fixed user weights associated to each design point (row) in \code{data}.
#' Also accepts a one-sided \link{formula} using the column names present in \code{data}. If \code{NULL},
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
#' @importClassesFrom Matrix dgCMatrix dgRMatrix dsCMatrix
#' @importFrom Matrix t diag chol
#' @importFrom reformulas mkReTrms
#' @examples
#'
#' ## Getting started
#'
#' ## unbalanced design grid
#' data <- data.frame(
#'   Treatment = c(-1, -1, 0, 0, 1, 1, -1, 0, 0, 1, -1, 1),
#'   Block = factor(rep(c("B1", "B2", "B3", "B4", "B5", "B6"), each = 2))
#' )
#' ## fixed effects quadratic model in Treatment ~ random intercepts by Block
#' design <- lme_design$new(
#'   formula = ~ Treatment + I(Treatment^2) ~ Block,
#'   eta = 10,       ## variance component ratio
#'   data = data
#' )
#'
#' ## approximate D-optimal design
#' design$optimize(criterion = "D", max_iter = 10)
#' ## integer design w/ 18 design points
#' design$round(m = 6, method = "optimal", seed = 1)
#'
#' ## Split-plot design
#'
#' ## Case study from [Otava & Mylona (2025), Section 4.1]
#'
#' data <- expand.grid(
#'   A = -1:1,
#'   B = -1:1,
#'   C = -1:1,
#'   D = -1:1,
#'   week = factor(1:7, levels = 1:7)
#' )
#'
#' design <- lme_design$new(
#'   formula = ~ (A + B + C + D)^2 + I(A^2) + I(B^2) + I(C^2) + I(D^2) ~ week,
#'   eta = 1,
#'   data = data
#' )
#' ## initial D-optimal approximate design
#' design$optimize(criterion = "D", max_iter = 5, return_value = FALSE)
#' if(FALSE) {
#'   ## parallel evaluation
#'   mirai::daemons(n = parallel::detectCores() / 2, seed = 1)
#'   ## D-optimal split-plot design w/ 3 points per whole plot
#'   exact_design <- design$round_split_plot(
#'    m = 3,
#'    hard_to_change = ~A,
#'    n_repeats = 100
#'    )
#'   mirai::daemons(0)
#' }
#'
#' ## Penalized optimal design
#'
#' design$update(
#'   cost = ~A,         ## increase cost of high levels
#'   alpha = 1          ## penalty scaling term
#' )
#' ## initial D-optimal approximate design
#' design$optimize(criterion = "D", max_iter = 5, return_value = FALSE)
#' if(FALSE) {
#'   mirai::daemons(n = parallel::detectCores() / 2, seed = 1)
#'   exact_design <- design$round_split_plot(
#'    m = 3,
#'    hard_to_change = ~A,
#'    n_repeats = 100
#'    )
#'   mirai::daemons(0)
#' }
#'
#' @export
lme_design <- R6Class(
  "lme_design",
  inherit = lm_design,
  public = list(
    formula = NULL,
    data = NULL,
    eta = NULL,
    cost = NULL,
    weights = NULL,
    design_weights = NULL,
    strata = NULL,
    alpha = 1,
    lambda = 0,
    upper = 1,

    #' @description
    #' Constructor. Initializes the experimental design to optimize with a linear mixed model,
    #' discrete design grid, and fixed variance component ratio.
    #' @param formula a two-part \link{formula} of the form \code{~ fix.eff ~ ran.eff}, where \code{fix.eff} is
    #' a linear model formula encoding the fixed effects design matrix and \code{ran.eff} is an expression
    #' encoding the random effects levels. \code{ran.eff} is evaluated as a factor, with the factor levels
    #' determining the random effects groups.
    #' @param eta numeric variance component ratio \eqn{\sigma_u^2 / \sigma_e^2}, where \eqn{\sigma_u^2} is the
    #' (single) variance component of the random effects and \eqn{\sigma_e^2} is the error variance.
    #' @param contrasts optional contrasts used when evaluating the fixed effects formula \code{fix.eff}.
    #' @param ... any additional arguments passed to \link{model.matrix} when evaluating the fixed effects formula \code{fix.eff}.
    initialize = function(formula, data, eta, weights = NULL, cost = NULL, contrasts = NULL, alpha = 1, lambda = 0, upper = 1, ...) {

      self$formula <- formula
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
      stopifnot(is.numeric(eta), length(eta) == 1L, eta >= 0)
      self$eta <- eta
      model_terms <- parse_design_matrices(
        self$formula,
        self$data,
        private$dots
      )
      private$.X <- model_terms$X
      private$Zt <- model_terms$Zt
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
    #' Inputs should be named \code{formula}, \code{data}, \code{eta}, \code{weights}, \code{cost}, \code{alpha}, \code{lambda} and/or \code{upper} respectively.
    #' @param ... named arguments \code{formula}, \code{data}, \code{eta}, \code{weights}, \code{cost}, \code{alpha}, \code{lambda} and/or \code{upper}.
    update = function(...) {

      dots <- list(...)

      if(is.element("formula", names(dots)) && !is.null(dots$formula)) {
        self$formula <- dots$formula
        ## cascade update data
        if(!is.element("data", names(dots))) {
          dots$data <- self$data
        }
      }
      if(is.element("eta", names(dots)) && !is.null(dots$eta)) {
        eta <- dots$eta
        stopifnot(is.numeric(eta), length(eta) == 1L, eta >= 0)
        self$eta <- eta
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
        ## update design matrices
        model_terms <- parse_design_matrices(
          self$formula,
          self$data,
          private$dots
        )
        private$.X <- model_terms$X
        private$Zt <- model_terms$Zt
        ## cascade update cost
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
    #' Optimizes the penalized approximate design based on the selected criterion by iterative
    #' reweighting of a convex approximation of the optimization problem.
    #' \code{"alias"}-optimality is included for completeness only, as the optimality criterion
    #' does not vary with the variance component ratio \eqn{\eta}.
    #' @param max_iter integer number of iterations sequential convex approximation. The default is 1,
    #' which returns the approximate design corresponding to a variance component ratio \eqn{\eta = 0}.
    #' @param tol numeric tolerance used to scale the minimum improvement threshold in the backtracking line
    #' search to determine a new set of approximate design weights.
    #' @param show_progress logical, display a progress bar during the rounding procedure (only if
    #' \code{max_iter > 1}). Defaults to \code{TRUE}.
    optimize = function(criterion = c("D", "A", "I", "G", "alias"), max_iter = 1, tol = .Machine$double.eps^(1/4),
                        scs_control = scs_control_dflt(), verbose = FALSE, return_value = TRUE, show_progress = TRUE) {
      stopifnot(
        is.logical(verbose), is.logical(show_progress), is.logical(return_value),
        is.numeric(max_iter), length(max_iter) == 1, max_iter >= 1,
        is.numeric(tol), length(tol) == 1, tol > 0
      )
      criterion <- match.arg(criterion, c("D", "A", "I", "G", "alias"))
      if(identical(criterion, "alias") || is.na(self$eta) || self$eta <= 0) {
        max_iter <- 1
      }
      if(max_iter > 1 && isTRUE(show_progress)) {
        pb <- progress_bar$new(
          format = "Iter.: :iter crit.: :crit) [:bar] :percent (:elapsed)",
          total = max_iter,
          clear = FALSE,
          show_after = 0
        )
        pb$tick(
          tokens = list(
            iter = 1L,
            crit = "N/A"
          )
        )
      }
      crit <- 0
      warm_start <- NULL
      design_weights <- NULL
      M <- NULL
      sca_terms <- list(X = private$.X)
      for(iter in seq_len(max_iter)) {
        ## optimize
        design <- optimize_impl(
          X = sca_terms$X,
          criterion = criterion,
          cost = self$cost,
          weights = self$weights,
          strata = self$strata,
          control = scs_control,
          alpha = self$alpha,
          lambda = self$lambda,
          upper = self$upper,
          start = warm_start,
          verbose = verbose,
          return_sol = isTRUE(iter > 1),
          v_sca = sca_terms$v
        )
        if(suppressWarnings(check_status(design))) {
          if(iter > 1) {
            Minv <- tryCatch(chol2inv(chol(M)), error = function(e) 0)  ## can M be singular?
            ## backtracking line search
            for(alpha in 0.5^(0:9)) {
              new_weights <- design_weights + alpha * (design$w - design_weights)
              M1 <- M_eta(w = pmax(new_weights, 0), X = private$.X, Z = t(private$Zt), eta = self$eta, lambda = self$lambda)
              M1_approx <- M_approx_impl(w = new_weights, sca_terms = sca_terms)
              if(!any(is.na(M1))) {
                crit_new <- switch(criterion,
                                   A = a_crit(w = new_weights, X = private$.X, M = M1, cost = self$cost, alpha = self$alpha),
                                   I = i_crit(w = new_weights, X = private$.X, M = M1, cost = self$cost, alpha = self$alpha),
                                   D = d_crit(w = new_weights, X = private$.X, M = M1, cost = self$cost, alpha = self$alpha),
                                   G = g_crit(w = new_weights, X = private$.X, M = M1, cost = self$cost, alpha = self$alpha),
                                   alias = alias_crit(w = new_weights, X = private$.X, M = M1, cost = self$cost, alpha = self$alpha)
                )
              } else {
                crit_new <- 0
              }
              ftol <- tol * alpha * sum(t(Minv) * (M1_approx - M))   ## Armijo condition
              if(crit_new > crit + ftol) break
            }
            if(crit_new > crit) {
              crit <- crit_new
              design_weights <- new_weights
              warm_start <- design$scs_solution     ## use only if sca_terms is defined
              M <- M1
            } else {
              break
            }
          } else {
            design_weights <- design$w
            M <- M_eta(w = pmax(design_weights, 0), X = private$.X, Z = t(private$Zt), eta = self$eta, lambda = self$lambda)
            crit <- switch(criterion,
                           A = a_crit(w = design_weights, X = private$.X, M = M, cost = self$cost, alpha = self$alpha),
                           I = i_crit(w = design_weights, X = private$.X, M = M, cost = self$cost, alpha = self$alpha),
                           D = d_crit(w = design_weights, X = private$.X, M = M, cost = self$cost, alpha = self$alpha),
                           G = g_crit(w = design_weights, X = private$.X, M = M, cost = self$cost, alpha = self$alpha),
                           alias = alias_crit(w = design_weights, X = private$.X, M = M, cost = self$cost, alpha = self$alpha)
            )
          }
          ## update sca matrices
          if(iter < max_iter) {
            sca_terms <- M_sca_terms(w = pmax(design_weights, 0), X = private$.X, Z = t(private$Zt), eta = self$eta)
          }
        } else {
          break
        }
        if(max_iter > 1 && isTRUE(show_progress) && iter < max_iter) {
          pb$tick(
            tokens = list(
              iter = iter + 1,
              crit = sprintf("%.3g", crit)
            )
          )
        }
      }
      self$design_weights <- structure(
        design_weights,
        criterion = list(
          crit.value = crit,
          crit = criterion
        )
      )
      if(isTRUE(return_value)) {
        return(self$design_weights)
      }
    },

    #' @description
    #' Optimizes the augmented penalized approximate design given a vector of initial design weights
    #' by iterative reweighting of a convex approximation of the optimization problem.
    #' \code{"alias"}-optimality is included for completeness only, as the optimality criterion
    #' does not vary with the variance component ratio \eqn{eta}.
    #' @param design_weights a numeric vector with initial design weights, must be the same length as the number of design points in \code{data}.
    #' Also accepts a one-sided \link{formula} made up of column names present in \code{data}.
    #' @param max_iter integer number of iterations sequential convex approximation. The default is 1,
    #' which returns the approximate design corresponding to a variance component ratio \eqn{\eta = 0}.
    #' @param gamma numeric contribution weight of the initial design in the augmented design problem, see also \sQuote{Augmented optimal design}
    #' above. Defaults to 0.5.
    #' @param tol numeric tolerance used to scale the minimum improvement threshold in the backtracking line
    #' search to determine a new set of approximate design weights.
    #' @param show_progress logical, display a progress bar during the rounding procedure (only if
    #' \code{max_iter > 1}). Defaults to \code{TRUE}.
    augment = function(design_weights, max_iter = 1, gamma = 0.5, criterion = c("D", "A", "I", "G", "alias"), tol = .Machine$double.eps^(1/4),
                       scs_control = scs_control_dflt(), verbose = FALSE, return_value = TRUE, show_progress = TRUE) {

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
        any(design_weights > 0),
        is.logical(verbose), is.logical(show_progress), is.logical(return_value),
        is.numeric(max_iter), length(max_iter) == 1, max_iter >= 1,
        is.numeric(tol), length(tol) == 1, tol > 0
      )
      gamma <- as.numeric(gamma)
      w0 <- design_weights / sum(design_weights)
      criterion <- match.arg(criterion, c("D", "A", "I", "G", "alias"))
      if(identical(criterion, "alias") || is.na(self$eta) || self$eta <= 0) {
        max_iter <- 1
      }
      if(max_iter > 1 && isTRUE(show_progress)) {
        pb <- progress_bar$new(
          format = "Iter.: :iter crit.: :crit) [:bar] :percent (:elapsed)",
          total = max_iter,
          clear = FALSE,
          show_after = 0
        )
        pb$tick(
          tokens = list(
            iter = 1L,
            crit = "N/A"
          )
        )
      }
      crit <- 0
      warm_start <- NULL
      design_weights <- NULL
      M <- NULL
      sca_terms <- list(X = private$.X)
      for(iter in seq_len(max_iter)) {
        ## optimize
        design <- optimize_impl(
          X = sca_terms$X,
          criterion = criterion,
          cost = self$cost,
          weights = self$weights,
          strata = self$strata,
          control = scs_control,
          alpha = self$alpha,
          lambda = self$lambda,
          upper = self$upper,
          gamma = gamma,
          w0 = w0,
          start = warm_start,
          verbose = verbose,
          return_sol = isTRUE(iter > 1),
          v_sca = sca_terms$v,
        )
        if(suppressWarnings(check_status(design))) {
          if(iter > 1) {
            Minv <- tryCatch(chol2inv(chol(M)), error = function(e) 0)  ## can M be singular?
            ## backtracking line search
            for(alpha in 0.5^(0:9)) {
              new_weights <- design_weights + alpha * (design$w - design_weights)
              M1 <- M_eta(w = pmax(new_weights, 0), X = private$.X, Z = t(private$Zt), eta = self$eta, lambda = self$lambda)
              M1_approx <- M_approx_impl(w = new_weights, sca_terms = sca_terms)
              if(!any(is.na(M1))) {
                crit_new <- switch(criterion,
                                   A = a_crit(w = new_weights, X = private$.X, M = M1, cost = self$cost, alpha = self$alpha),
                                   I = i_crit(w = new_weights, X = private$.X, M = M1, cost = self$cost, alpha = self$alpha),
                                   D = d_crit(w = new_weights, X = private$.X, M = M1, cost = self$cost, alpha = self$alpha),
                                   G = g_crit(w = new_weights, X = private$.X, M = M1, cost = self$cost, alpha = self$alpha),
                                   alias = alias_crit(w = new_weights, X = private$.X, M = M1, cost = self$cost, alpha = self$alpha)
                )
              } else {
                crit_new <- 0
              }
              ftol <- tol * alpha * sum(t(Minv) * (M1_approx - M))   ## Armijo condition
              if(crit_new > crit + ftol) break
            }
            if(crit_new > crit) {
              crit <- crit_new
              design_weights <- new_weights
              warm_start <- design$scs_solution     ## use only if sca_terms is defined
              M <- M1
            } else {
              break
            }
          } else {
            design_weights <- design$w
            M <- M_eta(w = pmax(design_weights, 0), X = private$.X, Z = t(private$Zt), eta = self$eta, lambda = self$lambda)
            crit <- switch(criterion,
                           A = a_crit(w = design_weights, X = private$.X, M = M, cost = self$cost, alpha = self$alpha),
                           I = i_crit(w = design_weights, X = private$.X, M = M, cost = self$cost, alpha = self$alpha),
                           D = d_crit(w = design_weights, X = private$.X, M = M, cost = self$cost, alpha = self$alpha),
                           G = g_crit(w = design_weights, X = private$.X, M = M, cost = self$cost, alpha = self$alpha),
                           alias = alias_crit(w = design_weights, X = private$.X, M = M, cost = self$cost, alpha = self$alpha)
            )
          }
          ## update sca matrices
          if(iter < max_iter) {
            sca_terms <- M_sca_terms(w = pmax(design_weights, 0), X = private$.X, Z = t(private$Zt), eta = self$eta)
          }
        } else {
          break
        }
        if(max_iter > 1 && isTRUE(show_progress) && iter < max_iter) {
          pb$tick(
            tokens = list(
              iter = iter + 1,
              crit = sprintf("%.3g", crit)
            )
          )
        }
      }
      self$design_weights <- structure(
        design_weights,
        criterion = list(
          crit.value = crit,
          crit = criterion
        )
      )
      if(isTRUE(return_value)) {
        return(self$design_weights)
      }
    },

    #' @description
    #' Round penalized approximate design to an exact (integer) design.
    #' The approximate design is initialized from the \code{design_weights} field, which is the optimal (approximate) design calculated by \code{$optimize()}.
    #' For method \code{"optimal"}, cost penalties and upper limit constraints present in \code{$cost} and \code{$upper}
    #' (if specified) are included in the criterion optimized by the point-exchange procedure.
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
    #' @param m integer number of design points included in the exact design.
    #' @param method the method used to calculate the exact design. The following choices are supported:
    #' \tabular{l}{
    #'  \code{"optimal"} optimality criterion point-exchange algorithm (default).\cr
    #'  \code{"efficient"} efficient rounding method of Pukelsheim & Reider (1992).
    #' }
    #' @param max_iter integer number of maximum iterations in point-exchange algorithms. Only used when
    #' \code{method = "optimal"}.
    #' @param n_repeats integer number of (random) repeats of the rounding/exchange algorithm. The exact design
    #' with the highest criterion value is returned. Defaults to 10.
    #' @param replicates logical, allow replicate design points in the exact design, defaults to \code{TRUE}.
    #' Can also be an integer number, limiting the maximum number of allowed replicates per design point.
    #' This argument overrides the \code{$upper} value constraint if existing.
    #' @param tol numeric tolerance used as lower threshold to continue iterations in point-exchange algorithms.
    #' @param seed (optional) random seed used to generate jitter to break ties in rounding the initial exact design.
    #' Not used in the case of (\code{mirai}) asynchronous evaluation, which requires the random seed to be set via the
    #' \code{seed} argument in \code{mirai::daemons}.
    #' @param design_weights (optional) vector of design weights used as initial approximate design, with length equal
    #' to the number of available design points. If \code{TRUE} (default), the approximate design is initialized from the
    #' \code{design_weights} field. If \code{FALSE}, the initial approximate design is generated by random sampling.
    #' @param augment_design (optional) integer design to augment, with length equal to the number of available design points.
    #' @param show_progress logical, display a progress bar during the rounding procedure (only if
    #' \code{n_repeats > 1}). Defaults to \code{TRUE}.
    #' @param ... additional arguments to tune algorithm/solver settings:
    #' \itemize{
    #' \item \itemize{
    #' \item \code{criterion}, character string specifying the optimality criterion for the rounding procedure in case
    #' \code{design_weights} has no \code{"criterion"} attribute, or to override the criterion specified in the \code{"criterion"} attribute.
    #' If specified, must be one of \code{"D"}, \code{"A"}, \code{"I"}, \code{"G"} or \code{"alias"}.}
    #' }
    round = function(m, method = c("optimal", "efficient"), max_iter = 100, n_repeats = 10, replicates = TRUE,
                     tol = .Machine$double.eps^(1/3), seed = NULL, design_weights = TRUE, augment_design =  NULL, show_progress = TRUE, ...) {
      ## check arguments
      method <- match.arg(method, c("optimal", "efficient"))
      dots <- list(...)
      if(isTRUE(design_weights)) {
        design_weights <- self$design_weights
        if(is.null(design_weights)) {
          stop("no initial design available to round. Call `$optimize()` first to initialize an approximate optimal design.")
        }
      } else if(isFALSE(design_weights)) {
        design_weights <- NULL
      } else {
        design_weights <- design_weights / sum(design_weights)
      }
      if(method == "efficient" && is.null(design_weights)) {
        stop("method 'efficient' can only be used if approximate design weights are available.")
      }
      if(!is.null(augment_design)) {
        stopifnot(
          is.numeric(augment_design), length(augment_design) == nrow(private$.X)
        )
        augment_design <- round(augment_design)
      }
      stopifnot(
        is.null(design_weights) || is.numeric(design_weights), is.null(design_weights) || length(design_weights) == nrow(private$.X),
        is.numeric(m), length(m) == 1, m > 0,
        is.numeric(tol), length(tol) == 1, tol > 0,
        is.numeric(max_iter), length(max_iter) == 1, max_iter > 0
      )
      if(!is.null(dots$criterion)) {
        crit <- match.arg(dots$criterion, c("I", "A", "D", "G", "alias"))
      } else if(!is.null(attr(design_weights, "criterion"))) {
        crit0 <- attr(design_weights, "criterion")
        crit <- match.arg(crit0$crit, c("I", "A", "D", "G", "alias"))
        attr(crit, "crit0") <- crit0
      } else {
        stop("no 'criterion' attribute found in 'design_weights'. Specify the 'criterion' argument to set an optimality criterion for the rounded design.")
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
          stopifnot(is.numeric(replicates), length(replicates) == 1, replicates >= 1)
          replicates <- min(replicates, nrow(X))
        }
        upper <- replicates / m
      } else {
        upper <- self$upper
      }
      if(upper < 1 && nrow(X) * max(floor(m * upper), 1) < m) {
        warning("upper limit constraint in `$upper` cannot be satisfied, reduce the number of design points 'm'.")
      }
      ## round strata proportions
      if(!is.null(self$strata)) {
        if(identical(method, "efficient")) {
          warning("stratification constraints are not optimized when using method 'efficient'.")
        }
        strata <- self$strata
        if(m < sum(strata$proportions > 0)) {
          warning("'m' is smaller than the number of individual strata, increase 'm' to (minimally) cover all strata")
        }
        strata$proportions <- round_puk(strata$proportions, m)
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
      if(!is.null(design_weights)) {
        supp_tol <- ifelse(crit == "alias" && method != "efficient", 0, (.Machine$double.eps)^(1/4))
        supp <- design_weights > supp_tol
        if(!is.null(augment_design)) {
          supp <- supp & (augment_design < .Machine$double.eps)
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
          w <- round_eta_impl(
            design_weights = design_weights,
            X = X,
            Zt = private$Zt,
            eta = self$eta,
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
            alpha = self$alpha
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
          round_eta_impl = round_eta_impl,
          design_weights = design_weights,
          X = X,
          Zt = private$Zt,
          eta = self$eta,
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
          alpha = self$alpha
        )
        m_list <- lapply(1:n_repeats, function(n) {
          mirai::mirai(
            .expr = {
              round_eta_impl(
                design_weights, X, Zt, eta, m, n_supp, supp_tol, supp, tol, max_iter,
                crit, method, augment_design, upper, strata, cost, lambda, alpha
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
      return(w_opt)
    },

    #' @description
    #' Round penalized approximate design to an exact (integer) split-plot design.
    #' This rounding procedure enforces fixed levels of the defined \code{hard_to_change} factors within each \emph{whole plot}, where the whole plots
    #' are determined by the random effects factors in the model \code{$formula}.
    #' Stratification constraints imposed with \code{$stratify()} are not taken into account as these are overridden by the split-plot design constraints.
    #' The approximate design is initialized from the \code{design_weights} field, which is the optimal (approximate) design calculated by \code{$optimize()}.
    #' To override the initial approximate (optimal) design, a vector of design weights can be passed using the \code{design_weights} argument.
    #' The vector must be the same length as the number of design points in \code{data}.
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
    #'
    #' @param m integer whole plot size. The number of design points per level of the whole plot factor(s). The total number of design points included in
    #' the split-plot design is \code{m} times the number of whole plot levels.
    #' @param hard_to_change hard-to-change \link{factor} variable of the same length as the number of rows in \code{data}. Also accepts
    #' a one-sided \link{formula} made up of column names present in \code{data}, in which case the factor levels are defined by the \link{interaction}
    #' between terms in the \link{model.frame} obtained from the \link{formula}.
    #' @param max_iter integer number of maximum iterations in point-exchange algorithms.
    #' @param n_repeats integer number of (random) repeats of the rounding/exchange algorithm. The exact design
    #' with the highest criterion value is returned. Defaults to 10.
    #' @param tol numeric tolerance used a lower threshold to continue iterations in point-exchange algorithms.
    #' @param replicates logical, allow replicate design points in the exact design, defaults to \code{TRUE}.
    #' Can also be an integer number, limiting the maximum number of allowed replicates per design point.
    #' This argument overrides the \code{$upper} value constraint if existing.
    #' @param seed (optional) random seed used to generate jitter to break ties in rounding the initial exact design.
    #' Not used in the case of (\code{mirai}) asynchronous evaluation, which requires the random seed to be set via the
    #' \code{seed} argument in \code{mirai::daemons}.
    #' @param design_weights (optional) vector of design weights used as initial approximate design.
    #' @param augment_design (optional) integer design to augment, with length equal to the number of available design points.
    #' @param show_progress logical, display a progress bar monitoring the repeats of the rounding/exchange procedure. Defaults to \code{TRUE}.
    #' @param ... additional arguments to tune algorithm settings:
    #' \itemize{
    #' \item \itemize{
    #' \item \code{criterion}, character string specifying the optimality criterion for the rounding procedure in case
    #' \code{design_weights} has no \code{"criterion"} attribute, or to override the criterion specified in the \code{"criterion"} attribute.
    #' If specified, must be one of \code{"D"}, \code{"A"}, \code{"I"}, or \code{"G"}.
    #' }
    #' }
    round_split_plot = function(m, hard_to_change, max_iter = 100, n_repeats = 10, replicates = TRUE, tol = .Machine$double.eps^(1/3), seed = NULL,
                                design_weights = NULL, augment_design = NULL, show_progress = TRUE, ...) {
      ## check arguments
      dots <- list(...)
      if(!is.null(self$strata)) {
        warning("stratification constraints are present, but are not taken into account when rounding to a split-plot design.")
        strata0 <- self$strata
        on.exit(self$strata <- strata0, add = TRUE, after = FALSE)
      }
      if(is.null(design_weights)) {
        design_weights <- self$design_weights
        if(is.null(design_weights)) {
          stop("No initial approximate design available to round. Call `$optimize()` first to initialize an approximate optimal design.")
        }
      } else {
        design_weights <- design_weights / sum(design_weights)
      }
      if(!is.null(augment_design)) {
        stopifnot(
          is.numeric(augment_design), length(augment_design) == length(design_weights)
        )
        augment_design <- round(augment_design)
      }
      stopifnot(
        is.numeric(design_weights), length(design_weights) == nrow(private$.X),
        is.numeric(m), length(m) == 1, m > 0,
        is.numeric(tol), length(tol) == 1, tol > 0,
        is.numeric(max_iter), length(max_iter) == 1, max_iter > 0
      )
      if(!is.null(dots$criterion)) {
        crit <- match.arg(dots$criterion, c("I", "A", "D", "G", "alias"))
      } else if(!is.null(attr(design_weights, "criterion"))) {
        crit0 <- attr(design_weights, "criterion")
        crit <- match.arg(crit0$crit, c("I", "A", "D", "G", "alias"))
        attr(crit, "crit0") <- crit0
      } else {
        stop("no 'criterion' attribute found in 'design_weights'. Specify the 'criterion' argument to set an optimality criterion for the rounded design.")
      }
      mw <- as.integer(nrow(private$Zt))
      wsize <- as.integer(m)
      m <- mw * wsize
      max_iter <- as.integer(max_iter)
      n_repeats <- as.integer(n_repeats)
      if(!is.null(self$weights)) {
        X <- sqrt(self$weights) * private$.X
      } else {
        X <- private$.X
      }
      ## determine upper limit
      upper <- 1
      if(!isTRUE(replicates)) {
        if(isFALSE(replicates)) {
          replicates <- 1
        } else {
          stopifnot(is.numeric(replicates), length(replicates) == 1, replicates >= 1)
          replicates <- min(replicates, nrow(X))
        }
        upper <- replicates / m
      } else {
        upper <- self$upper
      }
      if(upper < 1 && nrow(X) * max(floor(m * upper), 1) < m) {
        warning("upper limit constraint in `$upper` cannot be satisfied, reduce the number of design points 'm'.")
      }
      ## approximate strata proportions
      Zind <- as.integer(crossprod(private$Zt > 0, seq_len(nrow(private$Zt))))
      if(!missing(hard_to_change)) {
        if(inherits(hard_to_change, "formula")) {
          mf <- do.call(model.frame, args = c(list(formula = hard_to_change, data = self$data), private$dots))
          strata <- do.call(interaction, args = c(as.list(mf), drop = TRUE))
        } else {
          strata <- hard_to_change
        }
        stopifnot(
          is.factor(strata), nlevels(strata) > 0, length(strata) == nrow(self$data)
        )
        strata <- list(
          strata = strata,
          proportions = tapply(design_weights, strata, sum)
        )
      } else {
        hard_to_change <- NULL
        strata <- list(
          strata = factor(Zind),
          proportions = rep(wsize, times = max(Zind))
        )
      }
      w_opt <- NULL
      crit_opt <- -Inf
      strata_opt <- Inf
      supp_tol <- ifelse(crit == "alias", 0, (.Machine$double.eps)^(1/4))
      if(isTRUE(show_progress) && n_repeats > 1) {
        pb <- progress_bar$new(
          format = sprintf("Iter: :n (%s-crit: :crit) [:bar] :percent (:elapsed)", crit),
          total = n_repeats,
          clear = FALSE,
          show_after = 0
        )
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
          init <- init_split_plot_strata(
            design_weights = design_weights,
            X = X,
            Zind = Zind,
            m = m,
            wsize = wsize,
            supp_tol = supp_tol,
            hard_to_change = hard_to_change,
            strata = strata,
            augment_design = augment_design,
            upper = upper
          )
          if(is.null(init)) {
            next
          }
          w <- round_eta_impl(
            design_weights = init$design_weights,
            X = init$X,
            Zt = private$Zt[, init$strata$indices, drop = FALSE],
            eta = self$eta,
            m = m,
            n_supp = sum(init$supp),
            supp_tol = supp_tol,
            supp = init$supp,
            tol = tol,
            max_iter = max_iter,
            crit = crit,
            method = "optimal",
            augment_design = init$augment_design,
            upper = upper,
            strata = init$strata$strata,
            cost = self$cost[init$strata$indices],
            lambda = self$lambda,
            alpha = self$alpha,
            split_plot = TRUE
          )
          ## update criterion value
          strata_diff <- sum(attr(w, "strata") %% wsize)
          if(strata_diff > strata_opt) {
            crit_val <- -Inf
          } else {
            strata_opt <- strata_diff
            w0 <- numeric(length(design_weights))
            w0[init$strata$indices] <- w
            w <- update_crit_eta(w0 = design_weights, w1 = w0, X = X, Z = t(private$Zt), eta = self$eta, Zind = Zind, crit = crit,
                                 cost = self$cost, lambda = self$lambda, alpha = self$alpha, strata = strata)
            crit_val <- attr(w, "criterion")$crit.value
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
          init_split_plot_strata = init_split_plot_strata,
          round_eta_impl = round_eta_impl,
          design_weights = design_weights,
          X = X,
          Zt = private$Zt,
          Zind = Zind,
          eta = self$eta,
          m = m,
          wsize = wsize,
          supp_tol = supp_tol,
          hard_to_change = hard_to_change,
          tol = tol,
          max_iter = max_iter,
          crit = crit,
          augment_design = augment_design,
          upper = upper,
          strata = strata,
          cost = as.vector(self$cost),
          lambda = self$lambda,
          alpha = self$alpha
        )
        m_list <- lapply(1:n_repeats, function(n) {
          mirai::mirai(
            .expr = {
              init <- init_split_plot_strata(
                design_weights = design_weights,
                X = X,
                Zind = Zind,
                m = m,
                wsize = wsize,
                supp_tol = supp_tol,
                hard_to_change = hard_to_change,
                strata = strata,
                augment_design = augment_design,
                upper = upper
              )
              if(is.null(init)) {
                return(NULL)
              }
              w <- round_eta_impl(
                design_weights = init$design_weights,
                X = init$X,
                Zt = Zt[, init$strata$indices, drop = FALSE],
                eta = eta,
                m = m,
                n_supp = sum(init$supp),
                supp_tol = supp_tol,
                supp = init$supp,
                tol = tol,
                max_iter = max_iter,
                crit = crit,
                method = "optimal",
                augment_design = init$augment_design,
                upper = upper,
                strata = init$strata$strata,
                cost = cost[init$strata$indices],
                lambda = lambda,
                alpha = alpha,
                split_plot = TRUE
              )
              return(list(w = w, indices = init$strata$indices))
            }
          )
        })
        while (length(m_list) > 0) {
          m_id <- mirai::race_mirai(m_list)
          m_data <- m_list[[m_id]]$data
          m_list[[m_id]] <- NULL
          if(!is.null(m_data)) {
            strata_diff <- sum(attr(m_data$w, "strata") %% wsize)
            if(strata_diff > strata_opt) {
              crit_val <- -Inf
            } else {
              strata_opt <- strata_diff
              w0 <- numeric(length(design_weights))
              w0[m_data$indices] <- m_data$w
              w <- update_crit_eta(w0 = design_weights, w1 = w0, X = X, Z = t(private$Zt), eta = self$eta, Zind = Zind, crit = crit,
                                   cost = self$cost, lambda = self$lambda, alpha = self$alpha, strata = strata)
              crit_val <- attr(w, "criterion")$crit.value
            }
            if(crit_val > crit_opt) {
              w_opt <- w
              crit_opt <- crit_val
            }
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
      return(w_opt)
    },

    #' @description
    #' Evaluate optimality criterion of a penalized approximate or rounded (exact) design.
    #' If cost penalties are present in the \code{$cost} field, the penalized optimality criterion as described in the \sQuote{Optimality criteria} section of
    #' \link{lm_design} are returned.
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
      M <- M_eta(w = pmax(w, 0), X = X, Z = t(private$Zt), eta = self$eta, lambda = lambda)
      crit <- switch(criterion,
                     A = a_crit(w = w, X = X, M = M, cost = cost, alpha = alpha),
                     I = i_crit(w = w, X = X, M = M, cost = cost, alpha = alpha),
                     D = d_crit(w = w, X = X, M = M, cost = cost, alpha = alpha),
                     G = g_crit(w = w, X = X, M = M, cost = cost, alpha = alpha),
                     alias = alias_crit(w = w, X = X, M = M, cost = cost, alpha = alpha)
      )
      return(crit)
    },

    #' @description
    #' Calculates the variance-covariance matrix of the model parameters based on an approximate or rounded (exact) design \eqn{(w_1, \ldots, w_n)} according to:
    #' \deqn{\Sigma \ = \ (X^T\Omega^{1/2}V^{-1}\Omega^{1/2}X)^{-1}}
    #' Here \eqn{\Omega^{1/2} = \mathrm{diag}(\sqrt{w_1v_1}, \ldots, \sqrt{w_nv_n})} and \eqn{V^{-1}} is the inverse of the block-diagonal covariance matrix \eqn{V = \mathrm{diag}[V_1,\ldots,V_S]},
    #' consisting of the level-specific blocks \eqn{V_s = \sigma_e^2 \left( I_{n_s} + \eta \boldsymbol{1}_{n_s} \boldsymbol{1}_{n_s}^T \right) },
    #' with \eqn{n_s} is the number of observations for random effects level \eqn{s}.
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
      M <- M_eta(w = w * m, X = X, Z = t(private$Zt), eta = self$eta, scale = FALSE)
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
    #' @param newdata (optional) data frame or matrix with design points (rows) for which to return the standard error of prediction and which
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
        frm <- as.formula(self$formula[[2]])
        xpred <- do.call(model.matrix, args = c(list(object = frm, data = newdata), private$dots))
      } else {
        xpred <- private$.X
      }
      Sigma <- self$vcov(m, design_weights, sigma2)
      ses <- sqrt(rowSums(xpred * (xpred %*% Sigma)))
      return(ses)
    },

    #' @description
    #' Calculates the correlation matrix of the design matrix columns based on an approximate or rounded (exact) design \eqn{(w_1, \ldots, w_n)} according to:
    #' \deqn{R \ = \ D^{-1/2}(X_c^T \Omega^{1/2} V^{-1} \Omega^{1/2} X_c)D^{-1/2}}
    #' where \eqn{\Omega = \mathrm{diag}(\sqrt{w_1v_1}, \ldots, \sqrt{w_nv_n})} and \eqn{V^{-1}} is the inverse of the block-diagonal covariance matrix \eqn{V = \mathrm{diag}[V_1,\ldots,V_S]}.
    #' \eqn{X_c} is the design matrix \eqn{X} centered by the (weighted) column means:
    #' \deqn{\frac{\sum_{i=1}^n \omega_i X_{ij}}{\sum_{i=1}^n \omega_i}, \quad j = 1,\ldots, p}
    #' \eqn{D^{-1/2}} is \eqn{\mathrm{diag}(1/\sqrt{D_1},\ldots, 1/\sqrt{D_p})}, with \eqn{D_1,\ldots,D_p} the diagonal elements of \eqn{X_c^T \Omega^{1/2} V^{-1} \Omega^{1/2} X_c}.
    #' If \code{center = FALSE}, the uncentered correlation matrix is returned, which corresponds to the scale-invariant (unit diagonal) version of the information matrix \eqn{M(w | X, v)}.
    #' The design weights are initialized from the \code{design_weights} field, which is the optimal (approximate) design calculated by \code{$optimize()}.
    #' To override the default approximate (optimal) design weights, a vector of design weights can be passed using the \code{design_weights} argument.
    #' The vector must be the same length as the number of design points in \code{data}.
    #' @param design_weights (optional) vector of design weights for which to evaluate the information matrix.
    #' @param center logical, whether to center the design columns when evaluating the correlation matrix, defaults to \code{TRUE}.
    #' @param unscaled logical value, whether to include \eqn{V^{-1}} in the calculation of \eqn{R}. If \code{TRUE}, the ordinary
    #' correlation matrix (without random effects) \eqn{R \ = \ D^{-1/2}(X_c^T \Omega X_c)D^{-1/2}} is returned.
    corr = function(design_weights = NULL, center = TRUE, unscaled = FALSE) {
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
        is.numeric(w), length(w) == nrow(private$.X), any(w > 0)
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
      } else {
        Xc <- X
      }
      M <- M_eta(w = omega, X = Xc, Z = t(private$Zt), eta = ifelse(!isTRUE(unscaled), self$eta, 0), scale = FALSE)
      D <- sqrt(pmax(diag(M), .Machine$double.eps))
      R <- M / outer(D, D)
      rownames(R) <- colnames(R) <- colnames(X)
      return(R)
    }
  ),

  active = list(
    X = function(value) {
      if (missing(value)) {
        private$.X
      } else {
        stop("To modify the fixed effects design matrix, update the data or model formula with `$update()`.", call. = FALSE)
      }
    },
    Z = function(value) {
      if(missing(value)) {
        t(as.matrix(private$Zt))
      } else {
        stop("To modify the random effects design matrix, update the data or model formula with `$update()`.", call. = FALSE)
      }
    }
  ),

  private = list(
    dots = list(),
    .X = NULL,
    Zt = NULL,
    contrasts = NULL
  )

)

parse_design_matrices <- function(formula, data, kwargs) {
  ## fixed effects
  frm <- as.formula(formula[[2]])
  X <- do.call(model.matrix, args = c(list(object = frm, data = data), kwargs))
  ## random effects
  re_values <- tryCatch(eval(formula[[3]], envir = as.list(data)))
  if(inherits(re_values, "error")) {
    stop(sprintf("failed to evaluate random effects epxression: %s", re_values$message))
  }
  re_nm <- "._ran.eff"
  if(is.element(re_nm, names(data))) {
    ## probably never happens
    nms_uniq <- make.unique(c(names(data), re_nm))
    re_nm <- nms_uniq[length(nms_uniq)]
  }
  frm <- update(frm, call("~", call("+", quote(.), call("+", 1, as.name(re_nm)))))
  environment(frm) <- environment(formula)
  data_frm <- as.list(data)
  data_frm[[re_nm]] <- re_values
  mf <- do.call(model.frame, args = c(list(formula = frm, data = data_frm), drop.unused.levels = TRUE))
  attr(mf, "formula") <- frm
  ## create random effects terms
  re_terms <- mkReTrms(list(call("|", 1, as.name(re_nm))), mf)
  return(list(X = X, Zt = re_terms$Zt))
}

M_sca_terms <- function(w, X, Z, eta) {
  ## initialize
  Zind <- as.integer((Z > 0) %*% seq_len(ncol(Z)))
  wX <- w * X
  zeta <- 1 / (1 / eta + colSums(w * Z))
  Xs <- array(dim = dim(X))
  v <- matrix(nrow = ncol(X), ncol = ncol(Z))
  ## populate matrices
  for(j in 1:ncol(Z)) {
    ids <- which(Zind == j)
    if(any(ids)) {
      Sj <- colSums(wX[ids, , drop = FALSE])
      v[, j] <- zeta[j] / sqrt(eta) * Sj
      for(m in ids) {
        Xs[m, ] <- X[m, ] - sqrt(eta) * v[, j]
      }
    }
  }
  return(list(X = Xs, v = v))
}

M_approx_impl <- function(w, sca_terms) {
  w <- pmax(w, 0)
  M <- crossprod(sqrt(w) * sca_terms$X)
  M <- M + tcrossprod(sca_terms$v)
  return(M)
}

sample_strata_impl <- function(strata, Zind, wsize) {
  prop <- sample.int(
    n = length(strata$proportions),
    size = max(Zind),
    replace = TRUE,
    prob = pmax(strata$proportions, 0)
  )
  prop <- tabulate(prop, nbins = nlevels(strata$strata))
  blocks <- sample(rep(levels(strata$strata), times = prop), size = sum(prop), replace = FALSE)
  blocks <- factor(blocks, levels = levels(strata$strata))
  indices <- which(strata$strata == blocks[Zind])
  strata <- list(strata = interaction(strata$strata[indices], Zind[indices], drop = TRUE))
  strata$proportions <- rep(wsize, nlevels(strata$strata))
  return(list(strata = strata, indices = indices))
}

init_split_plot_strata <- function(design_weights, X, Zind, m, wsize, supp_tol, hard_to_change, strata,
                                   augment_design = NULL, upper = 1) {

  if(!is.null(hard_to_change)) {
    strata_n <- sample_strata_impl(
      strata = strata,
      Zind = Zind,
      wsize = wsize
    )
  } else {
    strata_n <- list(
      strata = strata,
      indices = seq_len(nrow(X))
    )
  }
  X_n <- X[strata_n$indices, , drop = FALSE]
  design_weights_n <- design_weights[strata_n$indices]
  if(upper < 1 && nrow(X_n) * max(floor(m * upper), 1) < m) {
    return(NULL)
  }
  ## initial rounding
  supp <- design_weights_n > supp_tol
  if(!is.null(augment_design)) {
    augment_design_n <- augment_design[strata_n$indices]
    supp <- supp & (augment_design_n < .Machine$double.eps)
  } else {
    augment_design_n <- NULL
  }
  return(
    list(
      design_weights = design_weights_n,
      X = X_n,
      strata = strata_n,
      augment_design = augment_design_n,
      supp = supp
    )
  )
}

round_eta_impl <- function(design_weights, X, Zt, eta, m, n_supp, supp_tol, supp, tol, max_iter, crit, method, augment_design = NULL,
                           upper = 1, strata = NULL, cost = NULL, lambda = 0, alpha = 1, split_plot = FALSE) {

  ## initial rounding
  if(isTRUE(split_plot) || is.null(design_weights) || m < n_supp || crit == "alias") {
    ## incomplete support
    w <- round_random(design_weights, nrow(X), m, tol = supp_tol, prob_min = ifelse(crit == "alias", 1 / length(design_weights), 0),
                      augment = augment_design, upper = upper)
  } else {
    ## complete support
    w <- round_puk(design_weights, m, augment = augment_design, upper = upper)
  }
  ## exchange optimization
  if(method == "optimal") {
    ## unweighted exchange (fast)
    w <- exchange_impl(w, X, ifelse(crit == "G", "D", crit), supp = supp, strata = strata, max_iter = max_iter,
                       tol = tol, augment = augment_design, cost = cost, lambda = lambda, alpha = alpha,
                       upper = upper, add_crit = FALSE)

    ## weighted exchange (slow)
    if(!is.null(augment_design)) {
      w <- w - augment_design
    }
    if(!is.null(attr(crit, "crit0"))) {
      w <- structure(w, criterion = attr(crit, "crit0"))
      attr(crit, "crit0") <- NULL
    }
    w <- exchange_impl(w, X, crit, strata = strata, max_iter = max_iter, tol = tol, augment = augment_design,
                       cost = cost, lambda = lambda, alpha = alpha, upper = upper, eta = eta, Zt = Zt, add_crit = !split_plot)
  } else if(method == "efficient") {
    if(!is.null(augment_design)) {
      w <- w + augment_design
      w <- update_crit_eta(design_weights, w, X = X, Z = t(Zt), eta = eta, crit = crit, cost = cost,
                           lambda = lambda, alpha = alpha, strata = strata)
      attr(w, "criterion")$crit.rel.eff <- NA
    } else {
      w <- update_crit_eta(design_weights, w, X = X, Z = t(Zt), eta = eta, crit = crit, cost = cost,
                           lambda = lambda, alpha = alpha, strata = strata)
    }
  }

  return(w)
}
