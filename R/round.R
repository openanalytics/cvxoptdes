round_random <- function(w, n, m, tol = (.Machine$double.eps)^(1/4), prob_min = 0, augment = NULL, upper = 1) {
  if(!is.null(w)) {
    supp <- (w > tol)
  } else {
    w <- rep(1 / n, n)
    supp <- rep(TRUE, n)
  }
  n_supp <- sum(supp)
  w_supp <- w[supp]
  if(upper < 1) {
    m_upper <- max(floor(m * upper), 1)
    if(m_upper < m / n_supp) {
      m_upper <- ceiling(m / n_supp)
    }
    if(!is.null(augment)) {
      m_upper <- pmax(rep(m_upper, n_supp) - augment[supp], 0)
      prob <- rep(pmax(w_supp, prob_min), times = m_upper)
      if(sum(prob > 0) < m) {
        prob <- rep(pmax(w_supp, 1 / n), times = m_upper)
      }
      wloc <- sample(rep(1:n_supp, times = m_upper), size = m, replace = FALSE, prob = prob)
    } else {
      prob <- rep(pmax(w_supp, prob_min), each = m_upper)
      if(sum(prob > 0) < m) {
        prob <- rep(pmax(w_supp, 1 / n), each = m_upper)
      }
      wloc <- sample(rep(1:n_supp, each = m_upper), size = m, replace = FALSE, prob = prob)
    }
  } else {
    prob <- pmax(w_supp, prob_min)
    if(sum(prob > 0) < m) {
      prob <- pmax(w_supp, 1 / n)
    }
    wloc <- sample.int(n_supp, size = m, replace = isTRUE(m > n_supp), prob = prob)
  }
  w_round <- tabulate(wloc, nbins = n_supp)
  if(n_supp < n) {
    w[!supp] <- 0
  }
  w[supp] <- w_round
  return(w)
}

round_puk <- function(w, m, tol = (.Machine$double.eps)^(1/4), augment = NULL, upper = 1) {
  n <- length(w)
  if(upper < 1) {
    supp <- rep(TRUE, n)
  } else {
    supp <- (w > tol)
  }
  n_supp <- sum(supp)
  ## break ties randomly
  w_supp <- w[supp] + tol * runif(n_supp)
  w_round <- pmax(ceiling(w_supp * (m - n_supp / 2)), 1)
  if(upper < 1) {
    w_round <- pmin(w_round, max(floor(m * upper), 1))
    if(!is.null(augment)) {
      w_round <- w_round - augment
    }
  }
  if(sum(w_round) != m) {
    ## remove points
    while(sum(w_round) > m) {
      i_max <- which.max((w_round - 1 + tol) / w_supp)
      w_round[i_max] <- w_round[i_max] - 1
    }
    ## add points
    while(sum(w_round) < m) {
      w_in <- pmax(w_round, tol)
      if(upper < 1) {
        if(!is.null(augment)) {
          w_in[(w_round + augment + 1) > max(m * upper, 1)] <- Inf
        } else {
          w_in[(w_round + 1) > max(m * upper, 1)] <- Inf
        }
      }
      i_min <- which.min(w_in / w_supp)
      w_round[i_min] <- w_round[i_min] + 1
    }
  }
  w[!supp] <- 0
  w[supp] <- w_round
  return(w)
}

# round_norm <- function(w, m, X, tol = (.Machine$double.eps)^(1/4), max_iter = 100, augment = NULL, upper = 1) {
# 	## target matrix
# 	n <- length(w)
# 	if(upper < 1) {
# 		tol <- -1
# 		w <- pmax(w, 0)
# 	}
# 	if(!is.null(augment)) {
# 		w1 <- (m * w + augment) / (m + sum(augment))
# 		m1 <- m + sum(augment)
# 		supp <- (w > tol)
# 		supp1 <- (w1 > tol)
# 		X_supp <- X[supp, , drop = FALSE]
# 		X_supp1 <- X[supp1, , drop = FALSE]
# 	} else {
# 		w1 <- w
# 		m1 <- m
# 		supp <- supp1 <- (w > tol)
# 		X_supp <- X_supp1 <- X[supp, , drop = FALSE]
# 	}
# 	n_supp <- sum(supp)
# 	invert <- TRUE
# 	wsqrt <- sqrt(pmax(w1[supp1], 0) / sum(w1[supp1]) * m1)
# 	M0 <- crossprod(wsqrt * X_supp1)
# 	Minv0 <- tryCatch(chol2inv(chol(M0)), error = function(e) e)
# 	if(inherits(Minv0, "error")) {
# 		## target information instead of inverse information
# 		Minv0 <- M0
# 		invert <- FALSE
# 		warning("Information matrix is not positive definite, minimizing the norm with respect to the information instead of the covariance matrix")
# 	}
# 	## initial exact design
# 	Minv <- NULL
# 	ntrial <- 0  ## limit sample attempts
# 	if(!is.null(augment)) {
# 		wloc0 <- cumsum(supp1)[rep(1:n, times = augment)]
# 	} else {
# 		wloc0 <- NULL
# 	}
# 	while(is.null(Minv) && ntrial < 50) {
# 		loc <- setdiff(1:n_supp, wloc0)
# 		if(length(loc) < m) {
# 			loc <- 1:n_supp
# 		}
# 		wloc <- sample(loc, size = m, replace = FALSE, prob = w[supp][loc])
# 		M <- crossprod(X_supp1[c(wloc, wloc0), , drop = FALSE])
# 		if(invert) {
# 			Minv <- tryCatch(chol2inv(chol(M)), error = function(e) NULL)
# 		} else {
# 			Minv <- M
# 		}
# 		ntrial <- ntrial + 1
# 	}
# 	if(is.null(Minv)) {
# 		## return highest weights
# 		if(!is.null(augment)) {
# 			w[augment > 0] <- 0
# 		}
# 		wloc <- order(w, decreasing = TRUE)[1:m]
# 		w[] <- 0
# 		w[wloc] <- 1
# 		return(w)
# 	}
# 	## initialize
# 	w_round <- rep(0, n_supp)
# 	w_round[wloc] <- 1
# 	norm0 <- norm(Minv - Minv0, type = "F")
# 	xi <- matrix(nrow = max_iter, ncol = 2)
# 	## exchange algorithm
# 	for(iter in 1:max_iter) {
# 		## delta function: add x_i, remove x_j
# 		delta <- matrix(0, nrow = m, ncol = n_supp)
# 		for(i in 1:m) {
# 			if(invert) {
# 				Minvi <- rank1_update(Minv, X_supp[wloc[i], ], sgn = -1)
# 				if(!is.null(Minvi)) {
# 					for(j in 1:n_supp) {
# 						Minvij <- rank1_update(Minvi, X_supp[j, ], sgn = 1)
# 						if(!is.null(Minvij)) {
# 							delta[i, j] <- norm0 / norm(Minvij - Minv0, type = "F")
# 						}
# 					}
# 				}
# 			} else {
# 				Mi <- Minv - tcrossprod(X_supp[wloc[i], ])
# 				for(j in 1:n_supp) {
# 					Mij <- Mi + tcrossprod(X_supp[j, ])
# 					delta[i, j] <- norm0 / norm(Mij - Minv0, type = "F")
# 				}
# 			}
# 		}
# 		if(upper < 1) {
# 			delta[, (w_round + 1) > max(m * upper, 1)] <- 0
# 		}
# 		delta_argmax <- which.max(delta)
# 		if(delta[delta_argmax] >  1 + .Machine$double.eps^0.25) {
# 			## exchange design point
# 			mod_m <- delta_argmax %% m
# 			xout <- wloc[ifelse(mod_m == 0, m, mod_m)]
# 			xin <- (delta_argmax - 1) %/% m + 1
# 			w_round[xout] <- w_round[xout] - 1
# 			w_round[xin] <- w_round[xin] + 1
# 			if(iter > 1) {
# 				## stop early when cycling
# 				if(any(xi[1:(iter - 1), 1] == xin & xi[1:(iter - 1), 2] == xout)) {
# 					break
# 				}
# 			}
# 			## update variables
# 			wloc <- rep(1:n_supp, times = w_round)
# 			if(invert) {
# 				Minv1 <- double_rank1_update(Minv, X_supp[xin, ], X_supp[xout, ])
# 				if(is.null(Minv1)) {
# 					## revert change
# 					w_round[xout] <- w_round[xout] + 1
# 					w_round[xin] <- w_round[xin] - 1
# 					break
# 				} else {
# 					Minv <- Minv1
# 				}
# 			} else {
# 				Minv <- Minv - tcrossprod(X_supp[xout, ]) + tcrossprod(X_supp[xin, ])
# 			}
# 			norm0 <- norm(Minv - Minv0, type = "F")
# 			xi[iter, ] <- c(xin, xout)
# 		} else {
# 			break
# 		}
# 	}
# 	w[!supp] <- 0
# 	w[supp] <- w_round
# 	return(w)
# }

lower_tri_idx <- function(i, j, k) {
  if(i > j) {
    k * (j - 1) + i - (j * (j + 1)) / 2
  } else if (j > i) {
    k * (i - 1) + j - (i * (i + 1)) / 2
  } else {
    NA_integer_
  }
}

round_highs <- function(w, m, X, orthogonal, tol = .Machine$double.eps^(1/3), zeros = NULL, Xzero = NULL, colmeans = NULL,
                        augment = NULL, control = list(), seed = NULL, crit = NA, n_supp = NULL, write_mps = NULL) {

  ## initialize
  if(is.null(w)) {
    w <- rep(1 / nrow(X), times = nrow(X))
  }
  n <- length(w)
  w0 <- m * w
  p <- ncol(X)
  highs_available <- .Call("R_highs_available", NULL, PACKAGE = "cvxoptdes")
  ## reduce problem size
  if(!is.null(n_supp)) {
    u <- range(w)
    is_binary <- all(w == u[1] | w == u[2])
    if(!is_binary && length(w) > n_supp) {
      supp <- (w > quantile(w, 1 - n_supp / length(w)))
    } else {
      supp <- rep(TRUE, length(w))
    }
  } else {
    supp <- rep(TRUE, length(w))
  }
  if(!is.null(augment)) {
    supp <- supp | (augment > 0)
    m1 <- m + sum(augment)
    w0 <- w0 + augment
  } else {
    m1 <- m
  }
  n_supp <- sum(supp)
  X_supp <- X[supp, , drop = FALSE]
  stopifnot("method 'orthogonal' (currently) requires m < n, with n the total number of design points." = m1 <= n_supp)

  ## highs optimization
  if(highs_available || requireNamespace("highs", quietly = TRUE)) {
    ## MILP constraints
    col_idx <- combn(p, 2)
    col_idx <- col_idx[, !is.infinite(orthogonal[t(col_idx)])]
    b <- abs(orthogonal[t(col_idx)])
    L <- length(b)
    X0 <- X_supp
    if(any(b > sqrt(.Machine$double.eps))) {
      X0 <- center_scale(w = w0[supp], X = X_supp, center = FALSE, scale = TRUE)
    }
    X0[abs(X0) < sqrt(.Machine$double.eps)] <- 0
    X0_sp <- as(X0, "dgCMatrix")
    A0 <- t(X0_sp[, col_idx[1, ], drop = FALSE] * X0_sp[, col_idx[2, ], drop = FALSE])
    ## optimization constraints
    ## A0 == |X[ij]| <= rho[ij]
    if(is.null(Xzero)) {
      ## A[, i] * A[, i] > eps
      b1_lwr <- rep(1e-4, p)
      b1_upr <- rep(Inf, p)
      A1 <- t(X0_sp * X0_sp)
    } else {
      b1_lwr <- b1_upr <- numeric(0)
      A1 <- Matrix(0, nrow = 0, ncol = n_supp)
    }
    ## sum(w) == m
    A2 <- Matrix(rep(1, times = n_supp), nrow = 1, ncol = n_supp, sparse = TRUE)
    ## zero entries
    n_zeros <- 0
    zeromin <- zeromax <- numeric(0)
    if(!is.null(Xzero)) {
      A3 <- as(t(Xzero[supp, , drop = FALSE]), "dgCMatrix")
      if(!is.null(zeros)) {
        if(is.factor(zeros)) {
          ## equal nr zeros across columns
          has_zero <- !is.na(zeros) & (rowSums(A3 != 0) > 0)
          A3 <- A3[has_zero, , drop = FALSE]
          zeros <- zeros[has_zero]
          zero_lvls <- tapply(zeros, zeros, length)
          zero_lvls <- zero_lvls[zero_lvls > 1]
          if(length(zero_lvls)) {
            n_zeros <- sum(zero_lvls)
            ## at least one zero per column
            zeromin <- rep(1, n_zeros)
            zeromax <- rep(m1 - 1, n_zeros)
            zero_idx <- split(seq_along(zeros), f = zeros)
            A3 <- cbind(A3, Diagonal(n = n_zeros, x = -1))
            A3_eq <- Matrix(0, nrow = 0, ncol = n_zeros)
            for(idx in zero_idx) {
              A3_eq_lvl <- Matrix(0, nrow = length(idx) - 1, ncol = n_zeros)
              A3_eq_lvl[cbind(seq_along(idx[-1]), idx[-length(idx)])] <- 1
              A3_eq_lvl[cbind(seq_along(idx[-1]), idx[-1])] <- -1
              A3_eq <- rbind(A3_eq, A3_eq_lvl)
            }
            A3 <- rbind(A3, cbind(Matrix(0, nrow = nrow(A3_eq), ncol = n_supp), A3_eq))
            zero_lwr <- zero_upr <- rep(0, nrow(A3))
            ## widen other matrices
            A0 <- cbind(A0, Matrix(0, nrow = nrow(A0), ncol = n_zeros))
            A1 <- cbind(A1, Matrix(0, nrow = nrow(A1), ncol = n_zeros))
            A2 <- cbind(A2, Matrix(0, nrow = nrow(A2), ncol = n_zeros))
          } else {
            ## at least one zero per column
            has_zero <- rowSums(A3 != 0) > 0
            zero_lwr <- rep(1, sum(has_zero))
            zero_upr <- rep(m1 - 1, sum(has_zero))
            A3 <- A3[has_zero, , drop = FALSE]
          }
        } else {
          ## fixed nr zeros per column
          has_zero <- !is.na(zeros) & (rowSums(A3 != 0) > 0)
          zero_lwr <- zero_upr <- round(m1 * zeros)
          zero_lwr <- zero_lwr[has_zero]
          zero_upr <- zero_upr[has_zero]
          A3 <- A3[has_zero, , drop = FALSE]
        }
      } else {
        ## at least one zero per column
        has_zero <- rowSums(A3 != 0) > 0
        zero_lwr <- rep(1, sum(has_zero))
        zero_upr <- rep(m1 - 1, sum(has_zero))
        A3 <- A3[has_zero, , drop = FALSE]
      }
    } else {
      has_zero <- rep(FALSE, p)
      zero_lwr <- zero_upr <- numeric(0)
      A3 <- Matrix(0, nrow = 0, ncol = n_supp)
    }
    if(!is.null(colmeans) && any(!is.na(colmeans))) {
      A4 <- as(t(X_supp[, !is.na(colmeans), drop = FALSE]), "dgCMatrix")
      if(n_zeros > 0) {
        A4 <- cbind(A4, Matrix(0, nrow = nrow(A4), ncol = n_zeros))
      }
      mean_lwr <- mean_upr <- colmeans[!is.na(colmeans)]
    } else {
      mean_lwr <- mean_upr <- numeric(0)
      A4 <- Matrix(0, nrow = 0, ncol = n_supp + n_zeros)
    }

    ## linear objective
    if(any(b > sqrt(.Machine$double.eps))) {
      n_aux <- sum(b > sqrt(.Machine$double.eps))
      col_idx_pos <- col_idx[, b > sqrt(.Machine$double.eps), drop = FALSE]
      A5 <- t(X0_sp[, col_idx_pos[1, ], drop = FALSE] * X0_sp[, col_idx_pos[2, ], drop = FALSE])
      zmax <- rowSums(abs(A5))
      if(n_zeros > 0) {
        A5 <- cbind(A5, Matrix(0, nrow = nrow(A5), ncol = n_zeros))
      }
      A5_neg <- cbind(A5, Diagonal(n = n_aux, x = -1))
      A5_pos <- cbind(A5, Diagonal(n = n_aux))
      A5 <- rbind(A5_neg, A5_pos)
      rho_lwr <- c(-2 * zmax, rep(0, n_aux))
      rho_upr <- c(rep(0, n_aux), 2 * zmax)
      if(length(unique(w0[supp])) > 1) {
        L <- c(abs(1 - w0[supp]) - abs(w0[supp]), rep(0, n_zeros), rep(1, n_aux))
      } else {
        L <- rep(c(0, 0, 1), times = c(n_supp, n_zeros, n_aux))
      }
      ## widen other matrices
      A0 <- cbind(A0, Matrix(0, nrow = nrow(A0), ncol = n_aux))
      A1 <- cbind(A1, Matrix(0, nrow = nrow(A1), ncol = n_aux))
      A2 <- cbind(A2, Matrix(0, nrow = nrow(A2), ncol = n_aux))
      A3 <- cbind(A3, Matrix(0, nrow = nrow(A3), ncol = n_aux))
      A4 <- cbind(A4, Matrix(0, nrow = nrow(A4), ncol = n_aux))
    } else {
      n_aux <- 0
      A5 <- Matrix(0, nrow = 0, ncol = n_supp + n_zeros)
      rho_lwr <- numeric(0)
      rho_upr <- numeric(0)
      zmax <- numeric(0)
      L <- c(abs(1 - w0[supp]) - abs(w0[supp]), rep(0, n_zeros))
    }
    ## control parameters
    hi_control <- list(threads = 1L, time_limit = 180, log_to_console = TRUE,
                       mip_detect_symmetry = TRUE, mip_heuristic_effort = 1.0, presolve = "on",
                       mip_feasibility_tolerance = 1e-5, mip_lifting_for_probing = 3L,
                       mip_heuristic_run_zi_round = TRUE, mip_heuristic_run_shifting = TRUE)
    if(is.list(control)) {
      hi_control <- modifyList(hi_control, control)
      hi_control <- Filter(Negate(is.na), hi_control)
    }
    if(!is.null(seed)) {
      hi_control$random_seed <- as.integer(seed)
    }
    ## solve MILP
    if(highs_available) {
      ## split control parameters
      ctrl_dbl <- Filter(function(x) is.numeric(x) && !is.integer(x), hi_control)
      if(length(ctrl_dbl)) {
        ctrl_dbl <- unlist(ctrl_dbl, recursive = FALSE)
      } else {
        ctrl_dbl <- NULL
      }
      ctrl_int <- Filter(is.integer, hi_control)
      if(length(ctrl_int)) {
        ctrl_int <- unlist(ctrl_int, recursive = FALSE)
      } else {
        ctrl_int <- NULL
      }
      ctrl_bool <- Filter(is.logical, hi_control)
      if(length(ctrl_bool)) {
        ctrl_bool <- unlist(ctrl_bool, recursive = FALSE)
      } else {
        ctrl_bool <- NULL
      }
      ctrl_str <- Filter(is.character, hi_control)
      if(length(ctrl_str)) {
        ctrl_str <- unlist(ctrl_str, recursive = FALSE)
      } else {
        ctrl_str <- NULL
      }
      if(is.null(augment)) {
        hi_sol <- .Call(
          "R_highs_solve",
          L,                     ## linearized costs
          c(rep(0, n_supp), zeromin, rep(0, n_aux)),
          c(rep(1, n_supp), zeromax, zmax),   ## upr
          rbind(A0, A1, A2, A3, A4, A5),     ## A
          as.numeric(c(-b, b1_lwr, m1, zero_lwr, mean_lwr, rho_lwr)),   ## lhs
          as.numeric(c(b, b1_upr, m1, zero_upr, mean_upr, rho_upr)),    ## rhs
          rep(c(2L, 1L), times = c(n_supp, n_zeros + n_aux)),   ## types
          FALSE,                     ## maximize
          0.,                        ## offset
          ctrl_dbl,
          ctrl_int,
          ctrl_bool,
          ctrl_str,
          write_mps,                 ## mps file path
          PACKAGE = "cvxoptdes"
        )
      } else {
        A <- rbind(A0, A1, A2, A3, A4, A5)
        sol0 <- as.vector(A[, 1:n_supp, drop = FALSE] %*% augment[supp])
        hi_sol <- .Call(
          "R_highs_solve",
          L,          ## linearized costs
          c(rep(0, n_supp), zeromin, rep(0, n_aux)),
          c(pmax(rep(1., n_supp) - augment[supp], 0.), zeromax, zmax),
          A,
          as.numeric(c(-b, b1_lwr, m1, zero_lwr, mean_lwr, rho_lwr) - sol0),  ## lhs
          as.numeric(c(b, b1_upr, m1, zero_upr, mean_upr, rho_upr) - sol0),   ## rhs
          rep(c(2L, 1L), times = c(n_supp, n_aux + n_zeros)),   ## types
          FALSE,                     ## maximum
          0.,                        ## offset
          ctrl_dbl,
          ctrl_int,
          ctrl_bool,
          ctrl_str,
          write_mps,                 ## mps file path
          PACKAGE = "cvxoptdes"
        )
      }
      if(is.element(hi_sol$hi_model_status, c(7L, 11L, 12L, 13L, 14L, 16L, 18L, 19L))) {
        if(any(is.na(hi_sol$sol))) {
          warning(sprintf("HiGHS solver failed to find a feasible solution within the specified time limit."))
          return(NULL)
        } else {
          w0[supp] <- round(hi_sol$sol[1:n_supp])
        }
      } else {
        warning(sprintf("HiGHS solver failed with status code: %d", hi_sol$hi_model_status))
      }
    } else {
      ## initialize model
      if(is.null(augment)) {
        hi_mod <- highs::highs_model(
          L = L,  ## linearized costs
          lower = c(rep(0, n_supp), zeromin, rep(0, n_aux)),
          upper = c(rep(1, n_supp), zeromax, zmax),
          A = rbind(A0, A1, A2, A3, A4, A5),
          lhs = c(-b, b1_lwr, m1, zero_lwr, mean_lwr, rho_lwr),
          rhs = c(b, b1_upr, m1, zero_upr, mean_upr, rho_upr),
          types = rep(c(2L, 1L), times = c(n_supp, n_aux + n_zeros)),   ## w_i \in {0, 1}
          maximum = FALSE,
          offset = 0
        )
      } else {
        A <- rbind(A0, A1, A2, A3, A4, A5)
        sol0 <- as.vector(A[, 1:n_supp, drop = FALSE] %*% augment[supp])
        hi_mod <- highs::highs_model(
          L = L,  ## linearized costs
          lower = c(rep(0, n_supp), zeromin, rep(0, n_aux)),
          upper = c(pmax(rep(1., n_supp) - augment[supp], 0.), zeromax, zmax),
          A = A,
          lhs = c(-b, b1_lwr, m1, zero_lwr, mean_lwr, rho_lwr) - sol0,
          rhs = c(b, b1_upr, m1, zero_upr, mean_upr, rho_upr) - sol0,
          types = rep(c(2L, 1L), times = c(n_supp, n_aux + n_zeros)),   ## w_i \in {0, 1}
          maximum = FALSE,
          offset = 0
        )
      }
      hi_control <- do.call(highs::highs_control, args = hi_control)
      hi_solver <- highs::hi_new_solver(hi_mod)
      highs::hi_solver_set_options(hi_solver, hi_control)
      highs::hi_solver_run(hi_solver)
      if(!is.null(write_mps)) {
        highs::hi_solver_write_model(hi_solver, write_mps)
      }
      if(is.element(highs::hi_solver_status(hi_solver), c(7L, 11L, 12L, 13L, 14L, 16L))) {
        hi_sol <- highs::hi_solver_get_solution(hi_solver)
        w0[supp] <- round(hi_sol$col_value[1:n_supp])
        if(sum(w0[supp]) < 1) {
          warning(sprintf("HiGHS solver failed to find a feasible solution within the specified time limit."))
          return(NULL)
        }
      } else {
        warning(sprintf("HiGHS solver failed with status: %s", highs::hi_solver_status_message(hi_solver)))
      }
    }
  } else {
    warning("HiGHS solver not available, install the 'highs' package to solve for an exact feasible solution.")
  }
  w[!supp] <- 0
  w[supp] <- w0[supp]
  return(w)
}

center_scale <- function(w, X, center = TRUE, scale = TRUE) {
  X0 <- X
  if(center) {
    wmean <- colSums(w / sum(w) * X)
    X0 <- sweep(X, 2, wmean, "-")
  }
  if(scale) {
    D <- sqrt(pmax(colSums(w * X0^2), 0))
    if(any(D < .Machine$double.eps)) {
      D[D < .Machine$double.eps] <- 1   ## zero variance
    }
    X0 <- sweep(X0, 2, D, "/")
  }
  return(X0)
}


