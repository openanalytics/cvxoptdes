dist_exchange <- function(w, X, strata = NULL, max_iter = 100, tol = .Machine$double.eps^(1/3),
                          type = c("dist-sum", "dist-min"), augment = NULL, orthogonal = NULL, upper = 1) {
  m <- sum(w)
  n <- nrow(X)
  wmax <- w
  blocks <- NULL
  if(!is.null(strata)) {
    blocks <- split(seq_len(n), f = strata$strata)
    sizes <- tapply(w, strata$strata, sum) - strata$proportions
    if(any(is.na(sizes))) {
      sizes[is.na(sizes)] <- 0
    }
  }
  wloc <- rep(1:n, times = w)
  wloc0 <- NULL
  if(!is.null(augment)) {
    wloc0 <- rep(1:n, times = augment)
    m <- length(wloc) + length(wloc0)
  }
  X1 <- X[c(wloc, wloc0), , drop = FALSE]
  xi <- matrix(nrow = max_iter, ncol = 2)
  ## pairwise distances
  alldist <- as.matrix(dist(X1, method = "euclidean"))
  newdist <- .Call(dist_impl, X1, X, as.integer(c(wloc, wloc0)), PACKAGE = "cvxoptdes")
  for(iter in 1:max_iter) {
    if(type == "dist-sum") {
      ## maximize sum of pairwise distances
      ## d(xi, xj) = -sum(alldist[i, ]) + sum(newdist[-i, j])
      delta <- t(array(colSums(newdist), dim = rev(dim(newdist))))
      delta <- delta - colSums(alldist) - newdist
    } else {
      ## maximize minimum pairwise distance
      ## minimum pairwise distance after removing x_i
      alldist[row(alldist) == col(alldist)] <- Inf
      allmin <- apply(alldist, 1, min)
      loomin <- .Call(delta_min_impl, as.matrix(allmin, ncol = 1), tol, PACKAGE = "cvxoptdes")
      delta1 <- array(loomin, dim = dim(newdist))
      delta1[cbind(seq_along(c(wloc, wloc0)), c(wloc, wloc0))] <- 0
      ## minimum pairwise distance after adding x_j
      delta2 <- .Call(delta_min_impl, newdist, 0.0, PACKAGE = "cvxoptdes")
      delta <- pmin(delta1, delta2) - min(allmin)
    }
    if(!is.null(augment)) {
      delta[(length(wloc) + 1):m, ] <- -Inf
    }
    if(upper < 1) {
      if(!is.null(augment)) {
        delta[, (w + augment + 1) > max(m * upper, 1)] <- -Inf
      } else {
        delta[, (w + 1) > max(m * upper, 1)] <- -Inf
      }
    }
    if(!is.null(blocks)) {
      delta1 <- array(-Inf, dim(delta))
      if(any(sizes != 0)) {
        ## switch between incorrect blocks
        inblocks <- unlist(blocks[sizes < 0], recursive = FALSE, use.names = FALSE)
        for(i in which(sizes > 0)) {
          indices <- as.matrix(expand.grid(na.omit(match(blocks[[i]], wloc)), inblocks))
          delta1[indices] <- delta[indices]
        }
      } else {
        ## switch within blocks only
        for(i in seq_along(blocks)) {
          indices <- as.matrix(expand.grid(na.omit(match(blocks[[i]], wloc)), blocks[[i]]))
          delta1[indices] <- delta[indices]
        }
      }
      delta <- delta1
    }
    delta_argmax <- which.max(delta)
    if(delta[delta_argmax] > tol || (!is.null(blocks) && any(sizes != 0))) {
      mod_m <- delta_argmax %% m
      xout_loc <- ifelse(mod_m == 0, m, mod_m)
      xout <- wloc[xout_loc]
      xin <- (delta_argmax - 1) %/% m + 1
      w[xout] <- w[xout] - 1
      w[xin] <- w[xin] + 1

      if(iter > 1 && (is.null(blocks) || all(sizes == 0))) {
        ## stop early when cycling
        if(any(xi[1:(iter - 1), 1] == xin & xi[1:(iter - 1), 2] == xout)) {
          break
        }
      }
      xi[iter, ] <- c(xin, xout)
      if(!is.null(blocks) && any(sizes != 0)) {
        blockout <- as.integer(strata$strata[xout])
        blockin <- as.integer(strata$strata[xin])
        sizes[blockout] <- sizes[blockout] - 1
        sizes[blockin] <- sizes[blockin] + 1
      }
      ## update distance matrices
      wloc <- rep(1:n, times = w)
      X1 <- X[c(wloc, wloc0), , drop = FALSE]
      xin_dist <- .Call(dist_impl, X[xin, , drop = FALSE], X1, -1L, PACKAGE = "cvxoptdes")
      xin_loc <- which(wloc == xin)[1] - 1L
      alldist <- alldist[-xout_loc, -xout_loc]
      alldist <- append_impl(alldist, xin_dist[-(xin_loc + 1)], after = xin_loc, margin = 1L)
      alldist <- append_impl(alldist, xin_dist, after = xin_loc, margin = 2L)
      xin_newdist <- .Call(dist_impl, X[xin, , drop = FALSE], X, as.integer(xin), PACKAGE = "cvxoptdes")
      newdist <- newdist[-xout_loc, ]
      newdist <- append_impl(newdist, xin_newdist, after = xin_loc, margin = 1L)
    } else {
      break
    }
  }
  if(!is.null(augment)) {
    w <- w + augment
  }
  ## update criterion
  alldist <- as.matrix(dist(X1, method = "euclidean"))
  if(type == "dist-sum") {
    attr(w, "criterion") <- list(
      crit = type,
      crit.value = sum(alldist[lower.tri(alldist, diag = FALSE)]),
      crit.rel.eff = NA
    )
  } else {
    attr(w, "criterion") <- list(
      crit = type,
      crit.value = min(alldist[lower.tri(alldist, diag = FALSE)]),
      crit.rel.eff = NA
    )
  }
  if(!is.null(strata)) {
    attr(w, "strata") <- tapply(w, strata$strata, sum)
  }
  return(w)
}

exchange_impl <- function(w, X, crit = c("I", "A", "D", "G", "alias"), supp = rep(TRUE, nrow(X)), strata = NULL, max_iter = 100, tol = .Machine$double.eps^(1/3),
                          augment = NULL, cost = NULL, lambda = 0, alpha = 1, upper = 1, ortho = NULL, add_crit = TRUE, eta = NA, Zt = NULL) {
  w1 <- w0 <- w
  n <- n_supp <- nrow(X)
  m <- m1 <- sum(w)
  p <- ncol(X)
  X_supp <- X
  cost_supp <- cost
  critval <- critval_m <- NA
  Minv <- blocks <- sizes <- NULL
  if(!is.null(strata)) {
    blocks <- split(seq_len(n), f = strata$strata)
    sizes <- tapply(w, strata$strata, sum) - strata$proportions
    if(any(is.na(sizes))) {
      sizes[is.na(sizes)] <- 0
    }
  } else if(is.null(augment) && (is.na(eta) || eta <= 0)) {
    w <- w1 <- w[supp]
    X_supp <- X[supp, , drop = FALSE]
    n_supp <- nrow(X_supp)
    if(!is.null(cost)) {
      cost_supp <- cost[supp]
    }
  }
  ## initial matrices
  wloc <- wloc1 <- rep(1:n_supp, times = w)
  XtX <- Z <- Zind <- NULL
  if(!is.null(augment)) {
    w1 <- w + augment
    wloc1 <- rep(1:n, times = w1)
    m1 <- sum(w1)
  }
  if(crit == "I") {
    XtX <- crossprod(X_supp)
  }
  if(is.na(eta) || eta <= 0) {
    M <- crossprod(X_supp[wloc1, , drop = FALSE])
    if(lambda > 0) {
      M <- M + diag(m1 * lambda, nrow = nrow(M))
    }
    if(crit != "alias") {
      Minv <- tryCatch(chol2inv(chol(M)), error = function(e) e)
      if(inherits(Minv, "error")) {
        if(n_supp < n) {
          w1 <- w0
        }
        if(add_crit) {
          attr(w1, "criterion") <- list(
            crit = crit,
            crit.value = 0,
            crit.rel.eff = 0
          )
        }
        return(w1)
      }
    }
    if(!is.null(cost)) {
      ## criterion without cost penalty
      critval <- switch(crit,
                        "I" = i_crit(X = X_supp, M = M / m1),
                        "A" = a_crit(M = M / m1),
                        "D" = d_crit(M = M / m1, ln = TRUE),
                        "G" = g_crit(X = X_supp, M = M / m1),
                        "alias" = alias_crit(M = M / m1)
      )
    }
  } else {
    Z <- as.matrix(t(Zt))
    Zind <- c(Z %*% seq_len(ncol(Z)))
    ## unweighted information matrix
    M <- M_eta(X = X[wloc1, , drop = FALSE], Z = Z[wloc1, , drop = FALSE], Zind = Zind[wloc1], eta = eta, lambda = lambda)
    critval_m <- switch(crit,
                        "I" = i_crit(X = X, M = M),
                        "A" = a_crit(M = M),
                        "D" = d_crit(M = M, ln = TRUE),
                        "G" = g_crit(X = X, M = M),
                        "alias" = alias_crit(M = M)
    )
    if(!is.null(cost)) {
      ## weighted information matrix
      M_scl <- M_eta(w = w1, X = X, Z = Z, Zind = Zind, eta = eta, scale = TRUE, lambda = lambda)
      critval <- switch(crit,
                        "I" = i_crit(X = X, M = M_scl),
                        "A" = a_crit(M = M_scl),
                        "D" = d_crit(M = M_scl, ln = TRUE),
                        "G" = g_crit(X = X, M = M_scl),
                        "alias" = alias_crit(M = M_scl)
      )
    }
  }

  ## optimize w/ limited support
  result <- exchange_impl_loop(w, wloc, m, X_supp, M, Minv, critval, critval_m, XtX = XtX, crit = crit, max_iter = max_iter, strata = strata,
                              blocks = blocks, sizes = sizes, tol = tol, augment = augment, cost = cost_supp, lambda = lambda, alpha = alpha,
                              upper = upper, ortho = ortho, eta = eta, Z = Z, Zind = Zind)

  for(nm in names(result)) {
    assign(nm, result[[nm]])
  }
  if(n_supp < n) {
    w <- numeric(n)
    w[supp] <- result$w
    wloc <- rep(1:n, times = w)
    if(crit == "I") {
      XtX <- crossprod(X)
    }
    M <- crossprod(X[wloc, , drop = FALSE])
    if(lambda > 0) {
      M <- M + diag(m1 * lambda, nrow = nrow(M))
    }
    if(crit != "alias") {
      Minv <- tryCatch(chol2inv(chol(M)), error = function(e) e)
      if(inherits(Minv, "error")) {
        if(add_crit) {
          attr(w, "criterion") <- list(
            crit = crit,
            crit.value = 0,
            crit.rel.eff = 0
          )
        }
        return(w)
      }
    }
    ## optimize w/ full support
    result <- exchange_impl_loop(w, wloc, m, X, M, Minv, result$critval, result$critval_m, XtX = XtX, crit = crit, max_iter = max_iter,
                              tol = tol, cost = cost, alpha = alpha, upper = upper, ortho = ortho)
    w <- result$w
  }
  if(!is.null(augment)) {
    w <- w + augment
    if(add_crit) {
      attr(w, "criterion") <- NULL
    }
  }
  ## update criterion
  if(add_crit) {
    if(is.na(eta) || eta <= 0) {
        w <- update_crit(w0 = w0, w1 = w, X = X, crit = crit, cost = cost, lambda = lambda, alpha = alpha, strata = strata, ortho = ortho)
    } else {
        w <- update_crit_eta(w0 = w0, w1 = w, X = X, Z = Z, eta = eta, Zind = Zind, M = NULL, crit = crit, cost = cost, lambda = lambda, alpha = alpha, strata = strata)
    }
  }
  return(w)
}

exchange_impl_loop <- function(w, wloc, m, X, M, Minv, critval, critval_m, XtX = NULL, crit = c("I", "A", "D", "G", "alias"), max_iter = 100, strata = NULL,
                              blocks = NULL, sizes = NULL, tol = .Machine$double.eps^(1/3), augment = NULL, cost = NULL, lambda = 0, alpha = 1, upper = 1,
                              ortho = NULL, eta = NA, Z = NULL, Zind = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  m1 <- sum(w)
  xi <- matrix(nrow = max_iter, ncol = 2)
  for(iter in 1:max_iter) {
    ## delta function
    if(is.na(eta) || eta <= 0) {
      if(crit == "D") {
        delta <- .Call(d_delta, as.integer(wloc), X, Minv, PACKAGE = "cvxoptdes")
      } else if(crit == "A") {
        delta <- .Call(a_delta, as.integer(wloc), X, Minv, PACKAGE = "cvxoptdes")
      } else if(crit == "I") {
        delta <- .Call(i_delta, as.integer(wloc), X, XtX, M, Minv, PACKAGE = "cvxoptdes")
      } else if(crit == "G") {
        delta <- .Call(g_delta, as.integer(wloc), X, M, Minv, PACKAGE = "cvxoptdes")
      } else if(crit== "alias") {
        delta <- .Call(alias_delta, as.integer(wloc), X, M, ortho, PACKAGE = "cvxoptdes")
      }
    } else {
      ## unweighted criterion
      if(crit == "D") {
        crit0 <- max(critval_m * p, log(.Machine$double.eps^0.5))
      } else if(crit == "A") {
        crit0 <- 1 / pmax(critval_m, .Machine$double.eps^0.5)
      } else if(crit == "I") {
        crit0 <- n / pmax(critval_m, .Machine$double.eps^0.5)
      } else if(crit == "G") {
        crit0 <- critval_m
      } else if(crit == "alias") {
        crit0 <- 1 - critval_m
      }
      delta <- exchange_delta_eta(wloc, X, Z, Zind, eta, crit0, crit, XtX)
    }
    delta0 <- delta
    if(!is.null(cost)) {
      ## Phi(M1) - Phi(M0)
      if(crit == "D") {
        delta_c <- log1p(pmax(delta, -1 + .Machine$double.eps^0.5)) / p
      } else if(crit == "A") {
        delta_c <- (1 / (1 / critval - m1 * delta) - critval)
      } else if(crit == "I") {
        delta_c <- (1 / (1 / critval - m1 / n * delta) - critval)
      } else if(crit == "G") {
        delta_c <- delta / m1
      } else if(crit == "alias") {
        delta_c <- delta
      }
      ## (Phi(M1) - cost(w1)) - (Phi(M0) - cost(w0))
      delta <- delta_c - cost_delta(wloc, cost, alpha / m)
    }
    if(upper < 1) {
      if(!is.null(augment)) {
        delta[, (w + augment + 1) > max(m * upper, 1)] <- -Inf
      } else {
        delta[, (w + 1) > max(m * upper, 1)] <- -Inf
      }
    }
    if(!is.null(blocks) && is.matrix(delta)) {
      delta1 <- array(-Inf, dim(delta))
      if(any(sizes != 0)) {
        ## switch between incorrect blocks
        inblocks <- unlist(blocks[sizes < 0], recursive = FALSE, use.names = FALSE)
        for(i in which(sizes > 0)) {
          indices <- as.matrix(expand.grid(na.omit(match(blocks[[i]], wloc)), inblocks))
          delta1[indices] <- delta[indices]
        }
      } else {
        ## switch within blocks only
        for(i in seq_along(blocks)) {
          indices <- as.matrix(expand.grid(na.omit(match(blocks[[i]], wloc)), blocks[[i]]))
          delta1[indices] <- delta[indices]
        }
      }
      delta <- delta1
    }
    delta_max <- max(delta, na.rm = TRUE)
    ## exchange design point
    if(delta_max > tol || (!is.null(blocks) && any(sizes != 0))) {
      delta_argmax <- which.max(delta)
      mod_m <- delta_argmax %% m
      xout <- wloc[ifelse(mod_m == 0, m, mod_m)]
      xin <- (delta_argmax - 1) %/% m + 1
      w[xout] <- w[xout] - 1
      w[xin] <- w[xin] + 1
      if(!is.null(blocks) && any(sizes != 0)) {
        blockout <- as.integer(strata$strata[xout])
        blockin <- as.integer(strata$strata[xin])
        sizes[blockout] <- sizes[blockout] - 1
        sizes[blockin] <- sizes[blockin] + 1
      }
      if(iter > 1 && (is.null(blocks) || all(sizes == 0))) {
        ## stop early when cycling
        if(any(xi[1:(iter - 1), 1] == xin & xi[1:(iter - 1), 2] == xout)) {
          break
        }
      }
      ## update variables
      wloc <- rep(1:n, times = w)
      if(is.na(eta) || eta <= 0) {
        M <- M - tcrossprod(X[xout, ]) + tcrossprod(X[xin, ])
        if(crit != "alias") {
          Minv <- double_rank1_update(Minv = Minv, u = X[xin, ], v = X[xout, ])
          if(is.null(Minv)) {
            ## revert last update
            w[xout] <- w[xout] + 1
            w[xin] <- w[xin] - 1
            break
          }
        }
        if(!is.null(cost)) {
          if(crit == "D") {
            critval <- critval + log1p(delta0[delta_argmax]) / p
          } else if(crit == "A") {
            critval <- 1 / (1 / critval - m1 * delta0[delta_argmax])
          } else if(crit == "I") {
            critval <- 1 / (1 / critval - m1 / n * delta0[delta_argmax])
          } else if(crit == "G"){
            critval <- critval + delta0[delta_argmax] / m1
          } else if(crit == "alias") {
            critval <- critval + delta0[delta_argmax]
          }
        }
      } else if(eta > 0) {
        if(!is.null(augment)) {
          w1 <- w + augment
          wloc1 <- rep(1:n, times = w1)
        } else {
          w1 <- w
          wloc1 <- wloc
        }
        ## unweighted information matrix
        M <- M_eta(X = X[wloc1, , drop = FALSE], Z = Z[wloc1, , drop = FALSE], Zind = Zind[wloc1], eta = eta, lambda = lambda)
        critval_m <- switch(crit,
                            "I" = i_crit(X = X, M = M),
                            "A" = a_crit(M = M),
                            "D" = d_crit(M = M, ln = TRUE),
                            "G" = g_crit(X = X, M = M),
                            "alias" = alias_crit(M = M)
        )
        if(!is.null(cost)) {
          ## weighted information matrix
          M_scl <- M_eta(w = w1, X = X, Z = Z, Zind = Zind, eta = eta, scale = TRUE, lambda = lambda)
          critval <- switch(crit,
                            "I" = i_crit(X = X, M = M_scl),
                            "A" = a_crit(M = M_scl),
                            "D" = d_crit(M = M_scl, ln = TRUE),
                            "G" = g_crit(X = X, M = M_scl),
                            "alias" = alias_crit(M = M_scl)
          )
        }
      }
      xi[iter, ] <- c(xin, xout)
    } else {
      break
    }
  }
  return(list(w = w, critval = critval, critval_m = critval_m))
}

double_rank1_update <- function(Minv, u, v) {
  Minv <- .Call(rank1_update, Minv, u, sgn = 1.0, PACKAGE = "cvxoptdes")
  if(!is.null(Minv)) {
    Minv <- .Call(rank1_update, Minv, v, sgn = -1.0, PACKAGE = "cvxoptdes")
  }
  return(Minv)
}

cost_delta <- function(wloc, cost, alpha = 1) {
  ## delta function: add xi, remove xj
  ## delta(xi, xj) = alpha * (c(xi) - c(xj))
  delta <- matrix(alpha * cost, byrow = TRUE, nrow = length(wloc), ncol = length(cost))
  delta <- delta - (alpha * cost[wloc])
  return(delta)
}

exchange_delta_eta <- function(wloc, X, Z, Zind, eta, crit0, crit, XtX = NULL) {

  ## initialize
  Xloc <- X[wloc, , drop = FALSE]
  Zloc <- Z[wloc, , drop = FALSE]
  Zind1 <- Zind[wloc]
  zseq <- apply(Zloc, 2, cumsum)[cbind(1:nrow(Zloc), Zind1)]
  ## group weights
  ni <- colSums(Zloc)
  ndelta <- matrix(c(0, 1, -1), nrow = 3L, ncol = length(ni))
  if(is.infinite(eta)) {
    zeta <- t(1 / (ni + t(ndelta)))
  } else {
    zeta <- t(1 / (1 / eta + ni + t(ndelta)))
  }
  ## group entities
  Vli <- vector(mode = "list", length = 3 * ncol(Z))
  Xi <- vector(mode = "list", length = ncol(Z))
  Ai <- vector(mode = "list", length = ncol(Z))
  A <- diag(0, nrow = ncol(X))
  for(i in 1:ncol(Z)) {
    ids <- (Zind1 == i)
    if(any(ids)) {
      Xi[[i]] <- Xloc[ids, , drop = FALSE]
      for(l in 1:3) {
        nli <- ni[i] + ndelta[l, i]
        if(nli > 0) {
          Vli[[(i - 1) * 3 + l]] <- diag(nli) - matrix(zeta[l, i], nrow = nli, ncol = nli)
        }
      }
      Ai[[i]] <- crossprod(Xi[[i]], Vli[[(i - 1) * 3 + 1]] %*% Xi[[i]])
      A <- A + Ai[[i]]
    } else {
      Vli[[(i - 1) * 3 + 2]] <- diag(1L) - matrix(zeta[2, i], nrow = 1L, ncol = 1L)
    }
  }
  ## initial delta matrix
  crit1 <- match(crit, c("A", "I", "D", "G", "alias")) - 1L
  delta <- .Call(delta_eta, as.integer(wloc), X, A, as.integer(zseq), crit1, crit0, as.integer(Zind), Vli, Xi, Ai, PACKAGE = "cvxoptdes")
  ## update delta matrix
  if (crit == "D") {
    delta <- exp(delta - crit0) - 1
  } else if(crit == "I") {
    delta <- crit0 - delta
  } else if(crit == "G") {
    delta <- delta - crit0
  } else {
    delta <- crit0 - delta
  }
  return(delta)
}

append_impl <- function(x, values, after = nrow(x), margin = 1L) {
  if(margin != 1L) {
    x <- t(x)
  }
  lengx <- nrow(x)
  if (after < 1) {
    x1 <- rbind(values, x, deparse.level = 0)
  } else if(after >= lengx) {
    x1 <- rbind(x, values, deparse.level = 0)
  } else {
    x1 <- rbind(x[1L:after, , drop = FALSE], values, x[(after + 1L):lengx, , drop = FALSE], deparse.level = 0)
  }
  if(margin != 1L) {
    x1 <- t(x1)
  }
  return(x1)
}

# d_delta <- function(wloc, X, Minv) {
# 	if(is.null(Minv)) {
# 		return(matrix(-Inf, nrow = length(wloc), ncol = nrow(X)))
# 	}
# 	## variance terms
# 	uMu <- rowSums(X * (X %*% Minv))
# 	vMv <- uMu[wloc]
# 	## covariance terms
# 	vMu <- (X[wloc, , drop = FALSE] %*% Minv) %*% t(X)
# 	## delta function: add xi, remove xj
# 	## delta(xi, xj) = d(xi) - d(xj) - d(xi)d(xj) + d(xj, xi)^2
# 	delta <- t(uMu - t(vMv + vMv %o% uMu)) + vMu^2
# 	return(delta)
# }

# a_delta <- function(wloc, X, Minv) {
# 	if(is.null(Minv)) {
# 		return(matrix(-Inf, nrow = length(wloc), ncol = nrow(X)))
# 	}
# 	## inverses
# 	Minv2 <- crossprod(Minv)
# 	## variance terms
# 	uMu <- rowSums(X * (X %*% Minv))
# 	vMv <- uMu[wloc]
# 	if(any(vMv == 1)) {
# 		vMv[vMv == 1] <- 1 - .Machine$double.eps^0.5
# 	}
# 	uM2u <- rowSums(X * (X %*% Minv2))
# 	vM2v <- uM2u[wloc]
# 	## covariance terms
# 	vMu <- (X[wloc, , drop = FALSE] %*% Minv) %*% t(X)
# 	vM2u <- (X[wloc, , drop = FALSE] %*% Minv2) %*% t(X)
# 	## delta function: add xi, remove xj
# 	## u2 = xj.M2.xj
# 	## v2 = xi.M2.xi
# 	## uv = xi.M2.xj
# 	## s = xi.M1.xi + (xj.M1.xi)^2 / (1 - xj.M1.xj)
# 	## w2 = v2 + 2 * xj.M1.xi / (1 - xj.M1.xj) * (uv) + (xj.M1.xi)^2 / (1 - xj.M1.xj)^2 * u2
# 	## delta(xi, xj) = u2 / (1 - xj.M1.xj) - w^2 / (1 + s)
# 	smat <- matrix(uMu, nrow = nrow(vMu), ncol = length(uMu), byrow = TRUE) + vMu^2 / (1 - vMv)
# 	wnorm <- matrix(uM2u, nrow = nrow(vMu), ncol = length(uM2u), byrow = TRUE)
# 	wnorm <- wnorm + 2 * vMu / (1 - vMv) * vM2u + vMu^2 / (1 - vMv)^2 * vM2v
# 	delta <- -(vM2v / (1 - vMv) - wnorm / (1 + smat))
# 	return(delta)
# }

# i_delta <- function(wloc, X, XtX, M, Minv) {
# 	uMinv <- solve(M, t(X))
# 	if(is.null(uMinv) || is.null(Minv)) {
# 		return(matrix(-Inf, nrow = length(wloc), ncol = nrow(X)))
# 	}
# 	## inverses
# 	vMinv <- uMinv[, wloc, drop = FALSE]
# 	## variance terms
# 	uMu <- rowSums(X * (X %*% Minv))
# 	vMv <- uMu[wloc]
# 	if(any(vMv == 1)) {
# 		vMv[vMv == 1] <- 1 - .Machine$double.eps^0.5
# 	}
# 	uHu <- colSums(uMinv * (XtX %*% uMinv))
# 	vHv <- uHu[wloc]
# 	## covariance terms
# 	vMu <- (X[wloc, , drop = FALSE] %*% Minv) %*% t(X)
# 	vHu <-  t(vMinv) %*% XtX %*% uMinv
# 	## delta function: add xi, remove xj
# 	## uHu = xj.M2.H.M2.xj
# 	## vHv = xi.M2.A.M2.xi
# 	## uHv = xi.M2.A.M2.xj
# 	## s = xi.M1.xi + (xj.M1.xi)^2 / (1 - xj.M1.xj)
# 	## wHw = vHv + 2 * xj.M1.xi / (1 - xj.M1.xj) * uHv + (xj.M1.xi)^2 / (1 - xj.M1.xj)^2 * uHu
# 	## delta(xi, xj) = u.H.u / (1 - xj.M1.xj) - w.H.w / (1 + s)
# 	smat <- matrix(uMu, nrow = nrow(vMu), ncol = length(uMu), byrow = TRUE) + vMu^2 / (1 - vMv)
# 	wnorm <- matrix(uHu, nrow = nrow(vMu), ncol = length(uHu), byrow = TRUE)
# 	wnorm <- wnorm + 2 * vMu / (1 - vMv) * vHu + vMu^2 / (1 - vMv)^2 * vHv
# 	delta <- -(vHv / (1 - vMv) - wnorm / (1 + smat))
# 	return(delta)
# }

# g_delta <- function(wloc, X, M) {
# 	uMinv <- t(solve(M, t(X)))
# 	if(is.null(uMinv)) {
# 		return(matrix(-Inf, nrow = length(wloc), ncol = nrow(X)))
# 	}
# 	n <- nrow(X)
# 	m <- length(wloc)
# 	## covariance terms
# 	zMz <- tcrossprod(X, uMinv)
# 	uMu <- diag(zMz)
# 	vMv <- uMu[wloc]
# 	if(any(vMv == 1)) {
# 		vMv[vMv == 1] <- 1 - .Machine$double.eps^0.5
# 	}
# 	vMu <- zMz[wloc, , drop = FALSE]
# 	## delta function: add xi, remove xj
# 	## s = xi.M1.xi + (xj.M1.xi)^2 / (1 - xj.M1.xj)
# 	## z.w =  xi.M1.z + xj.M1.xi / (1 - xj.M1.xj) * xj.M1.z
# 	## delta(xi, xj) = z.M1.z + xj.M1.z / (1 - xj.M1.xj) - (z.w)^2 / (1 + s)
# 	smat <- matrix(uMu, nrow = m, ncol = n, byrow = TRUE) + vMu^2 / (1 - vMv) ## (xj by xi)
# 	zw2 <- vMu / (1 - vMv)  ## (xj by xi)
# 	delta <- matrix(nrow = m, ncol = n)
# 	## loop over xi
# 	for(i in 1:n) {
# 		zw <- matrix(zMz[i, ], nrow = m, ncol = n, byrow = TRUE) + zw2[, i] * vMu  ## (xj by z)
# 		wnorm <- zw^2 / (1 + smat[, i])   ## (xj by z)
# 		delta[, i] <- 1 / apply(smat - wnorm, 1, max)
# 	}
# 	return(delta - 1 / max(uMu))
# }

# alias_delta <- function(wloc, X, M) {
# 	## initial correlations
# 	p <- ncol(X)
# 	n <- nrow(X)
# 	D <- sqrt(pmax(diag(M), .Machine$double.eps))
# 	R <- M / outer(D, D)
# 	corr0 <- mean(abs(R[lower.tri(R)]))
# 	X0j <- apply(X, 1, function(x) tcrossprod(x))
# 	dim(X0j) <- c(p, p, n)
# 	## delta function: add xi, remove xj
# 	delta <- matrix(0, nrow = length(wloc), ncol = n)
# 	for(i in seq_along(wloc)) {
# 		for(j in 1:n) {
# 			M1 <- M - X0j[, , wloc[i]] + X0j[, , j]
# 			if(all(diag(M1) > .Machine$double.eps)) {
# 				D <- sqrt(diag(M1))
# 				R1 <- M1 / outer(D, D)
# 				corr1 <- mean(abs(R1[lower.tri(R1)]))
# 				delta[i, j] <- corr0 / corr1 - 1
# 			}
# 		}
# 	}
# 	return(delta)
# }
