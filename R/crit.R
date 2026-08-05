update_crit <- function(w0 = NULL, w1, X, m = length(w1), M = NULL, crit = NA, cost = NULL, lambda = 0, alpha = 0, strata = NULL, ortho = NULL) {
	if(!is.null(attr(w0, "criterion"))) {
		if(is.na(crit) || identical(crit, attr(w0, "criterion")$crit)) {
			crit0 <- attr(w0, "criterion")
		} else {
			crit0 <- list(crit = crit, crit.value = NA)
		}
		crit1 <- switch(crit0$crit,
				A = a_crit(w1, X, m, M, cost = cost, lambda = lambda, alpha = alpha),
				I = i_crit(w1, X, m, M, cost = cost, lambda = lambda, alpha = alpha),
				D = d_crit(w1, X, m, M, ln = FALSE, cost = cost, lambda = lambda, alpha = alpha),
				G = g_crit(w1, X, m, M, cost = cost, lambda = lambda, alpha = alpha),
				alias = alias_crit(w1, X, m, M, ortho, cost = cost, alpha = alpha)
		)
		attr(w1, "criterion") <- list(
				crit = crit0$crit,
				crit.value = crit1,
				crit.rel.eff = ifelse(crit1 >= 0 && isTRUE(crit0$crit.value > 0), crit1 / crit0$crit.value, NA)
		)
	} else if (!is.na(crit)) {
		crit1 <- switch(crit,
				A = a_crit(w1, X, m, M, cost = cost, lambda = lambda, alpha = alpha),
				I = i_crit(w1, X, m, M, cost = cost, lambda = lambda, alpha = alpha),
				D = d_crit(w1, X, m, M, ln = FALSE, cost = cost, lambda = lambda, alpha = alpha),
				G = g_crit(w1, X, m, M, cost = cost, lambda = lambda, alpha = alpha),
				alias = alias_crit(w1, X, m, M, ortho, cost = cost, alpha = alpha)
		)
		attr(w1, "criterion") <- list(
				crit = crit,
				crit.value = crit1,
				crit.rel.eff = NA
		)
	}
	if(!is.null(strata)) {
		attr(w1, "strata") <- tapply(w1, strata$strata, sum)
	}
	return(w1)
}

update_crit_eta <- function(w0 = NULL, w1, X, Z, eta, Zind = NULL, M = NULL, crit = NA, cost = NULL, lambda = 0, alpha = 0, strata = NULL) {
	## weighted information matrix
	## M = XW^{1/2}V^{-1}W^{1/2}X'
	if(is.null(M)) {
		M <- M_eta(w = pmax(w1, 0), X = X, Z = Z, Zind = Zind, eta = eta, scale = TRUE)
	}
	if(!is.null(attr(w0, "criterion"))) {
		if(is.na(crit) || identical(crit, attr(w0, "criterion")$crit)) {
			crit0 <- attr(w0, "criterion")
		} else {
			crit0 <- list(crit = crit, crit.value = NA)
		}
		crit1 <- switch(crit0$crit,
				A = a_crit(w = w1, X = X, M = M, cost = cost, lambda = lambda, alpha = alpha),
				I = i_crit(w = w1, X = X, M = M, cost = cost, lambda = lambda, alpha = alpha),
				D = d_crit(w = w1, X = X, M = M, cost = cost, lambda = lambda, alpha = alpha),
				G = g_crit(w = w1, X = X, M = M, cost = cost, lambda = lambda, alpha = alpha),
				alias = alias_crit(w = w1, X = X, M = M, cost = cost, alpha = alpha)
		)
		attr(w1, "criterion") <- list(
				crit = crit0$crit,
				crit.value = crit1,
				crit.rel.eff = ifelse(crit1 >= 0 && isTRUE(crit0$crit.value > 0), crit1 / crit0$crit.value, NA)
		)
	} else if (!is.na(crit)) {
		crit1 <- switch(crit,
				A = a_crit(w = w1, X = X, M = M, cost = cost, lambda = lambda, alpha = alpha),
				I = i_crit(w = w1, X = X, M = M, cost = cost, lambda = lambda, alpha = alpha),
				D = d_crit(w = w1, X = X, M = M, cost = cost, lambda = lambda, alpha = alpha),
				G = g_crit(w = w1, X = X, M = M, cost = cost, lambda = lambda, alpha = alpha),
				alias = alias_crit(w = w1, X = X, M = M, cost = cost, alpha = alpha)
		)
		attr(w1, "criterion") <- list(
				crit = crit,
				crit.value = crit1,
				crit.rel.eff = NA
		)
	}
	if(!is.null(strata)) {
		attr(w1, "strata") <- tapply(w1, strata$strata, sum)
	}
	return(w1)
}

d_crit <- function(w, X, m = length(w), M = NULL, ln = FALSE, cost = NULL, lambda = 0, alpha = 0) {
	if(is.null(M)) {
		if(m < length(w)) {
			X <- X[w > 0, , drop = FALSE]
			w1 <- w[w > 0]
		} else {
		  w1 <- w
		}
		wsqrt <- sqrt(pmax(w1, 0) / sum(w1))
		M <- crossprod(wsqrt * X)
	}
	if(lambda > 0) {
		M <- M + diag(lambda, nrow = nrow(M))
	}
	crit <- .Call(crit_impl, M, NULL, NULL, 2L, FALSE, PACKAGE = "cvxoptdes")
	if(!is.null(cost) && alpha > 0) {
		crit <- crit - alpha * sum(cost * w) / sum(w)
	}
	if(ln) {
		return(crit)
	} else {
		return(exp(crit))
	}
}

a_crit <- function(w, X, m = length(w), M = NULL, cost = NULL, lambda = 0, alpha = 0) {
	if(is.null(M)) {
		if(m < length(w)) {
			X <- X[w > 0, , drop = FALSE]
			w1 <- w[w > 0]
		} else {
		  w1 <- w
		}
		wsqrt <- sqrt(pmax(w1, 0) / sum(w1))
		M <- crossprod(wsqrt * X)
	}
	if(lambda > 0) {
		M <- M + diag(lambda, nrow = nrow(M))
	}
	crit <- .Call(crit_impl, M, NULL, NULL, 0L, TRUE, PACKAGE = "cvxoptdes")
	if(!is.null(cost) && alpha > 0) {
		crit <- crit - alpha * sum(cost * w) / sum(w)
	}
	return(crit)
}

i_crit <- function(w, X, m = length(w), M = NULL, cost = NULL, lambda = 0, alpha = 0) {
	if(is.null(M)) {
		if(m < length(w)) {
			X1 <- X[w > 0, , drop = FALSE]
			w1 <- w[w > 0]
		} else {
			X1 <- X
			w1 <- w
		}
		wsqrt <- sqrt(pmax(w1, 0) / sum(w1))
		M <- crossprod(wsqrt * X1)
	}
	if(lambda > 0) {
		M <- M + diag(lambda, nrow = nrow(M))
	}
	XtX <- crossprod(X)
	crit <- .Call(crit_impl, M, X, XtX, 1L, TRUE, PACKAGE = "cvxoptdes")
	if(!is.null(cost) && alpha > 0) {
		crit <- crit - alpha * sum(cost * w) / sum(w)
	}
	return(crit)
}

g_crit <- function(w, X, m = length(w), M = NULL, cost = NULL, lambda = 0, alpha = 0) {
	if(is.null(M)) {
		if(m < length(w)) {
			X1 <- X[w > 0, , drop = FALSE]
			w1 <- w[w > 0]
		} else {
			X1 <- X
			w1 <- w
		}
		wsqrt <- sqrt(pmax(w1, 0) / sum(w1))
		M <- crossprod(wsqrt * X1)
	}
	if(lambda > 0) {
		M <- M + diag(lambda, nrow = nrow(M))
	}
	crit <- .Call(crit_impl, M, X, NULL, 3L, FALSE, PACKAGE = "cvxoptdes")
	if(!is.null(cost) && alpha > 0) {
		crit <- crit - alpha * sum(cost * w) / sum(w)
	}
	return(crit)
}

alias_crit <- function(w, X, m = length(w), M = NULL, ortho = NULL, cost = NULL, alpha = 0) {
	if(is.null(M)) {
		if(m < length(w)) {
			X1 <- X[w > 0, , drop = FALSE]
			w1 <- w[w > 0]
		} else {
			X1 <- X
			w1 <- w
		}
		wsqrt <- sqrt(pmax(w1, 0) / sum(w1))
		M <- crossprod(wsqrt * X1)
	}
	crit <- .Call(crit_impl, M, NULL, ortho, 4L, TRUE, PACKAGE = "cvxoptdes")
	if(!is.null(cost) && alpha > 0) {
		crit <- crit - alpha * sum(cost * w) / sum(w)
	}
	return(crit)
}

M_eta <- function(w, X, Z, Zind = NULL, eta, scale = FALSE, lambda = 0) {
	if(!is.matrix(Z)) {
		Z <- as.matrix(Z)
	}
	if(is.null(Zind)) {
		Zind <- c((Z > 0) %*% seq_len(ncol(Z)))
	}
	if(!missing(w)) {
		if(scale) {
			w_scl <- w / sum(w)
		} else {
			w_scl <- w
			lambda <- lambda * sum(w)
		}
		if(is.infinite(eta)) {
			zeta <- 1 /  colSums(w_scl * Z)
		} else {
			zeta <- 1 / (1 / eta + colSums(w_scl * Z))
		}
		Z <- sqrt(w_scl) * Z
		X <- sqrt(w_scl) * X
	} else {
		if(is.infinite(eta)) {
			zeta <- 1 / colSums(Z)
		} else {
			zeta <- 1 / (1 / eta + colSums(Z))
		}
	}
	A <- diag(0, nrow = ncol(X))
	for(i in 1:ncol(Z)) {
		ids <- (Zind == i)
		if(any(ids)) {
			ni <- sum(ids)
			Xi <- X[ids, , drop = FALSE]
			Vi <- diag(ni) - tcrossprod(sqrt(zeta[i]) * Z[ids, i])
			A <- A + crossprod(Xi, Vi %*% Xi)
		}
	}
	if(lambda > 0) {
	  A <- A + diag(lambda, nrow = nrow(A))
	}
	return(A)
}
