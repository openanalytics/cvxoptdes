require(cvxoptdes)

start <- proc.time()
cat("Start of unit_tests_nlm_design.R\n")

## helper functions
dotest <- function(itest, observed, expected) {
	if(!identical(observed, expected)) stop(sprintf("Unit test nlm-%s failed", itest), call. = FALSE)
}
dotest_tol <- function(itest, observed, expected, tol = (.Machine$double.eps)^0.25) {
	if(any(abs(observed - expected) > tol)) stop(sprintf("Unit test nlm-%s failed", itest), call. = FALSE)
}
dotest_ineq <- function(itest, left, right, tol = (.Machine$double.eps)^0.25) {
	if(any(left > right + tol)) stop(sprintf("Unit test nlm-%s failed", itest), call. = FALSE)
}

## 1.x optimize, round methods
data <- do.call(expand.grid, replicate(3, (-1):1, simplify = FALSE))
colnames(data) <- c("x1", "x2", "x3")

obj <- nlm_design$new(
		formula = ~(b1 * x1 + b2 * x2 + b3 * x3 + b12 * x1 * x2 + b13 * x1 * x3 + b23 * x2 * x3)^2,
		theta = list(b1 = 1, b2 = 1, b3 = 1, b12 = .1, b13 = .1, b23 = .1),
		data = data
)

### d--optimal
w_a0 <- structure(rep(1/nrow(data), nrow(data)), criterion = list(crit = "A", crit.value = 0.5654714))

w_a <- obj$optimize(criterion = "A")
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("1.1.1", a_val, 1.180618)

w_a <- obj$round(m = 8, method = "efficient", seed = 1, show_progress = FALSE)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("1.1.2", a_val, 1.070862)

w_a <- obj$round(m = 8, method = "optimal", design_weights = w_a0, seed = 1, n_repeats = 1)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("1.1.3", a_val, 1.092137)

w_a <- obj$round(m = 8, method = "dist-sum", seed = 1, show_progress = FALSE)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("1.1.4", a_val, 374.6052)

w_a <- obj$round(m = 8, method = "dist-min", design_weights = w_a0, seed = 1, n_repeats = 1)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("1.1.5", a_val, 7.967434)

### ridge penalty
obj$update(lambda = 1e-2)

w_a <- obj$optimize(criterion = "A")
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("1.2.1", a_val, 1.182555)

w_a <- obj$round(m = 5, method = "optimal", seed = 1, show_progress = FALSE)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("1.2.2", a_val, 0.009942405)

obj$update(lambda = 0)

## 2.x augment method

### a-optimal
w_a01 <- structure(c(1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 1, 0, 1, 0), criterion = list(crit = "A", crit.value = 0))

obj$update(data = transform(data, augment = w_a01))

w_a <- obj$augment(design_weights = ~augment, gamma = 0.5, criterion = "A")
dotest_tol("2.1.1", obj$crit(w_a + w_a01 / 4, "A"), 1.160885)

w_a <- obj$round(m = 4, method = "optimal", seed = 1, augment_design = w_a01, show_progress = FALSE)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("2.1.2", a_val, 1.092137)

### ridge penalty
obj$update(lambda = 1e-2)

w_a <- obj$augment(design_weights = w_a01, gamma = 0.5, criterion = "A")
dotest_tol("2.2.1", obj$crit(w_a + w_a01 / 4, "A"), 1.162896)

w_a <- obj$round(m = 1, method = "optimal", seed = 1, augment_design = w_a01, show_progress = FALSE)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("2.2.2", a_val, 0.009942405)

obj$update(lambda = 0)

## 3.x cost penalty

cost <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0)
obj$update(cost = cost)

### a-optimal
w_a <- obj$optimize(criterion = "A")
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("3.1.1", a_val, 0.9551377)

w_a <- obj$round(m = 8, method = "optimal", seed = 1, show_progress = FALSE)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("3.1.2", a_val, 0.8939291)

### augment
w_a02 <- structure(c(1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0), criterion = list(crit = "A", crit.value = 0))

w_a <- obj$augment(design_weights = w_a02, gamma = 0.5, criterion = "A")
dotest_tol("3.2.1", obj$crit(w_a + w_a02 / 3, "A"), 0.8785763)

w_a <- obj$round(m = 5, method = "optimal", seed = 1, augment_design = w_a02, show_progress = FALSE)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("3.2.2", a_val, 0.8939291)

### weights
obj$update(
  cost = NULL,
  weights = 1 / (100 * pmax(cost, 0.01))
)

w_a <- obj$optimize(criterion = "A")
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("3.3.1", a_val, 0.9551391)

w_a <- obj$round(m = 8, method = "optimal", seed = 1, criterion = "A", show_progress = FALSE)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("3.3.2", a_val, 0.8939291)

w_a <- obj$augment(design_weights = w_a02, gamma = 0.5, criterion = "A")
dotest_tol("3.3.3", obj$crit(w_a + w_a02 / 3, "A"), 0.8785763)

w_a <- obj$round(m = 5, method = "optimal", seed = 1, augment_design = w_a02, show_progress = FALSE)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("3.2.2", a_val, 0.8939291)

## 4.x upper limit

obj <- nlm_design$new(
  formula = ~(b1 * x1 + b2 * x2 + b3 * x3 + b12 * x1 * x2 + b13 * x1 * x3 + b23 * x2 * x3)^2,
  theta = list(b1 = 1, b2 = 1, b3 = 1, b12 = .1, b13 = .1, b23 = .1),
  data = data,
  upper = 0.1
)

### a-optimal
w_a <- obj$optimize(criterion = "A")
dotest_tol("4.1.1", max(w_a), 0.1)

w_a <- obj$round(m = 20, method = "optimal", seed = 1, replicates = FALSE, show_progress = FALSE)
dotest_tol("4.1.2", max(w_a), 1)

w_a <- obj$round(m = 17, method = "optimal", seed = 1, replicates = 1, augment_design = w_a02, show_progress = FALSE)
dotest_tol("4.1.3", max(w_a), 1)

## 5.x stratified design
obj$update(upper = 1)

obj$stratify(
		strata = ~x1,
		proportions = c(1/4, 1/2, 1/4)
)

w_a <- obj$optimize(criterion = "A")
dotest_tol("5.1.1", tapply(w_a, obj$strata$strata, sum), c("-1" = 1/4, "0" = 1/2, "1" = 1/4))

w_a <- obj$round(m = 8, method = "optimal", seed = 1, show_progress = FALSE)
dotest_tol("5.1.2", attr(w_a, "strata"), c("-1" = 2, "0" = 4, "1" = 2))

w_a <- obj$round(m = 8, method = "dist-sum", seed = 1, show_progress = FALSE)
dotest_tol("5.1.3", attr(w_a, "strata"), c("-1" = 2, "0" = 4, "1" = 2))

w_a <- obj$round(m = 8, method = "dist-min", seed = 1, show_progress = FALSE)
dotest_tol("5.1.4", attr(w_a, "strata"), c("-1" = 2, "0" = 4, "1" = 2))

## 6.x sparse design

obj$stratify(NULL)
obj$update(alpha = 1)

entropy <- function(w, eps = .Machine$double.eps) -sum((pmax(w, 0) + eps) * log(pmax(w, 0) + eps))

w_a <- obj$sparsify(max_iter = 2, criterion = "A", reweight_fun = "entropy", show_progress = FALSE)
dotest_ineq("6.1.1", entropy(w_a), entropy(rep(1/nrow(obj$X), nrow(obj$X))))

## 7.x other methods

### crit
w_a <- obj$optimize(criterion = "A")
dotest_tol("7.1.1", obj$crit(criterion = "A"), 1.180618)
dotest_tol("7.1.2", obj$crit(w_a, criterion = "D"), 8.008704)
dotest_tol("7.1.3", obj$crit(w_a, criterion = "I"), 0.3152352)
dotest_tol("7.1.4", obj$crit(w_a, criterion = "G"), 0.06952269)
dotest_tol("7.1.5", obj$crit(w_a, criterion = "alias"), 0.7986691)

## sep
ses <- obj$sep(m = 8, sigma2 = 1)
dotest_tol("7.2.1", range(ses), c(0.000000, 1.340886))
dotest("7.2.2", length(ses), 27L)

ses1 <- obj$sep(newdata = obj$data, design_weights = w_a, sigma2 = 1)
dotest_tol("7.2.3", ses1, sqrt(8) * ses)

### vcov
w_round <- obj$round(m = 8, method = "optimal", seed = 1, show_progress = FALSE)
cov_w <- obj$vcov(m = 8, design_weights = w_a, sigma2 = 1)

dotest_tol("7.3.1", diag(cov_w), c(rep(0.01703948, 3), rep(0.01825287, 3)))
dotest_tol("7.3.2", max(abs(cov_w[lower.tri(cov_w)])), 0.005583771)

### corr
corr_w <- obj$corr(design_weights = w_a)
dotest_tol("7.4.1", max(abs(corr_w[lower.tri(corr_w)])), 0.2793079)

cat(sprintf("End of unit_tests_nlm_design.R [elapsed: %.2fs]\n", (proc.time() - start)[3]))

