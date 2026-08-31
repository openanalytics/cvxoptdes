require(cvxoptdes)

start <- proc.time()
cat("Start of unit_tests_glm_gam_design.R\n")

## helper functions
dotest <- function(itest, observed, expected) {
  if(!identical(observed, expected)) stop(sprintf("Unit test glm-gam-%s failed", itest), call. = FALSE)
}
dotest_tol <- function(itest, observed, expected, tol = (.Machine$double.eps)^0.25) {
  if(any(abs(observed - expected) > tol)) stop(sprintf("Unit test glm-gam-%s failed", itest), call. = FALSE)
}
dotest_ineq <- function(itest, left, right, tol = (.Machine$double.eps)^0.25) {
  if(any(left > right + tol)) stop(sprintf("Unit test glm-gam-%s failed", itest), call. = FALSE)
}

## 1.x optimize, round methods
data <- do.call(expand.grid, replicate(3, (-1):1, simplify = FALSE))
colnames(data) <- c("x1", "x2", "x3")

obj_glm <- glm_design$new(
  formula = ~(x1 + x2 + x3)^2 - 1,
  data = data,
  family = gaussian(link = "identity"),
  theta = c(1, 1, 1, .1, .1, .1)
)

obj_gam <- gam_design$new(
  formula = ~(x1 + x2 + x3)^2 - 1,
  data = data,
  family = gaussian(link = "identity")
)

### d-optimal
w_d0 <- structure(rep(1/nrow(data), nrow(data)), criterion = list(crit = "A", crit.value = 0.5443311))

w_d <- obj_glm$optimize(criterion = "D")
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("1.1.1", d_val, 1)

w_d <- obj_gam$optimize(criterion = "D")
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("1.1.2", d_val, 1)

w_d <- obj_glm$round(m = 10, method = "efficient", seed = 1, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("1.1.3", d_val, 0.9640569)

w_d <- obj_gam$round(m = 10, method = "efficient", seed = 1, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("1.1.4", d_val,  0.9640569)

w_d <- obj_glm$round(m = 10, method = "optimal", seed = 1, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("1.1.5", d_val,  0.9640569)

w_d <- obj_glm$round(m = 10, method = "optimal", seed = 1, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("1.1.6", d_val,  0.9640569)

w_d <- obj_glm$round(m = 10, method = "dist-sum", design_weights = w_d0, seed = 1, n_repeats = 1)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("1.1.7", d_val, 159.1384)

w_d <- obj_gam$round(m = 10, method = "dist-min", seed = 1, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("1.1.8", d_val, 2)

### ridge penalty
obj_glm$update(lambda = 1e-2)

w_d <- obj_glm$optimize(criterion = "D")
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("1.2.1", d_val, 0.99)

w_d <- obj_glm$round(m = 5, method = "optimal", seed = 1, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_ineq("1.2.2", 0.4901072, d_val)

obj_glm$update(lambda = 0)

## 2.x augment method

w_d01 <- structure(c(1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 1), criterion = list(crit = "D", crit.value = 0))

obj_glm$update(data = transform(data, augment = w_d01))
obj_gam$update(data = transform(data, augment = w_d01))

w_d <- obj_glm$augment(design_weights = w_d01, gamma = 0.5, criterion = "D")
dotest_tol("2.1.1", obj_glm$crit(w_d + w_d01 / 5, "D"), 1)

w_d <- obj_gam$augment(design_weights = w_d01, gamma = 0.5, criterion = "D")
dotest_tol("2.1.2", obj_gam$crit(w_d + w_d01 / 5, "D"), 1)

w_d <- obj_glm$round(m = 5, method = "optimal", seed = 1, augment_design = w_d01, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("2.1.3", d_val, 0.9640569)

w_d <- obj_gam$round(m = 5, method = "optimal", seed = 1, augment_design = w_d01, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("2.1.4", d_val,  0.9640569)

## 3.x cost penalty
cost <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1)
obj_glm$update(cost = cost)
obj_gam$update(cost = cost)

w_d <- obj_glm$optimize(criterion = "D")
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("3.1.1", d_val, 0.7217672)

w_d <- obj_gam$optimize(criterion = "D")
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("3.1.2", d_val, 0.7217672)

w_d <- obj_glm$round(m = 10, method = "optimal", seed = 1, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("3.1.3", d_val, 0.6970324)

w_d <- obj_gam$round(m = 10, method = "optimal", seed = 1, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("3.1.4", d_val, 0.6970324)

### augment
w_d02 <- structure(c(1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0), criterion = list(crit = "D", crit.value = 0))

w_d <- obj_glm$round(m = 6, method = "optimal", seed = 1, augment_design = w_d02, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("3.2.1", d_val, 0.6970324)

### weights
obj_glm$update(
  cost = NULL,
  weights = 1 / (100 * pmax(cost, 0.01))
)
obj_gam$update(
  cost = NULL,
  weights = 1 / (100 * pmax(cost, 0.01))
)

w_d <- obj_glm$optimize(criterion = "D")
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("3.3.1", d_val, 0.7190719)

w_d <- obj_gam$optimize(criterion = "D")
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("3.3.2", d_val, 0.7190719)

w_d <- obj_glm$round(m = 10, method = "optimal", seed = 1, criterion = "D", show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("3.3.3", d_val, 0.6970324)

w_d <- obj_gam$round(m = 10, method = "optimal", seed = 1, criterion = "D", show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("3.3.4", d_val, 0.6970324)

w_d <- obj_glm$round(m = 6, method = "optimal", seed = 1, design_weights = w_d, augment_design = w_d02, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("3.3.5", d_val, 0.6970324)

w_d <- obj_gam$round(m = 6, method = "optimal", seed = 1, design_weights = w_d, augment_design = w_d02, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("3.3.6", d_val, 0.6970324)

## 4.x upper limit

obj_glm$update(
  weights = NULL,
  upper = 0.1
)

w_d <- obj_glm$optimize(criterion = "D")
dotest_tol("4.1.1", max(w_d), 0.1)

w_d <- obj_glm$round(m = 20, method = "optimal", seed = 1, replicates = FALSE, show_progress = FALSE)
dotest_tol("4.1.2", max(w_d), 1)

w_d <- obj_glm$round(m = 16, method = "optimal", seed = 1, replicates = 1, augment_design = w_d02, show_progress = FALSE)
dotest_tol("4.1.3", max(w_d), 1)

## 5.x stratified design
obj_glm <- glm_design$new(
  formula = ~(x1 + x2 + x3)^2 - 1,
  data = data,
  family = gaussian(link = "identity"),
  theta = c(1, 1, 1, .1, .1, .1)
)

obj_gam <- gam_design$new(
  formula = ~(x1 + x2 + x3)^2 - 1,
  data = data,
  family = gaussian(link = "identity")
)

obj_glm$stratify(
  strata = ~x1,
  proportions = c(1/4, 1/2, 1/4)
)
obj_gam$stratify(
  strata = ~x1,
  proportions = c(1/4, 1/2, 1/4)
)

w_d <- obj_glm$optimize(criterion = "D")
dotest_tol("5.1.1", tapply(w_d, obj_glm$strata$strata, sum), c("-1" = 1/4, "0" = 1/2, "1" = 1/4))

w_d <- obj_gam$optimize(criterion = "D")
dotest_tol("5.1.2", tapply(w_d, obj_gam$strata$strata, sum), c("-1" = 1/4, "0" = 1/2, "1" = 1/4))

w_d <- obj_glm$round(m = 8, method = "optimal", seed = 1, show_progress = FALSE)
dotest_tol("5.1.3", attr(w_d, "strata"), c("-1" = 2, "0" = 4, "1" = 2))

w_d <- obj_gam$round(m = 8, method = "optimal", seed = 1, show_progress = FALSE)
dotest_tol("5.1.4", attr(w_d, "strata"), c("-1" = 2, "0" = 4, "1" = 2))

w_d <- obj_glm$round(m = 8, method = "dist-sum", seed = 1, show_progress = FALSE)
dotest_tol("5.1.5", attr(w_d, "strata"), c("-1" = 2, "0" = 4, "1" = 2))

w_d <- obj_gam$round(m = 8, method = "dist-min", seed = 1, show_progress = FALSE)
dotest_tol("5.1.6", attr(w_d, "strata"), c("-1" = 2, "0" = 4, "1" = 2))

## 6.x sparse design

obj_glm$stratify(NULL)
obj_gam$stratify(NULL)

entropy <- function(w, eps = .Machine$double.eps) -sum((pmax(w, 0) + eps) * log(pmax(w, 0) + eps))

w_a <- obj_glm$sparsify(max_iter = 2, criterion = "A", reweight_fun = "entropy", show_progress = FALSE)
dotest_ineq("6.1.1", entropy(w_a), entropy(rep(1/nrow(obj_glm$X), nrow(obj_glm$X))))

w_a <- obj_gam$sparsify(max_iter = 2, criterion = "A", reweight_fun = "entropy", show_progress = FALSE)
dotest_ineq("6.1.2", entropy(w_a), entropy(rep(1/nrow(obj_gam$X), nrow(obj_gam$X))))

## 7.x other methods

### crit
w_d <- obj_glm$optimize(criterion = "D")
w_d <- obj_gam$optimize(criterion = "D")
dotest_tol("7.1.1", obj_glm$crit(criterion = "D"), 1)
dotest_tol("7.1.2", obj_gam$crit(criterion = "D"), 1)

## sep
ses <- obj_glm$sep(m = 10, sigma2 = 1)
dotest_tol("7.2.1", range(ses), c(0, 0.7745971))
dotest("7.2.2", length(ses), 27L)

ses1 <- obj_glm$sep(newdata = obj_gam$data, design_weights = w_d, sigma2 = 1)
dotest_tol("7.2.3", ses1, sqrt(10) * ses)

ses <- obj_gam$sep(m = 10, sigma2 = 1, type = "response")
dotest_tol("7.2.4", range(ses), c(0, 0.7745971))
dotest("7.2.5", length(ses), 27L)

ses1 <- obj_gam$sep(newdata = obj_gam$data, design_weights = w_d, sigma2 = 1, type = "response")
dotest_tol("7.2.6", ses1, sqrt(10) * ses)

### vcov
w_round <- obj_glm$round(m = 10, method = "optimal", seed = 1, show_progress = FALSE)
dotest_tol("7.3.1", obj_glm$vcov(m = 10, sigma2 = 1), diag(0.1, 6))
dotest_tol("7.3.2", diag(obj_gam$vcov(design_weights = w_round, sigma2 = 1)), rep(0.1071429, 6))
dotest_tol("7.3.3", sum(abs(obj_gam$vcov(design_weights = w_round, sigma2 = 1))), 0.8571429)

### corr
dotest_tol("7.4.1", obj_glm$corr(), diag(6))
dotest_tol("7.4.2", sum(abs(obj_glm$corr(design_weights = w_round, center = FALSE))), 8.4)
dotest_tol("7.4.3", sum(abs(obj_gam$corr(design_weights = w_round, center = TRUE))), 8.2)

### subset
design <- obj_glm$subset(w_round)
dotest("7.5.1", dim(design), c(10L, 3L))
dotest("7.5.2", colnames(design), c("x1", "x2", "x3"))

cat(sprintf("End of unit_tests_gam_glm_design.R [elapsed: %.2fs]\n", (proc.time() - start)[3]))

