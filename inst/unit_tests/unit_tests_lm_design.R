require(cvxoptdes)

start <- proc.time()
cat("Start of unit_tests_lm_design.R\n")

## helper functions
dotest <- function(itest, observed, expected) {
  if(!identical(observed, expected)) stop(sprintf("Unit test lm-%s failed", itest), call. = FALSE)
}
dotest_tol <- function(itest, observed, expected, tol = (.Machine$double.eps)^0.25) {
  if(any(abs(observed - expected) > tol)) stop(sprintf("Unit test lm-%s failed", itest), call. = FALSE)
}
dotest_ineq <- function(itest, left, right, tol = (.Machine$double.eps)^0.25) {
  if(any(left > right + tol)) stop(sprintf("Unit test lm-%s failed", itest), call. = FALSE)
}

## 1.x optimize, round methods
data <- do.call(expand.grid, replicate(3, (-1):1, simplify = FALSE))
colnames(data) <- c("x1", "x2", "x3")

obj <- lm_design$new(
  formula = ~(x1 + x2 + x3)^2 - 1,
  data = data
)

### d-optimal
w_d0 <- structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, .1, 0, .2, 0, 0, 0, .1,
                    0, .2, 0, .2, 0, 0, 0, .1, 0, .1), criterion = list(crit = "D", crit.value = 0.4))

w_d <- obj$optimize(criterion = "D")
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("1.1.1", d_val, 1)

w_d <- obj$round(m = 10, method = "efficient", seed = 1, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_ineq("1.1.2", 0.9607496, d_val)

w_d <- obj$round(m = 10, method = "optimal", seed = 1, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("1.1.3", d_val, 0.9640569)

w_d <- obj$round(m = 10, method = "dist-sum", design_weights = w_d0, seed = 1, n_repeats = 1)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("1.1.4", d_val, 159.1384)

w_d <- obj$round(m = 10, method = "dist-min", design_weights = w_d0, seed = 1, n_repeats = 1)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("1.1.5", d_val, 2)

### a-optimal
w_a0 <- structure(w_d0, criterion = list(crit = "A", crit.value = 0.03636364))

w_a <- obj$optimize(criterion = "A")
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("1.2.1", a_val, 0.166666)

w_a <- obj$round(m = 10, method = "efficient", seed = 1, show_progress = FALSE)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("1.2.2", a_val, 0.1555556)

w_a <- obj$round(m = 10, method = "optimal", design_weights = w_a0, seed = 1, show_progress = FALSE)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("1.2.3", a_val, 0.1555556)

w_a <- obj$round(m = 10, method = "dist-sum", design_weights = w_a0, seed = 1, n_repeats = 1)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("1.2.4", a_val, 159.1384)

w_a <- obj$round(m = 10, method = "dist-min", design_weights = w_a0, seed = 1, n_repeats = 1)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("1.2.5", a_val, 2)

### i-optimal
w_i0 <- structure(w_d0, criterion = list(crit = "I", crit.value = 0.06728972))

w_i <- obj$optimize(criterion = "I")
i_val <- attr(w_i, "criterion")$crit.value
dotest_tol("1.3.1", i_val, 0.3)

w_i <- obj$round(m = 10, method = "efficient", seed = 1, show_progress = FALSE)
i_val <- attr(w_i, "criterion")$crit.value
dotest_tol("1.3.2", i_val, 0.28)

w_i <- obj$round(m = 10, method = "optimal", design_weights = w_i0, seed = 1, show_progress = FALSE)
i_val <- attr(w_i, "criterion")$crit.value
dotest_tol("1.3.3", i_val, 0.28)

w_i <- obj$round(m = 10, method = "dist-sum", design_weights = w_i0, seed = 1, n_repeats = 1)
i_val <- attr(w_i, "criterion")$crit.value
dotest_tol("1.3.4", i_val, 159.13847)

w_i <- obj$round(m = 10, method = "dist-min", design_weights = w_i0, seed = 1, n_repeats = 1)
i_val <- attr(w_i, "criterion")$crit.value
dotest_tol("1.3.5", i_val, 2)

### g-optimal
w_g0 <- structure(w_d0, criterion = list(crit = "G", crit.value = 0.02105263))

w_g <- obj$optimize(criterion = "G")
g_val <- attr(w_g, "criterion")$crit.value
dotest_tol("1.4.1", g_val, 0.166667)

w_g <- obj$round(m = 10, method = "efficient", seed = 1, show_progress = FALSE)
g_val <- attr(w_g, "criterion")$crit.value
dotest_tol("1.4.2", g_val, 0.14)

w_g <- obj$round(m = 10, method = "optimal", design_weights = w_g0, seed = 1, show_progress = FALSE)
g_val <- attr(w_g, "criterion")$crit.value
dotest_tol("1.4.3", g_val, 0.14)

w_g <- obj$round(m = 10, method = "dist-sum", design_weights = w_g0, seed = 1, n_repeats = 1)
g_val <- attr(w_g, "criterion")$crit.value
dotest_tol("1.4.4", g_val, 159.1384)

w_g <- obj$round(m = 10, method = "dist-min", seed = 1, n_repeats = 1)
g_val <- attr(w_g, "criterion")$crit.value
dotest_tol("1.4.5", g_val, 2)

### alias-optimal
w_alias0 <-  structure(w_d0, criterion = list(crit = "alias", crit.value = 0.8015955))

w_alias <- obj$optimize(criterion = "alias")
alias_val <- attr(w_alias, "criterion")$crit.value
dotest_tol("1.5.1", alias_val, 1)

w_alias <- obj$round(m = 10, method = "optimal", seed = 1, n_repeats = 100, show_progress = FALSE)
alias_val <- attr(w_alias, "criterion")$crit.value
dotest_tol("1.5.3", alias_val, 1)

### ridge penalty
obj$update(lambda = 1e-2)

w_d <- obj$optimize(criterion = "D")
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("1.6.1", d_val, 0.99)

w_d <- obj$round(m = 5, method = "optimal", seed = 1, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("1.6.2", d_val, 0.4901072)

w_i <- obj$optimize(criterion = "I")
i_val <- attr(w_i, "criterion")$crit.value
dotest_tol("1.6.3", i_val, 0.303)

w_i <- obj$round(m = 5, method = "optimal", seed = 1, show_progress = FALSE)
i_val <- attr(w_i, "criterion")$crit.value
dotest_tol("1.6.4", i_val, 0.0206107)

w_a <- obj$optimize(criterion = "A")
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("1.6.5", a_val, 0.1683347)

w_a <- obj$round(m = 5, method = "optimal", seed = 1, show_progress = FALSE)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("1.6.6", a_val, 0.009475355)

w_g <- obj$optimize(criterion = "G")
g_val <- attr(w_g, "criterion")$crit.value
dotest_tol("1.6.7", g_val, 0.166667)

w_g <- obj$round(m = 5, method = "optimal", seed = 1, show_progress = FALSE)
g_val <- attr(w_g, "criterion")$crit.value
dotest_tol("1.6.8", g_val, 0.004741531)

w_alias <- obj$optimize(criterion = "alias")
alias_val <- attr(w_alias, "criterion")$crit.value
dotest_tol("1.6.9", alias_val, 1)

## 2.x augment method

obj <- lm_design$new(
  formula = ~(x1 + x2 + x3)^2 - 1,
  data = data
)

### d-optimal
w_d01 <- structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                     0, 1, 0, 1, 0, 0, 0, 1, 0, 1), criterion = list(crit = "D", crit.value = 0))

w_d <- obj$augment(design_weights = w_d01, gamma = 0.5, criterion = "D")
dotest_tol("2.1.1", obj$crit(w_d + w_d01 / 5, "D"), 0.9159439)

w_d <- obj$round(m = 4, method = "efficient", seed = 1, augment_design = w_d01, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("2.1.2", d_val, 0.9065106)

w_d <- obj$round(m = 4, method = "optimal", seed = 1, augment_design = w_d01, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("2.1.3", d_val, 0.9065106)

w_d <- obj$round(m = 4, method = "dist-sum", seed = 1, augment_design = w_d01, n_repeats = 1)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("2.1.4", d_val, 124.3699)

w_d <- obj$round(m = 5, method = "dist-min", seed = 1, augment_design = w_d01, n_repeats = 1)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("2.1.5", d_val, 2)

### i-optimal
w_i01 <- structure(w_d01, criterion = list(crit = "I", crit.value = 0))

data_w <- transform(data, weights = w_i01)
obj$update(data = data_w)

w_i <- obj$augment(design_weights = ~weights, gamma = 0.5, criterion = "I")
dotest_tol("2.2.1", obj$crit(w_i + data_w$weights / 5, "I"), 0.2755102)

w_i <- obj$round(m = 4, method = "optimal", seed = 1, augment_design = w_i01, show_progress = FALSE)
i_val <- attr(w_i, "criterion")$crit.value
dotest_tol("2.2.2", i_val, 0.2727273)

### a-optimal
w_a01 <- structure(w_d01, criterion = list(crit = "A", crit.value = 0))

w_a <- obj$augment(design_weights = w_a01, gamma = 0.5, criterion = "A")
dotest_tol("2.3.1", obj$crit(w_a + w_a01 / 5, "A"), 0.1525424)

w_a <- obj$round(m = 4, method = "optimal", seed = 1, augment_design = w_a01, show_progress = FALSE)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("2.3.2", a_val, 0.1509434)

### g-optimal
w_g01 <- structure(w_d01, criterion = list(crit = "G", crit.value = 0))

w_g <- obj$augment(design_weights = w_g01, gamma = 0.5, criterion = "G")
dotest_tol("2.4.1", obj$crit(w_g + w_g01 / 5, "G"),  0.1525441)

w_g <- obj$round(m = 4, method = "optimal", seed = 1, augment_design = w_g01, show_progress = FALSE)
g_val <- attr(w_g, "criterion")$crit.value
dotest_tol("2.4.2", g_val, 0.1509434)

### alias-optimal
w_alias01 <- structure(w_d01, criterion = list(crit = "alias", crit.value = 0.8737049))

w_alias <- obj$augment(design_weights = w_alias01, gamma = 0.5, criterion = "alias")
dotest_tol("2.5.1", obj$crit(w_alias + w_alias01 / 5, "alias"), 1)

w_alias <- obj$round(m = 4, method = "optimal", seed = 1, augment_design = w_alias01, show_progress = FALSE)
alias_val <- attr(w_alias, "criterion")$crit.value
dotest_tol("2.5.2", alias_val, 1)

### ridge penalty

obj$update(lambda = 1e-2)

w_d <- obj$augment(design_weights = w_d01, gamma = 0.5, criterion = "D")
dotest_tol("2.6.1", obj$crit(w_d + w_d01 / 5, "D"), 0.9259536)

w_d <- obj$round(m = 4, method = "optimal", seed = 1, augment_design = w_d01, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("2.6.2", d_val, 0.9165199)

w_i <- obj$augment(design_weights = w_i01, gamma = 0.5, criterion = "I")
dotest_tol("2.6.3", obj$crit(w_i + data_w$weights / 5, "I"),  0.2785151)

w_i <- obj$round(m = 4, method = "optimal", seed = 1, augment_design = w_i01, show_progress = FALSE)
i_val <- attr(w_i, "criterion")$crit.value
dotest_tol("2.6.4", i_val, 0.2757334)

w_a <- obj$augment(design_weights = w_a01, gamma = 0.5, criterion = "A")
dotest_tol("2.6.5", obj$crit(w_a + w_a01 / 5, "A"), 0.1542114)

w_a <- obj$round(m = 4, method = "optimal", seed = 1, augment_design = w_a01, show_progress = FALSE)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("2.6.6", a_val, 0.152613)

w_g <- obj$augment(design_weights = w_g01, gamma = 0.5, criterion = "G")
dotest_tol("2.6.7", obj$crit(w_g + w_g01 / 5, "G"),  0.1542118)

w_g <- obj$round(m = 4, method = "optimal", seed = 1, augment_design = w_g01, show_progress = FALSE)
g_val <- attr(w_g, "criterion")$crit.value
dotest_tol("2.6.8", g_val, 0.152613)

## 3.x cost penalty

data_cost <- transform(data, cost = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                                      0, 1, 0, 1, 0, 0, 0, 1, 0, 1))

obj <- lm_design$new(
  formula = ~(x1 + x2 + x3)^2 - 1,
  data = data_cost,
  cost = ~cost,
  alpha = 1
)

w_d <- obj$optimize(criterion = "D")
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("3.1.1", d_val, 0.7217662)

w_d <- obj$round(m = 10, method = "optimal", seed = 1, criterion = "D", show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("3.1.2", d_val, 0.6970324)

w_a <- obj$optimize(criterion = "A")
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("3.1.3", a_val, 0.1111112)

w_a <- obj$round(m = 10, method = "optimal", seed = 1, criterion = "A", show_progress = FALSE)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("3.1.4", a_val, 0.105)

w_i <- obj$optimize(criterion = "I")
i_val <- attr(w_i, "criterion")$crit.value
dotest_tol("3.1.5", i_val, 0.2047207)

w_i <- obj$round(m = 10, method = "optimal", seed = 1, criterion = "I", show_progress = FALSE)
i_val <- attr(w_i, "criterion")$crit.value
dotest_tol("3.1.6", i_val, 0.1935979)

w_g <- obj$optimize(criterion = "G")
g_val <- attr(w_g, "criterion")$crit.value
dotest_tol("3.1.7", g_val, 0.07344086)

w_g <- obj$round(m = 10, method = "optimal", seed = 1, criterion = "G", show_progress = FALSE)
g_val <- attr(w_g, "criterion")$crit.value
dotest_tol("3.1.8", g_val, 0.09090909)

w_alias <- obj$optimize(criterion = "alias")
alias_val <- obj$crit(w_alias, "alias")
dotest_tol("3.1.9", alias_val, 0.8722714)

w_alias <- obj$round(m = 10, method = "optimal", seed = 1, criterion = "alias", show_progress = FALSE)
alias_val <- attr(w_alias, "criterion")$crit.value
dotest_tol("3.1.10", alias_val, 0.9777777)

### augment
w_d02 <- structure(c(1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0), criterion = list(crit = "D", crit.value = 0))

w_d <- obj$round(m = 6, method = "optimal", seed = 1, design_weights = w_d, augment_design = w_d02, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("3.2.1", d_val, 0.6970324)

w_a <- obj$round(m = 6, method = "optimal", seed = 1, design_weights = w_a, augment_design = w_d02, show_progress = FALSE)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("3.2.2", a_val, 0.105)

w_i <- obj$round(m = 6, method = "optimal", seed = 1, design_weights = w_i, augment_design = w_d02, n_repeats = 25, show_progress = FALSE)
i_val <- attr(w_i, "criterion")$crit.value
dotest_ineq("3.2.3", 0.193598, i_val)

w_g <- obj$round(m = 6, method = "optimal", seed = 1, design_weights = w_g, augment_design = w_d02, show_progress = FALSE)
g_val <- attr(w_g, "criterion")$crit.value
dotest_tol("3.2.4", g_val, 0.09090909)

## 3.x weights
obj$update(
  cost = NULL,
  weights = ~ 1 / (100 * pmax(cost, 0.01))
)

w_d <- obj$optimize(criterion = "D")
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("3.3.1", d_val, 0.7190719)

w_d <- obj$round(m = 10, method = "optimal", seed = 1, criterion = "D", show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("3.3.2", d_val, 0.6970324)

w_a <- obj$optimize(criterion = "A")
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("3.3.3", a_val, 0.1111112)

w_a <- obj$round(m = 10, method = "optimal", seed = 1, criterion = "A", show_progress = FALSE)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("3.3.4", a_val, 0.105)

w_i <- obj$optimize(criterion = "I")
i_val <- attr(w_i, "criterion")$crit.value
dotest_tol("3.3.5", i_val, 0.3147453)

w_i <- obj$round(m = 10, method = "optimal", seed = 1, criterion = "I", show_progress = FALSE)
i_val <- attr(w_i, "criterion")$crit.value
dotest_tol("3.3.6", i_val, 0.2959766)

w_g <- obj$optimize(criterion = "G")
g_val <- attr(w_g, "criterion")$crit.value
dotest_tol("3.3.7", g_val, 0.1666666)

w_g <- obj$round(m = 10, method = "optimal", seed = 1, criterion = "G", show_progress = FALSE)
g_val <- attr(w_g, "criterion")$crit.value
dotest_tol("3.3.8", g_val, 0.1263158)

w_alias <- obj$optimize(criterion = "alias")
alias_val <- obj$crit(w_alias, "alias")
dotest_tol("3.3.9", alias_val, 0.8722714)

w_alias <- obj$round(m = 10, method = "optimal", seed = 1, criterion = "alias", show_progress = FALSE)
alias_val <- attr(w_alias, "criterion")$crit.value
dotest_ineq("3.3.10", 0.95, alias_val)

### augment
w_d <- obj$round(m = 6, method = "optimal", seed = 1, design_weights = w_d, augment_design = w_d02, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("3.4.1", d_val, 0.6970324)

w_a <- obj$round(m = 6, method = "optimal", seed = 1, design_weights = w_a, augment_design = w_d02, show_progress = FALSE)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("3.4.2", a_val, 0.105)

w_i <- obj$round(m = 6, method = "optimal", seed = 1, design_weights = w_i, augment_design = w_d02, show_progress = FALSE)
i_val <- attr(w_i, "criterion")$crit.value
dotest_tol("3.4.3", i_val, 0.2959766)

w_g <- obj$round(m = 6, method = "optimal", seed = 1, design_weights = w_g, augment_design = w_d02, show_progress = FALSE)
g_val <- attr(w_g, "criterion")$crit.value
dotest_tol("3.4.4", g_val, 0.1263158)

w_alias <- obj$round(m = 6, method = "optimal", seed = 1, design_weights = w_alias - w_d02, augment_design = w_d02, show_progress = FALSE)
alias_val <- attr(w_alias, "criterion")$crit.value
dotest_ineq("3.4.5", 0.9588524, alias_val)

## 4.x upper limit
obj <- lm_design$new(
  formula = ~(x1 + x2 + x3)^2 - 1,
  data = data,
  upper = 0.1
)

### d-optimal
w_d <- obj$optimize(criterion = "D")
dotest_tol("4.1.1", max(w_d), 0.1)

w_d <- obj$round(m = 19, method = "optimal", seed = 1, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("4.1.2", d_val, 0.701435)

w_d <- obj$round(m = 20, method = "optimal", seed = 1, replicates = FALSE, show_progress = FALSE)
dotest_tol("4.1.3", max(w_d), 1)

w_d <- obj$round(m = 19, method = "optimal", seed = 1, replicates = 2, show_progress = FALSE)
dotest_tol("4.1.4", max(w_d), 2)

w_d <- obj$round(m = 19, method = "dist-sum", seed = 1, show_progress = FALSE)
dotest_tol("4.1.5", max(w_d), 1)

w_d <- obj$round(m = 14, method = "optimal", seed = 1, replicates = 1, augment_design = w_d01, show_progress = FALSE)
dotest_tol("4.1.6", max(w_d), 1)

### a-optimal
w_a <- obj$optimize(criterion = "A")
dotest_tol("4.2.1", max(w_a), 0.1)

w_a <- obj$round(m = 10, method = "optimal", seed = 1, show_progress = FALSE)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("4.2.2", a_val, 0.1463415)

### i-optimal
w_i <- obj$optimize(criterion = "I")
dotest_tol("4.2.3", max(w_i), 0.1)

w_i <- obj$round(m = 10, method = "optimal", seed = 1, show_progress = FALSE)
i_val <- attr(w_i, "criterion")$crit.value
dotest_tol("4.2.4", i_val, 0.2653563)

### g-optimal
w_g <- obj$optimize(criterion = "G")
dotest_tol("4.2.5", max(w_g), 0.1)

w_g <- obj$round(m = 10, method = "optimal", seed = 1, show_progress = FALSE)
g_val <- attr(w_g, "criterion")$crit.value
dotest_tol("4.2.6", g_val, 0.1384615)

### alias-optimal
w_alias <- obj$optimize(criterion = "alias")
dotest_tol("4.2.7", max(w_alias), 0.1)

w_alias <- obj$round(m = 10, method = "optimal", seed = 1, n_repeats = 100, show_progress = FALSE)
alias_val <- attr(w_alias, "criterion")$crit.value
dotest_tol("4.2.8", alias_val, 1)

## 5.x stratified design
obj <- lm_design$new(
  formula = ~(x1 + x2 + x3)^2 - 1,
  data = data
)

obj$stratify(
  strata = ~x1,
  proportions = c(1/4, 1/2, 1/4)
)

w_d <- obj$optimize(criterion = "D")
dotest_tol("5.1.1", tapply(w_d, obj$strata$strata, sum), c("-1" = 1/4, "0" = 1/2, "1" = 1/4))

w_d <- obj$round(m = 12, method = "optimal", seed = 1, show_progress = FALSE)
dotest_tol("5.1.2", attr(w_d, "strata"), c("-1" = 3, "0" = 6, "1" = 3))

w_a <- obj$optimize(criterion = "A")
dotest_tol("5.2.1", tapply(w_a, obj$strata$strata, sum), c("-1" = 1/4, "0" = 1/2, "1" = 1/4))

w_a <-obj$round(m = 11, method = "optimal", seed = 1, show_progress = FALSE)
dotest_tol("5.2.2", attr(w_a, "strata"), c("-1" = 3, "0" = 5, "1" = 3))

w_i <- obj$optimize(criterion = "I")
dotest_tol("5.2.3", attr(w_i, "strata"), c("-1" = 1/4, "0" = 1/2, "1" = 1/4))

w_i <- obj$round(m = 12, method = "optimal", seed = 1, show_progress = FALSE)
dotest_tol("5.2.4", attr(w_i, "strata"), c("-1" = 3, "0" = 6, "1" = 3))

w_g <- obj$optimize(criterion = "G")
dotest_tol("5.2.5", tapply(w_g, obj$strata$strata, sum), c("-1" = 1/4, "0" = 1/2, "1" = 1/4))

w_g <- obj$round(m = 15, method = "optimal", seed = 1, show_progress = FALSE)
dotest_tol("5.2.6", attr(w_g, "strata"), c("-1" = 4, "0" = 7, "1" = 4))

w_alias <- obj$optimize(criterion = "alias")
dotest_tol("5.2.7", tapply(w_alias, obj$strata$strata, sum), c("-1" = 1/4, "0" = 1/2, "1" = 1/4))

w_alias <- obj$round(m = 12, method = "optimal", seed = 1, show_progress = FALSE)
dotest_tol("5.2.8", attr(w_alias, "strata"), c("-1" = 3, "0" = 6, "1" = 3))

## 6.x orthogonal design

obj <- lm_design$new(
  formula = ~0 + (x1 + x2 + x3)^2 + I(x1^2) + I(x2^2) + I(x3^2),
  data = data,
  orthogonal = list(
    x1 + x2 + x3 ~ . ~ 0
  ),
  zeros = list(
    x1 + x2 + x3 + I(x1^2) + I(x2^2) + I(x3^2) ~ 0.25,
    (x1 + x2 + x3)^2 - (x1 + x2 + x3) ~ 0.5
  ),
  col_means = x1 + x2 + x3 ~ 0
)

### d-optimal
w_d <- obj$optimize(criterion = "D")
dotest_ineq("6.1.1", obj$corr(w_d), obj$orthogonal)
dotest_tol("6.1.2", apply(obj$X, 2, function(x) sum(w_d[x == 0])), c(rep(0.25, 6), rep(0.5, 3)))
dotest_tol("6.1.3", colSums(w_d * obj$X[, 1:3]), rep(0, 3))

if(requireNamespace("highs", quietly = TRUE)) {

  w_d <- obj$round(m = 8, method = "orthogonal", highs_control = list(log_to_console = FALSE))
  dotest_ineq("6.2.1", obj$corr(w_d), obj$orthogonal)
  dotest_tol("6.2.2", apply(obj$X, 2, function(x) sum(w_d[x == 0])), c(rep(2, 6), rep(4, 3)))
  dotest_tol("6.2.3", colSums(w_d * obj$X[, 1:3]), rep(0, 3))

  w_d03 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0)

  w_d <- obj$round(m = 4, method = "orthogonal", augment_design = w_d03, highs_control = list(log_to_console = FALSE))
  dotest_ineq("6.2.4", obj$corr(w_d), obj$orthogonal)
  dotest_tol("6.2.5", apply(obj$X, 2, function(x) sum(w_d[x == 0])), c(rep(2, 6), rep(4, 3)))
  dotest_tol("6.2.6", colSums(w_d * obj$X[, 1:3]), rep(0, 3))

  w_d <- obj$round(m = 4, method = "orthogonal", augment_design = w_d03,
                   design_weights = c(0, 0, 0.25, 0.25, 0, 0, 0, 0.25, 0, 0.25, 0, 0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 0, 0, 0, 0), highs_control = list(log_to_console = FALSE))
  dotest_ineq("6.2.7", obj$corr(w_d), obj$orthogonal)
  dotest_tol("6.2.8", apply(obj$X, 2, function(x) sum(w_d[x == 0])), c(rep(2, 6), rep(4, 3)))
  dotest_tol("6.2.9", colSums(w_d * obj$X[, 1:3]), rep(0, 3))

} else {
  cat("Skipping orthogonal rounding unit tests, HiGHS solver not available.\n")
}

### optimality criteria
w_a <- obj$optimize(criterion = "A")
dotest_ineq("6.3.1", obj$corr(w_a), obj$orthogonal)
dotest_tol("6.3.2", apply(obj$X, 2, function(x) sum(w_a[x == 0])), c(rep(0.25, 6), rep(0.5, 3)))
dotest_tol("6.3.3", colSums(w_a * obj$X[, 1:3]), rep(0, 3))

w_i <- obj$optimize(criterion = "I")
dotest_ineq("6.3.4", obj$corr(w_i), obj$orthogonal)
dotest_tol("6.3.5", apply(obj$X, 2, function(x) sum(w_i[x == 0])), c(rep(0.25, 6), rep(0.5, 3)))
dotest_tol("6.3.6", colSums(w_i * obj$X[, 1:3]), rep(0, 3))

w_g <- obj$optimize(criterion = "G")
dotest_ineq("6.3.7", obj$corr(w_g), obj$orthogonal)
dotest_tol("6.3.8", apply(obj$X, 2, function(x) sum(w_g[x == 0])), c(rep(0.25, 6), rep(0.5, 3)))
dotest_tol("6.3.9", colSums(w_g * obj$X[, 1:3]), rep(0, 3))

w_alias <- obj$optimize(criterion = "alias")
dotest_ineq("6.3.10", obj$corr(w_alias), obj$orthogonal)
dotest_tol("6.3.11", apply(obj$X, 2, function(x) sum(w_alias[x == 0])), c(rep(0.25, 6), rep(0.5, 3)))
dotest_tol("6.3.12", colSums(w_alias * obj$X[, 1:3]), rep(0, 3))

### custom dot product constraints
obj <- obj$update(
  orthogonal = x1 + x2 + x3 ~ . ~ 0.1,
  zeros = NULL,
  col_means = NULL
)

w_d <- obj$optimize(criterion = "D")
dotest_ineq("6.4.1", obj$corr(w_d), obj$orthogonal)

w_a <- obj$optimize(criterion = "A")
dotest_ineq("6.4.2", obj$corr(w_a), obj$orthogonal)

w_i <- obj$optimize(criterion = "I")
dotest_ineq("6.4.3", obj$corr(w_i), obj$orthogonal)

w_g <- obj$optimize(criterion = "G")
dotest_ineq("6.4.5", obj$corr(w_g), obj$orthogonal)

w_alias <- obj$optimize(criterion = "alias")
dotest_ineq("6.4.6", obj$corr(w_alias), obj$orthogonal)

## 7.x sparse design

obj <- lm_design$new(
  formula = ~(x1 + x2 + x3)^2 - 1,
  data = data,
  alpha = 50
)

entropy <- function(w, eps = .Machine$double.eps) -sum((pmax(w, 0) + eps) * log(pmax(w, 0) + eps))
logsum <- function(w, eps = .Machine$double.eps) sum(log(pmax(w, 0) + sqrt(eps)))
powerp <- function(w, eps = .Machine$double.eps, p = 0.5) sum((pmax(w, 0) + eps)^p)

### alias-optimal
w_alias <- obj$sparsify(max_iter = 2, criterion = "alias", reweight_fun = "entropy", show_progress = FALSE)
dotest_ineq("7.1.1", entropy(w_alias), entropy(rep(1/nrow(obj$X), nrow(obj$X))))

w_alias <- obj$sparsify(max_iter = 2, criterion = "alias", reweight_fun = "log-sum", show_progress = FALSE)
dotest_ineq("7.1.2", logsum(w_alias), logsum(rep(1/nrow(obj$X), nrow(obj$X))))

w_alias <- obj$sparsify(max_iter = 2, criterion = "alias", reweight_fun = "power", p = 0.5, show_progress = FALSE)
dotest_ineq("7.1.3", powerp(w_alias), powerp(rep(1/nrow(obj$X), nrow(obj$X))))

## approximate omars
obj$update(
  formula = ~0 + (x1 + x2 + x3)^2 + I(x1^2) + I(x2^2) + I(x3^2),
  orthogonal = list(
    x1 + x2 + x3 ~ . ~ 0
  ),
  zeros = list(
    x1 + x2 + x3 + I(x1^2) + I(x2^2) + I(x3^2) ~ 0.25,
    (x1 + x2 + x3)^2 - (x1 + x2 + x3) ~ 0.5
  ),
  col_means = x1 + x2 + x3 ~ 0,
  alpha = 100
)

w_alias <- obj$sparsify(max_iter = 10, criterion = "alias", reweight_fun = "entropy", recenter = TRUE, show_progress = FALSE)
dotest_tol("7.2.1", sum(w_alias > .Machine$double.eps^(1/4)), 8L)
dotest_ineq("7.2.2", obj$corr(w_alias), obj$orthogonal)

## 8.x other methods

obj <- lm_design$new(
  formula = ~(x1 + x2 + x3)^2 - 1,
  data = data
)

### crit
w_d <- obj$optimize(criterion = "D")
dotest_tol("8.1.1", obj$crit(criterion = "D"), 1.0)
dotest_tol("8.1.2", obj$crit(w_d, criterion = "A"), 0.1666667)
dotest_tol("8.1.3", obj$crit(w_d, criterion = "I"), 0.3)
dotest_tol("8.1.4", obj$crit(w_d, criterion = "G"), 0.1666677)
dotest_tol("8.1.5", obj$crit(w_d, criterion = "alias"), 1.0)

### sep
ses <- obj$sep(m = 10, sigma2 = 1)
dotest_tol("8.2.1", range(ses), c(0, 0.7745971))
dotest("8.2.2", length(ses), 27L)

ses1 <- obj$sep(newdata = obj$data, design_weights = w_d, sigma2 = 1)
dotest_tol("8.2.3", ses1, sqrt(10) * ses)

### vcov
w_round <- obj$round(m = 10, method = "optimal", seed = 1, show_progress = FALSE)
dotest_tol("8.3.1", obj$vcov(m = 10, sigma2 = 1), diag(0.1, 6))
dotest_tol("8.3.2", diag(obj$vcov(design_weights = w_round, sigma2 = 1)), rep(0.1071429, 6))
dotest_tol("8.3.3", sum(obj$vcov(design_weights = w_round, sigma2 = 1)), 0.7142857)

### corr
dotest_tol("8.4.1", obj$corr(), diag(6))
dotest_tol("8.4.2", sum(abs(obj$corr(design_weights = w_round, center = FALSE))), 8.4)
dotest_tol("8.4.3", sum(abs(obj$corr(design_weights = w_round, center = TRUE))), 8.2)

cat(sprintf("End of unit_tests_lm_design.R [elapsed: %.2fs]\n", (proc.time() - start)[3]))
