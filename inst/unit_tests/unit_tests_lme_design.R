require(cvxoptdes)

start <- proc.time()
cat("Start of unit_tests_lme_design.R\n")

## helper functions
dotest <- function(itest, observed, expected) {
	if(!identical(observed, expected)) stop(sprintf("Unit test lme-%s failed", itest), call. = FALSE)
}
dotest_tol <- function(itest, observed, expected, tol = (.Machine$double.eps)^0.25) {
	if(any(abs(observed - expected) > tol)) stop(sprintf("Unit test lme-%s failed", itest), call. = FALSE)
}
dotest_ineq <- function(itest, left, right, tol = (.Machine$double.eps)^0.25) {
	if(any(left > right + tol)) stop(sprintf("Unit test lme-%s failed", itest), call. = FALSE)
}

## 1.x optimize, round
data <- data.frame(
  treatment = c(-1, -1, 0, 0, 1, 1, -1, 0, 0, 1, -1, 1),
  block = factor(rep(c("B1", "B2", "B3", "B4", "B5", "B6"), each = 2))
)

## w/ random effects
obj <- lme_design$new(
		formula = ~ treatment + I(treatment^2) ~ block,
		data = data,
		eta = 10
)

## w/o random effects
obj0 <- lme_design$new(
		formula = ~ treatment + I(treatment^2) ~ block,
		data = data,
		eta = 0
)

### d-optimal
w_d0 <- structure(rep(c(1/6, 0), each = 6), criterion = list(crit = "D", crit.value = 0.2814851))

w_d <- obj$optimize(criterion = "D", max_iter = 5, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_ineq("1.1.1", 0.3107272, d_val)

w_d <- obj0$optimize(criterion = "D")
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("1.1.2", d_val, 0.5291337)

w_d <- obj0$round(m = 12, method = "efficient", seed = 1, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("1.1.3", d_val, 0.5291337)

w_d <- obj$round(m = 12, method = "optimal", seed = 1, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("1.1.4", d_val, 0.2978544)

w_d <- obj0$round(m = 12, method = "optimal", design_weights = w_d0, seed = 1, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("1.1.5", d_val, 0.5291337)

### a-optimal
w_a0 <- structure(w_d0, criterion = list(crit = "A", crit.value = 0.02564103))

w_a <- obj$optimize(criterion = "A", max_iter = 10, show_progress = FALSE)
a_val <- attr(w_a, "criterion")$crit.value
dotest_ineq("1.2.1", a_val, 0.0834)

w_a <- obj0$optimize(criterion = "A")
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("1.2.2", a_val, 0.125)

w_a <- obj$round(m = 12, method = "optimal", seed = 1, show_progress = FALSE)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("1.2.3", a_val, 0.07804441)

w_a <- obj0$round(m = 12, method = "optimal", design_weights = w_a0, seed = 1, show_progress = FALSE)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("1.2.4", a_val, 0.125)

### i-optimal
w_i0 <- structure(w_d0, criterion = list(crit = "I", crit.value = 0.07692308))

w_i <- obj$optimize(criterion = "I", max_iter = 5, show_progress = FALSE)
i_val <- attr(w_i, "criterion")$crit.value
dotest_tol("1.3.1", i_val, 0.1814919)

w_i <- obj0$optimize(criterion = "I")
i_val <- attr(w_i, "criterion")$crit.value
dotest_tol("1.3.2", i_val, 0.333333)

w_i <- obj$round(m = 12, method = "optimal", seed = 1, show_progress = FALSE)
i_val <- attr(w_i, "criterion")$crit.value
dotest_tol("1.3.3", i_val, 0.1773161)

w_i <- obj0$round(m = 12, method = "optimal", design_weights = w_i0, seed = 1, show_progress = FALSE)
i_val <- attr(w_i, "criterion")$crit.value
dotest_tol("1.3.4", i_val, 0.333333)

### g-optimal
w_g0 <- structure(w_d0, criterion = list(crit = "G", crit.value = 0.07692308))

w_g <- obj$optimize(criterion = "G", max_iter = 5, show_progress = FALSE)
g_val <- attr(w_g, "criterion")$crit.value
dotest_tol("1.4.1", g_val, 0.1814919)

w_g <- obj0$optimize(criterion = "G")
g_val <- attr(w_g, "criterion")$crit.value
dotest_tol("1.4.2", g_val, 0.333333)

w_g <- obj$round(m = 12, method = "optimal", seed = 1, show_progress = FALSE)
g_val <- attr(w_g, "criterion")$crit.value
dotest_tol("1.4.3", g_val, 0.1681034)

w_g <- obj0$round(m = 12, method = "optimal", design_weights = w_g0, seed = 1, n_repeats = 1)
g_val <- attr(w_g, "criterion")$crit.value
dotest_tol("1.4.4", g_val, 0.333333)

### alias-optimal
w_alias0 <-  structure(w_d0, criterion = list(crit = "alias", crit.value = 0.7278345))

w_alias <- obj$optimize(criterion = "alias")
alias_val <- attr(w_alias, "criterion")$crit.value
dotest_tol("1.5.1", alias_val, 0.699103)

w_alias <- obj0$optimize(criterion = "alias")
alias_val <- attr(w_alias, "criterion")$crit.value
dotest_tol("1.5.2", alias_val, 0.6806489)

w_alias <- obj$round(m = 12, method = "optimal", seed = 1, show_progress = FALSE)
alias_val <- attr(w_alias, "criterion")$crit.value
dotest_tol("1.5.3", alias_val, 0.7524062)

w_alias <- obj0$round(m = 12, method = "optimal", seed = 1, show_progress = FALSE)
alias_val <- attr(w_alias, "criterion")$crit.value
dotest_tol("1.5.4", alias_val, 0.7278345)

### ridge penalty
obj$update(lambda = 1e-2)

w_d <- obj$optimize(criterion = "D", max_iter = 5, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("1.6.1", d_val, 0.3240983)

w_d <- obj$round(m = 5, method = "optimal", seed = 1, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("1.6.2", d_val, 0.2565065)

## 2.x split-plot

data <- expand.grid(
  A = c(0, 1),
  B = c(0, 1),
  C = c(0, 1),
  week = factor(1:3)
)

obj <- lme_design$new(
  formula = ~ (A + B + C)^2 ~ week,
  data = data,
  eta = 10
)

w_d <- obj$optimize(criterion = "D", show_progress = FALSE)

w_d <- obj$round_split_plot(m = 4, seed = 1, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("2.1.1", d_val, 0.1297256)
dotest_tol("2.1.2", attr(w_d, "strata") %% 4, rep(0, 3))

w_d <- obj$round_split_plot(m = 4, hard_to_change = ~A, seed = 1, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("2.1.3", d_val, 0.1052086)
dotest_tol("2.1.4", attr(w_d, "strata") %% 4, rep(0, 2))

## 3.x augment method

data <- data.frame(
  treatment = c(-1, -1, 0, 0, 1, 1, -1, 0, 0, 1, -1, 1),
  block = factor(rep(c("B1", "B2", "B3", "B4", "B5", "B6"), each = 2))
)

## w/ random effects
obj <- lme_design$new(
  formula = ~ treatment + I(treatment^2) ~ block,
  data = data,
  eta = 10
)

## w/o random effects
obj0 <- lme_design$new(
  formula = ~ treatment + I(treatment^2) ~ block,
  data = data,
  eta = 0
)

### d-optimal
w_d01 <- structure(rep(c(0, 1), each = 6), criterion = list(crit = "D", crit.value = 0.2814851))

w_d <- obj$augment(design_weights = w_d01, max_iter = 5, gamma = 0.5, criterion = "D", show_progress = FALSE)
dotest_ineq("3.1.1", 0.3086112,  obj$crit(w_d + w_d01 / 6, "D"))

w_d <- obj0$augment(design_weights = w_d01, gamma = 0.5, criterion = "D")
dotest_tol("3.1.2", obj0$crit(w_d + w_d01 / 6, "D"), 0.5291337)

w_d <- obj$round(m = 6, method = "optimal", seed = 1, augment_design = w_d01, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("3.1.3", d_val, 0.2814851)

w_d <- obj0$round(m = 6, method = "optimal", seed = 1, augment_design = w_d01, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("3.1.4", d_val, 0.5291337)

### i-optimal
w_i01 <- structure(w_d01, criterion = list(crit = "I", crit.value = 0.1468531))

data_w <- transform(data, weights = w_i01)
obj$update(data = data_w)
obj0$update(data = data_w)

w_i <- obj$augment(design_weights = ~weights, max_iter = 10, gamma = 0.5, criterion = "I", show_progress = FALSE)
dotest_ineq("3.2.1", 0.1721879, obj$crit(w_i + data_w$weights / 6, "I"))

w_i <- obj0$augment(design_weights = ~weights, gamma = 0.5, criterion = "I")
dotest_tol("3.2.2", obj0$crit(w_i + data_w$weights / 6, "I"), 0.333333)

w_i <- obj$round(m = 6, method = "optimal", seed = 1, augment_design = w_i01, show_progress = FALSE)
i_val <- attr(w_i, "criterion")$crit.value
dotest_tol("3.2.3", i_val, 0.1580991)

w_i <- obj0$round(m = 6, method = "optimal", seed = 1, augment_design = w_i01, show_progress = FALSE)
i_val <- attr(w_i, "criterion")$crit.value
dotest_tol("3.2.4", i_val, 0.3333333)

### a-optimal
w_a01 <- structure(w_d01, criterion = list(crit = "A", crit.value = 0.07023411))

w_a <- obj$augment(design_weights = w_a01, max_iter = 10, gamma = 0.5, criterion = "A", show_progress = FALSE)
dotest_ineq("3.3.1", 0.08159331, obj$crit(w_a + w_a01 / 6, "A"))

w_a <- obj0$augment(design_weights = w_a01, gamma = 0.5, criterion = "A")
dotest_tol("3.3.2", obj0$crit(w_a + w_a01 / 6, "A"), 0.125)

w_a <- obj$round(m = 6, method = "optimal", seed = 1, augment_design = w_a01, show_progress = FALSE)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("3.3.3", a_val, 0.0794177)

w_a <- obj0$round(m = 6, method = "optimal", seed = 1, augment_design = w_a01, show_progress = FALSE)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("3.3.4", a_val, 0.125)

### g-optimal
w_g01 <- structure(w_d01, criterion = list(crit = "G", crit.value = 0.1468531))

w_g <- obj$augment(design_weights = w_g01, max_iter = 5, gamma = 0.5, criterion = "G", show_progress = FALSE)
dotest_tol("3.4.1", obj$crit(w_g + w_g01 / 6, "G"),  0.1752357)

w_g <- obj$round(m = 6, method = "optimal", seed = 1, augment_design = w_g01, show_progress = FALSE)
g_val <- attr(w_g, "criterion")$crit.value
dotest_tol("3.4.2", g_val, 0.151211)

### alias-optimal
w_alias01 <- structure(w_d01, criterion = list(crit = "alias", crit.value = 0.7989924))

w_alias <- obj$augment(design_weights = w_alias01, gamma = 0.5, criterion = "alias")
dotest_tol("3.5.1", obj$crit(w_alias + w_alias01 / 6, "alias"), 0.7281696)

## 4.x cost penalty

data_cost <- transform(data, cost = rep(c(0, 1), times = 6))

obj <- lme_design$new(
  formula = ~ treatment + I(treatment^2) ~ block,
  data = data_cost,
  eta = 10,
  cost = ~cost,
  alpha = 1
)

obj0 <- lme_design$new(
  formula = ~ treatment + I(treatment^2) ~ block,
  data = data_cost,
  eta = 0,
  cost = ~cost,
  alpha = 1
)

### d-optimal
w_d <- obj$optimize(criterion = "D", max_iter = 5, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_ineq("4.1.1", 0.2299328, d_val)

w_d <- obj0$optimize(criterion = "D")
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("4.1.2", d_val, 0.5291335)

w_d <- obj$round(m = 12, method = "optimal", seed = 1, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("4.1.3", d_val, 0.2294685)

w_d <- obj0$round(m = 12, method = "optimal", seed = 1, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("4.1.4", d_val, 0.5291335)

### a-optimal
obj$update(cost = data_cost$cost)
obj0$update(cost = data_cost$cost)

w_a <- obj$optimize(criterion = "A", max_iter = 2, show_progress = FALSE)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("4.2.1", a_val, 0.04054074)

w_a <- obj0$optimize(criterion = "A")
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("4.2.2", a_val, 0.125)

w_a <- obj$round(m = 12, method = "optimal", seed = 1, show_progress = FALSE)
a_val <- attr(w_a, "criterion")$crit.value
dotest_tol("4.2.3", a_val, -0.04127726)

### i-optimal
w_i <- obj$optimize(criterion = "I", max_iter = 2, show_progress = FALSE)
i_val <- attr(w_i, "criterion")$crit.value
dotest_tol("4.3.1", i_val, 0.1097556)

w_i <- obj0$optimize(criterion = "I")
i_val <- attr(w_i, "criterion")$crit.value
dotest_tol("4.3.2", i_val, 0.333333)

w_i <- obj$round(m = 12, method = "optimal", seed = 1, show_progress = FALSE)
i_val <- attr(w_i, "criterion")$crit.value
dotest_tol("4.3.3", i_val, 0.1090175)

### alias-optimal
w_alias <- obj$optimize(criterion = "alias")
alias_val <- attr(w_alias, "criterion")$crit.value
dotest_tol("4.4.1", alias_val, 0.459041)

w_alias <- obj0$optimize(criterion = "alias")
alias_val <- attr(w_alias, "criterion")$crit.value
dotest_tol("4.4.2", alias_val, 0.6806492)

### augment
w_d <- obj$augment(design_weights = w_d01, max_iter = 5, gamma = 0.5, criterion = "D", show_progress = FALSE)
dotest_ineq("4.5.1", 0.2214764, obj$crit(w_d + w_d01 / 6, "D"))

w_d <- obj0$augment(design_weights = w_d01, gamma = 0.5, criterion = "D")
dotest_tol("4.5.2", obj0$crit(w_d + w_d01 / 6, "D"), 0.4120897)

w_d <- obj$round(m = 6, method = "optimal", seed = 1, augment_design = w_d01, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("4.5.3", d_val, 0.2180303)

w_d <- obj0$round(m = 6, method = "optimal", seed = 1, augment_design = w_d01, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("4.5.4", d_val, 0.4120897)

### weights
obj$update(
  cost = NULL,
  weights = ~ 1 / (100 * pmax(cost, 0.01))
)

obj0$update(
  cost = NULL,
  weights = ~ 1 / (100 * pmax(cost, 0.01))
)

w_d <- obj$optimize(criterion = "D", max_iter = 5, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_ineq("4.6.1", 0.1824448, d_val)

w_d <- obj0$optimize(criterion = "D")
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("4.6.2", d_val, 0.5291337)

w_d <- obj$round(m = 12, method = "optimal", seed = 1, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("4.6.3", d_val, 0.1631391)

w_d <- obj0$round(m = 12, method = "optimal", seed = 1, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("4.6.4", d_val, 0.5291335)

w_d <- obj$augment(design_weights = w_d01, max_iter = 10, gamma = 0.5, criterion = "D", show_progress = FALSE)
dotest_ineq("4.6.5", obj$crit(w_d + w_d01 / 6, "D"), 0.1862355)

w_d <- obj0$augment(design_weights = w_d01, gamma = 0.5, criterion = "D")
dotest_tol("4.6.6", obj0$crit(w_d + w_d01 / 6, "D"), 0.3981732)

w_d <- obj$round(m = 6, method = "optimal", seed = 1, augment_design = w_d01, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("4.6.7", d_val, 0.1774354)

w_d <- obj0$round(m = 6, method = "optimal", seed = 1, augment_design = w_d01, show_progress = FALSE)
d_val <- attr(w_d, "criterion")$crit.value
dotest_tol("4.6.8", d_val, 0.3981716)

## 5.x upper limit
obj <- lme_design$new(
  formula = ~ treatment + I(treatment^2) ~ block,
  data = data,
  eta = 10,
  upper = 0.1
)

### d-optimal
w_d <- obj$optimize(criterion = "D", max_iter = 5, show_progress = FALSE)
dotest_tol("5.1.1", max(w_d), 0.1)

w_d <- obj$round(m = 12, method = "optimal", seed = 1, show_progress= FALSE)
dotest_tol("5.1.2", max(w_d), 1)

w_d <- obj$round(m = 12, method = "optimal", seed = 1, replicates = FALSE, show_progress = FALSE)
dotest_tol("5.1.3", max(w_d), 1)

w_d <- obj$round(m = 6, method = "optimal", seed = 1, replicates = 1, augment_design = w_d01, show_progress = FALSE)
dotest_tol("5.1.4", max(w_d), 1)

w_a <- obj$optimize(criterion = "A", max_iter = 2, show_progress = FALSE)
dotest_tol("5.2.1", max(w_a), 0.1)

w_i <- obj$optimize(criterion = "I", max_iter = 2, show_progress = FALSE)
dotest_tol("5.2.2", max(w_i), 0.1)

w_g <- obj$optimize(criterion = "G", max_iter = 2, show_progress = FALSE)
dotest_tol("5.2.3", max(w_g), 0.1)

w_alias <- obj$optimize(criterion = "alias")
dotest_tol("5.2.4", max(w_alias), 0.1)

## 6.x stratified design
data <- expand.grid(
  A = c(0, 1),
  B = c(0, 1),
  C = c(0, 1),
  week = factor(1:3)
)

obj <- lme_design$new(
  formula = ~ (A + B + C)^2 ~ week,
  data = data,
  eta = 10
)
obj0 <- lme_design$new(
  formula = ~ (A + B + C)^2 ~ week,
  data = data,
  eta = 0
)

obj$stratify(
		strata = ~week,
		proportions = c(2/3, 0, 1/3)
)
obj0$stratify(
		strata = ~week,
		proportions = c(2/3, 0, 1/3)
)

### d-optimal
w_d <- obj$optimize(criterion = "D", max_iter = 2, show_progress = FALSE)
dotest_tol("6.1.1", tapply(w_d, obj$strata$strata, sum), c("-1" = 2/3, "0" = 0, "1" = 1/3))

w_d <- obj0$optimize(criterion = "D")
dotest_tol("6.1.2", tapply(w_d, obj$strata$strata, sum), c("-1" = 2/3, "0" = 0, "1" = 1/3))

w_d <- obj$round(m = 9, method = "optimal", seed = 1, show_progress = FALSE)
dotest_tol("6.1.3", attr(w_d, "strata"), c("-1" = 6, "0" = 0, "1" = 3))

w_d <- obj0$round(m = 9, method = "optimal", seed = 1, show_progress = FALSE)
dotest_tol("6.1.4", attr(w_d, "strata"), c("-1" = 6, "0" = 0, "1" = 3))

w_a <- obj$optimize(criterion = "A", max_iter = 2, show_progress = FALSE)
dotest_tol("6.2.1", tapply(w_a, obj$strata$strata, sum), c("-1" = 2/3, "0" = 0, "1" = 1/3))

w_i <- obj$optimize(criterion = "I", max_iter = 2, show_progress = FALSE)
dotest_tol("6.2.2", tapply(w_i, obj$strata$strata, sum), c("-1" = 2/3, "0" = 0, "1" = 1/3))

w_g <- obj$optimize(criterion = "G", max_iter = 2, show_progress = FALSE)
dotest_tol("6.2.3", tapply(w_g, obj$strata$strata, sum), c("-1" = 2/3, "0" = 0, "1" = 1/3))

w_alias <- obj$optimize(criterion = "alias")
dotest_tol("6.2.4", tapply(w_alias, obj$strata$strata, sum), c("-1" = 2/3, "0" = 0, "1" = 1/3))

## 7.x other methods

data <- data.frame(
  treatment = c(-1, -1, 0, 0, 1, 1, -1, 0, 0, 1, -1, 1),
  block = factor(rep(c("B1", "B2", "B3", "B4", "B5", "B6"), each = 2))
)

obj <- lme_design$new(
  formula = ~ treatment + I(treatment^2) ~ block,
  data = data,
  eta = 10
)

### crit
w_d <- obj$optimize(criterion = "D", max_iter = 5, show_progress = FALSE)
dotest_ineq("7.1.1", 0.3086754, obj$crit(criterion = "D"))
dotest_ineq("7.1.2", 0.07416199, obj$crit(w_d, criterion = "A"))
dotest_ineq("7.1.3", 0.1799158, obj$crit(w_d, criterion = "I"))
dotest_ineq("7.1.4", 0.1787619, obj$crit(w_d, criterion = "G"))
dotest_ineq("7.1.5", 0.76848, obj$crit(w_d, criterion = "alias"))

### sep
ses <- obj$sep()
dotest_ineq("8.2.1", 2.351431, min(ses))
dotest("8.2.2", length(ses), 12L)

ses1 <- obj$sep(m = 1, newdata = obj$data, design_weights = w_d, sigma2 = 1)
dotest_tol("8.2.3", ses1, ses)

### vcov
w_round <- obj$round(m = 12, method = "optimal", seed = 1, show_progress = FALSE)
dotest_tol("8.3.1", diag(obj$vcov(design_weights = w_round, sigma2 = 1)), c(2.4281305, 0.2363732, 0.4929453))
dotest_tol("8.3.2", sum(abs(obj$vcov(design_weights = w_round, sigma2 = 1))), 3.939636)

### corr
dotest_tol("8.4.1", obj$corr(design_weights = w_round, center = FALSE), diag(2))
dotest_tol("8.4.2", obj$corr(design_weights = w_round, center = TRUE), diag(2))

### subset
design <- obj$subset(w_round)
dotest("8.5.1", dim(design), c(12L, 2L))
dotest("8.5.2", colnames(design), c("treatment", "block"))

cat(sprintf("End of unit_tests_lme_design.R [elapsed: %.2fs]\n", (proc.time() - start)[3]))

