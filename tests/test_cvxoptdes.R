## source R-scripts in inst/unit_tests

## async evaluation
# options(warn = -1)
# mirai::daemons(n = 1, seed = 1)

## test scripts
source(system.file("unit_tests", "unit_tests_lm_design.R", package = "cvxoptdes"), local = new.env())
source(system.file("unit_tests", "unit_tests_nlm_design.R", package = "cvxoptdes"), local = new.env())
source(system.file("unit_tests", "unit_tests_lme_design.R", package = "cvxoptdes"), local = new.env())
source(system.file("unit_tests", "unit_tests_glm_gam_design.R", package = "cvxoptdes"), local = new.env())

# mirai::daemons(n = 0)
