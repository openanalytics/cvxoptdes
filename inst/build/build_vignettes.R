## build package vignettes
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1) {
	rmd_orig <- list.files(path = "inst/vignettes_orig", pattern = "\\.Rmd$", full.names = TRUE)
} else {
	rmd_orig <- file.path("inst/vignettes_orig", paste0(strsplit(args[[1]], ",\\s+")[[1]], ".Rmd"))
}

for(file in rmd_orig) {
	knitr::knit(
		input = file,
		output = file.path("vignettes", basename(file))
	)
}
