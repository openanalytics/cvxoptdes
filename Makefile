PKGNAME=cvxoptdes
PKGVERS=$(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)

.PHONY: doc vignette vignette-mkl readme manual globals test test-mkl covr build build-no-vignettes \
    install install-mkl check check-cran check-no-vignettes packamon docker docker-check docker-test docker-manual

all: doc check clean

doc:
	Rscript -e "roxygen2::roxygenize(\".\")"

vignette:
	Rscript --no-save inst/build/build_vignettes.R $(VIGNETTES)

vignette-mkl:
	Rscript-mkl --no-save inst/build/build_vignettes.R $(VIGNETTES)

readme:
	Rscript -e "rmarkdown::render('README.Rmd')"

manual: doc
	Rscript -e "devtools::build_manual(\".\")"

globals:
	Rscript -e "checkglobals::checkglobals(pkg=\".\")"

test:
	Rscript --vanilla tests/test_$(PKGNAME).R

test-mkl:
	Rscript-mkl --no-save tests/test_$(PKGNAME).R

covr:
	Rscript -e "covr::report(covr::package_coverage(path=\".\"), file = \"covr/coverage.html\", browse = TRUE)"

build:
	R CMD build .

build-no-vignettes:
	R CMD build . --no-build-vignettes

install:
	R CMD INSTALL --preclean .

install-mkl:
	. /opt/intel/oneapi/setvars.sh && \
	R-mkl CMD INSTALL --preclean .

check: build
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz

check-cran: build
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz --as-cran

check-no-vignettes: build-no-vignettes
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz

packamon:
	Rscript -e 'packamon::writeDockerfile("$(PWD)", dockerFilePath="Dockerfile", template="template.Dockerfile", overwrite=TRUE)'

docker:
	docker build -f Dockerfile --tag $(PKGNAME):$(PKGVERS) .

docker-check:
	docker run -it --rm $(PKGNAME):$(PKGVERS) R CMD check $(PKGNAME)_$(PKGVERS).tar.gz --no-manual

docker-test:
	docker run -it --rm $(PKGNAME):$(PKGVERS) Rscript -e "lapply(list.files(system.file(\"unit_tests\", package = \"$(PKGNAME)\"), full.names=TRUE), source, local=TRUE)"

docker-manual:
	docker run -it --user $(shell id -u):$(shell id -g) --rm -v $(PWD):/$(PKGNAME) -w /$(PKGNAME) rd2pdf:latest R CMD Rd2pdf . --force --no-preview

clean:
	$(RM) -r $(PKGNAME).Rcheck/
	$(RM) -f *.gcno
	$(RM) -rf covr
