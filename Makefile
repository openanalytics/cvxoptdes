PKGNAME=cvxoptdes
PKGVERS=$(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)

.PHONY: doc vignette vignette-mkl readme manual globals test test-mkl covr build install install-mkl check check-valgrind \
    packamon docker docker-check docker-test docker-asan docker-rchk pkgdocs-build pkgdocs-hugo pkgdocs-hugo-serve pkgdocs-server

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
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz --no-manual --ignore-vignettes

check-valgrind: build
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz --use-valgrind

packamon:
	Rscript -e 'packamon::writeDockerfile("$(CURDIR)", dockerFilePath="Dockerfile", template="template.Dockerfile", overwrite=TRUE)'

docker:
	docker build -f Dockerfile --tag $(PKGNAME):$(PKGVERS) .

docker-check:
	docker run -it --rm $(PKGNAME):$(PKGVERS) R CMD check $(PKGNAME)_$(PKGVERS).tar.gz --no-manual

docker-test:
	docker run -it --rm $(PKGNAME):$(PKGVERS) Rscript -e "lapply(list.files(system.file(\"unit_tests\", package = \"$(PKGNAME)\"), full.names=TRUE), source, local=TRUE)"

docker-manual:
	docker run -it --user $(shell id -u):$(shell id -g) --rm -v $(CURDIR):/$(PKGNAME) -w /$(PKGNAME) rd2pdf:latest R CMD Rd2pdf . --force --no-preview

docker-rchk: build-no-vignettes
	$(RM) -r $(PKGNAME).Rcheck/ && mkdir -p $(PKGNAME).Rcheck && \
	cp $(PKGNAME)_$(PKGVERS).tar.gz $(PKGNAME).Rcheck/. && \
	docker run --rm -v $(CURDIR)/$(PKGNAME).Rcheck:/rchk/packages cvxoptdes-rchk:latest \
	/rchk/packages/$(PKGNAME)_$(PKGVERS).tar.gz
	
docker-asan:
	docker run --rm -v $(CURDIR):/$(PKGNAME) -w /$(PKGNAME) cvxoptdes-asan:latest bash -c \
	"R CMD INSTALL --preclean . && \
	Rscript --no-save tests/test_$(PKGNAME).R"
		
pkgdocs-build:
	Rscript --no-save inst/build/build_docs.R $(CURDIR) && \
	cd docs && \
	hugo mod get github.com/google/docsy/theme@v0.16.0 && hugo mod tidy && \
	npm install --save-dev @docsy/theme

pkgdocs-hugo:
	hugo build --source $(CURDIR)/docs

pkgdocs-hugo-docker:
	docker run --rm -it --user $(shell id -u):$(shell id -g) --entrypoint sh -w /src -v $(CURDIR)/docs:/src \
	floryn90/hugo:ext-alpine -c "npm ci --cache /tmp/.npm-cache && hugo --gc --minify"

pkgdocs-serve:
	python3 -m http.server --directory $(CURDIR)/docs/public 1313

pkgdocs-hugo-serve: pkgdocs-hugo
	python3 -m http.server --directory $(CURDIR)/docs/public 1313

clean:
	$(RM) -r $(PKGNAME).Rcheck/
	$(RM) -f *.gcno
	$(RM) -rf covr
