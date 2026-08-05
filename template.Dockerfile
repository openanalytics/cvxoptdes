#include packamon.disclaimer

#include packamon.from

LABEL maintainer="joris.chau@openanalytics.eu"

ARG PKGNAME="cvxoptdes"
ENV PKGNAME=${PKGNAME}

##include packamon.system-dependencies

# System libraries (incl. system requirements for R packages)
RUN apt-get update && apt-get install --no-install-recommends -y \
    libcurl4-openssl-dev \
    libicu-dev \
    libssl-dev \
    libxml2-dev \
    make \
	cmake \
    pandoc \
    && rm -rf /var/lib/apt/lists/*
	
#include packamon.r-repos

#include packamon.r-dependencies

RUN R -q -e "options(warn = 2); \
    install.packages('checkglobals')"

##include packamon.local-r-dependencies

##include packamon.runtime-settings

## install local package(s)
RUN mkdir -p /tmp
WORKDIR /tmp
COPY . /tmp/${PKGNAME}
RUN R -q -e "install.packages(pkgbuild::build('${PKGNAME}'), repos = NULL, dependencies = FALSE)"
