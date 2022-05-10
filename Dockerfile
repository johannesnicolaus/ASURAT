FROM rocker/r-base:4.1.3

# set work directory
WORKDIR /opt

# set argument for git access token
ARG GIT_ACCESS_TOKEN
ENV GITHUB_PAT ${GIT_ACCESS_TOKEN}

# Install libraries
RUN apt update && \
    apt install -y libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libfftw3-dev \
    libtiff-dev \
    libgsl-dev \
    git \
    procps

# set git credentials
RUN git config --global credential.helper \
  '!f() { echo username=author; echo "password=${GIT_ACCESS_TOKEN}"; };f'

# Install renv
ENV RENV_VERSION 0.15.4
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))" && \
  R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

# Add ENV variables for Renv cache (not symlink)
RUN mkdir /opt/renvcache && \
 echo RENV_PATHS_CACHE='/opt/renvcache' >> $(R RHOME)/etc/Renviron.site

# set the renv project directory as the renv default dir
ENV RENV_DIR='/opt/ASURAT'

# clone repo, install packages and move symlinks
RUN git clone https://github.com/keita-iida/ASURAT.git && \
  cd ASURAT && \
  R -e 'renv::restore()' && \
  R CMD INSTALL . && \
  chmod +x /opt/ASURAT/inst/*.R && \
  mv $(cd $RENV_DIR && Rscript -e "cat(paste0(renv::paths\$library()))") /opt/Rlibsymlinks && \
  echo "R_LIBS=/opt/Rlibsymlinks" >> $(R RHOME)/etc/Renviron.site

# set environment variables for R scripts within package
ENV PATH=/opt/ASURAT/inst:$PATH