FROM rocker/geospatial:4
COPY . /workdir

RUN Rscript -e "remotes::install_github('ebird/ebirdst')"
RUN Rscript -e "install.packages(c('fields','rnaturalearth'), repos = 'http://cran.rstudio.com')"
RUN Rscript -e "remotes::install_github('ropensci/rnaturalearthhires')"
RUN Rscript -e "install.packages(c('styler','DT', 'viridis'), repos = 'http://cran.rstudio.com')"

WORKDIR /workdir
