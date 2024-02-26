FROM r-base 
RUN apt-get update && apt-get install -y --no-install-recommends \
    sudo \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    && rm -rf /var/lib/apt/lists/*
RUN install.r shiny
RUN R -e "install.packages(c('devtools','ggbeeswarm','dygraphs','comprehenr','pals','viridis','ggplot2','tidyr','shiny','dplyr','htmlwidgets','digest','bit','flexdashboard','ggthemes','highcharter','scatterD3'),dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "tinytex::install_tinytex()"
#Run R -e "devtools::install_github("rstudio/crosstalk")"
#RUN R -e "devtools::install_github("jcheng5/d3scatter")"
#RUN R -e "devtools::install_github("hadley/shinySignals"") 
COPY . ~/Desktop/test/test1/test2/LoxCodeR2022
WORKDIR ./LoxCodeR2022/
RUN ls
RUN R CMD build LoxCodeR2022
