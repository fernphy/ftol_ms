FROM rocker/verse:4.1.3

ARG DEBIAN_FRONTEND=noninteractive

############################
### Install APT packages ###
############################

# libpoppler-cpp-dev for pdftools

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    libpoppler-cpp-dev \
  && apt-get clean

####################################
### Install R packages with renv ###
####################################

# Create directory for renv project library
RUN mkdir renv

# Modify Rprofile.site so renv uses /renv for project library
RUN echo 'Sys.setenv(RENV_PATHS_LIBRARY = "/renv")' >> /usr/local/lib/R/etc/Rprofile.site

# Initialize a 'dummy' project and restore the renv library.
# Since the library path is specified as above, the library will be restored to /renv
RUN mkdir /tmp/project

COPY ./renv.lock /tmp/project

WORKDIR /tmp/project

# Don't use cache (the symlinks won't work from Rstudio server)
RUN Rscript -e 'install.packages("renv"); renv::consent(provided = TRUE); renv::settings$use.cache(FALSE); renv::init(bare = TRUE); renv::restore()'

#############################
### Other custom software ###
#############################

ENV APPS_HOME=/apps
RUN mkdir $APPS_HOME
WORKDIR $APPS_HOME

### gnparser ###
ENV APP_NAME=gnparser
ENV VERSION=1.4.0
ENV DEST=$APPS_HOME/$APP_NAME/$VERSION
RUN wget https://github.com/gnames/gnparser/releases/download/v$VERSION/gnparser-v$VERSION-linux.tar.gz \
  && tar xf $APP_NAME-v$VERSION-linux.tar.gz \
  && rm $APP_NAME-v$VERSION-linux.tar.gz \
  && mv "$APP_NAME" /usr/local/bin/

WORKDIR /home/rstudio/
