#
# dockerfile for qbound paper
#
# docker build --rm -t shabbychef/qbound .
#
# docker run -it --rm --volume $(pwd):/srv:rw fromo-crancheck
#
# Created: 2018.05.08
# Copyright: Steven E. Pav, 2018
# Author: Steven E. Pav
# Comments: Steven E. Pav

#####################################################
# preamble# FOLDUP
FROM shabbychef/crancheck
MAINTAINER Steven E. Pav, shabbychef@gmail.com
# UNFOLD

ENV DOCKERFILE_REFRESHED_AT 2018.05.08
# see http://crosbymichael.com/dockerfile-best-practices.html
#RUN echo "deb http://archive.ubuntu.com/ubuntu precise main universe" > /etc/apt/sources.list

RUN (apt-get clean -y ; \
 apt-get update -y -qq; \
 apt-get update --fix-missing );

#RUN (apt-get dist-upgrade -y ; \
#apt-get update -qq ; \
RUN (DEBIAN_FRONTEND=noninteractive DEBCONF_NONINTERACTIVE_SEEN=true apt-get install -q -y --no-install-recommends \ 
  libgs9 texlive-base texlive-binaries \
	libcupsimage2 libcups2 curl wget \
	qpdf pandoc ghostscript \
	texlive-latex-extra texlive-latex-base texlive-fonts-recommended texlive-fonts-extra latexmk ; \
	sync ; \
	mkdir -p /usr/local/lib/R/site-library ; \
	chmod -R 777 /usr/local/lib/R/site-library ; \
	sync ; \
	/usr/local/bin/install2.r knitr devtools doFuture doRNG dplyr ggplot2 hypergeo knitr LambertW microbenchmark quantmod SharpeR tidyr xtable ; \
	/usr/local/bin/r -l 'devtools' -e 'options(unzip="internal");install_github("shabbychef/aqfb_data");' ) 

ADD Makefile /srv/ 
ADD qbound.Rnw /srv/ 
ADD qbound.R /srv/ 
ADD *.bib /srv/ 
ADD *.sty /srv/ 
ADD sp100lr.rda /srv/ 

RUN mkdir -p /srv/cache /srv/figure ;

#RUN groupadd -g 1000 spav && useradd -g spav -u 1000 spav;
#USER spav

#####################################################
# entry and cmd# FOLDUP
# these are the default, but remind you that you might want to use /usr/bin/R instead?
# always use array syntax:
# fix these.
#ENTRYPOINT ["/usr/local/bin/R","CMD","check","--as-cran","--output=/tmp"]
ENTRYPOINT ["/bin/bash"]

# ENTRYPOINT and CMD are better together:
CMD ["-i"]
# UNFOLD

#for vim modeline: (do not edit)
# vim:nu:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=Dockerfile:ft=Dockerfile:fo=croql
