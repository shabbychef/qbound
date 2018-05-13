
######################
######################
# makefile generated 
# dude makefile 
# created by s.e.pav 
# $Id: Makefile 89 2006-01-26 20:59:08Z spav $
######################
######################

############### FLAGS ###############

# these don't work at times
# in that case set these by hand?

LATEX       := $(shell which latex)
BIBTEX      := $(shell which bibtex)
PDFLATEX    := $(shell which pdflatex)
MAKEINDEX   := $(shell which makeindex)
PAGER   		:= $(shell which less)
ASPELL  		:= $(shell which aspell)

RLIB         = /usr/lib64/R

TEXINPADD    = .:$(RLIB)/share/texmf/tex/latex

PRETEX       = TEXINPUTS=$(TEXINPADD):$$TEXINPUTS
PREBIB       = BSTINPUTS=$(TEXINPADD):$$BSTINPUTS \
               BIBINPUTS=$(TEXINPADD):$$BIBINPUTS 

PREIDX       = INDEXSTYLE=$(TEXINPADD):$$INDEXSTYLE

#undoes psfrag for pdf
UNPSFRAG		 = perl $(HOME)/sys/bin/unpsfrag.pl
#unroll commands
DETEXIFY		 = perl $(HOME)/sys/perl/detexify.pl

SCREEN_SIZE  = normal
#include	$(HOME)/sys/etc/.Makefile.local

PROJECT      = qbound
TEX_SOURCE   = $(PROJECT).tex
BIB_SOURCE   = $(PROJECT).bib
DVI_TARGET   = $(PROJECT).dvi
PDF_TARGET   = $(PROJECT).pdf
BBLS         = $(PROJECT).bbl

#SAVE
# tracked projects
PROJECTS     = $(PROJECT) 
#UNSAVE
# add on dependencies (subchapters of qbound)
R_DEPS 			 = 
TEX_EXTRAS   = SharpeR.sty SharpeR.bib 
# nonlocal dependencies
STY_FILES    = 

#aspell
ASPELL_FLAGS = 

# for running in docker this gets passed to the knit
# and controls speed of build. larger is faster.
RUNTIME_PARAM 		?= 1

GID 							?= $$UID
DOCKER_RUN_FLAGS 		= --user $$UID:$(GID)
DOCKER_ENV 				 = -e FOO_ENV='foo' -e RUNTIME_PARAM=$(RUNTIME_PARAM)

DOCKER 						?= $(shell which docker)
DOCKER_IMG 				 = .docker_img
DOCKER_NAME 			 = $(USER)/qbound

RESULTS_D 				 = output
DOWNSTREAM_D 			 = ../../qbound

############## DEFAULT ##############

default : all

############## MARKERS ##############

.PHONY   : 
.SUFFIXES: .tex .bib .dvi .ps .pdf .eps
.PRECIOUS: %.pdf 

############ BUILD RULES ############

.PHONY   : help 

# this will have to change b/c of inclusion file names...
help:  ## generate this help message
	@grep -h -P '^(([^\s]+\s+)*([^\s]+))\s*:.*?##\s*.*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'


.PHONY   : all

# an easy target
all : $(PROJECT).pdf  ## build the document by knitting source code

.PHONY   : doc

doc : $(RESULTS_D)/$(PROJECT).pdf  | $(RESULTS_D) ## build the document by knitting source code, for use in docker

%.tex : %.Rnw $(R_DEPS)
		Rscript -e 'require(knitr);knit("$<")'

%.R : %.Rnw
		Rscript -e 'require(knitr);knit("$<",tangle=TRUE)'

%.pdf : %.tex
	latexmk -f -bibtex -pdf -pdflatex="$(PDFLATEX)" -use-make $<
		
$(RESULTS_D) : 
	mkdir -p $@
	chmod 755 $@

$(RESULTS_D)/%.pdf : %.pdf
	cp $< $@
	echo "in docker, it runs as root, so chmod"
	chmod 777 $@

# tex extras
%.bbl : %.bib
	$(PREBIB) $(BIBTEX) $*

%.bbl : %.aux
	$(PREBIB) $(BIBTEX) $*

%.ind : %.idx
	$(PREIDX) $(MAKEINDEX) $*

# check a document
%.chk : %.dup %.spell

# check spelling
%.spell : %.tex
	$(ASPELL) $(ASPELL_FLAGS) --dont-tex-check-comments -t -l < $< | sort | uniq | $(PAGER)

# check duplicate words
%.dup : %.tex
	perl -an -F/\\s+/ -e 'BEGIN { $$last = q[]; $$line = 0; $$prevline = q[];}\
	$$line++;$$first = 1;\
	foreach $$word (@F) {\
	if ($$word eq $$last) {\
	if ($$first) { print qq[duplicate $$word, lines ],($$line-1),qq[-$$line:\n$$prevline$$_]; }\
	else { print qq[duplicate $$word, line $$line:\n$$_]; } }\
	$$last = $$word; $$first = 0; } \
	$$prevline = $$_;' < $< | $(PAGER)

############# CLEAN UP ##############

# clean up
%.clean : 
		-rm -f $*.aux $*.log $*.dvi $*.bbl $*.blg $*.toc $*.ilg $*.ind
		-rm -f $*.out $*.idx $*.lot $*.lof $*.brf $*.nav $*.snm
%.realclean : %.clean
		-rm -f $*.ps $*.pdf
		-rm -f $*-[0-9][0-9][0-9]*.eps $*-[0-9][0-9][0-9]*.pdf

############### RULES ###############


release.tex: qbound.tex
	perl -pe 's{figure/}{};' < $< > $@

release : release.tex
	mv release.tex qbound.tex
	@-echo "upload to arxiv"

# check it

spell: $(PROJECT).spell 

# clean up
clean: $(patsubst %,%.clean,$(PROJECTS))
realclean: $(patsubst %,%.realclean,$(PROJECTS))
		-rm -f Rplots.pdf

cleancache: 
		echo "killing knitr cache! ack!"
		-rm -rf cache

superclean: realclean cleancache

.tags :
	nice -n 18 ctags -f .tmp_tags --recurse --language-force=R --fields=+i `find . -regextype posix-egrep -regex '.*.R(nw)?'`;
	mv .tmp_tags $@

######################
# docker

.PRECIOUS: $(DOCKER_IMG)

docker_img : $(DOCKER_IMG) ## build the docker image

DOCKER_DEPS 			= Makefile $(PROJECT).R $(PROJECT).Rnw $(wildcard *.bib) $(wildcard *.sty) sp100lr.rda 

$(DOCKER_IMG) : Dockerfile $(DOCKER_DEPS)
	tar -czh . | $(DOCKER) build --rm -t $(DOCKER_NAME) -
	touch $@

.PHONY : docker_doc

docker_doc : $(DOCKER_IMG)  ## use docker to build the document
	$(DOCKER) run -it --rm $(DOCKER_ENV) -v $$(pwd)/$(RESULTS_D):/srv/$(RESULTS_D):rw --entrypoint="make" $(DOCKER_NAME) \
		"doc"

docker_R : $(DOCKER_IMG)  ## use docker to run R in this context.
	$(DOCKER) run -it --rm $(DOCKER_ENV) -v $$(pwd)/$(RESULTS_D):/srv/$(RESULTS_D):rw --entrypoint="R" $(DOCKER_NAME) 

docker_shell : $(DOCKER_IMG)  ## use docker to run bash in this context.
	$(DOCKER) run -it --rm $(DOCKER_ENV) -v $$(pwd)/$(RESULTS_D):/srv/$(RESULTS_D):rw --entrypoint="/bin/bash" $(DOCKER_NAME) "-i"

######################
# copy downstream

DOWNSTREAM_FILES 				 = $(foreach fname,README.md Dockerfile $(DOCKER_DEPS),$(DOWNSTREAM_D)/$(fname))

$(DOWNSTREAM_FILES) : $(DOWNSTREAM_D)/% : %
	cp $< $@

.PHONY : downstream_files

downstream_files : $(DOWNSTREAM_FILES)  ## copy files from this repo to the public, downstream repo.
		

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=149:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:tags=tags;:syn=make:ft=make:ai:si:cin:nu:fo=croqt:cino=p0t0c5(0:
