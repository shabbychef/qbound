#!/bin/bash -e
#
# copied from https://raw.githubusercontent.com/Strauman/travis-latexbuild/master/travis-texbuild.sh
# cf 
# https://tex.stackexchange.com/a/450256/2530
# 
# Created: 2019.11.13
# Author: Steven E. Pav
# Comments: Steven E. Pav

SCRIPT_PATH="$( cd "$(dirname "$0")" ; pwd -P )"
if [[ "$TRAVIS" != "true" ]]; then
  TRAVIS_BUILD_DIR=$SCRIPT_PATH
  TRAVIS_BUILD_NUMBER="0000"
fi
export CONFIG_FILE="$TRAVIS_BUILD_DIR/.travis.yml"
# Get the tex-scheme config option
export texscheme=$(awk -F "=" '/tex-scheme/ {print $2}' "$CONFIG_FILE")
export pushtype=$(awk -F "=" '/push-type/ {print $2}' "$CONFIG_FILE")
# export pushtype="branch"
DOCKER_IMAGE="strauman/travis-latexbuild:$texscheme"

# Copyright 2015 Brent Longborough
# Part of gitinfo2 package Version 2
# Release 2.0.7 2015-11-22
# Please read gitinfo2.pdf for licencing and other details
# -----------------------------------------------------
# Post-{commit,checkout,merge} hook for the gitinfo2 package
gitinfotwo() {
#
# Get the first tag found in the history from the current HEAD
	FIRSTTAG=$(git describe --tags --always --dirty='-*' 2>/dev/null)
	# Get the first tag in history that looks like a Release
	RELTAG=$(git describe --tags --long --always --dirty='-*' --match '[0-9]*.*' 2>/dev/null)
	# Hoover up the metadata
	git --no-pager log -1 --date=short --decorate=short \
			--pretty=format:"\usepackage[%
					shash={%h},
					lhash={%H},
					authname={%an},
					authemail={%ae},
					authsdate={%ad},
					authidate={%ai},
					authudate={%at},
					commname={%cn},
					commemail={%ce},
					commsdate={%cd},
					commidate={%ci},
					commudate={%ct},
					refnames={%d},
					firsttagdescribe={$FIRSTTAG},
					reltag={$RELTAG}
			]{gitexinfo}" HEAD > .git/gitHeadInfo.gin
}
setup_git() {
  if [ "$TRAVIS" == "true" ]; then
    echo "Testing on travis-ci...";
    git config --global user.email "travis@travis-ci.org"
    git config --global user.name "Travis CI"
  else
    echo "Testing locally...";
  fi
  git checkout --orphan "travis-$TRAVIS_BUILD_NUMBER"
  git reset
	git status
	if [[ $(git ls-files) ]]; then
		echo "will git ls_files: "
		git ls-files
		echo "done"
		echo "will git rm"
  	git rm -f --ignore-unmatch --cached $(git ls-files)
		echo "done with rm"
	fi
	echo "done with setup"
}

commit_pdfs() {
  #echo `ls tests`
  git add -f "*.pdf"
	echo $(git status)
  git commit --message "Travis build: $TRAVIS_BUILD_NUMBER"
}

upload_files() {
	echo "will add remote"
  git remote add origin-login "https://${GH_TOKEN}@github.com/$TRAVIS_REPO_SLUG.git"
	echo "will push"
  git push --set-upstream -f origin-login "travis-$TRAVIS_BUILD_NUMBER"
  echo "PUSHED PDFS TO BRANCH travis-$TRAVIS_BUILD_NUMBER"
}
# Only execute if branch doesn't start with travis-
if [[ "$TRAVIS_BRANCH" == travis-* ]]; then
  echo "On a travis branch ($TRAVIS_BRANCH). Not doing anything."
else
	gitinfotwo;
  # Now pull the appropriate docker
  docker pull $DOCKER_IMAGE
  # Run the docker and on the files
  docker run --mount src="$TRAVIS_BUILD_DIR/",target=/repo,type=bind $DOCKER_IMAGE
  exitCode=$?
  if [ $exitCode -ne 0 ]; then
    exit $exitCode;
  fi
  if [ "$pushtype" == "branch" ]; then
    echo "PUSHTYPE: $pushtype";
    setup_git;
		echo "will commit pdfs?"
    commit_pdfs;
		echo "done commit pdfs."
    [[ "$TRAVIS" == "true" ]] && upload_files;
  elif [ "$pushtype" == "release" ]; then
    echo "Push type release not supported yet."
  fi
fi

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=sh:ft=sh:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
