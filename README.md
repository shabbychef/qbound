
# Bounds on Portfolio Quality

[![Build Status](https://travis-ci.org/shabbychef/qbound.svg?branch=master)](https://travis-ci.org/shabbychef/qbound)

This repo holds the code to rebuild the paper, "Bounds on Portfolio Quality." 
There are a few ways to build the paper:

1. Build it locally on your computer via `knitr` and `R`. This will require
the following packages:
`knitr`, `remotes`, `future.apply`, `dplyr`, `ggplot2`, `hypergeo`, `knitr`, `LambertW`, 
`quantmod`, `SharpeR`, `tidyr`, `xtable`, `viridis`, 
and `tsrsa`, which is available from github via
```r
remotes::install_github("shabbychef/tsrsa/rpkg")
```
Then you can build via
```r
knitr::knit('qbound.Rnw')
```

2. The doc can be built via docker, which constructs an environment containing the proper
	 packages. This is most simply achieved via a makefile target:
```bash
$ make docker_doc
```
The doc will be deposited in `output/qbound.pdf`.

3. The doc can be built via docker with the docker image pulled from docker hub:
 ```bash
$ docker pull shabbychef/qbound
$ mkdir ./output
$ docker run -it --rm -v $(pwd)/output:/srv/output:rw --entrypoint="make" shabbychef/qbound "doc"
```
The doc will be deposited in `output/qbound.pdf`. You can control the build
speed/resolution tradeoff with the environment variable `RUNTIME_PARAM` as
follows:

 ```bash
$ docker pull shabbychef/qbound
$ mkdir ./output
$ docker run -it --rm -v $(pwd)/output:/srv/output:rw -e RUNTIME_PARAM=50 --entrypoint="make" shabbychef/qbound "doc"
```

## Testing

The simulations in this document can take several hours, but can take advantage of multicore
machines. To speed up the build with fewer simulations, either adjust the
variable `RUNTIME_PARAM` in the `qbound.Rnw` to default to a number bigger than 1, something like 50, say,
to get an approximately 50x speedup in simulation runtime.  This can be
controlled via the environment variable of the same name:

```bash
$ RUNTIME_PARAM=50 make docker_doc
```

## Published Version

The version that I had sent to an MDPI journal can be made via a separate target:

```bash
$ make docker_mdpi
# or
$ docker pull shabbychef/qbound
$ mkdir ./output
$ docker run -it --rm -v $(pwd)/output:/srv/output:rw -e RUNTIME_PARAM=50 --entrypoint="make" shabbychef/qbound "mdpi"
```

Don't ask why I never paid to have the article published in that journal.

### Screencast

Watch me stumble around building the doc:
[![asciicast](https://asciinema.org/a/nianUC5GtAQgJeO7eB2Zu3FyJ.png)](https://asciinema.org/a/nianUC5GtAQgJeO7eB2Zu3FyJ).

