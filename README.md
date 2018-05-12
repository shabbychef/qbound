
# Bounds on Portfolio Quality

This repo holds the code to rebuild the paper, "Bounds on Portfolio Quality." 
There are a few ways to build the paper:

1. Build it locally on your computer via `knitr` and `R`. This will require
the following packages:
`knitr`, `devtools`, `doFuture`, `doRNG`, `dplyr`, `ggplot2`, `hypergeo`, `knitr`, `LambertW`, 
`quantmod`, `SharpeR`, `tidyr`, `xtable`, 
and `aqfb.data`, which is available from github via
```r
install_github("shabbychef/aqfb_data")
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
The doc will be deposited in `output/qbound.pdf`. If you go this route, you can
log in to bash first and modify the `qbound.Rnw` file to change the
`RUNTIME_PARAM` to something larger than 1. That would look as follows:

 ```bash
$ docker pull shabbychef/qbound
$ mkdir ./output
$ docker run -it --rm -v $(pwd)/output:/srv/output:rw --entrypoint="/bin/bash" shabbychef/qbound "-i"
root@ddc0bdae65d9:/srv#  echo "you are in docker now; edit the file and run"
root@ddc0bdae65d9:/srv#  vi qbound.Rnw
root@ddc0bdae65d9:/srv#  make doc
root@ddc0bdae65d9:/srv#  exit
$ ls output/
```


## Testing

The simulations in this document can take several hours, but can take advantage of multicore
machines. To speed up the build with fewer simulations, adjust the parameter
`RUNTIME_PARAM` in the `qbound.Rnw` to a number bigger than 1, something like 50, say,
to get an approximately 50x speedup in simulation runtime. 


