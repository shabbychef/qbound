## ----'preamble', include=FALSE, warning=FALSE, message=FALSE-------------
library(knitr)

# set the knitr options ... for everyone!
# if you unset this, then vignette build bonks. oh, joy.
#opts_knit$set(progress=TRUE)
opts_knit$set(eval.after='fig.cap')
# for a package vignette, you do want to echo.
# opts_chunk$set(echo=FALSE,warning=FALSE,message=FALSE)
opts_chunk$set(warning=FALSE,message=FALSE)
#opts_chunk$set(results="asis")
opts_chunk$set(cache=TRUE,cache.path="cache/qbound")

#opts_chunk$set(fig.path="figure/",dev=c("pdf","cairo_ps"))
opts_chunk$set(fig.path="figure/qbound",dev=c("pdf"))
opts_chunk$set(fig.width=5,fig.height=4,dpi=64)

# doing this means that png files are made of figures;
# the savings is small, and it looks like shit:
#opts_chunk$set(fig.path="figure/",dev=c("png","pdf","cairo_ps"))
#opts_chunk$set(fig.width=4,fig.height=4)
# for figures? this is sweave-specific?
#opts_knit$set(eps=TRUE)

# this would be for figures:
#opts_chunk$set(out.width='.8\\textwidth')
# for text wrapping:
options(width=64,digits=2)
opts_chunk$set(size="small")
opts_chunk$set(tidy=TRUE,tidy.opts=list(width.cutoff=50,keep.blank.line=TRUE))

compile.time <- Sys.time()

# from the environment

# only recompute if FORCE_RECOMPUTE=True w/out case match.
FORCE_RECOMPUTE <- 
	(toupper(Sys.getenv('FORCE_RECOMPUTE',unset='False')) == "TRUE")


library(quantmod)
options("getSymbols.warning4.0"=FALSE)

library(SharpeR)

gen_norm <- rnorm
lseq <- function(from,to,length.out) { 
	exp(seq(log(from),log(to),length.out = length.out))
}

## ----'runtime', include=FALSE, warning=FALSE, message=FALSE, cache=TRUE----
# compiler flags!

# you can manipulate the 'runtime' by setting this variable.
# by convention the value 1 should take the longest time,
# while larger values should result in quicker execution
# (for low fidelity simulations)

RUNTIME_PARAM <- 1

## ----'prelook',eval=TRUE,echo=FALSE,cache=TRUE---------------------------
qbound <- function(p, n, zeta.s, use.ope=1) {
	# the upper bound on quality of the theorem.
	n.zet <- sqrt(n) * zeta.s
	retval <- sqrt(use.ope) * (n.zet * zeta.s) / sqrt(p - 1 + n.zet^2)
}
ex.n <- 5
ex.z <- 1.0
ex.p <- 10
myubound <- qbound(ex.p, ex.n, ex.z, 1)

## ----'haircut_study_setup',eval=TRUE,echo=FALSE,cache=TRUE---------------
suppressMessages({
	library(SharpeR)
  library(dplyr)
  library(tidyr)
})

# speed tests:
#X <- matrix(rnorm(100),nrow=1)
#norm(X,'2') - sqrt(X %*% t(X))
#microbenchmark('norm'={ norm(X,'2') },'ogram'={ sqrt(X %*% t(X)) })
# note that colMeans is afaster than apply
#library(microbenchmark)
#X <- matrix(rnorm(10000),ncol=100)
#microbenchmark('colmeans'={ colMeans(X,na.rm=TRUE) },'apply'={ as.vector(apply(X,MARGIN=2,mean,na.rm=TRUE)) })
# simple markowitz.
simple.marko <- function(rets) {
	mu.hat <- colMeans(rets,na.rm=TRUE)
	Sig.hat <- cov(rets)
	w.opt <- solve(Sig.hat,mu.hat)
	retval <- list('mu'=mu.hat,'sig'=Sig.hat,'w'=w.opt)
	return(retval)
}

# function that generates 0 mean, unit sd random variables as a vector#FOLDUP
genz.gauss <- rnorm
sq3 <- sqrt(3)
genz.uni <- function(n) {
	runif(n,min=-sq3,max=sq3)
}
make_genz.t <- function(df) {
	true.v <- df / (df - 2)
	adjust <- sqrt(1/true.v)
	function(n) { adjust * rt(n,df=df) }
}
make_genz.tukeyh <- function(h) {
	require(LambertW,quietly=TRUE)
	Gauss_input <- create_LambertW_input("normal", beta=c(0,1))
	params <- list(delta = c(h))
	LW.Gauss <- create_LambertW_output(LambertW.input = Gauss_input, theta = params)
	#get the moments of this distribution
	moms <- mLambertW(theta=list(beta=c(0,1),delta = h,gamma = 0,alpha = 1),distname=c("normal"))
	function(n) {
		Z <- LW.Gauss$r(n)
		Z <- (1 / moms$sd) * (Z - moms$mean)
		Z
	}
}
make_genz.lambertw <- function(gam) {
	require(LambertW,quietly=TRUE)
	Gauss_input = create_LambertW_input("normal", beta=c(0,1))
	params = list(delta = c(0), gamma=c(gam), alpha = 1)
	LW.Gauss = create_LambertW_output(LambertW.input = Gauss_input, theta = params)
	#get the moments of this distribution
	moms <- mLambertW(theta=list(beta=c(0,1),delta = 0,gamma = gam, alpha = 1), distname=c("normal"))
	function(n) {
		Z <- LW.Gauss$r(n)
		Z <- (1 / moms$sd) * (Z - moms$mean)
		Z
	}
}
#UNFOLD

# make multivariate pop. & sample w/ given zeta.star
gen.pop <- function(n,p,zeta.s,genf=genz.gauss,...) {
	true.mu <- matrix(rnorm(p),ncol=p)
	#generate an SPD population covariance. a hack.
	# maybe Bartlett would be faster
	#xser <- matrix(rnorm(p*(p + 50)),ncol=p)
	#true.Sig <- t(xser) %*% xser
	# since simple markowitz has the same performance invariant
	# to rotation, we can assume the true Sigma is diagonal
	true.Sig <- diag(p)
	pre.sr <- sqrt(true.mu %*% solve(true.Sig,t(true.mu)))
	#scale down the sample mean to match the zeta.s
	true.mu <- (zeta.s/pre.sr[1]) * true.mu 
  Z <- genf(n*p)
	Z <- matrix(Z,nrow=n,ncol=p)
	# true.Sig is eye, so no need to chol it. save some time
	#Z <- Z %*% chol(true.Sig)
	X <- Z + matrix(rep(true.mu,n),nrow=n,byrow=TRUE)
	retval <- list('X'=X,'mu'=true.mu,'sig'=true.Sig,'SNR'=zeta.s)
	return(retval)
}
# a single simulation
sample.haircut <- function(n,p,zeta.s,genf=genz.gauss,...) {
	popX <- gen.pop(n=n,p=p,zeta.s=zeta.s,genf=genf,...)
	smeas <- simple.marko(popX$X)
	# popX$sig is eye
	#ssnr <- (t(smeas$w) %*% t(popX$mu)) / sqrt(t(smeas$w) %*% popX$sig %*% smeas$w)
	ssnr <- (t(smeas$w) %*% t(popX$mu)) / sqrt(t(smeas$w) %*% smeas$w)
	as.numeric(ssnr)
}

# test
#sample.haircut(1000,10,zeta.s=0.1)

# set everything up
ope <- 253

n.sim <- ceiling(1e6/max(1,RUNTIME_PARAM))
n.per <- ceiling(1e3/max(1,RUNTIME_PARAM))

n.stok <- 6
n.yr <- 4
n.obs <- ceiling(ope * n.yr)
zeta.s <- 1.25 / sqrt(ope)   # optimal SNR, in daily units

t.df <- 4
tuk.h <- 0.15
lam.gam <- -0.2

repsim <- function(nsamp, nobs, nstok, zeta.s, genf, rseed=NULL) {
	if (!is.null(rseed)) {
		set.seed(rseed)
	}
	experiments <- replicate(nsamp, sample.haircut(n=nobs,p=nstok,zeta.s=zeta.s,genf=genf))
	data_frame(ssnr=sqrt(ope) * experiments,
						 n=nobs,p=nstok,zeta=zeta.s)
}
# testing
#repsim(100,genf=genz.gauss)
#repsim(nsamp=100,zeta=0.1,nstok=5,nobs=100,genf=genz.gauss)

suppressMessages({
  library(doRNG)
	registerDoRNG()
  # https://cran.r-project.org/web/packages/doFuture/vignettes/doFuture.html
  library(doFuture)
  registerDoFuture()
})

manysim <- function(nsamp, genf, nobs=n.obs, nstok=n.stok, zeta=zeta.s, nper=n.per, rseed=1234) {
  require(doFuture)
  registerDoFuture()
  plan(multiprocess)

	quals <- future_lapply(1:ceiling(nsamp/nper),
												 FUN=function(...) {
													 repsim(nsamp=nper,genf=genf,nobs=nobs,nstok=nstok,zeta=zeta)
												 },
												 future.seed=rseed) %>%
		bind_rows()

	# might have to down select
	quals <- quals[1:nsamp,]
	return(quals)
}

t0 <- Sys.time()

# gaussian
quals.gaussian <- manysim(n.sim,genf=genz.gauss)$ssnr

# uniform
quals.uni <- manysim(n.sim,genf=genz.uni)$ssnr

# t(4)
quals.t <- manysim(n.sim,genf=make_genz.t(t.df))$ssnr

# tukey(0.1)
quals.tukey <- manysim(n.sim,genf=make_genz.tukeyh(tuk.h))$ssnr

# Lambert W(0.1)
quals.lambert <- manysim(n.sim,genf=make_genz.lambertw(lam.gam))$ssnr

## ----'hcut_qtile',echo=FALSE,cache=FALSE---------------------------------
# print(summary(quals.gaussian))
# haircut approximation in the equation above
qqual <- function(p, df1, df2, zeta.s, use.ope=1, lower.tail=TRUE) {
	# incorporates the annualization factor
	# uses the sin(atan( ... )) formulation:
	atant <- atan((1/sqrt(df1-1)) * 
		qt(p,df=df1-1,ncp=sqrt(df2)*zeta.s,lower.tail=lower.tail))
	# a slightly better approximation is:
	# retval <- 1 - sin(atant - 0.0184 * zeta.s * sqrt(df1 - 1))
	retval <- sqrt(use.ope) * zeta.s * sin(atant)

	# use the beta formulation, but oh fuck, no negatives 
	# on this one. shit.
#	retval <- zeta.s * sqrt(use.ope * qbeta(p, shape1=0.5, shape2=0.5 * (df1-1),
#		ncp = df2 * zeta.s^2, lower.tail=lower.tail))

}
pqual <- function(q, df1, df2, zeta.s, use.ope=1, lower.tail=TRUE, log.p=FALSE) {
	sinthet <- pmin(pmax(q / (zeta.s * sqrt(use.ope)),-1),1)
	fooq <- sqrt(df1-1) * tan(asin(sinthet))
	retval <- pt(fooq,df=df1-1,ncp=sqrt(df2)*zeta.s,
				lower.tail=lower.tail,log.p=log.p)
	#retval <- min(max(retval,0),1)
	#retval
}
qbound <- function(p, n, zeta.s, use.ope=1) {
	# the upper bound on quality of the theorem.
	n.zet <- sqrt(n) * zeta.s
	retval <- sqrt(use.ope) * (n.zet * zeta.s) / sqrt(p - 1 + n.zet^2)
}
#apx.medv <- qqual(0.5,n.stok,n.obs,zeta.s,use.ope=ope)
#corr.digs <- function(true.v,est.v) {
#	err <- est.v - true.v
#	true.digs <- ceiling(log10(abs(true.v)))
#	err.digs <- trunc(log10(abs(err)))
#	n.signif <- true.digs - err.digs
#}
#n.digs <- corr.digs(median(quals.gaussian),apx.medv)

## ----'kstestit',eval=TRUE,cache=TRUE,echo=FALSE--------------------------

x <- quals.gaussian
x <- x[!is.na(x)]
sim.ks <- ks.test(x, y=pqual, df1=n.stok, df2=n.obs,
                    zeta.s=zeta.s, use.ope=ope)

ksp <- signif(sim.ks$p.value,digits=8)
kss.gaussian <- signif(sim.ks$statistic,digits=4)

get.stat <- function(x.quals) {
	x.quals <- x.quals[!is.na(x.quals)]
	res.ks <- ks.test(x.quals, y=pqual, df1=n.stok, df2=n.obs,
											zeta.s=zeta.s, use.ope=ope)
	res.stat <- signif(res.ks$statistic,digits=4)
}

kss.uni <- get.stat(quals.uni)
kss.t <- get.stat(quals.t)
kss.tukey <- get.stat(quals.tukey)
kss.lambert <- get.stat(quals.lambert)

## ----'damnoptions1',eval=TRUE,cache=FALSE,echo=FALSE---------------------
options(digits=3)

## ----'damnoptions2',eval=TRUE,cache=FALSE,echo=FALSE---------------------
options(digits=2)

## ----'haircutting',echo=FALSE,cache=TRUE,fig.width=5.00,fig.height=3.0,dpi=450,fig.cap=paste0("Q-Q plot of $10^{",round(log10(n.sim)),"}$ simulated \\txtQual values versus \\apxref{qual_dist} is shown. Units are `annual', \\ie \\yrtomhalf. Since the number of samples is very large, only a subset of $10^{",round(log10(max.ps)),"}$ points, uniformly selected by sample quantile, are plotted.")----
# qqplot;

# via qqplot;
#qqplot(qqual(ppoints(length(quals.gaussian)),n.stok,n.obs,zeta.s,use.ope=ope),quals.gaussian,
#			 xlab = "Theoretical Approximate Quantiles", ylab = "Sample Quantiles")
#qqline(quals.gaussian,datax=FALSE,distribution = function(p) {
# qqual(p,n.stok,n.obs,zeta.s,use.ope=ope) },
#			 col=2)

suppressMessages({
  library(ggplot2)
})

max.ps <- 1e4
if (length(quals.gaussian) > max.ps) {
	pvs <- ppoints(max.ps)
	subquals <- quantile(quals.gaussian, probs=pvs)
} else {
	subquals <- quals.gaussian
	pvs <- ppoints(subquals)
}

# or, see: 
# http://stackoverflow.com/a/4939257/164611
xydf <- as.data.frame(qqplot(x=qqual(pvs,n.stok,n.obs,zeta.s,use.ope=ope),
	y=subquals,plot=F))
ph <- ggplot(data=xydf,mapping=aes(x=x, y=y)) +
  geom_point() +
  geom_abline(intercept=0, slope=1, color="red") + 
  labs(x="Theoretical Approximate Quantiles", y="Sample Quantiles")

#ph <- ph + geom_smooth(method="lm", se=FALSE)	
print(ph)


## ----'qtiles',echo=FALSE,cache=TRUE,results='asis'-----------------------

qs <- c(0.1,0.25,0.5,0.75,0.95,0.975,0.99,0.995)
qs <- sort(1 - qs)
foo.tab <- data.frame(Gaussian=quantile(quals.gaussian,qs),
	uniform=quantile(quals.uni,qs),
	t=quantile(quals.t,qs),
	tukey=quantile(quals.tukey,qs),
	lambert=quantile(quals.lambert,qs),
	approximation=qqual(qs,n.stok,n.obs,zeta.s,use.ope=ope))
colnames(foo.tab) <- c("normal",
	"unif.",
	sprintf("t(%d)",t.df),
	sprintf("Tukey(%.2f)",tuk.h),
	sprintf("Lam.W(%.1f)",lam.gam),
	"approx.")
#rownames(foo.tab) <- 100 * qs
foo.tab <- cbind(data.frame("q-tile"=1*qs),foo.tab)

	#align=paste0(paste0(rep("l",dim(foo.tab)[2]-1),collapse=""),"|","l"),
	#align=c(rep("l",dim(foo.tab)[2]),"|","l"),
	#align=rep("c",dim(foo.tab)[2]+1),
library(xtable)
xres <- xtable(foo.tab,
	digits=c(3,3,rep(4,dim(foo.tab)[2]-1)),
	label="tab:qtiles",
	align=c("c","r|",rep("c",dim(foo.tab)[2]-2),"|c"),
	caption=paste0("Empirical quantiles of portfolio \\txtQual",
		" from $10^{",log10(n.sim),"}$",
		" simulations of ",n.obs," days of ",n.stok," assets, with",
		" maximal \\txtSR of $",zeta.s * sqrt(ope),"\\yrtomhalf$ are",
		" given, along with the approximate quantiles from",
		" \\apxref{qual_dist}. Units of \\txtQual are `annual',",
		" \\ie \\yrtomhalf."))
#print(xres,include.rownames=TRUE)
print(xres,include.rownames=FALSE)

## ----'ks_table',echo=FALSE,cache=TRUE,results='asis'---------------------

foo.tab <- data.frame(Gaussian=kss.gaussian,
	uniform=kss.uni,
	t=kss.t,
	tukey=kss.tukey,
	lambert=kss.lambert)
colnames(foo.tab) <- c("normal",
	"unif.",
	sprintf("t(%d)",t.df),
	sprintf("Tukey(%.2f)",tuk.h),
	sprintf("Lam.W(%.1f)",lam.gam))

	#align=paste0(paste0(rep("l",dim(foo.tab)[2]-1),collapse=""),"|","l"),
	#align=c(rep("l",dim(foo.tab)[2]),"|","l"),
	#align=rep("c",dim(foo.tab)[2]+1),
library(xtable)
xres <- xtable(foo.tab,
	digits=c(3,rep(4,dim(foo.tab)[2])),
	label="tab:ks_stats",
	align=c("c",rep("c",dim(foo.tab)[2])),
	caption=paste0("\\txtKS statistic comparing the empirical",
		" CDF to that of \\apxref{qual_dist}",
		" over $10^{",log10(n.sim),"}$",
		" simulations of ",n.obs," days of ",n.stok," assets, with",
		" maximal \\txtSR of $",zeta.s * sqrt(ope),"\\yrtomhalf$ are",
		" given for the different returns distributions."))
#print(xres,include.rownames=TRUE)
print(xres,include.rownames=FALSE)

## ----'means',echo=FALSE,cache=TRUE,results='asis'------------------------

ubound <- qbound(n.stok, n.obs, zeta.s, ope)

require(hypergeo)

teff <- 0.5 * n.obs * zeta.s^2
U <- c(n.stok/2,3/2)
L <- c((2+n.stok)/2,1/2)
ub2 <- ope * zeta.s^2 * exp(- teff + sum(lgamma(U)) - sum(lgamma(L)))
ub2 <- ub2 * hypergeo::genhypergeo(U, L, teff)

aggf <- function(x) {
	retv <- c(mean(x,na.rm=TRUE),mean(x^2,na.rm=TRUE))
}

foo.tab <- data.frame(Gaussian=aggf(quals.gaussian),
	uniform=aggf(quals.uni),
	t=aggf(quals.t),
	tukey=aggf(quals.tukey),
	lambert=aggf(quals.lambert),
	bound=c(ubound,ub2))
colnames(foo.tab) <- c("normal",
	"unif.",
	sprintf("t(%d)",t.df),
	sprintf("Tukey(%.2f)",tuk.h),
	sprintf("Lam.W(%.1f)",lam.gam),
	"bound")

mean.tab <- foo.tab[1,]
meansq.tab <- foo.tab[2,]
colnames(meansq.tab)[length(colnames(meansq.tab))] <- "approx."

library(xtable)
xres <- xtable(mean.tab,
	digits=c(1,1,rep(3,dim(mean.tab)[2]-1)),
	label="tab:means",
	align=c("c",rep("c",dim(mean.tab)[2]-1),"|c"),
	caption=paste0("Empirical mean portfolio \\txtQual",
		" from $10^{",log10(n.sim),"}$",
		" simulations of ",n.obs," days of ",n.stok," assets, with",
		" maximal \\txtSR of $",zeta.s * sqrt(ope),"\\yrtomhalf$ are",
		" given, along with the upper bound from",
		" \\theoremref{qual_bound}. Units of \\txtQual are `annual',",
		" \\ie \\yrtomhalf."))
print(xres,include.rownames=FALSE)

xres2 <- xtable(meansq.tab,
	digits=c(1,1,rep(3,dim(meansq.tab)[2]-1)),
	label="tab:meanssq",
	align=c("c",rep("c",dim(meansq.tab)[2]-1),"|c"),
	caption=paste0("Empirical mean of \\emph{squared} portfolio \\txtQual",
		" from $10^{",log10(n.sim),"}$",
		" simulations of ",n.obs," days of ",n.stok," assets, with",
		" maximal \\txtSR of $",zeta.s * sqrt(ope),"\\yrtomhalf$ are",
		" given, along with the approximate value from \\eqnref{hypergeo_extwo}.",
		" Units of squared \\txtQual are `annual',",
		" \\ie \\yrto{-1}."))
print(xres2,include.rownames=FALSE)

mu.gap <- (mean.tab$bound - mean.tab$normal) / (mean.tab$bound)

## ----'quality_heatmaps',eval=TRUE,echo=FALSE,cache=TRUE------------------
suppressMessages({
  require(SharpeR)
})

ope <- 253

# set everything up
n.sim3 <- ceiling(1e5/max(1,RUNTIME_PARAM))
n.per3 <- ceiling(5e3/max(1,RUNTIME_PARAM))

params <- tidyr::crossing(tibble::tribble(~nyr,0.5,1,2,4,8),
                          tibble::tribble(~p,2,4,8,16),
                          tibble::tibble(zeta=sqrt(c(0.125,0.25,0.5,1,2)) / sqrt(ope))) %>%
  mutate(nday=ceiling(ope * nyr)) 

registerDoRNG()
registerDoFuture()

sims <- params %>% 
  group_by(nyr,p,zeta,nday) %>%
    summarize(sims=list(manysim(nsamp=n.sim3,nper=n.per3,nobs=nday,nstok=p,zeta=zeta,genf=genz.gauss))) %>%
  ungroup() %>%
  unnest()

# also do ks.test ... 
ksfunc <- function(sims,df1,df2,zeta,use.ope=ope) {
  htest <- ks.test(sims,y=pqual, df1=df1, df2=df2, zeta.s=zeta,use.ope=use.ope)
  data_frame(stat=htest$statistic,pvalue=htest$p.value)
}

ks.stats <- sims %>%
  group_by(nyr,p,zeta,nday) %>%
    summarize(stats=list(ksfunc(ssnr,df1=p,df2=nday,zeta=zeta))) %>%
  ungroup() %>%
  unnest() %>%
  rename(nstock=p,value=stat)

## ----'ksplots',eval=TRUE,echo=FALSE,cache=FALSE,fig.width=5.75,fig.height=4.00,dpi=450,fig.cap=paste0("The \\txtKS statistic for \\apxref{qual_dist} over $10^{",round(log10(n.sim3)),"}$ simulations of Gaussian returns is plotted versus $\\psnropt \\wrapParens{\\nlatf-1}/\\sqrt{\\ssiz}$, with \\psnropt in annualized terms, and \\ssiz measured in years.  There is one line for each combination of $\\ssiz$ and $\\psnropt$. The line color corresponds to the `effect size', $\\sqrt{\\ssiz}\\psnropt$, which is unitless."),eval.after='fig.cap'----

ks.dat <- ks.stats %>%
  mutate(zeta.an = sqrt(ope) * zeta) %>%
  mutate(effsize = sqrt(nyr) * zeta.an) %>%
  mutate(pbyn = zeta.an * (nstock - 1) / sqrt(nyr)) %>%
  mutate(effect.size = factor(round(effsize,digits=3))) %>%
  mutate(nyr = factor(nyr)) %>%
  mutate(zeta.an = factor(zeta.an))

suppressMessages({
  library(ggplot2)
})

#ph <- ph + facet_grid(.~zeta.an)
ph <- ks.dat %>%
  ggplot(aes(x=pbyn,y=value,group=interaction(nyr,zeta.an))) +
  geom_line(aes(colour=effect.size)) +
  labs(x=expression(zeta["*"] * (p-1) / sqrt(n)),y="KS statistic") +
  guides(colour=guide_legend(title=expression(sqrt(n)*zeta["*"]))) +
  scale_x_log10()
#ph <- ph + scale_y_log10()
print(ph)


## ----'ksheat',eval=TRUE,echo=FALSE,cache=FALSE,fig.width=6.00,fig.height=4.25,dpi=450,fig.cap=paste0("The \\txtKS statistic for \\apxref{qual_dist} over $10^{",round(log10(n.sim3)),"}$ simulations of Gaussian returns is indicated, by color, versus \\nlatf, and the `total effect size,' $\\sqrt{\\ssiz}\\psnropt$, which is a unitless quantity. Different facets are for different values of \\ssiz (in years)."),eval.after='fig.cap'----

ks.dat <- ks.stats %>%
  mutate(zeta.an = sqrt(ope) * zeta) %>%
  mutate(effsize = round(sqrt(nyr) * zeta.an,digits=3)) %>%
  mutate(effsize = factor(effsize)) %>%
  mutate(nstock = factor(nstock)) %>%
  mutate(nyr = factor(nyr)) 

suppressMessages({
  library(ggplot2)
})

ph <- ks.dat %>%
  rename(years=nyr) %>%
  ggplot(aes(x=nstock,y=effsize)) +
  geom_tile(aes(fill=value)) +
  facet_grid(.~years,labeller=label_both) +
  scale_fill_gradient2(limits=c(min(ks.dat$value),max(ks.dat$value)),
                       low="steelblue",high="red") +
  labs(x="# assets",y="effect size") +
  guides(fill=guide_legend(title="KS stat"))

#ph <- ph + scale_y_log10()
print(ph)


## ----'forexample',echo=FALSE,cache=FALSE---------------------------------
# reset these, just in case...
ope <- 253
n.stok <- 6
n.yr <- 4
n.obs <- ceiling(ope * n.yr)
zeta.s <- 1.25 / sqrt(ope)   # optimal SNR, in daily units

big.n.stok <- 24
big.zeta <- 1.60 / sqrt(ope)   # optimal SNR, in daily units
big.qb <- qbound(big.n.stok, n.obs, big.zeta, ope)
lil.qb <- qbound(n.stok, n.obs, zeta.s, ope)

## ----'grow_bound',eval=TRUE,echo=FALSE,cache=TRUE,fig.width=5.00,fig.height=3.25,dpi=450,fig.cap=paste0("The upper bound of \\theoremref{qual_bound} is plotted versus \\nlatf for different scaling laws for \\psnropt.  These scaling laws correspond to $\\psnropt = \\psnr[0]\\nlatf^{\\gamma}$, with $\\gamma$ taking values ",sgammas.summary,".  The constant terms, \\psnr[0], are adjusted so that $\\psnropt =",sqrt(ope) * zeta.s,"\\yrtomhalf$ for $\\nlatf = ",n.stok,"$ for all the lines. The bound uses $\\ssiz=",n.obs,"$, corresponding to ",n.yr," years of daily observations."),eval.after='fig.cap'----

bound.experiment <- function(pow,n.obs,zeta.s,
                             n.stok=n.stok,ope=ope,plims=c(2,250)) {
  require(dplyr,quietly=TRUE)
  require(tidyr,quietly=TRUE)
	all.ps <- unique(round(exp(seq(log(min(plims)),
																	log(max(plims)),
																	length.out=140))))
	zeta.0 <- zeta.s / (n.stok ^ pow)
	all.zeta <- zeta.0 * (all.ps ^ pow)
	all.bnd <- sqrt(ope * n.obs) * all.zeta^2 / sqrt(all.ps - 1 + n.obs * all.zeta^2)

	#population.max=sqrt(ope) * all.zeta,
	foo.df <- data_frame(p=all.ps,
                       pow=rep(pow,length(all.bnd)),
                       bound=all.bnd)
  retv <- foo.df %>% tidyr::gather(key=Sharpe,value=value,-p,-pow)
}

#	plims = c(2,250)
#	melt.df <- NULL
#	for (ppp in pow) {
#		addon.df <- bound.experiment(ppp,n.obs,zeta.s,n.stok=n.stok,ope=ope,plims=plims)
#		if (is.null(melt.df))
#			melt.df <- addon.df
#		else
#			melt.df <- rbind(melt.df,addon.df)
#	}
#	#as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
#	melt.df$pow <- as.numeric(melt.df$pow)
#
#	require(ggplot2)
#	ph <- ggplot(data=melt.df,aes(x=p,y=value,group=pow,colour=pow))
#	ph <- ph + geom_line()
#	ph <- ph + labs(x="# assets",
#									y="Signal-Noise ratio (annualized)")
#	ph <- ph + scale_x_log10(limits = c(1,max(plims)))
#	ph <- ph + scale_y_log10(limits = c(0.5,4),
#		breaks=2^seq(from=-1,to=2,by=0.5))
#	print(ph)

bound.plot <- function(pow,n.obs,zeta.s,plims=c(2,250),...) {
	melt.df <- NULL
	for (ppp in pow) {
		addon.df <- bound.experiment(ppp,n.obs,zeta.s,plims=plims,...)
		if (is.null(melt.df))
			melt.df <- addon.df
		else
			melt.df <- rbind(melt.df,addon.df)
	}
	melt.df$gamma <- factor(melt.df$pow)

	require(ggplot2,quietly=TRUE)
	ph <- ggplot(data=melt.df,aes(x=p,y=value,group=gamma,colour=gamma))
	ph <- ph + geom_line()
	ph <- ph + labs(x="# assets",
									y="Signal-Noise ratio (annualized)")
	ph <- ph + guides(colour=guide_legend(title=expression(gamma)))
	ph <- ph + scale_x_log10(limits = c(1,max(plims)))
	ph <- ph + scale_y_log10(limits = c(0.5,4),
		breaks=2^seq(from=-1,to=2,by=0.5))
	return(ph)
}

#pgammas <- c(0.15,0.20,0.25,0.30)
#pgammas <- c(0.19,0.22,0.25,0.28)
#pgammas <- seq(from=0.19,to=0.28,by=0.03)
pgammas <- seq(from=0.15,to=0.35,by=0.05)

sgammas <- sort(pgammas)
sgammas.txt <- paste0(paste0(sgammas[1:(length(sgammas)-1)],collapse=', '),", and ",sgammas[length(sgammas)],sep='')
sgammas.summary <- paste0("between ",sgammas[1]," and
",sgammas[length(sgammas)])

print(bound.plot(pgammas,n.obs=n.obs,zeta.s=zeta.s,
	n.stok=n.stok,ope=ope))


## ----'grow_sqrt',eval=TRUE,echo=FALSE,cache=TRUE,fig.width=5.00,fig.height=3.25,dpi=450,fig.cap=paste0("Some quantiles of the \\txtQual of the \\txtMP, under \\apxref{qual_dist}, are plotted versus \\nlatf for different scaling laws for \\psnropt. The ",length(sgammas)," panels represent different values of $\\gamma$, \\viz ", sgammas.txt, ". The bound uses $\\ssiz=",n.obs,"$, corresponding to ",n.yr," years of daily observations."),eval.after='fig.cap'----

scale.experiment <- function(pow,n.obs,zeta.s,
                             n.stok=n.stok,ope=ope,plims=c(2,250)) {
  require(dplyr,quietly=TRUE)
  require(tidyr,quietly=TRUE)
	all.ps <- unique(round(exp(seq(log(min(plims)),
																	log(max(plims)),
																	length.out=140))))
	zeta.0 <- zeta.s / (n.stok ^ pow)
	all.zeta <- zeta.0 * (all.ps ^ pow)
	all.75 <- qqual(0.75,all.ps,n.obs,all.zeta,use.ope=ope)
	all.50 <- qqual(0.5,all.ps,n.obs,all.zeta,use.ope=ope)
	all.25 <- qqual(0.25,all.ps,n.obs,all.zeta,use.ope=ope)
	foo.df <- data_frame(p=all.ps,
		`pop. max`=sqrt(ope) * all.zeta,
		pow=pow,
		`.75 q`=as.numeric(all.75),
		`.50 q`=as.numeric(all.50), 
		`.25 q`=as.numeric(all.25))
  retv <- foo.df %>%
    tidyr::gather(key=Sharpe,value=value,-p,-pow)
}

scale.plot <- function(pow,n.obs,zeta.s,plims=c(2,250),...) {
	melt.df <- NULL
	for (ppp in pow) {
		addon.df <- scale.experiment(ppp,n.obs,zeta.s,plims=plims,...)
		if (is.null(melt.df))
			melt.df <- addon.df
		else
			melt.df <- rbind(melt.df,addon.df)
	}
	#melt.df$pow <- factor(melt.df$pow)

	require(ggplot2,quietly=TRUE)
	require(dplyr,quietly=TRUE)
	ph <- melt.df %>%
    rename(`gamma`=pow) %>%
    ggplot(aes(x=p,y=value,colour=Sharpe)) +
      geom_line() +
      facet_grid(.~`gamma`,labeller=label_both) +
      labs(x="# assets",
           y="Signal-Noise ratio (annualized)") +
      guides(colour=guide_legend(title="SNR")) +
      scale_x_log10(limits = c(1,max(plims))) +
      scale_y_log10(limits = c(0.5,4), 
                    breaks=2^seq(from=-1,to=2,by=0.5))
	return(ph)
}

#pgammas <- c(0.15,0.20,0.25,0.30)
#pgammas <- c(0.19,0.22,0.25,0.28)
#pgammas <- seq(from=0.19,to=0.28,by=0.03)
pgammas <- seq(from=0.21,to=0.29,by=0.04)

sgammas <- sort(pgammas)
sgammas.txt <- paste0(paste0(sgammas[1:(length(sgammas)-1)],collapse=', '),", and ",sgammas[length(sgammas)],sep='')

print(scale.plot(pgammas,n.obs=n.obs,zeta.s=zeta.s,
	n.stok=n.stok,ope=ope))


## ----'qlook',echo=FALSE,cache=FALSE--------------------------------------
if (!require(aqfb.data,quietly=TRUE) && require(devtools)) {
  # get the 10 industry data
  devtools::install_github('shabbychef/aqfb_data')
}
# now look at the data
library(aqfb.data)
# need mind10?
data(mind10)
ff11.xts <- mind10
require(SharpeR)

ff11.xts <- na.omit(ff11.xts)

obn.ret <- function(x.xts) {
	mur.xts <- xts((1/dim(x.xts)[2]) * rowSums(x.xts),order.by=time(x.xts))
}

mur.xts <- obn.ret(ff11.xts)

sr1 <- SharpeR::as.sr(mur.xts)
sr2 <- SharpeR::as.sropt(ff11.xts)

ret.TEO <- time(ff11.xts)
ret.TEO.0 <- as.Date(ret.TEO[1])
ret.TEO.f <- as.Date(ret.TEO[length(ret.TEO)])

dgu.yr <- 5
dgqb <- qbound(dim(ff11.xts)[2], dgu.yr, sr2$sropt, use.ope=1)
# upper tail prob:
dfp <- pqual(sr1$sr, dim(ff11.xts)[2], dgu.yr, sr2$sropt, use.ope=1, lower.tail=FALSE)

## ----'spantest',echo=FALSE,cache=TRUE------------------------------------

G <- matrix(1,nrow=1,ncol=dim(ff11.xts)[2])
dsr2 <- SharpeR::as.del_sropt(ff11.xts,G)
dKRS <- inference(dsr2,type='KRS')
dMLE <- inference(dsr2,type='MLE')
dUNB <- inference(dsr2,type='unbiased')

myci <- confint(dsr2)

asqb <- qbound(dim(ff11.xts)[2]-1, dgu.yr, sqrt(dKRS), use.ope=1)


## ----'sp100_divs',cache=FALSE,eval=TRUE,echo=FALSE-----------------------
load('sp100lr.rda')

asof <- 'March 21, 2014'

sp1.TEO <- time(sub.lr)
sp1.TEO.0 <- as.Date(sp1.TEO[1])
sp1.TEO.f <- as.Date(sp1.TEO[length(sp1.TEO)])


## ----'sp100_grow',eval=TRUE,echo=FALSE,cache=TRUE,fig.width=5.00,fig.height=3.25,dpi=450,fig.cap=paste0("Growth of estimated \\psnropt versus \\nlatf for the S\\&P 100 Index names, in alphabetical order, showing the `Apple effect.'"),eval.after='fig.cap'----

foo.df <- data.frame(df=seq(1:length(KRSs)),KRS=KRSs,
		MLE=MLEs,meanKRS=rowMeans(buncho.KRSs))

require(ggplot2)
ph <- ggplot(data=foo.df,aes(x=df,y=KRS))
ph <- ph + geom_line()
ph <- ph + labs(x="# assets",
								y=expression(zeta["*"]))
ph <- ph + scale_x_log10()
#ph <- ph + scale_y_log10()
print(ph)


## ----'sp100_box',eval=TRUE,echo=FALSE,cache=TRUE,fig.width=5.00,fig.height=3.25,dpi=450,fig.cap=paste0("Growth of estimated \\psnropt versus \\nlatf for the S\\&P 100 Index names is shown over ",dim(buncho.KRSs)[2]," permutations of the stocks. There is effectively \\emph{no} diversification benefit here beyond an equal weight portfolio."),eval.after='fig.cap'----

foo.df <- data.frame(df=rep(1:(dim(buncho.KRSs)[1]),dim(buncho.KRSs)[2]),
	KRS=as.vector(t(buncho.KRSs)))
foo.df <- foo.df[foo.df$df < 11,]

require(ggplot2)
ph <- ggplot(data=foo.df,aes(x=factor(df),y=KRS))
ph <- ph + geom_boxplot()
ph <- ph + labs(x="# assets",
								y=expression(zeta["*"]))
#ph <- ph + scale_x_log10()
#ph <- ph + scale_y_log10()
print(ph)


