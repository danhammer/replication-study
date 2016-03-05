################################################################
### Code for paper:
### "A New Approach to Estimating the Production Function for Housing"
### by Dennis Epple, Brett Gordon, and Holger Sieg, forthcoming AER.
###
### This code is used to estimate the production function.
###
### Written by Brett Gordon, June 2008
### Comments should be sent to brg2114@columbia.edu.
################################################################

library(MASS)
library(splines)
library(deSolve)
library(lokern)
library(xtable)
library(stargazer)

#######################################################
### Read in the data
#######################################################

residential = read.csv("Pittsburgh_post1995.csv", header=TRUE)
residential.N = nrow(residential)

commercial = read.csv("Pittsburgh_Commercial_Properties.tsv", header=TRUE, sep="\t")
commercial.N = nrow(commercial)

###################
## SUMMARY TABLE ##
###################

.summary <- function (names, resid=TRUE) {
	# Returns the summary table in LaTeX
	if (resid==TRUE) {
		data <- residential
	} else {
		data <- commercial
	}

	df <- data.frame()

	for (name in names) {
		x <- data[[name]]
		d <- data.frame(
			"Mean"   = mean(x),
			"Median" = median(x),
			"SD"     = sd(x),
			"Min."   = min(x),
			"Max."   = max(x)
		)
		df <- rbind(df, d)
	}
	row.names(df) <- names
	return(xtable(df))
}

print.xtable(
	.summary(c("v", "pland", "lotarea", "tcog")),
	file="residential.tex"
)

print.xtable(
	.summary(c("v", "pland", "lotarea"), resid=FALSE),
	file="commercial.tex"
)

#############
## FIGURES ##
#############

attach(residential)

# Scatter plot that appeared in the original code but not in the published
# article
png(filename="scatter.png", height=350, width=600)
res <- lm(log(pland) ~ log(v))
plot(
	log(v), log(pland), 
	pch=20, cex=0.5,
	col=rgb(0, 0, 0, alpha=0.11),
	xlab="log(housing value per unit land)",
	ylab="log(land prices)",
	bty="n", 
	family="serif"

)
abline(res, col=rgb(0, 0, 1, alpha=0.5))
dev.off()

png(filename="density.png", height=350, width=600)
par(mfrow=c(1,2))
plot(
	density(v),
	main="",
	xlim=c(0,80),
	ylab="Density",
	xlab="housing value per unit land",
	bty="n",
	family="serif"
)
plot(
	density(pland),
	main="",
	xlim=c(0,20),
	ylab="Density",
	xlab="land prices",
	bty="n",
	family="serif"
)
dev.off()


############################
## OLS REGRESSION RESULTS ##
############################

lm.loglin = lm(log(pland) ~ log(v))
lm.lin = lm(pland ~ v - 1)
lm.quad = lm(pland ~ v + I(v^2) - 1)
lm.cub = lm(pland ~ v + I(v^2) + I(v^3) - 1)

stargazer(
	lm.loglin, lm.lin, lm.quad, lm.cub, 
	title="Regression Results", align=TRUE,
	omit.stat=c("LL","ser","f")
)

### Non-parametric estimation of r(v)

gkern.deriv = glkerns(
	v, pland, 
	deriv=1,
	n.out=300, 
	hetero=TRUE,
	bandwidth=sd(v)*1.1
)

gkern = glkerns(
	v, pland, 
	deriv=0,
	n.out=500, 
	hetero=TRUE,
	bandwidth=sd(v)*1.1
)

gkern.interp = interpSpline(
	gkern$x.out, gkern$est,
	na.action=na.omit,
	bSpline=TRUE
)

gkern.log = glkerns(
	log(v), log(pland), 
	deriv=0,
	n.out=300, 
	hetero=TRUE
)

gkern.log.interp = interpSpline(
	gkern.log$x.out, gkern.log$est,
	na.action=na.omit,
	bSpline=TRUE
)

gkern.log.fitted = predict(
	gkern.log.interp, log(v)
)

regs = data.frame(
	"Log-Linear"=lm.loglin$fitted.values, 
	"Linear"=log(lm.lin$fitted.values),
	"Quadratic"=log(lm.quad$fitted.values),
	"Cubic"=log(lm.cub$fitted.values), 
	"Log Kernel"=gkern.log.fitted$y
)


#######################################################
### Derive supply functions for each
#######################################################

vseq.n = 1000
vseq = seq(1,100,length.out=vseq.n)

### Calculate supply functions for quad and cubic cases
r = lm.quad$coefficients
p.quad = (vseq^r[1]) * exp(2*r[2]*(vseq-1))
s.quad = vseq/p.quad

r = lm.cub$coefficients
c = 2*r[2] + 1.5*r[3]
p.cub = (vseq^r[1]) * exp(2*r[2]*(vseq-1) + 1.5* r[3] * (vseq^2 - 1))
s.cub = vseq/p.cub

# Now create price grid given this upperbound
lb = min(range(p.quad)[1],range(p.cub)[1])
ub = min(range(p.quad)[2],range(p.cub)[2])
pgrid = seq(lb,ub,length.out=vseq.n)


### Calculate supply functions for linear and log-linear case s
alpha = as.numeric(
	exp(lm.loglin$coefficients[1]) * lm.loglin$coefficients[2]
)
beta = as.numeric(
	lm.loglin$coefficients[2] - 1
)

# Determine upper bound on p, since beta < 0 implies (c+log(p)) < 0
ec = exp(-alpha/beta)

### Linear case
k = lm.lin$coefficients[1]
s.lin = pgrid^((1-k)/k)

### Log-linear case
s.loglin = (1/pgrid)*((1 + (beta/alpha)*log(pgrid))^(1/beta))

### Solve ODE for kernel supply function

rprime = interpSpline(
	gkern.deriv$x.out, gkern.deriv$est, 
	na.action=na.omit,
	bSpline=TRUE
)

calc.supply.deriv = function(t, y, params) {
	return(list((y/t)*((1/(predict(rprime, y*t))$y)-1)))
}

xstart = c(1)
pgrid.kern = c(1, seq(1.001, ub, length.out=vseq.n-1))

ode.out = lsoda(
	xstart, pgrid.kern, calc.supply.deriv, 
	parms=1, 
	rtol=1e-5, 
	atol=1e-7
)

s.kern = data.frame(s=ode.out[ , 2], p=ode.out[ , 1])

s.funcs = data.frame(
	"Log-Linear"=s.loglin,
	"Linear"=s.lin,
	"Quadratic"=s.quad,
	"Cubic"=s.cub, 
	"Kernel"=s.kern$s
)

png(filename="supply.png", height=350, width=600)
op <- par(family = "serif")
plot(
	log(s.funcs$Log.Linear), log(pgrid),
	type="l",
	bty="n",
	xlab="log(Supply)",
	ylab="log(Price)",
	family="serif"
)
lines(log(s.funcs$Kernel), log(pgrid), lty=3)
lines(log(s.funcs$Linear), log(pgrid), lty=2)
legend(
	"topleft",
	c("log-linear", "linear", "kernel"), 
	lty=1:3,
	bty="n"
)
par(op)
dev.off()


#######################################################
### Production Function Estimation
#######################################################

v.lin.hat = s.lin*pgrid
v.loglin.hat = s.loglin*pgrid
v.quad.hat = s.quad*p.quad
v.cub.hat = s.cub*p.cub
v.kernel = s.kern$p*s.kern$s

lm.lin.fitted = lm.lin$coefficients[1]*v.lin.hat
lm.loglin.fitted = exp(lm.loglin$coefficients[1] + lm.loglin$coefficients[2]*log(v.loglin.hat))
lm.quad.fitted = lm.quad$coefficients[1]*v.quad.hat + lm.quad$coefficients[2]*v.quad.hat^2
lm.cub.fitted = lm.cub$coefficients[1]*v.cub.hat + lm.cub$coefficients[2]*v.cub.hat^2 + lm.cub$coefficients[3]*v.cub.hat^3
kernel.fitted = predict(gkern.interp,v.kernel)$y

m.lin = v.lin.hat-lm.lin.fitted 
m.loglin = v.loglin.hat-lm.loglin.fitted
m.quad = v.quad.hat-lm.quad.fitted
m.cub = v.cub.hat-lm.cub.fitted
m.kernel = v.kernel - kernel.fitted

##Normalize the production function for the kernel case
m.kernel = m.kernel + (1-m.kernel[1])

prod.funcs = data.frame("Log-Linear"=s.loglin,"Linear"=s.lin,"Quadratic"=s.quad, "Cubic"=s.cub, "Kernel"=s.kern$s)

m.grids = cbind(m.loglin, m.lin, m.quad, m.cub,m.kernel)

matplot(m.grids,prod.funcs, col=1:ncol(prod.funcs), lty=1:ncol(prod.funcs),lwd=c(2,2,2,2), type="l", 
	  main="Production Functions", xlab="m = (M/L)", ylab="Output", xlim=c(0,80), ylim=c(0,42))
legend(60, 12, names(prod.funcs), col=1:ncol(prod.funcs), lty=1:ncol(prod.funcs))

png(filename="production.png", height=350, width=600)
op <- par(family = "serif")
plot(
	log(prod.funcs$Log.Linear), log(m.lin),
	type="l",
	bty="n",
	xlab="log(m)",
	ylab="log(q)",
	family="serif"
)
lines(log(prod.funcs$Kernel), log(m.lin), lty=3)
lines(log(prod.funcs$Linear), log(m.lin), lty=2)
legend(
	"topleft",
	c("log-linear", "linear", "kernel"), 
	lty=1:3,
	bty="n"
)
par(op)
dev.off()


#######################################################
### Calculate supply elasticities
#######################################################

require(splines)
s.loglin.elas = calc.elas(pgrid,s.loglin)
s.lin.elas = calc.elas(pgrid,s.lin)
s.quad.elas = calc.elas(pgrid,s.quad)
s.cub.elas = calc.elas(pgrid,s.cub)
s.kern.elas = calc.elas(s.kern$p,s.kern$s)

s.elas = list("Linear"=s.lin.elas, "Log-linear"=s.loglin.elas, "Quad"=s.quad.elas,
		  "Cubic"=s.cub.elas, "Kernel"=s.kern.elas)
sapply(s.elas,mean)
sapply(s.elas,range)
sapply(s.elas,sd)

### Calculate demand elasticities

p.inv = function(v,s.interp){	
	pstar = rep(NA,length(v))
	s.inv = function(p,v,sfunc){ ((predict(sfunc,p))$y*p - v) }
	for(i in 1:length(v)){
		pstar[i] = (uniroot(s.inv,c(0.0001,100),v=v[i],sfunc=s.interp))$root
	}
	return(pstar)
}

s.lin.interp = interpSpline(pgrid,s.lin,na.action=na.omit,bSpline=TRUE)
s.loglin.interp = interpSpline(pgrid,s.loglin,na.action=na.omit,bSpline=TRUE)
s.quad.interp = interpSpline(p.quad,s.quad,na.action=na.omit,bSpline=TRUE)
s.cub.interp = interpSpline(p.cub,s.cub,na.action=na.omit,bSpline=TRUE)
s.kern.interp = interpSpline(s.kern$p,s.kern$s,na.action=na.omit,bSpline=TRUE)

p.lin.inv = p.inv(v,s.lin.interp)
p.loglin.inv = p.inv(v,s.loglin.interp)
p.quad.inv = p.inv(v,s.quad.interp)
p.cub.inv = p.inv(v,s.cub.interp)
p.kern.inv = p.inv(v,s.kern.interp)

elas.lin = calc.average.elas(p.lin.inv,s.lin.interp)
elas.loglin = calc.average.elas(p.loglin.inv,s.loglin.interp)
elas.quad = calc.average.elas(p.quad.inv,s.quad.interp)
elas.cub = calc.average.elas(p.cub.inv,s.cub.interp)
elas.kern = calc.average.elas(p.kern.inv,s.kern.interp)

p.loglin.inv.median = p.inv(median(v),s.loglin.interp)
elas.loglin.median = calc.average.elas(p.loglin.inv.median,s.loglin.interp)

#######################################################
### Calculate the elasticity of substitution for land vs. non-land inputs
#######################################################

n = 50
L = lotarea
m = v - exp(lm.loglin$coefficients[1] + lm.loglin$coefficients[2]*log(v))
M = m*L
Lseq = seq(min(L),max(L),length.out=n)
Mseq = seq(min(M),max(M),length.out=n)

#######################################################
### Functions
#######################################################

pfunc.CD = function(k,Lseq,Mseq){
	Q = matrix(NA,nrow=length(Lseq),ncol=length(Mseq))
	A = 1/((1-k)^(1-k))
	for(i in 1:(length(Lseq))){ for(k in 1:(length(Mseq))){
		Q[i,k] = A*(Lseq[i]^k)*(Mseq[k]^(1-k))
	}}
	return(Q)
}

Q.CD = pfunc.CD(lm.lin$coefficients[1],Lseq,Mseq)
persp(Lseq,Mseq,Q.CD,theta=-20,phi=45,xlim=c(500000,1200000),ylim=c(20000,300000))

pfunc.full = function(Lseq,Mseq,prod.interp){
	Q = matrix(NA,nrow=length(Lseq),ncol=length(Mseq))
	for(i in 1:(length(Lseq))){
		Q[i,] = pfunc.full2(Lseq[i],Mseq,prod.interp)
	}
	return(Q)
}

pfunc.full2 = function(lval,mval,prod.interp){
	return(lval*(predict(prod.interp,mval/lval)$y))
}

pfunc.full.L = function(lval, mval,prod.interp, ord = 1){
	temp = function(x,mval){ return(pfunc.full2(x,mval,prod.interp))}
	return(fdiff(lval, temp, h=NULL, order=ord, accur=4, mval))
}

pfunc.full.M = function(mval,lval,prod.interp, ord = 1){
	temp = function(x,lval){ return( pfunc.full2(lval,x,prod.interp))}
	return(fdiff(mval, temp, h=NULL, order=ord, accur=4, lval))
}

pfunc.full.LM = function(Lval, Mval, prod.interp){
	h = (.Machine$double.eps)^(1/3)*abs(min(Lval,Mval))

	f.xy = pfunc.full2(Lval+h, Mval+h, prod.interp)
	f.x_y = pfunc.full2(Lval+h, Mval-h, prod.interp)
	f._xy = pfunc.full2(Lval-h, Mval+h, prod.interp)
	f._x_y = pfunc.full2(Lval-h, Mval-h, prod.interp)

	return(((f.xy-f.x_y)-(f._xy-f._x_y))/(4*h^2))
}

calc.elas.sub = function(Lseq, Mseq, prod.interp){
	n.L = length(Lseq)
	n.M = length(Mseq)
	Q.elas.sub = matrix(NA,nrow=n.L, ncol=n.M)
	for(i in 1:n.L){
		for(k in 1:n.M){
			F.L = pfunc.full.L(Lseq[i],Mseq[k],prod.interp, ord=1)
			F.M = pfunc.full.M(Mseq[k],Lseq[i],prod.interp, ord=1)
			F.LL = pfunc.full.L(Lseq[i],Mseq[k],prod.interp, ord=2)
			F.MM = pfunc.full.M(Mseq[k],Lseq[i],prod.interp, ord=2)
			F.LM = pfunc.full.LM(Lseq[i],Mseq[k],prod.interp)

			numer = F.L*F.M*(F.L*Lseq[i] + F.M*Mseq[k])
			denom = Lseq[i]*Mseq[k]*((2*F.L*F.M*F.LM)  - ((F.L^2)*F.MM) - ((F.M^2)*F.LL))
	
			Q.elas.sub[i,k] = numer/denom
		}
	}
	return(Q.elas.sub)
}

require(splines)
prod.lin.interp = interpSpline(m.lin,s.lin,na.action=na.omit,bSpline=TRUE)
prod.loglin.interp = interpSpline(m.loglin,s.loglin,na.action=na.omit,bSpline=TRUE)
prod.quad.interp = interpSpline(m.quad,s.quad,na.action=na.omit,bSpline=TRUE)
prod.cub.interp = interpSpline(m.cub,s.cub,na.action=na.omit,bSpline=TRUE)

require(R.basic)

Q = pfunc.full(Lseq, Mseq,prod.loglin.interp)
Q.vec = as.vector(Q)
L.vec = as.vector(rep(Lseq,n))
M.vec = rep(NA, length(Q.vec))
for(i in 1:n){ for(k in 1:n){
	M.vec[(i-1)*n + k] = Mseq[i]
}}
plot3d(L.vec,M.vec,Q.vec)

Q.elas.loglin = calc.elas.sub(Lseq, Mseq, prod.loglin.interp)
Q.elas.lin = calc.elas.sub(Lseq, Mseq, prod.lin.interp)
Q.elas.quad = calc.elas.sub(Lseq, Mseq, prod.quad.interp)
Q.elas.cub = calc.elas.sub(Lseq, Mseq, prod.cub.interp)



##############################################
### Simulate supply and production functions to calcalate standard errors

sup.loglin = function(grid,parms){
	alpha = as.numeric(exp(parms[1])*parms[2])
	beta = as.numeric(parms[2] - 1)
	return((1/grid)*((1 + (beta/alpha)*log(grid))^(1/beta)))
}

prod.loglin = function(q,s){
	v.hat = s.loglin*pgrid
	lm.loglin.fitted = exp(lm.loglin$coefficients[1] + lm.loglin$coefficients[2]*log(v.hat))
	m.loglin = v.hat-lm.loglin.fitted
	return(list(m=m.loglin,q=s))
}

sim.sup.prod = function(mod,grid,sfunc,pfunc,nsim){
	r = mod$coefficients[,"Estimate"]
	sd = mod$coefficients[,"Std. Error"]
	r.sim = matrix(NA,nrow=nsim,ncol=length(r))
	for(i in 1:length(r)){
		r.sim[,i] = rnorm(nsim,mean=r[i],sd=sd[i])
	}
	
	s.sd = rep(NA,length(grid))
	s.mat = matrix(NA,nrow=length(grid),ncol=nsim)
	for(k in 1:length(grid)){
		for(i in 1:nsim){
			s.mat[k,i] = sfunc(grid[k],r.sim[i,])		
		}
		s.sd[k] = sd(s.mat[k,])
	}

	p.sd = rep(NA,length(grid))
	p.mat = matrix(NA,nrow=length(grid),ncol=nsim)
	for(i in 1:nsim){
		## Calculate the prod function for this supply function
		pfunc.temp = pfunc(grid,s.mat[,i])
		p.mat[,i] = pfunc.temp$q
	}
	for(k in 1:length(grid)) p.sd[k] = sd(p.mat[k,])
	
	return(list(s.sd=s.sd,p.sd=p.sd))
}

sim.funcs = sim.sup.prod(summary(lm.loglin),pgrid,sup.loglin,prod.loglin,100)

s.loglin.ub = s.loglin + 2*sim.funcs$s.sd
s.loglin.lb = s.loglin - 2*sim.funcs$s.sd
prod.loglin.lb = s.loglin - 2*sim.funcs$p.sd
prod.loglin.ub = s.loglin + 2*sim.funcs$p.sd

# Supply and production plots
plot(s.loglin,pgrid,type="l",main="Loglin Supply Function with \n 95% Confidence Band",xlab="Supply",ylab="Price", xlim=c(0,40))
lines(s.loglin.lb,pgrid,col="red",lty=2)
lines(s.loglin.ub,pgrid,col="red",lty=2)

plot(m.loglin,s.loglin, type="l",main="Loglin Production Function with \n 95% Confidence Band",xlab="m",ylab="q", xlim=c(0,80))
lines(m.loglin, prod.loglin.ub, col="red",lty=2)
lines(m.loglin, prod.loglin.lb, col="red",lty=2)

# Supply and production plots in log-log scale
plot(log(s.loglin),log(pgrid),type="l",main="Log-linear Supply Function with \n 95% Confidence Band",
					xlab="log(Supply)",ylab="log(Price)", xlim=c(0,log(40)), lwd=2)
lines(log(s.loglin.lb),log(pgrid),col="black",lty=2)
lines(log(s.loglin.ub),log(pgrid),col="black",lty=2)

plot(log(m.loglin),log(s.loglin), type="l",main="Log-linear Production Function with \n 95% Confidence Band",
					xlab="log(m)",ylab="log(q)", xlim=c(0,log(80)), lwd=2)
lines(log(m.loglin), log(prod.loglin.ub), col="black",lty=2)
lines(log(m.loglin), log(prod.loglin.lb), col="black",lty=2)



###########################################
# A function to compute highly accurate first- and second-order derivatives
# From Fornberg and Sloan (Acta Numerica, 1994, p. 203-267; Table 1, page 213)
fdiff <- function(x, fun, h=NULL, order=1, accur=4, ...) {
	macheps <- .Machine$double.eps
	
	if (order==1) {
		if(is.null(h)) h <- macheps^(1/3)* abs(x)
		ifelse (accur==2, w <- c(-1/2,1/2), w <- c(1/12,-2/3, 2/3,-1/12))
		ifelse (accur==2, xd <- x + h*c(-1,1), xd <- x + h*c(-2,-1,1,2))
		return(sum(w*fun(xd,...))/h)
	}
	else if (order==2) {
		if(is.null(h)) h <- macheps^(1/4)* abs(x)
		ifelse (accur==2, w <- c(1,-2,1), w <- c(-1/12,4/3,-5/2,4/3,-1/12))
		ifelse (accur==2, xd <- x + h*c(-1,0,1), xd <- x + h*c(-2,-1,0,1,2))
		return(sum(w*fun(xd,...))/h^2)
	}
}

# y is a data frame, x is a matrix
plot.with.legend = function(x,y,col=1:ncol(x),lty=1:ncol(x),xlim=NULL,ylim=NULL,
			titles=names(y),lposx=(max(x)-max(x)/4),lposy=max(y)/3,main="",xl="",yl="",...){
	matplot(x,as.matrix(y),col=col, lty=lty,lwd=c(2,2,2,2), type="l",
		 main=main,xlab=xl,ylab=yl,xlim=xlim,ylim=ylim, ...)
	legend(lposx,lposy, titles, col=col, lty=lty)
}

### Calculate the derivative of a function over a given range of values
calc.deriv = function(grid,fun,interp.obj=NULL,high.density=FALSE){
      	require(splines)
	fun.interp = NULL
	if(! is.null(interp.obj)){
		fun.interp = interp.obj
	}
	else if(!high.density && length(grid)>200){
		interp.grid = rep(NA,200)
		interp.grid[1] = grid[1]
		interp.grid[200] = grid[200]
		interp.grid = sample(grid,198)
		fun.grid = fun[interp.grid]
		fun.interp = interpSpline(interp.grid,fun.grid,na.action=na.omit,bSpline=FALSE)
	}
	else{
		fun.interp = interpSpline(grid,fun,na.action=na.omit,bSpline=FALSE)
	}
	d = rep(NA,length(grid))
	for(i in 1:length(grid)){
		d[i] = fdiff(grid[i],function(x){(predict(fun.interp,x))$y})
	}
	return(list(x=grid,d=d))
}

calc.elas = function(grid,sfunc,bspline=TRUE){
	require(splines)
	s.interp = interpSpline(grid,sfunc,na.action=na.omit,bSpline=bspline)
	s.elas = rep(NA,length(grid))
	for(i in 1:length(grid)){
		s.elas[i] = (fdiff(grid[i],function(x){(predict(s.interp,x))$y}))*(grid[i]/sfunc[i])
	}
	return(s.elas)
}

### pa - actual set of prices corresponding to the data
### s - supply function
calc.average.elas = function(pa,s.interp,bspline=TRUE){
	require(splines)
	s.pred = predict(s.interp,pa)
	s.elas = rep(NA,length(pa))
	for(i in 1:length(pa)){
		s.elas[i] = (fdiff(pa[i],function(x){(predict(s.interp,x))$y}))*(pa[i]/s.pred$y[i])
	}
	return(s.elas)
}
