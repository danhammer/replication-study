################################################################
### Code for paper:
### "A New Approach to Estimating the Production Function for Housing"
### by Dennis Epple, Brett Gordon, and Holger Sieg, forthcoming AER.
###
### This code uses the implementation of the 
### estimator in Hausman, Newey, Ichimura, and Powell (Journal of Econometrics, 1991)
### and applies it to the data set in the paper on Pittsburgh housing.
###
### Written by Brett Gordon, June 2008
### Comments should be sent to brg2114@columbia.edu.
################################################################


source("HNIP_estimator.R")

INST.HAS.CONSTANT = TRUE
BASE.EQ.CONSTANT = FALSE

data = read.csv("Pittsburgh_post1995.csv", header=TRUE)
N = nrow(data)

attach(data)

hnip.estimates <- function(K) {

	p.cont = K+1  ### dimensions of the instruments
	p.dummy = length(unique(muni)) - 1   ##
	P = p.cont + p.dummy  #### Dimension of the instruments

	### Create design matrix of dummy variables based on municipalities
	Q.d = factor(muni)
	Q.df = model.matrix(~ Q.d - 1)

	### Estimate with the dummies as IVs
	beta.est = HNP.estimator(
		K, p.cont, tcog,
		(as.matrix(Q.df))[,2:(p.dummy+1)], 
		v, pland
	)

	return(beta.est)
}

print(hnip.estimates(2))
print(hnip.estimates(3))

#############################################
### Bootstrap for standard errors
#############################################
B = 1000
beta.hat = matrix(NA, nrow=B, ncol=K + 1)  ## add one for chi.stat

for(b in 1:B){

	index.b = sample(N, replace=T)
	v.b = v[index.b]
	pland.b = pland[index.b]
	muni.b = muni[index.b]
	tcog.b = tcog[index.b]

	p.dummy = length(unique(muni.b)) - 1   ##
	P = p.cont + p.dummy  #### Dimension of the instruments

	Q.d.b = factor(muni.b)
	Q.df.b = model.matrix(~ Q.d.b - 1)

	beta.hat[b,] = HNP.estimator(
		K, p.cont, tcog.b,
		(as.matrix(Q.df.b))[,2:(p.dummy+1)],
		v.b, pland.b
	)

	if(b %% 50 == 0){
		print(paste("b = ", b)) ## print iterations
	}
}

est = sapply(
	1:K, 
	function(i){ 
		return(
			round(c(mean(beta.hat[,i]), sqrt(var(beta.hat[,i])) ), digits=8)) 
	}
)

png(filename="hist.png", height=400, width=400)
hist(
	beta.hat[ , 1], 
	main="",
	col="gray",
	xlab="",
	breaks=20,
	border="white",
	ylab="",
	family="serif"
)
dev.off()


detach(data)










