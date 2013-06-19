dnormmix <-
function(r, mixture, burnin=1:50, log=FALSE) {
	mix.mean <- apply(mixture[["means"]][-burnin, ], 2, mean)
	mix.prec <- apply(mixture[["precs"]][-burnin, ], 2, mean)
	p <- apply(mixture[["P"]][-burnin, ], 2, mean)
	K <- length(p)
	## Find likelihood for each component
	lik.comp <- function(r, comp) {
		p[comp]*dnorm(r, mean=mix.mean[comp], sd=1/sqrt(mix.prec[comp]))
	}
	## Create likelihood array with likelihoods of from each component
	liks <- sapply(1:K, lik.comp, r=r)

	d <- rowSums(liks)
	if(log) {
		d <- log(d)
	}
	return(d)
}
