loglik.normmix <-
function(r, mixture, K, burnin=1:50) {
	loglike <- dnormmix(r, mixture, K, burnin=1:50, log=TRUE)
	return(sum(loglike))
}
