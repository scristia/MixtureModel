loglik.normmix <-
function(r, mixture) {
	loglike <- dnormmix(r, mixture, log=TRUE)
	return(sum(loglike))
}
