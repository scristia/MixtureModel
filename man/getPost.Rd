\name{getPost}
\alias{getPost}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
gets posterior distribution and performs model selection
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getPost(r, kmin = 2, kmax = 5, S = 200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{kmin}{
%%     ~~Describe \code{kmin} here~~
}
  \item{kmax}{
%%     ~~Describe \code{kmax} here~~
}
  \item{S}{
%%     ~~Describe \code{S} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (r, kmin = 2, kmax = 5, S = 200) 
{
    r <- r[!is.na(r)]
    n <- length(r)
    loglik <- c()
    bic <- c()
    kpost <- list()
    for (k in kmin:kmax) {
        nn <- rep(length(r)/k, k)
        alpha <- rep(1, k)
        mus <- kmeans(r, centers = k, nstart = 15)$centers
        mus <- sort(mus)
        print(k)
        print(mus)
        post <- gibbs(r = r, nn = nn, K = k, alpha = alpha, mu0 = mus, 
            S = S)
        loglik[k - kmin + 1] <- loglik.normmix(r, post)
        bic[k - kmin + 1] <- -2 * loglik[k - kmin + 1] + (3 * 
            k - 1) * log(n)
        kpost[[k - kmin + 1]] <- post
    }
    results <- data.frame(K = seq(kmin, kmax), `log-likelihood` = loglik, 
        BIC = bic)
    print(results)
    return(kpost)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
