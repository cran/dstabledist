
#####################################################################################################################################
#Simulations
### rdstable(n,alpha,lambda)


#' @title The discrete stable distribution:  random generation
#' @description Generates random variates from a discrete stable distribution \code{DS(alpha,lambda)}.
#' @param n number of random values to return.
#' @param alpha tail index parameter \code{alpha} in the interval= \eqn{(0, 1]}
#' @param lambda positive location parameter \code{lambda>0}
#' @importFrom stabledist rstable
#' @importFrom stats rpois dpois
#' @importFrom Rdpack reprompt
#' @return returns random variates from \code{DS(alpha,lambda)}. A warning is displayed for invalid parameter values.
#' @export
#' @references
#' \insertRef{DEVROYE1993349}{dstabledist}


#' @examples
#' rdstable(10,alpha=1,lambda=1) #this is Poisson
#' rdstable(10,alpha=0.5,lambda=1) # heavier tail more prone to extremes
#' rdstable(10,alpha=0.1,lambda=1) # heavier tail more prone to extremes

rdstable<-function (n,alpha,lambda=1){
  if (n<1 || alpha<=0 || alpha>1 || lambda<=0)
  {return('Invalid arguments: n should be a positive integer, n>=1, alpha should be between 0 and 1, alpha in (0,1], and  lambda should be positive, lambda>0.')}
  if(alpha==1) {return(stats::rpois(n,lambda))}
  #define scale function
  g_scale <- function(a) {
    iu <- complex(real=0, imaginary=1)
    return(abs(1 - iu * tan(pi * a / 2)) ^ (-1 / a))
  }
  #define continuous stable scale gamma g
  g <- g_scale(alpha)
  # positive stable random variable
  ps <- stabledist::rstable(n=n,alpha=alpha,beta=1,gamma=g,delta=0,pm=1)
  # simulate discrete stable vector
  DS <- rep(NA,n)
  l <- ps*(lambda^(1/alpha))
  for (i in 1:length(l)){DS[i]=rpois(1,l[i])}
  return(DS)}








