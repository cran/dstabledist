#####################################################################################################################################
# probability generating function (pgf)
### pgdstable(z,alpha,lambda)


#' @title The discrete stable distribution: probability generating function
#' @description Computes probability generating function of a discrete stable distribution \code{DS(alpha,lambda)}.
#' @param z argument of probability generating function, \code{z} in the interval= \eqn{[-1, 1]}.
#' @param alpha tail index parameter \code{alpha} in the interval= \eqn{(0, 1]}.
#' @param lambda positive location parameter \code{lambda>0}.
#' @return Returns value of probability generating function of \code{DS(alpha,lambda)}. A warning is displayed for invalid parameter values.
#' @export
#' @examples
#' pgdstable(c(-1,0,1),0.5,1)
#' pgdstable(c(-1,0,1),1,1) #This is Poisson
#'  curve(pgdstable(x,1,lambda=1), c(-1,1),col=1,ylab='prob. gen. fun.',xlab='z')
#'  curve(pgdstable(x,0.5,lambda=1), c(-1,1),col=2,add=TRUE)
#'  curve(pgdstable(x,0.2,lambda=1), c(-1,1),col=4,add=TRUE)
#' legend('topleft',legend=c(1,0.5,0.1), col=c(1,2,4), lty = 1, title='alpha')
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{Steutel_1979}{dstabledist}

#Probability generating function (pgf)
### pgdstable(z,alpha,lambda): function for probability generating function (pgf) of discrete stable distribution DS(alpha,lambda)
pgdstable<-function (z,alpha,lambda=1){
  if (!all(abs(z)<=1) || alpha<=0 || alpha>1 || lambda<=0)
    {return('Invalid arguments. Required:  all(abs(z)<=1),  alpha in (0,1] and lambda>0')}
  arg<- (1-z)
  arg2<-arg^alpha
  pgf<-exp(-lambda*(arg2))
  return(pgf)
}






