#####################################################################################################################################
# probability mass function (pmf)
### ddstable(k,alpha,lambda)

#####################################################################################################################################

#' @title The discrete stable distribution: formulae for the probabilities (density)
#' @description Computes the value of the formulae for the probabilities (density) of a discrete stable distribution \code{DS(alpha,lambda)}, by combining the explicit and fast asymptotic formulae.
#' @param x a vector of non-negative integer quantiles, \code{k>=0}
#' @param alpha tail index parameter \code{alpha} in the interval= \eqn{(0, 1]}
#' @param lambda positive location parameter \code{lambda>0}
#' @importFrom stats dpois
#' @return Returns the value of the formulae for the probabilities (density) of \code{DS(alpha,lambda)}.
#' @export
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{CHRISTOPH1998243}{dstabledist}

#' @examples
#' ddstable(c(0,1,2,100),1,lambda=1)#This is Poisson with lambda=1
#' dpois(c(0,1,2,100),1)#Checking with dpois
#' ddstable(c(0,1,2,100),0.5,lambda=1) # tail is heavier
#' ddstable(c(0,1,2,3,6,100),0.5,lambda=3) # change in location

ddstable<-function(x,alpha,lambda){
   ddstable_exact<-function(k,alpha,lambda){ #k>=0,
    p=0
    p[1]=exp(-lambda)
    if(k==0) {return(exp(-lambda))}
    if(p[1]==0) {return(rep(0,k+1))}
   #  if (alpha==1) {return(dpois(k,lambda,log=FALSE));break} #can be switched to dpois
    v=0
    ds=1:k
    for (s in 1:k){v[s]=(alpha-ds[s]+1)/ds[s]}
    for (i in 0:(k-1)){  sum_m=0:i;
    for (m in 0:i){ sum_m[m+1]<-p[i-m+1]*(m+1)*((-1)^m)*prod(v[1:(m+1)])}
    abi1<-sum(sum_m)
    rm(sum_m)
    p[i+2]=(lambda/(i+1))*abi1
    #print(c(i+2,p[i+2])) #UNCOMMENT for large k if want to see workflow on screen
    }
    return(p)
  }
  ddstable_tail1<-function(df,alpha,lambda,m){
    P=0
      for (i in 1:length(df)){ #if(k==0) {return(exp(-lambda))}
      k=df[i]
      M<-vector('numeric',m) #m is number of sums, floor((alpha+1)/alpha), but recommended <=15
      for (j in 1:m){
        M[j]<-(((-1)^(j+1))/j)*(lambda^j)*sin(alpha*j*pi)*gamma(alpha*j+1)*k^(-alpha*j-1)
      }
      P[i]=(1/pi)*sum(M,na.rm=NA)#;print(P)
      {if (P[i]<0) {m=m-1; {if (m>=1) {P[i]<-ddstable_tail1(k,alpha,lambda,m)[1]} else {P[i]='NA';return(P)}}}}
      {if (P[i]>1) {m=m-1; {if (m>=1) {P[i]<-ddstable_tail1(k,alpha,lambda,m)[1]} else {P[i]='NA';return(P)}}}}
      {if (P[i]=='NA') {return(P)}}
    }
    return(P)
  }
  ddstable_tail<-function(df,alpha,lambda){
    if (alpha==1) {return(dpois(df,lambda,log=FALSE))} # switch to Poisson
    P=0
     m=floor((alpha+1)/alpha)
    # if(m==0) {return ('m>0!')}
    if (m>20) {m<-20}
      for (i in 1:length(df)){ #if(k==0) {return(exp(-lambda))}
      k=df[i]
      M<-vector('numeric',m) #m is number of sums, floor((alpha+1)/alpha), but recommended <=15
      for (j in 1:m){
        M[j]<-(((-1)^(j+1))/j)*(lambda^j)*sin(alpha*j*pi)*gamma(alpha*j+1)*k^(-alpha*j-1)
      }
      P[i]=(1/pi)*sum(M,na.rm=NA)#;print(P)
      {if (P[i]<0) {m=m-1; {if (m>=1) {P[i]<-ddstable_tail1(k,alpha,lambda,m)[1]} else {P[i]='NA';return(P)}}}}
      {if (P[i]>1) {m=m-1; {if (m>=1) {P[i]<-ddstable_tail1(k,alpha,lambda,m)[1]} else {P[i]='NA';return(P)}}}}
      {if (P[i]=='NA') {return(P)}}
    }
     print(m)
    return(P)
  }
  kv=x #  "kv=vector of x values"
  if (!all(kv>=0)||!all(kv==round(kv)))  {return('Invalid arguments. Required: k is a vector of non-negative integer quantiles')}
  if ( alpha<=0 || alpha>1 )
  {return('Invalid arguments. Required:  alpha  in (0,1] ')}
    if (lambda<=0)
  {return('Invalid arguments. Required: lambda>0')}
  tr=1500 #tr is threshold from exact to asymptotic.
  #ddstable_exact works for any k, but due to computational power it can be very slow for k~>1500.
  pk=0;  pk1=0; pk2=0;
  df1=kv[kv<=tr];
  df2=kv[kv>tr];
  if (length(df1)>0) {kk=max(df1);  pk_all=ddstable_exact(kk,alpha,lambda);  pk1=pk_all[df1+1]}
  if (length(df2)>0) {pk2=ddstable_tail(df2,alpha,lambda)}
  pk=c(pk1,pk2)
  if (length(df1)==0) pk=pk2;
  if (length(df2)==0) pk=pk1;
  return(pk)}









