#$Id: /R/local/age groups/mid.R 108 2005-11-17T12:12:56.154042Z mpiktas  $

#Distributions used for fitting
#Generalized gamma distribution density function. Note: this variant is without
#normalizing constant, so it is really not a pdf. Use dggamma and dggamma.t for the
#"real" pdf. For sampling from generalized gamma use rgg function.
#dgg <- function(x,theta) x^theta[1]*exp(theta[2]*x+theta[3]*x^2+theta[4]*x^3)
dgg <- function(x,pv) dggamma.t(x,c(0,pv))
dggamma <- function(x,eta,t1,t2,t3,t4) x^t1*exp(t2*x+t3*x^2+t4*x^3+eta)
dggamma.t <- function(x,pv)dggamma(x,pv[1],pv[2],pv[3],pv[4],pv[5])
dgg.d <- function(x,theta) {
 cbind(log(x),x,x^2,x^3)*dgg(x,theta)   
}
norm <- function(x,p=2)(sum(abs(x)^p))^(1/p)
#Utility functions for generalized gamma distribution
gg.cnt <- function(theta)theta[1]>-1 &theta[4]<0
#Calculates mean of generalized gamma distribution
gg.mn <- function(t,maxa) {
    xdgg <- function(x,pv)x*dggamma.t(x,pv)
    integrate(xdgg,0,maxa,pv=t)$val
}
gg.tftmg <- function(theta,maxa,...) {
    tr <- rgg(100000,theta,maxa)
    gn <- Gini(tr$var)
    
    ex <- gg.mn(c(tr$eta,theta),maxa)
    c(ex,gn)
}

sm.tft<-function(theta,...){
 
 ex <- theta[2]*sm.ex(theta)
# gn <- sm.gini(theta)

 xsm <- function(x,theta)  qsm(x,theta)
  pk <- function(xx,theta) {
        integrate(xsm,xx[1],xx[2],subdivisions=1000,theta=theta)$value
  }

   q<-c(seq(0,1,by=0.1))
#q<-c(0,0.2,0.4,0.6,0.8,1)


   n<-10 #deciles

   md<-c()
   for(i in 1:n){
     md<-c(md,pk(c(q[i],q[i+1]),theta)/0.1)
   }


c(md) #!!!!!!!!!!!!!!!!!!
}



#Singh-Maddala distribution
dsm <- function(x,theta)exp(log(theta[1]) +log(theta[3]) +log(x)*(theta[1] - 1) - theta[1]*log(theta[2]) -log(1 + (x/theta[2])^theta[1])*(1+theta[3]))
psm <- function(q,theta)1-(1+(q/theta[2])^theta[1])^(-theta[3])
qsm <- function(p,theta)theta[2]*((1-p)^(-1/theta[3])-1)^(1/theta[1])
rsm <- function(n,theta)qsm(runif(n),theta)
dsm.d <- function(x,theta) {
    tt <- exp(theta[1]*(log(x/theta[2])))

    d1 <- 1/theta[1]+log(x/theta[2])-(theta[3]+1)*log(x/theta[2])*(1-1/(1+tt))
    d2 <- -theta[1]/theta[2]+theta[1]*(theta[3]+1)/theta[2]*(1-1/(1+tt))
    d3 <- 1/theta[3]-log(1+tt)
    cbind(d1,d2,d3)*dsm(x,theta)
}


#Utility functions for Singh-Maddala distribution
#second parameter is scale parameter, so it can be calculated from other two and the
#mean of distribution
sm.ex <- function(t) {                                # b=EX/SM.ex(c(a,q)); EX=b*SM.ex(c(a,q))
    if (length(t)==2) t <- c(t[1],0,t[2])
    beta(1+1/t[1],t[3]-1/t[1])/beta(1,t[3]) 
}
#provide the triplet of SM parameters from mean and first and third parameters
sm.bchange <- function(EX,t) {
    c(t[1],EX/SM.ex(t),t[3])
}
#These constraints ensure, that SM will have mean.
sm.cnt <- function(theta) theta[1]*theta[3]>1

#Leave dots, for compatibility reasons                                        
sm.tftmg <- function(theta,...) {
    ex <- theta[2]*sm.ex(theta)
    gn <- sm.gini(theta)
    c(ex,gn)
}
sm.tftmm <- function(theta,...) {
    ex <- theta[2]*sm.ex(theta)
    md <- qsm(0.5,theta)
    c(ex,md)
}
sm.gini <- function(theta)if((theta[3]*theta[1])>1)1-beta(theta[3],2*theta[3]-1/theta[1])/beta(theta[3]-1/theta[1],2*theta[3]) else -1

sm.tftq<-function(theta,...){
  f<-function(x,theta) x*dsm(x,theta)
  q1<-integrate(f,0,qsm(0.2,theta),theta=theta)$val
  q2<-integrate(f,qsm(0.2,theta),qsm(0.4,theta),theta=theta)$val
  q3<-integrate(f,qsm(0.4,theta),qsm(0.6,theta),theta=theta)$val
  q4<-integrate(f,qsm(0.6,theta),qsm(0.8,theta),theta=theta)$val
  q5<-integrate(f,qsm(0.8,theta),Inf,theta=theta)$val

c(q1,q2,q3,q4,q5)
}


#Exponential and lognormal mixture distribution
delnm <- function(x,theta) theta[1]*dlnorm(x,theta[2],theta[3])+(1-theta[1])*dexp(x,1/theta[4])
pelnm <- function(q,theta) theta[1]*plnorm(q,theta[2],theta[3])+(1-theta[1])*pexp(q,1/theta[4])
relnm <- function(n,theta) {
    res <- sample(c(1,0),n,rep=TRUE,prob=c(theta[1],1-theta[1]))
    r1 <- res==1
    res[r1] <- rlnorm(sum(r1),theta[2],theta[3])
    res[!r1] <- rexp(sum(!r1),1/theta[4])
    res
}
delnm.d <- function(x,theta) {
    t1 <- log(x)-theta[2]
    t2 <- dlnorm(x,theta[2],theta[3])

    d1 <- dlnorm(x,theta[2],theta[3])-dexp(x,1/theta[4])

    d2 <- theta[1]*t1/theta[3]^2*t2
    d3 <- (-1/theta[3]+t1^2/theta[3]^3)*t2*theta[1]
    d4 <- (-1/theta[4]+x/theta[4]^2)*dexp(x,1/theta[4])*(1-theta[1])
    cbind(d1,d2,d3,d4)
}

elnm.cnt <- function(theta) theta[1]>0 & theta[1]<1 & theta[4]>0 &theta[3]>0
elnm.tftmg <- function(theta,...) {
    ex <- theta[1]*exp(theta[2]+theta[3]^2/2)+(1-theta[1])*theta[4]
    int.fun <- function(x,theta)(pelnm(x,theta)-1/2)*x*delnm(x,theta)
    gn <- integrate(int.fun,0,Inf,theta=theta)$val*2/ex
    c(ex,gn)
    
}

dgm <- function(x,theta)dgamma(x,theta[1],rate=theta[2])
pgm <- function(q,theta)pgamma(q,theta[1],rate=theta[2])
qgm <- function(p,theta)qgamma(p,theta[1],rate=theta[2])
rgm <- function(n,theta)rgamma(n,theta[1],rate=theta[2])

dgm.d <- function(x,theta) {
    d1 <- +log(theta[2])-digamma(theta[1])+log(x)
    d2 <- theta[1]/theta[2]-x
    cbind(d1,d2)*dgm(x,theta)
}

gm.cnt <- function(theta)sum(theta<=0)==0

#The function is here for compatibility reasons. Also some assumptions are being
#made taking in mind, that we are working with income data.

makeclasses <- function(a,cellno,tail=NULL) {
    #Classes are the intervals , according to which observations are grouped.
    #The union of classes should span the range of the observed variable and
    #different classes should be disjoint.
    #a - the information about intervals of the classes, can be:
    #   1. The two column array, where each row represents class. The left edges of 
    #class intervals are in first column, the right edges are in the second column.
    #   2. The vector with the left edges of the intervals
    #   3. The vector with the right edges of the intervals.
    #   4. The vector with the all edges of intervals. 
    #cellno - the number of cells (the number of classes)
    #tail - NULL, for default tails: left tail defaults to 0, right tail to Inf.
    #depending on the a, can be left tail for case 3, right tail for case 2 and
    #vector with both left and right tails for case 1.

    if(is.null(dim(a)) || min(dim(a))==1) {
        #If a is column or row vector transform it to simple vector.
        if(!is.null(dim(a))) a <- a[,,drop=TRUE]
        ni <- length(a)
        if(ni==cellno+1) {
            #We have case 4.
            classes <- cbind(a[-ni],a[-1])
            
            if(!is.null(tail)) {
                if(length(tail)==1) warning("Inapropriate tail, left and right tails should be specified. Ignoring the tail option")
                else {

                    if(diff(tail)<0) {
                        warning("Left tail is bigger than right tail, exchanging them")
                        tail <- range(tail)
                    }
                    if (tail[1]<classes[1,2]) classes[1,1] <- tail[1]
                else warning("Left tail is higher than the right edge of the first class interval. Changing nothing")
                if (tail[2]>classes[cellno,1]) classes[cellno,2] <- tail[2]
                else warning("Right tail is lower than the left edge of the last class interval. Changing nothing")
                }
            }
        }
        else {
            if(ni!=cellno) stop("The number of classes and cells should be the same")
            if(a[1]==0) {
                #We have case number 2
                if(!is.null(tail)) {
                    tail <- max(tail)
                    if(tail<a[ni]) {
                        warning("Right tail is lower than the left edge of the last class interval. Reverting to default")
                        tail <- Inf
                    }
                }
                else tail <- Inf
                classes <- cbind(a,c(a[-1],tail))
            }
            else {
                #We have case number 3
                if(!is.null(tail)) {
                    tail <- min(tail)
                    if(tail>a[1]) {
                        warning("Left tail is higher than the right edge of the first class interval. Reverting to default")
                        tail <- 0
                    }
                }
                else tail <- 0
                classes <- cbind(c(tail,a[-ni]),a)
            }
        }
    }
    else {
        #We have case number 1
        if(sum(range(dim(a))!=c(2,cellno))>0)stop("Either number of classes is incorrect or, the array has more than two collumns(rows)")
        if(dim(a)[2]!=2) {
            warning("There are more than two columns. Assuming that there are two rows")
            a <- t(a)
        }
        classes <- a
        if(!is.null(tail)) {

            if(length(tail)==1) warning("Inapropriate tail, left and right tails should be specified. Ignoring the tail option")
            else {
                if(diff(tail)<0) {
                    warning("Left tail is bigger than right tail, exchanging them")
                    tail <- range(tail)
                }
                if (tail[1]<classes[1,2]) classes[1,1] <- tail[1]
                else warning("Left tail is higher than the right edge of the first class interval. Changing nothing")
                if (tail[2]>classes[cellno,1]) classes[cellno,2] <- tail[2]
                else warning("Right tail is lower than the left edge of the last class interval. Changing nothing")
            }
        }
    }
    colnames(classes) <- NULL
    classes
}
#Frequency Based Maximum Likelihood.
#Given frequencies produces fit, for supplied density function
#Truncates the distribution, check if this is really apropriate for you!

fbml <- function(a,o,
                 p,pdf,
                 constr=function(theta)TRUE,
                 cdf=NULL,
                 normalize=TRUE,
                 tail=NULL,
                 method="BFGS",hessian=TRUE,
                 silent=FALSE,
                 ...) {
    #a - lower or upper bounds of intervals.
    #o - number of people in given interval.
    #p - initial values for theta.
    #pdf - probability density function we are fitting. The first argument should
    #be quantile, the second should be  a numeric vector containing all the
    #parameters of distribution.
    #constr - function returning TRUE for valid pdf parameters
    #cdf - cumulative density function. Not necessary. The same naming conventions
    #as for pdf applies.
    #normalize - if TRUE, the pdf is given with the precision of constant
    #tail - argument for makeclasses function, which is applied to argument a, before
    #calculations
    #method - parameter for optim function.
    #hessian - parameter for optim function.
    #silent - if FALSE, information about fit and fitting progress is produced.
    #... - additional parameters for optim

    #if(sum(abs(a)==Inf)>0) stop("Infinite length intervals are not accepted")
    if(sum(a<0)>0) stop("Negative incomes are not supported")
    ea <- as.list(match.call(call=sys.call(sys.parent(n=1))))
    #Get default arguments listed for future reference.
    ea$constr <- constr
    ea$method <- method
    ea$hessian <- hessian
    ea$silent <- silent
    ea$cdf <- eval(cdf)
    ea$pdf <- eval(pdf)
    ea$p <- p
   
    l <- length(o)
    classes <- makeclasses(a,l,tail)

    N <- sum(o)
    Ok <- o                                       

    PDF <- function(x,theta) pdf(x,theta)
    pk <- function(xx,theta) {
        if(is.null(cdf)) {
            res <- try(integrate(PDF,xx[1],xx[2],theta=theta,stop.on.error=F)$val,silent=TRUE) #The error we encounter is infinite integral, for debugging do silent=F
            if(class(res)=="try-error") Inf
            else res

        }
        else cdf(xx[2],theta)-cdf(xx[1],theta)
    }
    pteta <- function(theta) {
        e <- apply(classes,1,pk,theta=theta)
        if(normalize) e <- e/sum(e)
        e
    }
    C2 <- lgamma(N+1)-sum(lgamma(Ok+1))
    minuslogl <- function(theta) {
        if (!constr(theta)) Inf
        else -sum(Ok*log(pteta(theta)))-C2
    }
    if (!silent) {
        cat("\nSTART")
        cat("\n ",minuslogl(p),"\n")
    }

    fit <- try(optim(p,minuslogl,method=method,hessian=hessian,...))
    
    if(class(fit)=="try-error") {
        #The usual error is with BFGS method, retry with Nelder-Mead
        fit <- optim(p,minuslogl,...)
    }
    if (!silent) print(fit)
    E <- pteta(fit$p)*N
    sE <- E>1
    if(sum(!sE)>0)warning("Adjusting Chi2")
    Chi2 <- sum((E[sE]-Ok[sE])^2/E[sE])
    L0 <- minuslogl(p)
    L1 <- minuslogl(fit$p)
    #truncate distribution at max(a)
    rtail <- max(classes)

    if(is.null(cdf)) {
        normc <- integrate(PDF,0,rtail,theta=fit$p)$value
        CDF <- NULL
    }
    else {
        normc <- cdf(rtail,fit$p)
        CDF <- function(x) {
            res <- cdf(x,fit$p)/normc
            res[res>1] <- 1
            res
        }
          
    }
    PDF <- function(x) {
        res <- pdf(x,fit$p)/normc
        res[x>=rtail] <- 0
        res
    }
    
    output <- list(a=a,o=o,e=E,
                   coefficients=fit$p,fit=fit,L0=L0,L1=L1,Chi2=Chi2,
                   pdf=PDF,cdf=CDF,CALL=ea,classes=classes)
    class(output) <- "fbml"
    output
    
}

hist.fbml <- function(x,
                      new=TRUE,
                      pE=TRUE,
                      pChi2=TRUE) {

    #x - fbml output
    #new - if TRUE, produce new plot
    #pE - if TRUE plot expected values as red bars
    #pChi2 - if TRUE print Chi2 value as title

   
   o <- x$o
   cl <- x$classes
   
   ni <- length(o)
   a <- cl[,2]
   ma <- max(x$a)
   if(ma!=Inf && ma>cl[ni,1]) {
       a[ni] <- ma
   }
   else a[ni] <- cl[ni,1]+0.5*(cl[ni,1]-cl[1,1])
       
   
   
   y <- o/sum(o)
   y <- y/diff(c(0,a))
   my <- max(y)
   if (!is.null(x$pdf)) {
       xx <- seq(0,max(a),len=101)
       yy <- x$pdf(xx)
       myy <- max(yy)
       if(myy==Inf) myy <- my
       my <- max(my,myy)
   }
   
   if(new)x11()
   par(mar=c(2,2,2,0))
   plot.new()


   plot.window(range(a,0),range(my,0),main='diufcaie')
   axis(1,at=my.atx <- seq(0,max(a),len=10),lab=formatC(my.atx,dig=2))
   axis(2)
   if(pChi2) {
       ttl <- ""
       if(!is.null(x$Chi2))ttl <- round(x$Chi2,3)
       if(!is.null(x$aChi2))ttl <- paste(ttl,round(x$aChi2,3)," ")
       title(ttl)

   }
   if(pE) {
       if(!is.null(x$e)) {
           e <- x$e
           e <- e/sum(e)
           e <- e/diff(c(0,a))
           rect(c(0,a[-ni]),0,a,e,col='red',border=FALSE)
       }
   }
   if(!is.null(x$pdf)) lines(xx,yy)
   rect(c(0,a[-ni]),0,a,y)
   
}
#A wrapper function for fbml. Fits multiple groups at once. It is assumed that
#every group will have the same distribution, but with different parameters. It
#is also assumed that intervals are the same for all groups (as it is the case
#usually)
fbmlgroup <- function(a,o,
                 p,pdf,
                 constr=function(theta)TRUE,
                 cdf=NULL,
                 normalize=TRUE,
                 tail=NULL,
                 method="BFGS",hessian=TRUE,
                 silent=FALSE,
                 ...) {
    #a - lower or upper bounds of intervals. 
    #o - number of people in given interval. A matrix, where each column represents
    #different distribution.
    #p - initial values for theta. Can be a matrix with different initial values for
    #different groups
    #pdf - parameter for fbml.
    #constr - parameter for fbml.
    #cdf - parameter for fbml.
    #method - parameter for optim function.
    #hessian - parameter for optim function.
    #silent - if FALSE, information about fit and fitting progress is produced.
    #... - additional parameters for optim

    ng <- dim(o)[2]
    if(is.null(dim(p))) p <- array(p,dim=c(length(p),ng))
    else if(dim(p)[2]!=ng) {
        warning("Number of starting values does not equal to number of groups")
        p <- array(p,dim=c(dim(p)[1],ng))
    }
    o <- as.array(o)
    
    output <- list()
    for(i in 1:ng) output <- c(output,list(fbml(a,o[,i],p[,i],pdf,constr,cdf,normalize,tail,method,hessian,silent,...)))
    class(output) <- "fbmlgroup"
    output
}

hist.fbmlgroup <- function(x,mfr=c(4,4),pE=TRUE,pChi2=TRUE) {
    par(mfrow=mfr)
    invisible(lapply(x,hist,pE=TRUE,pChi2=TRUE,new=FALSE))
}

#Quah GMM, based on paper  "One Third of the World?”?s Growth
#and Inequality"  by Danny Quah. Fits distribution,  having its
#characteristics. The idea is that, the characteristcs really are functionals
#of underlying distribution, thus we can use GMM for fitting.

qgmm <- function(x,p,
                 tftheta,
                 constr=function(theta)TRUE,
                 omega=rep(1,length(x)),
                 method="BFGS", hessian=TRUE,
                 silent=FALSE,
                 ...) {
    #x - the characteristics of distribution. Usually mean and some other, Gini for
    #example.
    #p - starting values for optimization.
    #tftheta - the functionals of underlying distribution giving the same
    #characteristics x, but depending on theta.
    #constr - constraints for theta.
    #omega - weight matrix, only diagonal can be given.
    #method - optim parameter.
    #hessian - optim parameter.
    #silent - if FALSE, information about fit and fitting progress is produced.
    #... - additional parameters for optim 
    
    ea <- as.list(match.call(call=sys.call(sys.parent(n=1))))
    #Get default arguments listed for future reference.
    #ea$tftheta <- eval(tftheta)
    ea$constr <- eval(constr)
    ea$method <- method
    ea$hessian <- hessian
    ea$silent <- silent
    if(is.null(dim(omega))) omega <- diag(omega)
    ea$omega <- omega
    ea$p <- p
    och <- chol(omega)
    if(dim(omega)[1]!=length(x) | diff(dim(omega))!=0 ) stop("Wrong weights!")
    qfun <- function(theta,...) {
        if(constr(theta)) {
            tt <- tftheta(theta,...)-x
            sum((och%*%tt)^2)

        }
        else Inf
    }

    if(!silent) {
        cat("\n THE START\n")
        cat(qfun(p))
        cat("\n")
    }
    fit <- try(optim(p,qfun,method=method,hessian=hessian,...))
    
    if(class(fit)=="try-error") {
        #The usual error is with BFGS method, retry with Nelder-Mead
        fit <- optim(p,qfun,method="Nelder-Mead",hessian=hessian,...)
        ea$method <- "Nelder-Mead"
    }
    if (!silent) print(fit)
    output <- list(coefficients=fit$p,x=x,tftheta=eval(tftheta,...),fit=fit,CALL=ea)
    class(output) <- "qgmm"
    output
}

qgmmgroup <- function(x,p,
                 tftheta,
                 constr=function(theta)TRUE,
                 omega=rep(1,dim(x)[1]),
                 method="BFGS", hessian=TRUE,
                 silent=FALSE,
                 ...) {
    #x - the characteristics of distribution. Usually mean and some other, Gini for
    #example. Should be array (matrix) with group data in columns.
    #p - starting values for optimization, can be array.
    #rtail - right tail of distribution at which it is truncated.
    #tftheta - the functionals of underlying distribution giving the same
    #characteristics x, but depending on theta.
    #constr - constraints for theta.
    #omega - weight matrix, only diagonal can be given. For the moment, only
    #one weight matrix for all age groups can be given.
    #method - optim parameter.
    #hessian - optim parameter.
    #silent - if FALSE, information about fit and fitting progress is produced.
    #... - additional parameters for optim
    ng <- dim(x)[2]
    if(is.null(dim(p))) p <- array(p,dim=c(length(p),ng))
    else if(dim(p)[2]!=ng) {
        warning("Number of starting values does not equal to number of groups")
        p <- array(p,dim=c(dim(p)[1],ng))
    }
    output <- list()
    for(i in 1:ng) output <- c(output,list(qgmm(x[,i,drop=TRUE],p[,i,drop=TRUE],tftheta,constr,omega,method,hessian,silent,...)))
    class(output) <- "qgmmgroup"
    output
}

as.fbml.qgmm <- function(qout,a,o,pdf,cdf=NULL,tail=NULL) {


       
    l <- length(o)
    classes <- makeclasses(a,l,tail)

    N <- sum(o)  #total Number of observations
    Ok <- o                                       

    constr <- qout$CALL$constr
    PDF <- function(x,theta) pdf(x,theta)
    pk <- function(xx,theta) {
        if(is.null(cdf)) {
            res <- try(integrate(PDF,xx[1],xx[2],theta=theta,stop.on.error=F)$val,silent=TRUE) #The error we encounter is infinite integral, for debugging do silent=F
            if(class(res)=="try-error") Inf
            else res

        }
        else cdf(xx[2],theta)-cdf(xx[1],theta)
    }
    pteta <- function(theta) {
        e <- apply(classes,1,pk,theta=theta)
        e <- e/sum(e)
        e
    }
    C2 <- lgamma(N+1)-sum(lgamma(Ok+1))
    minuslogl <- function(theta) {
        if (!constr(theta)) Inf
        else -sum(Ok*log(pteta(theta)))-C2
    }
    

    L0 <- minuslogl(qout$CALL$p)
    L1 <- minuslogl(coef(qout))
    e <- pteta(coef(qout))*N
    Chi2 <- sum((e-o)^2/e)
    #if(max(a)!=qout$rtail)warning("Truncation points differ")
    rtail <- qout$rtail
    if(is.null(cdf)) {
        normc <- integrate(PDF,0,rtail,theta=coef(qout))$value
        CDF <- NULL
    }
    else {
        normc <- cdf(rtail,coef(qout))
        CDF <- function(x) {
            res <- cdf(x,coef(qout))/normc
            res[res>1] <- 1
            res
        }
          
    }
    PDF <- function(x) {
        res <- pdf(x,coef(qout))/normc
        res[x>=rtail] <- 0
        res
    }
    CALL <- qout$CALL
    CALL$pdf <- eval(pdf)
    CALL$cdf <- eval(cdf)
    output <- list(a=a,o=o,e=e,coefficients=coef(qout),fit=qout$fit,L0=L0,L1=L1,Chi2=Chi2,pdf=PDF,cdf=CDF,CALL=CALL,classes=classes)
    class(output) <- "fbml"
    output
}
as.fbmlgroup.qgmmgroup <- function(qout,a,o,pdf,cdf=NULL,tail=NULL) {
    output <- list()
    for(i in 1:dim(o)[2]) output <- c(output,list(as.fbml.qgmm(qout[[i]],a,o[,i],pdf,cdf,tail)))
    class(output) <- "fbmlgroup"
    output
}

obre <- function(a,o,
                 p,
                 pdf,grad,
                 gamma=1,c=Inf,
                 tail=NULL,
                 normalize=FALSE,
                 constraints=function(theta)TRUE,
                 silent=FALSE,
                 control=list()
                 ) {

    #a - lower or upper bounds of intervals.
    #o - number of people in given interval.
    #p - initial values for theta.
    #pdf - probability density function we are fitting. The first argument should
    #be quantile, the second should be  a numeric vector containing all the
    #parameters of distribution.
    #grad - gradient of probability density function with respect to parameters.
    #Should be function with the same arguments as pdf, but should return a matrix
    #with the number of rows the same as length of the first argument, and the number
    #of columns the same as the length of the second argument.
    #normalize - if TRUE, the pdf is given with the precision of constant
    #tail - argument for makeclasses function, which is applied to argument a, before
    #calculations
    #constraints - function returning TRUE for valid pdf parameters
    #eta - convergence criteria
    #silent - if FALSE, information about fit and fitting progress is produced.

    if(sum(a<0)>0) stop("Negative incomes are not supported")
    ea <- as.list(match.call(call=sys.call(sys.parent(n=1))))
    ea <- lapply(ea,eval)
   

    con <- list(trace=0,reltol=sqrt(.Machine$double.eps),abstol=sqrt(.Machine$double.eps),maxit=300)
    con[names(control)] <- control
    if(con$abstol!=con$reltol) {
        con$abstol <- min(con$abstol,con$reltol)
    }
    
    l <- length(o) #Number of intervals
    classes <- makeclasses(a,l,tail)

    dth <- length(p) #dimension of theta
    N <- sum(o)  #total Number of observations
    
    PDF <- function(x,theta) pdf(x,theta)

    gPDF <- function(x,theta)grad(x,theta)

    #integrates given function with parameters theta, on interval (xx[1],xx[2])
    intc <- function(xx,theta,fun) {
        res <- try(integrate(fun,xx[1],xx[2],theta=theta,stop.on.error=F)$val,silent=TRUE) #The error we encounter is infinite integral, for debugging do silent=F
        if(class(res)=="try-error") Inf
        else res
    }

    #returns probabilities of classes for candidate theta
    ktheta <- function(theta) {
        e <- apply(classes,1,intc,fun=PDF,theta=theta)
        e
    }
    #returns derivatives of probabilities of classes for candidate theta
    dktheta <- function(theta) {
        hifun <- function(i) {
            gri <- function(x,theta)gPDF(x,theta)[,i]
            apply(classes,1,intc,fun=gri,theta=theta)
        }
        sapply(1:dth,hifun)
    }

    pj <- o/sum(o)
    #Calculate initial A
    a <- rep(0,dth)
    kt <- ktheta(p)
    dkt <- dktheta(p)
    etta <- 1

    if(normalize) {
        etta <- sum(kt)
        deta <- apply(dkt,2,sum)
        kt <- kt/etta
        dkt <- (dkt-outer(kt,deta))/etta
        
    }
    if(abs(sum(dkt))>1e-5) stop("Sum of gradient of class probabilities did not equal zero. Check your gradient function or use normalize=TRUE")
    
    dlkt <- sweep(dkt,1,kt,"/")
    nullkt <- kt==0
    nullno <- sum(nullkt)

    if(nullno>0) {
        warning("Some of classes has null probabilities")
        dlkt[nullkt,] <- 0
    }

    A <- apply(dkt,1,function(xx)outer(xx,xx))
    A <- sweep(A,2,kt,"/")

    if (nullno>0) A[,nullkt] <- 0
    
    A <- apply(A,1,sum)
    #since A is symmetric, we do not worry about how it was converted to vector
    A <- matrix(A,ncol=dth)
    #A <- solve(A)

    theta <- 0
    delta <- p
    i <- 1
    
    while(norm(delta)>con$abstol && i<con$maxit) {
        theta <- theta+delta
        if(con$trace>0) {
            cat("\n Iteration ",i," convergence criteria: ", norm(delta), "\n ")
            if(con$trace>1)cat("Parameter values:" ,theta,"\n")
        
        }
        
        kt <- ktheta(theta)
        nullkt <- kt==0
        nullno <- sum(nullkt)
        #browser()
        dkt <- dktheta(theta)
    
        if(normalize) {
            etta <- sum(kt)
            deta <- apply(dkt,2,sum)
            kt <- kt/etta
            dkt <- (dkt-outer(kt,deta))/etta
        }
        if(abs(sum(dkt))>1e-5) warning("Sum of gradient of class probabilities did not equal zero. Check your gradient function or use normalize=TRUE")
        
        dlkt <- sweep(dkt,1,kt,"/")
        
        if (nullno>0) {
            warning("Some of the classes have null probabilities")
            dlkt[nullkt,] <- 0
        }
        ddif <- dkt-outer(kt,a)
        W <- apply(dlkt,1,function(xx)min(1,c/norm(solve(A,(xx-a)))))
        if (nullno>0) W[nullkt] <- 0

        AA <- apply(cbind(ddif*W,dlkt),1,function(xx)outer(xx[1:dth],xx[dth+1:dth]))
        AA <- apply(AA,1,sum)
        AA <- matrix(AA,ncol=dth,byrow=TRUE)

        aa <- sweep(dkt,1,W,"*")
        aa <- apply(aa,2,sum)
        aa <- aa/sum(kt*W)
        
        delta <- sweep(ddif,1,W*(pj/kt)^gamma,"*")
        if(nullno>0)delta[nullkt,] <- 0
        delta <- apply(delta,2,sum)
        delta <- solve(AA,delta)

        if(con$trace>2)cat("Gradient:",delta,"\n")
        if(!constraints(theta+delta))stop("Process did not converge, try choosing different starting values")
        a <- aa
        A <- AA
        
        i <- i + 1

    }
    skt <- N*kt>0.5
    if (sum(!skt)>0)warning("Adjusting Chi2")
    
    aChi2 <- N*sum((kt[skt]-pj[skt])^2/kt[skt])
    Chi2 <- N*sum((kt-pj)^2/kt)
    if(!silent) {
        cat("\n Convergence results:\n")
        cat("Parameter values: ", theta, "\n")
        cat("Convergence criteria: ",norm(delta),"\n")
        cat("Number of iterations: ",i,"\n")
        cat("Goodness of fit, Chi2: ",aChi2,"\n")
    }


    if(normalize) outpdf <- function(x)PDF(x,theta=theta)/etta
    else outpdf <- function(x)PDF(x,theta)

    
    output <- list(a=a,o=o,e=N*kt,
                   coefficients=theta,
                   Chi2=Chi2,
                   aChi2=aChi2,
                   pdf=outpdf,CALL=ea,classes=classes,
                   fitt=list(normc=etta,delta=delta,a=a,A=A,W=W,eta=con$abstol))
    class(output) <- "obre"
    output
}

hist.obre <- function(x,
                      new=TRUE,
                      pE=TRUE,
                      pChi2=TRUE) {
    hist.fbml(x,new=new,pE=pE,pChi2=pChi2)
}
