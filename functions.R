library(ggplot2)
library(reshape2)
library(forecast)
library(urca)
library(plyr)
library(TSA)
library(gridExtra)
library(TTR)
library(manipulate)

######################################

spalva="dodgerblue4"
col1="#FAA537"
col2="#5F88A1"
col3="#57854E"
colt="white"
col4="steelblue4"
ec = c("#8B6F66","#7C6D96","#0078BF","#206779","#57854E","#FDCF41","#E07A3F","#BD4F5C")

######################################



press=function(x,str,plot=F){
  
  oldx=x
  
  xform=x/mean(x)-1
  xform=xform/str
  
  x=(xform+1)*mean(x)
  
  if(plot==T){
    matplot(cbind(oldx,x),lty=1,pch=20,type="l",main=paste("Sqeezing strength:",str,sep=" "),ylab="")
    points(oldx,pch=20)
    points(x,pch=20,col=2)
  }
  
  return(x)
}


level=function(x,lvl,plot=F){
  
  oldx=x
  olvl=mean(oldx)
  
  xform=x/mean(x)
  x=xform*lvl
  
  if(plot==T){
    matplot(cbind(oldx,x),lty=1,pch=20,type="l",main=paste("Old lvl",round(olvl,2),">>> New lvl",lvl,sep=" "),ylab="")
    points(oldx,pch=20)
    points(x,pch=20,col=2)
  }
  
  return(x)
}


sm=function(x,df=6,plot=F){
  
  oldx=x
  
  x=smooth.spline(x, df=df)$y
  
  if(plot==T){
    matplot(cbind(oldx,x),lty=1,pch=20,type="l",main=paste("Smoothing strength:",df,sep=" "),ylab="")
    points(oldx,pch=20)
    points(x,pch=20,col=2)
  }
  
  return(x)
}

# ff=lininterpol(cps.eval,dfor)
# 
# x=cps.eval
# xframe=dfor

lininterpol=function(x,xframe){
  x=as.numeric(unlist(x))
  outx=c()
  
  for(k in 1:length(x)){
    
    minnr=max(which(x[k]>=xframe[,2]))
    maxnr=min(which(x[k]<=xframe[,2]))
    
    if(maxnr==Inf)maxnr=minnr
    
    s1=as.numeric(xframe[minnr,2])
    s2=as.numeric(xframe[maxnr,2])
    
    ratio=(x[k]-s1)/(s2-s1)
    if(is.na(ratio)) ratio=1
    if(ratio==Inf) ratio=1
    
    outx[k]=ratio*(xframe[maxnr,1]-xframe[minnr,1])+xframe[minnr,1]
  }
  return(outx)
}



over100=function(xframe){
  
  xframe=apply(xframe,2,function(x){
    on=c()
    while(any(x>100)){
      
      on=unique(c(on,which(x>100)))
      over=sum(x[on]-100)
      x[on]=100
      x[-on]=x[-on]+over/length(x[-on])
      
    }
    return(x)}) 
  
  return(xframe)
}


critbound=function(x,crit=1.5){
  (100/x)^(1/crit)/100
  (1/x)^(1/crit)/2
  
}


critbound=function(x,crit=3){
  kritai=seq(4,crit,length.out=100)
  xsai=seq(1,100,by=1)
  
  minnr=max(which(x>=xsai))
  smash=kritai[minnr]
  
  1/x/smash
}


pos=function(x){
  out=F
  pos=sum(diff(x)>=0)
  neg=sum(diff(x)<0)
  
  if(pos>=neg)out=T
  return(out)
}



# iternorm ----------------------------------------------------------------
# iteratyviam normavimui parasyta funkcija

iternorm=function(matout,ttn,plottest=F,plottest2=T){
  
  omatout=matout
  
  if(plottest){
    matplot(t(matout),type="l",lty=1)
    matplot(apply(matout,2,function(x)x/sum(x)),type="l",lty=1)
  }
  
  shemat=t(apply(matout,2,function(x)x/sum(x)))
  shematd=apply(shemat,2,function(x)x[2:26]/x[1:25])
  
  protype=(shematd*ttn[2:26]/ttn[1:25])
  opr=protype
  
  if(plottest){
    matplot(protype,type="l",main="Starting point",lty=1)
    abline(h=1,col=4,lty=2,lwd=1)
    #     matplot(shemat,type="l",main="Starting point",lty=1)
    #     abline(h=1,col=4,lty=2,lwd=1)
    
  }
  
  
  
  totsign=sign(diff(ttn))
  
  testos=c()
  for(k in 1:25){ 
    sg=totsign[k]
    if(sg==1 | sg==0){
      if(all(protype[k,]>=1)) testos[k]=T else testos[k]=F
    }
    
    if(sg==-1){
      if(all(protype[k,]<1)) testos[k]=T else testos[k]=F
    }
  }
  
  
  
  iter=1
  while(!all(testos)){
    #     while(any(protype<1)){
    
    les=(1:25)[!testos]
    
    for(j in les){
      
      x1=shemat[j,]
      x2=shemat[j+1,]
      
      or= x2/x1*ttn[j+1]/ttn[j]
      r=or
      i=1
      t=seq(0,1,by=0.01)
      combo=x2
      
      #       print(j)
      
      TTT=F
      if(totsign[j]>=0){if(any(r<1)) TTT=T}
      if(totsign[j]<0){if(any(r>1)) TTT=T}
      
      while(TTT){
        
        #     cat(t[i],"\n")
        
        combo=x1*t[i] + x2 * (1-t[i])
        
        r=combo/x1*ttn[j+1]/ttn[j]
        
        TTT=F
        if(totsign[j]>=0){if(any(r<1)) TTT=T}
        if(totsign[j]<0){if(any(r>1)) TTT=T}
        
        #         plot(r,ylim=range(or),type="l",main=paste("Iteration",as.character(i-1)),
        #              xlab=paste("lambda=",t[i]))
        #         abline(h=1,col=4,lty=2)
        if(i>length(t))  {cat("NEPAVYKO SUKONVERGUOTI ITERATYVIAI","\n"); break}
        i=i+1
      }
      
      shemat[j+1,]=combo
    }
    
    shematd=apply(shemat,2,function(x)x[2:26]/x[1:25])
    
    protype=(shematd*ttn[2:26]/ttn[1:25])
    
    if(plottest){
      matplot(protype,type="l",main=paste("Iteration",as.character(iter)),lty=1,ylim=range(opr))
      abline(h=1,col=4,lty=2,lwd=1)
    }
    
    testos=c()
    for(k in 1:25){ 
      sg=totsign[k]
      if(sg==1 | sg==0){
        if(all(protype[k,]>=1)) testos[k]=T else testos[k]=F
      }
      
      if(sg==-1){
        if(all(protype[k,]<1)) testos[k]=T else testos[k]=F
      }
    }
    
    iter=iter+1
  }
  
  
  matoutre=shemat
  for(i in 1:26){
    matoutre[i,]=shemat[i,]*ttn[i]*10
  }
  
  matoutre=t(matoutre)
  
  if(plottest2){  
    par(mfrow=c(2,2))
    matplot(t(normal.normal(matout,ttn)),type="l",lty=1,main="BEFORE (NORMAL.NORMAL)",xlab="Year",ylab="% device in use")
    lines(ttn,col="firebrick3",lwd=1.5,pch=8,type="o")
    matplot(t(matoutre),type="l",lty=1,main="AFTER",
            ylab=paste("Iterations Nr.",as.character(iter-1)),xlab="Year")
    lines(ttn,col="firebrick3",lwd=1.5,pch=8,type="o")
    
    matplot(normal.normal(matout,ttn),type="l",lty=1,main="BEFORE",xlab="Decile",ylab="% device in use")
    matplot(matoutre,type="l",lty=1,main="AFTER",
            ylab=paste("Iterations Nr.",as.character(iter-1)),xlab="Decile")
    par(mfrow=c(1,1))
  }
  
  return(matoutre)
}


# monotonic ---------------------------------------------------------------
# >= 

monotonic=function(matout,ttn){
  
  pos=which(diff(ttn)>=0)
  neg=which(diff(ttn)<0)
  
  out=matout
  for(i in 1:10){
    x=as.numeric(out[i,])
    for(j in pos[pos!=1]){
      if(x[j]<x[j-1]) x[j]=x[j-1]
    }
    out[i,]=x
  }
  
  for(i in 1:10){
    x=as.numeric(out[i,])
    for(j in neg[neg!=1]){
      if(x[j]>x[j-1]) x[j]=x[j-1]
    }
    out[i,]=x
  }
  
  
  return(out) 
}

# matplot(t(matout),type="o",pch=20)
# matplot(t(monotonic(matout,ttn)),type="l")


monotonic.income=function(yframe){
  
  out=yframe
  for(i in 1:10){
    x=as.numeric(out[i,])
    for(j in 2:dim(yframe)[2]){
      if(!is.na(x[j])){
        if(x[j]<x[j-1]) x[j]=x[j-1]
      }
    }
    out[i,]=x
  }
  
  return(out) 
}


# mv  := monotonic vector
mv <- function(x,way="auto",split=0){
  
  x <- as.numeric(unlist(x))
  ox <- x
  
  #auto
  if(way =="auto"){
    if(sum(diff(x)>0)>=sum(diff(x)<0)) way  <- "up" else way <- "down"
    #     print(paste("Method AUTO used:","the way ->",way,"<- chosen"))
  }
  
  #up
  if(way == "up"){
    for(i in 2:length(x)){
      if(x[i]<x[i-1]) x[i] <- x[i-1]
    }
    
    bad <- which(x > ox)
    bad <- sort(unique(c(bad,bad+split)))
    bad <- setdiff(bad,length(x)+1:(split+1))
    x[bad] <- NA
    
    if(is.na(x[length(x)])){
      x[length(x)] <- x[max(which(!is.na(x)))] 
    }
  }
  
  #down
  if(way == "down"){
    for(i in 2:length(x)){
      if(x[i]>x[i-1]) x[i] <- x[i-1]
    }
    
    bad <- which(x < ox)
    bad <- sort(unique(c(bad,bad+split)))
    bad <- setdiff(bad,length(x)+1:(split+1))
    x[bad] <- NA
    
    if(is.na(x[length(x)])){
      x[length(x)] <- x[max(which(!is.na(x)))] 
    }
  }
  
  x <- MASplineVector(x)
  
  return(x) 
}

# nutemptas total ---------------------------------------------------------

nutemptas.total=function(matout,tt){
  
  for(y in 2005:2030){
    ff=matout[y-2004]
    ff=ff*as.numeric(tt[as.character(y)])/mean(as.numeric(unlist(ff)))
    matout[y-2004]=ff
  }
  
  
  kryptis=apply(matout,1,function(x){
    min(which(diff(as.numeric(x))<0))-1
  })
  stablenr=min(kryptis[kryptis>5])
  if(stablenr==Inf)stablenr=2
  
  
  for(y in (2005+stablenr):2030){
    ff=matout[2005+stablenr-2004-1]
    ff=ff*as.numeric(tt[as.character(y)])/mean(as.numeric(unlist(ff)))
    matout[y-2004]=ff
  }
  
  matout[stablenr]=matout[(stablenr-1)]*0.5+matout[(stablenr+1)]*0.5
  matout[stablenr]=matout[stablenr]*as.numeric(tt[as.character(stablenr+2005-1)])/
    mean(as.numeric(unlist(matout[stablenr])))
  
  return(matout)
}

# nutemptas total ---------------------------------------------------------

nutemptas.total.brutal=function(matout,tt){
  
  matout[c(1:4,6:26)]=matout[,5]
  
  for(y in 2005:2030){
    ff=matout[y-2004]
    ff=ff*as.numeric(tt[as.character(y)])/mean(as.numeric(unlist(ff)))
    matout[y-2004]=ff
  }
  
  return(matout)
}

# normal.normal ---------------------------------------------------------

normal.normal=function(matout,ttn){
  
  for(i in 1:26){
    matout[,i]=level(matout[,i],lvl=ttn[i],plot=F)
  }
  
  return(matout)
}

# smooth.mat ---------------------------------------------------------

smooth.mat=function(matout,df){
  
  for(i in 1:10){
    matout[i,]=sm(as.numeric(matout[i,]),df=df,plot=F)
  }
  
  return(matout)
}

# fordec ------------------------------------------------------------------
# function, which counts everything, what you need (well, only all possesions)


fordec=function(xhist,yhist,yframe,totline,plottest=F,plottest3=F,xframe){
  
  d=data.frame(cbind(y=yhist,x=xhist))
  ld=d
  
  model.np <- npreg(x ~ y ,regtype = "lc",
                    gradients = TRUE,
                    data = d,ckertype="gaussian",ckerorder=2)
  #       model.np2 <- npreg(x ~ y ,regtype = "ll",
  #                          gradients = TRUE,
  #                          data = d,ckertype="gaussian",ckerorder=2)
  if(plottest){  
    plot(model.np)
    points(d,col=4,pch=20)
  }
  #       plot(model.np2)
  
  plot(yframe[,1],xframe[,1],xlim=range(yframe),ylim=range(xframe),pch=20)
  points(yframe[,2],xframe[,2],xlim=range(yframe),ylim=range(xframe),col=2,pch=20)
  points(yframe[,3],xframe[,3],xlim=range(yframe),ylim=range(xframe),col=3,pch=20)
  points(yframe[,4],xframe[,4],xlim=range(yframe),ylim=range(xframe),col=4,pch=20)
  points(yframe[,5],xframe[,5],xlim=range(yframe),ylim=range(xframe),col=i,pch=20)
  
  
  
  
  
  
  
  
  # cps.eval=data.frame(y=seq(min(d$y),max(d$y),by=100))
  cps.eval=data.frame(y=seq(min(d$y),max(d[-10,]$y),by=100))
  foremaster=predict(model.np, newdata = cps.eval)
  #   foremaster=sm(foremaster,df=8,plot=T)
  
  
  
  models=list()
  
  models[[1]]=try(nls(x ~ SSlogis(y, alpha, xmid, scale),data = ld),silent=T)
  models[[2]]=try(nls(x ~ SSweibull(y,Asym,Drop,lrc,pwr),data = ld),silent=T)
  models[[3]]=try(nls(x ~ SSasympOff(y, Asym, lrc, c0),data = ld),silent=T)
  models[[4]]=try(nls(x ~ SSasympOrig(y, Asym, lrc),data = ld),silent=T)
  models[[5]]=try(nls(x ~ SSbiexp(y, A1, lrc1, A2, lrc2),data = ld),silent=T)
  models[[6]]=try(nls(x ~ SSfpl(y, A, B, xmid, scal),data = ld),silent=T)
  models[[7]]=try(nls(x ~ SSgompertz(y, Asym, b2, b3),data = ld),silent=T)
  models[[8]]=try(nls(x ~ SSmicmen(y, Vm, K),data = ld),silent=T)
  models[[9]]=try(nls(x ~ SSasymp(y, Asym, R0, lrc),data = ld),silent=T)
  
  names(models)=c("SSlogis","SSweibull","SSasympOff","SSasympOrig",
                  "SSbiexp","SSfpl","SSgompertz","SSmicmen","SSasymp")
  
  
  alo=data.frame(y=seq(min(yframe,na.rm=T),max(yframe,na.rm=T),by=100))
  
  fitmat=data.frame(matrix(ncol=9,nrow=dim(alo)[1]))
  for(i in 1:9){
    if(class(models[[i]])=="nls") 
      fitmat[,i]=predict(models[[i]],newdata=alo)
    
  }
  
  ne=c(1:9)[apply(fitmat,2,function(x){all(!is.na(x))})]
  plot(model.np,lwd=2)
  
  plot(ld,xlim=range(yframe,na.rm=T),col=4,pch=20,ylim=range(fitmat,na.rm=T))
  for(i in ne){
    lines(as.numeric(unlist(alo)),fitmat[,i],col=i)  
  }
  legend("bottomright", names(models)[ne],col=ne,lty=1)
  
  
  mapes=c()
  for(j in 1:length(ne)){
    i=ne[j]
    mapes[j]=mean(abs(ld$x-predict(models[[i]]))/ld$x)*100   
  }
  names(mapes)=as.character(ne)
  best=as.numeric(names(mapes[which(mapes==min(mapes))]))
  
  plot(ld,xlim=range(yframe,na.rm=T),col=4,pch=20,ylim=range(fitmat,na.rm=T),main=names(models)[best])
  lines(as.numeric(unlist(alo)),fitmat[,best],col=best)  
  
  bestfit=as.numeric(fitmat[,best])
  
  
  # jei nesuranda jokio gero modelio, reik kazka daryt:
  if(sd(bestfit)/mean(bestfit)>5 | any(bestfit<(-10)) | any(bestfit>180)){
    cat("NEI VIENAS IS S KREIVES MODELIU NETINKA!","\n")
  }
  
  
  
  dd=d
  d=d[-10,]
  
  p1y=c(seq(min(yframe,na.rm=T),min(d$y),by=100),min(d$y)+100)
  p2y=c(max(d$y)-100,seq(max(d$y),max(yframe,na.rm=T),by=100))
  
  p1fit=predict(models[[best]],newdata=data.frame(y=p1y))
  p2fit=predict(models[[best]],newdata=data.frame(y=p2y))
  
  r1=p1fit[2:length(p1fit)]/p1fit[1:(length(p1fit)-1)]
  r2=p2fit[2:length(p2fit)]/p2fit[1:(length(p2fit)-1)]
  
  p1fitf=1/rev(cumprod(rev(r1)))*foremaster[1]
  p2fitf=cumprod(r2)*foremaster[length(foremaster)]
  
  miss=c(p1fitf,foremaster,p2fitf)
  
  dfor=cbind(x=miss,y=c(p1y[-length(p1y)],seq(min(d$y),max(d$y),by=100),p2y[-1]))
  
  if(plottest)plot(dfor[,2],dfor[,1],type="l")
  points(ld,col=4,pch=20)
  
  yframe=monotonic.income(yframe)
  
  matout=data.frame(matrix(NA,nrow=10,ncol=2030-2005+1))
  
  for(y in 2005:2030){
    cps.eval <- data.frame(y=yframe[as.character(y)])
    ff=lininterpol(cps.eval,dfor)
    matout[y-2004]=ff   
  }
  
  tempmatout=matout
  # normavimas --------------------------------------------------------------
  
  ttn=as.numeric(tt)
  
  #   matout=monotonic(matout,ttn)
  matout=iternorm(matout,ttn,plottest=T)
  
  qwe=nutemptas.total.brutal(tempmatout,tt)
  #   matplot(t(data.frame(over100(qwe))),type="o",lty=1,pch=20,cex=0.8,main="NUTEMPTAS TOTAL")
  
  matout=data.frame(over100(matout))
  
  
  par(mfrow=c(1,1))
  
  if(plottest){
    matplot(matout,type="o",lty=1,pch=20,cex=1)
    matplot(t(matout),type="o",lty=1,pch=20,cex=0.8)
    lines(as.numeric(apply(matout,2,mean)),col="firebrick3",lwd=1.5,pch=8,type="o")
  }
  
  
  
  if(plottest3){
    matplot(t(data.frame(over100(qwe))),type="o",lty=1,pch=20,cex=0.8,main="NUTEMPTAS TOTAL",ylim=range(matout))
    lines(ttn,col="firebrick3",lwd=1.5,pch=8,type="o")
    matplot(t(matout),type="o",lty=1,pch=20,cex=0.8,ylim=range(matout),main="Final (normalus metodas)")
    lines(ttn,col="firebrick3",lwd=1.5,pch=8,type="o")
    
  }
  
  matout=as.data.frame(matout)
  names(matout)=as.character(2005:2030)
  return(matout)
}

emibase=function(){
  cat("Loading EMI base package ...","\n")
  load("K:\\GMID Research\\SM-route\\EMI base\\EMI base.Rdata",envir = .GlobalEnv)
  cat("Succesfully loaded","\n")
}



titlepage <- function(title,sub="",sizes=15,date=Sys.time(),cola="skyblue4",author=""){
  
  df <- data.frame()
  
  
  ggplot(df)+theme_bw()+
    geom_rect(aes(xmin = -5, xmax = 5 , ymin = -0, ymax = 5),fill=cola)+
    geom_text(aes( x = -4.75, y = 0.25*11/8),
              hjust=0, vjust=0,
              label = title,
              size = sizes, colour = "white",fontface=1,family="AvantGarde")+
    geom_text(aes( x = 4.75, y = -4),
              hjust=1, vjust=1.2,
              label = date,
              size = sizes/2, colour = cola,fontface=1,family="AvantGarde")+
    geom_text(aes( x = 4.75, y = -4),
              hjust=1, vjust=-0.2,
              label = author,
              size = sizes/2, colour = cola,fontface=1,family="AvantGarde")+
    geom_text(aes( x = -4.75, y = -0.25*11/8),
              hjust=0, vjust=1,
              label = sub,
              size = sizes/1.5, colour = cola,fontface=1,family="AvantGarde")+
    scale_y_continuous(limits=c(-5, 5),expand=c(0,0))+
    scale_x_continuous(limits=c(-5, 5),expand=c(0,0))+
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.grid = element_blank(), 
          axis.title=element_blank(),
          panel.border = element_blank(),
          plot.margin= unit(c(-0.25,-0.25,-5.25,-5.25),"mm"))
}


thief <- function(x,c=3,plot=F,ratio=F,critvalue=0.05,abscritvalue = 0){
  
  dx <- diff(x)
  if(ratio) dx <- x[2:length(x)]/x[1:(length(x)-1)] - 1
  
  ranges <- c(mean(dx,na.rm=T) - c * sd(dx,na.rm=T),mean(dx,na.rm=T) + c * sd(dx,na.rm=T))
  if(plot) plotg(dx,abline=ranges)
  
  out <- F
  
  if((any( dx < ranges[1] | dx > ranges[2], na.rm=T )) & 
       max(abs(dx),na.rm=T)/mean(x,na.rm=T) >= critvalue & 
       mean(x,na.rm=T) >= abscritvalue) out <- T
  return(out)
  
}


thiefnr <- function(x,c=3,plot=F,ratio=F){
  
  dx <- diff(x)
  if(ratio) dx <- x[2:length(x)]/x[1:(length(x)-1)] - 1
  
  ranges <- c(mean(dx,na.rm=T) - c * sd(dx,na.rm=T),mean(dx,na.rm=T) + c * sd(dx,na.rm=T))
  if(plot) plotg(dx,abline=ranges)
  
  out <- which( dx < ranges[1] | dx > ranges[2]) + 1
  return(out)
  
}



rankg2 <- function(xframe,main="",ran=NULL,outgr=F){
  
  nn <- names(xframe)[1]
  names(xframe)[1] <- "Object"
  mxframe <- melt(xframe,id.vars=names(xframe)[1])
  mxframe <- ddply(mxframe,~variable,function(x){
    cbind(x,Place = length(x$value) - rank(x$value) + 1)
  })
  
  if(is.null(ran)) ran <- range(xframe[,2:3])
  g <- ggplot(mxframe,aes(x=variable,y=Place,col=Object,group=Object)) + geom_line(size=rel(1)) + 
    theme_bw() + theme_pb() + scale_colour_emi(name=nn) + ggtitle(main)  +
    scale_size_continuous(limits=ran) + geom_point(aes(label=Object,size=value)) + 
    geom_text(aes(label=Object),colour="grey60",vjust=1.5) +
    scale_y_reverse(breaks=1:dim(xframe)[1]) + 
    theme(panel.grid.minor = element_blank())
  
  if(outgr) return(g) else  print(g)
  
}


rankg <- function(xframe,main="",ran=NULL,outgr=F){
  
  nn <- names(xframe)[1]
  names(xframe)[1] <- "Object"
  xframe <- cbind(xframe,rank1 = length(xframe[,2]) - rank(xframe[,2]) + 1,
                  rank2 = length(xframe[,3]) - rank(xframe[,3]) + 1,mean = apply(xframe[,2:3],1,mean))
  
  
  if(is.null(ran)) ran <- range(xframe[,2:3])
  g <- ggplot(xframe,aes(x=rank1,y=rank2,col=Object,group=Object))  + 
    geom_abline(intercept = 0, slope = 1) + 
    geom_point() + 
    theme_bw() + theme_pb() + scale_colour_emi(name=nn) + ggtitle(main)  +
    scale_size_continuous(name="Average Value",limits=ran) + geom_point(aes(label=Object,size=mean))+
    xlab(paste("Place:",names(xframe)[2])) + ylab(paste("Place:",names(xframe)[3])) +
    geom_text(aes(label=Object),colour="grey60",vjust=1.5) + 
    scale_x_reverse(breaks=1:dim(xframe)[1]) + scale_y_reverse(breaks=1:dim(xframe)[1]) + 
    theme(panel.grid.minor = element_blank())
  
  if(outgr) return(g) else  print(g)
  
}


mutation <- function(x,y,alpha=1){
  mult <- seq(0,1,length.out=length(x))^alpha
  z <- x * (1-mult) + y * mult
  return(z)
}


gr <- function(x){
  c(NA,x[-1]/x[-length(x)] - 1)
}


daug <- function(x,gro){
  
  outx <- x
  nna <- range(which(!is.na(x)))
  
  if(nna[2]!=length(x)){
    outx[(nna[2]+1):length(outx)] <- outx[nna[2]] * cumprod(gro[(nna[2]+1):length(gro)] + 1)
  }
  
  if(nna[1]!=1){
    outx[1:(nna[1]-1)] <- outx[nna[1]] * rev(cumprod(rev(1/(gro[1:(nna[1]-1) + 1] + 1))))
  }
  return(outx)
  
}




regionalpos <- function(k,pframe,plot=F,starter = 10000000,itermax = 10){
  
  abq <- as.numeric(k[,c("a","b.cntUSD","q")])
  xframe <- pframe[,c("Average.income","est.Pos.norm")]
  
  #funs
  weights <- function(x){dsm(x,abq)}   
  posfun <- function(x){
    y <- splinefun(xframe,method = "monoH.FC")(x)
    y[y<0] <- 0
    y[y>100] <- 100
    return(y)
  }
  fpcompose <- function(x){weights(x)*posfun(x)}
  ###
  
  if(plot){
    plot(xframe,type="p",pch=20,xlim=c(0,200000),main=paste(k$CountryCode,k$City,k$Year,sep=" - "),col="firebrick3")
    curve(posfun,0,200000,add=T)
    par(new=TRUE)
    curve(weights,0,200000,xaxt="n",yaxt="n",xlab="",ylab="",col=4)
    axis(4)
    
    #     curve(fpcompose,0,200000,main=paste(k$CountryCode,k$City,k$Year,sep=" - "))
    #     points(xframe,pch=20)
  }
  
  divider <- 1
  out <- NA
  iter <- 0
  while((class(out)=="try-error" | is.na(out)) & iter < itermax){
    out <- try(integrate(fpcompose,lower = 0,upper = starter/divider)$value)
    divider <- divider * 5
    iter <- iter + 1
  }
  
  if(iter >= itermax){print("ALARM ALARM : MAX ITERATIONS REACHED")
                      break}
  #   
  #   if(!low){ 
  #     out <- try(integrate(fpcompose,lower = 0,upper = 100000000)$value)
  #     if(class(out)=="try-error") out <- integrate(fpcompose,lower = 0,upper = 10000000)$value
  #   } else {
  #     out <- try(integrate(fpcompose,lower = 0,upper = 1000000)$value)
  #     if(class(out)=="try-error") out <- integrate(fpcompose,lower = 0,upper = 100000)$value
  #   }
  
  return(out)
  
}



regionalposrest <- function(kof,pframe,plot=F,starter = 10000000){
  
  abq <- kof
  #   print(abq)
  
  xframe <- pframe
  
  #   funs
  weights <- function(x){dsm(x,abq)}   
  posfun <- function(x){
    y <- splinefun(xframe,method = "monoH.FC")(x)
    y[y<0] <- 0
    y[y>100] <- 100
    return(y)
  }
  fpcompose <- function(x){weights(x)*posfun(x)}
  ##
  
  if(plot){
    plot(xframe,type="p",pch=20,xlim=c(0,200000),col="firebrick3")
    curve(posfun,0,200000,add=T)
    par(new=TRUE)
    curve(weights,0,200000,xaxt="n",yaxt="n",xlab="",ylab="",col=4)
    axis(4)
  }
  
  if(any(abq<0) | any(is.na(abq))) {out <- -100} else {
    
    divider <- 1
    out <- NA
    while(class(out)=="try-error" | is.na(out)){
      out <- try(integrate(fpcompose,lower = 0,upper = starter/divider)$value)
      divider <- divider * 5
      if(divider>10000) {out <- 0;break}
    }
  }
  
  return(out)
  
}



stretch <- function(x,a,b){
  
  c <- (x[1] + x[length(x)])/2
  ratio <- ((b - a)/2) / (x[length(x)] - c)
  c <- as.numeric(c)
  xnew <- (x - c) * as.numeric(ratio) + (a + b)/2
  
  return(xnew)
}


posfunction<- function(x,pframe){
  y <- splinefun(pframe,method = "monoH.FC")(x)
  y[y<0] <- 0
  y[y>100] <- 100
  return(y)
}

posfunctionreverse<- function(x,pframe){
  if(all(pframe[,2]==mean(pframe[,2]))) pframe[1,2] <- pframe[1,2] - 0.000001
  y <- splinefun(pframe[,c(2,1)],method = "monoH.FC")(x)
  return(y)
}

pframetransform <- function(pframe){
  
  postotal <- mean(pframe$est.Pos.norm)
  incometotal <- mean(pframe$Average.income)
  
  incometotalfake <- posfunctionreverse(postotal,pframe)
  ratio <- incometotal/incometotalfake
  out <- pframe
  out[,1] <- out[,1] * ratio
  
  out <- list(out,ratio)
  
  return(out)
  
}


pframetransformup <- function(pframe){
  
  postotal <- mean(pframe$est.Pos.norm)
  incometotal <- mean(pframe$Average.income)
  
  postotalfake <- posfunction(incometotal,pframe)
  ratio <- postotal/postotalfake
  out <- pframe
  out[,2] <- out[,2] * ratio
  return(out)
  
}


b.converter <- function(a,q,inctot){
  inctot * a * gamma(q)/(gamma(1/a) * gamma(q - 1/a))
}

b.converter.reverse <- function(a,b,q){b * beta(1+1/a,q-1/a)/beta(1,q)}


test.all.conditions <- function(bigdata){
  
  ###########################################################################
  # testing conditions ------------------------------------------------------
  ###########################################################################
  
  if("ProductName" %in% names(bigdata)){
    names(bigdata)[names(bigdata)=="ProductName"] <- "Possession"
    namechanged <- T
  }
  
  CC <- unique(bigdata$CountryCode)
  
  for(cc in CC){
    #     print(cc)
    
    xframe <- bigdata[bigdata$CountryCode == cc,]
    
    pc <- xframe[xframe$Possession == "Possession of Personal Computer",years]
    int <- xframe[xframe$Possession == "Possession of Internet Enabled Computer",years]
    br <- xframe[xframe$Possession == "Possession of Broadband Internet Enabled Computer",years]
    
    ctv <- xframe[xframe$Possession == "Possession of Colour TV Set",years]
    sat <- xframe[xframe$Possession == "Possession of Satellite TV System",years]
    cbl <- xframe[xframe$Possession == "Possession of Cable TV",years]
    
    
    emptyPC <- all(is.na(pc))
    emptyINT <- all(is.na(int))
    emptyBR <- all(is.na(br))
    emptyCTV <- all(is.na(ctv))
    emptySAT <- all(is.na(sat))
    emptyCBL <- all(is.na(cbl))
    
    
    # PC INT BR ---------------------------------------------------------------
    
    #   matplotg(cbind(pc,int,br),text=T,stx=2005)
    nr <- max(c(nrow(ctv),nrow(sat),nrow(cbl),nrow(pc),nrow(int),nrow(br)))
    nc <- max(c(ncol(ctv),ncol(sat),ncol(cbl),ncol(pc),ncol(int),ncol(br)))
    
    if(emptyBR) br <- matrix(0,nrow=nr,ncol=nc)
    if(emptyINT) int <- br
    if(emptyPC) pc <- int
    
    pcint <- round(pc - int,2)    
    
    if(any(as.numeric(unlist(pcint)) < 0)){
      print(paste(cc,"int > pc",sep=" - "))
      #       matplotg(cbind(t(pc),t(int)),text=T,stx=2005,main=paste(cc,"INT > PC before"))
      int[pcint < 0] <- pc[pcint < 0]
      #       matplotg(cbind(t(pc),t(int)),text=T,stx=2005,main=paste(cc,"INT > PC after"))
      
      if(!emptyINT) bigdata[bigdata$CountryCode == cc & bigdata$Possession == "Possession of Internet Enabled Computer",years] <- int
    }
    
    intbr <- round(int - br,2)
    
    if(any(as.numeric(unlist(intbr)) < 0)){
      print(paste(cc,"br > int",sep=" - "))
      #             matplotg(cbind(t(int),t(br)),text=T,stx=2005,main=paste(cc,"BR > INT before"))
      br[intbr < 0] <- int[intbr < 0]
      #             matplotg(cbind(t(int),t(br)),text=T,stx=2005,main=paste(cc,"BR > INT after"))
      
      if(!emptyBR) bigdata[bigdata$CountryCode == cc & bigdata$Possession == "Possession of Broadband Internet Enabled Computer",years] <- br
    }
    
    # CTV SAT CBL -------------------------------------------------------------
    
    if(emptyCTV) ctv <- matrix(100,nrow=nr,ncol=nc)
    if(emptySAT) sat <- matrix(0,nrow=nr,ncol=nc)
    if(emptyCBL) cbl <- matrix(0,nrow=nr,ncol=nc)
    
    satcbl <- sat + cbl
    ctvsatcbl <- round(ctv - satcbl,2)
    satratio <- sat/satcbl
    
    if(any(as.numeric(unlist(ctvsatcbl)) < 0)){
      print(paste(cc,"sat + cbl > ctv",sep=" - "))
      #             matplotg(cbind(t(ctv),t(sat),t(cbl),satcbl=t(sat+cbl)),text=T,stx=2005,main=paste(cc,"SAT + CBL > CTV before"))
      sat[ctvsatcbl < 0] <- ctv[ctvsatcbl < 0] * satratio[ctvsatcbl < 0]
      cbl[ctvsatcbl < 0] <- ctv[ctvsatcbl < 0] * (1 - satratio[ctvsatcbl < 0])
      #             matplotg(cbind(t(ctv),t(sat),t(cbl),satcbl=t(sat+cbl)),text=T,stx=2005,main=paste(cc,"SAT + CBL > CTV after"))
      
      if(!emptySAT) bigdata[bigdata$CountryCode == cc & bigdata$Possession == "Possession of Satellite TV System",years] <- sat
      if(!emptyCBL) bigdata[bigdata$CountryCode == cc & bigdata$Possession == "Possession of Cable TV",years] <- cbl
    }
    
  }
  
  
  ###########################################################################
  # testing constraints < 100; > 0 ------------------------------------------
  ###########################################################################
  
  indexbad0 <- which(bigdata[,years] < 0,arr.ind = T)[,1]
  
  for(i in indexbad0){
    
    x <- bigdata[i,years]
    #     plotg(x,main=paste(bigdata[i,"City"],bigdata[i,"Possession"],"before",sep=" - "))
    x[x<0] <- 0
    #     plotg(x,main=paste(bigdata[i,"City"],bigdata[i,"Possession"],"after",sep=" - "))
    bigdata[i,years] <- x
    
  }
  
  indexbad100 <- which(bigdata[,years] > 100,arr.ind = T)[,1]
  
  for(i in indexbad100){
    
    x <- bigdata[i,years]
    #     plotg(x,main=paste(bigdata[i,"City"],bigdata[i,"Possession"],"before",sep=" - "))
    x[x>100] <- 100
    #     plotg(x,main=paste(bigdata[i,"City"],bigdata[i,"Possession"],"after",sep=" - "))
    bigdata[i,years] <- x
    
  }
  
  # if there WAS 100, there WILL BE 100
  indexwere100 <- which(bigdata[,years[-1:-5]] == 100,arr.ind = T)[,1]
  indexwere100 <- unique(indexwere100)
  
  for(i in indexwere100){
    
    print(i)
    x <- bigdata[i,years]
    nr <- min(which(x==100))
    #     plotg(x,main=paste(bigdata[i,"City"],bigdata[i,"Possession"],"before",sep=" - "))
    x[nr:length(x)] <- 100
    #     plotg(x,main=paste(bigdata[i,"City"],bigdata[i,"Possession"],"after",sep=" - "))
    bigdata[i,years] <- x
    
  }
  
  
  if(namechanged) {names(bigdata)[names(xframe)=="Possession"] <- "ProductName"}
  
  return(bigdata)
  
}



hatfun <- function(x,right = 1,stp=1,crit = 0.01,plot=F){
  
  out <- c()
  for(i in 2:(length(x)-right)){
    if((x[i] > x[i-1] & x[i] > (1 + crit) * x[i+right]) |
         (x[i] < x[i-1] & x[i] < (1 - crit) * x[i+right])) out[i] <- T else out[i] <- F
  }
  
  hatnr <- which(out)
  hatnr[hatnr<=stp] <- NA
  hatnr <- hatnr[!is.na(hatnr)]
  
  if(plot) {plot(x); points(hatnr,x[hatnr],col=2,pch=8)}
  
  return(hatnr)
  
}



distratio <- function(theta1,theta2,p=0.35,plot=F){
  
  divider <- (psm(qsm(0.5+p,theta1),theta2) - psm(qsm(0.5-p,theta1),theta2))
  
  retribution <- (integrate(f,qsm(0.5-p,theta1),qsm(0.5+p,theta1),theta=theta2)$value/divider)/
    (integrate(f,qsm(0.5-p,theta1),qsm(0.5+p,theta1),theta=theta1)$value/(p*2))
  
  if(plot){
    x <- seq(0,theta1[2]*2,length.out=1000)
    plot(x,dsm(x,theta1),type="l",main=paste0("p = ",round(p,2),", divider = ",round(divider,2),", ratio = ",round(retribution,2)))
    lines(x,dsm(x,theta2),type="l",col=4)
    abline(v=c(qsm(0.5+p,theta1),qsm(0.5-p,theta1)),col=3)
    abline(v=c(qsm(0.5,theta1)),col="darkgreen",lty=2)
    
  }
  
  return(retribution)
  
}

setp <- function(ptt,m=0.42,alpha=1){
  (m) * (1.02 - ptt/100)^(alpha)
}



hatfun.solve <- function(x,itermax=50,...){
  
  solve <- T
  iter <- 0
  while(solve & iter <= itermax){
    iter <- iter + 1
    hnr <- hatfun(x,...)
    if(length(hnr)>0){
      x[hnr] <- NA
      x <- MASplineVector(x)
    } else {solve <- F}
  }
  
  return(x)
}



thief.solve <- function(x,itermax=50,...){
  
  solve <- T
  iter <- 0
  while(solve & iter <= itermax){
    iter <- iter + 1
    hnr <- thiefnr(x,c...)
    if(length(hnr)>0){
      x[hnr] <- NA
      x <- MASplineVector(x)
    } else {solve <- F}
  }
  
  return(x)
}




displaysm <- function(theta,main=""){
  arg <- seq(0,theta[2]*4,by=10)
  plot(arg,dsm(arg,theta=theta),type="l",main=main)
}


durmatplot.simple <- function(xframe,totframe,legend=T,...){
  
  #   if(is.null(ylim)) ylim=range(xframe)
  xai <- as.numeric(rownames(xframe))
  matplot(xai,xframe,lty=1,type="l",col=emi_pal()(dim(xframe)[2]),xlab="Year",xlim=c(2003,2032),...)
  lines(totframe[,c(1,3)],lwd=3)
  grid(col="grey70")
  if(legend)legend("bottomright",col=emi_pal()(dim(xframe)[2]),lty=1,legend=names(xframe),cex=0.5)
  text(rep(2004,ncol(xframe)),xframe[1,],labels = names(xframe),col=emi_pal()(dim(xframe)[2]),cex=0.5)
  text(rep(2031,ncol(xframe)),xframe[nrow(xframe),],labels = names(xframe),col=emi_pal()(dim(xframe)[2]),cex=0.5)
  
  
}

durmatplot.simple.2 <- function(xframe,xframenew,totframe,main="",ylim=NULL){
  
  if(is.null(ylim)) ylim=range(xframe)
  xai <- as.numeric(rownames(xframe))
  matplot(xai,xframe,lty=1,type="l",main=main,col=emi_pal()(dim(xframe)[2]),xlab="Year",ylab="Possession %",ylim=ylim)
  grid(col="grey85")
  matplot(xai,xframenew,lty=2,type="l",main=main,col=emi_pal()(dim(xframe)[2]),xlab="Year",ylab="Possession %",ylim=ylim,add=T)
  
  lines(totframe[,c("Index","value")],lwd=3)
  legend("bottomright",col=emi_pal()(dim(xframe)[2]),lty=1,legend=names(xframe),cex=0.5)
  
}



# allpos graph ------------------------------------------------------------

allpos.graph <- function(allpos,dt,incallnormabq,date=F){
  
  unikas <- unique(allpos[,c("CountryCode","MeasureType","ProductName")])
  unikas <- arrange(unikas,MeasureType,CountryCode)
  
  filename <- "plots/allpos.pdf"
  if(date) filename <- paste0("plots/allpos ",Sys.Date(),".pdf")
  
  pdf(file=filename,width=13,height=10)
  titlepage("All Durable Goods",sub = "Forecasts",author="Povilas. B")
  
  for(i in 1:nrow(unikas)){
    #   for(i in 1:86){
    
    cc <- unikas[i,1]
    mt <- unikas[i,2]
    pn <- unikas[i,3]
    
    print(paste(cc,mt,pn,sep=" - "))
    
    xframe <- allpos[allpos$CountryCode==cc & allpos$MeasureType==mt & allpos$ProductName==pn,]
    nxframe <- xframe$Subsector
    xframe <- data.frame(t(xframe[,years]))
    names(xframe) <- nxframe
    totframe <- as.numeric(dt[dt$CountryCode==cc & dt$ProductName==pn,paste0("Y",years)])
    totframe <- data.frame(cbind(Index=years,variable="Total",value = totframe))
    totframe$value <- as.numeric(as.character(totframe$value))
    totframe$Index <- as.numeric(as.character(totframe$Index))
    
    layout(matrix(c(1,1,1,1,2,3,4,5), 2, 4, byrow = TRUE),heights = c(2,1))
    
    durmatplot.simple(xframe,totframe,main=paste(cc,mt,pn,sep=" - "))
    
    #   g1 <- matplotg(xframe,main=paste(cc,mt,pn,sep=" - "),stx=2005,outgr=T,
    #                  nrowslegend = 2,text=T,point=F,numbertext = F)
    #   g1 <- g1 + geom_line(data=totframe,col="black",size=rel(1.2))
    #   print(g1)
    
    #  ------------------------------------------------------------------------
    
    framepnabq <- incallnormabq[incallnormabq$CountryCode==cc & incallnormabq$Possession==pn & 
                                  incallnormabq$MeasureType==mt,]
    framepnabq$Possession <- NULL
    mframepnabq <- melt(framepnabq,id.vars=names(framepnabq)[1:4])
    names(mframepnabq)[5] <- "var"
    
    mframepnabq <- dcast(mframepnabq,CountryCode+MeasureType + variable + var ~ Subsector)
    
    kofcn <- kofrecall[kofrecall$CountryCode==cc & kofrecall$Possession==pn,]
    kofcn$CountryName <- b.converter.reverse(kofcn$a,kofcn$b,kofcn$q)
    names(kofcn)[length(kofcn)] <- "value"
    names(kofcn)[3] <- "variable"
    mkofcn <- melt(kofcn,id.vars=names(kofcn)[1:3])
    names(mkofcn)[4] <- "var"
    mkofcn$variable <- as.numeric(as.character(mkofcn$variable))
    
    durmatplot.simple(data.frame(mframepnabq[mframepnabq$var=="value",-1:-4],row.names = years,check.names=F),
                      mkofcn[mkofcn$var=="value",-1:-2],main="Average Income",legend=F,ylab="")
    
    durmatplot.simple(data.frame(mframepnabq[mframepnabq$var=="a",-1:-4],row.names = years,check.names=F),
                      mkofcn[mkofcn$var=="a",-1:-2],main="a",legend=F,ylab="")
    
    durmatplot.simple(data.frame(mframepnabq[mframepnabq$var=="b",-1:-4],row.names = years,check.names=F),
                      mkofcn[mkofcn$var=="b",-1:-2],main="b",legend=F,ylab="")
    
    durmatplot.simple(data.frame(mframepnabq[mframepnabq$var=="q",-1:-4],row.names = years,check.names=F),
                      mkofcn[mkofcn$var=="q",-1:-2],main="q",legend=F,ylab="")
    
    #   g2 <- ggplot(mframepnabq,aes(x=variable,y=value,colour=Subsector)) + theme_bw() + theme_pb()+
    #     facet_wrap(~var,nrow=1,scales = "free_y") + geom_line() + scale_colour_emi() + 
    #     theme(strip.background = element_rect(colour = "grey20", fill = "skyblue4"),
    #           strip.text = element_text(colour = "white")) + 
    #     geom_line(data=mkofcn,col="black",size=rel(1.2)) + theme(legend.position="none")
    #   
    #   
    #   g <- arrangeGrob(g1,g2,heights=c(3,1),nrow=2)
    #   print(g)
    
  }
  
  dev.off()
  
}




# over100 v2 --------------------------------------------------------------

over100.v2 <- function(xframe,border=1,divover=2,plot=F){
  
  goodmax <- max(xframe[xframe<100])
  if(length(which(xframe > 100,arr.ind = T)[,2])==0) minbadnr <- Inf else 
    minbadnr <- min(which(xframe > 100,arr.ind = T)[,2])
  oxframe <- xframe
  
  xframe <- t(apply(xframe,1,function(x){
    if(any(x > 100)){
      #       bcnr <- max(c(1,min(which(x > 100)) - border))
      bcnr <- max(c(1,minbadnr - border))
      
      if(as.character(x[bcnr]) == as.character(goodmax)){
        
        x[x > 100] <- as.numeric(x[bcnr])
      } else {
        
        x[x > 100] <- (x[x > 100] - 100)/divover + x[x > 100]
        if (mean(as.numeric(x[bcnr:length(x)]))!=100){
          x[bcnr:length(x)] <- stretch(x[bcnr:length(x)],x[bcnr],goodmax)
        }
      }
    }
    return(x)
  }))
  
  xframe[xframe>100] <- 100
  
  if(plot){
    #       matplotf(t(oxframe),abline=c(0,100),text=T,main="Before")
    matplotf(t(xframe),abline=c(0,100),text=T,main="After")
    
  }
  
  return(xframe)
  
}



ce.constantination <- function(bigdataframe,cpi){
  
  for(cc in unique(bigdataframe$CountryCode)){
    
    xframe <- bigdataframe[bigdataframe$CountryCode==cc,]
    
    xframe[,years] <- t(apply(xframe[,years],1,function(x){  # sitas yra tik del to, kad nera padarytas tenure!!!
      MASplineVector(x,k = 1.1)
    }))
    
    cpicc <- cpi[cpi$CountryCode==cc,years]
    cpicc <- cpicc/cpicc[,"2013"]*100
    
    xframe[,years] <- t(apply(xframe[,years],1,function(x){
      x/as.numeric(cpicc)*100
    }))
    
    bigdataframe[bigdataframe$CountryCode==cc,] <- xframe
    
  }
  
  return(bigdataframe)
}


converter <- function(x,cc,year){
  
  xframe <- framelist[[cc]]
  imat <- xframe[xframe$ProductName=="inc",-1:-3]
  emat <- xframe[xframe$ProductName=="exp",-1:-3]
  
  imean <- apply(imat,2,mean)
  emean <- apply(emat,2,mean)
  
  ratio <- as.numeric(imean[year-2004]/emean[year-2004])
  
  x * ratio
}



# testing.all.conditions.cities -------------------------------------------

test.all.conditions.cities <- function(bigdata,plot=F){
  
  ###########################################################################
  # testing conditions ------------------------------------------------------
  ###########################################################################
  
  if(plot){
    pdf(file=paste("plots/test all conditions cities",Sys.Date(),".pdf"),width=13,height=8)
    titlepage("Regions Conditions",sub = "",author="Povilas. B")
  }
  
  if("ProductName" %in% names(bigdata)){
    names(bigdata)[names(bigdata)=="ProductName"] <- "Possession"
    namechanged <- T
  }
  
  
  cclist <- unique(bigdata[,c("CityCode")])
  #   cclist <- cclist[cclist!="Total"]
  
  for(cc in cclist){
    #     print(cc)
    
    xframe <- bigdata[bigdata$CityCode == cc,]
    
    # pc <- opc <- xframe[xframe$Possession == "Possession of Personal Computer",years]
    int <- oint <- xframe[xframe$Possession == "Possession of Internet Enabled Computer",years]
    br <- obr <- xframe[xframe$Possession == "Possession of Broadband Internet Enabled Computer",years]
    
    ctv <- octv <- xframe[xframe$Possession == "Possession of Colour TV Set",years]
    sat <- osat <- xframe[xframe$Possession == "Possession of Satellite TV System",years]
    cbl <- ocbl <- xframe[xframe$Possession == "Possession of Cable TV",years]
    
    
    # emptyPC <- all(is.na(pc))
    emptyINT <- all(is.na(int))
    emptyBR <- all(is.na(br))
    emptyCTV <- all(is.na(ctv))
    emptySAT <- all(is.na(sat))
    emptyCBL <- all(is.na(cbl))
    
    
    # PC INT BR ---------------------------------------------------------------
    
    #   matplotg(cbind(pc,int,br),text=T,stx=2005)
    # nr <- max(c(nrow(ctv),nrow(sat),nrow(cbl),nrow(pc),nrow(int),nrow(br)))
    # nc <- max(c(ncol(ctv),ncol(sat),ncol(cbl),ncol(pc),ncol(int),ncol(br)))
    
    nr <- max(c(nrow(ctv),nrow(sat),nrow(cbl),nrow(int),nrow(br)))
    nc <- max(c(ncol(ctv),ncol(sat),ncol(cbl),ncol(int),ncol(br)))
    
    if(emptyBR) br <- matrix(0,nrow=nr,ncol=nc)
    if(emptyINT) int <- br
    # if(emptyPC) pc <- int
    
    PG <- F
    
    # pcint <- round(pc - int,2)    
    
    # if(any(as.numeric(unlist(pcint)) < 0)){
    #   print(paste(cc,"int > pc",sep=" - "))
    #   int[pcint < 0] <- pc[pcint < 0]      
    #   if(!emptyINT) bigdata[bigdata$CityCode == cc & bigdata$Possession == "Possession of Internet Enabled Computer",years] <- int
    #   PG <- T
    # }
    
    intbr <- round(int - br,2)
    
    if(any(as.numeric(unlist(intbr)) < 0)){
      print(paste(cc,"br > int",sep=" - "))
      br[intbr < 0] <- int[intbr < 0]
      if(!emptyBR) bigdata[bigdata$CityCode == cc & bigdata$Possession == "Possession of Broadband Internet Enabled Computer",years] <- br
      PG <- T
    }
    
    if(plot & PG){
      # matplotf(data.frame(cbind(pc=as.numeric(t(opc)),int=as.numeric(t(oint)),br=as.numeric(t(obr)))),text=T,stx=2005,main=paste(xframe$CountryCode[1],xframe$City[1],"before",sep=" - "),abline=c(0,100))
      # matplotf(data.frame(cbind(pc=as.numeric(t(pc)),int=as.numeric(t(int)),br=as.numeric(t(br)))),text=T,stx=2005,main=paste(xframe$CountryCode[1],xframe$City[1],"after",sep=" - "),abline=c(0,100))
      matplotf(data.frame(cbind(int=as.numeric(t(oint)),br=as.numeric(t(obr)))),text=T,stx=2005,main=paste(xframe$CountryCode[1],xframe$City[1],"before",sep=" - "),abline=c(0,100))
      matplotf(data.frame(cbind(int=as.numeric(t(int)),br=as.numeric(t(br)))),text=T,stx=2005,main=paste(xframe$CountryCode[1],xframe$City[1],"after",sep=" - "),abline=c(0,100))
    }
    
    
    # CTV SAT CBL -------------------------------------------------------------
    
    if(emptyCTV) ctv <- matrix(100,nrow=nr,ncol=nc)
    if(emptySAT) sat <- matrix(0,nrow=nr,ncol=nc)
    if(emptyCBL) cbl <- matrix(0,nrow=nr,ncol=nc)
    
    satcbl <- sat + cbl
    ctvsatcbl <- round(ctv - satcbl,2)
    satratio <- sat/satcbl
    
    if(any(as.numeric(unlist(ctvsatcbl)) < 0)){
      print(paste(cc,"sat + cbl > ctv",sep=" - "))
      #             matplotf(cbind(t(ctv),t(sat),t(cbl),satcbl=t(sat+cbl)),text=T,stx=2005,main=paste(cc,"SAT + CBL > CTV before"))
      sat[ctvsatcbl < 0] <- ctv[ctvsatcbl < 0] * satratio[ctvsatcbl < 0]
      cbl[ctvsatcbl < 0] <- ctv[ctvsatcbl < 0] * (1 - satratio[ctvsatcbl < 0])
      #             matplotf(cbind(t(ctv),t(sat),t(cbl),satcbl=t(sat+cbl)),text=T,stx=2005,main=paste(cc,"SAT + CBL > CTV after"))
      
      if(!emptySAT) bigdata[bigdata$CityCode == cc & bigdata$Possession == "Possession of Satellite TV System",years] <- sat
      if(!emptyCBL) bigdata[bigdata$CityCode == cc & bigdata$Possession == "Possession of Cable TV",years] <- cbl
    }
    
  }
  
  if(plot){dev.off()}
  
  ###########################################################################
  # testing constraints < 100; > 0 ------------------------------------------
  ###########################################################################
  
  indexbad0 <- unique(which(bigdata[,years] < 0,arr.ind = T)[,1])
  
  for(i in indexbad0){
    
    x <- bigdata[i,years]
    #     plotg(x,main=paste(bigdata[i,"City"],bigdata[i,"Possession"],"before",sep=" - "))
    x[x<0] <- 0
    #     plotg(x,main=paste(bigdata[i,"City"],bigdata[i,"Possession"],"after",sep=" - "))
    bigdata[i,years] <- x
    
  }
  
  indexbad100 <- unique(which(bigdata[,years] > 100,arr.ind = T)[,1])
  
  for(i in indexbad100){
    
    x <- bigdata[i,years]
    #         plotg(x,main=paste(bigdata[i,"City"],bigdata[i,"Possession"],"before",sep=" - "))
    x[x>100] <- 100
    #         plotg(x,main=paste(bigdata[i,"City"],bigdata[i,"Possession"],"after",sep=" - "))
    bigdata[i,years] <- x
    
  }
  
  # if there WAS 100, there WILL BE 100
  indexwere100 <- which(bigdata[,years[-1:-5]] == 100,arr.ind = T)[,1]
  indexwere100 <- unique(indexwere100)
  
  for(i in indexwere100){
    
    print(i)
    x <- bigdata[i,years]
    nr <- min(which(x==100))
    #         plotg(x,main=paste(bigdata[i,"City"],bigdata[i,"Possession"],"before",sep=" - "))
    x[nr:length(x)] <- 100
    #         plotg(x,main=paste(bigdata[i,"City"],bigdata[i,"Possession"],"after",sep=" - "))
    bigdata[i,years] <- x
    
  }
  
  
  #   if(namechanged) {names(bigdata)[names(xframe)=="Possession"] <- "ProductName"}
  
  return(bigdata)
  
}




avgdec.from.qbq <- function(theta){
  
  pk <- function(xx,theta){integrate(qsm,xx[1],xx[2],subdivisions=1000,theta=theta)$value}  
  
  q<-c(seq(0,1,by=0.1))
  fit <- c()
  for(i in 1:(length(q)-1)){
    fit <- c(fit,pk(c(q[i],q[i+1]),theta)/(q[i+1] - q[i]))
  }
  
  return(fit)
}




cities.possession.graphs <- function(bigdata,date=T){
  
  ucp <- unique(bigdata[,c("CityCode")])
  ucp <- ucp[ucp %in% tier$CityCode & !(ucp %in% bigdata$CountryCode)] # cia kad breztu tik miestus
  
  # ucp[780:810]
  # 16*ccc,h=13*ccc
  pdftitle <- paste0("plots/cities possession ",gsub(":","-",Sys.time()),".pdf")
  
  pdf(pdftitle,width=16 * 0.75,height=13 * 0.75)
  titlepage("Forecasted Regions and Cities", "Possessions",date=Sys.time(),size=15,author="Povilas Bockus")
  
  for(index in 1:length(ucp)){
    #   for(index in 1:10){
    
    print(ucp[index])
    
    cc <- as.character(ucp[index])
    #   pn <- as.character(ucp[index,"Possession"])
    
    xframe <- bigdata[bigdata$CityCode==cc,]
    
    cn <- xframe$CountryCode[1]  
    rc <- tier[tier$CityCode==cc,"RegionCode"][1]
    
    xframe <- cbind(xframe,RegionCode = NA)
    xframe$RegionCode <- rc
    maintitle <- paste(cn,paste0(tier[tier$CityCode==cc,"Region"]," (",rc,")"),paste0(xframe$City[1]," (",cc,")"),sep=" - ")
    
    totxframe <- bigdata[bigdata$CityCode==cn & bigdata$Possession %in% xframe$Possession,]
    regxframe <- bigdata[bigdata$CityCode==rc & bigdata$Possession %in% xframe$Possession,]
    
    #   group1 <- c(2495,2504,2505)
    #   group2 <- c(12981,2498,2515,59926,2502,2497,2514)
    
    for(groupnr in 1:6){
      
      if(groupnr==1) group <- c(2495,2504,2505)
      if(groupnr==2) group <- c(12981,2498,2515,59926,2502,2497,2514)
      if(groupnr==3) group <- c(2496,2499,12980,12987)
      if(groupnr==4) group <- c(2511,  68800,	2506,	59927,	87218	,152648)
      if(groupnr==5) group <- c(2512,  2513,	2516,	2500,	12988,	2494)
      if(groupnr==6) group <- c(12982,  2501,	2503,	2508)
      
      print(group)
      
      # pagrindinis grafikas
      
      xframe.g <- xframe[xframe$ProductID %in% group,]
      totxframe.g <- totxframe[totxframe$ProductID %in% group,]
      regxframe.g <- regxframe[regxframe$ProductID %in% group,]
      
      xframe.g <- data.frame(t(xframe.g[,years]))
      totxframe.g <- data.frame(t(totxframe.g[,years]))
      regxframe.g <- data.frame(t(regxframe.g[,years]))
      
      emptyregtest <- dim(regxframe.g)[2]==0
      
      names(xframe.g) <- names(totxframe.g) <- xframe[xframe$ProductID %in% group,]$Possession
      if(!emptyregtest) names(regxframe.g) <- xframe[xframe$ProductID %in% group,]$Possession
      
      layout(matrix(c(1,1,1,2,3,4), 2, 3, byrow = TRUE),heights = c(2.5,1))
      
      
      matplotf(totxframe.g,legend = F,abline=c(0,100,120),stx=2005,lty=1,lwd=0.75,turnoffpar = T,yaxt='n',main=maintitle)
      if(!emptyregtest) matplotf(regxframe.g,legend = F,abline=c(0,100),stx=2005,lty=5,add=T,lwd=1,turnoffpar = T)
      matplotf(xframe.g,legend = F,abline=c(0,100),stx=2005,add=T,lwd=3.25,turnoffpar = T)
      axis(2, at=seq(0,100,by=20), labels=as.character(seq(0,100,by=20)))
      
      rect(-1, 150,3000,105,col = "white")
      
      legend("topleft",col=emi_pal()(ncol(xframe.g)),lty=1,box.col="white",lwd=3.25,
             legend = names(xframe.g),cex=1,horiz=F,y.intersp = 0.8)
      
      if(!emptyregtest){
        legend("topright",lty=c(1,5,1),lwd=c(3.25,1,0.75),box.col="white",
               legend = c("City","Region","Country"),cex=1,y.intersp = 0.8)
      } else {
        legend("topright",lty=c(1,6),lwd=c(3.25,0.75),box.col="white",
               legend = c("City","Country"),cex=1,y.intersp = 0.8)  
      }
      
      cola <- emi_pal()(ncol(xframe.g))
      names(cola) <- xframe[xframe$ProductID %in% group,]$Possession
      
      # w(w=17,h=12)
      
      
      # originaliu duomenu rankingas
      
      years2 <- as.character(2005:2020)
      ccrank <- c()
      codeslength <- c()
      
      for(pn in unique(xframe[xframe$ProductID %in% group,]$Possession)){
        
        his <- hp[hp$CountryCode==cn & hp$ProductName==pn,]
        xis <- bigdata[bigdata$CountryCode==cn & bigdata$Possession==pn & bigdata$CityCode %in% tier$CityCode,]
        codes <- intersect(xis$CityCode,his$RegionCode)  
        codeslength[pn] <- length(codes)
        
        if(length(codes)==0) next
        if(!cc %in% codes) next
        
        his <- arrange(his[his$RegionCode %in% codes,c("RegionCode",years2)],RegionCode)
        xis <- arrange(xis[xis$CityCode %in% codes,c("CityCode",years2)],CityCode)
        
        xis[which(is.na(his),arr.ind = T)] <- NA
        
        rankmat <- data.frame(cbind(CityCode = xis[,1],Original = apply(his[,-1],1,function(x){mean(x,na.rm=T)}),
                                    Modelled = apply(xis[,-1],1,function(x){mean(x,na.rm=T)})),stringsAsFactors=F)
        rankmat[,2:3] <- apply(rankmat[,2:3],2,as.numeric)
        rankmat[,2:3] <- apply(rankmat[,2:3],2,rank)
        rankmat[,2:3] <- apply(rankmat[,2:3],2,function(x){nrow(rankmat) + 1 - x})
        
        ccrank <- rbind(ccrank,cbind(Possession = pn,rankmat[rankmat[,1]==cc,]))
        
      }
      
      if(is.null(ccrank)){
        plot(c(1,1),type="n",main="Orignal vs Modelled ranking",ylab="",xlab="",yaxt='n',xaxt='n')
        text(x=1.5,y=1,labels = "No Original Data for This City")
      } else {
        plot(ccrank[,c("Original","Modelled")],pch=20,col=cola[ccrank$Possession],
             xlim=rev(c(1,max(codeslength))),ylim=rev(c(1,max(codeslength))),main="Orignal vs Modelled ranking")
        abline(a=0,b=1,col="grey50")  
        points(ccrank[,c("Original","Modelled")],pch=20,col=cola[ccrank$Possession],cex=1.8,
               xlim=rev(c(1,max(codeslength))),ylim=rev(c(1,max(codeslength))))
      }
      
      
      # prognoziu rankingas
      
      ccrank2 <- c()
      codeslength2 <- c()
      
      for(pn in unique(xframe[xframe$ProductID %in% group,]$Possession)){
        
        xis <- bigdata[bigdata$CountryCode==cn & bigdata$Possession==pn & bigdata$CityCode %in% tier$CityCode,]
        codeslength2[pn] <- dim(xis)[1]
        
        rankmat <- data.frame(cbind("CityCode" = xis[,"CityCode"],"2014" = xis[,"2014"],"2030" = xis[,"2030"]),stringsAsFactors=F,check.names=F)
        rankmat[,2:3] <- apply(rankmat[,2:3],2,as.numeric)
        rankmat[,2:3] <- apply(rankmat[,2:3],2,rank)
        rankmat[,2:3] <- apply(rankmat[,2:3],2,function(x){nrow(rankmat) + 1 - x})
        
        ccrank2 <- rbind(ccrank2,cbind(Possession = pn,rankmat[rankmat[,1]==cc,]))
        
      }
      
      
      plot(ccrank2[,c("2014","2030")],pch=20,col=cola[ccrank2$Possession],
           xlim=rev(c(1,max(codeslength2))),ylim=rev(c(1,max(codeslength2))),main="2014 vs 2030 forecasts ranking")
      abline(a=0,b=1,col="grey50")  
      points(ccrank2[,c("2014","2030")],pch=20,col=cola[ccrank2$Possession],cex=1.8,
             xlim=rev(c(1,max(codeslength2))),ylim=rev(c(1,max(codeslength2))))  
      
      # miestu ir saliu income
      
      avicn <- avicity[avicity$CountryCode==cn,]
      avicn <- avicn[,c("CityCode","Year","AI")]
      avicn <- dcast(avicn,CityCode ~ Year,value.var="AI")
      avicn[,-1] <- t(apply(avicn[,-1],1,function(x){
        if(!all(is.na(x))) x <- MASplineVector(x)
        return(x)
      }))
      avicno <- avicn
      avicn <- data.frame(t(avicn[-1]))
      names(avicn) <- avicno[,1]
      
      inccolour <- "#00AED9"
      
      matplotf(avicn,legend = F,stx=2005,lwd=0.8,turnoffpar = T,main="Average Income",col="grey65")
      if(cn %in% names(avicn)) lines(2005:2030,as.numeric(avicn[,cn]),lwd=1.5,col=4,lty=1)
      if(rc %in% names(avicn)) lines(2005:2030,as.numeric(avicn[,rc]),lwd=1.5,col=4,lty=5)
      lines(2005:2030,as.numeric(avicn[,cc]),lwd=3,col=4)
      
    }
    
  }
  
  
  dev.off()
  
}



# cities.possession.graphs.grubus -----------------------------------------

cities.possession.graphs.grubus <- function(bigdata,date=T,filename="",info1="",info2=""){
  
  ucp <- unique(bigdata[,c("CountryCode")])
  #   ucp <- ucp[ucp %in% tier$CityCode & !(ucp %in% bigdata$CountryCode)] # cia kad breztu tik miestus
  
  # ucp[780:810]
  # 16*ccc,h=13*ccc
  
  pdftitle <- paste0(filename,gsub(":","-",Sys.time()),".pdf")
  
  pdf(pdftitle,width=16 * 0.75,height=13 * 0.75)
  print(titlepage(info1,info2,date=Sys.time(),size=15,author="Povilas Bockus"))
  
  for(index in 1:length(ucp)){
    #   for(index in 1:10){
    
    print(ucp[index])
    cn <- ucp[index]
    
    bigxframe <- bigdata[bigdata$CountryCode==cn,]
    pnlist <- unique(bigxframe$Possession)
    
    for(pn in pnlist){
      print(pn)
      
      xframe <- bigxframe[bigxframe$Possession==pn,]
      
      maintitle <- paste(cn,pn,sep=" - ")
      
      xframe.g <- data.frame(t(xframe[,years]))
      names(xframe.g) <- xframe$CityCode
      
      
      layout(matrix(c(1,1,1,2,3,4), 2, 3, byrow = TRUE),heights = c(2.5,1))
      
      notcitynames <- xframe$CityCode[!xframe$CityCode %in% tier$CityCode]
      citynames <- xframe$CityCode[xframe$CityCode %in% tier$CityCode]
      citycodes <- intersect(xframe[xframe$CityCode %in% citynames,"CityCode"],tier$CityCode)
      if(length(notcitynames)==1) notcitynames <- rep(notcitynames,2)
      if(length(citynames)==0) citynames <- setdiff(xframe$CityCode,cn)
      if(length(citycodes)==0) citycodes <- setdiff(xframe$CityCode,cn)
      if(length(citynames)==1) citynames <- rep(citynames,2)
      if(length(citycodes)==1) citycodes <- rep(citycodes,2)
      
      g1.xframe.g <- xframe.g
      names(g1.xframe.g) <- paste(xframe$City,xframe$CityCode,sep=":")
      
      matplotf(g1.xframe.g,legend = F,abline=c(0,100),stx=2005,lty=1,turnoffpar = T,main=maintitle,text = T,col="white")
      matplotf(xframe.g[,notcitynames],legend = F,stx=2005,lty=1,turnoffpar = T,add=T,col="grey70")
      matplotf(xframe.g[,citynames],legend = F,abline=c(0,100),stx=2005,lty=1,turnoffpar = T,main=maintitle,text = T,add=T,lwd=2)
      lines(2005:2030,xframe.g[,cn],lwd=4)
      
      
      
      cola <- emi_pal()(nrow(xframe))
      names(cola) <- xframe$CityCode
      
      # w(w=17,h=12)
      
      
      # originaliu duomenu rankingas
      
      years2 <- as.character(2005:2020)
      ccrank <- c()
      codeslength <- c()
      
      
      his <- hp[hp$CountryCode==cn & hp$ProductName==pn,]
      xis <- xframe
      codes <- intersect(xis$CityCode,his$RegionCode)  
      codeslength <- length(codes)
      
      if(length(codes)==0 | !cn %in% codes) ccrank <- NULL else {
        
        his <- arrange(his[his$RegionCode %in% codes,c("RegionCode",years2)],RegionCode)
        xis <- arrange(xis[xis$CityCode %in% codes,c("CityCode",years2)],CityCode)
        
        xis[which(is.na(his),arr.ind = T)] <- NA
        
        rankmat <- data.frame(cbind(CityCode = xis[,1],Original = apply(his[,-1],1,function(x){mean(x,na.rm=T)}),
                                    Modelled = apply(xis[,-1],1,function(x){mean(x,na.rm=T)})),stringsAsFactors=F)
        rankmat[,2:3] <- apply(rankmat[,2:3],2,as.numeric)
        rankmat[,2:3] <- apply(rankmat[,2:3],2,rank)
        rankmat[,2:3] <- apply(rankmat[,2:3],2,function(x){nrow(rankmat) + 1 - x})
        
        #       ccrank <- rbind(ccrank,cbind(Possession = pn,rankmat[rankmat[,1]==cc,]))
        ccrank <- rankmat
      }
      
      if(is.null(ccrank)){
        plot(c(1,1),type="n",main="Orignal vs Modelled ranking",ylab="",xlab="",yaxt='n',xaxt='n')
        text(x=1.5,y=1,labels = "No Original Data for This City")
      } else {
        plot(ccrank[,c("Original","Modelled")],pch=20,col="grey70",
             xlim=rev(c(1,max(codeslength))),ylim=rev(c(1,max(codeslength))),main="Orignal vs Modelled ranking")
        abline(a=0,b=1,col="grey40")  
        points(ccrank[,c("Original","Modelled")],pch=20,col="grey70",cex=1.8,
               xlim=rev(c(1,max(codeslength))),ylim=rev(c(1,max(codeslength))))
        points(ccrank[ccrank$CityCode %in% citycodes,c("Original","Modelled")],pch=20,col=emi_pal()(length(citycodes)),cex=1.8,
               xlim=rev(c(1,max(codeslength))),ylim=rev(c(1,max(codeslength))))
      }
      
      
      # prognoziu rankingas
      
      ccrank2 <- c()
      codeslength2 <- c()
      
      xis <- xframe
      codeslength2 <- dim(xis)[1]
      
      rankmat <- data.frame(cbind("CityCode" = xis[,"CityCode"],"2014" = xis[,"2014"],"2030" = xis[,"2030"]),stringsAsFactors=F,check.names=F)
      rankmat[,2:3] <- apply(rankmat[,2:3],2,as.numeric)
      rankmat[,2:3] <- apply(rankmat[,2:3],2,rank)
      rankmat[,2:3] <- apply(rankmat[,2:3],2,function(x){nrow(rankmat) + 1 - x})
      
      
      plot(rankmat[,c("2014","2030")],pch=20,col=cola[ccrank2$Possession],
           xlim=rev(c(1,max(codeslength2))),ylim=rev(c(1,max(codeslength2))),main="2014 vs 2030 forecasts ranking")
      abline(a=0,b=1,col="grey40")  
      points(rankmat[,c("2014","2030")],pch=20,col="grey70",cex=1.8,
             xlim=rev(c(1,max(codeslength2))),ylim=rev(c(1,max(codeslength2))))  
      points(rankmat[rankmat$CityCode %in% citycodes,c("2014","2030")],pch=20,col=emi_pal()(length(citycodes)),cex=1.8,
             xlim=rev(c(1,max(codeslength2))),ylim=rev(c(1,max(codeslength2))))
      
      # miestu ir saliu income
      
      avicn <- avicity[avicity$CountryCode==cn,]
      avicn <- avicn[,c("CityCode","Year","AI")]
      avicn <- dcast(avicn,CityCode ~ Year,value.var="AI")
      avicn[,-1] <- t(apply(avicn[,-1],1,function(x){
        if(!all(is.na(x))) x <- MASplineVector(x)
        return(x)
      }))
      avicno <- avicn
      avicn <- data.frame(t(avicn[-1]))
      names(avicn) <- avicno[,1]
      
      inccolour <- "#00AED9"
      
      matplotf(avicn,legend = F,stx=2005,lwd=0.8,turnoffpar = T,main="Average Income",col="grey65")
      matplotf(avicn[,citycodes],add=T,turnoffpar = T,
               lwd=3,stx=2005,legend=F,text=F)
      if(cn %in% names(avicn)) lines(2005:2030,as.numeric(avicn[,cn]),lwd=4,lty=1)
      
    }
    
  }
  
  
  dev.off()
  
}


# cities.possession.graphs.grubus -----------------------------------------

cities.possession.graphs.grubus.only.after.income <- function(bigdata,date=T,filename="",info1="",info2=""){
  
  #   obigdata <- bigdata
  #   oallcountries <- bigdata$CountryCode
  #   splitmap <- unique(splitmap[,c("CityCode","GraphGroup","GraphName")])
  #   splitmap$GraphGroup <- paste(splitmap$GraphGroup,splitmap$GraphName,sep=" - ")
  # bigdata <- merge(bigdata, splitmap[,c("CountryCode","GraphGroup")], by="CountryCode", all.x = T)
  # bigdata$CountryCode <- bigdata$GraphGroup 
  #   
  
  ucp <- unique(bigdata[,c("CountryCode")])
  #   ucp <- ucp[ucp %in% tier$CityCode & !(ucp %in% bigdata$CountryCode)] # cia kad breztu tik miestus
  
  # ucp[780:810]
  # 16*ccc,h=13*ccc
  
  pdftitle <- paste0(filename,gsub(":","-",Sys.time()),".pdf")
  
  pdf(pdftitle,width=16 * 0.75,height=13 * 0.75)
  print(titlepage(info1,info2,date=Sys.time(),size=15,author="Povilas Bockus"))
  
  for(index in 1:length(ucp)){
    #   for(index in 1:10){
    
    print(ucp[index])
    cn <- ucp[index]
    
    bigxframe <- bigdata[bigdata$CountryCode==cn,]
    pnlist <- unique(bigxframe$Possession)
    
    for(pn in pnlist){
      print(pn)
      
      xframe <- bigxframe[bigxframe$Possession==pn,]
      maintitle <- paste(cn,pn,sep=" - ")
      
      xframe.g <- data.frame(t(xframe[,years]))
      names(xframe.g) <- xframe$City
      
      layout(matrix(c(1,1,1,2,3,4), 2, 3, byrow = TRUE),heights = c(2.5,1))
      
      notcitynames <- xframe$City[!xframe$CityCode %in% tier$CityCode]
      citynames <- xframe$City[xframe$CityCode %in% tier$CityCode]
      citycodes <- intersect(xframe[xframe$City %in% citynames,"CityCode"],tier$CityCode)
      if(length(notcitynames)==1) notcitynames <- rep(notcitynames,2)
      if(length(citynames)==0) citynames <- setdiff(xframe$City,"Total")
      if(length(citycodes)==0) citycodes <- setdiff(xframe$CityCode,cn)
      if(length(citynames)==1) citynames <- rep(citynames,2)
      if(length(citycodes)==1) citycodes <- rep(citycodes,2)
      
      matplotf(xframe.g,legend = F,abline=c(0,100),stx=2005,lty=1,turnoffpar = T,main=maintitle,text = T,col="white")
      matplotf(xframe.g[,notcitynames],legend = F,stx=2005,lty=1,turnoffpar = T,add=T,col="grey70")
      matplotf(xframe.g[,citynames],legend = F,abline=c(0,100),stx=2005,lty=1,turnoffpar = T,main=maintitle,text = T,add=T,lwd=2)
      lines(2005:2030,xframe.g[,"Total"],lwd=4)
      
      
      
      cola <- emi_pal()(nrow(xframe))
      names(cola) <- xframe$CityCode
      
      # w(w=17,h=12)
      
      
      # originaliu duomenu rankingas
      
      plot(1, type="n", axes=F, xlab="", ylab="")
      
      #       years2 <- as.character(2005:2020)
      #       ccrank <- c()
      #       codeslength <- c()
      #       
      #       
      #       his <- hp[hp$CountryCode==cn & hp$ProductName==pn,]
      #       xis <- xframe
      #       codes <- intersect(xis$CityCode,his$RegionCode)  
      #       codeslength <- length(codes)
      #       
      #       if(length(codes)==0 | !cc %in% codes) ccrank <- NULL else {
      #         
      #         his <- arrange(his[his$RegionCode %in% codes,c("RegionCode",years2)],RegionCode)
      #         xis <- arrange(xis[xis$CityCode %in% codes,c("CityCode",years2)],CityCode)
      #         
      #         xis[which(is.na(his),arr.ind = T)] <- NA
      #         
      #         rankmat <- data.frame(cbind(CityCode = xis[,1],Original = apply(his[,-1],1,function(x){mean(x,na.rm=T)}),
      #                                     Modelled = apply(xis[,-1],1,function(x){mean(x,na.rm=T)})),stringsAsFactors=F)
      #         rankmat[,2:3] <- apply(rankmat[,2:3],2,as.numeric)
      #         rankmat[,2:3] <- apply(rankmat[,2:3],2,rank)
      #         rankmat[,2:3] <- apply(rankmat[,2:3],2,function(x){nrow(rankmat) + 1 - x})
      #         
      #         #       ccrank <- rbind(ccrank,cbind(Possession = pn,rankmat[rankmat[,1]==cc,]))
      #         ccrank <- rankmat
      #       }
      #       
      #       if(is.null(ccrank)){
      #         plot(c(1,1),type="n",main="Orignal vs Modelled ranking",ylab="",xlab="",yaxt='n',xaxt='n')
      #         text(x=1.5,y=1,labels = "No Original Data for This City")
      #       } else {
      #         plot(ccrank[,c("Original","Modelled")],pch=20,col="grey70",
      #              xlim=rev(c(1,max(codeslength))),ylim=rev(c(1,max(codeslength))),main="Orignal vs Modelled ranking")
      #         abline(a=0,b=1,col="grey40")  
      #         points(ccrank[,c("Original","Modelled")],pch=20,col="grey70",cex=1.8,
      #                xlim=rev(c(1,max(codeslength))),ylim=rev(c(1,max(codeslength))))
      #         points(ccrank[ccrank$CityCode %in% citycodes,c("Original","Modelled")],pch=20,col=emi_pal()(length(citycodes)),cex=1.8,
      #                xlim=rev(c(1,max(codeslength))),ylim=rev(c(1,max(codeslength))))
      #       }
      
      
      # prognoziu rankingas
      
      ccrank2 <- c()
      codeslength2 <- c()
      
      xis <- xframe
      codeslength2 <- dim(xis)[1]
      
      rankmat <- data.frame(cbind("CityCode" = xis[,"CityCode"],"2014" = xis[,"2014"],"2030" = xis[,"2030"]),stringsAsFactors=F,check.names=F)
      rankmat[,2:3] <- apply(rankmat[,2:3],2,as.numeric)
      rankmat[,2:3] <- apply(rankmat[,2:3],2,rank)
      rankmat[,2:3] <- apply(rankmat[,2:3],2,function(x){nrow(rankmat) + 1 - x})
      
      
      plot(rankmat[,c("2014","2030")],pch=20,col=cola[ccrank2$Possession],
           xlim=rev(c(1,max(codeslength2))),ylim=rev(c(1,max(codeslength2))),main="2014 vs 2030 forecasts ranking")
      abline(a=0,b=1,col="grey40")  
      points(rankmat[,c("2014","2030")],pch=20,col="grey70",cex=1.8,
             xlim=rev(c(1,max(codeslength2))),ylim=rev(c(1,max(codeslength2))))  
      points(rankmat[rankmat$CityCode %in% citycodes,c("2014","2030")],pch=20,col=emi_pal()(length(citycodes)),cex=1.8,
             xlim=rev(c(1,max(codeslength2))),ylim=rev(c(1,max(codeslength2))))
      
      # miestu ir saliu income
      
      avicn <- avicity[avicity$CountryCode==cn,]
      avicn <- avicn[,c("CityCode","Year","AI")]
      avicn <- dcast(avicn,CityCode ~ Year,value.var="AI")
      avicn[,-1] <- t(apply(avicn[,-1],1,function(x){
        if(!all(is.na(x))) x <- MASplineVector(x)
        return(x)
      }))
      avicno <- avicn
      avicn <- data.frame(t(avicn[-1]))
      names(avicn) <- avicno[,1]
      
      inccolour <- "#00AED9"
      
      matplotf(avicn,legend = F,stx=2005,lwd=0.8,turnoffpar = T,main="Average Income",col="grey65")
      matplotf(avicn[,citycodes],add=T,turnoffpar = T,
               lwd=3,stx=2005,legend=F,text=F)
      if(cn %in% names(avicn)) lines(2005:2030,as.numeric(avicn[,cn]),lwd=4,lty=1)
      
    }
    
  }
  
  
  dev.off()
  
}






# cities.possession.graphs.grubus -----------------------------------------

average.income.graphs.grubus <- function(bigdata,date=T,filename="",info1="",info2=""){
  
  ucp <- unique(bigdata[,c("CountryCode")])
  #   ucp <- ucp[ucp %in% tier$CityCode & !(ucp %in% bigdata$CountryCode)] # cia kad breztu tik miestus
  
  # ucp[780:810]
  # 16*ccc,h=13*ccc
  
  pdftitle <- paste0(filename,gsub(":","-",Sys.time()),".pdf")
  
  pdf(pdftitle,width=16 * 0.75,height=10 * 0.75)
  print(titlepage(info1,info2,date=Sys.time(),size=15,author="Povilas Bockus"))
  
  for(index in 1:length(ucp)){
    #   for(index in 1:10){
    
    print(ucp[index])
    cn <- ucp[index]
    
    bigxframe <- bigdata[bigdata$CountryCode==cn,]
    citycodes <- bigxframe$CityCode[bigxframe$CityCode %in% tier$CityCode]  
    
    avicn <- avicity[avicity$CountryCode==cn,]
    #     nnames <- unique(avicn$City)
    
    avicn <- avicn[,c("CityCode","City","Year","AI")]
    avicn <- dcast(avicn,CityCode + City ~ Year,value.var="AI")
    avicn[,-1:-2] <- t(apply(avicn[,-1:-2],1,function(x){
      if(!all(is.na(x))) x <- MASplineVector(x)
      return(x)
    }))
    avicno <- avicn
    avicn <- data.frame(t(avicn[-1:-2]))
    names(avicn) <- avicno[,2]
    avicn2 <- avicn
    names(avicn2) <- avicno[,1]
    
    inccolour <- "#00AED9"
    
    matplotf(avicn,legend = F,stx=2005,lwd=0.8,turnoffpar = T,main=paste0(cn," - Average Income"),col="grey65",text=T)
    if(ncol(data.frame(avicn2[,citycodes]))>1){
      matplotf(avicn2[,citycodes],add=T,turnoffpar = T,
               lwd=3,stx=2005,legend=F,text=F)
    }
    if(ncol(data.frame(avicn2[,citycodes]))==1){
      matplotf(data.frame(avicn2[,citycodes]),add=T,turnoffpar = T,
               lwd=3,stx=2005,legend=F,text=F)
    } 
    if(cn %in% names(avicn)) lines(2005:2030,as.numeric(avicn2[,cn]),lwd=4,lty=1)
    
  }
  
  
  dev.off()
  
}




# cities.possession.graphs.grubus -----------------------------------------

abq.summary.graphs.grubus <- function(kof,date=T,filename="",info1="",info2=""){
  
  ucp <- unique(kof[,c("CountryCode")])
  #   ucp <- ucp[ucp %in% tier$CityCode & !(ucp %in% bigdata$CountryCode)] # cia kad breztu tik miestus
  
  # ucp[780:810]
  # 16*ccc,h=13*ccc
  
  pdftitle <- paste0(filename,gsub(":","-",Sys.time()),".pdf")
  
  pdf(pdftitle,width=16 * 0.75,height=12 * 0.75)
  print(titlepage(info1,info2,date=Sys.time(),size=15,author="Povilas Bockus"))
  
  for(index in 1:length(ucp)){
    
    print(ucp[index])
    cn <- ucp[index]
    
    bigxframe <- kof[kof$CountryCode==cn,]
    
    layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
    
    # average 
    
    avg.vector <- ddply(bigxframe,~CityCode,function(x){
      x<<-x
      print(as.character(x$CityCode[1]))
      ttt <- b.converter.reverse(x[,"a"],x[,"b.cntUSD"],x[,"q"])
      out <- data.frame(cbind(x[,1:4],AI=ttt),check.names=F,stringsAsFactors=F)
      return(out)
    })
    
    avg.vector <- dcast(avg.vector,CityCode + City + CountryCode ~ Year,value.var="AI")
    avg.vector.o <- avg.vector
    avg.vector <- data.frame(t(avg.vector[,years]))
    names(avg.vector) <- avg.vector.o$City
    
    matplotf(avg.vector,legend = F,stx=2005,lwd=0.8,turnoffpar = T,main=paste0(cn," - Average Income"),text=T)
    if(cn %in% avg.vector.o$CityCode) lines(2005:2030,as.numeric(avg.vector.o[avg.vector.o$CityCode==cn,years]),lwd=4,lty=1)
    
    # median income
    
    median.vector <- ddply(bigxframe,~CityCode,function(x){
      x<<-x
      print(as.character(x$CityCode[1]))
      m <- c()
      for(y in 1:length(years)){m[y] <- qsm(0.5,as.numeric(x[y,c("a","b.cntUSD","q")]))}
      ttt <- m
      out <- data.frame(cbind(x[,1:4],Median=ttt),check.names=F,stringsAsFactors=F)
      return(out)
    })
    
    median.vector <- dcast(median.vector,CityCode + City + CountryCode ~ Year,value.var="Median")
    median.vector.o <- median.vector
    median.vector <- data.frame(t(median.vector[,years]))
    names(median.vector) <- median.vector.o$City
    
    matplotf(median.vector,legend = F,stx=2005,lwd=0.8,turnoffpar = T,main=paste0(cn," - Median Income"),text=T)
    if(cn %in% median.vector.o$CityCode) lines(2005:2030,as.numeric(median.vector.o[median.vector.o$CityCode==cn,years]),lwd=4,lty=1)
    
    
    
    # Mode income
    
    mode.vector <- ddply(bigxframe,~CityCode,function(x){
      x<<-x
      print(as.character(x$CityCode[1]))
      m <- c()
      for(y in 1:length(years)){m[y] <- modesm(as.numeric(x[y,c("a","b.cntUSD","q")]))}
      ttt <- m
      out <- data.frame(cbind(x[,1:4],Mode=ttt),check.names=F,stringsAsFactors=F)
      return(out)
    })
    
    mode.vector <- dcast(mode.vector,CityCode + City + CountryCode ~ Year,value.var="Mode")
    mode.vector.o <- mode.vector
    mode.vector <- data.frame(t(mode.vector[,years]))
    names(mode.vector) <- mode.vector.o$City
    
    matplotf(mode.vector,legend = F,stx=2005,lwd=0.8,turnoffpar = T,main=paste0(cn," - Mode Income"),text=T)
    if(cn %in% mode.vector.o$CityCode) lines(2005:2030,as.numeric(mode.vector.o[mode.vector.o$CityCode==cn,years]),lwd=4,lty=1)
    
    
    #all in one
    
    matplotf(avg.vector,legend = F,stx=2005,lwd=0.8,turnoffpar = T,main=paste0(cn," - All in One"),text=F,col="dodgerblue4",abline=range(rbind(avg.vector,median.vector,mode.vector)))
    if(cn %in% avg.vector.o$CityCode) lines(2005:2030,as.numeric(avg.vector.o[avg.vector.o$CityCode==cn,years]),lwd=5,lty=1,col="blue")
    
    matplotf(median.vector,legend = F,stx=2005,lwd=0.8,turnoffpar = T,main=paste0(cn," - Median Income"),text=F,col="firebrick3",add=T)
    if(cn %in% median.vector.o$CityCode) lines(2005:2030,as.numeric(median.vector.o[median.vector.o$CityCode==cn,years]),lwd=5,lty=1,col="red")
    
    matplotf(mode.vector,legend = F,stx=2005,lwd=0.8,turnoffpar = T,main=paste0(cn," - Mode Income"),text=F,col="goldenrod",add=T)
    if(cn %in% mode.vector.o$CityCode) lines(2005:2030,as.numeric(mode.vector.o[mode.vector.o$CityCode==cn,years]),lwd=5,lty=1,col="orange")
    
    legend("topleft",col=c("dodgerblue4","firebrick3","goldenrod"),lty=1,legend=c("Average","Median","Mode"),cex=0.5)
    
    abline(h=range(rbind(avg.vector,median.vector,mode.vector)),col="white")
    
  }
  dev.off()
  
}



# original.data.barplots --------------------------------------------------

original.data.barplots <- function(hp,avicity,date=T,filename="",info1="",info2=""){
  
  
  pdftitle <- paste0(filename,gsub(":","-",Sys.time()),".pdf")
  
  pdf(pdftitle,width=16 * 0.75,height=13 * 0.75)
  print(titlepage(info1,info2,date=Sys.time(),size=15,author="Povilas Bockus"))
  
  
  for(cn in sort(unique(hp$CountryCode))){
    print(cn)
    
    #     avis <- avicity[avicity$CountryCode==cn,]
    #     avis <- ddply(avis,~CityCode + City,function(x){mean(x$AI)})
    
    his <- hp[hp$CountryCode==cn,]
    ohis <- his
    #     notnayears <- years2[apply(his[,years2],2,function(x){!all(is.na(x))})]
    
    his <- cbind(his,meanvalue=apply(his[,years2],1,function(x){mean(x,na.rm=T)}))
    his <- his[,c("RegionCode","CountryName","ProductName","meanvalue")]
    
    #     avis <- arrange(avis,CityCode)
    his <- arrange(his,RegionCode)
    
    histemplate <- unique(his[,c("RegionCode","CountryName")])
    
    ###
    
    n <- length(unique(his$ProductName)) + 2
    if(n %in% 1:3) nnrow <- 1
    if(n %in% 4:6) nnrow <- 2
    if(n %in% 7:12) nnrow <- 3
    if(n %in% 13:20) nnrow <- 4
    if(n %in% 21:30) nnrow <- 5
    if(n %in% 31:42) nnrow <- 6
    if(n %in% 43:56) nnrow <- 7
    if(n > 56) nnrow <- 8
    
    
    #     par(mfrow=c(nnrow,ceiling(n/nnrow)))
    
    ###  
    
    #     colinc <- rep("firebrick2",dim(avis)[1])
    #     colinc[which(avis[,1]==cn)] <- "darkred"
    #     colinc[substr(avis$CityCode,1,1)=="R" & !is.na(as.numeric(substr(avis$CityCode,2,2)))] <- "darkorange"
    #     
    #     barplot(avis[,3],names.arg = avis[,2],col=colinc, las=3,main=paste(cn,"Average Income",sep=" - "),cex.names=0.4)
    #     abline(h=avis[,3][which(avis[,1]==cn)],col="darkred",lwd=2)
    #     abline(h=c(0),col="grey60",lty=2)
    #     
    
    for(pn in unique(his$ProductName)){
      
      notnayears <- years2[apply(ohis[ohis$ProductName==pn,years2],2,function(x){!all(is.na(x))})]
      
      hhis <- his[his$ProductName==pn,c("RegionCode","meanvalue")]
      #       hhisx <- avis
      hhisx <- histemplate
      
      hhisx <- merge(hhisx,hhis,by.x="RegionCode",by.y="RegionCode",all.x=T)
      astosttotal <- mean(as.numeric(rposout[rposout$CountryCode==cn & rposout$CityCode==cn & rposout$Possession==pn,notnayears]))
      hhisx <- data.frame(rbind(c(paste(cn,"Astos"),paste(cn,"Astos"),astosttotal),hhisx))
      hhisx$meanvalue <- as.numeric(hhisx$meanvalue)
      
      colinc2 <- rep("#FAA537",dim(hhisx)[1]+1)
      colinc2[which(hhisx[,1]==cn)] <- "darkblue"
      colinc2[which(hhisx[,1]==paste(cn,"Astos"))] <- "dodgerblue2"
      colinc2[substr(hhisx$RegionCode,1,1)=="R" & !is.na(as.numeric(substr(hhisx$RegionCode,2,2)))] <- "chartreuse3"
      colinc2[which(hhisx$meanvalue > 100 | hhisx$meanvalue < 0)] <- "red"
      
      ylimas <- c(min(c(0,hhisx$meanvalue),na.rm=T),max(c(100,hhisx$meanvalue),na.rm=T))
      
      par(mar=c(9, 4, 4, 2) + 0.1)
      barplot(hhisx[,"meanvalue"],names.arg = hhisx[,2],col=colinc2, las=3,main=paste(cn,pn,sep=" - "),ylim=ylimas,cex.names=0.8)
      abline(h=hhisx[,"meanvalue"][which(hhisx[,1]==paste(cn,"Astos"))],col="dodgerblue2",lwd=2)
      abline(h=hhisx[,"meanvalue"][which(hhisx[,1]==cn)],col="darkblue",lwd=2)
      abline(h=c(0,100),col="grey60",lty=2)
      
    }
    
    #     plot.new()
    #     legend("topleft",legend=paste(hhisx[,1],hhisx[,2],sep=" - "),cex=1)
    
    
  }
  dev.off()
  
}



lastinput <- function(dirname){
  files <- dir(dirname)
  out <- files[which.max(as.numeric(gsub("\\D","",files)))]
  cat(paste0("Loading: ",paste0(dirname,"/",out)),"\n")
  return(paste0(dirname,"/",out))
}


modesm <- function(theta,n=1000){ 
  arg <- seq(0,theta[2]*4,length.out=n)
  return(arg[which.max(dsm(arg,theta))]) 
}


posfun.reverse <- function(x){
  y <- splinefun(pframe,method = "monoH.FC")(x)
  return(y)
}






mono.pos <- function(out,tot,...){
  
  ratio <- out/tot
  ratio <- sm(ratio,df=length(ratio)*0.7,plot=F)
  out <- ratio * tot
  
  way <- NULL
  if(sum(diff(tot)>0)>=sum(diff(tot)<0)) way  <- "up" else way <- "down"
  
  if(way == "down"){
    break.point <- max(which(diff(tot) > 0) + 1)
    if(break.point!=-Inf){
      out[-1:-break.point] <- mv(out[-1:-break.point],way=way,...)
    } else {
      out <- mv(out,way=way,...)
    }
  } 
  
  if(way == "up"){
    break.point <- max(which(diff(tot) < 0) + 1)
    if(break.point!=-Inf){
      out[-1:-break.point] <- mv(out[-1:-break.point],way=way,...)
    } else {
      out <- mv(out,way=way,...)
    }
  } 
  
  return(out)
}





# cities.possession.graphs.grubus -----------------------------------------

cities.possession.graphs.grubus.tik.miestai.split <- function(bigdata,splitmap,date=T,filename="",info1="",info2=""){
  
  obigdata <- bigdata
  bigdata <- bigdata[bigdata$CityCode %in% splitmap$CityCode,]
  
  
  ucp <- as.character(unique(splitmap[,c("GraphGroup")]))
  #   ucp <- ucp[ucp %in% tier$CityCode & !(ucp %in% bigdata$CountryCode)] # cia kad breztu tik miestus
  
  # ucp[780:810]
  # 16*ccc,h=13*ccc
  
  pdftitle <- paste0(filename,gsub(":","-",Sys.time()),".pdf")
  
  pdf(pdftitle,width=16 * 0.8,height=13 * 0.7)
  print(titlepage(info1,info2,date=Sys.time(),size=15,author="Povilas Bockus"))
  
  for(index in 1:length(ucp)){
    #   for(index in 1:10){
    
    print(ucp[index])
    cn <- ucp[index]
    
    citiessplit <- as.character(splitmap[splitmap$GraphGroup %in% cn,"CityCode"])
    grname <- as.character(splitmap[splitmap$GraphGroup %in% cn,"GraphName"][1])
    
    bigxframe <- bigdata[bigdata$CityCode %in% citiessplit,]
    bigxframeall <- bigdata[bigdata$CountryCode == bigxframe$CountryCode[1],]
    
    pnlist <- unique(bigxframe$Possession)
    
    for(pn in pnlist){
      print(pn)
      
      xframe <- bigxframe[bigxframe$Possession==pn,]
      
      maintitle <- paste(cn,grname,pn,sep=" - ")
      
      xframe.g <- data.frame(t(xframe[,years]))
      names(xframe.g) <- xframe$City
      
      
      xframeall <-  bigxframeall[bigxframeall$Possession==pn,]     
      matplotf(t(xframeall[,years]),legend = F,abline=c(0,100),stx=2005,lty=1,turnoffpar = T,main=maintitle,text = F,col="grey90")
      
      matplotf(xframe.g,legend = F,abline=c(0,100),stx=2005,lty=1,turnoffpar = T,main=maintitle,text = T,add=T)
      lines(2005:2030,xframe.g[,"Total"],lwd=4)    
      
    }
    
  }
  
  
  dev.off()
  
}




add.cities.which.are.equal.to.countries <- function(cc,citycode,cityname,pd,rposout){
  
  xframe <- pd[pd$CountryCode==cc,]
  
  xframe <- ddply(xframe,~ CountryCode + Possession,function(x){
    x <<- x
    # })
    cat(paste0("Calculating Country Totals: ",x$CountryCode[1]),"\n")
    out <- x[1,]
    out[,c("CityCode","City")] <- rep("Total",2)
    out[years] <- ddply(pd[pd$CountryCode==out$CountryCode & pd$Possession==out$Possession,c("Year","est.Pos.norm")],
                        ~Year,function(y){mean(y[,"est.Pos.norm"])})$V1
    out$CityCode <- out$CountryCode
    return(out)
  })
  
  xframe <- data.frame(cbind(Possession = xframe$Possession,CountryCode = xframe$CountryCode, CityCode = citycode,
                             City = cityname, xframe[,years],ProductID = xframe$ProductID),check.names=F)
  
  xframe <- xframe[xframe$ProductID %in% unique(rposout$ProductID),]
  return(xframe)
  
}



# compare2data.frames -----------------------------------------------------

compare2data.frames <- function(xframe.new,xframe.old,datayears = as.character(2005:2030)){
  
  id.name.new <- names(xframe.new)[-which(names(xframe.new) %in% datayears)]
  id.name.old <- names(xframe.old)[-which(names(xframe.old) %in% datayears)]
  
  xframe.new <- cbind(xframe.new,key = apply(xframe.new[,id.name.new],1,function(x){paste(x,collapse = "-")}))
  xframe.old <- cbind(xframe.old,key = apply(xframe.old[,id.name.old],1,function(x){paste(x,collapse = "-")}))
  
  if(all(xframe.new$key %in% xframe.old$key) & all(xframe.old$key %in% xframe.new$key)){
    print("Data Frames can be compared")
  }
  
  xframe.new <- arrange(xframe.new,key)
  xframe.old <- arrange(xframe.old,key)
  
  ratioframe <- xframe.new
  
  ratioframe[,datayears] <- round(xframe.new[,datayears]/xframe.old[,datayears],5)
  changedframe <- which(ratioframe[,datayears]!=1,arr.ind = T)
  changedframe <- unique(changedframe[,1])
  
  outframe <- rbind(cbind(Version = "New",xframe.new[changedframe,]),
        cbind(Version = "Old",xframe.old[changedframe,]))
  
  return(outframe)

}


####city orginalaus tasko (ratio regiono ir city) radimas kai turime salies possessions by decile ir income bei regiono orginalu taska
## 1. issibreziam total salies kreive (x asis - income by decile, y-possessions by decile), metu, kuriu orginalu taska mes turim
## 2. pasiimame regional orginalu taska ir randame "jo pajamas" ant 1. kreives
## 3. randam ratio tarp is abqUSD isskaiciuotu tu metu vidutiniu regiono ir miesto pajamu
## 4. ta ratio pritaikom fiktyviom regiono pajamom (2. punktas) ir gaunam fiktyvias miesto pajamas
## 5. pagal fiktyvias miesto pajamas ir total kreive (1. punktas) randam "orginalu city taska", 
##    prie kurio veliau pritraukinesim nauja is income gauta possessions kreive
## 6. ratio.out - (tam tikru metu pos city (suskaiciuotas naudojantis regional org.))/regional orig
fromIncomePossessionCurve <- function(city, cc, pn , rposout, hp, abqUSD, final.pbd, hp.cc, inc.cc){
  
#   hp.cc <- hp[hp$CountryCode==cc & hp$ProductName==pn,]
#   inc.cc <- rposout[rposout$CountryCode==cc & rposout$Possession==pn,]
  year <- apply(hp.cc[,as.character(2000:2020)], 2, function(x){any(!is.na(x)) })
  year <- names(year[which(year==T)])
  year <- year[length(year)]
  
  if (length(year)==0){
    ratio.out=NA
    return.out(ratio.out)
    break
  }
  
  #is income sumodeliuoti regioniniai duomenys pritempiami prie orginaliu
  out <- NULL
  for (i in unique(inc.cc$CityCode)){
    
    s <- subset(inc.cc, CityCode==i)
    if (substr(i, 1, 2) %in% c("R0", "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9")){
      
      hpsh <- hp.cc[hp.cc$RegionCode==i,]
      hpsh <- melt(hpsh, id=c("RegionCode", "Region/CityName", "CountryName", "CountryCode", "ProductID", "ProductName", "Unit"))
      
      if (any(!is.na(hpsh$value))){
        point <- nrow(hpsh[!is.na(hpsh$value),])
        value <- hpsh[!is.na(hpsh$value),][point,]$value
        skirt <- value/eval(parse(text=paste0("s$`", year, "`")))
        s[,as.character(2005:2030)] <- s[,as.character(2005:2030)]*skirt
        out <- rbind(s, out)
      }
    }
  }
  
  if (length(out)!=0){
    
    income.pos <- final.pbd[final.pbd$CountryCode==cc & final.pbd$Possession==pn & final.pbd$Year==year,]
    income.pos <- income.pos[order(income.pos$Decile.No),]
    income.pos <- income.pos[, c("est.Pos.norm", "Average.income")]
    
    # plot(income.pos$Average.income, income.pos$est.Pos.norm, type="o")
    # points(income.pos$Average.income, income.pos$est.Pos.norm, col="red")
    
    #suskaciuojam nauja city lygi ir prie jo pritempiam is miestu income gautas possession kreives
    new.city <- NULL
    
    s <- subset(inc.cc, CityCode==city)
    
    reg <- as.character(city.codes[city.codes$CityCode==city,]$RegionCode)
    reg <- substr(reg, 1,4)
    pos.reg <- out[out$CityCode==reg,year]
    inc.reg <- mean(sm.tft(as.numeric(abqUSD[abqUSD$CityCode==reg & abqUSD$Year==year, c("a", "b.cntUSD", "q")])))
    inc.city <- mean(sm.tft(as.numeric(abqUSD[abqUSD$CityCode==city & abqUSD$Year==year, c("a", "b.cntUSD", "q")])))
    ratio <- inc.city/inc.reg
    
    if (round(ratio, 3)==1){
      
      s[,as.character(2005:2030)] <- out[out$CityCode==reg,as.character(2005:2030)]
      ratio.out <- 1
    }else{
      
      income.pos1 <- rbind(income.pos, data.frame(est.Pos.norm=c(pos.reg), Average.income=c(NA)))
      income.pos1 <- income.pos1[order(income.pos1$est.Pos.norm),]
      
      if (is.na(income.pos1$Average.income[nrow(income.pos1)]) | is.na(income.pos1$Average.income[1]) |
          (!(all(income.pos$est.Pos.norm==sort(income.pos$est.Pos.norm)))  & 
          (!all(income.pos$est.Pos.norm==sort(income.pos$est.Pos.norm, decreasing = T))))){
        
        ratio.out <- NA
      }else{
        
        x1 <- income.pos1[which(is.na(income.pos1$Average.income))-1,]$Average.income
        x2 <- income.pos1[which(is.na(income.pos1$Average.income))+1,]$Average.income
        y1 <- income.pos1[which(is.na(income.pos1$Average.income))-1,]$est.Pos.norm
        y2 <- income.pos1[which(is.na(income.pos1$Average.income))+1,]$est.Pos.norm
        k <- (y1-y2)/(x1-x2)
        b <- y2-x2*(y1-y2)/(x1-x2)
        
        inc.reg.f <- (pos.reg-b)/k
        inc.city.f <- ratio*inc.reg.f
        
        if (min(income.pos$Average.income)>inc.city.f){
          
          pos.city <- min(income.pos$est.Pos.norm)
        }else{
        
          income.pos2 <- rbind(income.pos, data.frame(est.Pos.norm=c(NA), Average.income=c(inc.city.f)))
          income.pos2 <- income.pos2[order(income.pos2$Average.income),]
          x1 <- income.pos2[which(is.na(income.pos2$est.Pos.norm))-1,]$Average.income
          x2 <- income.pos2[which(is.na(income.pos2$est.Pos.norm))+1,]$Average.income
          y1 <- income.pos2[which(is.na(income.pos2$est.Pos.norm))-1,]$est.Pos.norm
          y2 <- income.pos2[which(is.na(income.pos2$est.Pos.norm))+1,]$est.Pos.norm
          
          if (is.na(x2)){
            x2 <- income.pos2[which(is.na(income.pos2$est.Pos.norm))-2,]$Average.income
            y2 <- income.pos2[which(is.na(income.pos2$est.Pos.norm))-2,]$est.Pos.norm
          }
          
          k <- (y1-y2)/(x1-x2)
          b <- y2-x2*(y1-y2)/(x1-x2)
          
          #pos.city naujas pos tam tikriems metams pritraukus prie regiono orginaliu
          pos.city <- b+k*inc.city.f
        }
        
        ratio.out <- pos.city/pos.reg 
      }
    }
  }else{
    
    ratio.out <- NA
  }
  
  return(ratio.out)
}


# sarasas, tokiu atveju, kur pavertus absoliuciomis reiksmemis 
# susumuoti miestu possession of durables virsija salies possession of durables

list.cities.surpass.country <- function(hh, rposout){
  
  rposout <- rposout[!(substr(rposout$CityCode, 1, 2) %in% paste0("R", c(0:9))),]
  
  hh <- hh[,c("CityCode", "ParentCountryCode", paste0("Y", 2005:2030))]
  hh <- rename(hh, c(ParentCountryCode="CountryCode"))
  hh[hh$CityCode=="",]$CityCode <- hh[hh$CityCode=="",]$CountryCode
  hh <- melt(hh,  id=c("CityCode", "CountryCode"))
  hh <- rename(hh, c(variable="Year", value="hh"))
  hh$Year <- gsub("Y", "", hh$Year)
  
  rposout <- melt(rposout, id=c("CountryCode", "CityCode", "City", "ProductID", "Possession"))
  rposout <- rename(rposout, c(variable="Year", value="poss"))
  
  rposout <- merge(rposout, hh, by=c("CityCode", "CountryCode", "Year"), all.x=T)
  head(rposout[is.na(rposout$hh),]) ##has to be empty
  
  rposout$hh.poss <- rposout$poss*rposout$hh/100
  
  n.rposout <- rposout[rposout$City=="Total",]
  c.rposout <- rposout[rposout$City!="Total",]
  
  n.rposout <- rename(n.rposout, c(poss="n.poss", hh="n.hh", hh.poss="n.hh.poss"))
  n.rposout <- n.rposout[,c("CountryCode", "Year", "ProductID", "n.poss", "n.hh", "n.hh.poss")]
  
  c.rposout <- merge(c.rposout, n.rposout, by=c("CountryCode", "ProductID", "Year"))
  
  out <- ddply(c.rposout, .(CountryCode, ProductID, Year, Possession), summarize, 
               city.sum= sum(hh.poss), n.hh.poss=mean(n.hh.poss))
  
  return (unique(out[out$city.sum>out$n.hh.poss, c("CountryCode", "Possession", "ProductID")]))
  
}


##taisymnas

fixing.cities.surpass.country <- function(hh, rposout, plot=T){
  
  year <- as.character(2005:2030)
  
  rposout.out <- rposout
  rposout <- rposout[!(substr(rposout$CityCode, 1, 2) %in% paste0("R", c(0:9))),]
  
  hh <- hh[,c("CityCode", "ParentCountryCode", paste0("Y", 2005:2030))]
  hh <- rename(hh, c(ParentCountryCode="CountryCode"))
  hh[hh$CityCode=="",]$CityCode <- hh[hh$CityCode=="",]$CountryCode
  hh <- melt(hh,  id=c("CityCode", "CountryCode"))
  hh <- rename(hh, c(variable="Year", value="hh"))
  hh$Year <- gsub("Y", "", hh$Year)
  
  rposout <- melt(rposout, id=c("CountryCode", "CityCode", "City", "ProductID", "Possession"))
  rposout <- rename(rposout, c(variable="Year", value="poss"))
  
  rposout <- merge(rposout, hh, by=c("CityCode", "CountryCode", "Year"), all.x=T)
  head(rposout[is.na(rposout$hh),]) ##has to be empty
  
  rposout$hh.poss <- rposout$poss*rposout$hh/100
  
  n.rposout <- rposout[rposout$City=="Total",]
  c.rposout <- rposout[rposout$City!="Total",]
  
  n.rposout <- rename(n.rposout, c(poss="n.poss", hh="n.hh", hh.poss="n.hh.poss"))
  n.rposout <- n.rposout[,c("CountryCode", "Year", "ProductID", "n.poss", "n.hh", "n.hh.poss")]
  
  c.rposout <- merge(c.rposout, n.rposout, by=c("CountryCode", "ProductID", "Year"))
  
  out <- ddply(c.rposout, .(CountryCode, ProductID, Year, Possession), summarize, 
               city.sum= sum(hh.poss), n.hh.poss=mean(n.hh.poss))
  
  bad.list <- unique(out[out$city.sum>out$n.hh.poss, c("CountryCode", "Possession", "ProductID")])
  
  corrected <- NULL
  for (i in 1:nrow(bad.list)){
    
    s <- c.rposout[c.rposout$CountryCode==bad.list[i,"CountryCode"] & c.rposout$ProductID==bad.list[i,"ProductID"],]
    
    ss <- ddply(s, .(CountryCode, ProductID, Year), summarize,
          city.hh=sum(hh), city.hh.poss=sum(hh.poss))  
   
    s <- merge(s, ss, by=c("CountryCode", "ProductID", "Year"))
    s$hh.left <- s$n.hh-s$city.hh
    s$hh.poss.left <- s$hh.left*(s$n.poss/2)/100
    s$new.poss <- (s$hh.poss*s$n.hh.poss)/(s$city.hh.poss+s$hh.poss.left)/s$hh*100
    
    if (length(s[s$new.poss>s$poss,]$new.poss)!=0){
      
      s[s$new.poss>s$poss,]$new.poss <- s[s$new.poss>s$poss,]$poss
    }
    corrected <- rbind(s, corrected)
    
    for (kk in unique(s$CityCode)){
      
      s.kk <- s[s$CityCode==kk,]
      s.kk <- s.kk[order(s.kk$Year),]
      rposout.out[rposout.out$CountryCode==bad.list[i,"CountryCode"] & 
                    rposout.out$ProductID==bad.list[i,"ProductID"] &
                    rposout.out$CityCode==kk,year] <- s.kk$new.poss
      
    }
    
  }
  
  if (plot){
    
    for.graph <- corrected[,c("CountryCode", "ProductID", "Possession", "CityCode", "City", "Year", 
                              "n.poss", "poss", "new.poss")]
    
    filename="plots/5.A. Forecasted possessions after norming (where was necessary)/Forecasted possessions after norming_"
    
    pdftitle <- paste0(filename,gsub(":","-",Sys.time()),".pdf")
    
    pdf(pdftitle,width=16 * 0.75,height=13 * 0.75)
    print(titlepage("Forecasted possessions after norming",date=Sys.time(),size=15,author="Euromonitor"))
    par(ps = 10, lab = c(14,10,10))
    
    for (i in unique(for.graph$CountryCode)){
      
      sg <- for.graph[for.graph$CountryCode==i,]
      
      for (j in unique(sg$Possession)){
      
        sgg <- sg[sg$Possession==j,]
        
        c.sgg <- dcast(sgg[,c("City", "Year", "poss")], City~Year, value.var = c("poss"))
        n.sgg <- dcast(sgg[,c("City", "Year", "n.poss")], City~Year, value.var = c("n.poss"))[1,year]
        new.sgg <- dcast(sgg[,c("City", "Year", "new.poss")], City~Year, value.var = c("new.poss"))
        
        matplot(c(2005:2030), t(new.sgg[,year]), lty=1, type="l", lwd = 1, col=emi_pal()(nrow(new.sgg)),
                main=paste(i, j, sep=" - "), xlab=c("Year"), ylab=c("Possession"), xlim=c(2004.5, 2030.5),
                ylim=c(1, 100))
        matlines(c(2005:2030), t(c.sgg[,year]), lty=2, col=emi_pal()(nrow(c.sgg)))
        lines(c(2005:2030), t(n.sgg[,year]), col="black", lwd=3)
        grid()
        text(rep(2030.5, nrow(new.sgg)), t(new.sgg[,year])[26,], labels=new.sgg[,1], cex=0.8, col=emi_pal()(nrow(new.sgg)))
        text(rep(2004.5, nrow(new.sgg)), t(new.sgg[,year])[1,], labels=new.sgg[,1], cex=0.8, col=emi_pal()(nrow(new.sgg)))
      }
    }
    
    dev.off()
  }
  
  return(rposout.out)
}


lastinput.zip <- function(dirname){
  files <- dir(dirname)
  out <- files[which.max(as.numeric(gsub("\\D","",files)))]
  out2 <- gsub(".zip", ".csv", out)
  #   cat(paste0("Loading: ",paste0(dirname,"/",out)),"\n")
  paste0("unz(","'",dirname, "/", out, "'", ",","'", out2, "')")
  return(paste0("unz(","'",dirname, "/", out, "'", ",","'", out2, "')"))
}


issue.smooth <- function(bt, rposout){
  issue.sm <- bt[bt$Issue=="smooth",]
  
  for(index in 1:nrow(issue.sm)){
    
    x <- issue.sm[index,]
    cn <- x$CountryCode[1]
    pn <- x$ProductName[1]
    city <- x$City[1]
    
    if(x$City == "all"){
      
      xframe <- rposout[rposout$CountryCode==x$CountryCode & rposout$Possession==x$ProductName,]
      xframeratio <- xframe
      tot <- as.numeric(xframe[xframe$City == "Total",years])
      xframeratio[,years] <- t(apply(xframe[,years],1,function(x){x/tot}))
      
      matplotf(t(xframe[years]),main=paste(cn,pn,"BEFORE",sep=" - "),legend=F)
      matplotf(t(xframeratio[years]),main=paste(cn,pn,"BEFORE",sep=" - "),legend=F)
      
      xframeratio[,years] <- t(apply(xframeratio[,years],1,function(x){sm(x,df=length(x)*0.5)}))
      xframe[,years] <- t(apply(xframeratio[,years],1,function(x){x * tot}))
      
      matplotf(t(xframeratio[years]),main=paste(cn,pn,"AFTER",sep=" - "),legend=F)
      matplotf(t(xframe[years]),main=paste(cn,pn,"AFTER",sep=" - "),legend=F)
      
      rposout[rposout$CountryCode==x$CountryCode & rposout$Possession==x$ProductName,] <- xframe
      
    } else {
      
      xframe <- rposout[rposout$CountryCode==x$CountryCode & rposout$Possession==x$ProductName & rposout$City %in% c("Total",x$City),]
      
      xframeratio <- xframe
      tot <- as.numeric(xframe[xframe$City == "Total",years])
      xframeratio[,years] <- t(apply(xframe[,years],1,function(x){x/tot}))
      
      matplotf(t(xframe[years]),main=paste(cn,pn,"BEFORE",sep=" - "),legend=F,abline=c(0,100))
      #     matplotf(t(xframeratio[years]),main=paste(cn,pn,"BEFORE",sep=" - "),legend=F)
      
      xframeratio[,years] <- t(apply(xframeratio[,years],1,function(x){sm(x,df=length(x)*0.5)}))
      xframe[,years] <- t(apply(xframeratio[,years],1,function(x){x * tot}))
      
      #     matplotf(t(xframeratio[years]),main=paste(cn,pn,"AFTER",sep=" - "),legend=F)
      matplotf(t(xframe[years]),main=paste(cn,pn,"AFTER",sep=" - "),legend=F,abline=c(0,100))
      
      rposout[rposout$CountryCode==x$CountryCode & rposout$Possession==x$ProductName & rposout$City %in% c("Total",x$City),] <-  xframe
      
    }
    
  }
  
  return(rposout)
}


issue.stretch.f <- function(bt, rposout){
  
  issue.stretch <- bt[bt$Issue=="stretch",]
  
  for(index in 1:nrow(issue.stretch)){
    
    x <- issue.stretch[index,]
    cn <- x$CountryCode[1]
    pn <- x$ProductName[1]
    city <- x$City[1]
    arg <- x$arg
    
    citylist <- c()
    if(city=="all"){
      citylist <- unique(rposout[rposout$CountryCode==x$CountryCode,"City"])
    } else {
      citylist <- city
    }
    
    for(cityname in citylist){
      
      xframe <- rposout[rposout$CountryCode==x$CountryCode & rposout$Possession==x$ProductName & rposout$City==cityname ,]
      matplotf(t(xframe[years]),main=paste("BEFORE",cityname),abline=c(0,100))
      
      xframe[,years] <- t(apply(xframe[,years],1,function(x){stretch(as.numeric(x) + seq(0,0.1,length.out=length(as.numeric(x))),x[1],arg)}))
      matplotf(t(xframe[years]),main=paste("AFTER",cityname),abline=c(0,100)) 
      rposout[rposout$CountryCode==x$CountryCode & rposout$Possession==x$ProductName & rposout$City %in% unique(xframe$City),] <- xframe
    }
  }
}


issue.total.ratio.f <- function(bt,rposout){
  
  issue.total.ratio <- bt[bt$Issue=="total.ratio",]
  
  for(index in 1:nrow(issue.total.ratio)){
    
    x <- issue.total.ratio[index,]
    arg <- x$arg
    #   xframe <- rposout[rposout$CountryCode==x$CountryCode & rposout$Possession==x$ProductName & rposout$City==x$City ,]
    
    if(x$City == "all"){
      
      xframe <- rposout[rposout$CountryCode==x$CountryCode & rposout$Possession==x$ProductName ,]
      xframe <- xframe[!xframe$City %in% c(x$CountryCode,"Total"),]
      tot <- rposout[rposout$CountryCode==x$CountryCode & rposout$Possession==x$ProductName & rposout$City=="Total" ,]
      
      for(i in 1:dim(xframe)[1]){
        ratio <- xframe[i,as.character(arg)] / tot[,as.character(arg)]
        xframe[i,years] <- tot[,years] * ratio
        rposout[rposout$CountryCode==x$CountryCode & rposout$Possession==x$ProductName & rposout$City==xframe$City[i],] <- xframe[i,]
      }   
      
    } else {
      
      xframe <- rposout[rposout$CountryCode==x$CountryCode & rposout$Possession==x$ProductName & rposout$City==x$City ,]
      if(x$CountryCode=="IR" & x$City=="Isfahan") xframe <- xframe[1,]
      
      tot <- rposout[rposout$CountryCode==x$CountryCode & rposout$Possession==x$ProductName & rposout$City=="Total" ,]
      ratio <- xframe[,as.character(arg)] / tot[,as.character(arg)]
      
      matplotg(t(xframe[years]),main="BEFORE")
      
      xframe[,years] <- tot[,years] * ratio
      matplotg(t(xframe[years]),main="AFTER") 
      rposout[rposout$CountryCode==x$CountryCode & rposout$Possession==x$ProductName & rposout$CityCode==xframe$CityCode ,] <- xframe
      
    }
  }
  
  return(rposout)
  
}


graphs.Badthings.adjustion <- function(bt, rposout, s5.rposout, filename){
  
  pdftitle <- paste0(filename,gsub(":","-",Sys.time()),".pdf")
  
  pdf(pdftitle,width=16 * 0.75,height=13 * 0.5)
  print(titlepage("Bad Things","Possession after Adjustion",date=Sys.time(),size=15,author="Povilas Bockus"))
  
  # changed graphs
  # ucp <- unique(bt[,c(1,3,4)])
  # ucp <- unique(bt[,c(1,2,3,4)])
  ucp <- bt
  ucp <- ucp[ucp$Issue!="dont.make.fakes",]
  
  for(index in 1:nrow(ucp)){
    print(index)
    
    x <- ucp[index,]
    if(x$City=="all"){
      xnew <- rposout[rposout$CountryCode==x$CountryCode[1] & rposout$Possession==x$ProductName[1],]
      xold <- s5.rposout[s5.rposout$CountryCode==x$CountryCode[1] & s5.rposout$Possession==x$ProductName[1],]
    } else {
      xnew <- rposout[rposout$CountryCode==x$CountryCode[1] & rposout$Possession==x$ProductName[1] & rposout$City %in% c("Total",x$City[1]),]
      xold <- s5.rposout[s5.rposout$CountryCode==x$CountryCode[1] & s5.rposout$Possession==x$ProductName[1] & rposout$City %in% c("Total",x$City[1]),]
    }
    
    xnew2 <- data.frame(t(xnew[,years]))
    xold2 <- data.frame(t(xold[,years]))
    
    names(xnew2) <- xnew$City
    names(xold2) <- xold$City
    
    matplotf(xold2,text=T,legend=F,abline=c(0,100),stx=2005,
             main=paste(x$CountryCode[1],x$ProductName[1],x$Issue[1],x$City[1],"BEFORE",sep=" - "))
    lines(2005:2030,xold2$Total,lwd=3)
    
    matplotf(xnew2,text=T,legend=F,abline=c(0,100),stx=2005,
             main=paste(x$CountryCode[1],x$ProductName[1],x$Issue[1],x$City[1],"AFTER",sep=" - "))
    lines(2005:2030,xnew2$Total,lwd=3)
    
  }
  
  dev.off()
  
}


make.fakes <- function(hpregion, hp, dont.make, possession_codes){
  
  for(index in 1:length(hpregion)){
    
    rc <- as.character(hpregion[index])
    
    print(paste(rc,sep=" - "))
    
    his <- hp[hp$RegionCode == rc,] #paziurim, kokius turim tam miestui/regionui istorinius
    cc <- his$CountryCode[1]
    
    dont.make <- dont.make.fakes[dont.make.fakes$CountryCode==cc,]
    
    # pc int br klasteris (2506, 59927, 87218)
    
    # clustID <- c(2506,59927,87218)
    clustID <- c(59927,87218)
    
    if(any(clustID %in% unique(his$ProductID)) & 
       !all(clustID %in% unique(his$ProductID)) & !("Internet" %in% dont.make$ProductName)){
      
      print("SOMETHING MISSING IN original: `pc int br`")
      
      missID <- clustID[which(!clustID %in% unique(his$ProductID))]
      notmissID <- clustID[which(clustID %in% unique(his$ProductID))]
      
      #klausimas, kodėl tik iki 2020 metų
      fakeout <- data.frame(matrix(NA,nrow=length(missID),ncol=length(years2) + 1))
      fakeout[,1] <- missID
      names(fakeout) <- c("ProductID",years2)
      
      ### Čia skirtingi variantai, žiūrint kiek duomenų trūksta, prikuriama vis skirtingai
      
      if(length(missID)==2){   
        for(id in missID){
          ratio <- rposout[(rposout$ProductID %in% notmissID) & (rposout$CityCode %in% rc),years2]/
            rposout[(rposout$ProductID %in% id) & (rposout$CityCode %in% rc),years2]
          fakeout[fakeout$ProductID==id,years2] <- as.numeric(hp[hp$ProductID==notmissID & hp$RegionCode==rc,years2]/ratio)  
        }
      }
      
      ### jei trūksta tik PC (bet yra šalies PC ir interneto total tarp istorinių duoemenų)
      if(length(missID)==1 & missID==2506){
        if(dim(hp[hp$CountryCode == cc & hp$CountryName=="Total" & hp$ProductID==2506,])[1]!=0 & 
           dim(hp[hp$CountryCode == cc & hp$CountryName=="Total" & hp$ProductID==59927,])[1]!=0){
          #pasiima šalies santykį
          ratio <- hp[hp$CountryCode == cc & hp$CountryName=="Total" & hp$ProductID==59927,years2]/
            hp[hp$CountryCode == cc & hp$CountryName=="Total" & hp$ProductID==2506,years2]
          if(sum(!is.na(ratio)) == 1) ratio <- rep(ratio[!is.na(ratio)],length(ratio))
          ratio <- MASplineVector(ratio) #is prognozuoja ratio su moving average, 
          # jei jo truksta kuriems nors metams
          
          fakeout[fakeout$ProductID==missID,years2] <- as.numeric(his[his$ProductID==59927,years2])/ratio  
          
        } else {
          
          ratio <- rposout[rposout$ProductID==59927 & rposout$CityCode==rc,years2]/
            rposout[rposout$ProductID==missID & rposout$CityCode==rc,years2]
          fakeout[fakeout$ProductID==missID,years2] <- as.numeric(hp[hp$ProductID==59927 & hp$RegionCode==rc,years2]/ratio)  
        }
        
      }
      
      if(length(missID)==1 & missID==59927){
        
        if(dim(hp[hp$CountryCode == cc & hp$CountryName=="Total" & hp$ProductID==2506,])[1]!=0 & 
           dim(hp[hp$CountryCode == cc & hp$CountryName=="Total" & hp$ProductID==59927,])[1]!=0){
          
          ratio <- hp[hp$CountryCode == cc & hp$CountryName=="Total" & hp$ProductID==2506,years2]/
            hp[hp$CountryCode == cc & hp$CountryName=="Total" & hp$ProductID==59927,years2]
          if(sum(!is.na(ratio)) == 1) ratio <- rep(ratio[!is.na(ratio)],length(ratio))
          ratio <- MASplineVector(ratio)
          
          fakeout[fakeout$ProductID==missID,years2] <- as.numeric(his[his$ProductID==2506,years2])/ratio  
          
        } else {
          
          ratio <- rposout[rposout$ProductID==87218 & rposout$CityCode==rc,years2]/
            rposout[rposout$ProductID==missID & rposout$CityCode==rc,years2]
          fakeout[fakeout$ProductID==missID,years2] <- as.numeric(hp[hp$ProductID==87218 & hp$RegionCode==rc,years2]/ratio)  
        }
        
      }
      
      if(length(missID)==1 & missID==87218){
        
        if(dim(hp[hp$CountryCode == cc & hp$CountryName=="Total" & hp$ProductID==87218,])[1]!=0 & 
           dim(hp[hp$CountryCode == cc & hp$CountryName=="Total" & hp$ProductID==59927,])[1]!=0){
          
          ratio <- hp[hp$CountryCode == cc & hp$CountryName=="Total" & hp$ProductID==59927,years2]/
            hp[hp$CountryCode == cc & hp$CountryName=="Total" & hp$ProductID==87218,years2]
          if(sum(!is.na(ratio)) == 1) ratio <- rep(ratio[!is.na(ratio)],length(ratio))
          ratio <- MASplineVector(ratio)
          
          fakeout[fakeout$ProductID==missID,years2] <- as.numeric(his[his$ProductID==59927,years2])/ratio  
          
        } else {
          ratio <- rposout[rposout$ProductID==59927 & rposout$CityCode==rc,years2]/
            rposout[rposout$ProductID==missID & rposout$CityCode==rc,years2]
          fakeout[fakeout$ProductID==missID,years2] <- as.numeric(hp[hp$ProductID==59927 & hp$RegionCode==rc,years2]/ratio)  
          
        }
        
        
      }
      
      
      
      out <- his[rep(1,dim(fakeout)[1]),]
      out[,names(fakeout)] <- fakeout
      out$ProductName <- possession_codes[(match(fakeout$ProductID,possession_codes$ProductID)),"ProductName"]
      out[,as.character(2000:2004)] <- NA
      
      hp <- rbind(hp,out)
    }
    
    ### čia pabaiga tų skirtingų variantų, kurie susiję su PC
    
    ### tęsiam toliau kitų fake kūrimą - ctv stv cbl klasteris (2499, 12987, 12980)
    
    clustID2 <- c(2499,12987,12980)
    
    if(any(clustID2 %in% unique(his$ProductID)) & 
       !all(clustID2 %in% unique(his$ProductID)) &
       !("TV" %in% dont.make$ProductName)){
      
      print("SOMETHING MISSING IN: `ctv stv cbl`")
      
      missID <- clustID2[which(!clustID2 %in% unique(his$ProductID))]
      notmissID <- clustID2[which(clustID2 %in% unique(his$ProductID))]
      
      fakeout <- data.frame(matrix(NA,nrow=length(missID),ncol=length(years2) + 1))
      fakeout[,1] <- missID
      names(fakeout) <- c("ProductID",years2)
      
      ###
      
      if(length(missID)==2){   
        for(id in missID){
          if(dim(rposout[rposout$ProductID==id & rposout$CityCode==rc,years2])[1]==0) next
          ratio <- rposout[rposout$ProductID==notmissID & rposout$CityCode==rc,years2]/
            rposout[rposout$ProductID==id & rposout$CityCode==rc,years2]
          fakeout[fakeout$ProductID==id,years2] <- as.numeric(hp[hp$ProductID==notmissID & hp$RegionCode==rc,years2]/ratio)  
        }
      }
      
      if(length(missID)==1 & missID==2499){
        ratio <- rposout[rposout$ProductID==12980 & rposout$CityCode==rc,years2]/
          rposout[rposout$ProductID==missID & rposout$CityCode==rc,years2]
        fakeout[fakeout$ProductID==missID,years2] <- as.numeric(hp[hp$ProductID==12980 & hp$RegionCode==rc,years2]/ratio)  
      }
      
      if(length(missID)==1 & missID==12987){
        if(dim(rposout[rposout$ProductID==missID & rposout$CityCode==rc,years2])[1]==0) next
        ratio <- rposout[rposout$ProductID==12980 & rposout$CityCode==rc,years2]/
          rposout[rposout$ProductID==missID & rposout$CityCode==rc,years2]
        fakeout[fakeout$ProductID==missID,years2] <- as.numeric(hp[hp$ProductID==12980 & hp$RegionCode==rc,years2]/ratio)  
      }
      
      if(length(missID)==1 & missID==12980){
        if(dim(rposout[rposout$ProductID==missID & rposout$CityCode==rc,years2])[1]==0) next
        ratio <- rposout[rposout$ProductID==12987 & rposout$CityCode==rc,years2]/
          rposout[rposout$ProductID==missID & rposout$CityCode==rc,years2]
        fakeout[fakeout$ProductID==missID,years2] <- as.numeric(hp[hp$ProductID==12987 & hp$RegionCode==rc,years2]/ratio)  
      }
      
      ###
      
      out <- his[rep(1,dim(fakeout)[1]),]
      out[,names(fakeout)] <- fakeout
      out$ProductName <- possession_codes[(match(fakeout$ProductID,possession_codes$ProductID)),"ProductName"]
      out[,as.character(2000:2004)] <- NA
      
      hp <- rbind(hp,out)
    } 
    
  }
  
  return(hp)
}



