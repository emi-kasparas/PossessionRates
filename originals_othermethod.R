
# issue...originals.other.method ------------------------------------------
  
  for(index in 1:nrow(bt.other.method)){
    
    x <- bt.other.method[index,]
    cc <- x$CountryCode[1]
    pn <- x$ProductName[1]
    city <- x$City[1]
    
    print(paste(cc,city,pn,sep=" - "))
    
    his <- hp[hp$CountryCode == cc & hp$ProductName== pn,]
    his <- his[!duplicated(his[,1]),]
    his <- his[!apply(his[,years2],1,function(x){all(is.na(x))}),] 
    
    citylist <- c()
    if(city=="all"){
      citylist <- rposout[rposout$CountryCode==cc & rposout$Possession==pn,"City"]
    } else {
      citylist <- city
    }
    
    xframe <- s8rposout[s8rposout$CountryCode==cc & s8rposout$Possession==pn,]
    if(dim(xframe)[1]==0) next
    if(dim(his)[1]==0) next
    
    oxframe <- xframe
    
    # priforecastinamas totalas
    totalratio <- his[his$CountryName=="Total",years2]/xframe[xframe$City=="Total",years2]
    if(sum(!is.na(totalratio))==1) totalratio <- rep(totalratio[!is.na(totalratio)],length(totalratio))
    totalratio <- MASplineVector(totalratio)
    his[his$CountryName=="Total",years2] <- xframe[xframe$City=="Total",years2] * totalratio
    
    # suskaiciuojami santykiai ir jie sukalami ant rposout
    levelmat <- his
    
    
    # #multi
    # totaloriginal <- as.numeric(his[his$CountryName=="Total",years2])
    # totaldatabase <- as.numeric(xframe[xframe$City=="Total",years2])
    # levelmat[,years2] <- t(apply(levelmat[,years2],1,function(x){as.numeric(x)/totaloriginal * totaldatabase}))
    
    #adi
    totaloriginal <- as.numeric(his[his$CountryName=="Total",years2])
    totaldatabase <- as.numeric(xframe[xframe$City=="Total",years2])
    levelmat[,years2] <- t(apply(levelmat[,years2],1,function(x){as.numeric(x)-totaloriginal + totaldatabase}))
    
    levelmat[,years2][levelmat[,years2] < 0 & !is.na(levelmat[,years2] < 0)] <- 
      min(levelmat[,years2][levelmat[,years2] > 0 & !is.na(levelmat[,years2] < 0)])
    
    #   levelhp[levelhp$CountryCode == cc & levelhp$ProductName== pn,years2] <- levelmat[,years2]
    
    
    # ███████████████████████ pridedami fiktyvus originalus
    
    miss <- intersect(setdiff(xframe$CityCode,levelmat$RegionCode),tier$RegionCode)
    
    out <- c()
    #missing regions adjustion
    if(length(miss)!=0){ 
      for(mr in miss){
        print(paste("Missing Region:",mr))
        
        mc <- tier[tier$RegionCode == mr,"CityCode"]
        out <- levelmat[levelmat$CountryName=="Total",]
        out$RegionCode <- mr
        out[,c("Region/CityName")] <- NA
        out[,c("CountryName")] <- "Temp"
        
        if(any(mc %in% unique(levelmat$RegionCode))){
          out[,years2] <- as.numeric(apply(levelmat[levelmat$RegionCode %in% c(mc,cc),years2],2,mean))
          levelmat <- rbind(levelmat,out)
          
        } else next
        
        print(out)
        
      }
    }
    
    # pridedami fiktyvus originalus miestam, jei turim regionus
    
    setdiff(xframe$CityCode,levelmat$RegionCode)
    
    regs <- intersect(intersect(xframe$CityCode,tier$RegionCode),levelmat$RegionCode)
    
    for(rr in regs){
      
      cinr <- tier[tier$RegionCode == rr,"CityCode"]
      cinr <- intersect(cinr[!cinr %in% unique(levelmat$RegionCode)],xframe$CityCode)
      
      if(length(cinr)!= 0){
        for(i in 1:length(cinr)){
          print(paste("City changed because of Region: ", cinr[i],"because of",rr))
          x <- oxframe[oxframe$CityCode==cinr[i],years]
          xreg <- oxframe[oxframe$CityCode==rr,years]
          ratio <- x[,years]/xreg
          #         x <-  xframe[xframe$CityCode==rr,years] * ratio
          #         xframe[xframe$CityCode==cinr[i],years]  <- x
          out <- levelmat[levelmat$RegionCode==rr,]
          out$RegionCode <- cinr[i]
          out[,years2] <- out[,years2] * ratio[,years2]
          levelmat <- rbind(levelmat,out)
        }
      }
    }
    
    levelmat[,years2][levelmat[,years2]>100] <- 100
    
    
    
    p <- xframe[xframe$City=="Total",years]/100
    #███
    citymustbetchanged <- intersect(unique(levelmat$RegionCode)[unique(levelmat$RegionCode)!=cc],xframe$CityCode)
    citylistcodes <- xframe[xframe$City %in% citylist,"CityCode"]
    citymustbetchanged <- intersect(citymustbetchanged,citylistcodes)
    
    if(length(citymustbetchanged)==0) print("ALARM ALARM THERE IS NO CITY WITH THAT NAME")
    
    for(mcc in citymustbetchanged){
      
      x <- xframe[xframe$CityCode==mcc,]
      if(dim(x)[1]==0) next
      y <- levelmat[levelmat$RegionCode==mcc,]
      nnay <- years2[which(!is.na(y[,years2]))]
      ymean <- mean(as.numeric(y[,nnay]))
      xmean <- mean(as.numeric(x[,nnay]))
      
      
      #     if(dim(bt.other.method[bt.other.method$CountryCode==cc & bt.other.method$ProductName==pn,])[1]>0){
      ### 
      
      cat(paste("Other originals adjustion method used:",mcc),"\n")  
      
      nnay <- floor((as.numeric(nnay[1]) + as.numeric(nnay[length(nnay)]))/2)
      kofori <- kof[kof$CityCode==mcc,]
      kofcn <- kof[kof$CityCode==cc,]
      
      koforiini <- as.numeric(kofori[kofori$Year==nnay,c("a","b.cntUSD","q")])
      pmat <- pd[pd$CountryCode==cc & pd$Possession==pn,] 
      pframe <- pmat[pmat$Year==nnay,c("est.Pos.norm","Average.income")] 
      
      b.pre <- koforiini[2]
      koforiini[2] <- posfun.reverse(ymean)
      
      if(koforiini[2]<0) koforiini[2] <- min(pframe[,2])
      
      
      if(any(is.na(koforiini))) next
      
      minifun <- function(b){
        #       print(b)
        if(b<0) out <- 10^200 else {
          theta <- data.frame(matrix(c(koforiini[1],b,koforiini[3]),nrow=1))
          names(theta) <- c("a","b.cntUSD","q")
          out <- (regionalpos(theta,pframe,plot=F,starter = tail(pframe$Average.income,1) * 20) - ymean)^2
        }
        return(out)
      }
      
      # aaaa <- Sys.time()
      fitb <- optim(koforiini[2],fn = minifun,control = list(maxit=10,trace=0,REPORT=2),
                    method="L-BFGS-B")$par
      # Sys.time() - aaaa   
      
      ratio <- b.pre/fitb
      
      koforinew <- kofori
      koforinew$b.cntUSD <- koforinew$b.cntUSD/ratio
      
      #   print(paste0(mcc,"   ",fitb))
      out <- c()
      for(y in 2005:2030){
        #       print(y)
        k <- koforinew[koforinew$Year==y,]
        pframe <- pmat[pmat$Year==y,]
        out[y-2004] <- regionalpos(k,pframe,plot=F,starter = tail(pframe$Average.income,1) * 20)
      }
      
      x.out <- out
      
      tot <- as.numeric(xframe[xframe$CityCode==cc,years])
      
      if(all(diff(tot)<=0) | all(diff(tot)>=0)){
        x.out <- mono.pos(x.out,tot,split = 1)
      }
      
      x[,years] <- x.out
      
      xframe[xframe$CityCode==mcc,] <- x
      
    }
    
    # pramusimas over100.v2
    xframe[,years] <- over100.v2(xframe[,years],plot=F,divover=2,border = 1)
    #   xframe[,years] <- over100.v2(xframe[,years],plot=F,divover=2,border = 1)
    
    #   if(dim(xframe)[1]==0) next
    
    ####################################
    
    #   xframe[,years] <- over100.v2(xframe[,years],plot=F,divover=2,border = 1)
    
    #   notchanged <- setdiff(xframe$CityCode[which(apply(oxframe[,years] - xframe[,years],1,mean) == 0)],c(cc,"Total"))
    #   if(length(notchanged)!=0) print(paste("There left few 'CityCodes' which was not changed or was equal (strange):",notchanged))
    
    #███
    rposout[rposout$CountryCode==cc & rposout$Possession==pn & rposout$City%in% citylist,] <- xframe[xframe$City%in% citylist,]
    #███
    
    
  }
  

# manipulation ------------------------------------------------------------

bl <- c()

manipulate.fun <- function(rposout){
  
  require(manipulate)
  
  manipulate({
    # print("Check")
    xframe <- rposout[rposout$CountryCode==cc & rposout$Possession==pn,]
    
    if(dim(xframe)[1]==0) print("There is NO such Country x Possession") else {
      
      xframe.print <- data.frame(t(xframe[,years]))
      names(xframe.print) <- xframe$City
      
      matplotf(xframe.print,main=paste(cc,pn,sep=" - "),legend=F,point=p,abline=c(minabline,maxabline),text=text)
      if(total) lines(xframe.print[,"Total"],lwd=2)
      
      if(add) {print(paste0("Added to Bad List: ",cc," - ",pn))
        bl <<- rbind(bl,c(cc,pn))
      }
    }
  },
  
  cc = picker(as.list(unique(rposout$CountryCode)),label = "Country Code"),
  pn = picker(as.list(unique(rposout$Possession)),label = "Possession"),
  p = checkbox(F, "Points"),
  text = checkbox(F, "Text"),
  total = checkbox(F, "Total Line"),
  minabline = slider(min = 0,max=100,initial = 0,step = 1),
  maxabline = slider(min = 0,max=100,initial = 100,step = 1),
  add = button("Add to Bad List")
  )
  
}



manipulate(
plot(cars, xlim=c(x.min,x.max),axes=axes.menu), 
  x.min=slider(0,15), 
  x.max=slider(15,30),
axes.menu = checkbox(F,"Axes"))



