source("mid108(2).R")
source("functions.R")
library("grid")

#emibase()
load("input/EMI base.RData")


options(scipen = 7)

# INPUT -------------------------------------------------------------------
# inputas pasiima pacias naujausias versijas

# is possessions by decile csv file'o, padarome Rdata file'a,
# kuriame idedam naujus possessions by decile (is csv) ir Average income by decile is senesnio Rdata file'o.
# dirname1 <- "input/Possessions by Decile/final_poss_by_decile csv/"
# dirname2 <- "input/Possessions by Decile/final_poss_by_decile RData/"
# create.possessions.by.decile.RData(dirname1, dirname2)

load(lastinput("input/Possessions by Decile/final_poss_by_decile RData"))
pd <- arrange(final.pbd,Decile.No,Year,CountryName) # possessions by decile
kof <- read.csv(file=lastinput("input/Income abq"),check.names=F,stringsAsFactors=F)
possession_codes <- read.csv("input/possession_codes.csv", check.names=F, stringsAsFactors = F)
bt <- read.csv(file="input/bad things.csv",check.names=F,stringsAsFactors=F, sep=";")

# kodėl reikalingas šitas kof koregavimas? imi 1,3,4 stulpelį ir ten pridedi kažkokius minusus???
kof <- kof[which(!(duplicated(apply(kof[,c(1,3,4)],1,function(x){paste(x,collapse="-")})))),]
abqUSD <- kof

ccid <- read.csv(file="input/city codes id.csv",check.names=F,stringsAsFactors=F)

years <- as.character(2005:2030)

# cia ispiesiamas summary abq (bet ji galim ir uzkomentuot, jei ten viskas ok)
abq.summary.graphs.grubus(kof,filename="plots/1. Average income/abq summary ",info1="abq Summary",info2="Average, Median, Mode")


# mapping -----------------------------------------------------------------

splitmap <- read.csv(file="input/countries mapping splits.csv",check.names=F,stringsAsFactors=F)

splitmap <- ddply(splitmap,~ GraphGroup + GraphName,function(x){
  x<<-x
  out <- c(unique(x$CityCode),substr(unique(x$RegionCode),1,4),unique(x$CountryCode))
  out <- c(unique(x$CityCode),unique(x$CountryCode))
  out <- data.frame(cbind(CityCode = out,GraphGroup = x$GraphGroup[1],GraphName = x$GraphName[1]))
  return(out)
})

osplitmap <- splitmap

# PROCESS -----------------------------------------------------------------
# recalculating average decile income from countries abq ------------------

cnkof <- kof[kof$CityCode==kof$CountryCode,]
cnlist <- unique(cnkof$CountryCode)

# suskaičiuoja naujas vidutines pajamas miesto deciliuose, kad po to 
# būtų galima suskaičiuoti possessions
for(cn in unique(cnkof$CountryCode)){
  cat(paste0("Average deciles recalculation: ",cn),"\n")
  
  xframe <- cnkof[cnkof$CountryCode==cn,] # selecting one country
  
  for(y in 2005:2030){
    theta <- as.numeric(xframe[xframe$Year==y,c("a","b.cntUSD","q")])
    avgdec <- avgdec.from.qbq(theta)
    
    pd[pd$CountryCode==cn & pd$Year==y,"Average.income"] <- rep(avgdec,each=length(unique(pd[pd$CountryCode==cn & pd$Year==y,"Possession"])))     
  }
}

# average cities income recalculation -------------------------------------
# this code selects cities from kof, takes its a, b, q, 
# and gets city's average income for each year

avicity <- kof

avicity <- ddply(avicity,~CityCode,function(x){
  x<<-x
  print(as.character(x$CityCode[1]))
  ttt <- b.converter.reverse(x[,"a"],x[,"b.cntUSD"],x[,"q"])
  out <- data.frame(cbind(x[,1:4],AI=ttt),check.names=F,stringsAsFactors=F)
  return(out)
})

#  ------------------------------------------------------------------------

years <- as.character(2005:2030)
years2 <- as.character(2005:2020)

###################

# skaiciuojame viskam, kam turime income distributionus
ukofcncc <- unique(kof[,c("CountryCode","CityCode","City")]) 

upcn <- unique(pd[,c("CountryCode","Possession")])

# Pasalinami saliu duomenys (a,b,q), nes gerus pasiimam is Astos
# ukofcncc <- ukofcncc[ukofcncc$CityCode!=ukofcncc$CountryCode,]


# # sita reiks istrint kai noresim visus produktus, bet jei reik tik keliu, tai cia juos issiskiriam

# upcn <- upcn[upcn$Possession %in% c("Possession of Mobile Telephone"),]

# upcn <- upcn[upcn$Possession %in% c("Possession of Personal Computer","Possession of Internet Enabled Computer",
#                                     "Possession of Broadband Internet Enabled Computer"),]
# #


## pridedam possessions kiekvienam regionui/miestui
ukofcncc <- merge(ukofcncc,upcn,by="CountryCode",all.x=T)
## čia išmetam NA eilutes (gal jau ir nebereikia šito)
ukofcncc <- ukofcncc[!is.na(ukofcncc$Possession),]
## išmetam eilutes, kur šalies kodas = miesto kodas, t.y., total eilutes
ukofcncc <- ukofcncc[ukofcncc$CountryCode!=ukofcncc$CityCode,]

# sukuriamas sablonas outputui sukalti
rposout <- cbind(ukofcncc,ukofcncc[,rep(1,length(years))])
names(rposout)[-1:-4] <- years
rposout[years] <- NA

# pridedamos TOTAL eilutes, kuriose bus salies possesionai, kad grafikus normaliai piest.. Galima sita ir ismest
rposout <- ddply(rposout,~ CountryCode + Possession,function(x){
  x <<- x
  # })
  cat(paste0("Calculating Country Totals: ",x$CountryCode[1]),"\n")
  out <- x[1,]
  out[,c("CityCode","City")] <- rep("Total",2)
  out[years] <- ddply(pd[pd$CountryCode==out$CountryCode & pd$Possession==out$Possession,c("Year","est.Pos.norm")],
                      ~Year,function(y){mean(y[,"est.Pos.norm"])})$V1
  out$CityCode <- out$CountryCode
  return(rbind(x,out))
})


# FORECASTING -------------------------------------------------------------
# possessions forecasting using income distributions ----------------------

ukofcnccl1 <- as.character(unique(ukofcncc$CountryCode)) # tiesiog visos šalys
pd <- arrange(pd,Decile.No)

s3rposout <- rposout # saving some rations just for recovering if needed

for(indexl1 in 1:length(ukofcnccl1)){ # darome ciklą tiek kartų, kiek yra šalių
  
  cn <- as.character(ukofcnccl1[indexl1]) # imam šalį
  
  ukofcnccl2 <- ukofcncc[ukofcncc$CountryCode==cn,] # imam jos regionus ir possessionus
  kofl1 <- kof[kof$CountryCode==cn,] # šalies income distribution parametrai
  pmatl1 <- pd[pd$CountryCode==cn,] # šalies decilių possession duomenys
  
  # dabar jau sukam ciklą tos šalies viduje
  for(index in 1:nrow(ukofcnccl2)){  
    
    cc <- as.character(ukofcnccl2[index,"CityCode"]) # konkretus miestas
    pn <- as.character(ukofcnccl2[index,"Possession"]) # konkretus possession
    
    # pasižiūrim šalies total duomenis (kuriuos įsidėjom jau anksčiau)
    tot <- as.numeric(rposout[rposout$CityCode==cn & rposout$Possession==pn,years])
    
    cat(paste("Forecasting: ",cn,cc,ukofcnccl2[index,"City"],pn,sep=" - "),"\n")
    
    out <- c()
    kmat <- kofl1[kofl1$CityCode==cc,] # miesto ar regiono income distribution parametrai visiems metams
    pmat <- pmatl1[pmatl1$Possession==pn,] # šalies decilių possession visiems metams, bet konkretaus possession
    
    # gaunam koks miesto income 2005-2030, iš jo distribution
    cityincome <- b.converter.reverse(kmat$a,kmat$b.cntUSD,kmat$q)
    avicity[avicity$CityCode %in% c(cn),"AI"]
    
  
    for(y in 2005:2030){
      #       print(y)
      k <- kmat[kmat$Year==y,] # vienerių metų income distribution parametrai
      pframe <- pmat[pmat$Year==y,] # šalies possession pagal decilius, vieniems metams)
      
      #integruojam ir gaunam possession, kai ciklas prasisuka, tai visiems metams
      out[y-2004] <- regionalpos(k,pframe,plot=T,starter = tail(pframe$Average.income,1) * 20)
    }
    
    oout <- out
    
    # šita funkcija patikrina monotoniškumą - kad possession visąlaik augtų arba visąlaik kristų
    # iš principo patvarko duomenis labai nedaug
    out <- hatfun.solve(out,itermax=10,stp=9,plot=T)
    
    # monotoniskumas - čia tiems atvejais, kai total (t.y.šalis) visąlaik auga/krenta - 
    # tada griežčiau padaro, kad ir miestas visąlaik mažėtų/augtų.
    
    if (substr(cc, 1, 2)=="UZ"){ #monotoniskumas tik Uzbekistano atveju uzdetas
      if(all(diff(tot)<=0) | all(diff(tot)>=0)){
        out <- mono.pos(out,tot,split = 1)
      }
    }

    # matplotg(cbind(oout,out,tot),main=paste(cn,cc,ukofcnccl2[index,"City"],pn,sep=" - "))
    
    rposout[rposout$CountryCode==cn & rposout$CityCode==cc & rposout$Possession==pn,years] <- out
    
  }
}


# possessions adjustion using historical data -----------------------------

tier <- read.csv(file="input/Region city mapping file.csv",check.names=F,stringsAsFactors=F)
tier[,"RegionCode"] <- substr(tier[,"RegionCode"],1,4)
regnew <- read.csv(file="input/Region codes new.csv",check.names=F,stringsAsFactors=F)
tier <- merge(tier,regnew[c("RegionCode","RegionCodeOld")],all.x=T,by="RegionCode")


cities.possession.graphs.grubus.only.after.income(rposout,filename="plots/2. Forecasted possessions only from income/Forecasted possessions only from income ",
                                                  info1="Possessions of Durables",info2="Forecasts only after income method step")

# padaromi pritempimai prie originaliu duomenu ----------------------------

s2rposout <- rposout
# rposout <- s2rposout
qwe <- rposout

save(s2rposout,file=paste0("temp/","s2rposout only after income",gsub(":","-",Sys.time()),".Rdata"))
# 
# load("temp/s2rposout only after income2016-06-06 07-19-18.RData") # sita reikia uzsiloadint, jei norim daryt nuo tos stadijos, kai jau turim susimodeliave pagal income
# rposout <- s2rposout
# 

# input
# hp yra historical possessions

hp <- read.csv(file=lastinput("input/Original data"),check.names=F,stringsAsFactors=F)

# išmetam šalis, kurių galėjo būti originaliuose duomenyse, bet mums jų nereikia
hp <- hp[hp$CountryCode %in% intersect(unique(hp$CountryCode),unique(rposout$CountryCode)),]

# jei yra 2004 metų taškas, priskiriame jį 2005 metams
hp[apply(hp[,years2],1,function(x){all(is.na(x))}) & !is.na(hp[,"2004"]),"2005"] <- 
  hp[apply(hp[,years2],1,function(x){all(is.na(x))}) & !is.na(hp[,"2004"]),"2004"]

# kurie originalus islipa is veziu 
hp[unique(which(hp[,years2]>100,arr.ind = T)[,1]),]
hp[unique(which(hp[,years2]<0,arr.ind = T)[,1]),]
hp[duplicated(hp[,c("RegionCode","ProductName")]),]

## 
hp <- hp[!duplicated(hp[,c("RegionCode","ProductName")]),]

###

# isbreziami grafikai (sita irgi galim komentuot jei nenorim kad breztu)
# original.data.barplots(hp,avicity,filename = "plots/3. Original data barplots/Original data barplots ",info1 = "Original Data Barplots",info2 = "Blue - Total, Green - Region, Yellow - City")

##
# idsnames <- unique(hp[,c("ProductID","ProductName")])
pnii <- unique(pd[,c("Possession","ProductID")])

# hp <- hp[hp$CountryCode %in% c("CO"),] # 

#countries and poss for these not making fakes
dont.make.fakes <- bt[bt$Issue=="dont.make.fakes",]

### Kodėl ne visi regionai yra tarp sumodeliuotų?? nes
### po antro veiksmo sumažėja duomenų
hpregion <- unique(hp[,c("RegionCode")])
hpregion <- intersect(hpregion,rposout$CityCode)

# hpregion <- hpregion[!hpregion$RegionCode %in% unique(hp$CountryCode),]

# prikabinimas product ID
rposout <- merge(rposout,possession_codes,by.x="Possession",by.y="ProductName",all.x=T)
rposout <- arrange(rposout,CityCode,Possession,CountryCode)
rposout[is.na(rposout$ProductID) & rposout$Possession=="Possession of Black and White TV Set","ProductID"]  <-  2496

#isorderinama, kad is pradziu prikabintu totalus, jei nera, o tik tada miestam uzvarytu
hpregion <- c(hpregion[hpregion %in% rposout$CountryCode],hpregion[!hpregion %in% rposout$CountryCode])

# hpregion <- c("ID", paste0("ID0", 1:9), paste0("ID", 10:16), paste0("R", 330:362))

###########################################################################
# FIKTYVIU ORIGINALIU DUOMENU SUKURIMAS -----------------------------------
###########################################################################

# sukuria fiktyvius duomenis 1) internet ir broadband internet; 2) bei cable TV, satellite, colour TV 
# atvejais, kai yra bent vienas is grupes sukuriami fiktyvus duomenys, 
# pvz. jeigu miestui turim colour TV, tai sukriam sattellite ir cable TV.
 
hp <- make.fakes(hpregion, hp, dont.make, possession_codes)

original.data.barplots(hp,avicity,filename = "plots/3. Original data barplots/Original data barplots after fakes ",info1 = "Original Data Barplots",info2 = "Blue - Total, Green - Region, Yellow - City")

###########################################################################
# 'FIKTYVIŲ DUOMENŲ' PRITEMPIMAS PRIE ORIGINALAUS TAŠKO -------------------
###########################################################################

### Veiksmų eiliškumas:
# 1. Ima kiekvieną regioną atskirai (pirmas for ciklas)
# 2. 'priforecastinamas totalas' - pakoreguojame pagal C&C research total ir Cities research total skirtumą
# 3. 'suskaiciuojami santykiai ir jie sukalami ant rposout' - visi cities research duomenys pakoreguojami pagal C&C duomenis ADITYVIAI
# 4. Regionai, kuriems turimi tik miestai, užpildomi pagal to regiono miestų vidurkį.
# 5. Miestai, kuriems turime tik regioną, suskaiciuojam orginal taska rementis orginaliu region tasku
#    multiplikatyviai (cia taikome du metodus ir issirenkam viena kuris duoda artimesni durable share 
#    regionui)
# 6. 'pramusimas over100.v2' - čia pataisome 100% pramušimą su atskira funkcija
# 7. Grafikai

# rposout - visada yra mūsų visi gaminami duomenys, kažką pakoreguojam, bet paliekam ta patį pavadinimą. 

s8rposout <- rposout

# pasiliekame tik tuos istorinius duomenis, kuriuos esame primodeliavę.
hpregion <- unique(hp[,c("CountryCode","ProductName")])
hpregion <- hpregion[hpregion$CountryCode %in% rposout$CountryCode & hpregion$ProductName %in% rposout$Possession,]
hpregion <- arrange(hpregion,CountryCode)


levelhp <- hp

filename <- "plots/4. Forecasted possessions adjusted using originals/Forecasted possessions adjusted using originals "
pdftitle <- paste0(filename,gsub(":","-",Sys.time()),".pdf")

pdf(pdftitle,width=13,height=10)
print(titlepage("Possessions Forecasts", "Before and After Original data adjustion",date=Sys.time(),size=15,author="Povilas Bockus"))

# which(hpregion[,1]=="AT" & hpregion[,2]=="Possession of Personal Computer")
# hpregion <- hpregion[hpregion$CountryCode=="IN",]
kkk <- 1
for(index in 1:nrow(hpregion)){
  
  cc <- as.character(hpregion[index,"CountryCode"])
  pn <- hpregion[index,"ProductName"]
  
  hp.cc <- hp[hp$CountryCode==cc & hp$ProductName==pn,]
  inc.cc <- rposout[rposout$CountryCode==cc & rposout$Possession==pn,]
  
  print(paste(cc,pn,sep=" - "))
  
  # his <- hp[hp$CountryCode == cc & hp$ProductName== pn,]
  his <- hp.cc
  his <- his[!duplicated(his[,1]),] ### ŠITIE DUPLICATED TRYNIMAI YRA BLOGAI??? juos reikia sutvarkyt source, bet ne imti ir trinti
  his <- his[!apply(his[,years2],1,function(x){all(is.na(x))}),] 
  
  xframe <- rposout[rposout$CountryCode==cc & rposout$Possession==pn,]
  if(dim(xframe)[1]==0) next
  if(dim(his)[1]==0) next
  
  oxframe <- xframe
  
  # priforecastinamas totalas (t.y. dirbam tik su ŠALIES TOTAL - C&C VS Cities)
  ## gaunam ratio šalis originali (iš 'his') vs šalis sumodeliuota (iš xframe)
  totalratio <- his[his$CountryName=="Total",years2]/xframe[xframe$City=="Total",years2]
  ## užpildom, jei istorinių duomenų yra tik vieni metai
  if(sum(!is.na(totalratio))==1) totalratio <- rep(totalratio[!is.na(totalratio)],length(totalratio))
  ## jei yra bent dvi istorinės reikšmės, jis su MASplineVector kažką padaro
  totalratio <- MASplineVector(totalratio)
  ## čia kažkodėl užrašome pakoreguotus ŠALIES istorinius (toliau visi miestai pakoreguojami eilutėse 500-515)
  his[his$CountryName=="Total",years2] <- xframe[xframe$City=="Total",years2] * totalratio
  
  # suskaiciuojami santykiai ir jie sukalami ant rposout
  levelmat <- his
  
  ## jamam jau pakoreguotus šalies total original
  totaloriginal <- as.numeric(his[his$CountryName=="Total",years2])
  ## čia jamam išmodeliuotus (kurie turbūt sutampa su C&C research)
  totaldatabase <- as.numeric(xframe[xframe$City=="Total",years2])
  ## Pridedam kiekvienam regionui/miestui skirtumą, koks buvo rastas tarp C&C research ir Cities Research šalies 
  ## duomenų. Kai kur gali gautis neigiami skaičiai ???
  levelmat[,years2] <- t(apply(levelmat[,years2],1,function(x){as.numeric(x)-totaloriginal + totaldatabase}))
  
  #   levelhp[levelhp$CountryCode == cc & levelhp$ProductName== pn,years2] <- levelmat[,years2]
    
  # ███████████████████████ pridedami fiktyvus originalus - 
  
  # regionai, kuriems NĖRA originalių duomenų 
  miss <- intersect(setdiff(xframe$CityCode,levelmat$RegionCode),tier$RegionCode)
  
  out <- c()
  
  # missing regions adjustion - (KAI TURĖJOME KAŽKOKIUS MIESTŲ DUOMENIS)
  if(length(miss)!=0){ 
    for(mr in miss){
      print(paste("Missing Region:",mr))
      
      mc <- tier[tier$RegionCode == mr,"CityCode"]
      out <- levelmat[levelmat$CountryName=="Total",]
      out$RegionCode <- mr
      out[,c("Region/CityName")] <- NA
      out[,c("CountryName")] <- "Temp"
      
      # paima visų to regiono miestų vidurkį, ir priskiria jį regionui
      if(any(mc %in% unique(levelmat$RegionCode))){
        out[,years2] <- as.numeric(apply(levelmat[levelmat$RegionCode %in% c(mc,cc),years2],2,mean))
        levelmat <- rbind(levelmat,out)
        
      } else next
      
      print(out)
      
    }
  }
  
  # pridedami fiktyvus originalus miestam, jei turim regionus
  
  # regionai, kuriems turim originalių duomenų
  regs <- intersect(intersect(xframe$CityCode,tier$RegionCode),levelmat$RegionCode)
  
  for(rr in regs){
    
    ## miestai, kuriame yra tame regione
    cinr <- tier[tier$RegionCode == rr,"CityCode"]
    ## dar išmetam tuos miestus, kuriems turim originalių duomenų
    cinr <- intersect(cinr[!cinr %in% unique(levelmat$RegionCode)],xframe$CityCode)
    
    if(length(cinr)!= 0){
      
      for(i in 1:length(cinr)){
        
        x <- oxframe[oxframe$CityCode==cinr[i],years] # miesto possession iš pajamų
        xreg <- oxframe[oxframe$CityCode==rr,years] # regiono possession iš pajamų
        ### Čia multiplikatyviai prideda miestus
        
        ### čIA IŠ PRADŽIŲ DVIEM SKIRTINGAIS BŪDAIS SPĖJAME MIESTO POSSESSION ATSIŽVELGDAMI Į REGION
        
        ## pirmas būdas primityvus multiplikatyvus
        ratio1 <- x[,years]/xreg
        ratio2 <- ratio1
        
        ## antras būdas labiau komplikuotas, išsamiau aprašytas functions.R faile
        ratio22 <- fromIncomePossessionCurve(city=cinr[i], cc=cc, pn=pn , rposout=rposout,
                                            hp=hp, abqUSD=abqUSD, final.pbd=final.pbd, hp.cc, inc.cc)
        ratio2[,as.character(c(2005:2030))] <- ratio22
        
        ### ŠITAS IF IŠRENKA TĄ BŪDĄ, KAI MIESTO SUMODELIUOTI DUOMENYS YRA ARTIMESNI REGIONUI
        ### (T.Y. TOKS LYG KONSERVATYVESNIS POŽIŪRIS GAUNASI)
        
        if (all(is.na(ratio2))){
          
          ratio <- ratio1
        }else{
          
          if (abs(1-mean(as.numeric(ratio1))) > abs(1-mean(as.numeric(ratio2)))){
            ratio <- ratio2
          }else{
            ratio <- ratio1
          }
        }

        out <- levelmat[levelmat$RegionCode==rr,] # priskiriam ou region duomenis
        # original.point <-  as.numeric( out[,years2])[!is.na(as.numeric( out[,years2]))][length(as.numeric( out[,years2])[!is.na(as.numeric( out[,years2]))])]
        # pos.inc <- xreg[,names(out[,years2][!is.na(as.numeric( out[,years2]))][length(as.numeric( out[,years2])[!is.na(as.numeric( out[,years2]))])])]
        # pos.n <- xframe[xframe$City=="Total",names(out[,years2][!is.na(as.numeric( out[,years2]))][length(as.numeric( out[,years2])[!is.na(as.numeric( out[,years2]))])])]
        out$RegionCode <- cinr[i] # pakeiciam region code, kad out jau butu miesto duomenys
        out[,years2] <- out[,years2] * ratio[,years2] # out padauginam is region/city ratio (is pajamu)
        levelmat <- rbind(levelmat,out)
        
        # sink(paste0("orig only regions/orig only regions_", kkk, ".txt"))
        # print(paste0("City changed because of Region: ", cinr[i], "; ","because of ",
        #              rr, "; ", pn, "; ",  mean(as.numeric(ratio1)), "; ", ratio22, "; ", original.point,"; ",
        #              pos.inc, "; ", pos.n))
        # sink()
        # kkk <- kkk+1
      }
    }
    
  }
  
  levelmat[,years2][levelmat[,years2]>100] <- 100
  
    
  p <- xframe[xframe$City=="Total",years]/100
  
  #███
  ## Paima visus šalies miestus ir regionus -  ir tada dar pakoreguoja
  citymustbetchanged <- intersect(unique(levelmat$RegionCode)[unique(levelmat$RegionCode)!=cc],xframe$CityCode)
  
  for(mcc in citymustbetchanged){
    
    ### čia mūsų modeliuotas possession iki dabar
    x <- xframe[xframe$CityCode==mcc,] 
    ### jei dar nesumodeliuovom, praleidžiam
    if(dim(x)[1]==0) next
    ### čia grynai originalūs arba sukurti 'fiktyvūs' pagal regionus ankstesniame cikle
    y <- levelmat[levelmat$RegionCode==mcc,] 
    ### parenkame metus, bet tik tuos kuriems turime kažkokius duomenis (research ar fiktyvius)
    nnay <- years2[which(!is.na(y[,years2]))]
    ### originalūs duomenys
    ymean <- mean(as.numeric(y[,nnay]))
    ### istoriniai duomenys
    xmean <- mean(as.numeric(x[,nnay]))
    
    ### čia jamam adityvų skirtumą
    x.aditive <- x[,years] - xmean + ymean
    ### čia multiplikatyvų skirtumą
    x.multiplicative <- x[,years] * ymean/xmean
      
    ### svorių parinkimo nesupratau - p yra gautas anksčiau tiesiog kaip possession???
    ### jei possession mažas, tai labiau ima multiplikatyvų. Jei possession didelis
    ### tada labiau žiūri į adityvų.
    ### 
    ### o neturėtų šitas priklausyt labiau nuo to, kaip nutolęs originalus nuo 
    ### istorinių duomenų? Nes jei nenutolęs, tai ir problemos gal nėra???
    x.out <- p * x.aditive + (1-p) * x.multiplicative
      
#     }
#         matplotg(cbind(as.numeric(x.aditive),as.numeric(x.multiplicative),as.numeric(x.out),as.numeric(xframe[xframe$CityCode==mcc,years])),abline=c(0,100))
    #     matplotg(cbind(as.numeric(x[,years]),as.numeric(x.out),out),stx=2005)
#         matplotg(cbind(as.numeric(x[,years]),as.numeric(x.out),as.numeric(out)),abline=c(0,100))

    x[,years] <- x.out
    
    xframe[xframe$CityCode==mcc,] <- x
    
  }
  
  # pramusimas over100.v2
  xframe[,years] <- over100.v2(xframe[,years],plot=F,divover=2,border = 1)
  # xframe[,years] <- over100.v2(xframe[,years],plot=F,divover=2,border = 1)
  
  if(dim(xframe)[1]==0) next
  
  ####################################
  
  xframe[,years] <- over100.v2(xframe[,years],plot=F,divover=2,border = 1)
  
  notchanged <- setdiff(xframe$CityCode[which(apply(oxframe[,years] - xframe[,years],1,mean) == 0)],c(cc,"Total"))
  if(length(notchanged)!=0) print(paste("There left few 'CityCodes' which was not changed or was equal (strange):",notchanged))
  
  #███
  rposout[rposout$CountryCode==cc & rposout$Possession==pn,] <- xframe
  #███
  
  #graphs
  gdata <- data.frame(t(xframe[,years]))
  names(gdata) <- xframe$City
  gdata2 <- gdata
  gdata2[,] <- t(oxframe[,years])
  
  
  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
  
  matplotf(gdata2,stx=2005,legend=F,text=T,abline=c(0,100),turnoffpar=T,
           main=paste(xframe$CountryCode[1],xframe$Possession[1],"Before",sep=" - "))
  lines(years,as.numeric(xframe[xframe$CityCode==cc,years]),lwd=3)
  
  matplotf(gdata,stx=2005,legend=F,text=T,abline=c(0,100),turnoffpar=T,
           main=paste(xframe$CountryCode[1],xframe$Possession[1],"After",sep=" - "))
  lines(years,as.numeric(xframe[xframe$CityCode==cc,years]),lwd=3)
  
  matplotf(gdata2,stx=2005,legend=F,text=T,abline=c(0,100),col="firebrick3",turnoffpar=T,
           main=paste(xframe$CountryCode[1],xframe$Possession[1],sep=" - "))
  matplotf(gdata,stx=2005,legend=F,text=T,abline=c(0,100),col="dodgerblue3",turnoffpar=T,
           main=paste(xframe$CountryCode[1],xframe$Possession[1],sep=" - "),add=T,lwd=2)
  lines(years,as.numeric(xframe[xframe$CityCode==cc,years]),lwd=3)
  legend("bottomright",legend = c("Before Originals","After Originals","Total"),col=c("firebrick3","dodgerblue3","black"),lty=1,cex=0.6)
  

}


dev.off()

#regionai, is kuriu padaryti orginalus miestu possessionai (yra orginalus region taskas bet nera city)
#dvieju variantu ratio sujungiam i viena file'a

# region.list <- NULL
# for (i in dir("orig only regions/")){
# 
#   region1 <- read.table(paste0("orig only regions/", i))
#   region.list <- rbind(region.list, region1)
# }
# region.list <- region.list[order(region.list$V2),]
# 
# sink("orig only regions all.txt")
# print(region.list$V2)
# sink()

###########################################################################
# testing conditions (satellite + cabel < television, etc)-----------------
###########################################################################

s4rposout <- rposout
save(s4rposout,file=paste0("temp/","s4rposout after originals adjustion ",gsub(":","-",Sys.time()),".Rdata"))

# source("regions durables conditions testing.R")
rposout <- test.all.conditions.cities(rposout)

cities.possession.graphs.grubus(bigdata=rposout,filename="plots/5. Forecasted possessions after originals and conditions/Forecasted possessions after originals and conditions ",
                                info1="Possessions of Durables",info2="Forecasts after adjustion using originals and 1 time conditions")


save(rposout,file=paste0("temp/","s4rposout after conditions and original adjustion",gsub(":","-",Sys.time()),".Rdata"))
# load("temp/s4rposout after conditions and original adjustion2016-06-06 08-54-08.RData")


###########################################################################
# galutiniai pataisymai (monotoniskumai ir glodinimai) --------------------
###########################################################################

s5.rposout <- rposout

# Exceptionu nuskaitymas, kad imtu kitoki pritempimo buda atitinkamai saliai x produktui
bt.other.method <- bt[bt$Issue=="originals.other.method",]
source("originals_othermethod.R")

rposout <- issue.smooth(bt, rposout)
rposout <- issue.stretch.f(bt, rposout)
rposout <- issue.total.ratio.f(bt, rposout)

################
# grafikai po pataisymu
graphs.Badthings.adjustion(bt, rposout, s5.rposout, 
                           filename= "plots/6. Badthings adjustion/Badthings adjustion ")


###########################################################################
# fixing where sum(city.poss.hh)>n.poss.hh --------------------------------
###########################################################################

# tikrinamame, ar pavertus absoliuciomis reiksmemis susumuoti miestu possession of durables nevirsija salies
hh <- read.csv(file="input/number of households 2015-11-23.csv",check.names = F,stringsAsFactors = F)

# bad.list <- list.cities.surpass.country(hh, rposout)
rposout <- fixing.cities.surpass.country(hh, rposout, plot=T)

###########################################################################
# testing conditions ------------------------------------------------------
###########################################################################

s6.rposout <- rposout
rposout <- test.all.conditions.cities(rposout,plot=F)

list.cities.surpass.country(hh, rposout) #turi buti tuscia jeigu ne reiktu kazka pergalvoti!!!!!!

###########################################################################
# Final Graphs ------------------------------------------------------------ 
###########################################################################

cities.possession.graphs.grubus(rposout,filename="plots/7. Forecasted possessions after exeptions/Forecasted possessions after exeptions ",
                                info1="Possessions of Durables",info2="Forecasts after exeptions and conditions check")
cities.possession.graphs.grubus.tik.miestai.split(rposout,splitmap,filename="plots/8. Forecasted possessions only cities mapping/Forecasted possessions only cities mapping ",
                                                  info1="Possessions of Durables",info2="Forecasts after exeptions and conditions cities splitted")


#  ------------------------------------------------------------------------
# final graphs ------------------------------------------------------------
#  ------------------------------------------------------------------------

bigdata <- rposout
cities.possession.graphs(bigdata)
# cities.possession.graphs.grubus(bigdata)
  

# output all ------------------------------------------------------------------
### pridedami 3 miestai, kuriu nebuvo anksciau, bet jie sutampa su salimi
BH01 <- add.cities.which.are.equal.to.countries(cc = "BH",citycode = "BH01",cityname = "Manama",pd,rposout)
KW01 <- add.cities.which.are.equal.to.countries(cc = "KW",citycode = "KW01",cityname = "Kuwait City",pd,rposout)
QT01 <- add.cities.which.are.equal.to.countries(cc = "QT",citycode = "QT01",cityname = "Doha",pd,rposout)

rposout <- rbind.fill(rposout,BH01,KW01,QT01)

rposout <- rposout[,c("CountryCode","CityCode","City","ProductID","Possession",years)]
outputname <- paste0("output/rposout final ",gsub(":","-",Sys.time()),".csv")
write.csv(rposout,file=outputname,row.names=F)

# cities output -----------------------------------------------------------

allrposout <- rposout

citiesposout <- rposout[rposout$CityCode %in% ccid$CityCode,]
citiesposout <- merge(citiesposout,ccid[,c("CityCode","CityCodeID")],all.x=T,by.x="CityCode",by.y="CityCode")

outputname <- paste0("output/cities output final ",gsub(":","-",Sys.time()),".csv")
write.csv(citiesposout,file=outputname,row.names=F)


citiesposoutoutput <- citiesposout
citiesposoutoutput <- citiesposoutoutput[,c("ProductID","CityCodeID",years)]
names(citiesposoutoutput) <- c("ProductID","CountryCodeID",paste0("Y",years))

outputname <- paste0("output/cities output with id final ",gsub(":","-",Sys.time()),".csv")
write.csv(citiesposoutoutput,file=outputname,row.names=F)


##check countries and durables where all cities are below country level

# all <- read.csv("output/rposout final 2015-11-24 14-25-22.csv", stringsAsFactors = F, check.names=F)
# 
# for (cn in setdiff(unique(all$CountryCode), c("BH", "QT", "KW"))){
#   
#   s <- subset(all, CountryCode==cn)
#   
#   for (pi in unique(s$ProductID)){
#     
#     ss <- subset(s, ProductID==pi)
#     
#     ss.city <- ss[ss$CityCode!=ss$CountryCode & !(substr(ss$CityCode, 1, 2) %in% paste0("R", 0:9)),]
#     
#     ss.city <- apply(ss.city[,as.character(c(2005:2030))], 2, max)
#     
#     ss.country <- ss[ss$CityCode==ss$CountryCode, as.character(c(2005:2030))]
#     
#     if (all(as.numeric(ss.city)<as.numeric(ss.country))){
#       
#       print(paste(cn, pi, sep=" - "))
#     }
#   }
# }

# Compare 2 data frames (this may be commented) ---------------------------


# xframe.new <- read.csv(file="output/cities output final 2015-11-24 14-26-56.csv",check.names=F,stringsAsFactors=F)
#   
# dirname2 <- "K:/GMID Research/CITY/1 - Data/0 - Country data dump/CitiesData"
# eval(parse(text=paste("xframe.old <- read.csv(",lastinput.zip(dirname2),", stringsAsFactors=F)")))
# xframe.old <- xframe.old[xframe.old$ProductID %in% unique(xframe.new$ProductID),]
# 
# xframe.old <- xframe.old[, c("ProductID", "ProductName", "CityCountryCode", "CityCode", "CityName", 
#                      paste0("Y", c(2005:2030)))]
# xframe.old <- rename(xframe.old, c(CityName="City", CityCountryCode="CityCodeID", CityName="City",
#                                    ProductName="Possession"))
# names(xframe.old) <- gsub("^Y", "", names(xframe.old))
# 
# xframe.old <- merge(xframe.old, city.codes[, c("CityCode", "CountryCode")], all.x=T)
# xframe.old[is.na(xframe.old$CountryCode),] #has to be empty
# 
# xframe.old <- xframe.old[!(xframe.old$CountryCode=="SA" & xframe.old$ProductID==12980),]
# 
# xframe.old$version <- "old"
# xframe.new$version <- "new"
# 
# compare <- rbind(xframe.old, xframe.new)
# compare$id <- paste(compare$CityCode, compare$ProductID)
# 
# years <- as.character(c(2005:2030))
# 
# compare$changes <- NA
# 
# for (i in unique(compare$id)){
#   
#   old <- mean(as.numeric(compare[compare$id==i & compare$version=="old", years]))
#   new <- mean(as.numeric(compare[compare$id==i & compare$version=="new", years]))
#   
#   changed <- abs((new-old)/old)*100
#   compare[compare$id==i,]$changes <- changed
# }
# 
# 
# write.csv(compare, "output/chnages comparing with old poss.csv", row.names=F)



# changed.frame <- compare2data.frames(xframe.new,xframe.old)



