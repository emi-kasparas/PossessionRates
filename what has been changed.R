# old <- read.csv("output/cities output with id final 2015-04-16 15-49-54.csv", stringsAsFactors=F)
# d.new <- read.csv("output/cities output with id final 2015-05-19 13-24-19.csv", stringsAsFactors=F)
# city.codes <- read.csv("input/city codes id.csv", stringsAsFactors=F)
# city.codes <- rename(city.codes, c(CityCodeID="CountryCodeID"))
# p.codes <- read.csv("input/possession_codes.csv", stringsAsFactors=F)
# 
# setdiff.data.frame <-
#   function(A,B) A[ !duplicated( rbind(B,A) )[ -seq_len(nrow(B))] , ]
# 
# changed.new <- setdiff.data.frame(d.new, old)
# changed.old <- setdiff.data.frame(old, d.new)
# 
# changed.new$type <- "new"
# changed.old$type <- "old"
# 
# changed <- rbind(changed.new, changed.old)
# 
# changed <- merge(changed, city.codes, by=c("CountryCodeID"))
# changed <- merge(changed, p.codes, by=c("ProductID"))
# 
# changed <- changed[,c("type","CityCode", "CityName", "ProductName", "ProductID", "CountryCodeID", paste0("Y", c(2005:2030)))]
# 
# write.csv(changed, paste0("output/cities changed_", Sys.Date(),".csv"), row.names=F)

d.new <- read.csv("output/cities output final 2016-06-10 00-56-45.csv", stringsAsFactors=F)
d.old <- read.csv("input/PossData2016-05-03.csv", stringsAsFactors = F)

d.new <- melt(d.new, id=c("CityCode", "CountryCode", "City", "ProductID", "Possession", "CityCodeID"))
d.new <- rename(d.new, c(variable="Year", value="new", Possession="ProductName"))
d.new$Year <- sub("X", "", d.new$Year)
d.new$ProductName <- NULL

d.old <- d.old[,c("ProductID", "ProductName", "CityCode", paste0("Y", 2005:2030))]
d.old <- melt(d.old, id=c("ProductID", "ProductName", "CityCode"))
d.old$variable <- sub("Y", "", d.old$variable)
d.old <- rename(d.old, c(variable="Year", value="old"))

data <- merge(d.new, d.old, by=c("ProductID", "CityCode", "Year"))

data$diff.old.new <- abs((round(data$new,2)-round(data$old,2))/round(data$old,2)*100)

##output

data.out <- melt(data, id=c("ProductID", "CityCode", "Year", "CountryCode", "City",
                            "CityCodeID", "ProductName"))
data.out <- dcast(data.out, ProductID+CityCode+CountryCode+City+CityCodeID+ProductName+
                    variable~Year, value.var="value")

write.csv(data.out, paste0("output/Compare_old_and_new ", Sys.Date(), ".csv"), row.names = F)


#add totals (old and new) to data
load("input/Possessions by Decile/final_poss_by_decile RData/final_poss_by_decile 2016-06-03.RData")
tnew <- final.pbd

load("input/Possessions by Decile/final_poss_by_decile RData/final_poss_by_decile 2015-12-10.RData")
told <- final.pbd

tnew$pos <- "new"
told$pos <- "old"

tt <- rbind(tnew, told)
tt <- tt[,c("Year", "Decile.No", "CountryCode", "ProductID", "est.Pos.norm", "pos")]

tt <- ddply(tt, .(Year, CountryCode, ProductID, pos), summarise, t.poss=mean(est.Pos.norm))
tt <- dcast(tt, Year+CountryCode+ProductID~pos, value.var="t.poss")

tt <- rename(tt, c(new="t.new", old="t.old"))

data <- merge(data, tt, by=c("Year", "CountryCode", "ProductID"), all.x=T)


#duomenys,kurie pasikeite daugiau nei per 10proc
changes <- unique(data[data$diff.old.new>10,c("CityCode", "ProductID", "CountryCode", "ProductName", "City")])

changes$id <- paste0(changes$CityCode, changes$ProductID, changes$CountryCode,
                     changes$ProductName, changes$City)
data$id <- paste0(data$CityCode, data$ProductID, data$CountryCode,
                  data$ProductName, data$City)

data.diff10proc <- data[data$id %in% changes$id,]
data.diff10proc$id <- NULL
data$id <- NULL

##graphs

comparing.old.new.data(data, savedir = "plots/9. Comparing old and new/",
                       filename = "Comparing old and new poss ")

comparing.old.new.data(data.diff10proc, savedir = "plots/9. Comparing old and new/",
                       filename = "Comparing old and new poss diff 10proc ")


