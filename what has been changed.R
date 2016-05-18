old <- read.csv("output/cities output with id final 2015-04-16 15-49-54.csv", stringsAsFactors=F)
d.new <- read.csv("output/cities output with id final 2015-05-19 13-24-19.csv", stringsAsFactors=F)
city.codes <- read.csv("input/city codes id.csv", stringsAsFactors=F)
city.codes <- rename(city.codes, c(CityCodeID="CountryCodeID"))
p.codes <- read.csv("input/possession_codes.csv", stringsAsFactors=F)

setdiff.data.frame <-
  function(A,B) A[ !duplicated( rbind(B,A) )[ -seq_len(nrow(B))] , ]

changed.new <- setdiff.data.frame(d.new, old)
changed.old <- setdiff.data.frame(old, d.new)

changed.new$type <- "new"
changed.old$type <- "old"

changed <- rbind(changed.new, changed.old)

changed <- merge(changed, city.codes, by=c("CountryCodeID"))
changed <- merge(changed, p.codes, by=c("ProductID"))

changed <- changed[,c("type","CityCode", "CityName", "ProductName", "ProductID", "CountryCodeID", paste0("Y", c(2005:2030)))]

write.csv(changed, paste0("output/cities changed_", Sys.Date(),".csv"), row.names=F)
