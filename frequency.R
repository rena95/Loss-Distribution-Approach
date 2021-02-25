library("lubridate")
library("tidyr")
library("MASS")
# Import the data
# Data are in the following format i.e. LONG FORMAT
#  CODICE  ETLEV1   PTL        DATA
# 1 80478   1   426572.79 2014-08-21
# 2 728215  3   68087.18  2013-02-01
# 3 862286  4   347551.44 2008-05-07
# 4 4679476 1   5270.05  2011-04-18
# 5 1247    7   52112.14  2011-05-19
load("dati.RData")

# #Prepare the data - create a list where each element
# # is a dataframe with the number of events for ET
 prepara.frequenza<- function(data){
   dataFrequency <- list(dim = 7)
   for(i in 2:8){
    dataFrequency[[i-1]]<- cbind(data[1],data[i])
    colnames(dataFrequency[[i-1]]) <- c("DATA_RIL","LOSS")
   }
   names(dataFrequency)<- ets
   return(dataFrequency)
 }

 
dati_m <- dati %>% mutate(DATA_RIL_M = paste0(substr(DATA_RIL,1,8),"01"))
# Short format:
# in the following format i.e. SHORT FORMAT
#  DATA_RIL        ET1       ET2         ET3          ET4        ET5         ET6       ET7
# 1 2007-01-01  1580852.61  13125504    5115556.9  4610128     291992.90  33255.07  3379720
# 2 2007-02-01  47430.63  4700840     1168138.7    1901531    48842.28   56339.82  1895047
# 3 2007-03-01  562063.16 11638957  4199101.9     79446932    1018598.83   89886.24  5952838
dati_m <- dati_m %>% 
  group_by(ETLEV1, DATA_RIL_M) %>%
  summarise(n=n()) %>% 
  spread(ETLEV1,n) %>% 
  arrange(DATA_RIL_M) %>% data.frame  
  
dati_m<- prepara.frequenza(dati_m)
dati_m <- lapply(dati_m,na.omit)
# Estimation of lambda for monthly losses occurrences
freq<- lapply(dati_m,function(et){
  parameter <- fitdistr(et$LOSS, densfun = "Poisson")
  return(list(lambda=parameter))
})
#Save monthly losses
save(freq,file="frequency.RData")
