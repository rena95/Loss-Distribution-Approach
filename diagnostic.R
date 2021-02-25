# Graphs to set up the threshold used to split body and tail in the severity distribution for the seven ETs
rm(list=ls())
library(fExtremes,verbose=FALSE,warn.conflicts=FALSE)
library(evir,quietly=QUIETLY,verbose=FALSE,warn.conflicts=FALSE)
library("dplyr")
library("ggplot2")
library("tidyr")
# Import data in a list: each element contain the data of an ET
load("dati_list.RData")
try(expr, silent = TRUE)
ets <- paste0("ET",1:7)

################################################
######### COMPUTE THE xi with ML method########
################################################

 parameters.ml <- lapply(ets,function(et){
            U_x <- round(seq(0.90,0.999,length.out = 100),4)
					  lapply(U_x,function(p){
					  	U 		<- vector()
					  	U 		<- gpdFit(dati.list[[et]]$PTL,u=quantile(dati.list[[et]]$PTL,p,na.rm=TRUE),"mle")@parameter$u
					  	xi 		<- vector()
					  	xi 		<- gpdFit(dati.list[[et]]$PTL,u=quantile(dati.list[[et]]$PTL,p,na.rm=TRUE),"mle")@fit$par.ests[[1]]
					  	beta 	<- vector()
					  	beta 	<- gpdFit(dati.list[[et]]$PTL,u=quantile(dati.list[[et]]$PTL,p,na.rm=TRUE),"mle")@fit$par.ests[[2]]
					  	exceed 	<- vector()
					  	exceed  <- gpdFit(dati.list[[et]]$PTL,u=quantile(dati.list[[et]]$PTL,p,na.rm=TRUE),"mle")@data$exceedances
					  	return(list(U=U,xi=xi,beta=beta,exceed=exceed))
					  	 })
          })
 names(parameters.ml) <- ets
 ################################################
 ######### COMPUTE THE xi with PW method########
 ################################################	
 
 parameters.pw <- lapply(ets,function(et){
     U_x <- round(seq(0.90,0.999,length.out = 100),4)
     lapply(U_x,function(p){
       U 		<- vector()
       U 		<- gpdFit(dati.list[[et]]$PTL,u=quantile(dati.list[[et]]$PTL,p,na.rm=TRUE),"pwm")@parameter$u
       xi 		<- vector()
       xi 		<- gpdFit(dati.list[[et]]$PTL,u=quantile(dati.list[[et]]$PTL,p,na.rm=TRUE),"pwm")@fit$par.ests[[1]]
       beta 	<- vector()
       beta 	<- gpdFit(dati.list[[et]]$PTL,u=quantile(dati.list[[et]]$PTL,p,na.rm=TRUE),"pwm")@fit$par.ests[[2]]
       exceed 	<- vector()
       exceed  <- gpdFit(dati.list[[et]]$PTL,u=quantile(dati.list[[et]]$PTL,p,na.rm=TRUE),"pwm")@data$exceedances
       return(list(U=U,xi=xi,beta=beta,exceed=exceed))
             })
         })
 names(parameters.pw) <- ets

###############################################################################################
###DRAWN GRAPH OF Xi ESTIMATED WITH Hill METHOD AS FUNCTION OF 'U' FOR THE 7 EVENT TYPES###
###############################################################################################

if(do.plot){
  lapply(ets, function(et){
    dir <- paste0("../img/chapter_2/Hill_plot_",et,".png")
    png(filename = dir,width = 1500,height = 500)
    xi.hill <- evir::hill(dati.list[[et]]$PTL,start = 20,end=120, option = "xi")
    dev.off()
    
  })
}  


########################################################
###COMPUTE THE PARAMETERS FOR THE CHOOSEN THRESHOLD#############
########################################################
U <- c(6.7e5,2.3e5,2.45e5,1.3e6,0.43e5,0.95e5,4.5e5)
gpd_par <- list()

for(et in 1:length(ets)){
  xi 		  <- gpdFit(dati.list[[ets[et]]]$PTL,u=U[et],"pwm")@fit$par.ests[[1]]
  beta   	<- gpdFit(dati.list[[ets[et]]]$PTL,u=U[et],"pwm")@fit$par.ests[[2]]
  exceed 	<- vector()
  exceed  <- gpdFit(dati.list[[ets[et]]]$PTL,u=U[et],"pwm")@data$exceedances
  gpd_par[[et]]<- list(xi,beta,exceed )
  names(gpd_par[[et]]) <- c("xi","beta","excesses")
}
names(gpd_par) <- ets

save(gpd_par,file="gpd_parameters.RData")
