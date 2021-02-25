library(copula)
library(dplyr)
library(QRM)
library(evir)
library(goftest)
library(parallel)
library(ggplot2)
source("copula_lib.R")
source("fitting.R")
# load the dataset where observation's losses are monthly aggregated
# i.e. these are the observations used to estimate the parameters of the copulas
load("dati_pseudo_spread_monthly.RData")
# Dataset used in the silumation
load("dati_list.RData")
#Generalised Pareto distribution parameter - xi and beta
load("gpd_parameters.RData")
#load the fitted lambda
load("frequency.RData")
# load the t copula where correlation matrix with 12 df is estimated under ML
# load Clayton - Frank - Gumbel copula where alpha is estimated under ML
load("fitted_parameter_ml.RData")
############PARAMETERS############
#Selected threshold for the 7 ETs
U         <- c(6.7e5,2.3e5,2.45e5,1.3e6,0.43e5,0.95e5,4.5e5)
ets       <- paste0("ET",1:7)
nSim      <- 5000000
nBlock    <- 500
do.margin.short <- FALSE
do.fit          <- FALSE
do.plot         <- FALSE
set.seed(8052016) 

# Find the correlation matrix as Example 3 - Chapter 3
# If the third argument "type.delta" is set equal:
# "cov" the perturbations of Kendall tau come from a normal covariance matrix
# else  from  unif[0,b] random variables with b=0.5
copula <- copula(dati.pseudo.spread.monthly,10,"un")
names(copula) <- c("monthly_data","monte_Carlo_perturbation"
                  ,"emp_tail_copulae","distance","min_corr_matrix","df")



#Compute the margin
margin <- function(dati,lambda,k,xi,beta,n){
  tot_loss   <- vector()
  sim_loss   <- vector()
  n_ext      <- vector()
  n_ev <- rpois(n,lambda)
  #Number of events less than the threshold out all events
  F_u <- length(dati[dati<k])/length(dati)
  n_coda <- 0
  for(i in 1:n){
    if(n_ev[i]==0){next}
    else{
      U        <- runif(n_ev[i])
      for(j in 1:length(U)){
        if(U[j] < F_u){
           U[j] <- (U[j]/F_u)
           sim_loss[j] <- dati[ceiling(length(dati)* U[j])]#sim_loss <- dati[[et]]$PTL[length(dati[[et]]$PTL)* U]
        }
        else{
          U[j] <- ((U[j]-F_u)/(1-F_u))
          sim_loss[j]    <- k + evir::qgpd(U[j], xi = xi, beta = beta)
        }
      }
    tot_loss[i] <- sum(sim_loss)
    sim_loss <- vector()
     }
  }
  # Using our data, is.na(tot_loss) is never TRUE
  tot_loss[is.na(tot_loss)] <-0
  return(tot_loss)
}

# NOTE: The following part is really heavy and time-consuming
#-----------------------------------------------------------
# Compute the aggregate annual loss distribution for ET1
# 500 block of 10000 aggregated losses generated
    et <- 1
    lenBlock    <-  nSim/nBlock #=10000
    # My laptop has 4 cores
    cl <- makeCluster(4)
    # Only for window!
    # Pass the variable and export 
    # these variables to the other R processes 
    # in the cluster
    clusterExport(cl, varlist = c("et","dati.list","freq",
                                  "U","gpd_par","margin",
                                  "lenBlock"))
    # Export the package evir as the function evir::qgpd
    # it is used in the code
    clusterEvalQ(cl, library(evir))
    # Parallel working
    aggr_loss_ET1 <- parallel::parLapply(cl,1:nBlock, function(i){
        #  lambda is estimated on monthly observations, thus it must 
        # be multiplied by 12
        margin(dati.list[[et]]$PTL,freq[[et]]$lambda$estimate*12
               ,U[et],gpd_par[[et]]$xi,gpd_par[[et]]$beta,lenBlock)
      })

     stopCluster(cl)
     
  # Compute 500 Vars for the 500 blocks
    Vars_ET1  <- lapply(aggr_loss_ET1, function(i){
         i <- quantile(i,0.999,na.rm=TRUE)
       })
    Var_bar_ET1 <- mean(unlist(Vars_ET1)) 
    
  save(aggr_loss_ET1,file = "aggr_loss_ET1.RData")
  
  data_VaR <- data.frame(n=1:500,
                         VaR=unlist(Vars_ET1))
  ggplot( data_VaR, aes(x=VaR/1e6))+
    geom_histogram(aes(y=..density..), colour="black", fill="white",bins=50)+
    geom_density(alpha=.5, fill="#FF6666",adjust=2)+ 
    scale_x_continuous(limits = c(0, 2500))+
    labs(title=paste0("Asymptotic distribution for ET1"),x="VaR")+
    theme_classic()
  ggsave(filename =paste0("../img/chapter_4/ET1/CLT_VaR1.png"))
  
#-----------------------------------------------------
  # Compute the aggregate annual loss distribution for ET2
  # 500 block of 10000 aggregated losses generated
  et <- 2
  lenBlock    <-  nSim/nBlock #=10000
  # My laptop has 4 cores
  cl <- makeCluster(4)
  clusterExport(cl, varlist = c("et","dati.list","freq",
                                "U","gpd_par","margin",
                                "lenBlock"))
  clusterEvalQ(cl, library(evir))
  # Parallel working
  aggr_loss_ET2 <- parallel::parLapply(cl,1:nBlock, function(i){
    margin(dati.list[[et]]$PTL,freq[[et]]$lambda$estimate*12
           ,U[et],gpd_par[[et]]$xi,gpd_par[[et]]$beta,lenBlock)
  })
  
  stopCluster(cl)
  # Compute 500 Vars for the 500 blocks
  Vars_ET2  <- lapply(aggr_loss_ET2, function(i){
    i <- quantile(i,0.999,na.rm=TRUE)
  })
  Var_bar_ET2 <- mean(unlist(Vars_ET2)) 
  
  save(aggr_loss_ET2,file = "aggr_loss_ET2.RData")
  # Plot Vars ET2
  data_VaR <- data.frame(n=1:500,
                         VaR=unlist(Vars_ET2))
  ggplot( data_VaR, aes(x=VaR/1e6))+
    geom_histogram(aes(y=..density..), colour="black", fill="white",bins=50)+
    geom_density(alpha=.5, fill="#FF6666",adjust=2)+ 
    scale_x_continuous(limits = c(80, 120))+
    labs(title=paste0("Asymptotic distribution for ET2"),x="VaR")+
    theme_classic()
  ggsave(filename =paste0("../img/chapter_4/ET2/CLT_VaR2.png"))
  #-----------------------------------------------------
  # Compute the aggregate annual loss distribution for ET3
  # 500 block of 10000 aggregated losses generated
  et <- 3
  lenBlock    <-  nSim/nBlock #=10000
  # My laptop has 4 cores
  cl <- makeCluster(4)
  clusterExport(cl, varlist = c("et","dati.list","freq",
                                "U","gpd_par","margin",
                                "lenBlock"))
  clusterEvalQ(cl, library(evir))
  # Parallel working
  aggr_loss_ET3 <- parallel::parLapply(cl,1:nBlock, function(i){
    margin(dati.list[[et]]$PTL,freq[[et]]$lambda$estimate*12
           ,U[et],gpd_par[[et]]$xi,gpd_par[[et]]$beta,lenBlock)
  })
  
  stopCluster(cl)
  # Compute 500 Vars for the 500 blocks
  Vars_ET3  <- lapply(aggr_loss_ET3, function(i){
    i <- quantile(i,0.999,na.rm=TRUE)
  })
  Var_bar_ET3 <- mean(unlist(Vars_ET3)) 
  
  save(aggr_loss_ET3,file = "aggr_loss_ET3.RData")
  # Plot Vars ET3
  data_VaR <- data.frame(n=1:500,
                         VaR=unlist(Vars_ET3))
  ggplot( data_VaR, aes(x=VaR/1e6))+
    geom_histogram(aes(y=..density..), colour="black", fill="white",bins=50)+
    geom_density(alpha=.5, fill="#FF6666",adjust=2)+ 
    scale_x_continuous(limits = c(40, 55))+
    labs(title=paste0("Asymptotic distribution for ET3"),x="VaR")+
    theme_classic()
  ggsave(filename =paste0("../img/chapter_4/ET3/CLT_VaR3.png"))
  #-----------------------------------------------------
  # Compute the aggregate annual loss distribution for ET4
  # 500 block of 10000 aggregated losses generated
  et <- 4
  lenBlock    <-  nSim/nBlock #=10000
  # My laptop has 4 cores
  cl <- makeCluster(4)
  clusterExport(cl, varlist = c("et","dati.list","freq",
                                "U","gpd_par","margin",
                                "lenBlock"))
  clusterEvalQ(cl, library(evir))
  # Parallel working
  aggr_loss_ET4 <- parallel::parLapply(cl,1:nBlock, function(i){
    margin(dati.list[[et]]$PTL,freq[[et]]$lambda$estimate*12
           ,U[et],gpd_par[[et]]$xi,gpd_par[[et]]$beta,lenBlock)
  })
  
  stopCluster(cl)
  # Compute 500 Vars for the 500 blocks
  Vars_ET4  <- lapply(aggr_loss_ET4, function(i){
    i <- quantile(i,0.999,na.rm=TRUE)
  })
  Var_bar_ET4 <- mean(unlist(Vars_ET4)) 
  
  save(aggr_loss_ET4,file = "aggr_loss_ET4.RData")
  # Plot Vars ET4
  data_VaR <- data.frame(n=1:500,
                         VaR=unlist(Vars_ET4))
  ggplot( data_VaR, aes(x=VaR/1e6))+
    geom_histogram(aes(y=..density..), colour="black", fill="white",bins=50)+
    geom_density(alpha=.5, fill="#FF6666",adjust=2)+ 
    scale_x_continuous(limits = c(500, 2000))+
    labs(title=paste0("Asymptotic distribution for ET4"),x="VaR")+
    theme_classic()
  ggsave(filename =paste0("../img/chapter_4/ET4/CLT_VaR4.png"))

  #-----------------------------------------------------
  # Compute the aggregate annual loss distribution for ET5
  # 500 block of 10000 aggregated losses generated
  et <- 5
  lenBlock    <-  nSim/nBlock #=10000
  # My laptop has 4 cores
  cl <- makeCluster(4)
  clusterExport(cl, varlist = c("et","dati.list","freq",
                                "U","gpd_par","margin",
                                "lenBlock"))
  clusterEvalQ(cl, library(evir))
  # Parallel working
  aggr_loss_ET5 <- parallel::parLapply(cl,1:nBlock, function(i){
    margin(dati.list[[et]]$PTL,freq[[et]]$lambda$estimate*12
           ,U[et],gpd_par[[et]]$xi,gpd_par[[et]]$beta,lenBlock)
  })
  
  stopCluster(cl)
  # Compute 500 Vars for the 500 blocks
  Vars_ET5  <- lapply(aggr_loss_ET5, function(i){
    i <- quantile(i,0.999,na.rm=TRUE)
  })
  Var_bar_ET5 <- mean(unlist(Vars_ET5)) 
  
  save(aggr_loss_ET5,file = "aggr_loss_ET5.RData")  
  # Plot Vars ET5
  data_VaR <- data.frame(n=1:500,
                         VaR=unlist(Vars_ET5))
  ggplot( data_VaR, aes(x=VaR/1e6))+
    geom_histogram(aes(y=..density..), colour="black", fill="white",bins=50)+
    geom_density(alpha=.5, fill="#FF6666",adjust=2)+ 
    scale_x_continuous(limits = c(9, 12))+
    labs(title=paste0("Asymptotic distribution for ET5"),x="VaR")+
    theme_classic()
  ggsave(filename =paste0("../img/chapter_4/ET5/CLT_VaR5.png"))
  
  
  #-----------------------------------------------------
  # Compute the aggregate annual loss distribution for ET6
  # 500 block of 10000 aggregated losses generated
  et <- 6
  lenBlock    <-  nSim/nBlock #=10000
  # My laptop has 4 cores
  cl <- makeCluster(4)
  clusterExport(cl, varlist = c("et","dati.list","freq",
                                "U","gpd_par","margin",
                                "lenBlock"))
  clusterEvalQ(cl, library(evir))
  # Parallel working
  aggr_loss_ET6<- parallel::parLapply(cl,1:nBlock, function(i){
    margin(dati.list[[et]]$PTL,freq[[et]]$lambda$estimate*12
           ,U[et],gpd_par[[et]]$xi,gpd_par[[et]]$beta,lenBlock)
  })
  
  stopCluster(cl)
  # Compute 500 Vars for the 500 blocks
  Vars_ET6  <- lapply(aggr_loss_ET6, function(i){
    i <- quantile(i,0.999,na.rm=TRUE)
  })
  Var_bar_ET6 <- mean(unlist(Vars_ET6)) 
  
  save(aggr_loss_ET6,file = "aggr_loss_ET6.RData")  
  # Plot Vars ET6
  data_VaR <- data.frame(n=1:500,
                         VaR=unlist(Vars_ET6))
  ggplot( data_VaR, aes(x=VaR/1e6))+
    geom_histogram(aes(y=..density..), colour="black", fill="white",bins=50)+
    geom_density(alpha=.5, fill="#FF6666",adjust=2)+ 
    scale_x_continuous(limits = c(30, 180))+
    labs(title=paste0("Asymptotic distribution for ET6"),x="VaR")+
    theme_classic()
  ggsave(filename =paste0("../img/chapter_4/ET6/CLT_VaR6.png"))
  
  #-----------------------------------------------------
  # Compute the aggregate annual loss distribution for ET7
  # 500 block of 10000 aggregated losses generated
  et <- 7
  lenBlock    <-  nSim/nBlock #=10000
  # My laptop has 4 cores
  cl <- makeCluster(4)
  clusterExport(cl, varlist = c("et","dati.list","freq",
                                "U","gpd_par","margin",
                                "lenBlock"))
  clusterEvalQ(cl, library(evir))
  # Parallel working
  aggr_loss_ET7<- parallel::parLapply(cl,1:nBlock, function(i){
    margin(dati.list[[et]]$PTL,freq[[et]]$lambda$estimate*12
           ,U[et],gpd_par[[et]]$xi,gpd_par[[et]]$beta,lenBlock)
  })
  
  stopCluster(cl)
  # Compute 500 Vars for the 500 blocks
  Vars_ET7 <- lapply(aggr_loss_ET7, function(i){
    i <- quantile(i,0.999,na.rm=TRUE)
  })
  Var_bar_ET7 <- mean(unlist(Vars_ET7)) 
  
  save(aggr_loss_ET7,file = "aggr_loss_ET7.RData")    
  # Plot Vars ET7
  data_VaR <- data.frame(n=1:500,
                         VaR=unlist(Vars_ET7))
  ggplot( data_VaR, aes(x=VaR/1e6))+
    geom_histogram(aes(y=..density..), colour="black", fill="white",bins=50)+
    geom_density(alpha=.5, fill="#FF6666",adjust=2)+ 
    scale_x_continuous(limits = c(400, 900))+
    labs(title=paste0("Asymptotic distribution for ET7"),x="VaR")+
    theme_classic()
  ggsave(filename =paste0("../img/chapter_4/ET7/CLT_VaR7.png"))

# Version Shortcode without parallel computing
if(do.margin.short){
  # the element in sim_margin[[i]][[j]] is what 
  # we called S^i_j - see chapter 4.1 
  sim_margin <- lapply(1:7, function(et){
    lenBlock    <-  nSim/nBlock
    m <- list()
    VaR    <- vector()
    lapply(1:nBlock, function(i){
      #  lambda is estimated on monthly observations, thus it must 
      # be multiplied by 12
      margin(dati.list[[et]]$PTL,freq[[et]]$lambda$estimate*12
             ,U[et],gpd_par[[et]]$xi,gpd_par[[et]]$beta,n=lenBlock)
    })
  }) 
}
#-----------------------------------------------------------

#Kolmogorov-Smirnov test for normality of the VaRs_block
km_Vars<- list()
for(i in 1:7){
  #extract the p_value of the Kolmogorov-Smirnov 
  km_Vars[[i]] <- ks.test(unlist(get(paste0("Vars_ET",i)) ),"pnorm",
          mean(unlist(get(paste0("Vars_ET",i)))),
          sd(unlist(get(paste0("Vars_ET",i)))))[[2]]
} 
km_test <- data.frame(ET=paste0("ET",1:7),
                      km_pvalue=round(unlist(km_Vars),3)) 


# Join the aggregate annual loss distribution
# in one list
sim_margin <- lapply(1:7,function(et){
  get(paste0("aggr_loss_ET",et))
})
names(sim_margin) <- paste0("ET",1:7)

# Join all the Vars in a list 
Var.perfect_corr <- lapply(1:7,function(et){
  get(paste0("Vars_ET",et))
})
names(Var.perfect_corr) <- paste0("ET",1:7)  

# Join all the Economic Capital for the 7 ETs in a list 
capital_charge.perf_cor <- lapply(1:7,function(et){
  get(paste0("Var_bar_ET",et))
})
names(capital_charge.perf_cor) <- paste0("ET",1:7)  




if(do.fit){
  # Check the fit of the distributions in "sim_margin " by QQ plots:
  # empirical dfs vs simulated body-tail df
  # empirical dfs vs lognormal df
  # empirical dfs vs normal    df
  # empirical dfs vs exponential df
  lapply(1:7, function(et){
    fit_margin <- fit_QQ(dati.list[[et]]$PTL,U[et],gpd_par[[et]]$xi,gpd_par[[et]]$beta,1000) 
    png(filename = paste0("../img/chapter_4/ET",et,"/qq_edfVSsim.png"))
    qqplot(fit_margin[[1]],fit_margin[[2]],ylab="Empirical df quantiles",
           xlab="Body-tail df quantiles")
    abline(a=0,b=1)
    dev.off()
    png(filename =  paste0("../img/chapter_4/ET",et,"/qq_edfVSlog.png"))
    qqplot(fit_margin[[1]],fit_margin[[3]],ylab="Empirical df quantiles",
           xlab="Log-normal df quantiles")
    abline(a=0,b=1)
    dev.off() 
    png(filename =  paste0("../img/chapter_4/ET",et,"/qq_edfVSnorm.png"))
    qqplot(fit_margin[[1]],fit_margin[[4]],ylab="Empirical df quantiles",
           xlab="Normal df quantiles")
    abline(a=0,b=1)
    dev.off()
    png(filename =  paste0("../img/chapter_4/ET",et,"/qq_edfVSexp.png"))
    qqplot(fit_margin[[1]],fit_margin[[5]],ylab="Empirical df quantiles",
           xlab="Exponential df quantiles")
    abline(a=0,b=1)
    dev.off()
  })
  # Anderson-Darling test for the seven ETs
  A.D <- lapply(1:7, function(et){
    dati <- dati.list[[et]]$PTL-U[et]
    sel  <- which(dati>0)
    dati <- dati[sel]
    goftest::ad.test(dati,"pgpd",xi=gpd_par[[et]]$xi,
                     beta=gpd_par[[et]]$beta)
  })
  # Kolmogorov-Smirnov test for the seven ETs
  K.S  <- lapply(1:7, function(et){
    dati <- dati.list[[et]]$PTL-U[et]
    sel  <- which(dati>0)
    dati <- dati[sel]
    goftest::ks.test(dati,"pgpd",xi=gpd_par[[et]]$xi,
                     beta=gpd_par[[et]]$beta)
  })
  # Summarize the statistics and p-value for Anderson-Darling
  # and Kolmogorov-Smirnov 
  distr_free_test <- data.frame(
    ET        = paste0("ET",1:7),
    km_test   = sapply(1:7, function(et){K.S[[et]][[1]]}),
    km_pvalue = sapply(1:7, function(et){K.S[[et]][[2]]}),
    ad_test   = sapply(1:7, function(et){A.D[[et]][[1]]}),
    ad_pvalue = sapply(1:7, function(et){A.D[[et]][[2]]})
    
  )
  # Sensitivity of the VaR vs 50 quantiles from 0.99 to 0.999
  # Note that we refer to the loss distribution and not to
  # frequency/severity aggregated loss distribution
  for(et in 1:7){
   distr <- sapply(seq(0.99,0.999,length.out = 50),function(p){
    quantile(dati.list[[et]]$PTL,p)})

   VaR  <- data.frame(q=seq(0.99,0.999,length.out = 50), VaR=(distr/1000))

  ggplot(data=VaR, aes(x=q, y=VaR, group=1))+
    geom_line()+
    #geom_errorbar(aes(ymin=VaR-SD, ymax=VaR+SD), width=.1) +
    scale_x_continuous(breaks=round(seq(0.99,0.999,length.out=20),3))+
    scale_y_continuous(breaks=round(seq(min(VaR$VaR),max(VaR$VaR),length.out=10),0))+
    theme(axis.text.x = element_text(angle=45))+
    labs(x="Confidence level",y = paste0("Quantile ET",et," in k/euro"))+
    theme_bw()
    ggsave(filename =paste0("../img/chapter_4/ET",et,"/VaR_severity_path_ET_",et,".png"),width=7,height=5) 
    }
}



if(do.plot){
  for (et in 1:7){
  # drawn how the VaR evolves over different confidence intervals
  # inizialize the matrix
  var.path <- matrix(nrow = 100,ncol = 50)
  # compute the VaR for each confidence level
  # the i-th row of var.path containes the 50 quantiles
  # of the i-th block
  for (i in 1:100) {
    var.path[i,]<- unlist(quantile(sim_margin[[et]][[i]],seq(0.99,0.999,length.out=50)))
  }
  # take the column average of var.path[i,] so as to apply the central limit theorem
  VaR  <- data.frame(q=seq(0.99,0.999,length.out=50),VaR=(apply(var.path, MARGIN=2, FUN=mean)/1e6))
  ggplot(data=VaR, aes(x=q, y=VaR, group=1))+
    geom_line()+
    scale_x_continuous(breaks=round(seq(0.99,0.999,length.out=10),3))+
    scale_y_continuous(breaks=round(seq(min(VaR$VaR),max(VaR$VaR),length.out=20),0))+
    theme(axis.text.x = element_text(angle=45))+
    labs(x="Confidence level",y = paste0("Aggregate annual loss quantiles for ET",et," in mln/euro"))+
    theme_bw()
    ggsave(filename =paste0("../img/chapter_4/ET",et,"/VaR_path_ET_",et,".png"),width=7,height=5)
  }
}


# Compute the computational error as shown
# in section 4.1
CI_builder <- function(VaR,crit_value){
  sd   <- sd(VaR)
  lim_CIM <- crit_value*sd/sqrt(nBlock)
    }
CE <- lapply(VaRs, function(i){
          CI_builder(unlist(i),1.96)/mean(unlist(i))
})



U <- function(lenBlock,type.dependence){
    if(type.dependence=="t"){
      # take the off-diagonal elements of the correlation matrix estimated 
      # under minimum distance method in chapter 3
      P = P2p(copula$min_corr_matrix)
      t_cop    <- tCopula(param=P,dim=7,dispstr = "un", df = copula$df)
      U        <- rCopula(lenBlock,t_cop)
    }
    # take the off-diagonal elements of the correlation matrix estimated 
    # under ML method with 12 degrees of fredom i.e. best t-copula under AIC
    # see table 9, chapter 3
    else if(type.dependence=="t_ml"){
      P = P2p(fitted_parameter$t_fitted_corr_ml)
      t_cop    <- tCopula(param=P,dim=7,dispstr = "un", df = 12)
      U        <- rCopula(lenBlock,t_cop)
  }
    else if(type.dependence=="clayton"){
      alpha = fitted_parameter$clayton_fitted_alpha_ml
      clayton_cop <- claytonCopula(param=alpha,dim=7)
      U           <- rCopula(lenBlock,clayton_cop)
    }
  else if(type.dependence=="frank"){
    alpha = fitted_parameter$frank_fitted_alpha_ml
    frank_cop <- frankCopula(param=alpha,dim=7)
    U           <- rCopula(lenBlock,frank_cop)
  }
  else if(type.dependence=="gumbel"){
    alpha = fitted_parameter$gumbel_fitted_alpha_ml
    gumbel_cop <- gumbelCopula(param=alpha,dim=7)
    U           <- rCopula(lenBlock,gumbel_cop)
  }
  else if(type.dependence=="uncorr"){
      U        <- runif(lenBlock*7)
      U        <- matrix(U,ncol = 7,nrow=lenBlock)
    }
    #else if(type.dependence=="perfcorr"){
      #U        <- runif(lenBlock)
      #U        <- matrix(rep(U,7),ncol = 7,nrow=lenBlock)
    #}
  return(U)
}


capital_charge <- function(margin,nSim,nBlock,type.dependence){
  # Sort the aggregate simulated losses for the 7 ETs
  margin <-lapply(margin, function(et){
    lapply(et, function(i){
      i <- sort(i,decreasing = FALSE, na.last = FALSE)
    })
  })

   lenBlock <-  nSim/nBlock
   VaR_U    <- vector()
  # Extract a (0,1) values from a copula or uniform distribution and 
  # save it in a list where each element is a matrix with 7 columns
  # and "lenBlock" rows
  # i.e. the matrix has on the rows values in (0,1)
  # that are able to capture the correlation among ETs
   U        <- lapply(1:nBlock,function(sim){U(lenBlock,type.dependence)})
   # Inizialize the empty list where each element is a matrix that
   # will take the simulated and diversified losses i.e. the simulated
   # losses obtained by the function "margin" multiplied by the element in U
   losses   <-  vector(mode = "list", length = nBlock)
   losses   <- lapply(losses, function(mat){mat<-matrix(ncol = 7,nrow = lenBlock)})
   # First loop: take the number of block
    # Second loop: take the number of columns of the matrixs in the list i.e. the ETs
     # Third lopp: take the number of the simulated losses in each block
   for(sim in 1:nBlock){
      for(et in 1:7){
        for(i in 1:lenBlock){
       # Note that margin[[et]][[i]] is a vector (we called its elements S^i_j - see chapter 4.1),
       # therefore I pick up the elemet
       # "ceiling(lenBlock * U[i,et])" of this sort vector
        losses[[sim]][i,et] <- margin[[et]][[sim]][ceiling(lenBlock * U[[sim]][i,et])]
      }
    }
    # Sum the columns of the i-th row of "losses[[sim]]" matrix, in such a way the aggregate 
    # annual and diversified simulated loss is obtained
    # In other words we are performing the sum_{j=1}^7 S^i_j
    aggreg_loss <- apply( losses[[sim]],1,sum)
    # Take the quantile 0.999 of such aggregate and diversified simulated losses
    VaR_U[sim] <- quantile(aggreg_loss,0.999,na.rm=TRUE)
   }
  # By central limit theory the final VaR will be just the average of the VaRs
  # computed from each block
  VaR <- mean(VaR_U)
  return(list(VaR_U,VaR))
}



type_copula <- c("t","t_ml","clayton","frank","gumbel","uncorr")

# compute the capital charge under various copulas
capital_charge_all_copulas <- lapply(type_copula, function(i){
  capital_charge(sim_margin,nSim=nSim,nBlock = nBlock,i)})

names(capital_charge_all_copulas) <- type_copula

capital_charge_result <- list(capital_charge_perf_corr = capital_charge_all_copulas)

#save(capital_charge_result,file="capital_charge_result.RData")


br <- seq(min(capital_charge[[1]]),max(capital_charge[[1]]),length.out = 50)
hist(capital_charge[[1]],breaks = br)