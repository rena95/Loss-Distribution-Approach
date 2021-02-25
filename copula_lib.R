load("dati_pseudo.RData")
# Data are in the following format i.e. LONG FORMAT
#  CODICE ETLEV1   PTL        DATA
# 1	80478	1	426572.79	2014-08-21
# 2	728215	3	68087.18	2013-02-01
# 3	862286	4	347551.44	2008-05-07
# 4	4679476	1	25270.05	2011-04-18
# 5	1247	7	52112.14	2011-05-19
#   .............................
# Monthly aggregation of the losses
# mlosses return a data.frame (n X k) where k are the ETs
mlosses <- function(dati){
  library("lubridate")
  dati <- dati %>% dplyr::mutate(DATA_RIL = ymd(paste0(year(DATA_RIL),
  	"-",month(DATA_RIL),"-","1")))
  
  dati <- dati %>% dplyr::group_by(ETLEV1,DATA_RIL) %>% 
  dplyr::summarise(PTL = sum(PTL))

  dati <- dati %>% tidyr::spread(ETLEV1,PTL) %>% 
  dplyr::arrange(DATA_RIL)

  dati[is.na(dati)]   <- 0
  colnames(dati)[2:8] <- paste0("ET",1:7)
  return(dati <- as.data.frame(dati))
}
# mlosses transform data uploaded from load("dati_pseudo.RData") 
# in the following format i.e. SHORT FORMAT
#  DATA_RIL        ET1       ET2         ET3          ET4           ET5         ET6       ET7
# 1	2007-01-01	1580852.61	13125504    5115556.9	 4610128	291992.90	 33255.07	 3379720
# 2	2007-02-01	47430.63	4700840	    1168138.7	 1901531	48842.28	 56339.82	 1895047
# 3	2007-03-01	562063.16	11638957	4199101.9	 79446932	1018598.83   89886.24	 5952838
# .................................................................................


#####-----Estimation of t-copula parameters-----#####
#Lindskog estimator 
r_star   <- function(tau){
              sin(pi/2*tau)
            }

#Embrechts estimator 
lambda_u <- function(r_star,v){
              2*pt(-sqrt((v+1)*((1-r_star)/(1+r_star))),v+1)
            }
#Euclidean distance
euclidean_distance <- function(p,q){
  sqrt(sum((p - q)^2))
}

# Empirical tail copulae
# Reference "Nonparametric estimation of tail dependence"
# SCHMIDT,STADTMULLER (2006)
emp_tail_copulae <- function(dati,b, min.lambda=20){
  # Note that "dati" must be a data.frame in short format!
  # Number of ET
  nET  <- ncol(dati)
  # Number of losses i.e. number of monthly aggregated losses
  # we have 12 years, therefore there are at most 132 observations for each ET
  n  <- nrow(dati)
  # Rank value of the data matrix
  rank <- apply(dati, 2, rank)
  # Begin Empirical tail copulae
  # We build the estimator at page 7 of Schmidt, Stadtmuller- equation 13
  lambda_u_matrix <- matrix(nrow =nET,ncol = nET )
  tmp             <- vector()
  for(i in 1:(nET-1)){
    for(j in (i+1):nET){
      # m is defined at pag 15 Frahm, Junker 
      # when this function is called "b" is set equal 10
      m  <- round(sqrt(n-2*b))
      for(s in 1:n){
        tmp[s] <- sum(rank[,i]>n-s & rank[,j]>n-s)/s
      }
      # this is the plateu-finding algorithm as defined by
      # Frahm, Junker,Schmidt - Estimating the tail dependence 
      # coefficient: properties and pitfall - pag 15
      smooth_tmp <- vector()
      for(l in 1:(n-b)){
      		# smooth tmp by the box kernel then
      		# the empirical tail copulae estimator is induced by 
      		# homogeneity
      		 smooth_tmp[l] <- sum(tmp[l:(l+b)])/b
      }
      smooth_tmp.sd   <- 2*sd(smooth_tmp)
      n.m <- n-m
      lambda_hat <- vector()
      for(s in 1:n.m){
        #p_k
        lambda_hat[s] <- sum(abs(smooth_tmp[s] - smooth_tmp[(s+1):(s+m-1)])) 
      }
      TF <- lambda_hat <= smooth_tmp.sd
      if (any(TF)){min <-max(min.lambda, which(TF)[1])
      lambda_u_matrix[i,j] <- lambda_u_matrix[j,i] <- lambda_u <- mean(smooth_tmp[min:(min+m-1)])
	  }
    }
  }
  lambda_u_matrix[is.na(lambda_u_matrix)] <- 0
  diag(lambda_u_matrix)                   <- 1
  return(lambda_u_matrix)
}



objective_df <- function(v, data, P){
	-sum(dcopula.t(data, v, P, log = TRUE))
}


# This function is used to estimate the correlation matrix R* 
# and the degrees of fredom for the t-copula 
# as explained in the example 3 - section 3.7
copula <- function(dati,b,type.delta, min.lambda=20, guess_df = 12){
    library(QRM)
	cop                                <- list()
	cop$monthly                        <- dati
	n.classi      					   <- ncol(cop$monthly[2:8])
	cop$monthly [is.na(cop$monthly  )] <- 0
	
	# Emprical cumulative distribution function for each ET
	cop$ecd  <- apply(cop$monthly[2:8], 2 ,edf, adjust = TRUE)
	
	#compute r_star con Lindskog estimator
	tau           <- Kendall(cop$monthly [2:8])
	r_star        <- r_star(tau)
	#Number of Monte Carlo simulation
	nMC                 <- 10000
	#compute the d_f given r_star - Note QRM::fit.tcopula() can also be used
	minim <- nlminb(guess_df,objective_df, data = cop$ecd , P = r_star )
	v     <- round(minim[[1]])
	# Step 1 - Example 3 - Chapter 3
	# Perform Monte Carlo simulation creating "nMC" random Matrixes from:
	# 1) a normal with correlation matrix=cop$dev_st  
	# 2) from uniform random variables in [0,b], here  b=0.5
	
	if (type.delta=="cov"){
	#set up the correlation matrix from which extract the random matrix
	  cop$dev_st       <- matrix(0.4,nrow=7,ncol=7)
	  diag(cop$dev_st) <- 0
	  cop$perturb      <- array(rnorm(nMC*n.classi*n.classi, mean=0.0,
                                 sd=cop$dev_st),
                                dim=c(n.classi,n.classi,nMC))
	  #Force the simulated matrix to be symmetric					   
	  cop$perturb      <- array(apply(cop$perturb ,3, function(P) {
	  	as.matrix(Matrix::forceSymmetric(P))}),
	  dim=c(n.classi,n.classi,nMC))
	}
	
	else{
	  #Example 3 - Step 1: generate nMC unif[0,b] random variables
	  unif_b        <-  runif(nMC*n.classi*n.classi,0,0.5)
	  cop$perturb   <-  array(unif_b,dim=c(n.classi,n.classi,nMC))
	  cop$perturb    <- array(apply(cop$perturb  ,3, function(P) {
	  	as.matrix(Matrix::forceSymmetric(P))}),
	                          dim=c(n.classi,n.classi,nMC))
	  # set the diagonal element equal 0
	  for (i in 1:nMC) {diag(cop$perturb[,,i]) <- 0}
	  #diag(cop$perturb) <- 0
	}
	#emp_tail_copulae 
	cop$emp <- emp_tail_copulae(cop$monthly[2:8],b)
	#Create the array 7 X 7 X N of montecarlo perturbation
	matr.new.tail     <- array(dim=c(n.classi,n.classi,nMC))
	distance          <- vector()
	matr.new.corr.min <- matrix(nrow = 7, ncol = 7)
	
	# New tau matrix after first perturbation
	matr.new        <- pmax(tau+cop$perturb[,,1],tau)
	
	#compute the positive linear correlation matrix r_star given the perturbation
	# to get positive matrix (using Matrix::nearPD), this is the eigenvalue alghorithm described 
	# in chapter 3 to define positive the matrix
	matr.new.corr <- pmin(as.matrix(Matrix::nearPD(r_star(matr.new))$mat),1)
	
	#compute the Embrechts tail matrix
	matr.new.tail[,,1] <- lambda_u(matr.new.corr,v)
	
	#compute the euclidean distance
	distance[1]   <- distance.min <- euclidean_distance(matr.new.tail[,,1],cop$emp)
	
	#Try the matrix for which the distance between cop$emp (emp_tail_copulae)
	# and matr.new.tail is minimum
	for (i in 2:dim(cop$perturb)[3]){
	  #take the greatest number between r_star and r_star+pert
		matr.new      <- pmax(tau+cop$perturb[,,i],tau)
		#compute the positive linear correlation matrix r_star given the perturbation
		matr.new.corr <- pmin(as.matrix(Matrix::nearPD(r_star(matr.new))$mat),1)
		#compute the Embrechts tail matrix - Step 4 - Example 3
		matr.new.tail[,,i] <- lambda_u(matr.new.corr,v)
		distance[i] <- euclidean_distance(matr.new.tail[,,i],cop$emp)
		# Step 6 and 7 - Example 3
		if (distance[i] < distance.min){
			distance.min           <- distance[i]
			matr.new.corr.min      <- matr.new.corr
		}
	
			}
			#Select matr.new.corr related to matr.new.tail.min
			return(list(cop$monthly  ,cop$perturb, cop$emp, distance , matr.new.corr.min,v))

}




