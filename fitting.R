fit_QQ <- function(dati,k,xi,beta,n){

        x_bd    <- vector()
        x_log   <- vector()
        x_norm  <- vector()
        x_exp   <- vector() 
        sim_bd  <- vector()
        sim_ed  <- vector()
        sim_log <- vector()
        sim_norm<- vector()
        sim_exp <- vector()
        U       <- runif(n)
        # run "n" pseudo observation from 
        # body-tail distribution i.e. empirical-GPD distribution
        for(i in 1:n){
          x_bd[i]    <- body_tail(dati,k,xi,beta,n,U[i])
        }
        # run "n" pseudo observation from lognormal distribution
         log_par <- MASS::fitdistr(dati,"lognormal")
        for(i in 1:n){
          x_log[i]    <- stats::qlnorm(U[i],log_par$estimate[[1]],log_par$estimate[[2]])  
        }
         # run "n" pseudo observation from normal distribution
        norm_par <- MASS::fitdistr(dati,"normal")
         for(i in 1:n){
           x_norm[i]    <- stats::qnorm(U[i],norm_par$estimate[[1]],norm_par$estimate[[2]])  
         }
        # run "n" pseudo observation from exponential distribution
        exp_par <- MASS::fitdistr(dati,"exponential")
        for(i in 1:n){
          x_exp[i]    <- stats::qexp(U[i],exp_par$estimate[[1]])  
        }
        # compute the empirical cumulative dfs functions
        p_ed       <- stats::ecdf(dati)
        p_bd       <- stats::ecdf(x_bd)
        p_log      <- stats::ecdf(x_log)
        p_norm      <- stats::ecdf(x_norm)
        p_exp      <- stats::ecdf(x_exp)
           for(i in 1:n){
             # simulate from empirical cumulative dfs functions
             sim_ed [i]    <- p_ed(dati[ceiling(length(dati)*U[i])])
             sim_bd [i]    <- p_bd(dati[ceiling(length(dati)*U[i])])
             sim_log[i]    <- p_log(dati[ceiling(length(dati)*U[i])])
             sim_norm[i]   <- p_norm(dati[ceiling(length(dati)*U[i])])
             sim_exp[i]    <- p_exp(dati[ceiling(length(dati)*U[i])])
           }
          return(list(sim_ed=sim_ed,sim_bd=sim_bd,sim_log=sim_log,
                      sim_norm=sim_norm,sim_exp=sim_exp))
}
body_tail <- function(dati,k,xi,beta,n,U){
    sim_loss   <- vector()
    #Number of events less than the trheshold out all events
    F_u <- length(dati[dati<k])/length(dati)
    n_coda <- 0
        for(j in 1:length(U)){
          if(U[j] < F_u){
            U[j] <- (U[j]/F_u)
            sim_loss[j] <- dati[ceiling(length(dati)* U[j])]
          }
          else{
            U[j] <- ((U[j]-F_u)/(1-F_u))
            sim_loss[j]    <- k + qgpd(U[j], xi = xi, beta = beta)
          }
        }
    return(sim_loss)
    }
