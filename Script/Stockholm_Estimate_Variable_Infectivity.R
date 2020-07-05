


# Remove workspace
rm(list = ls(all=TRUE))

heading <- "
#------------------------------------------------------------------------------------
#
# FILE:         Stockholm_Estimate_Variable_Infectivity
# CONTACT:      Lisa Brouwers
# EMAIL:        analysenheten@folkhalsomyndigheten.se
# AFFILIATION:  FD-AN
# PROJECT:      Modelling COVID-19
#
#
# Created:      2020-03-15 
# Updated:      2020-06-22
# R version:    3.5.2
#
# What the script does: This script estimates the parameters of the infectivity described in the report 
#                       'Estimates of the the number of infected individuals during the covid-19 outbreak 
#                       in the Stockholm region, Västra Götaland region, Skåne region, and Dalarna region, Sweden.'
#                       The main focus of this script is the region of Stockholm, and it produces the main results for
#                       this region. 
#                       
#                       
#                   
#                   
#------------------------------------------------------------------------------------
\n"





#-------------#
# Disposition #
#-------------#

# 1. Some important notes on notations
# 2. Paths (You need to create your own, see section)
# 3. Libraries
# 4. Read incidence data per region
# 5. Read Population data
# 6. Choice of parameter values and region
# 7. Visualise reported cases/incidence for the four regions 
# 8. Basic functions and set.up
# 9. Main estimation function
# 10. Estimation
# 11. Initial inspection of results
# 12. Calculate results to save in tables. 
#     - Save estimated parameters, their SE, and CI 
#     - Bootsrap CIs
#     - Save estimated R0, Re, infectivity and their CI  
#     - Save estimated peak days and their CI  
#  13. Plot figures
#  14. Analysing increased contacts




#-----------------------------------------------------------------------------------
# 1. Some important notes on notations 
#-----------------------------------------------------------------------------------


#### Compartments ###
##
## S := Susceptible 
## E := Exposed, infected but not infectious. 
## I_symp := Infected case that becomes reported
## I_asymp := Infected case that not becomes reported, i.e. unreported case.
## R_1 := First recovered/immune compartment. Individuals in this compartment are still assumed to test positive on PCR-tests for SARS-CoV-2.
## R_2 := Second recovered/immune compartment. Individuals in this compartnebt do not test positive on PCR-tests.
##
####

#### Rates and parameters specified ###
##
## p_lower_inf := how much lower infectivity unreported cases have. In this analysis consistently set to 1 (equal infectivity).
## t_b := 16 march 2020. Day Swedes were recommended to work from home.
## eta := rate of leaving incubation period. Denoted rho in report.
## gamma_D := rate of leaving the infectious compartment, i.e. 1/gamma_D = mean time infectious
## gamma_pos :=  rate of leaving the first recovered compartment.
##
####

#### Parameters estimated  ###
##
## p_symp := fraction reported cases. Denoted p_r in report. Sometimes denoted pb in script. 
## p_asymp = 1 - p_symp := fraction unreported cases. Denoted p_u in report.
##
## The time dependent infectivity beta takes three parameters except time t,
## beta = beta(t, theta, p, epsilon). 
## 
## IMPORTANT
## The parameter p is denoted delta in the report!
## 
####

#### Misc
##
## Gloria 2 
## Study that measured the prevalence of active or near time infection (by PCR-test) 
## in Stockholm between 17th March and 3rd April. 
## Estimate = 2.5% had infection.
##
## Gloria 3
## Same kind of study as Gloria 2 but between 21st April and 24th April
## Estimate = 2.3% had infection.
##
####


#-----------------------------------------------------------------------------------
# 2. Paths 
#-----------------------------------------------------------------------------------


# Input your own path instead of "YOUR_PATH" where folder is located in project.path (end with "/").
# e.g. project.path 	     <- "C:/Users/Modelling/"
#
# Also, create a folder within your main folder named "Data" where you save the four txt-files :
# Data_Dalarna_no_pers_aldre_prim_2020-06-05.txt
# Data_skane_no_pers_aldre_prim_2020-06-05.txt
# Data_Stockholm_no_pers_aldre_prim_2020-06-05.txt
# Data_VGR_no_pers_aldre_prim_2020-06-05.txt
#
# Create a folder within the main folder named Results, with two folders named Figures and Tables
#


project.path  <- YOUR_PATH

# 
data.path     <- paste(project.path, "Data", sep="") 
figure.path   <- paste(project.path, "Results/Figures", sep="") 
table.path    <- paste(project.path, "Results/Tables", sep="") 



#-----------------------------------------------------------------------------------
# 3. Libraries
#-----------------------------------------------------------------------------------

library(reshape2)
library(readxl)
library(writexl)
library(openxlsx)     # to write tables in excel
library(RColorBrewer)
require(rootSolve)    # to load function multiroot that finds roots to system of equations
library(Rcpp)
require(deSolve)
library(RColorBrewer)


library("coda")
library("adaptMCMC")


#-----------------------------------------------------------------------------------
# 4. Read incidence data per region
#-----------------------------------------------------------------------------------



Stockholm_Data  <- read.table(file=paste(data.path, "/Data_Stockholm_no_pers_aldre_prim_2020-06-05.txt",sep=""), sep = " ", header=TRUE)
Dalarna_Data    <- read.table(file=paste(data.path, "/Data_Dalarna_no_pers_aldre_prim_2020-06-05.txt",sep=""), sep = " ", header=TRUE)
VGR_Data        <- read.table(file=paste(data.path, "/Data_VGR_no_pers_aldre_prim_2020-06-05.txt",sep=""), sep = " ", header=TRUE)
Skane_Data      <- read.table(file=paste(data.path, "/Data_skane_no_pers_aldre_prim_2020-06-05.txt",sep=""), sep = " ", header=TRUE)


#-----------------------------------------------------------------------------------
# 5. Read Population data
#-----------------------------------------------------------------------------------

N_Stockholm       <- 2374550
N_Dalarna         <- 287795
N_Skane           <- 1376659
N_Vastra_Gotaland <- 1724529

N_region <- c("Stockholm" = N_Stockholm, 
              "Dalarna" = N_Dalarna,
              "Skane" = N_Skane,
              "Vastra Gotaland" = N_Vastra_Gotaland
)


#-----------------------------------------------------------------------------------
# 6. Choice of parameter values and region
#-----------------------------------------------------------------------------------



## Region

REGION     <- "Stockholm"
Fixed_data <- Stockholm_Data


## Parameters

# eta1 is our chosen value for eta, the mean time in the latency stage, the exposed compartment.
# gammaD1 is our chosen value for gammaD, the mean time being infectious is set to 5 days

gammaD1 <- 1/5 
eta1 <- 1/5.1


Days_pos <- 5 # Number of days, beyond the infectivity period, someone test positivte on PCR-test
# Since infectivity period 5 days, assigning Days_pos = 5 => 10 days total testing positive

gammaPos <- 1/Days_pos # inverse of how many days as recovered you test positive


SEED_W_ONE <- FALSE # The initial values of number of infected. Should be false if we seed with unreported cases as well.
p_lower_inf_use <- 1 # Factor of the infectivity of unreported cases in comparison to reported cases. Here they are equal.



#-----------------------------------------------------------------------------------
# 7. Visualise reported cases/incidence for the four regions 
#-----------------------------------------------------------------------------------


dayatyear_feb_july <- c(32, 61, 92, 122, 153, 153+30)
NameDateFebJuly <- as.Date(dayatyear_feb_july,origin ="2019-12-31")
NameDateFebJuly <- c("Feb", "Mar", "Apr", "May", "Jun", "Jul")


x_Stockholm <- as.numeric(as.Date(Stockholm_Data$Datum)) - as.numeric(as.Date("2019-12-31"))
x_Dalarna   <- as.numeric(as.Date(Dalarna_Data$Datum)) - as.numeric(as.Date("2019-12-31"))
x_VGR       <- as.numeric(as.Date(VGR_Data$Datum)) - as.numeric(as.Date("2019-12-31"))
x_Skane     <- as.numeric(as.Date(Skane_Data$Datum)) - as.numeric(as.Date("2019-12-31"))


sum(Stockholm_Data$Incidens)
100 * sum(Stockholm_Data$Incidens)/N_Stockholm

sum(VGR_Data$Incidens)
100 * sum(VGR_Data$Incidens)/N_Vastra_Gotaland

sum(Skane_Data$Incidens)
100 * sum(Skane_Data$Incidens)/N_Skane

sum(Dalarna_Data$Incidens)
100 * sum(Dalarna_Data$Incidens)/N_Dalarna



par(mfrow = c(2,2), mar = c(6.1, 4.1, 4.1, 2.1)) #

# Dalarna

  plot(x_Dalarna, Dalarna_Data$Incidens, type = "o", pch = 16,
       ylab = "Reported cases", xlab ="",xlim = c(range(dayatyear_feb_july)),
       ylim = c(0, max(Dalarna_Data$Incidens) + 5),
       xaxt ='n', main = "Daily number of reported cases: Dalarna")

  grid(nx=NA, ny=NULL)
  abline(v = dayatyear_feb_july, col = "lightgray", lty = "dotted", lwd = par("lwd"))
  axis(side = 1, at = dayatyear_feb_july, label = NameDateFebJuly,las=2)
  lines(x_Dalarna, Dalarna_Data$Incidens, type = "o", pch = 16)

# Skåne

  plot(x_Skane, Skane_Data$Incidens, type = "o", pch = 16,
       ylab = "Reported cases", xlab ="",xlim = c(range(dayatyear_feb_july)),
       ylim = c(0, max(Skane_Data$Incidens) + 5),
       xaxt ='n', main = "Daily number of reported cases: Skåne")

  grid(nx = NA, ny = NULL)
  abline(v = dayatyear_feb_july, col = "lightgray", lty = "dotted", lwd = par("lwd"))

  axis(side = 1, at = dayatyear_feb_july, label = NameDateFebJuly, las = 2)
  lines(x_Skane, Skane_Data$Incidens, type = "o", pch = 16)


# Stockholm

  plot(x_Stockholm, Stockholm_Data$Incidens, type = "o", pch = 16,
       ylab = "Reported cases", xlab = "", xlim = c(range(dayatyear_feb_july)),
       ylim = c(0, max(Stockholm_Data$Incidens) + 5),
       xaxt ='n', main = "Daily number of reported cases: Stockholm")

  grid(nx = NA, ny = NULL)
  abline(v = dayatyear_feb_july, col = "lightgray", lty = "dotted", lwd = par("lwd"))

  axis(side = 1, at = dayatyear_feb_july, label = NameDateFebJuly, las = 2)
  lines(x_Stockholm, Stockholm_Data$Incidens, type = "o", pch = 16)

# VGR

  plot(x_VGR, VGR_Data$Incidens, type = "o", pch = 16,
       ylab = "Reported cases", xlab = "", xlim = c(range(dayatyear_feb_july)),
       ylim = c(0, max(VGR_Data$Incidens) + 10),
       xaxt ='n', main = "Daily number of reported cases: VGR")

  grid(nx = NA, ny = NULL)
  abline(v = dayatyear_feb_july, col = "lightgray", lty = "dotted", lwd = par("lwd"))

  axis(side = 1, at = dayatyear_feb_july, label = NameDateFebJuly,las = 2)
  lines(x_VGR, VGR_Data$Incidens, type = "o", pch = 16)


  
  
#-----------------------------------------------------------------------------------
# 8. Basic functions and set.up
#-----------------------------------------------------------------------------------


# Rounding functions for plotting
roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

roundUp <- function(x) 10^ceiling(log10(x))

# Colours for plotting
Blues <- brewer.pal(n = 9, name = "Blues")


# Function for quantilies
CRI_90_low <- function(x){
  return(quantile(x, probs = 0.9))
}
CRI_90_up <- function(x){
  return(quantile(x, probs = 0.9))
}

CRI_95_low <- function(x){
  return(quantile(x, probs = 0.025))
}
CRI_95_up <- function(x){
  return(quantile(x, probs = 0.975))
}


# Functions for transformation of parameters
logit <- function(p) {
  log(p) - log(1 - p)
}


expit <- function(x) { 
  1 / (1 + exp(-x)) 
} 





#------------------#
#    C++ Model     #
#------------------#

# read the model and ODE solver
source(paste(project.path,"Script/SEIR_model_Two_R_C++.R",sep=""))



#-----------------------------------------------------------------------------------
# 9. Main estimation function
#-----------------------------------------------------------------------------------




Atol = 1e-8
Rtol = 1e-10


# This is the main function estimating the parameters of the time-dependent infectivity 
# Its most important outputs 
#   - Opt_par: The estimated parameters, p, theta, epsilon, and p_symp
#     p, theta, and epsilon are the parameters governing the time-dependent infectivity and 
#     p_symp the fraction of new cases that becomes reported. 
#
#   - SEIR_model: the model that calculates the number of individuals in each compartments at different 
#                 time points
#
#   - Initial_values: The estimated initial values for SEIR_model, if seed_w_one = FALSE. 
#
# The parameter p_symp can only be estimated based on Stockholm data, since we from this region 
# have the results from Gloria 2 and Gloria 3
#



# In optimisation, iter is how many iterations/guesses



Estimate_function_Stockholm <- function(Region = "Stockholm", 
                                        p_lower_inf = 1, 
                                        gamma_pos = gammaPos, 
                                        gammaD = gammaD1, 
                                        eta = eta1, 
                                        iter = 100,
                                        seed_w_one = FALSE){
  
  if(gamma_pos == Inf){gamma_pos <- 0}

  N <- N_Stockholm

  Incidence <- Stockholm_Data$Incidens
  Datum     <- as.Date(Stockholm_Data$Datum)
  Day       <- as.numeric(Datum) - as.numeric(as.Date("2019-12-31"))
  

  dayatmonth <- c(1,31,29,31,30,31,30,31,31,30,31,30,31)
  dayatyear <- cumsum(dayatmonth)
  Namedate <- as.Date(dayatyear, origin = "2019-12-31")
  
  
  # p is delta in the report.
  # p_b is p_r in the report. b stands for the swedish word "bekräftad" = reported/verified
  Opt_par_names_t <- c("logit_p", "epsilon", "log_theta", "logit_p_b")  

  GuessesT <- function(){ # guesses for the optimisation
    
    #transformed 
    u_p <- logit(runif(1, 0.05, 0.6)) #
    u_e <- runif(1,-0.6, 0)    #
    u_t <- log(runif(1, 0, 10))    #
    if(p_lower_inf >= 0.5){ u_t <- log(runif(1, 0, 2)) }
    if(p_lower_inf >= 0.8){ u_t <- log(runif(1, 0, 1)) }
    u_pb <- logit(runif(1, 0, 0.06)) # prob reported
    return(c(u_p, u_e, u_t, u_pb))
    
  }
  

  t_today <- as.numeric(as.Date("2020-03-16")) - as.numeric(as.Date("2019-12-31")) # day Swedes were recommended to work from home
  # time-dependent infectivity beta 
  beta.peak.free <- function(t,  p, epsilon,theta){
    res <- ((1 - p)/(1 + exp(epsilon*(-(t - t_today)))) + p)*theta
    return(res)
  }
  
  # Time-dependent basic reproduction number
  Basic_repr <- function(t,p,epsilon,theta,gamma,p_symp){ 
    res <- p_symp * beta.peak.free(t,p,epsilon,theta)/gamma + (1 - p_symp)*p_lower_inf*beta.peak.free(t,p,epsilon,theta)/gamma
    return(res)
    
  }
  
  
  model <- function(time, state, parameters) {
    ## The vector state contains:
    ## S       <- state[1] # susceptibles
    ## E       <- state[2] # latent/exposed but not infectious
    ## I_symp  <- state[3] # infected who get reported
    ## I_asymp <- state[4] # infected who remain non-reported
    ## R1       <- state[5] # recovered/immune but testing positive
    ## R2       <- state[6] # recovered/immune but NOT testing positive
    par <- as.list(c(state, parameters))
    with(par, {
      dS <- -beta.peak.free(time,p,epsilon,theta) * S * I_symp/N - p_lower_inf*beta.peak.free(time,p,epsilon,theta) * S * I_asymp/N
      dE <-  beta.peak.free(time,p,epsilon,theta) * S * I_symp/N + p_lower_inf*beta.peak.free(time,p,epsilon,theta) * S * I_asymp/N - eta*E
      dI_symp <- p_symp * eta * E       - gammaD * I_symp
      dI_asymp <- (1 - p_symp)* eta * E - gammaD * I_asymp
      dR1 <- gammaD * (I_symp + I_asymp) - gamma_pos * R1
      dR2 <- gamma_pos * R1
      dx <- c(dS, dE, dI_symp, dI_asymp, dR1, dR2)
      list(dx)
    }
    )
  }
   
  
  RSS <- function(parameters) {
    ## The vector parameters contains the 4 parameters to estimate:
    ## parameters[1] = p
    ## parameters[2] = epsilon 
    ## parameters[3] = theta 
    ## parameters[4] = p_symp/pb 
    
    if(seed_w_one == TRUE){
      init <- c(S = N - Incidence[1], 
                E = 0, 
                I_symp = Incidence[1], 
                I_asymp = 0, 
                R1 = 0, 
                R2 = 0)
    }else{
      init <- c(S = N - Incidence[1]*(1 + (1-expit(parameters[4]))/expit(parameters[4])), 
                E = 0, 
                I_symp = Incidence[1], 
                I_asymp = Incidence[1]*(1-expit(parameters[4]))/expit(parameters[4]), 
                R1 = 0, 
                R2 = 0)
    }
    
    
    ## Create vector named parms to be sent to model function SEIRmodel2
    parms = c(p = expit(parameters[1]), epsilon = parameters[2], theta = exp(parameters[3]),
              t_today = t_today, 
              N = N, 
              p_lower_inf = p_lower_inf, 
              eta = eta, 
              p_symp = expit(parameters[4]), 
              gamma_pos = gamma_pos, 
              gammaD = gammaD,
              t_start = Day[1],
              t_end = Day[length(Day)])
    
    Dummy_infectivity <- beta.peak.free(t = c(0: 700), p = parms[1], epsilon = parms[2],  theta = parms[3])
    ## if the infectivity is negative, throw away guess
    if(min(Dummy_infectivity) < 0 ){
      res <- 10^12
      #print("negative infectivity")
      return(res)
    }else{
      
      #out <- ode(y = init, times = Day, func = SEIRmodel, parms = parms, atol = 1e-8, rtol = 1e-10)
      out <- SEIRmodel2(state = init, parms = parms)
      
      fit_S <- out[ , 2]
      fit_E <- out[ , 3]
      fit_I_symp <- out[ , 4]
      fit_I_asymp <- out[ , 5]
      
      fit_R1 <- out[ , 6]
      fit_R2 <- out[ , 7]
      fit_R <- fit_R1 + fit_R2
      
      
      ## Model's number of individuals who are able to test positive on PCR test
      if(gamma_pos == 0){ # if only infectious individuals test positive
        Infectious    <- fit_I_symp + fit_I_asymp 
        InfectiousF   <-  Infectious[40:47] # Gloria 2
        InfectiousFG3 <-  Infectious[65:68] # Gloria 3
        
      }else{
        # including the R1 in those who test positive
        Infectious    <- fit_I_symp + fit_I_asymp + fit_R1
        InfectiousF   <-  Infectious[40:47] # Gloria 2
        InfectiousFG3 <-  Infectious[65:68] # Gloria 3
        
      }
      
      
      fitted_incidence  <- expit(parameters[4]) * fit_E[1:length(Incidence)] * eta
      logMLEcases <- - (length(Incidence)) * log(sum((Incidence - fitted_incidence)^2)) / 2 
      
      
      binom_est     <- round(mean(InfectiousF))
      binom_estG3   <- round(mean(InfectiousFG3))
      
      #gloria 2 days [40:47]
      logMLE_SthlmWeek1 <- dbinom(binom_est, size = N, prob = 0.025, log = TRUE)
      
      # gloria 3, days[65:68]
      logMLE_SthlmWeek2 <- dbinom(binom_estG3, size = N, prob = 0.023, log = TRUE)
      
      #return(-logMLEcases - logMLE_SthlmWeek1)
      return(-logMLEcases - logMLE_SthlmWeek1 - logMLE_SthlmWeek2)
      
      
    }
  }
  
  
  
  #print("Optimisation initialised")
  Guess <- GuessesT()
  Opt <- 0
  #conl <- list(maxit = 1000)
  conl <- list(maxit = 1000, abstol = 1e-8, reltol = 1e-10) 
  
  Opt <- optim(Guess, RSS, control = list(conl), hessian = TRUE)
   
  
  while(Opt$convergence>0){
    Guess <- GuessesT()
    Opt <- optim(Guess, RSS, control = list(conl), hessian = TRUE)  
  }
   
  
  
  i = 1 
  while(i <= iter){
    Guess <- GuessesT()
    
    conl <- list(maxit = 1000, parscale = c(1, 1, 0.01, 10 ), abstol = 1e-8, reltol = 1e-10) 
    
    Opt2 <- optim(Guess, RSS, control = list(conl), hessian = TRUE)
    
    if((Opt2$convergence == 0)){
      
      # if( ((i/iter)*100) %% 10 == 0){
      # How_much_done <- (i/iter)*100
      # print(paste(How_much_done,"% done", sep=""))
      # }
      
      i = i + 1
      if(Opt2$value < Opt$value){Opt <- Opt2}
    }
    
  }
  
  Opt_par_transformed <- Opt$par 
  names(Opt_par_transformed) <- Opt_par_names_t 
  
  Opt_par <- c(p = unname(expit(Opt_par_transformed["logit_p"])), 
               epsilon = unname(Opt_par_transformed["epsilon"]), 
               theta = unname(exp(Opt_par_transformed["log_theta"])),
               p_symp = unname(expit(Opt_par_transformed["logit_p_b"]))
  )
  
  if(seed_w_one==TRUE){
    init <- c(S = N - Incidence[1], 
              E = 0, 
              I_symp = Incidence[1], 
              I_asymp = 0, 
              R1 = 0, 
              R2 = 0)
  }else{
    init <- c(S = N - Incidence[1]*(1 + (1-unname(Opt_par[4]))/unname(Opt_par[4])), 
              E = 0, 
              I_symp = Incidence[1], 
              I_asymp = Incidence[1] * (1-unname(Opt_par[4]))/unname(Opt_par[4]), 
              R1 = 0, 
              R2 = 0)
  }
  
  return(list(Observed_incidence = Incidence, Population_size = N, Day = Day, dayatyear = dayatyear, 
              Namedate = Namedate, Optimisation = Opt, Opt_par = Opt_par, Seasonal_forcing = beta.peak.free, 
              Basic_reproduction = Basic_repr, Initial_values = init, SEIR_model = model))
  
  
  
}


#-----------------------------------------------------------------------------------
# 10. Estimation
#-----------------------------------------------------------------------------------


gammaD <- gammaD1
eta <- eta1


set.seed(483958913)

  
Est_par_model <- Estimate_function_Stockholm(Region = REGION, 
                                             iter = 150, 
                                             gamma_pos = gammaPos, 
                                             p_lower_inf = p_lower_inf_use, 
                                             seed_w_one = SEED_W_ONE)



# Days of incidence based on studied region
Day <- Est_par_model$Day
Last_date_used <- as.Date(Day[length(Day)], origin = "2019-12-31")


N   <- Est_par_model$Population_size
dayatyear <- Est_par_model$dayatyear
Namedate  <- Est_par_model$Namedate 


#Observed incidence based on studied region
Observed_incidence <-  Est_par_model$Observed_incidence 
Est <- Est_par_model$Optimisation


options(digits=3)

## Estimated parameters.
Opt_par <- Est_par_model$Opt_par


# functions based on model scenario
Basic_repr <- Est_par_model$Basic_reproduction
beta       <- Est_par_model$Seasonal_forcing
SEIR_model <- Est_par_model$SEIR_model

# initial values based on model scenario
init       <- Est_par_model$Initial_values



# if logLikelihood
LogL_MLE <- -Est$value
RSS_value <- exp(-(2/length(Observed_incidence)) * LogL_MLE)


H <- Est$hessian
sdParams <- sqrt(diag(solve(H)))

options("scipen"=100, "digits"=4)
#default options("scipen"=0, "digits"=7)

  

#-----------------------------------------------------------------------------------
# 11. Initial inspection of results
#-----------------------------------------------------------------------------------


#----------------------------------------------------------
# Look at the results. 
# The observed vs estimated nr reported cases. 
# And estimated infectivity and basic reproductive number
#----------------------------------------------------------
  
  
  
  t <- (Day[1]):(Day[length(Day)]+223) # time in days
  fit <- data.frame(ode(y = init, times = t, func = SEIR_model , parms = Opt_par))
  fit_S <- fit[ , 2]
  fit_E <- fit[ , 3]
  fit_I_symp <- fit[ , 4]
  fit_I_asymp <- fit[ , 5]
  fit_R1 <- fit[ , 6]
  fit_R2 <- fit[ , 7]
  
  fit_R <- fit_R1 + fit_R2
  fit_I <- fit_I_symp + fit_I_asymp
  fit_cum_inf <- N - fit_S
  fit_antibodies <- fit_cum_inf - fit_E
  
  fitted_incidence  <- Opt_par[4] * fit_E * eta
  fitted_incidence_non_report  <- (1 - Opt_par[4]) * fit_E * eta
  Incidence_both <- fit_E * eta
  
  
  # tot 
  fit_cum_inf[length(fit_cum_inf)]/N
  
  # Gloria 2, point-estimate 2.5%
  
  if(Days_pos == 0){
    Infectious  <- fit_I_symp + fit_I_asymp #+ fit_E
    InfectiousF <-  Infectious[40:47]
    mean(InfectiousF/N)
  }else{
    # including the R1 in those who test positive
    Infectious <- fit_I_symp + fit_I_asymp + fit_R1
    InfectiousF   <-  Infectious[40:47]
    mean(InfectiousF/N)
  }
    
  

  



# Look at the observed vs estimated nr reported cases. And estimated infectivity and basic reproductive number

  
dayatyear_feb_jan <- c(32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336, 367)
NameDateFebJan <- as.Date(dayatyear_feb_jan,origin ="2019-12-31")
  
  

NameNumber <- paste("/",REGION,"_", Last_date_used,"_",
                    "Days_pos_", Days_pos + 5,
                    "_Reported_cases.pdf", sep="")

## Plot the reported and fitted cases, uncomment pdf and dev.off to save fig.

#pdf(paste(figure.path, NameNumber, sep=""), width=9, height=5 )
  
  par(mfrow = c(1,1), mar = c(6.1, 4.1, 4.1, 2.1)) # 
  
  plot(Day,Observed_incidence,type="o", ylab = "Reported cases", xlab ="", xlim = c(t[1],366), xaxt ='n', 
       main = paste("Estimated and simulated number of reported cases: ", REGION, sep= ""))
  #grid
  grid(nx = NA, ny = NULL)
  abline(v = dayatyear_feb_jan, col = "lightgray", lty = "dotted", lwd = par("lwd"))
  #
  
  lines(t,fitted_incidence,  lwd = 2, col = "blue")
  axis(side = 1, at = dayatyear_feb_jan, label = NameDateFebJan, las = 2)
  legend("topright", c("Observed reported cases", "Fitted reported cases"), lwd=c(1,2), pch = c(1,NA), lty = c(1,1), col = c("black","blue"))

#dev.off()






t_Re   <- c(Day[1]:200)
Re     <- Basic_repr(t_Re, p = Opt_par[1], epsilon = Opt_par[2],  theta = Opt_par[3], gamma = gammaD, p_symp = Opt_par[4]) * fit_S[t_Re] / N
infect <- beta(t_Re, p = Opt_par[1], epsilon = Opt_par[2],  theta = Opt_par[3])

contact_reduced = infect[1]/infect[which(t_Re == Day[length(Day)])]



NameNumber <- paste("/",REGION, "_", Last_date_used,
                    "_Re_Infectivity_", 
                    "Days_pos_", Days_pos + 5, 
                    ".pdf",sep="")

  
#pdf(paste(figure.path, NameNumber, sep=""), width=9*1.25, height=4.5*1.3 )
  
  par(mfrow = c(1,2), mar = c(6.1, 4.1, 6.1, 5.1)) # 
  plot(t_Re, Re, 
       ylim=c(0, roundUpNice(max(Re))),
       type="l", ylab="R(t)",lwd=2, 
       main = paste("Estimated reproductive number \n ", REGION, sep = "") , 
       xlab="", xaxt ='n')
  
  grid(nx=NA, ny=NULL)
  abline(v=dayatyear_feb_jan,col = "lightgray", lty = "dotted", lwd = par("lwd"))
  
  lines(t_Re, Re, lwd = 2)
  abline(v=Day[length(Day)], lty = 2)
  axis(side = 1, at = dayatyear_feb_jan, label = NameDateFebJan, las = 2)
  
  
  plot(t_Re, infect,type="l", ylab="Infectivity",lwd=2, 
       ylim = c(0, roundUpNice(max(infect))),
       main = paste("Estimated infectivity \n ", REGION, sep = "") , 
       xlab = "", xaxt='n')
  
  grid(nx=NA, ny=NULL)
  abline(v = dayatyear_feb_jan, col = "lightgray", lty = "dotted", lwd = par("lwd"))
  lines(t_Re, infect, lwd = 2)
  
  abline(v=Day[length(Day)], lty = 2)
  
  axis(side = 1, at = dayatyear_feb_jan, label = NameDateFebJan, las = 2)

#dev.off()


  


#-----------------------------------------------------------------------------------
# 12. Calculate results to save in tables. 
#     - Save estimated parameters, their SE, and CI 
#     - Bootsrap CIs
#     - Save estimated R0, Re, infectivity and their CI  
#     - Save estimated peak days and their CI  
#-----------------------------------------------------------------------------------

  


CI_level_05 <- 0.025

## CIs of estimated parameters
p_high        <- expit(qnorm(1-CI_level_05, mean = Est$par[1], sd = sdParams[1], lower.tail = TRUE, log.p = FALSE))
epsilon_high  <- qnorm(1-CI_level_05, mean = Est$par[2], sd = sdParams[2], lower.tail = TRUE, log.p = FALSE)
theta_high    <- exp(qnorm(1-CI_level_05, mean = Est$par[3], sd = sdParams[3], lower.tail = TRUE, log.p = FALSE))
pb_high       <- expit(qnorm(1-CI_level_05, mean = Est$par[4], sd = sdParams[4], lower.tail = TRUE, log.p = FALSE))

p_low         <- expit(qnorm(CI_level_05, mean = Est$par[1], sd = sdParams[1], lower.tail = TRUE, log.p = FALSE))
epsilon_low   <- qnorm(CI_level_05, mean = Est$par[2], sd = sdParams[2], lower.tail = TRUE, log.p = FALSE)
theta_low     <- exp(qnorm(CI_level_05, mean = Est$par[3], sd = sdParams[3], lower.tail = TRUE, log.p = FALSE))
pb_low        <- expit(qnorm(CI_level_05, mean = Est$par[4], sd = sdParams[4], lower.tail = TRUE, log.p = FALSE))

## For bootstrap CI's
p.v        <- expit(rnorm(1000, mean = Est$par[1], sd = sdParams[1]))
epsilon.v  <- rnorm(1000,mean = Est$par[2], sd = sdParams[2])
theta.v    <- exp(rnorm(1000,mean = Est$par[3], sd = sdParams[3]))
pb.v       <- expit(rnorm(1000, mean = Est$par[4], sd = sdParams[4]))


R0.v.Dag1 <- c()
R0.v.DagSista <- c()
for(i in 1:length(p.v)){
  
  R0.v.Dag1[i]      <- Basic_repr(Day[1],p = p.v[i], epsilon = epsilon.v[i], theta = theta.v[i], gamma = gammaD, p_symp = pb.v[i]) 
  R0.v.DagSista[i]  <- Basic_repr(Day[length(Day)],p = p.v[i], epsilon = epsilon.v[i], theta = theta.v[i], gamma = gammaD, p_symp = pb.v[i]) 
}



R0_low  <- c(CRI_95_low(R0.v.Dag1), CRI_95_low(R0.v.DagSista))
R0_high <- c(CRI_95_up(R0.v.Dag1), CRI_95_up(R0.v.DagSista))

R0_Mean <- Basic_repr(Day, p = Opt_par[1], epsilon = Opt_par[2],  theta = Opt_par[3], gamma = gammaD, p_symp = Opt_par[4])



#############################################
## save estimated parameters and their SE  ##
#############################################


res_param <- c(round(c(5 + Days_pos),digits=4), 
               round(mean(InfectiousF/N), digits = 5), 
               round(c(LogL_MLE, RSS_value, Opt_par[1], sdParams[1], Opt_par[2], sdParams[2], Opt_par[3], sdParams[3], 1-Opt_par[4], sdParams[4]), digits = 4))

names(res_param) <- c("Days pos","27 mars - 3 april", "LogL", "RSS" ,"delta", "logit s.e.", "epsilon", "s.e.", "theta", "log s.e.",  "p_o", "logit s.e.")
CIp       <- paste("[",round(p_low,digits=3), ", ",round(p_high,digits=3),"]",sep="")
CIepsilon <- paste("[",round(epsilon_low,digits=3),", ",round(epsilon_high,digits=3),"]", sep="")
CItheta   <- paste("[",round(theta_low,digits=3),", ", round(theta_high,digits=3),"]", sep="")
CIp_o     <- paste("[",round(1-pb_high,digits=4),", ", round(1-pb_low,digits=4),"]", sep="")

CI_param <- c("", "", "", "", CIp, "", CIepsilon, "", CItheta, "", CIp_o, "")

MAT_para <- matrix(c(res_param,CI_param), ncol = 12, nrow = 2, byrow = TRUE)

df.res <- as.data.frame(MAT_para)
colnames(df.res) <- names(res_param)


XL_file_name <- paste(table.path,"/",REGION, "_", Last_date_used, 
                        "_Days_pos_", Days_pos + 5,"_parameters.xlsx", sep ="") 


write.xlsx(df.res, XL_file_name)




##########################
## Create bootstrap CIs ##
##########################

t <- (Day[1]):(Day[length(Day)]+4*31) # time in days
fit         <- data.frame(ode(y = init, times = t, func = SEIR_model , parms = Opt_par))
fit_S       <- fit[ , 2]
fit_E       <- fit[ , 3]
fit_I_symp  <- fit[ , 4]
fit_I_asymp <- fit[ , 5]
fit_R       <- fit[ , 6]
fit_I       <- fit_I_symp + fit_I_asymp
fit_I_E     <- fit_E + fit_I
fit_cum_inf <- N - fit_S
fit_antibodies <- fit_cum_inf - fit_E



fitted_incidence  <- Opt_par[4] * fit_E * eta
fitted_incidence_non_report  <- (1 - Opt_par[4]) * fit_E * eta

fitted_incidence_tot <- fit_E * eta






Cum_Inf_inc <- data.frame(Cumulative_infected = fit_cum_inf, 
                          Incidence_reported = fitted_incidence, 
                          Incidence_non_reported = fitted_incidence_non_report,
                          Incidence_tot = fitted_incidence_tot)

fit_Cum_Inf_inc <- cbind(Date = as.Date(t, origin = "2019-12-31"), fit, Cum_Inf_inc)


file_name_X <- paste(table.path,"/",REGION,"_",Last_date_used,"_",
                     "Days_pos_", Days_pos + 5, 
                     "_Fitted_SEIR_and_Incidence_with_time.xlsx", sep ="")


## if you want to save complete time-series of number in each compartment + incidence.
## then uncomment the row below. 

# write.xlsx(fit_Cum_Inf_inc, file_name_X )




fit_S.v       <- matrix(rep(NA, times = length(t) * length(p.v)), ncol = length(t))
fit_E.v       <- matrix(rep(NA, times = length(t) * length(p.v)), ncol = length(t))
fit_I_symp.v  <- matrix(rep(NA, times = length(t) * length(p.v)), ncol = length(t))
fit_I_asymp.v <- matrix(rep(NA, times = length(t) * length(p.v)), ncol = length(t))
Fit_I.v       <- matrix(rep(NA, times = length(t) * length(p.v)), ncol = length(t))
fit_cum_inf.v <- matrix(rep(NA, times = length(t) * length(p.v)), ncol = length(t))

fitted_incidence.v            <- matrix(rep(NA, times = length(t) * length(p.v)), ncol = length(t))
fitted_incidence_nonreport.v  <- matrix(rep(NA, times = length(t) * length(p.v)), ncol = length(t))
fitted_incidence_tot.v        <- matrix(rep(NA, times = length(t) * length(p.v)), ncol = length(t))
effective_reprod.v            <- matrix(rep(NA, times = length(t) * length(p.v)), ncol = length(t))



for(i in 1:length(p.v)){
  Opt_parDummy = c(p = p.v[i], epsilon = epsilon.v[i], theta = theta.v[i],  p_symp = pb.v[i])
  
  fitDummy <- data.frame(ode(y = init, times = t, func = SEIR_model , parms = Opt_parDummy))
  fit_S.v[i,] <- fitDummy[ , 2]
  fit_E.v[i,] <- fitDummy[ , 3]
  fit_I_symp.v[i,] <-  fitDummy[ , 4]
  fit_I_asymp.v[i,] <- fitDummy[ , 5]
  Fit_I.v[i,] <-  fitDummy[ , 4] + fitDummy[ , 5]
  fit_cum_inf.v[i,] <- N - fitDummy[ , 2]
  
  fitted_incidence.v[i,]            <- Opt_par[4] *  fitDummy[ , 3] * eta
  fitted_incidence_nonreport.v[i,]  <- (1 - Opt_par[4]) *  fitDummy[ , 3] * eta
  fitted_incidence_tot.v[i,]        <- fitDummy[ , 3] * eta
 
  effective_reprod.v[i,] <- Basic_repr(fitDummy[,1], p = Opt_parDummy[1], epsilon = Opt_parDummy[2], theta = Opt_parDummy[3], gamma = gammaD, p_symp = Opt_parDummy[4]) * fitDummy[ , 2] /N
}




cum_inf_mean       <- apply(fit_cum_inf.v, MARGIN = 2, FUN = mean) 
cum_inf_median     <- apply(fit_cum_inf.v, MARGIN = 2, FUN = median) 
cum_inf_95_up_CRI  <- apply(fit_cum_inf.v, MARGIN = 2, FUN = CRI_95_up) 
cum_inf_95_low_CRI <- apply(fit_cum_inf.v, MARGIN = 2, FUN = CRI_95_low) 


fit_antibodies_mean       <- apply(fit_cum_inf.v - fit_E.v, MARGIN = 2, FUN = mean) 
fit_antibodies_median     <- apply(fit_cum_inf.v - fit_E.v, MARGIN = 2, FUN = median) 
fit_antibodies_95_up_CRI  <- apply(fit_cum_inf.v - fit_E.v, MARGIN = 2, FUN = CRI_95_up) 
fit_antibodies_95_low_CRI <- apply(fit_cum_inf.v - fit_E.v, MARGIN = 2, FUN = CRI_95_low) 


fit_I_mean       <- apply(Fit_I.v, MARGIN = 2, FUN = mean) 
fit_I_median     <- apply(Fit_I.v, MARGIN = 2, FUN = median) 
fit_I_95_up_CRI  <- apply(Fit_I.v, MARGIN = 2, FUN = CRI_95_up) 
fit_I_95_low_CRI <- apply(Fit_I.v, MARGIN = 2, FUN = CRI_95_low) 

fitted_Incidence_mean       <- apply(fitted_incidence.v, MARGIN = 2, FUN = mean) 
fitted_Incidence_median     <- apply(fitted_incidence.v, MARGIN = 2, FUN = median) 
fitted_Incidence_95_up_CRI  <- apply(fitted_incidence.v, MARGIN = 2, FUN = CRI_95_up) 
fitted_Incidence_95_low_CRI <- apply(fitted_incidence.v, MARGIN = 2, FUN = CRI_95_low) 

fitted_Incidence_tot_mean       <- apply(fitted_incidence_tot.v, MARGIN = 2, FUN = mean) 
fitted_Incidence_tot_median     <- apply(fitted_incidence_tot.v, MARGIN = 2, FUN = median) 
fitted_Incidence_tot_95_up_CRI  <- apply(fitted_incidence_tot.v, MARGIN = 2, FUN = CRI_95_up) 
fitted_Incidence_tot_95_low_CRI <- apply(fitted_incidence_tot.v, MARGIN = 2, FUN = CRI_95_low) 



effective_reprod_mean       <- apply(effective_reprod.v, MARGIN = 2, FUN = mean) 
effective_reprod_median     <- apply(effective_reprod.v, MARGIN = 2, FUN = median) 
effective_reprod_95_up_CRI  <- apply(effective_reprod.v, MARGIN = 2, FUN = CRI_95_up) 
effective_reprod_95_low_CRI <- apply(effective_reprod.v, MARGIN = 2, FUN = CRI_95_low) 






maxdagen <- as.Date(t[which(fit_I==max(fit_I))],origin = "2019-12-31" )

minDag <- as.Date(t[which(fit_I_95_low_CRI==max(fit_I_95_low_CRI))],origin = "2019-12-31" )
maxDag <- as.Date(t[which(fit_I_95_up_CRI==max(fit_I_95_up_CRI))],origin = "2019-12-31" )

maxdagen_incidens <- as.Date(t[which(fitted_incidence_tot == max(fitted_incidence_tot))], origin = "2019-12-31" )
minDagIncidens    <- as.Date(t[which(fitted_Incidence_tot_95_low_CRI == max(fitted_Incidence_tot_95_low_CRI))], origin = "2019-12-31" )
maxDagIncidens    <- as.Date(t[which(fitted_Incidence_tot_95_up_CRI == max(fitted_Incidence_tot_95_up_CRI))], origin = "2019-12-31" )







######################################################
## save estimated R0, Re, infectivity and their CI  ##
######################################################


# save estimated R0 and their uncertainty

res_R0        <- round(c(infect[1], infect[which(t_Re == Day[length(Day)])] , contact_reduced, 
                         infect[which(t_Re == 162)], 
                         R0_Mean[1], R0_Mean[length(Day)],effective_reprod_mean[length(Day)]), digits = 3)

names(res_R0) <- c("infectivity(start)","infectivity(June 5)","infectivity reduced",  "infectivity 10/6", "R0(start)", "R0(June 5)", "Re(June 5)")

CI_R01        <- paste("[",round(R0_low[1],digits=3), ", ", round(R0_high[1],digits=3),"]", sep="")
CI_ROend      <- paste("[",round(R0_low[2],digits=3), ", ", round(R0_high[2],digits=3),"]", sep="")
CI_Reend      <- paste("[",round(effective_reprod_95_low_CRI[length(Day)],digits=3), ", ", round(effective_reprod_95_up_CRI[length(Day)],digits=3),"]",sep="")
CIR0          <- c("","","","", CI_R01, CI_ROend,CI_Reend)

MAT_R0    <- matrix(c(res_R0,CIR0), ncol = 7, nrow = 2,byrow=TRUE)
df.resR0  <- as.data.frame(MAT_R0)
colnames(df.resR0) <- names(res_R0)


XL_R0_file_name <- paste(table.path,"/",REGION, "_", Last_date_used, "_",
                           "Days_pos_", Days_pos + 5, 
                           "_R0_Re_Infectivity", ".xlsx", sep ="")





write.xlsx(df.resR0, XL_R0_file_name )










###########################################
## save estimated peak days and their CI ##
###########################################




# 2020-04-11 = day 102
# 2020-05-01 = day 122
# 2020-06-01 = day 153
# 2020-07-01 = day 183
# 2020-08-01 = day 245
# 2020-09-01 = day 245



fit_cum_low   <- cum_inf_95_low_CRI
fit_cum_high  <- cum_inf_95_up_CRI


res_days <-c(round(fit_cum_inf[which(t==122)],digits=0), round(fit_cum_inf[which(t==122)]/N,digits=3),  
             round(fit_cum_inf[which(t==153)],digits=0), round(fit_cum_inf[which(t==153)]/N,digits=3),  
             round(fit_cum_inf[which(t==183)],digits=0), round(fit_cum_inf[which(t==183)]/N,digits=3),  
             round(fit_cum_inf[which(t==245)],digits=0), round(fit_cum_inf[which(t==245)]/N,digits=3),  
             round(fit_cum_inf[which(t==245)],digits=0), round(fit_cum_inf[which(t==245)]/N,digits=3))

res_maxprevalens <- round(max(fit_I),digits=0)
res_days <- c(res_days, as.character(maxdagen), res_maxprevalens)

res_maxincidens <- round(max(fitted_incidence_tot), digits = 0)
res_days <- c(res_days, as.character(maxdagen_incidens),res_maxincidens)


names(res_days)   <- c("2020-05-01", "", "2020-06-01","","2020-07-01","","2020-08-01","","2020-09-01","", 
                       "Peak-day prevalence", "Prevalence at peak-day", "peak-day incidence", "Incidence at peak-day")



CIdag1maj         <- paste("[",round(fit_cum_low[which(t==122)],digits=0),", ",round(fit_cum_high[which(t==122)],digits=0),"]",sep="")
CIdag1majProc     <- paste("[",round(fit_cum_low[which(t==122)]/N,digits=3),", ",round(fit_cum_high[which(t==122)]/N,digits=3),"]",sep="")

CIdag1juni        <- paste("[",round(fit_cum_low[which(t==153)],digits=0),", ",round(fit_cum_high[which(t==153)],digits=0),"]",sep="")
CIdag1juniProc    <- paste("[",round(fit_cum_low[which(t==153)]/N,digits=3),", ",round(fit_cum_high[which(t==153)]/N,digits=3),"]",sep="")

CIdag1juli        <- paste("[",round(fit_cum_low[which(t==183)],digits=0),", ",round(fit_cum_high[which(t==183)],digits=0),"]",sep="")
CIdag1juliProc    <- paste("[",round(fit_cum_low[which(t==183)]/N,digits=3),", ",round(fit_cum_high[which(t==183)]/N,digits=3),"]",sep="")

CIdag1aug         <- paste("[",round(fit_cum_low[which(t==214)],digits=0),", ",round(fit_cum_high[which(t==214)],digits=0),"]",sep="")
CIdag1augProc     <- paste("[",round(fit_cum_low[which(t==214)]/N,digits=3),", ",round(fit_cum_high[which(t==214)]/N,digits=3),"]",sep="")

CIdag1sep         <- paste("[",round(fit_cum_low[which(t==245)],digits=0),", ",round(fit_cum_high[which(t==245)],digits=0),"]",sep="")
CIdag1epProc      <- paste("[",round(fit_cum_low[which(t==245)]/N,digits=3),", ",round(fit_cum_high[which(t==245)]/N,digits=3),"]",sep="")

CImaxdag          <- paste("[",as.character(as.Date(minDag, origin = "2019-12-31")),", ",as.character(as.Date(maxDag, origin = "2019-12-31")),"]",sep="")
CIprevalensMaxdag <- paste("[",round(max(fit_I_95_low_CRI),digits=0),", ",round(max(fit_I_95_up_CRI),digits=0),"]",sep="")
CIincidems_maxdag <- paste("[",as.character(as.Date(minDagIncidens, origin = "2019-12-31")),", ",as.character(as.Date(maxDagIncidens, origin = "2019-12-31")),"]",sep="")
CIincidensMaxdag  <- paste("[",round(max(fitted_Incidence_tot_95_low_CRI),digits=0),", ",round(max(fitted_Incidence_tot_95_up_CRI),digits=0),"]",sep="")



CI_dag <- c(CIdag1maj,CIdag1majProc, 
            CIdag1juni, CIdag1juniProc, 
            CIdag1juli, CIdag1juliProc, 
            CIdag1aug, CIdag1augProc,
            CIdag1sep, CIdag1epProc,
            CImaxdag, CIprevalensMaxdag, 
            CIincidems_maxdag, CIincidensMaxdag)

MAT_dag <- matrix(c(res_days,CI_dag),ncol = 14 ,nrow=2,byrow=TRUE)

df.dag <- as.data.frame(MAT_dag)
colnames(df.dag) <- names(res_days)




XL_file_name <- paste(table.path,"/", REGION, "_", Last_date_used,"_",
                        "Days_pos_", Days_pos + 5,
                        "_Estimated_peak_days", ".xlsx", sep ="")





write.xlsx(df.dag, XL_file_name )





#-----------------------------------------------------------------------------------
#  13. Plot figures
#-----------------------------------------------------------------------------------

## Figure with I + E to compare with antibody studies. One then needs to look 2-3 weeks before
## week 14 - week 22


Week_14_to_22_Mondays <- c(43, 50, 57, 64, 71, 78, 85, 92, 99)
W14_to_W22 <- 43:99

Week_14_to_22_Wednesdays <- Week_14_to_22_Mondays + 2
as.Date(t[Week_14_to_22_Mondays], origin = "2019-12-31")

fit_antibodies[Week_14_to_22_Mondays]/N
fit_antibodies_95_low_CRI[Week_14_to_22_Mondays]/N
fit_antibodies_95_up_CRI[Week_14_to_22_Mondays]/N





NameDate_Week_14_to_22_Mondays <- paste0("W",c(14:22))


NameNumber <- paste("/", REGION, "_", Last_date_used,
                      "_Days_pos_", Days_pos + 5,
                    "_AntiB_CI", ".pdf",sep="")




CEX = 1.5
CEX.AXIS = 1.3


pdf(paste(figure.path, NameNumber, sep = ""), width = 7.5, height = 6 )

  
  par(mfrow = c(1,1), mar = c(6.3, 4.1, 4.1, 2.1)) # 
  # bottom, left, top, and right. The default is c(5.1, 4.1, 4.1, 2.1).
  
  MainTitle <- paste("Fitted SEIR model covid-19 ", REGION ,
                     "\n Fraction possible to undergo seroconversion",sep="")
  
  # 2 weeks for sero to appear
  plot(W14_to_W22, fit_antibodies[W14_to_W22]/N, type = "l", col = "black", ylab = "I + R", lwd = 2,
       xlab = "", ylim = c(0, roundUpNice(max(fit_antibodies_95_up_CRI[W14_to_W22]/N))), xaxt = 'n', main = MainTitle)
  
  
  grid(nx=NA, ny=NULL)
  
  COL = as.vector(col2rgb(Blues[4]))
  COL = rgb(red = COL[1], green = COL[2], blue = COL[3], max = 255, alpha = 150)
  
  abline(v = Week_14_to_22_Mondays, col = "lightgray", lty = "dotted", lwd = par("lwd"))
  
  polygon(x = c(seq_along(W14_to_W22) + W14_to_W22[1]-1, rev(seq_along(W14_to_W22)+ W14_to_W22[1]-1)),
          y = c(fit_antibodies_95_low_CRI[W14_to_W22]/N, rev(fit_antibodies_95_up_CRI[W14_to_W22]/N)),
          lty = 0, col = COL)
  
  lines(W14_to_W22, fit_antibodies[W14_to_W22]/N, col = "black", lwd = 2)
  
  axis(side = 1, at = Week_14_to_22_Mondays, label = NameDate_Week_14_to_22_Mondays, las = 2)
  
  
dev.off()








## ~4 months ahead in time

t <- (Day[1]):(Day[length(Day)] + 4*31) # time in days
fit <- data.frame(ode(y = init, times = t, func = SEIR_model, parms = Opt_par))
fit_S <- fit[ , 2]
fit_E <- fit[ , 3]
fit_I_symp <- fit[ , 4]
fit_I_asymp <- fit[ , 5]
fit_I <- fit_I_symp + fit_I_asymp
fit_I_E <- fit_E + fit_I
fit_cum_inf <- N - fit_S



fitted_incidence  <- Opt_par[4] * fit_E * eta

fitted_incidence_low <- fitted_Incidence_95_low_CRI[1:length(t)]
fitted_incidence_high <- fitted_Incidence_95_up_CRI[1:length(t)]

fitted_I_high <- fit_I_95_up_CRI[1:length(t)]
fitted_I_low <- fit_I_95_low_CRI[1:length(t)]









# 1 march, 1 april, ..., 1 october
dayatyear_march_september <- c(61, 61 + 31, 
                               61 + 31 + 30, 
                               61 + 31 + 30 + 31,
                               61 + 31 + 30 + 31 + 30, 
                               61 + 31 + 30 + 31 + 30 + 31,
                               61 + 31 + 30 + 31 + 30 + 31 + 31,
                               61 + 31 + 30 + 31 + 30 + 31 + 31 + 30)
NameDateMarchSeptember <- as.Date(dayatyear_march_september, origin = "2019-12-31")





NameNumber <- paste("/", REGION,"_", Last_date_used, "_",
                    "Days_pos_", Days_pos + 5, 
                    "_Infected_recovered_CI.pdf", sep = "")




pdf(paste(figure.path, NameNumber, sep=""), width=7.5, height=6 )
  
  par(mfrow = c(1,1), mar = c(6.3, 4.1, 4.1, 2.1)) # 
  # bottom, left, top, and right. The default is c(5.1, 4.1, 4.1, 2.1).
  
  MainTitle <- paste("Fitted SEIR model covid-19 ", REGION ,
                     "\n Cumulative fraction of infected: E + I + R",sep="")
  
  
  # 2 weeks for sero to appear
  plot(fit$time, fit_cum_inf[1:length(t)]/N, type="l", col="black", ylab = "E + I + R" ,lwd=2,
       xlab="",
       ylim = c(0,0.3),
       #ylim = c(0,roundUpNice(max(fit_cum_high[1:length(t)]/N))), 
       xaxt='n', main = MainTitle)
  
  grid(nx=NA, ny=NULL)
  
  COL = as.vector(col2rgb(Blues[4]))
  COL = rgb(red = COL[1], green = COL[2], blue = COL[3], max = 255, alpha = 150)
  
  abline(v=dayatyear_march_september, col = "lightgray", lty = "dotted", lwd = par("lwd"))
  
  polygon(x = c(seq_along(fit$time) + fit$time[1]-1, rev(seq_along(fit$time) + fit$time[1] - 1)),
          y = c(fit_cum_low[1:length(t)]/N, rev(fit_cum_high[1:length(t)]/N)),
          lty = 0, col = COL)
  
  lines(fit$time, fit_cum_inf[1:length(t)]/N, col = "black", lwd = 2)
  
  axis(side = 1, at = dayatyear_march_september, label = NameDateMarchSeptember, las = 2)



dev.off()







NameNumber <- paste("/",REGION,"_",Last_date_used, "_",
                    "Days_pos_", Days_pos + 5,
                    "_Reported_and_total_prevalence_CI.pdf",sep="")


pdf(paste(figure.path, NameNumber, sep=""), width=9*1.5, height=4.5*1.3 )
  
  
  par(mfrow = c(1,2), mar = c(6.1, 4.1, 6.1, 5.1)) # 
  
  ## Prognosis until May
  
  plot(fit$time, fitted_incidence, type="l", col="red", ylab = "Incidence", main = "" ,lwd=2,
       xlab="",ylim=c(0,roundUpNice(max(fitted_incidence_high,Observed_incidence))),xaxt='n')
  
  
  grid(nx=NA, ny=NULL)
  
  
  abline(v=dayatyear_march_september,col = "lightgray", lty = "dotted", lwd = par("lwd"))
  
  lines(fit$time, fitted_incidence,col="red",lwd=2)
  points(Day, Observed_incidence)
  
  lines(fit$time, fitted_incidence_low, lty=2)
  lines(fit$time, fitted_incidence_high, lty=2)
  
  axis(side = 1, at = dayatyear_march_september, label = NameDateMarchSeptember,las=2)
  
  
  ## next plot

  plot(fit$time, fit[ , 4]+fit[ , 5], ylim = c(0, roundUpNice(max(fitted_I_high)) ) , type="l", col="red", 
       ylab = "Number of infectious",lwd=2,
       xlab="",xaxt='n')
  
  grid(nx=NA, ny=NULL)
  abline(v=dayatyear_march_september,col = "lightgray", lty = "dotted", lwd = par("lwd"))
  lines(fit$time, fit[ , 4]+fit[ , 5],col="red",lwd=2)
  
  lines(fit$time, fitted_I_high, lty=2)
  lines(fit$time, fitted_I_low, lty=2)
  
  
  axis(side = 1, at = dayatyear_march_september, label = NameDateMarchSeptember,las=2)
  
  
  MainTitle <- paste("Fitted SEIR model covid-19: ", REGION, sep="")
  
  title(MainTitle, outer = TRUE, line = -2.5)
  
  
  title("Estimated and observed number of daily reported cases", outer = TRUE, line = -5, adj = 0.085)
  title("Number of infectious at the same time", outer = TRUE, line = -5, adj = 0.825)

dev.off()









#-----------------------------------------------------------------------------------
# 14. Analysing increased contacts
#-----------------------------------------------------------------------------------

# In the main analysis we used case data until 5th June and estimated the parameters of the 
# time-dependent infectivity. Here, we allow this estimated infectivity to increase after 9th June 
# (first day of increase is 10th June). The increase continues (linearly) until 31th August 
# when it reaches its maximum, the infectivity then stays at this maximal value.



gamma_pos   <- 1/Days_pos # inverse of how many days as recovered you test positive
if(Days_pos == 0){gamma_pos = 0}


p           <- unname(Opt_par[1])
epsilon     <- unname(Opt_par[2])
theta       <- unname(Opt_par[3])
p_symp      <- unname(Opt_par[4])



## No contact increase as base 


t <- c(Day[1]:365)
fit <- data.frame(ode(y = init, times = t, func = SEIR_model , parms = Opt_par))
fit_S <- fit[ , 2]
fit_E <- fit[ , 3]
fit_I_symp <- fit[ , 4]
fit_I_asymp <- fit[ , 5]
fit_I <- fit_I_symp + fit_I_asymp
fit_I_E <- fit_E + fit_I
fit_cum_inf <- N - fit_S



fitted_incidence_reported_same  <- p_symp * fit_E * eta
fitted_incidence_non_reported_same  <- (1 - p_symp) * fit_E * eta
Incidence_both_same <- fit_E * eta




## Contact increase


## initiate lists to save incidence results in
fitted_incidence_reported_inc_contact <- list()
fitted_incidence_non_reported_inc_contact <- list()
Incidence_both_inc_contact <- list()





## Contact increase from Day_rise to 31th August 
## We start the contact increase from 10th June = day 162 
## This means that the first day with an actual increase is 10th June. 
## E.g. infectivity 8th June = 0.1932, 9th June = 0.1932, 10 juni 0.1937



Day_rise <- 162


base_infectivity <- unname(beta(Day_rise, p = Opt_par[1], epsilon = Opt_par[2],  theta = Opt_par[3]))

JanJune <- 1:(Day_rise - 2)

June_to_max_slow <- (Day_rise-1):244 # to 31th August 
Max_slow_dec     <- 245:366 


min_cont <- 1 
# choose which factors to increase the contacts/infectivity with
# here, 20%, 40%, 60%, 80%, and 100%
Increase_contact <- seq(1.2, 2, by = 0.2)


for(i in 1:length(Increase_contact)){
  
  max_cont <- Increase_contact[i]
  
  Val_JanJune   <- rep(min_cont, length(JanJune))
  
  Val_JuneToMax_Slow <- seq(min_cont, max_cont, by = (max_cont-min_cont)/(length(June_to_max_slow)-1))
  Val_MaxToDec_Slow  <- rep(max_cont, length(Max_slow_dec))
  

  Values_slow_rise <- c(Val_JanJune, Val_JuneToMax_Slow, Val_MaxToDec_Slow)
  Values <- Values_slow_rise
   
  
  
  beta_G <- function(t){ 
    
    res <- base_infectivity * Values[t]

    return(res) 
  }
  
  
  t_today <- as.numeric(as.Date("2020-03-16")) - as.numeric(as.Date("2019-12-31")) # dagen för jobba hemma
  
  beta.peak.free <- function(t,p, epsilon, theta){
    
    res <- (t < Day_rise)*((1-p)/(1+exp(epsilon*(-(t-t_today)))) + p)* theta + 
      (t>=Day_rise)*((1-p)/(1+exp(epsilon*(-(Day_rise-t_today)))) + p)* theta * Values[t]
    
    return(res)
  }
  
  

  
  # rate of going from latency to infectious eta
  seir.model.asymptomatics <- function(time, state, parameters) {
    # S       <- state[1] # susceptibles
    # E       <- state[2] # latent/exposed but not infectious
    # I_symp  <- state[3] # infected who get reported
    # I_asymp <- state[4] # infected who remain non-reported
    # R1       <- state[5] # recovered/immune but testing positive
    # R2       <- state[6] # recovered/immune but NOT testing positive
    par <- as.list(c(state, parameters))
    with(par, {
      dS <- -beta.peak.free(t = time, p = p, epsilon = epsilon, theta = theta) * S * I_symp/N - p_lower_inf_use*beta.peak.free(t = time, p = p, epsilon = epsilon, theta = theta) * S * I_asymp/N
      dE <-  beta.peak.free(t = time, p = p, epsilon = epsilon, theta = theta)  * S * I_symp/N + p_lower_inf_use*beta.peak.free(t = time, p=p, epsilon = epsilon, theta = theta)  * S * I_asymp/N - eta*E
      dI_symp <- p_symp * eta * E       - gammaD * I_symp
      dI_asymp <- (1 - p_symp)* eta * E - gammaD * I_asymp
      dR1 <- gammaD * (I_symp + I_asymp) - gamma_pos * R1
      dR2 <- gamma_pos * R1
      dx <- c(dS, dE, dI_symp, dI_asymp, dR1, dR2)
      list(dx)
    }
    )
  }
  
  
  
  
  if(SEED_W_ONE==TRUE){
    init <- c(S = N - Incidence[1] ,
              E = 0, 
              I_symp = Incidence[1], 
              I_asymp = 0, 
              R1 = 0, 
              R2 = 0)
  }else{
    init <- c(S = N - Observed_incidence[1]*(1 + (1-p_symp)/p_symp) , 
              E = 0, 
              I_symp = Observed_incidence[1], 
              I_asymp = Observed_incidence[1]*(1-p_symp)/p_symp , 
              R1 = 0, 
              R2 = 0)
  }
  
  
  
  fit <- data.frame(ode(y = init, times = t, func = seir.model.asymptomatics , parms = Opt_par))
  
  
  fit_E <- fit[ , 3]
  fitted_incidence_reported_inc_contact[[i]]      <- p_symp * fit_E * eta
  fitted_incidence_non_reported_inc_contact[[i]]  <- (1 - p_symp) * fit_E * eta
  Incidence_both_inc_contact[[i]] <- fit_E * eta
  
 
  
  Incidence_df <- data.frame(time = fit[,1], 
                             date = as.Date(fit[,1], origin = "2019-12-31"), 
                             Incidence = fit_E * eta, 
                             Incidence_report = p_symp * fit_E * eta, 
                             Incidence_nonreport = (1 - p_symp) * fit_E * eta)
  
  
  
  
  file_name_incidenceX <- paste(table.path,"/",REGION,"_",Last_date_used,"_",
                                "Days_pos_", Days_pos + 5, 
                                "_Incidence_Increased_contacts_by_", 
                                round((max_cont-1)*100,digits=0),
                                "percent", ".xlsx", sep ="")
  
  ## if you want to save reported incidence, non-reported incidence and the combined incidence
  ## uncomment below
  
  # write.xlsx(Incidence_df, file_name_incidenceX )
  
  
  
  
}

COL  <- rev(brewer.pal(n = 9, name = "Blues"))
COL <- COL[c(1,2,4,6,7)]




First_feb   <- as.numeric(as.Date("2020-02-01"))-as.numeric(as.Date("2019-12-31"))
First_march <- as.numeric(as.Date("2020-03-01"))-as.numeric(as.Date("2019-12-31"))

dayatyear_feb_jan <- c(First_feb, First_march, 92, 122, 
                               153, 183 , 214, 245,275, 306,336, 367)

NameDateFebJan <- as.Date(dayatyear_feb_jan,origin ="2019-12-31")




NameNumber <- paste("/",REGION , "_", Last_date_used, "_",
                    "Days_pos_", 5 + Days_pos, 
                    "_Reported_incidence_contact_increase", ".pdf", sep="")



pdf(paste(figure.path, NameNumber, sep=""), width=10, height=5)
  
  
  # Look at the estimated reported cases and fitted
  #bottom, left, top, right c(5.1, 4.1, 4.1, 2.1)
  par(mfrow = c(1,1), mar = c(6.1, 4.1, 4.1, 14.1)) # 
  yMax <- roundUpNice(max( fitted_incidence_reported_inc_contact[[length(Increase_contact)]],Observed_incidence)+100)
  #yMax <- 300
  plot(Day,Observed_incidence,type="p", ylab = "Reported cases", xlab ="", 
       xlim = range(dayatyear_feb_jan), 
       xaxt ='n', 
       main = paste("Estimated and simulated number of reported cases: ", REGION, sep=""),
       ylim = c(0, yMax))
  
  grid(nx=NA, ny=NULL)
  #grid
  abline(v=dayatyear_feb_jan,col = "lightgray", lty = "dotted", lwd = par("lwd"))
  #
  June10 <- which(t == Day_rise)
  lines(t,fitted_incidence_reported_same,  lwd = 2, col = "black")
  
  for(i in 1:length(Increase_contact)){
    plot_incidence <-  fitted_incidence_reported_inc_contact[[i]]
    lines(t[(June10 + 1):length(t)], plot_incidence[(June10 + 1):length(t)],  lwd = 2, col = COL[i])  
  }
  
  
  
  axis(side = 1, at = dayatyear_feb_jan, label = NameDateFebJan, las = 2)
  
  Legend_Name <- paste("Contact increase ", (Increase_contact-1)*100, "%", sep = "")
  
  op <- par(xpd=TRUE) 
  
  legend(390,yMax, c("Observed reported cases", "Fitted cases no change", Legend_Name), 
         lwd=c(1,2,2,2,2,2,2), pch = c(1,NA,NA,NA,NA,NA,NA), lty =c(NA,1,1,1,1,1,1), col = c("black","black",COL))
  
  
  par(op) 

dev.off()






