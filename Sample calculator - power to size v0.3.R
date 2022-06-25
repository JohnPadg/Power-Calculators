##############################################################
#
#   Padgham's Unofficial WEBEXPO R scripts
#
#   Power to Sample Size Analysis
#  
#   'Give me your Power, I'll give you sample numbers'
#
#   V0.3 May 2022  
#
#############################################################


start_time <- Sys.time()

###### library ####

library(rjags)

library(ggplot2)

require(gridExtra)

library(here)

library(MASS)

library(tolerance)

##### WEBEXPO sourcing programs ####

## data preparation / generation

source(here("RANDOM SAMPLE GENERATION", "webexpo.seg.randomgeneration.R"))

source(here("DATA PREPARATION", "SEG ANALYSIS", "webexpo.seg.dataprep.R"))

# JAGS models

source(here("JAGS MODELS", "SEG ANALYSIS", "webexpo.seg.mainbayesian.R"))

source(here("JAGS MODELS", "SEG ANALYSIS", "webexpo.seg.informedvarbayesian.R"))

source(here("JAGS MODELS", "SEG ANALYSIS", "webexpo.seg.informedvarbayesian.models.R"))

source(here("JAGS MODELS", "SEG ANALYSIS", "webexpo.seg.uninformativebayesian.R"))

source(here("JAGS MODELS", "SEG ANALYSIS", "webexpo.seg.uninformativebayesian.models.R"))

source(here("JAGS MODELS", "SEG ANALYSIS", "webexpo.seg.riskbandbayesian.R"))

source(here("JAGS MODELS", "SEG ANALYSIS", "webexpo.seg.riskbandbayesian.models.R"))


# Data interpretation

source(here("RESULT INTERPRETATION", "SEG ANALYSIS", "webexpo.seg.interpretation.R"))

source(here("RESULT INTERPRETATION", "SEG ANALYSIS", "webexpo.seg.summary.R"))




####### DEFINE VALUES OF INTEREST ##############


#Relevant OEL
input.OEL <- 100

# Criteria Percentile
target_perc <- 95 #Usually 95

# Credible interval
probacred <- 40 # 40 = 70% upper credible limit

#Desired power of sample plan
desired.power <- 0.8 #arbitrary - between 0 & 1

#Starting sample size
init_samp_size <- 3 #minimum 1 - higher saves time


#Number of simulated experiments to test power of each sample size
nSimulatedDataSets <- 500 #arbitrary large - more is slower but better


#Report power for each sample size?
verbose = TRUE #True - report all, False - report only successful sample size


#Data input
sample.historic <- c("10","21","34","24","33","13","22","21","17")





##################  Create Hypotheses ################


          # Bayesian calculation of historic data using JAGS models
    
          mcmc.historic <- Webexpo.seg.globalbayesian.jags( data.sample = sample.historic ,
                                                     is.lognormal = TRUE , 
                                                     error.type = "none" ,
                                                     oel = input.OEL, prior.model = "informedvar" )
            
          # Identify hypothesis
          historic.chain <- exp( mcmc.historic$mu.chain + qnorm( target_perc / 100 ) * mcmc.historic$sigma.chain )
          historic.utl <- quantile( historic.chain , 1 - ( 100 - probacred ) / 200 )
          
          if(historic.utl<input.OEL) {
            hypothesis.compliant = TRUE # Is the P95 70%UCL estimated to be less than OEL? 
          } else {
            hypothesis.compliant = FALSE
          }



          # Record mu & sigma parameter distributions to sample randomly from in tests
          
          distribution.historic.mu <- mcmc.historic$mu.chain

                                                 
          distribution.historic.sigma <- mcmc.historic$sigma.chain

        
          
################## GENERATE SAMPLES FROM PARAMETERS + RUN TEST  ################
          
          # Set up records of the tests below
          list.simulated.P95 <- data.frame(matrix(ncol=1,nrow=nSimulatedDataSets))
          
          verbose.description <- data.frame(matrix(ncol=2,nrow=0))
          colnames(verbose.description) <- c("Sample Size","Power")
          
          
          
          # Initialize loop for incrementing `sample_size`:
          sample_size <- init_samp_size
          not_powerful_enough = TRUE
          too_many_samples = FALSE
          
          
          
          # Increment `sample_size` until desired power is achieved:
          
          while(not_powerful_enough) {
            
            list.simulated.musig <- data.frame(generated.Mu = replicate(nSimulatedDataSets,sample(distribution.historic.mu,1,replace=TRUE)),
                                               generated.sigma = replicate(nSimulatedDataSets,sample(distribution.historic.sigma,1,replace=TRUE)))
            
            list.simulated.samples <- data.frame(samples = replicate(nSimulatedDataSets,
                                                                     rlnorm(sample_size,list.simulated.musig$generated.Mu,list.simulated.musig$generated.sigma)))
            
            list.simulated.tests <- apply(list.simulated.samples,2,log)
            
            #Calculate UTL70,95 function
            UTL_function <- function(x){
              y<-normtol.int(x,alpha=0.3,P=0.95,side=1)
              return(y[,5])
            }
            
            list.simulated.P95 <- apply(list.simulated.tests,MARGIN = 2,FUN = UTL_function)
            list.simulated.P95 <- lapply(list.simulated.P95,exp)
            
            #Calculate Power based on if goal is reached
            
            if (hypothesis.compliant==TRUE){  
              power.level <- sum(list.simulated.P95 < input.OEL)/length(list.simulated.P95)
              
              verbose.description[nrow(verbose.description)+1,] <- c(paste("For sample size = ",sample_size),paste("- power = ", power.level)) 
            }
            else{
              power.level <- sum(list.simulated.P95 > input.OEL)/length(list.simulated.P95)
              
              verbose.description[nrow(verbose.description)+1,] <- c(paste("For sample size = ",sample_size),paste("- power = ", power.level))
            }
            
            
            # Stop loop if desired power is obtain, or +1 sample size and repeat
            if (power.level > desired.power) {  
              not_powerful_enough = FALSE
            } 
            else {
              sample_size <- sample_size + 1
            }
            
            #Emergency Stop
            if (sample_size > 50) {
              not_powerful_enough = FALSE
              too_many_samples = TRUE
            }
            
          } 
          #End if desired power is achieved
          
######### REPORT POWER LEVEL ##############                  
          
          if (verbose==TRUE){
          print(verbose.description)
            }

          if (verbose==FALSE){
            print(tail(verbose.description))
          }
          
          if (hypothesis.compliant==TRUE){
            print("This power is based on a hypothesis that the UTL < OEL")
          } else {
            print("This power is based on a hypothesis that the UTL > OEL")
          }
          
          if (too_many_samples==TRUE){
            print("Sample size is greater than 50... the program was stopped short")
          }
            
end_time <- Sys.time()            
            
total_time <- end_time - start_time  
print(total_time)
          