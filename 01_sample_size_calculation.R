#rm(list=ls())

# load necessary libraries for this simulation
library(survival)       # for survival analysis
library(tidyverse)      # for data manipulation and analysis
library(multidplyr)     # for parallel computing
library(parallel)       # for finding number of cores

#===============================================================================

data_and_analysis <- function(hr = .80,        # default hazard ratio treatment vs. control
                              n0 = 100,        # default sample size for control
                              n1 = 50,         # default sample size for treatment
                              dropr = .05,     # default dropout rate in percentage
                              meddrop = 12,    # default median dropout time in months
                              etime = 24,      # default enrollment time in months
                              atime = 42,      # default final analysis time in months
                              medsurv0 = 36)   # default median survival time in months for control
{
  n = n0 + n1
  
  # calculate rate for the control and treatment from median survival time
  lambda0 = log(2)/medsurv0
  medsurv1 = medsurv0/hr
  lambda1 = log(2)/medsurv1
  
  t0 = runif(n)*etime        # simulate entry time
  td = if(dropr==0){         # simulate dropout time
    rep(Inf, n)
  } else {
    rexp(n, -log(1-dropr)/meddrop)
  } 
  t1 = c(rexp(n1, lambda1), rexp(n0, lambda0))            # simulate end time
  arm = c(rep('trt', n1), rep('ctr', n0))                 # create arm covariate
  status = ifelse(t1+t0>=atime, 0, ifelse(t1<=td, 1, 0))  # assign patient status 
  time = as.numeric(pmin(atime-t0, pmin(t1, td)))         # obtain survival time
  
  dat <- bind_cols(t0=t0, td=td, t1=t1, arm=arm)       
  
  # fit survival model and obtain summary
  fit = survival::coxph(Surv(time, status) ~ arm, data = dat) %>% summary()
  
  # estimate true HR and 80% and 90% upper confidence interval
  estHR = exp(fit$coefficients[1])
  uCi80 = exp(fit$coefficients[1] + qnorm(.90)*fit$coefficients[3])
  uCi90 = exp(fit$coefficients[1] + qnorm(.95)*fit$coefficients[3])
  
  # logical output (1, 0) if estimated true HR > .80
  pHR = as.numeric(estHR > .80)
  
  # logical output (1, 0) if the upper limit of CI below 1.0 
  puCi80 = as.numeric(uCi80 < 1.0)
  puCi90 = as.numeric(uCi90 < 1.0)
  
  # logical output (1, 0) if null hypothesis (H0: HR=1) is rejected
  pPval80 = as.numeric(fit$coefficients[5] < .20)
  pPval90 = as.numeric(fit$coefficients[5] < .10)
  
  # keeping simulated data as well as statistics in a nested structure
  return(list(data = tibble(arm=arm, status=status, time=time),
              stat = tibble(trueHR=hr, n0=n0, n1=n1, n=n, estHR=estHR, pHR=pHR, 
                          puCi80=puCi80, puCi90=puCi90,
                          pPval80=pPval80, pPval90=pPval90)))
}



#===============================================================================

power_by_HR_N <- function(iter,                  # number of iterations
                          hr = c(.65, .70, .80), 
                          n0 = 350,  
                          n1 = c(30, 60, 90),   
                          dropr = .05,      
                          meddrop = 12,    
                          etime = 24,          
                          atime = 42,        
                          medsurv0 = 36) 
  
{
  # find combinations of hr, n0, and n1; apply function to obtain statistics,
  # then only extract statistics (stat) from the results
  comboFit <- expand_grid(hr, n0, n1) %>%
                  mutate(results=pmap(., data_and_analysis)) 
  return(comboFit)
}


#===============================================================================
# setup multiple core for parallel computation
cl <- detectCores() # automatically detect the number of cores
cl
cluster <- new_cluster(cl)

# include R-libraries and functions to be used to the cluster
cluster_library(cluster, c("tidyverse", "survival", "dplyr"))
cluster_copy(cluster, c('power_by_HR_N', 'data_and_analysis'))


# apply functions using multiple cores
set.seed(1212021)

# Start the clock!
ptm <- proc.time()

iter=10
d <- bind_cols(iterVal = 1:iter)

iter_result <- d %>% 
  group_by(iterVal) %>%
  partition(cluster) %>%
  mutate(output = map(iterVal, 
                      ~power_by_HR_N(iter=., 
                                     hr = c(.65, .90, 1.0, 1.1, 1.5),           
                                     n0 = c(250, 350),  
                                     n1 = c(100, 200, 300)))) %>% collect()

proc.time() - ptm
# End the clock!


# summarize the output
final_result <- iter_result %>% 
                  unnest(cols=output) %>%
                  .$results %>%
                  map_dfr('stat') %>%
                  group_by(trueHR, n0, n1, n) %>%
                  summarise_all(mean, na.rm=TRUE)


#===============================================================================
# apply functions using single core
# single core----------
set.seed(1212021)
# Start the clock!
ptm <- proc.time()
iter=5000
d <- bind_cols(iterVal = 1:iter)

iter_resultx <- d %>% 
                  group_by(iterVal) %>%
                  mutate(output = map(iterVal, 
                                      ~power_by_HR_N(iter=., 
                                                     hr = c(.65, .90, 1.0, 1.1, 1.5),           
                                                     n0 = c(250, 350),  
                                                     n1 = c(100, 200, 300)))) 
proc.time() - ptm
# End the clock!


# summarize the output
final_resultx <- iter_result %>% 
                  unnest(cols=output) %>%
                  .$results %>%
                  map_dfr('stat') %>%
                  group_by(trueHR, n0, n1, n) %>%
                  summarise_all(mean, na.rm=TRUE)

#===============================================================================




