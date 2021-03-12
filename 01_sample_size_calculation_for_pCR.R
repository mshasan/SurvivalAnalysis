## ----setup, include=FALSE---------------------------------------------------------------------------------
#knitr::opts_chunk$set(echo = TRUE)


## ----loadlibraries, message=FALSE, warning=FALSE----------------------------------------------------------
#rm(list=ls())
# setwd('/home/mhasan10/BLC2002')

# load necessary libraries for this simulation
library(tidyverse)      # for data manipulation and analysis
library(multidplyr)     # for parallel computing
library(parallel)       # for finding the number of cores
library(DT)             # for interactive tables


## ----functionSimulateSingleDataset, message=FALSE, warning=FALSE------------------------------------------
#===============================================================================

data_and_analysis <- function(pCR = 0.30,   # default true pCR rate
                              p0 = 0.20,    # default thresholds comparing estimated pCR rate with
                              n = 60,       # default sample size 
                              alpha = 0.20) # False positive rate      
{
  b1 = b2 = 0.5    # parameters of beta prior 
  m = 0:n          # observed number of participants 

  #Obtaining the lest number of observed participants needed to obtain >(1-alpha/2)        
  #probability of having estimated pCR rate is greater than the given threshold p0 
  m_th = bind_cols(m = m, p = pbeta(p0, b1+m, b2+n-m, lower.tail = FALSE)) %>% 
          mutate(pgt = p>(1-alpha/2)) %>%
          filter(pgt==TRUE) %>%
          slice_head() %>% 
          select('m') %>%
          as.numeric()
  
  # out1-out4 below are different ways of obtaining the logical outcomes,
  # consequently obtaining probabilities using all iterations of simulation
  
  # Checking (TRUE/FALSE) whether sum of samples (0, 1) exceeds m_th number
  out1 = sample(x=c(0, 1), size=n, prob=c(1-pCR, pCR), replace=TRUE) %>%
            sum() >= m_th
  
  # generate binomial sample, calculate P(estPCR>p0|pCR), then compare whether
  # this is greater than the (1-alpha/2)
  x = rbinom(n=1, size=n, prob=pCR) 
  out2 = pbeta(p0, b1+x, b2+n-x, lower.tail = FALSE) > (1-alpha/2)
  
  # conduct one-sample proportion tests for H0:phat<=p0 vs. H1:phat>p0,
  # then compare pvalue with Type I error (alpha) 
  out3 = prop.test(x=x, n=n, conf.level=1-alpha/2, p=p0, 
                   alternative="greater", correct=TRUE)$p.value < alpha/2

  # perform normal-approximation similar to out2
  a = b1+x
  b = b2+n-x
  out4 = pnorm(p0, mean=a/(a+b), sd=sqrt((a*b)/((a+b+1)*(a+b)^2)), 
               lower.tail = FALSE) > (1-alpha/2)
  
  # futility probability P(phat<=p0) is greater than (1-alpha/2)
  futility = pbeta(p0, b1+x, b2+n-x) > (1-alpha/2)

  # keeping statistics in a nested structure
  return(stat = tibble(reqSubj = m_th,
                       reqSubjBinom = x,
                       sampleBased=out1,
                       binomBased=out2,
                       ptestBased=out3,
                       normalBased=out4,
                       futility = futility))
}


## ----functionSimulateMultipleDatasets, message=FALSE, warning=FALSE---------------------------------------
#===============================================================================

power_by_combo <- function(iter,                  # number of iterations
                         pCR = c(.15, .30, .40, .55),
                         p0 = c(.20, .25, .35, .40, .45),
                         n = c(30, 50, 60, 70, 100),
                         alpha = c(.2, .4)) 
  
{
  # find combinations of pCR, n, and p0;then 
  # apply function to obtain statistics
  comboFit <- expand_grid(pCR, n, p0, alpha) %>%
                  mutate(results=pmap(., data_and_analysis)) 
  return(comboFit)
}


## ----setupMultiCore, message=FALSE, warning=FALSE---------------------------------------------------------
#===============================================================================
# setup multiple core for parallel computation
cl <- detectCores() # automatically detect the number of cores
cl
cluster <- new_cluster(cl)

# include R-libraries and functions to be used to the cluster
cluster_library(cluster, c("tidyverse", "dplyr"))
cluster_copy(cluster, c('power_by_combo', 'data_and_analysis'))


## ----applyMultiCore, message=FALSE, warning=FALSE---------------------------------------------------------
# apply functions using multiple cores
set.seed(616)

# Start the clock!
ptm <- proc.time()

iter=5000
d <- bind_cols(iterVal = 1:iter)

iter_result <- d %>% 
  group_by(iterVal) %>%
  partition(cluster) %>%
  mutate(output = map(iterVal, 
                      ~power_by_combo(iter=., 
                                     pCR = c(.15, .30, .40, .55),
                                     p0 = c(.20, .25, .35, .40, .45),
                                     n = c(30, 50, 60, 70, 100),
                                     alpha = c(.2, .4)))) %>% collect()

ptm_mc <- proc.time() - ptm
ptm_mc
# End the clock!


## ----iterationResults, message=FALSE, warning=FALSE-------------------------------------------------------
iter_result


## ----summaryResults, message=FALSE, warning=FALSE---------------------------------------------------------
# summarize the output
final_result <- iter_result %>% 
                  unnest(cols=output) %>%
                  unnest(cols=results) %>%
                  group_by(pCR, n, p0, alpha) %>%
                  summarise_all(mean, na.rm=TRUE)

# Create interactive table
final_result %>%
  select(-'iterVal') %>% 
  mutate_if(is.numeric, round, 2) %>%
  mutate(reqSubjBinom=ceiling(reqSubjBinom)) %>%
  DT::datatable(extensions = 'Buttons',
            options = list(dom = 'Blfrtip',
                           buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                           lengthMenu = list(c(10, 25, 50, -1),
                                             c(10, 25, 50, "All"))))





