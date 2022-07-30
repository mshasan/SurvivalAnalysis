# rm(list=ls())
# getwd()

library(readxl)
library(xlsx)
library(dplyr)
library(tidyverse)      # for data manipulation and analysis
library(lubridate)
library(scales)
library(multidplyr)     # for parallel computing
library(parallel)       # for finding the number of cores

dataXls <- read_excel("aaaa18_20-Jul-2022.xlsx")     #subject enrolled timeline
#View(dataXls)

date <- as.character(dataXls[6, -c(1:4)])
cumEnroll <- as.numeric(dataXls[29, -c(1:4)])
WklyEnroll <- cumEnroll - lag(cumEnroll)
WklyEnroll[is.na(WklyEnroll)] <- 0


wks <- seq(1, length(WklyEnroll), 1) # weeks
mnth <- wks*7/30.4375     

set.seed(616)

df <- tibble(Weeks = wks,
             months = mnth,
             cumEnroll = cumEnroll,
             WklyEnroll = WklyEnroll) 

df0 <- df %>% filter(WklyEnroll==0)
#dim(df0)
dfn <- df %>% filter(WklyEnroll!=0)
#dim(dfn) 
     
dfn2 <-   tidyr::uncount(dfn, WklyEnroll[WklyEnroll!=0], .remove = FALSE)

df2x <- df0 %>% bind_rows(dfn2) %>% arrange(Weeks)
#View(df2x)


fun_event <- function(iter, 
                      df2 = df2x,
                      hr = .65,
                      dropr = .05,
                      meddrop = 12,
                      medsurv0 = 36,
                      cutoff = 1:20) #cutoff weeks
{
  lambda0 = log(2)/medsurv0
  medsurv1 = medsurv0/hr
  lambda1 = log(2)/medsurv1
  
  mnth <- (cutoff*7)/30.4375
  
  df3 <- df2 %>%
    mutate(trtGrp = sample(1:2, n(), replace = TRUE),      #trt=1 is A
           t0 = months,                                                 # simulate entry time)
           td = rexp(n(), -log(1-dropr)/meddrop),                            # dropout time
           t1 = (trtGrp==1)*rexp(n(), lambda1) + (trtGrp==2)*rexp(n(), lambda0),
           event = ifelse(outer(X=t1+t0, Y=mnth, FUN='>'), 0, ifelse(t1<=td, 1, 0))) %>%
    colSums()
  
  return(df3)
}

#fun_event()

# setup multiple core for parallel computation
cl <- detectCores() # automatically detect the number of cores
cl
cluster <- new_cluster(cl)

# include R-libraries and functions to be used to the cluster
cluster_library(cluster, c("tidyverse", "dplyr"))
cluster_copy(cluster, c('fun_event', 'df2x'))


## ----applyMultiCore, message=FALSE, warning=FALSE----------------------------------------------------------
# apply functions using multiple cores
set.seed(616)

# Start the clock!
ptm <- proc.time()

iter = 10000
d <- bind_cols(iterVal = 1:iter)

iter_result <- d %>% 
  group_by(iterVal) %>%
  partition(cluster) %>%
  mutate(output = map(iterVal, 
                      ~fun_event(iter=., df2=df2x, cutoff=1:400))) %>% collect()

ptm_mc <- proc.time() - ptm
ptm_mc
# End the clock!


## ----iterationResults, message=FALSE, warning=FALSE--------------------------------------------------------
iter_result


## ----summaryResults, message=FALSE, warning=FALSE----------------------------------------------------------
# summarize the output
final_result <- iter_result %>% 
  unnest_wider(col = output) %>%
  pivot_longer(cols = starts_with('event'),
               names_to = 'eventWk',
               values_to = 'nEvent') 

final_result2 <- final_result%>%
                    group_by(eventWk) %>%
                    summarise_at('nEvent', list(p025 = ~quantile(., probs = 0.025),
                                                median = median, 
                                                p975 = ~quantile(., probs = 0.975)))

final_result3 <- final_result2 %>%
                     mutate(eWks = as.numeric(str_extract(eventWk, "[0-9]+"))) %>%
                     arrange(eWks)

#matplot(final_result4$Weeks, cbind(final_result4[,c(3, 6:8)]), type='l', xlab='Weeks', ylab='Events')

final_result4 <- df2x %>% 
                   right_join(final_result3, by=c("Weeks" = "eWks")) %>%
                   mutate(dates = lubridate::mdy("12-11-2022") + lubridate::weeks(Weeks- 1))
#View(final_result4)

IAdata = final_result4[which.max(final_result4$median==171),]
PAdata = final_result4[which.max(final_result4$median==245),]

# get time/event plots
timeLineData <- final_result4 %>%
                      pivot_longer(cols = c(cumEnroll, p025, median, p975),
                                   names_to = 'enrEvnt',
                                   values_to = 'nEnrEvnt') %>%
                      mutate(lineType = case_when(enrEvnt == 'cumEnroll' ~ '1',
                                                  enrEvnt == 'p025' ~ '2',
                                                  enrEvnt == 'median' ~ '3',
                                                  enrEvnt == 'p975' ~ '2')) 



timeLinePlots <- timeLineData %>% 
                      ggplot(aes(x = dates, group = lineType)) + 
                      geom_point(aes(y = nEnrEvnt, colour = lineType), size=.5) +
                      scale_y_continuous(breaks = seq(0, 750, 100)) +
                      scale_x_date(date_breaks = "12 week",
                                   expand = c(.0297, 0),
                                   labels = date_format("%d-%b-%Y")) +
                      geom_segment(aes(x = IAdata$dates , y = 0, xend = IAdata$dates, yend = 300, col = 'red')) +
                      geom_segment(aes(x = PAdata$dates , y = 0, xend = PAdata$dates, yend = 400, col = 'red')) +
                      geom_text(aes(x = IAdata$dates, y = 330, 
                                    label = paste0('Interim Analysis (', format(x=IAdata$dates, format="%d-%b-%Y"), ')'))) + 
                      geom_text(aes(x = PAdata$dates, y = 430, 
                                    label = paste0('Primary Analysis (', format(x=PAdata$dates, format="%d-%b-%Y"), ')'))) + 
                      geom_text(aes(x = as_date('2028-07-16'), y = 275, label = 'Events with 95% Confidence Interval')) +
                      geom_text(aes(x = as_date('2023-07-23'), y = 500, label = 'Enrollments')) +
                      labs(x='Dates',
                           y = 'Number of Enrolled / Events') +
                      ggtitle('Timeline of projected enrollment and events') + 
                      theme_bw() +
                      theme(axis.text.x = element_text(angle = 45, hjust = 1),
                            legend.position='none', 
                            legend.title = element_blank(),
                            plot.title = element_text(hjust = 0.5)) 







