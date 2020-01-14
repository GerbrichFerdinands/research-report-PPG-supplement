## -----------------------------------------------------------------------------
set.seed(1234)
source("Manuscript figures and matrices/printMat.R")

library(tidyverse)
library(PPGtools)
library(spam)
library(xtable)
library(cowplot)
library(ggthemes)

# color palette
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# path to store figures and matrices
mypath <- "Manuscript figures and matrices/figures/"


## -------------------------------------------------------------------------------------
library(PPGtools)
# prepare data
raw_signal <- prepInput(rec, "Green", tstart = 10, tstop = 80)

# Plot the red, green, and blue colour channel of the raw signal 
rgb <- 
  ggplot(data.frame(rec), aes(x = time)) +
  geom_line(aes(y=Red), colour = cbp2[7]) +
  geom_line(aes(y=Green), colour = cbp2[4]) +
  geom_line(aes(y=Blue), colour = cbp2[6]) + 
  labs(x = "Time (s)", y = "", title = "Raw PPG signal") +
  theme_cowplot() + 
  theme(legend.position = "none")

# save figure - raw data  ---------------------------------------------------------
postscript(paste0(mypath, "rawsignal.eps"))
rgb
dev.off()

# plot only green raw signal 
onlyg <-
  ggplot(data.frame(rec), aes(x = time)) +
  geom_line(aes(y=Green), colour = cbp2[1]) +
  labs(x = "Time (s)", y = "", title = "Raw PPG signal, green channel") +
  scale_x_continuous(limits = c(10, 80)) + 
  scale_y_continuous(limits = c(.15, .25)) + 
  theme_cowplot() + 
  theme(legend.position = "none")

# save ffigure  only green signal ----------------------------------------------
postscript(paste0(mypath, "greensignal.eps"))
onlyg
dev.off()

## -----------------------------------------------------------------------------
# first order differences in timesteps:
tdiff <- data.frame(diff1 = diff(raw_signal$time))

# summary statistics on interval widths
ts <- tdiff %>% 
  summarise(n = length(diff1)+1, # number of timepoints 
            min = min(diff1), 
            median = median(diff1), 
            mean = mean(diff1), 
            max = max(diff1), 
            sd = sd(diff1), 
            range = max-min)
ts 
ts$max/ts$min # minimum interval is 12 times shorter then maximum interval 


## -------------------------------------------------------------------------------------
# plot distribution of time interval widths in the signal
timedist <-
  ggplot(tdiff, aes(diff1, colour = diff1)) +
  geom_histogram(bins = 50, colour=NA, fill=cbp2[1]) + 
  labs(x = "Interval width", y = "") +
  scale_x_continuous(breaks=seq(0.01,0.09, 0.01)) +
  scale_y_continuous(
    # don't expand y scale at the lower end
    expand = expand_scale(mult = c(0, 0.05))) + 
  theme_cowplot()
timedist

# save figure - time hist  ---------------------------------------------------------
postscript(paste0(mypath, "timehist.eps"))
timedist
dev.off()




## -------------------------------------------------------------------------------------
# sampling rate of the smartphone camera (frames per second)
m <- length(raw_signal$time) # number of frames
sec <- tail(raw_signal$time, n=1) - head(raw_signal$time, n=1) # number of seconds
fps <- m/sec
fps # frames per second
1/fps # expected interval width at a constand sampling rate

# scale raw signal
raw_signal$time_scaled <- raw_signal$time *fps


# for four different lambda values 
  lambda <- 10^(seq(-1, 5, 2))
  lambda <- matrix(lambda, dimnames = list(paste0('lambda_', 1:length(lambda)), NULL))

# for three order differences 
  for(d in 1:3){
      
    # 1. FILTER THE DATA
      
      # smooth assuming equal timesteps (ES)
      zES <- smoothWE(raw_signal = raw_signal, lambda = lambda , d = d, uni = TRUE)
      #pdates <- data.frame(raw_signal, zES) 
      # plot for smoothing equal timesteps (pes)
      pes1 <-
        plotLambda(raw_signal, zES, "") + 
        scale_x_continuous(limits = c(30, 40)) + 
        scale_y_continuous(limits = c(0.15, 0.225)) +
        theme(legend.position = "none")
      # save plot to .eps file 
      postscript(paste0(mypath, "pes1d", d, ".eps"))
      plot(pes1)
      dev.off()
      
      # smooth assuming unequal time steps (UES)
      zUES <- smoothWE(raw_signal = raw_signal, lambda = lambda , d = d, uni = FALSE)
      #pdatues <- data.frame(raw_signal, zUES)
      # plot for smoothing assuming unequal timesteps (pues)
      pues1 <- 
        plotLambda(raw_signal, zUES, "") + 
        scale_x_continuous(limits = c(30, 40)) + 
        scale_y_continuous(limits = c(0.15, 0.225))+
        theme(legend.position = "none")
      # save plots to .eps file 
      postscript(paste0(mypath, "pues1d", d, ".eps"))
      plot(pues1)
      dev.off()
   
      
     # 2. DETREND THE RAW SIGNAL WITH A LARGE LAMBDA VALUE  
     detrendES <- data.frame(raw_signal, 
                        ES = raw_signal$Green-zES[,"l4"],
                        UES = raw_signal$Green-zUES[,"l4"]) 
     
     # assuming equal timesteps (ES)
     dtES <-
        ggplot(detrendES, aes(x=time, y = ES)) +
        geom_line() +
        labs(x = "Time (s)", 
        y = "signal") +
        scale_x_continuous(limits = c(30, 40)) + 
        theme_cowplot() 
      
      # save plot
      postscript(paste0(mypath, "detrended_es_d", d, ".eps"))
      plot(dtES)
      dev.off()
      
      # assuming unequal time steps (UES)
      dtUES <-
        ggplot(detrendES, aes(x=time, y = UES)) +
        geom_line() +
        labs(x = "Time (s)", 
        y = "signal") +
        scale_x_continuous(limits = c(30, 40)) + 
        theme_cowplot() 
      
      # save plot
      postscript(paste0(mypath, "detrended_ues_d", d, ".eps"))
      plot(dtUES)
      dev.off()

}


## -------------------------------------------------------------------------------------
# set up grid of lambda values 
lambda <- 10^(seq(-3, 1, by = .01))
lambda <- matrix(lambda, dimnames = list(paste0('lambda_', 1:length(lambda)), NULL))
plot(lambda)

# empty frame to save  cross-validation standard errors
cve <- data.frame(lambda = lambda, 
                  scv_e1 = 0, 
                  scv_e2 = 0, 
                  scv_e3 = 0, 
                  scv_ue1 = 0,
                  scv_ue2 = 0,
                  scv_ue3 = 0)

# indicators:
steps <- c(TRUE, FALSE) # equal and unequal timesteps
I <- matrix(1:6, ncol=2, nrow=3) # index for cve dataframe 

for(s in 1:2){ # for equal or unequal timesteps
    U <- steps[s] 
    for(d in 1:3){ # for 3 order differences
      i <- I[d,s]
      
      # estimate trend
      trend <- smoothWE(raw_signal = raw_signal, lambda = 10^5 , d = d, uni = U)
      # detrend the raw signal
      detrended_signal <- raw_signal
      detrended_signal$Green <- detrended_signal$Green - trend[,1]

      # compute cross-validation standard error over the whole lambda grid
      for(l in 1:length(lambda)){
          z <- smoothWE(raw_signal=detrended_signal, 
                        lambda = lambda[l,], 
                        d=d, uni = U, cv = TRUE)
          cve[l,i+1] <- z$cve
      }
    }
}


# plot the results of the grid search in the 6 (2*3) conditions. 
for(i in 1:6){
  nm <- names(cve)[i+1] # condition
  grides <- 
    ggplot(data = cve, aes(x=lambda, y = cve[,i+1])) +
    labs(title = "", x = "lambda", y = "cross-validation standard error") +
    theme_bw() + 
    geom_line() 
    
    grides
      
    # save plot
    postscript(paste0(mypath, nm, ".eps"))
    plot(grides)
    dev.off()
}

## -------------------------------------------------------------------------------------


