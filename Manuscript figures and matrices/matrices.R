# Gerbrich Ferdinands ----------------------------------------------------------
# Research report matrices producer --------------------------------------------

# path to store figures and matrices
mypath <- "Manuscript figures and matrices/matrices/"
source("Manuscript figures and matrices/printMat.R")

# ------------------------------------------------------------------------------
# example: signal of length m=5 
m <- 5

# first order differences matrix
d1 <- diff(diag(m), differences = 1)
D1 <- xtable(d1, digits = 0, align=rep("",ncol(d1)+1))
print(D1, 
      floating=FALSE, 
      tabular.environment="bmatrix", 
      hline.after=NULL, 
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      file = paste0(mypath, "/mat_d_1.txt"))

# ------------------------------------------------------------------------------
# second order differences matrix
d2 <- diff(diag(m), differences = 2)
D2 <- xtable(d2, digits = 0, align=rep("",ncol(d2)+1)) 
print(D2, 
      floating=FALSE, 
      tabular.environment="bmatrix", 
      hline.after=NULL, 
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      file = paste0(mypath, "mat_d_2.txt"))

# ------------------------------------------------------------------------------
library(PPGtools)
raw_signal <- prepInput(rec, "Green", tstart = 30, tstop = 40)
sec <- tail(raw_signal$time, n=1) - head(raw_signal$time, n=1)# number of seconds
fps <- length(raw_signal$time)/sec # frames p/s

Y <- raw_signal[25:29,]
m <- length(Y[,"time"])
y <- matrix(Y[,"Green"])
lambda <- 1
t <- matrix(Y[,"time"])

# time steps matrix 
steps <- xtable(t(t), 
                digits = 2, 
                align=rep("",ncol(t(t))+1))

print(steps, floating=FALSE, 
      tabular.environment="bmatrix", 
      hline.after=NULL, 
      include.rownames=FALSE, 
      include.colnames=FALSE,
      file = paste0(mypath, "mat_timesteps.txt"))

## -------------------------------------------------------------------------------------
# VD matrix
d <- 1 
# E = identity matrix
E <- diag(m)
# D = differences matrix
D <- diff(E, differences = d)
# V = m-1 by m-1 matrix with 1/delta x on its diagonal.
V <- diag.spam(m-d)
V <- as.matrix(V)
#V <- diag.spam(m-1)
diag(V) <- 1/diff(t, differences = d)

VD <- V%*%D

vd <- xtable(VD,digits = 2, align=rep("",ncol(VD)+1)) # We repeat empty string 6 times
print(vd, 
      floating=FALSE, 
      tabular.environment="bmatrix", 
      hline.after=NULL, 
      include.rownames=FALSE, 
      include.colnames=FALSE,
      file = paste0(mypath, "mat_VD.txt"))

## -------------------------------------------------------------------------------------
# VD matrix with scaled time steps 
scaled_steps <- t %*% fps
sV <- diag.spam(m-d)
diag(sV) <- 1/diff(scaled_steps, differences = d)
sV <- as.matrix(sV)
sVD <- sV %*% D

scaled_VD <- xtable(sVD,digits = 2, align=rep("",ncol(sVD)+1)) # We repeat empty string 6 times
print(scaled_VD, 
      floating=FALSE, 
      tabular.environment="bmatrix", 
      hline.after=NULL, 
      include.rownames=FALSE, 
      include.colnames=FALSE,
      file = paste0(mypath, "mat_VD_scaled.txt"))

## -----------------------------------------------------------------------------