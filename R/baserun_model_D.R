source('load_files.R')
### Initialize modedl parameters
source('getParameters.R')
library(TMB)
library(ggplot2)

parameters <- getParameters()
par <- baseparameters.trimmed(parameters)

nyear <- length(parameters$years)

# Run the operating model 
OM <- runOM(parameters)


nsurvey <- 1 # Number of surveys 
nsize <- 20 # Number of survey bins 

# Calculate survey bins 
nsurvey.bin <- seq(1,parameters$wInf*2, length.out = parameters$nsize)
nsurvey.means <- (nsurvey.bin[2:20]+nsurvey.bin[1:19])/2

# Misc. calculations for model 
w <- par$w
Amat <- matrix(0,length(par$w))
idx <- 2:length(par$w)
Z0 <- par$mu0prefactor*par$wInf^par$mu0exponent

Amat[idx] <- -par$gg[(idx-1)]*1/par$dw[idx]
phiMat <- (1+(w/(par$wInf*par$alphaMature))^(-10))^(-1)
co <- 0.1 
phiMat[w<(co*par$alphaMature*par$wInf)] <- 0


#
n <- matrix(NA,(nsize-1)*nyear)
sd.surv <- matrix(NA,nsize-1,nyear)
yearidx <- seq(1,nyear*(nsize-1), by = (nsize-1))

for (i in 1:nyear){
  tmpn <- OM$survey[i,] # Order survey for tmb 
  n[yearidx[i]:(yearidx[i]+nsize-2)] <- tmpn
  
}


data <- list(      dt = 1,
                   wInf = par$wInf,
                   Zbase = Z0,
                   phiMat = phiMat,
                   iTimeMax = nyear,
                   mbins = nsurvey.bin,
                   yearidx = yearidx-1,
                   wlength = length(par$w),
                   nlength = nsize-1,
                   w = par$w,
                   dw = par$dw,
                   Mpredin = (parameters$a*par$A*w^(par$n-1)),
                   gg = par$gg,
                   alphaEgg = par$alpha.egg,
                   Sequ = par$Sequ,
                   Amat = Amat,
                   SDlog = log(0.1),
                   logSDR = log(0.6),
                   widxM = 50,
                   widxF = which.max(par$w),
                   a0 = log(1),
                   # Finally add the data
                   Catch = OM$Catch.obs,
                   survey = n,
                   nobs = length(n),
                   catchability = parameters$q,
                   logSDzero = log(0.1)
)
U <- matrix(0, 3 , data$iTimeMax)

parameters.tmb = list(Rmaxlog = log(parameters$Rmax),
                  #SDlog = log(0.2),#,#,
                  SDlogsurv = log(0.01),
                  logSDM = log(0.2),
                  logSDR = log(0.5),
                  logSDF = log(0.2),
                  lognF = log(0.05),
                  logmuF = log(10),
                  lognFs = log(0.2),
                  logmuFs = log(2),
                  logQ = log(1e-5),
                  lFstart = log(1),
                  #                  logphi = log(1),
                  U = U)#



compile("baserun_model_D.cpp")
dyn.load(dynlib("baserun_model_D"))

obj <-MakeADFun(data,parameters.tmb,DLL="baserun_model_D", random = 'U')

lower <- obj$par-Inf

# Parameters related to selectivity 
lower[names(lower) == 'lognF'] <- log(0.02)
lower[names(lower) == 'logmuF'] <- log(1)
lower[names(lower) == 'lognFs'] <- log(0.001)
lower[names(lower) == 'logmuFs'] <- log(1)
upper <- obj$par+Inf

# Restrain the selectvitiy parameters in realistic values
upper[names(upper) == 'lognF'] <- log(0.5)
upper[names(upper) == 'logmuF'] <- log(100)
upper[names(upper) == 'lognFs'] <- log(0.4)
upper[names(upper) == 'logmuFs'] <- log(100)


system.time(opt<-nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper))

rep<-sdreport(obj)
rep

# Uncertainty 
sdrep <- summary(rep)
rep.values<-rownames(sdrep)

unique(rep.values)

getUncertainty <- function(name){
  
  df <- data.frame(name = sdrep[rep.values == name,1])
  df$SE <- sdrep[rep.values == name,2]
  df$min <- df$name-2*df$SE
  df$max <- df$name+2*df$SE
  
  if(dim(df)[1] == nyear){
    df$year <- parameters$years
  }
  
  return(df)  
}

Catch <- getUncertainty('Catchest')
Fest <- getUncertainty('Fsave')
Bio <- getUncertainty('Bio')
Msave <- getUncertainty('Msave')
Rsave <- getUncertainty('R')

source('plotFit.R')
source('plotSE.R')
M.obs <- OM$M2[,data$widxM]
plotFit(F)

