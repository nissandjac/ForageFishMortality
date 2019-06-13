# Create a baserun using just the single species model, but with error terms 
runOM <- function(parameters,seed = 3844335){
# Run the single species simulation and estimate the natural mortality + fishing mortality from survey and catch
source('load_files_OM.R')
  
A <- h*0.6*0.4 # Growth constant 

param <- baseparameters(wInf,Rmax,h,a,A, kappa)
w <- makegrid(param)[[1]]
param$fishing <- "Trawl"
param$F0 <- 0
param$Mpred <- a*A*w^(-0.25)

# Create a random walk mortality 
set.seed(seed)   #3844334

nyear <- length(years)

Fforage <- matrix(NA,nyear)
Fforage[1] <- 1# Initial fishing 
param$Mpred <- a*A*w^(-0.25)
phi <- matrix(NA,nyear)
phi[1] <- 1
M <- NA
Merr <- 0
errF <- 1

for (i in 2:nyear){
  errF[i] <- exp(rnorm(n = 1, mean = 0, sd = Fdev))
  Fforage[i] <- Fforage[i-1]*errF[i]

  Merr[i] <- rnorm(n = 1, mean = 0, sd = Mdev)
  phi[i] <- phi[i-1]*exp(Merr[i])
  Mpred <- a*param$A*w^(param$n-1)*phi[i]
  M[i] <- Mpred[50]

}

F0 <- Fforage
param$dt <- 1
param$tEnd <- nyear
#

SF.year <- IterateSpectrum.modified3(param,F0,Merr)
w.app <- t(matrix(rep(SF.year$w,nyear),nrow = length(SF.year$w)))
dw.app<- t(matrix(rep(SF.year$dw,nyear),nrow = length(SF.year$w)))
Catch <- rowSums(SF.year$Fin*(SF.year$N*w.app*dw.app))


nsurvey.bin <- seq(1,wInf*2, length.out = nsize)

nsurvey.means <- (nsurvey.bin[2:nsize]+nsurvey.bin[1:(nsize-1)])/2
nsurvey.dw <- nsurvey.bin[2:nsize]-nsurvey.bin[1:(nsize-1)]
size <- SF.year$w# 10^(seq(log10(1),log10(6000),length.out = nsize)) # Survey weight bins

surv.sel <- (1+(SF.year$w/(0.01 * wInf))^-2)^-1 # Survey selectivity
surv.sel[SF.year$w < 1] <- 0 # No fish smaller than 1 g

N.survey <- array(NA,dim = c(nyear,nsize-1)) # There always a bin less

# Set up the measuring vector

cuts <- cut(SF.year$w,nsurvey.bin)

for (i in 1:nyear){
  
    err <-  exp(rnorm(n = nsize-1,mean = 0, sd = surv.sd)) #  Error within that survey
    x <- tapply(SF.year$N[i,]*q*surv.sel, cuts, sum)
    x[is.na(x)] <- 0
    
    N.survey[i,] <-  as.numeric(x)*err # Add observation error to the survey
}

surv.sel.app <- t(matrix(rep(surv.sel,nyear),nrow = length(SF.year$w)))

err.catch <- exp(rnorm(n = nyear,mean = 0, sd = sd.catch))
Obs.Y <- Catch*err.catch # Add observation error to the catch
 

OM.mortality <- list(M2 = SF.year$M2, Catch.obs = Obs.Y, Catch.true = Catch,
                     survey = N.survey, SF = SF.year, Merr = Merr, Ferr = errF, pars = param)

# # Consider adding recruitment deviation and growth deviation
return(OM.mortality)
}

