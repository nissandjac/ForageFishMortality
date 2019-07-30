# Run the sprat from the North Sea (closely resembles our model orgaNism
# load('NS_survey_spectra.Rdata')
source('load_files.R')
source('gradient.R')
source('optimize_sizebased.R')
source('run_optimize_sizebased.R')
source('baseparameters_trimmed.R')
source('C:/Users/Nis/Dropbox/UW/Fluctuations paper/sprat/IterateSpectrumModified_pt3_TMB.R')

library(ggplot2)
library(dplyr)

load('data/sprat_datras.Rdata')

a.weight <- 0.0058	# From fishbase 
b.weight <- 3.08
Linf <- 14.7
Winf <- a.weight*Linf^b.weight # Looks like a fair number considering data
nsize <- 12# Number of survey bins 

K <- 0.5
wmat <- a.weight*10^b.weight # 
myF <- 8
kappa <- 1e8 # Used to initialize

Wmax <- max(df.sprat$weight[df.sprat$CPUE > 0])
eta.mat <- wmat/Winf
eta.F <- 5/Winf # Pretty high L50

# Metabolic params 
n <- 3/4

w <- unique(df.sprat$weight)

A <- b.weight*a.weight^(1-n)*K*Linf^(b.weight*(1-n))

h <- A/(0.6*0.4)

a <- 0.35#0.35*1/(b.weight*K*eta.mat^(1-n)) # Prior natural mortality

Rmax <- 1e8
par <- baseparameters.trimmed(Winf,h,a,A,Rmax,eta.mat)
w <- par$w

par$fishing <- "Trawl"
par$F0 <- 1
par$Mpred <- a*A*w^(n-1)
par$nF <- eta.F
par$tEnd <- 100
par$Rmax <- Rmax
par$nSpecies <- 1
par$dt <- 1

S <- IterateSpectrum.modified3.tmb(par,S = NA, F0 = rep(0.3, par$tEnd), Mdev = rep(0, par$tEnd))

# Calc Fmsy 
source('YieldCalc.R')
Fvec <- seq(0,2, length.out = 30)
Yield <- matrix(NA,30)
for (i in 1:length(Fvec)){
  S <- IterateSpectrum.modified3(par, S= NA, F0 = rep(Fvec[i], par$tEnd), Mdev = rep(0, par$tEnd))
  N <- S$N[par$tEnd/par$dt,]
  Yield[i] <- sum(N*S$w*S$dw*S$Fin[par$tEnd/par$dt,])
  
}

plot(Fvec,Yield, type = 'l')

# Official assessment uses year 1974 as the first year 
year <- unique(df.sprat$year)
year <- year[year <2013]
nyear <- length(year)

df.sprat <- df.sprat[df.sprat$year < 2013,] # Not anymore catch data

# Find the survey selectivity
sprat.mean <- df.sprat %>%
  group_by(weight) %>%
  summarise(N = mean(CPUE, na.rm = T))

sel.survey <- (1+(S$w/(0.1 * Winf))^-5)^-1
sel.survey[S$w < 0.1] <- 0

bins <- exp(seq(log(1),log(Winf[1]*1.1), length.out = nsize))
m.bins <- (bins[2:nsize]+bins[1:(nsize-1)])/2
dw <- (bins[2:nsize]-bins[1:(nsize-1)])

cuts <- cut(sort(unique(df.sprat$weight)),bins, right = F,labels = F)
cuts[is.na(cuts)] <- max(cuts, na.rm =T)
cuts[sort(unique(df.sprat$weight)) < bins[1]] <- NA

#N.est.bin <- tapply(N.est, cuts, sum)

# Load sprat data from the SMS model 
Sprat.obs <- read.csv('Sprat_sms.csv')
sprat.ass <- read.table('sprat_assessment.txt', header = T)
sprat.ass.2 <- read.csv('Sprat_assessment_output.csv')

Landings.obs <- c(sprat.ass$Landings, sprat.ass$Landings[length(sprat.ass$Landings)])

Catch.obs <- Sprat.obs %>%
  group_by(Year) %>%
  summarise(Catch = sum(Yield, na.rm = T),
            Biomass = sum(Bio, na.rm = T),
            SSB = sum(SSB))
#Catch.obs$Catch <- Catch.obs$Catch*0.25 # It's in quarters


M.obs <- Sprat.obs %>%
  group_by(Year) %>%
  summarise(M = sum(M2))
# Read the sprat mortality 
Sprat.M <- read.csv('spratNS_mortality.csv')
M.obs <- data.frame(M = Sprat.M[,5], Year = Sprat.M$year) # 12 g fish
# plot(Sprat.M[,2], ylim = c(0,2.5), type = 'l')
# lines(Sprat.M[,3])
# lines(Sprat.M[,4])
# lines(Sprat.M[,5])5
# Find weights for reference 
Sprat.weight <- read.csv('spratNS_weight.csv')
Sprat.weight <- colMeans(Sprat.weight[,2:5])*1000 # in G


F.obs <- Sprat.obs$F[Sprat.obs$Age == 3]


sprat.optim <- list(survey = df.sprat, catch = Catch.obs)
##Z <- Z0 + Fin + par$Mpred


par$dt <- 1 # Decrease later 
par$tEnd <- nyear

Amat <- matrix(0,length(w))
Z0 <- par$mu0prefactor*par$wInf^par$mu0exponent

Amat[2:length(w)] <- -par$gg[(2:length(par$w)-1)]*par$dt/par$dw[2:length(par$w)]

phiMat <- (1+(w/(par$wInf*par$alphaMature))^(-10))^(-1)
co <- 0.1 
phiMat[w<(co*par$alphaMature*par$wInf)] <- 0

Fsel <- (1+(par$w/(par$nF * par$wInf))^-par$myF)^-1
Fsel[par$w < 0.01] <- 0
yearidx <- seq(1,nyear*(nsize-1), by = (nsize-1))

df.plot <- df.sprat %>%
  group_by(year) %>%
  summarise(Survey.bio = sum(surveybiomass))

plot(1:nyear,df.plot$Survey.bio/mean(df.plot$Survey.bio), type = 'l')
lines(2:(nyear+1),sprat.ass$SSB/mean(sprat.ass$SSB), col = 'red')
#year <- 
#sprat.ass$SSB <- sprat.ass$SSB*1000000 # Convert to gram


# Approximate initial Rmax 
sum(2e12*par$Sequ*par$w*par$dw)/median(sprat.ass$SSB)
# Apprximate initial Q 
median(df.plot$Survey.bio)/median(sprat.ass$SSB)

# Calculate the 'true' Mdev 

Mdev.t <- matrix(NA, nyear)
Mtest <- matrix(NA, nyear, length(par$w))

Mdev.t[1] <- M.obs$M[1]/(a*par$A*w^(par$n-1))[50]

Mtest[1,] <- (a*par$A*w^(par$n-1))*Mdev.t[1]

for(i in 2:nyear){
  
  Mdev.t[i] <-  M.obs$M[i]/Mtest[i-1,50]
  Mtest[i,] <- Mtest[i-1,]*Mdev.t[i]
}

#plot(Mtest[,50], type = 'l')
data <-list(      dt = 1,
                  wInf = par$wInf,
                  Zbase = Z0,
                  phiMat = phiMat,
                  iTimeMax = nyear,
                  mbins = bins,
                  yearidx = yearidx-1,
                  wlength = length(par$w),
                  nlength = nsize-1,
                  w = par$w,
                  dw = par$dw,
                  Mpredin = (a*par$A*w^(par$n-1)),
                  gg = par$gg,
                  alphaEgg = par$alpha.egg,
                  Sequ = par$Sequ,
                  surveysel = sel.survey,
                  Amat = Amat,
                  widxM = 50,
                  widxF = 52,    
                  # MDEV = log(Mdev.t),
                  #  logphi = log(1),
                  #                   logSDzero = log(0.5),
                  
                  # Finally add the data
                  # SDlog = log(0.35),
                  logSDR = log(0.6),
                  Catch = sprat.ass$Landings*1000,
                  survey = sprat.optim$survey$CPUE,
                  nobs = length(sprat.optim$survey$CPUE)
                  # Q = log(1e-6)
                  # #lognF = log(0.15),
                  # Fstart = 1.3,
                  
                  
)

U <- matrix(0, 3 , data$iTimeMax)
parameters = list(Rmaxlog = log(1e10),
                  SDlog = log(0.4),#,#,
                  SDlogsurv = log(0.7),
                  #  logSDR = log(0.6),
                  logSDM = log(0.1),
                  logSDF = log(0.1),
                  lognF = log(0.01),
                  Q = log(1e-3),
                  #    logphi = log(1),
                  lFstart = log(1),
                  # logR0 = log(0.2),
                  lognFs = log(0.01),
                  logmuF = log(10),
                  logmuFs = log(6),
                  #  lalphaEgg = log(par$alpha.egg),
                  
                  U = U)#

library(TMB)
compile("runSprat2_TMB.cpp")

dyn.load(dynlib("runSprat2_TMB"))

obj <-MakeADFun(data,parameters,DLL="runSprat2_TMB", random = 'U')

lower <- obj$par-Inf
# lower['lognF'] <- log(0.0001)
# lower['lognFs'] <- log(0.0001)
lower['logmuF'] <- log(1)
lower['logmuFs'] <- log(1)
# #lower['lFstart'] <- log(0.5)
# #lower['logphi'] <- log(0.5)
# 
upper <- obj$par+Inf
#  # 
upper['lognF'] <- log(1)
upper['lognFs'] <- log(1)
#upper['SDlog'] <- log(0.1)
#upper['logSDM'] <- log(0.2)
#  #upper['R0'] <- log(2)
# upper['logmuF'] <- log(50)
# upper['logmuFs'] <- log(50)
#upper['logSDR'] <- log(0.6) # COrresponding to the real SR relationship
# upper['Rmaxlog'] <- log(1e2)
#upper['Q'] <- log(1e-4)
#  
system.time(opt<-nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper)) # If error one of the random effects is unused

rep<-sdreport(obj)

rep

#Uncertainty 
sdrep <- summary(rep)
rep.values<-rownames(sdrep)


unique(rep.values)

getUncertainty <- function(name){
  
  df <- data.frame(name = sdrep[rep.values == name,1])
  df$SE <- sdrep[rep.values == name,2]
  df$min <- df$name-2*df$SE
  df$max <- df$name+2*df$SE
  
  if(dim(df)[1] == nyear){
    df$year <- year
  }
  
  return(df)  
}

Catch <- getUncertainty('Catchest')
Fest <- getUncertainty('Fsave')
Bio <- getUncertainty('Bio')
Msave <- getUncertainty('Msave')
Fsel <- getUncertainty('Fsel')
Rsave <- getUncertainty('R')
SSB <- getUncertainty('SSB')
Nest <- getUncertainty('Nest')
Nest$weight <- rep(m.bins, nyear)
Nsave <- getUncertainty('Nsave')

Nsave$w <- rep(data$w, nyear)

source('plotFit.R')

plotFit(plotexp = 0)

dev.off()
plot(sprat.ass$Year,sprat.ass$Landings*1000, type ='l', lwd = 2)
#plot(Catch.obs$Year,Catch.obs$Catch)
lines(year,Catch$name, col = 'red')
lines(year,Catch$min, lty = 16, col = 'red')
lines(year,Catch$max, lty = 16, col = 'red')

exp(rep$par.fixed)

# Plot size fit for appendix 
#jpeg(filename = )
dev.off()

par(mfrow = c(ceiling(sqrt(nsize)),ceiling(sqrt(nsize))), mar = c(1,3,0,0), oma = c(4,4,1,1))

for (i in 1:data$nlength){
  Ntmp <- Nest$name[Nest$weight == sprat.mean$weight[i]]
  Ntmp.min <- Nest$min[Nest$weight == sprat.mean$weight[i]]
  Ntmp.max <- Nest$max[Nest$weight == sprat.mean$weight[i]]
  
  stmp <- data$survey[yearidx+i-1]
  
  yl = c(min(Ntmp[Ntmp > 0],stmp[stmp >0]), max(Ntmp,stmp))
  
  plot(year,Ntmp, type = 'l',  ylim = yl, log = 'y', lwd = 2)
  
  lines(year, Ntmp.min, lty = 2)
  lines(year,Ntmp.max, lty = 2)
  points(year,stmp, col = 'red')
  lines(year,stmp, col = alpha('red', alpha = 0.2), lty = 16)
  
}

library(scales)
# size spectrum each year 
png(filename = 'C:/Users/Nis/Dropbox/UW/Fluctuations paper/Multispecies model/single species model/single species operating model/Pub figures/sizespectrum.png', 
    width = 16*0.8, height = 12*0.8, 
    unit = 'cm', res = 400)

yidx <- c(yearidx, yearidx[length(yearidx)]+length(sprat.mean$weight))
par(mfrow = c(6,7), mar = c(1,1,1,1), oma = c(3,3,0.5,0.5))

for (i in 1:nyear){
  Ntmp <- Nest$name[yidx[i]:(yidx[i+1]-1)]
  Ntemp.min <- Nest$min[yidx[i]:(yidx[i+1]-1)]
  Ntemp.max <- Nest$max[yidx[i]:(yidx[i+1]-1)]
  
  stmp <- data$survey[yidx[i]:(yidx[i+1]-1)]
  
  yl = c(min(Ntmp[Ntmp > 0],stmp[stmp >0])*0.8, max(Ntmp,stmp)*1.5)
  
  plot(sprat.mean$weight,Ntmp, type = 'l',  ylim = yl, log = 'xy', lwd = 1.5)
  points(sprat.mean$weight,stmp, col = 'red')
  # lines(sprat.mean$weight,Ntemp.min, lty = 2)
  # lines(sprat.mean$weight,Ntemp.max, lty = 2)
  
  polygon(c(sprat.mean$weight, rev(sprat.mean$weight)), 
          c(Ntemp.min, rev(Ntemp.max)), 
          border = NA, col = alpha("gray", alpha = 0.5))
  
  
  if (i == 22){
    mtext(text = 'size spectrum (#/g)',side = 2, padj = -2.4)
    
  }
  
  if(i == nyear){
    mtext(text = 'weight (g)', side  = 1, padj = 2.5)
  }
  
}
#

dev.off()  

par(mar = c(1,1,1,1), oma = c(3,3,1,1))
plot(sprat.ass$Year,sprat.ass$Landings*1000/1e6, lwd = 2, type = 'o', col = 'red')
#plot(Catch.obs$Year,Catch.obs$Catch)
lines(year,Catch$name/1e6, col = 'black', lwd = 2)
polygon(c(sprat.ass$Year, rev(sprat.ass$Year)), 
        c(Catch$min/1e6, rev(Catch$max/1e6)), 
        border = NA, col = alpha("gray", alpha = 0.5))
mtext(text = 'Landings (1000 t)', side =2 , padj = -3.2)
mtext(text = 'Year', side = 1, padj = 3)


#plotFit(F)

# Plot the deviations and see how they are
idx <- seq(1,(length(rep$par.random)-2), by = 3)

Fdev <- rep$par.random[idx]
Mdev <- rep$par.random[idx+1]
Rdev <- rep$par.random[idx+2]

hist(exp(Mdev))
hist(exp(Fdev))
hist(exp(Rdev))


# Recalculate mortality 
Mchange <- matrix(NA, nyear)
Mchange[1] <- data$Mpredin[50]*exp(Mdev[1])

for (i in 2:nyear){
  Mchange[i] <- Mchange[i-1]*exp(Mdev[i])
}

plot(Mchange, type = 'l')







