# Run the sprat from the North Sea 

source('load_data.R')

# Plug in some global parameters and values from fishbase 

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

a <- 0.35 # Prior natural mortality

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

# Official assessment uses year 1974 as the first year 
year <- unique(df.sprat$year)
year <- year[year <2013]
nyear <- length(year)

df.sprat <- df.sprat[df.sprat$year < 2013,] # Not anymore catch data

# Find the survey selectivity
sprat.mean <- df.sprat %>%
  group_by(weight) %>%
  summarise(N = mean(CPUE, na.rm = T))

bins <- exp(seq(log(1),log(Winf[1]*1.1), length.out = nsize))
m.bins <- (bins[2:nsize]+bins[1:(nsize-1)])/2
dw <- (bins[2:nsize]-bins[1:(nsize-1)])



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
M.obs <- data.frame(M = Sprat.M[,5], Year = Sprat.M$year) # 12 g fish

Sprat.weight <- colMeans(Sprat.weight[,2:5])*1000 # in G


F.obs <- Sprat.obs$F[Sprat.obs$Age == 3]


sprat.optim <- list(survey = df.sprat, catch = Catch.obs)


# Prepare some calculations for TMB 

par$dt <- 1 # Change for 
par$tEnd <- nyear

Amat <- matrix(0,length(w))
Z0 <- par$mu0prefactor*par$wInf^par$mu0exponent

Amat[2:length(w)] <- -par$gg[(2:length(par$w)-1)]*par$dt/par$dw[2:length(par$w)]

phiMat <- (1+(w/(par$wInf*par$alphaMature))^(-10))^(-1)
co <- 0.1 
phiMat[w<(co*par$alphaMature*par$wInf)] <- 0

yearidx <- seq(1,nyear*(nsize-1), by = (nsize-1))


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
                  Amat = Amat,
                  widxM = 50,
                  widxF = 52,   
                  logSDR = log(0.6),
                  Catch = sprat.ass$Landings*1000,
                  survey = sprat.optim$survey$CPUE,
                  nobs = length(sprat.optim$survey$CPUE)
)

U <- matrix(0, 3 , data$iTimeMax)
parameters = list(Rmaxlog = log(1e10),
                  SDlog = log(0.4),#,#,
                  SDlogsurv = log(0.7),
                  logSDM = log(0.1),
                  logSDF = log(0.1),
                  lognF = log(0.01),
                  Q = log(1e-3),
                  lFstart = log(1),
                  lognFs = log(0.01),
                  logmuF = log(10),
                  logmuFs = log(6),
                  U = U)#

library(TMB)
compile("runSprat_TMB.cpp")

dyn.load(dynlib("runSprat_TMB"))

obj <-MakeADFun(data,parameters,DLL="runSprat_TMB", random = 'U')

lower <- obj$par-Inf

lower['logmuF'] <- log(1)
lower['logmuFs'] <- log(1)

# 
upper <- obj$par+Inf
#  # 
upper['lognF'] <- log(1)
upper['lognFs'] <- log(1)

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
