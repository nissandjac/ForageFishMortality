baseparameters.trimmed  = function(wInf,h,a,A,Rmax,alphaMature){ 
  source('makegrid.R')
  source('calcequ.R')
  param <- list()
  
  epsA <- 0.8
  param$epsA <- epsA
  fc <- 0.2
  a <- a # From kkkalis 2016
  param$a <- a
  #A <- h*0.6*0.4 # 0.6 = par$alpha, 0.4 = par$F0est-fc
  param$A <- A
  param$h <- A/(0.6*0.4)
  
  param$dt <- 0.2;        # Time step 0.1
  param$tEnd <- 40;       # End time # originally 40
  param$wInf <- wInf#logspace(log10(10), log10(1000000), nSpecies)';
  param$mu0prefactor <- 0 #2      # Natural mortality
  param$mu0exponent <- -1/4
  param$n <- 3/4; # 3/4 
  param$w0 <- 0.001                     # Weigth to recruit at
  param$alphaMature <- alphaMature# Fraction of wInf to mature at 0$25
  param$eRepro <- 0.1               # Eff$ of gonad production - 0$1, CC 0$05
  w <- matrix()
  
  i <- 1
  w[1] <- param$w0[[1]]
  
  while (w[i] < 2*param$wInf){
    w[i+1] <- w[i] + 0.2*w[i]
    i <- i + 1
  }
  
  #dw <- gradient(w);
  k <- log(1+0.2);
  dw <- k*param$w0*exp((0:(length(w)-1))*k)
  Z0 <- param$mu0prefactor * param$wInf^param$mu0exponent
  
  param$w <- w
  param$dw <- dw
  #Sequ <- calcequ(param, w);
  lambda <- 2 - param$n + 0.8
  Sequ <-  wInf^(-lambda+param$n+a)*w^(-param$n - a)*
    (1 - (w/wInf)^(1-param$n))^((a+
                                   wInf^(1-param$n)*Z0/h)/(1-param$n))
  Sequ[is.nan((Sequ))] <- 0
  param$Sequ <- Sequ
  
  n <- param$n;
  q <- param$q; 
  
  param$Ninit <- Sequ*Rmax # Why?
  param$Rmax <- Rmax
  # param$Ninit[w>param$alphaMature*param$wInf] <- 0
  
  param$myF <- 10
  #param$aF <- 0$2;
  param$nF <- 0.05
  param$alpha.egg <- (A*param$eRepro*(1-param$epsA)*param$wInf^(1-param$n))/w[1]
  param$psi <- (1+(param$w/(param$wInf*param$alphaMature))^(-10))^(-1)
  param$psi[param$w<(0.1*param$alphaMature*wInf)] <- 0
  param$gg <- param$A*w^param$n*(1-(param$w/param$wInf)^(1-param$n)*(param$epsA+(1-param$epsA)*param$psi))
  param$gg[param$w > param$wInf] <- 0
  
  
  # Pars usually calculated within the loop 
  
  return(param)
  # Save the pars to run 
  
}