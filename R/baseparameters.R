baseparameters  = function(wInf,Rmax,h,a,A,kappa){
  source('makegrid.R')
  source('calcequ.R')
  param <- list()
  
  epsA <- 0.8
  param$epsA <- epsA
  fc <- 0.2
  a <- a # From kkkalis 2016
  
  #A <- h*0.6*0.4 # 0.6 = par$alpha, 0.4 = par$F0est-fc
  param$A <- A
  beta <- 100;
  param$h <- A/(0.6*0.4)
  
  

  
  nSpecies <- length(wInf);
  param$nSpecies <- nSpecies;   # no$ of species
  param$fGridexp <- 0.2;     # Expansion of grid
  param$nGridPP <- 50;
  param$dt <- 0.2;        # Time step 0.1
  param$tEnd <- 40;       # End time # originally 40
  param$iSave <- 1;        # how often to save results (in iterations)
  
  param$f0est <- 0.6; # equilibrium feeding level, for which h-bar was estimated
  #   kappa <- 0$005;  # The value of kappa that I want to have at equilibrium
  # Used together with f0est for the calculation of gamma
  param$f0 <- 0.6; # Feeding level used for calculation of equ$ spectra$
  
  # Species-specific:
  param$wInf <- t(wInf);#logspace(log10(10), log10(1000000), nSpecies)';
  param$mu0prefactor <- 1 #2      # Natural mortality
  param$mu0exponent <- -1/4
  param$muS0 <- 0;                        # Prefactor for starvation mortality
  param$wMax <- 2*max(param$wInf);
  
  # Growth parameters:
  param$alpha <- 0.6;                    # Assimilation efficiency
  param$n <- 3/4; # 3/4 
  # Scaling of intake      
  param$k <- epsA*A*wInf^(param$n-1)
  param$ks <- fc*param$alpha*h         # Activity
  param$p <- 3/4 #0.75; #param$n;                    # Scaling of std$ metabolism
  #param$ks <- 0.08*param$h; # 2.4;        # std$ metabolism (and activity)
  #param$ks <- 0
  # 
  # if exist('VBk','var')
  # hCalc <- (3*VBk)$/(wInf$^(-1/3)*param$alpha*param$f0)+(param$ks/(param$alpha*param$f0)); 
  # else 
  
  # end
  
  
  param$fc <- param$ks/(param$alpha*param$h);
  #fprintf('Critical feeding level: #f\n',param$fc);
  
  # Encounter:
  param$q <- 0.8;                        # Scaling of search volume
  param$beta <- beta;                     # Predation/prey mass ratio
  param$sigma <- 1.3;                            # Width of size preference (1$3)
  lambda <- 2+param$q-param$n;
  alphae <- sqrt(2*pi)*param$sigma*param$beta^(lambda-2)*exp((lambda-2)^2*param$sigma^2/2);
  param$gamma <- (param$f0est*param$h / (alphae*kappa*(1-param$f0est)));
  param$v <- 1*matrix(1,nSpecies)  # Vulnerabilities
  
  
  # Recruitment:
  param$w0 <- 0.001                     # Weigth to recruit at
  param$R0 <- 0/param$wInf                   # "Background" recruitment flux
  param$rho <- param$wInf^(0)                # Egg survival
  param$rho <- param$rho/param$rho[1]
  param$alphaMature <- 0.25          # Fraction of wInf to mature at 0$25
  param$eRepro <- 0.1                # Eff$ of gonad production - 0$1, CC 0$05
  
  # Primary production
  param$typeResource <- 1 # semi-chemostat
  param$rR <- 4
  param$kappaR <- kappa;# *(1/param$f0est - 1) $/ (1$/param$f0 -1); #0$04;#kappa; #Sequ$alphap*Sequ$kappa; # Adjusted upward to give expected f_0
  param$PPmin <- 0.001 
  param$kR <- -2-param$q+param$n; 
  param$lR <- 0.25;
  param$wRcut <- 20;                  # Cut off of background spectrum
  
  
  w <- makegrid(param)[[1]];
  param$kappaPP <- kappa;
  Z0 <- param$mu0prefactor * param$wInf^param$mu0exponent
  
  
  Sequ.t <-  wInf^(-lambda+param$n+a)*w^(-param$n - a)*
    (1 - (w/wInf)^(1-param$n))^((a+wInf^(1-param$n)*as.numeric(Z0)/h)/(1-param$n))
  Sequ.t[is.nan((Sequ.t))] <- 0
  
  param$Sequ <- Sequ.t
  param$Ninit <- Sequ.t*Rmax # Why?
  
  
  # Sequ <- calcequ(param, w);
  # param$a <- Sequ$a;
  # 
  # 
  #param$kappaPP <- kappa/param$rPP*(param$rPP +  kappa*kappa*alphae*alphap/(alphae*kappa+param$h));
  #
  # Fix recruitment:
  param$nRecruitmentType <- 2; # Beverton-Holt
  
  # discount <- 0.01;
  # param$fN0 <- Sequ$Nequ[,1] * discount; # Note, the initial level is discounted with the# expected reduction in the
  # backgroun spectrum
  
  n <- param$n;
  q <- param$q; 
  # a <- param$a[1];
  # param$kappaRmax <- param$kappaPP*1e4
  # 
  # 
  # 
  # param$Rmax <- (param$alpha*param$f0est*param$h*param$w0^param$n-param$ks*param$w0^param$p)*param$fN0
  # 
  # tmpA = param$wInf[1];
  # tmpB = (log10(param$wInf[param$nSpecies])-log10(param$wInf[1]))/(param$nSpecies-1)
  # dW = tmpB*tmpA*10^(tmpB*((1:param$nSpecies)-1))
  # 
  # 
  # if (length(param$h) == 1){
  #   param$Rmax = as.numeric(param$kappaRmax*(param$alpha*param$f0*param$h*param$w0^param$n-param$ks*
  #                                              param$w0^param$p)*param$wInf^(2*param$n-param$q-3+param$a[1])*dW)
  # }else{
  #   param$Rmax = as.numeric(param$kappaRmax*(param$alpha*param$f0*param$h*param$w0^param$n-param$ks*
  #                                              param$w0^param$p)*param$wInf^(2*param$n-param$q-3+param$a)*dW)
  # }
  # 
  # param$fRandomRecruitment <- 0;
  # 
  # Fishing:
  param$F0 <- 0 * matrix(1,nSpecies); # Overall fishing pressure on
  # ind$ species$
  
  
  # Other parameters:
  param$bExtremeSwitching <- 0
  param$bVerbose <- 1
  
  param$expPrice <- 0*0.25;
  param$facCost <- 0*0.0025;
  param$SSBlimit <- 0.0; # limit where there is an extra cost of SSB falls below this
  param$a <- a
  
  # -------------------------------------------------
  # Fishing is  specified as a "trawl" selectivity pattern through
  # F0 and wFstart$ An alternive is to define a function which returns the
  # fishing mortality as a function of (param,iSpecies,w)$
  # -------------------------------------------------
  # ----------------------------------------------------
  # Params for trawling function
  # ----------------------------------------------------- # All from Andersen & Rice 2010
  
  param$c <- 0.01 
  param$myF <- 10
  #param$aF <- 0$2;
  param$nF <- 0.05
  param$gSigma <- 1.5 # 
  param$Rmax <- Rmax
 
  # Example of a function which specifices a trawl fishing pattern:
  # param$funcFishing <- @(param,iSpecies,w) 0*w(w>0$05*param$wInf(iSpecies))
  return(param)
}