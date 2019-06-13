#
# Iterate the spectrum$ 
# Input: param: parameters, see baseparameters$m
#        S: optional parameter with the results of a previous run which is
#           used for initial conditions$
# Output: S: results, see the end of this script$
#
# -------------------------------------------------------------------------
#
#

IterateSpectrum <- function(param, F0, Mdev){
 #
  S <- list()
  
  # Fix fishing mortality and natural mortality to dt 
  
  # tlength <- 1:(param$tEnd/param$dt)
  # 
  # if(length(tlength) > length(F0)){
  #   
  #   F0 <- rep(F0, each = 1/param$dt)
  # }
  # 
  
  source("makegrid.r")
  source('gradient.r')
  source('fishing.r')
  # ---------------------------------------------------------
  # Set up grid for weights (w)
  # ---------------------------------------------------------
  
  wl <- makegrid(param);
  w <- wl[[1]]
  dw <- wl[[2]]
  
  nGrid <- length(w);
  Mpred <- param$Mpred
  
  # wPP is the primary spectrum grid
  #--------------------
  # Set up the species:
  # ---------------------------------------------------------
  wend <- length(w)
  #
  # First distribute the parameters over all species and calculate
  # intermediate variable:
 
  # Preallocate some stuff 

  
  Winf <- as.numeric(param$wInf)
  wMature <- Winf*param$alphaMature;
  Z0 <- param$mu0prefactor * Winf^param$mu0exponent; # Changed from Z0 to mu0
  
  # Set up the background spectrum 
  
  
  #
  # Initial spectrum:
  #  Define species specific h and gamma 
  
  IntakeMax <- param$h * w^param$n;
  
  #
  # Allocation to reproductive mass:
  #
  
  # tmp <- w/param$wInf;
  # psi <- tmp^(1-param$n)*1/(1+(tmp/param$alphaMature)^(-10));     # cut off before maturation
  # # param to mature cutoff
  # 
  # psi[tmp>1] <- 1
  
  psi <- (1+(w/(Winf*param$alphaMature))^(-10))^(-1)
  co <- 0.1 
  psi[w<(co*param$alphaMature*Winf)] <- 0
  param$F0 <- F0[1]
  
  if (length(param$fishing) > 0){
    type = param$fishing
    Fin <- fishing(param,1,w,type)
  }
  
  
  
  #-----------------------------------------------------------------------------
  # If there was a previous run, override some of the defaults given above 
  # -----------------------------------------------------------------------------
  
  N <- param$Sequ*param$Rmax
  # Allocate the variables:
  A <- matrix(0,1,nGrid);
  B <- matrix(0,1,nGrid);
  S.m <- matrix(0,1,nGrid);
  R <- matrix(0,nGrid,1);
  
  
  #---------------------------------------------------------
  # Initialize:
  # ---------------------------------------------------------
  iTimeMax <- param$tEnd/param$dt
  nSave <- floor(param$tEnd/(param$dt*param$iSave))
  Nsave <- array(0,dim = c(nSave, nGrid))
  NtotSave <- matrix(0,nSave, nGrid)
  Rsave <- matrix(0,nSave)
  M2save <- array(0,dim = c(nSave, nGrid))
  Zsave <- array(0,dim = c(nSave,nGrid))
  #Rpsave <- matrix(0,nSave)
  Biomass <- matrix(0,nSave)
  Rtotsave <- matrix(0,nSave,1)
  gSave <- array(0, dim = c(nSave, nGrid))
  MsSave <- array(0, dim = c(nSave, nGrid))
  Fsave <- array(0, dim = c(nSave, nGrid))
  SSB <- matrix(0,nSave)
  
  dt <- param$dt
  # -------------------------------------------------------------------------
  # Calculate total spectrum  (just for the first iteration):
  # -------------------------------------------------------------------------

  Ntot <- sum(N);
  idx <- 2:nGrid

  # -------------------------------------------------------------------------
  # Initialize fishing mortality and natural mortality 
  # -------------------------------------------------------------------------
  
  Z <- matrix(NA,nGrid,iTimeMax)
  Mpred <- matrix(NA,nGrid,iTimeMax)
  Ftot <- matrix(NA,nGrid,iTimeMax)
  a <- param$a
  
  Mpred[,1] <- a*param$A*w^(param$n-1)
  Ftot[,1] <- Fin
  Z[,1] <- param$mu0prefactor*Winf^param$mu0exponent + Ftot[,1] + Mpred[,1]
  # Set up Mortality 
  phi <- 1
  for (i in 2:length(Mdev)){
    phi[i] <- phi[i-1]*exp(Mdev[i])

  }  

  #  if(length(tlength) > length(Mdev)){
  #  phi <- rep(phi, each = 1/param$dt)
  # }

 #  if(length(tlength) > length(Mdev)){
 #  phi <- rep(1, iTimeMax)
 #  ix <- seq(1, iTimeMax, by = 1/param$dt)
 #  phi[ix] <- phitemp
 # }

  for (i in 2:iTimeMax){

    param$F0 <- F0[i]
    type = param$fishing
    Ftot[,i] <- fishing(param,1,w,type)   
    
    Mpred[,i] = a*param$A*w^(param$n-1)*phi[i]
    
    Z[,i] = param$mu0prefactor*Winf^param$mu0exponent + Ftot[,i] + Mpred[,i]
    
      }
  
  
  #rec.err <- rlnorm(1:param$tEnd, meanlog = 0, sdlog =0.6)
  
  # if(length(rec.err) < iTimeMax){
  #   rec.err <- rep(rec.err, each = 1/param$dt)
  # }
  
  #---------------------------------------------------------
  # main loop:
  # ---------------------------------------------------------
  for (iTime in 1:iTimeMax){
    
    g.err <- rnorm(1,0, 0)
    
    gg <- param$A*w^param$n*(1-(w/Winf)^(1-param$n)*(param$epsA+(1-param$epsA)*psi))
    gg[w > Winf] <- 0
    
    gg <- gg*exp(g.err) # Natural variation on growth
    #
    #
    # Set up matrix for derivative:
    #
    A[idx] <- -gg[(idx-1)]*dt/dw[idx]
    B[idx] <- 1 + gg[idx]*dt/dw[idx] + Z[idx,iTime]*dt
    S.m[idx] <- N[idx]
    
    #
    # BC at upstream end (recruitment)
    #
    alpha.egg <- (param$A*param$eRepro*(1-param$epsA)*param$wInf^(1-param$n))/w[1]
    SSB.temp <- sum((1+(w/(Winf*param$alphaMature))^(-10))^(-1)*N*w*dw)  # Egg production (mass/time)
    
    Rp <- alpha.egg*SSB.temp
    # Beverton-Holt
    # Add a recruitment error term
    rec.err <- rlnorm(1, meanlog = 0, sdlog =0.6)
    
    R <- (param$Rmax*Rp/(param$Rmax+Rp))*rec.err
    #R <- SSB.temp/(param$alpha.egg+param$beta.R*SSB.temp)
    # Add physiological recruitment 
    #R <- Rp
    
    B[1] <- 1 + gg[1]*dt/dw[1] + Z[1,iTime]*dt
    N[1] <- (N[1]+R*dt/dw[1])/B[1]
    
    #
    # Invert matrix
    #
    for (j in 2:nGrid){
      N[j] <- (S.m[j]-A[j]*N[j-1])/B[j]
    }
    
    # Set up fishing for the following year

    
    #
    # save results:
    #
    
    if ((iTime %% param$iSave) == 0) {
      iSave <- iTime/param$iSave;
      #eSave(iSave,i,:) <- e(i,:);
      gSave[iSave,] <- gg
      Nsave[iSave,] <- N
      Rsave[iSave] <- R; # recruitment
      Biomass[iSave] <- sum(N*w*dw);
      M2save[iSave,] <- Mpred[,iSave]
      SSB[iSave] <- SSB.temp
     # Rpsave[iSave] <- Rp; # egg production
      Rtotsave[iSave] <- Rtotsave[iSave] + R;
      Fsave[iSave,] <- Ftot[,iSave]
      
    }
      # Things that change every year

    #
    
    # Calculate total spectrum:
    #
    Ntot <- N
    #
    # save results:
    #
    if ((iTime %% param$iSave) == 0){ 
      iSave <- iTime/param$iSave;
      NtotSave[iSave,] <- Ntot
      Zsave[iSave,] <- Z[,iSave]
    }
  }
  # --------------------------------------------------------------
  # --------------------------------------------------------------
  rm(S)
  S <- list()
  #
  # Numerical parameters:
  #
  S$t <- seq(param$dt,param$tEnd, by = param$dt*param$iSave) # The time steps where the
  # solution is saved
  #
  # Grid:
  #
 # S$nSpecies <- nSpecies;# No$ of species (trait classes)$
  S$w <- w;              # Individual weight
  S$dw <- dw;            # Difference between weight classes
  #
  # Species specific rates calculated directly from param:
  #
  S$Fin <- Fsave              # Fishing mortality (wInf,weight)$
  #
  # Reproduction & recruitment:
  #
  S$R <- Rsave;          # Recruitment flux (time,wInf) measured in numbers/time
  #
  # Spectra:
  #
  S$Biomass <- Biomass;  # Total biomass (time,wInf)
  S$SSB <- SSB
  S$Ntot <- NtotSave;    # Community spectrum exclusing the resource spectrum (time,weight)
  S$N <- Nsave;          # Species spectra (time,wInf,weight)
  #
  S$M2 <- M2save
  #                            end
  #     
  return(S)
}
