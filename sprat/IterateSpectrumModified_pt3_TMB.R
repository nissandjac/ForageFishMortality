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

IterateSpectrum.modified3.tmb<- function(param, S, F0, Mdev){
  
  if (length(S) == 1) {
    S <- list()
  }
  source("makegrid.r")
  source('gradient.r')
  source('fishing.r')
  # ---------------------------------------------------------
  # Set up grid for weights (w)
  # ---------------------------------------------------------
  # 
  # wl <- makegrid(param);
  # w <- wl[[1]]
  # dw <- wl[[2]]
  w <- param$w
  dw <- param$dw
  
  nSave<- 1
  
  
  nGrid <- length(w);
  Mpred <- param$Mpred
  
  # wPP is the primary spectrum grid
  
  # Copy the spectrum in the x-direction if needed:

  #---------------------------------------------------------
  # Set up the species:
  # ---------------------------------------------------------
  nSpecies <- param$nSpecies;
  wend <- length(w)
  #
  # First distribute the parameters over all species and calculate
  # intermediate variable:
  #
  onez <- matrix(1,nSpecies)
  
  # Preallocate some stuff 
  
  psi <- matrix(0,nSpecies,length(w))
  e <- matrix(0,nSpecies,nGrid)   
  gg <- matrix(0,nSpecies,nGrid)
  Fin <- matrix(0,nGrid)
  # Maximum intake:
  
  # Bioenergetics

  Winf <- as.numeric(param$wInf);
  wMature <- Winf*param$alphaMature;
  Z0 <- param$mu0prefactor * Winf^param$mu0exponent; # Changed from Z0 to mu0
  
  # Set up the background spectrum 
  
  
  if (length(S) == 0){
    if (length(param$Ninit) > 0){
      N <- param$Ninit
    }
  }
  
  #
  # Initial spectrum:
  #  Define species specific h and gamma 
  
  param$psi <- (1+(param$w/(param$wInf*param$alphaMature))^(-10))^(-1)
  param$psi[param$w<(0.1*param$alphaMature*Winf)] <- 0
  
  
  
  
  param$F0 <- F0[1]
  
  if (length(param$fishing) > 0){
    type = param$fishing
    Fin <- fishing(param,1,w,type)
  }
  
  
  #-----------------------------------------------------------------------------
  # If there was a previous run, override some of the defaults given above 
  # -----------------------------------------------------------------------------
  if (length(S) > 1){ 
    nPP <- S$nPP[dim(S$nPP)[1],  ]
    N <- S$N[dim(S$N)[1],  ]
  }
  
  # Allocate the variables:
  A <- matrix(0,nSpecies,nGrid);
  B <- matrix(0,nSpecies,nGrid);
  S.m <- matrix(0,nSpecies,nGrid);
  R <- matrix(0,nGrid,1);
  
  
  #---------------------------------------------------------
  # Initialize:
  # ---------------------------------------------------------
  iTimeMax <- param$tEnd/param$dt
  nSave <- floor(param$tEnd/(param$dt))
  Nsave <- array(0,dim = c(nSave, nGrid))
  NtotSave <- matrix(0,nSave, nGrid)
  Rsave <- matrix(0,nSave)
  M2save <- array(0,dim = c(nSave, nGrid))
  Zsave <- array(0,dim = c(nSave,nGrid))
  Rpsave <- matrix(0,nSave)
  Biomass <- matrix(0,nSave)
  Catch <- matrix(0,nSave)
  
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
  phi <- matrix(1,nSave)
  phi[1] <- 1
  
  Mpred <- a*param$A*w^(par$n-1)
  
  
  #---------------------------------------------------------
  # main loop:
  # ---------------------------------------------------------
  for (iTime in 1:iTimeMax){
    
    
    param$F0 <- F0[iTime]
    type = param$fishing
    Fin <- fishing(param,1,w,type)
    
    gg <- param$A*w^param$n*(1-(w/Winf)^(1-param$n)*(param$epsA+(1-param$epsA)*psi))
    gg[w > Winf] <- 0
    
    
    param$gg <- param$A*w^param$n*(1-(param$w/param$wInf)^(1-param$n)*(param$epsA+(1-param$epsA)*param$psi))
    
    #
    # Total mortality:
    # 
    Z <- Z0 + Fin + Mpred
    #
    # Set up matrix for derivative:
    #
    A[idx] <- -gg[(idx-1)]*dt/dw[idx]
    B[idx] <- 1 + gg[idx]*dt/dw[idx] + Z[idx]*dt
    S.m[idx] <- N[idx]
    
    #
    # BC at upstream end (recruitment)
    #
    alpha.egg <- (param$A*param$eRepro*(1-param$epsA)*param$wInf^(1-param$n))/w[1]
    SSB.temp <- sum((1+(w/(param$wInf*param$alphaMature))^(-10))^(-1)*N*w*dw)  # Egg production (mass/time)
    
    Rp <- alpha.egg*SSB.temp
    # Beverton-Holt
    # Add a recruitment error term
    rec.err <- rlnorm(1, meanlog = 0, sdlog =1)
    
    
    R <- (param$Rmax*Rp/(param$Rmax+Rp))#*rec.err
    #R <- SSB.temp/(param$alpha.egg+param$beta.R*SSB.temp)
    # Add physiological recruitment 
    #R <- Rp
    
    B[1] <- 1 + gg[1]*dt/dw[1] + Z[1]*dt
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
    
    if ((iTime %% 1) == 0) {
      iSave <- iTime/1;
      #eSave(iSave,i,:) <- e(i,:);
      gSave[iSave,] <- gg
      Nsave[iSave,] <- N
      Rsave[iSave] <- R; # recruitment
      Biomass[iSave] <- sum(N*w*dw);
      Catch[iSave] <- sum(Fin*N*w*dw)
      M2save[iSave,] <- Mpred
      SSB[iSave] <- SSB.temp
       Rpsave[iSave] <- Rp; # egg production
      Rtotsave[iSave] <- Rtotsave[iSave] + R;
      Fsave[iSave,] <- Fin
      
    }
    # Things that change every year
    
    #
    
    if (iTime < iTimeMax){
      
      phi[iTime+1] <- phi[iTime]*exp(Mdev[iTime+1])
      Mpred <- a*param$A*w^(param$n-1)*phi[iTime+1]
      
    }
   
    
    # Calculate total spectrum:
    #
    Ntot <- N
    #
    # save results:
    #
    if ((iTime %% 1) == 0){ 
      iSave <- iTime/1;
      NtotSave[iSave,] <- Ntot
      Zsave[iSave,] <- Z
    }
  }
  # --------------------------------------------------------------
  # --------------------------------------------------------------
  rm(S)
  S <- list()
  #
  # Numerical parameters:
  #
  S$t <- seq(param$dt,param$tEnd, by = param$dt*1) # The time steps where the
  # solution is saved
  #
  # Grid:
  #
  S$nSpecies <- nSpecies;# No$ of species (trait classes)$
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
  S$Rp <- Rpsave
  #
  # Spectra:
  #
  S$Biomass <- Biomass;  # Total biomass (time,wInf)
  S$Catch <- Catch
  S$SSB <- SSB
  S$Ntot <- NtotSave;    # Community spectrum exclusing the resource spectrum (time,weight)
  S$N <- Nsave;          # Species spectra (time,wInf,weight)
  #
  S$M2 <- M2save
  #                            end
  #     
  return(S)
}
