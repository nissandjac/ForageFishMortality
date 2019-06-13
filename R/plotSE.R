plotSE <- function(plotExp = FALSE, fig.name = 'SE.jpg'){
  
  if (plotExp == TRUE){
    jpeg(filename = fig.name, width = 16, height = 12, units = 'cm', res = 1200)
    
  }
  
  par(mfrow = c(3,2), mar = c(4,4,1.5,1.5))
  
  w.idx <- 50
  F.idx <- length(OM.mortality$SF$w)
  SE.mortality <- (Msave$name-M.obs)/M.obs
  SE.m.min <- (Msave$min-M.obs)/M.obs
  SE.m.max <- (Msave$max-M.obs)/M.obs
  yl <- c(-2,2)
  
  plot(year,SE.mortality,ylab =  'M2 bias', type = 'l', lty = 1, lwd = 2, ylim = yl, col = 'red')
  lines(year,rep(0,length(year)), col = 'black')
  lines(year,SE.m.min, col = 'red')
  lines(year,SE.m.max, col = 'red')
  polygon(c(year, rev(year)), c(SE.m.max, rev(SE.m.min)), 
          border = NA, col = "#FF000050")
  text(2010,yl[2]*0.95, 'a')
  
  
  SE.F <- (Fest$name-OM.mortality$SF$Fin[,F.idx])/OM.mortality$SF$Fin[,F.idx]
  SE.F.min <-  (Fest$max-OM.mortality$SF$Fin[,F.idx])/OM.mortality$SF$Fin[,F.idx]
  SE.F.max<-  (Fest$min-OM.mortality$SF$Fin[,F.idx])/OM.mortality$SF$Fin[,F.idx]
  
  plot(year,SE.F,ylab =  'F bias', type = 'l', lty = 1, lwd = 2, col = 'red', ylim = yl)
  lines(year,rep(0,length(year)))
  lines(year,SE.F.min, col = 'red')
  lines(year,SE.F.max, col = 'red')
  polygon(c(year, rev(year)), c(SE.F.max, rev(SE.F.min)), 
          border = NA, col = "#FF000050")
  text(2010,yl[2]*0.95, 'b')
  
  
  SE <- (Bio$name-OM.mortality$SF$Biomass)/OM.mortality$SF$Biomass
  SE.B.min <- (Bio$min-OM.mortality$SF$Biomass)/OM.mortality$SF$Biomass
  SE.B.max <- (Bio$max-OM.mortality$SF$Biomass)/OM.mortality$SF$Biomass

  plot(year,SE, ylim = yl, type  = 'l', lwd = 2, col = 'red')
  lines(year,rep(0,nyear))
  lines(year,SE.B.min, col = 'red')
  lines(year,SE.B.max, col = 'red')
  polygon(c(year, rev(year)), c(SE.B.max, rev(SE.B.min)), 
          border = NA, col = "#FF000050")
  text(2010,yl[2]*0.95, 'c')
  
  
  SE.catch <- (Catch$name-OM.mortality$Catch.true)/OM.mortality$Catch.true
  SE.c.min <- (Catch$min-OM.mortality$Catch.true)/OM.mortality$Catch.true
  SE.c.max <- (Catch$max-OM.mortality$Catch.true)/OM.mortality$Catch.true
  
  
  plot(year,SE.catch, ylim = yl, lwd = 2, col = 'red', type = 'l')
  lines(year,rep(0,nyear))
  lines(year,SE.c.min, col = 'red')
  lines(year,SE.c.max, col = 'red')
  polygon(c(year, rev(year)), c(SE.c.max, rev(SE.c.min)), 
          border = NA, col = "#FF000050")
  
  text(2010,yl[2]*0.95, 'd')
  
  if(plotExp == TRUE){
   dev.off() 
  }
  SE.R <- (Rsave$name-OM.mortality$SF$R)/OM.mortality$SF$R
  SE.R.min <- (Rsave$min-OM.mortality$SF$R)/OM.mortality$SF$R
  SE.R.max <- (Rsave$max-OM.mortality$SF$R)/OM.mortality$SF$R
  
  
  plot(year,SE.R, ylim = yl, lwd = 2, col = 'red', type = 'l')
  lines(year,rep(0,nyear))
  lines(year,SE.R.min, col = 'red')
  lines(year,SE.R.max, col = 'red')
  polygon(c(year, rev(year)), c(SE.R.max, rev(SE.R.min)), 
          border = NA, col = "#FF000050")
  
  text(2010,yl[2]*0.95, 'd')
  
  if(plotExp == TRUE){
    dev.off() 
  }
  
   
print(paste('SE mortality = ',100*median(SE.mortality)))  
print(paste('SE F = ',100*median(SE.F)))  
print(paste('SE SSB = ',100*median(SE)))  
print(paste('SE Catch =', 100*median(SE.catch)))
print(paste('SE R = ',100*median(SE.R)))  
  
  
  
}