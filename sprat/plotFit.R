plotFit <- function(plotexp = FALSE){
  
  if(plotexp == TRUE){
    
    # wd <- "C:/Users/Nis/Dropbox/UW/Fluctuations paper/Multispecies model/single species model/single species operating model/Pub figures/Figure6.jpg"
    # jpeg(filename = wd, width = 16, height = 12, units = 'cm', res = 1200)
    # 
    wd <- "C:/Users/Nis/Dropbox/UW/Fluctuations paper/Documents/ICES estimate mortality/ICES JMS submission/Resubmission/Resubmitted/Figure6.pdf"
    
    cairo_pdf(filename = wd, width = 17*0.393701, height = 12*0.393701)
    
  }
  
  par(mfrow = c(2,2), mar = c(1,3.5,1,1), oma = c(2.5,0.1,0,0))
  
  yl = c(min(Msave$min,M.obs$M), max(Msave$name,M.obs$M))
  yl <- c(0.3,1.2)
  
  M.obs <- M.obs[M.obs$Year<2013,]
  plot(year,Msave$name, ylim = yl, type = 'l', col = 'red4', lwd = 2, 
       ylab = '', xlab = '')
  lines(M.obs$Year,M.obs$M, lwd =2, type = 'o')
  polygon(c(year, rev(year)), c(Msave$max, rev(Msave$min)), 
          border = NA, col = alpha('red', alpha = 0.2))
  text(2010,yl[2]*0.95, 'a')
  mtext(text = expression(paste(mu['M'],' (yr'^{-1},')')),side = 2, padj =-1.5)
  
  yl = c(min(Fest$name,sprat.ass$F0)-0.15, max(sprat.ass$F0,Fest$name)+0.15)
  #yl <- c(0.2,3)
  
  plot(year,Fest$name, type = 'l', col = 'red4', ylim = yl, lwd  = 2,
       ylab = '', xlab = '')
  lines(sprat.ass$Year,sprat.ass$F0, lwd = 2, type = 'o')
  polygon(c(year, rev(year)), c(Fest$max, rev(Fest$min)), 
          border = NA, col = alpha('red', alpha = 0.2))
  text(2010,yl[2]*0.95, 'b')
  mtext(text = expression(paste('F',' (yr'^{-1},')')), side = 2, padj = -1.5)
  
  
  yl = c(min(Bio$name,sprat.ass$SSB*100), max(Bio$name,sprat.ass$SSB*1000))/(1e6)
  #yl <- c(0.2,3)
  
  plot(year,Bio$name/(1e6), type = 'l', ylim = yl, 
       xlab = '', ylab = '', col = alpha('red4', alpha = 1), lwd = 2, log = '')#lines(SF$Biomass)
  lines(sprat.ass$Year,sprat.ass$SSB*1000/(1e6), col = alpha('black', alpha = 0.8), lwd = 2, type = 'o')
  # lines(year,Bio$min/mean(Bio$name))
  # lines(year,Bio$max/mean(Bio$name))
  polygon(c(year, rev(year)), c(Bio$max/(1e6), rev(Bio$min/(1e6))), 
          border = NA,  col = alpha('red', alpha = 0.2))
  text(2010,yl[2]*0.95, 'c')
  mtext(text = 'Biomass (1000 t)', side = 2, padj = -2.7)
  mtext(text = 'Year', side = 1, padj = 2.8)
  # 
  # yl = c(min(Catch$name,Catch.obs$Catch), max(Catch$name,Catch.obs$Catch))
  # 
  # plot(year,Catch$name, ylim = yl, type = 'l', col = 'red', ylab = 'Catch (tonnes)', lwd =2)
  # lines(Catch.obs$Year,Catch.obs$Catch, lwd = 2)
  # polygon(c(year, rev(year)), c(Catch$max, rev(Catch$min)),
  #         border = NA, col = "#FF000050")
  # text(2010,yl[2]*0.95, 'd')


  yl = c(min(Rsave$name/mean(Rsave$name),sprat.ass$Recruitment/mean(sprat.ass$Recruitment))*0.5, 
         max(Rsave$name/mean(Rsave$name),sprat.ass$Recruitment/mean(sprat.ass$Recruitment))*1.4)
  if(yl[1] <0){
    yl[1] <- 0
  }
  
  plot(year,Rsave$name/mean(Rsave$name), ylim = yl, type = 'l', col = 'red4', ylab = '', lwd =2, log = '',
       xlab= '')
  lines(sprat.ass$Year,sprat.ass$Recruitment/mean(sprat.ass$Recruitment), lwd = 2, type = 'o')
  polygon(c(year, rev(year)), c(Rsave$max/mean(Rsave$name), rev(Rsave$min/mean(Rsave$name))),
          border = NA, col = alpha('red', alpha = 0.2))
  text(2010,yl[2]*0.95, 'd')
  mtext('Relative recruitment', side = 2, padj = -2.7)
  mtext(text = 'Year', side = 1, padj = 2.8)
  
  if(plotexp == TRUE){
   dev.off() 
  }
  
}