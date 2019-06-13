plotFit = function(plotExp = FALSE, fig.name = 'est_figs.jpg'){


if (plotExp == TRUE){
  jpeg(filename = fig.name, width = 16, height = 12, units = 'cm', res = 300)
  
}

year <- parameters$years

#time.yr <- seq(year[1]+OM$pars$dt, year[length(year)]+1, by = OM$pars$dt)-1
time.yr <- year

par(mfrow = c(2,2), mar = c(4,4,1,1))

yl = c(min(Msave$name,M.obs)-0.05, max(Msave$name,M.obs)+0.1)

plot(year,Msave$name, type = 'l', col = alpha('red', alpha= 1), xlab = 'Year', ylab = 'M', lwd =2, ylim = yl)

points(time.yr,M.obs, type = 'o', col = alpha('black', alpha = 0.6))
lines(year,Msave$max, col = alpha('red', 0.2))
polygon(c(year, rev(year)), c(Msave$max, rev(Msave$min)), 
        border = NA, col = "#FF000050")
lines(year,Msave$min, col = alpha('red', 0.2))
lines(year,Msave$max, col = alpha('red', 0.2))
text(2010,yl[2]*0.95, 'a')

# 
# 
yl = c(min(Fest$name,OM$SF$Fin[,63])-0.1, max(Fest$name,OM$SF$Fin[,63])+0.1)
# if(max(yl) > 2){
#   yl[2] <- 2
# }

plot(year,Fest$name , type = 'l', col = 'red', xlab = 'Year', ylab = 'F', lwd = 2, ylim = yl)
points(time.yr,OM$SF$Fin[,63], col = alpha('black', alpha= 0.6), type = 'o')
lines(year,Fest$max, col = alpha('red', 0.2))
polygon(c(year, rev(year)), c(Fest$max, rev(Fest$min)), 
        border = NA, col = "#FF000050")
lines(year,Fest$min, col = alpha('red', 0.2))
text(2010,yl[2]*0.95, 'b')

yl <- c(min(Bio$name,OM$SF$Biomass)*0.8, max(Bio$name,OM$SF$Biomass)*1.2)

plot(year,Bio$name, type = 'l',  xlab = 'Year', ylab = 'Biomass', col = 'red', lwd = 2, ylim = yl, log = '')
points(time.yr, OM$SF$Biomass, type = 'o', col = alpha('black', alpha= 0.6))
lines(year,Bio$min, col = 'red')
lines(year,Bio$max, col = 'red')
polygon(c(year, rev(year)), c(Bio$min, rev(Bio$max)), 
        border = NA, col = "#FF000050")
text(2010,yl[2]*0.95, 'c')

# plot(year,Catch$name, col = 'red', type = 'l', lwd = 2, xlab = 'Year', ylab = 'Catch', ylim = yl)
# points(year,OM$Catch.true, type = 'o', col = alpha('black', alpha= 0.6))
# lines(year, Catch$min, col = 'red')
# lines(year,Catch$max, col = 'red')
# polygon(c(year, rev(year)), c(Catch$max, rev(Catch$min)), 
#         border = NA, col = "#FF000050")
# 
# lines(rep(year[15],100), seq(-1,1e10,length.out = 100), lty = 15)
# lines(rep(year[40],100), seq(-1,1e10,length.out = 100), lty = 15)
# text(2010,yl[2]*0.95, 'd')
# # 
# Plot recruitment rather than cathc 
yl <- c(min(Rsave$name,OM$SF$R)*0.8, max(Rsave$name,OM$SF$R)*1.1)
plot(year,Rsave$name, col = 'red', type = 'l', lwd = 2, xlab = 'Year', ylab = 'Recruitment', ylim = yl, log = '')
points(time.yr,OM$SF$R, type = 'o', col = alpha('black', alpha= 0.6))
lines(year, Rsave$min, col = 'red')
lines(year,Rsave$max, col = 'red')
polygon(c(year, rev(year)), c(Rsave$max, rev(Rsave$min)),
        border = NA, col = "#FF000050")

text(2010,yl[2]*0.95, 'd')
#
text(2010,yl[2]*0.95, 'd')
# 

if (plotExp == TRUE){
dev.off()
}




}

