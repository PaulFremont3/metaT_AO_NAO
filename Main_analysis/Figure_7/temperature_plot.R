

env_arctic <- read.table('env_arctic_3.txt', header=T)
env_arctic <- env_arctic[order(env_arctic$T, decreasing = T),]
env_arctic <- env_arctic[!(env_arctic$Station %in% c('143', '146', '149', '191')),]

pdf('Temperature_plot.pdf', width=14)
par(mar=c(10.1, 10.1, 4.1, 2.1))
plot(1:length(env_arctic$Station), env_arctic$T, type='l', lwd=10, col='red', xaxt='n', lty=1, xlab='', ylab='T (°C)', cex.lab=3, cex.axis=3)
axis(side = 1, at = c(1:length(env_arctic$Station)), labels = paste(env_arctic$Station, 'SRF', sep=''), las=2, cex.axis=3)
dev.off()
