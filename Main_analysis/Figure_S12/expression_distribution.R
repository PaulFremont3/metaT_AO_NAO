setwd("~/expression_distribution")

dat=read.table('phate_training_expression_GGZZ_0_0_0.txt')
env=read.table('../chl_arctic/env_arctic_4.txt')
env =env[!(c(env$Station %in% c(143, 146,149,191))),]

new_dat=NULL
means=NULL
pdf('distribs_expr.pdf')
for (i in 1:20){
  dt=dat[,i]
  me=mean(dat[dat[,i]!=0,i])
  sde=sd(dat[dat[,i]!=0,i])
  means=c(means, me)
  hist(dat[dat[,i]!=0,i], breaks=1000, main=env$Station[i])
  mtext(paste(round(me, 1), round(sde,1), sep=' '), 3)
  me_b=mean(dat[dat[,i]!=0,i]-me)
  print(sum(dat[dat[,i]!=0,i]-me))
  hist(dat[dat[,i]!=0,i]-me, breaks=1000, main=env$Station[i])
  mtext(paste(round(me_b, 1), round(sde,1), sep=' '), 3)
  
  dt[dt==0]=NA
  dt=dt-me
  dt[is.na(dt)]=0
  new_dat=cbind(new_dat, dt-me)
}
dev.off()

sds=function(x){
  sd(x[x!=0])
}


h=apply(dat,1, sds)
hnew=apply(new_dat,1, sds)

me_sd=mean(h)
sd_sd=sd(h)
me_sdn=mean(hnew)
sd_sdn=sd(hnew)
pdf('sds_bias_and_unbiased.pdf')
hist(h, breaks=1000)
mtext(paste(round(me_sd, 1), round(sd_sd,1), sep=' '), 3)
hist(hnew, breaks = 1000)
mtext(paste(round(me_sdn, 1), round(sd_sdn,1), sep=' '), 3)
dev.off()

pdf('bias-plot.pdf')
plot(1:length(env$Station), means, xaxt='n', pch=19)
axis(1,at = 1:length(env$Station), labels = env$Station, las=2)
abline(h = 0, lty=2)
dev.off()