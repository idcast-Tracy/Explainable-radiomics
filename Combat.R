# -------------- ComBat ---------------- #
#dev.off()
cat("\014"); rm(list = ls()); options(warn = -1); options(digits=3); cc0=0
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)); library(sva)

# ---------------------- Data upload -----------------------------
filename0 = c('Hospital_I_T1.csv', 'Hospital_II_T1.csv',  'Hospital_III_T1.csv')
Groupname = c('Hospital I', 'Hospital II', 'Hospital III'); Batch = c(1, 2, 3)
cc0 = 19; mydata = data.frame(); for (f in c(1:length(filename0))) {ref.data = read.csv(filename0[f]); names(ref.data)[1]<-"p_ID"; ref.data[,1] = as.factor(ref.data[,1]); mydata = rbind(mydata, cbind(Groupname[f], Batch[f], ref.data))}; names(mydata)[1:2] = c('group', 'batch')

mydata = data.frame(mydata[order(mydata$batch), ])
Hospital = unique(mydata$group); Hospital
Color = c("#0bb6c3", "#b394e5", "#efd15a", "#dc5c27", "#fcbea5","#d9af09", "#f17f4f" )

ref.batch = 2; cc = ifelse(cc0>0, cc0+2, cc); V = c()
total_x = mydata[, -c(1:cc)]; combat_x = ComBat(dat=t(total_x), batch=mydata$batch, mod=NULL, ref.batch=ref.batch, par.prior=TRUE, prior.plots=FALSE); total_x <- t(combat_x)
ComBat.data = cbind(mydata[, c(1:cc)], total_x)
if (length(filename0)==1){write.csv(ComBat.data, paste0('ComBat_', filename0), row.names = F)}
if (length(filename0)>=2){for (j in c(1:length(Batch))) {
    write.csv(ComBat.data[which(ComBat.data$batch==j), -c(1:2)], paste0('ComBat_', filename0[j]), row.names = F)}}

library(ggplot2); library(factoextra); library(FactoMineR)
# ---- PCA: Before correction -----------
Bef.Data = mydata
Bef.pca = PCA(Bef.Data[,-c(1:cc)], graph=FALSE)
# ---- PCA: After correction -----------
Af.Data  = ComBat.data
Af.pca  = PCA( Af.Data[,-c(1:cc)], graph=FALSE)

lim = 60
pdf('Before_ComBat.pdf', width=4, height = 4)
fviz_pca_ind(Bef.pca, geom.ind="point", pointsize=3, pointshape=21, 
             fill.ind = as.factor(Bef.Data$group),
             palette = Color, addEllipses=TRUE, title = 'Before batch calibration',
             legend.title="")+ 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
  theme(legend.position='none', plot.title = element_text(hjust = 0.5))+
  xlim(-lim, lim) + ylim(-lim, lim) #+
dev.off()

pdf('After_ComBat.pdf', width=5, height = 4)
fviz_pca_ind(Af.pca, geom.ind="point", pointsize=3, pointshape=21, 
             fill.ind = as.factor(Af.Data$group),
             palette = Color, addEllipses=TRUE, title = 'After batch calibration',legend.title="Group")+ 
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
  theme(plot.title = element_text(hjust = 0.5))+
  xlim(-lim, lim) + ylim(-lim, lim)
dev.off()
print('finish')

