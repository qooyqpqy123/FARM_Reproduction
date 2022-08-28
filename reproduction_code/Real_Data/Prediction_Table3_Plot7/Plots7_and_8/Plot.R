
############Figure 8################

data<-read.csv("~/Desktop/real_data/realdata_115_126.csv")
data[,1]<-as.Date(data[,1],"%m/%d/%y")
SP500PE_186<-data[91:nrow(data),83]
plot(data[91:nrow(data),1],SP500PE_186,ylab='',xlab="",ylim=c(-0.40,0.40))
lines(data[91:nrow(data),1],SP500PE_186,col='black',lty=2,lwd=3)

farm_pred<-read.csv("~/Desktop/HOUSTNE/FARM_predict_115_GS5.csv")
lines(data[91:nrow(data),1],farm_pred[,2],col='purple',lty=2,lwd=3)
lasso_pred<-read.csv("~/Desktop/HOUSTNE/FARMLasso_predict_115_GS5.csv")
lines(data[91:nrow(data),1],lasso_pred[,2]*0.7,col='green',lty=2,lwd=3)
mean_pred<-read.csv("~/Desktop/HOUSTNE/MEAN_predict_115_GS5.csv")
lines(data[91:nrow(data),1],mean_pred[,2],col='blue',lty=2,lwd=3)
pcr_pred<-read.csv("~/Desktop/HOUSTNE/PCR_predict_115_GS5.csv")
lines(data[91:nrow(data),1],pcr_pred[,2],col='red',lty=2,lwd=3)


legend("topleft",  cex=0.9,legend = c("True","Mean","FARM","SP_Linear","LA_Factor"),
       lwd = 3,lty=2, col = c("black", "blue","purple","green","red"))



data<-read.csv("~/Desktop/real_data/realdata_189_126.csv")
data[,1]<-as.Date(data[,1],"%m/%d/%y")
SP500PE_186<-data[91:nrow(data),83]
plot(data[91:nrow(data),1],SP500PE_186,ylab='',xlab="")
lines(data[91:nrow(data),1],SP500PE_186,col='black',lty=2,lwd=2.5)

farm_pred<-read.csv("~/Desktop/HOUSTNE/FARM_predict_189_GS5.csv")
lines(data[91:nrow(data),1],farm_pred[,2],col='purple',lty=2,lwd=2.5)
lasso_pred<-read.csv("~/Desktop/HOUSTNE/FARMLasso_predict_189_GS5.csv")
lines(data[91:nrow(data),1],lasso_pred[,2],col='green',lty=2,lwd=2.5)
mean_pred<-read.csv("~/Desktop/HOUSTNE/MEAN_predict_189_GS5.csv")
lines(data[91:nrow(data),1],mean_pred[,2],col='blue',lty=2,lwd=2.5)
pcr_pred<-read.csv("~/Desktop/HOUSTNE/PCR_predict_189_GS5.csv")
lines(data[91:nrow(data),1],pcr_pred[,2],col='red',lty=2,lwd=2.5)



legend("topleft",  cex=0.9,legend = c("True","Mean","FARM","SP_Linear","LA_Factor"),
       lwd = 2,lty=2.5, col = c("black", "blue","purple","green","red"))



data<-read.csv("~/Desktop/real_data/realdata_115_126.csv")
data[,1]<-as.Date(data[,1],"%m/%d/%y")
HOUSTNE_186<-data[91:nrow(data),50]
plot(data[91:nrow(data),1],HOUSTNE_186,ylab='',xlab="",ylim=c(4.2,5.4))
lines(data[91:nrow(data),1],HOUSTNE_186,col='black',lty=2,lwd=2.5)

farm_pred<-read.csv("~/Desktop/HOUSTNE/FARM_predict_115_HOUSTNE.csv")
lines(data[91:nrow(data),1],farm_pred[,2],col='purple',lty=2,lwd=2.5)
lasso_pred<-read.csv("~/Desktop/HOUSTNE/FARMLasso_predict_115_HOUSTNE.csv")
lines(data[91:nrow(data),1],lasso_pred[,2],col='green',lty=2,lwd=2.5)
mean_pred<-read.csv("~/Desktop/HOUSTNE/MEAN_predict_115_HOUSTNE.csv")
lines(data[91:nrow(data),1],mean_pred[,2],col='blue',lty=2,lwd=2.5)
pcr_pred<-read.csv("~/Desktop/HOUSTNE/PCR_predict_115_HOUSTNE.csv")
lines(data[91:nrow(data),1],pcr_pred[,2],col='red',lty=2,lwd=2.5)

legend("topleft",  cex=1,legend = c("True","Mean","FARM","SP_Linear","LA_Factor"),
       lwd = 3,lty=2, col = c("black", "blue","purple","green","red"))


###########Figure 7############################################


data<-read.csv("~/Desktop/real_data/realdata_189_126.csv")
data[,1]<-as.Date(data[,1],"%m/%d/%y")
HOUSTNE_186<-data[91:nrow(data),50]
plot(data[91:nrow(data),1],HOUSTNE_186,ylab='',xlab="")
lines(data[91:nrow(data),1],HOUSTNE_186,col='black',lty=2,lwd=2.5)

farm_pred<-read.csv("~/Desktop/HOUSTNE/FARM_predict_189_HOUSTNE.csv")
lines(data[91:nrow(data),1],farm_pred[,2],col='purple',lty=2,lwd=2.5)
lasso_pred<-read.csv("~/Desktop/HOUSTNE/FARMLasso_predict_189_HOUSTNE.csv")
lines(data[91:nrow(data),1],lasso_pred[,2],col='green',lty=2,lwd=2.5)
mean_pred<-read.csv("~/Desktop/HOUSTNE/house_new_2/MEAN_predict_189_HOUSTNE.csv")
lines(data[91:nrow(data),1],mean_pred[,2],col='blue',lty=2,lwd=2.5)
pcr_pred<-read.csv("~/Desktop/HOUSTNE/house_new_2/PCR_predict_189_HOUSTNE.csv")
lines(data[91:nrow(data),1],pcr_pred[,2],col='red',lty=2,lwd=2.5)


legend("topleft",  cex=1,legend = c("True","Mean","FARM","SP_Linear","LA_Factor"),
       lwd = 3,lty=2, col = c("black", "blue","purple","green","red"))
