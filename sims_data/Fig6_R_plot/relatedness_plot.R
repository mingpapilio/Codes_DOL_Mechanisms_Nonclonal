raw<-read.csv("out_porp_10.txt",sep="\t")
raw<- as.matrix(raw)
#
rep= 10
gen= as.numeric(length(raw[,1]))
# 
length.log<- gen* 0.2
temp.log<- matrix(NA, length.log, 2)
data_porp<- matrix(NA, 4, 7)
data_cord<- matrix(NA, 4, 7)
data_cord[,1]= c(0.125, 0.25, 0.5, 1.0)
data_porp[,1]= data_cord[,1]
colnames(data_cord)<-c("R", "mu_n08e0.8", "sd", "mu_n08e1.0", "sd", "mu_n32e0.8","sd")
colnames(data_porp)<- colnames(data_cord)
#
raw<-read.csv("out_porp_10.txt",sep="\t")
raw<- as.matrix(raw)
for(i in 1: length.log){
  temp.log[i, 1]= mean(raw[i+ 0.8*gen, 2:11])
  temp.log[i, 2]= sd(raw[i+ 0.8*gen, 2:11])
}
#
data_porp[4,2]<- mean(temp.log[,1])
data_porp[4,3]<- mean(temp.log[,2])
#
raw<-read.csv("out_s_10.txt",sep="\t")
raw<- as.matrix(raw)
for(i in 1: length.log){
  temp.log[i, 1]= mean(raw[i+ 0.8*gen, 2:11])
  temp.log[i, 2]= sd(raw[i+ 0.8*gen, 2:11])
}
#
data_cord[4,2]<- mean(temp.log[,1])
data_cord[4,3]<- mean(temp.log[,2])
#
library(Hmisc)
#
dev.new()
plot(data_porp[,1], data_porp[,2], ylim=c(0,1), xlim=c(0,1), type="l", col= "skyblue",
     xlab="Relatedness", ylab="Coordination level")
errbar(data_porp[,1], data_porp[,2], data_porp[,2]+data_porp[,3], data_porp[,2]-data_porp[,3], add=T)
points(data_porp[,1], data_porp[,4], col="orange", type="l")
errbar(data_porp[,1], data_porp[,4], data_porp[,4]+data_porp[,5], data_porp[,4]-data_porp[,5], add=T)
##
dev.new()
plot(data_porp[,1], data_porp[,2], ylim=c(0,1), xlim=c(0,1), type="l", col= "skyblue",
     xlab="Relatedness", ylab="Coordination level")
errbar(data_porp[,1], data_porp[,2], data_porp[,2]+data_porp[,3], data_porp[,2]-data_porp[,3], add=T)
points(data_porp[,1], data_porp[,6], col="orange", type="l")
errbar(data_porp[,1], data_porp[,6], data_porp[,6]+data_porp[,7], data_porp[,6]-data_porp[,7], add=T)
##
dev.new()
plot(data_porp[,1], data_porp[,4], ylim=c(0,1), xlim=c(0,1), type="l", col= "skyblue",
     xlab="Relatedness", ylab="Coordination level")
errbar(data_porp[,1], data_porp[,4], data_porp[,4]+data_porp[,5], data_porp[,4]-data_porp[,5], add=T)
points(data_porp[,1], data_porp[,8], col="orange", type="l")
errbar(data_porp[,1], data_porp[,8], data_porp[,8]+data_porp[,9], data_porp[,8]-data_porp[,9], add=T)
##
dev.new()
plot(data_porp[,1], data_porp[,6], ylim=c(0,1), xlim=c(0,1), type="l", col= "skyblue",
     xlab="Relatedness", ylab="Coordination level")
errbar(data_porp[,1], data_porp[,6], data_porp[,6]+data_porp[,7], data_porp[,6]-data_porp[,7], add=T)
points(data_porp[,1], data_porp[,8], col="orange", type="l")
errbar(data_porp[,1], data_porp[,8], data_porp[,8]+data_porp[,9], data_porp[,8]-data_porp[,9], add=T)
##

fwrite(data_porp, file="data_porp.txt")
fwrite(data_cord, file="data_cord.txt")
