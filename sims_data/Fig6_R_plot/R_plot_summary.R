library(Hmisc)
library(data.table)
data_porp<- fread("data_porp.txt")
data_cord<- fread("data_cord.txt")
#
dev.new()
plot(data_porp[,R], data_porp[,avg_porp_n], ylim=c(0,1), xlim=c(0,1), type="l", col= "skyblue",
     xlab="Relatedness", ylab="Proportion of helpers")
errbar(data_porp[,R], data_porp[,avg_porp_n], data_porp[,avg_porp_n]+data_porp[,asd_porp_n]/sqrt(10), data_porp[,avg_porp_n]-data_porp[,asd_porp_n]/sqrt(10), add=T)
points(data_porp[,R], data_porp[,avg_porp], type="l")
errbar(data_porp[,R], data_porp[,avg_porp], data_porp[,avg_porp]+data_porp[,asd_porp]/sqrt(10), data_porp[,avg_porp]-data_porp[,asd_porp]/sqrt(10), add=T)
points(data_porp[,R], data_porp[,avg_porp_e], col="orange", type="l")
errbar(data_porp[,R], data_porp[,avg_porp_e], data_porp[,avg_porp_e]+data_porp[,asd_porp_e]/sqrt(10), data_porp[,avg_porp_e]-data_porp[,asd_porp_e]/sqrt(10), add=T)
#
dev.new()
plot(data_cord[,R], data_cord[,avg_s_n], ylim=c(0,1), xlim=c(0,1), type="l", col= "skyblue",
     xlab="Relatedness", ylab="Proportion of helpers")
errbar(data_cord[,R], data_cord[,avg_s_n], data_cord[,avg_s_n]+data_cord[,asd_s_n]/sqrt(10), data_cord[,avg_s_n]-data_cord[,asd_s_n]/sqrt(10), add=T)
points(data_cord[,R], data_cord[,avg_s], type="l")
errbar(data_cord[,R], data_cord[,avg_s], data_cord[,avg_s]+data_cord[,asd_s]/sqrt(10), data_cord[,avg_s]-data_cord[,asd_s]/sqrt(10), add=T)
points(data_cord[,R], data_cord[,avg_s_e], col="orange", type="l")
errbar(data_cord[,R], data_cord[,avg_s_e], data_cord[,avg_s_e]+data_cord[,asd_s_e]/sqrt(10), data_cord[,avg_s_e]-data_cord[,asd_s_e]/sqrt(10), add=T)
