raw<-read.csv("out_s_10.txt",sep="\t")
#
rep= 10
gen= as.numeric(length(raw[,1]))
#
select<- seq(10, gen, by= 10)
raw_sel<- raw[select,]
#
dev.new()
plot(raw_sel[,1], raw_sel[,2], type="l", col="lightgrey", ylim=c(0,1.0), xlab="generation", ylab="level or proportion")
for(i in 2: rep){
  points(raw_sel[,1], raw_sel[,i+1], type="l", col="lightgrey")
}
mean= rep(NA, length(raw_sel[,1]))
sd= rep(NA, length(mean))
for(i in 1: length(mean)){
  mean[i]= mean(as.matrix(raw_sel[i,2:rep+1]))
  sd[i]= sd(as.matrix(raw_sel[i,2:rep+1]))
}
###
raw<-read.csv("out_porp_10.txt",sep="\t")
raw_sel<- raw[select,]
for(i in 2: rep){
  points(raw_sel[,1], raw_sel[,i+1], type="l", col="pink")
}
mean_cord= rep(NA, length(mean))
sd= rep(NA, length(mean))
for(i in 1: length(mean)){
  mean_cord[i]= mean(as.matrix(raw_sel[i,2:rep+1]))
  sd[i]= sd(as.matrix(raw_sel[i,2:rep+1]))
}
points(raw_sel[,1], mean, type="l",col="black", lwd=3)
points(raw_sel[,1], mean_cord, type="l",col="black", lty=3, lwd=3)
