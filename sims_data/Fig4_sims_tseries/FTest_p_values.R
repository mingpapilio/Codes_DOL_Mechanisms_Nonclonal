raw_s<-read.csv("out_s_10.txt",sep="\t")
raw_p<-read.csv("out_porp_10.txt",sep="\t")
#
rep= 10
gen= as.numeric(length(raw_s[,1]))
#
select<- seq(0.9*gen+10, gen, by= 10)
s_sel<- raw_s[select,]
p_sel<- raw_p[select,]
## Stat test between averages of each trial
trial_s_avg<- colMeans(s_sel[,(2:11)])
trial_p_avg<- colMeans(p_sel[,(2:11)])
shapiro.test(trial_s_avg)$p.value ## larger than 0.05 suggests a good fit
shapiro.test(trial_p_avg)$p.value ## larger than 0.05 suggests a good fit
var.test(trial_p_avg,trial_s_avg)$p.value ## smaller than 0.05 suggests significantly different
##
## Stat test through time =====
var.test(s_sel[,2], p_sel[,2])$p.value ## smaller than 0.05 suggests significantly different
var.test(s_sel[,3], p_sel[,3])$p.value ## smaller than 0.05 suggests significantly different
var.test(s_sel[,4], p_sel[,4])$p.value ## smaller than 0.05 suggests significantly different
var.test(s_sel[,5], p_sel[,5])$p.value ## smaller than 0.05 suggests significantly different
var.test(s_sel[,6], p_sel[,6])$p.value ## smaller than 0.05 suggests significantly different
var.test(s_sel[,7], p_sel[,7])$p.value ## smaller than 0.05 suggests significantly different
var.test(s_sel[,8], p_sel[,8])$p.value ## smaller than 0.05 suggests significantly different
var.test(s_sel[,9], p_sel[,9])$p.value ## smaller than 0.05 suggests significantly different
var.test(s_sel[,10], p_sel[,10])$p.value ## smaller than 0.05 suggests significantly different
var.test(s_sel[,11], p_sel[,11])$p.value ## smaller than 0.05 suggests significantly different

## =====
library(lawstat)
levene.test(trial_p_avg,trial_s_avg)
bartlett.test(trial_p_avg,trial_s_avg)
