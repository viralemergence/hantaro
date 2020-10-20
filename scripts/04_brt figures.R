## hantaro 04: processing rodent hantavirus BRT
## danbeck@ou.edu

## clean environment & plots
rm(list=ls()) 
graphics.off()

## libraries
library(tidyr)
library(ggplot2)

## load files
setwd("~/Desktop/hantaro/data/clean files")
pcr_brts=readRDS("pcr brts.rds")
comp_brts=readRDS("comp brts.rds")

## average predictions: PCR
pcr_apreds=lapply(pcr_brts,function(x) x$predict)
pcr_apreds=do.call(rbind,pcr_apreds)

## aggregate
pcr_apreds=data.frame(aggregate(pred~treename,data=pcr_apreds,mean),
                      aggregate(cpred~treename,data=pcr_apreds,mean)['cpred'], ## holding wos constant
                      aggregate(hPCR~treename,data=pcr_apreds,prod)["hPCR"],
                      aggregate(competence~treename,data=pcr_apreds,prod)["competence"])

## type
pcr_apreds$type='PCR'

## average predictions: competence
comp_apreds=lapply(comp_brts,function(x) x$predict)
comp_apreds=do.call(rbind,comp_apreds)

## aggregate
comp_apreds=data.frame(aggregate(pred~treename,data=comp_apreds,mean),
                       aggregate(cpred~treename,data=comp_apreds,mean)['cpred'], ## holding wos constant
                       aggregate(hPCR~treename,data=comp_apreds,prod)["hPCR"],
                       aggregate(competence~treename,data=comp_apreds,prod)["competence"])

## type
comp_apreds$type='competence'

## apreds
apreds=rbind.data.frame(pcr_apreds,comp_apreds)

## long to wide
apreds2=spread(apreds[c('treename','type','cpred')],type,cpred)
comp_apreds$comp=comp_apreds$competence

## merge
apreds2=merge(apreds2,comp_apreds[c("treename","hPCR","comp")],by="treename")

## fix names
names(apreds2)=c("treename","pred_comp","pred_pcr","PCR","competence")

## save
preds=apreds2

## clean
rm(apreds,apreds2,comp_apreds,pcr_apreds)

## load raw data
setwd("~/Desktop/hantaro/data/clean files")
data=read.csv('hantaro cleaned response and traits.csv')

## classify true negatives
data$type=ifelse(data$studies>0 & data$hPCR==0 & data$competence==0,"true negative","other")

## merge into predictions
preds=merge(preds,data[c("treename","type","studies")],by="treename")
rm(data)

## write file
setwd("~/Desktop/hantaro/data/clean files")
write.csv(preds,"hantaro predictions.csv")