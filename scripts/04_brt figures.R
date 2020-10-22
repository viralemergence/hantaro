## hantaro 04: processing rodent hantavirus BRT
## danbeck@ou.edu

## clean environment & plots
rm(list=ls()) 
graphics.off()

## libraries
library(tidyr)
library(ggplot2)
library(sciplot)

## load files
setwd("~/Desktop/hantaro/data/clean files")
pcr_brts=readRDS("pcr brts.rds")
comp_brts=readRDS("comp brts.rds")

## get net AUC
mean(c(sapply(pcr_brts,function(x) x$testAUC),sapply(comp_brts,function(x) x$testAUC)))
se(c(sapply(pcr_brts,function(x) x$testAUC),sapply(comp_brts,function(x) x$testAUC)))

## independent auc
mean(sapply(pcr_brts,function(x) x$testAUC))
se(sapply(pcr_brts,function(x) x$testAUC))
mean(sapply(comp_brts,function(x) x$testAUC))
se(sapply(comp_brts,function(x) x$testAUC))

## t test
n=length(sapply(pcr_brts,function(x) x$testAUC))
adata=data.frame(auc=c(sapply(pcr_brts,function(x) x$testAUC),
                       sapply(comp_brts,function(x) x$testAUC)),
                 response=c(rep('infection',n),rep('competence',n)))
rm(n)

## factor
adata$response=factor(adata$response,levels=c('infection','competence'))

## t test
t.test(auc~response,data=adata,
       alternative='two.sided',
       var.equal=F)                 

## relative importance for PCR
vinf=lapply(pcr_brts,function(x) x$rinf)
pcr_vinf=do.call(rbind,vinf)

## relative importance for competence
vinf=lapply(comp_brts,function(x) x$rinf)
comp_vinf=do.call(rbind,vinf)

## aggregate mean
vdata_pcr=data.frame(aggregate(rel.inf~var,data=pcr_vinf,mean),
                     aggregate(rel.inf~var,data=pcr_vinf,se)["rel.inf"])
names(vdata_pcr)=c("var","rel.inf","rse")
vdata_pcr=vdata_pcr[order(vdata_pcr$rel.inf,decreasing=T),]

## aggregate mean
vdata_comp=data.frame(aggregate(rel.inf~var,data=comp_vinf,mean),
                      aggregate(rel.inf~var,data=comp_vinf,se)["rel.inf"])
names(vdata_comp)=c("var","rel.inf","rse")
vdata_comp=vdata_comp[order(vdata_comp$rel.inf,decreasing=T),]

## clean
rm(pcr_vinf,comp_vinf,vinf)

## rank
vdata_pcr$pcr_rank=1:nrow(vdata_pcr)

## comp
vdata_comp$comp_rank=1:nrow(vdata_comp)

## rel inf
vdata_pcr$pcr_imp=vdata_pcr$rel.inf/100
vdata_comp$comp_imp=vdata_comp$rel.inf/100

## combine ranks
ranks=merge(vdata_pcr[c("var","pcr_rank","pcr_imp")],
            vdata_comp[c("var","comp_rank","comp_imp")],
            by="var")
rm(vdata_comp,vdata_pcr)

## correlate
cor.test(ranks$pcr_rank,ranks$comp_rank,method="spearman")

## trim to non-zero
ranks2=ranks[-which(ranks$pcr_imp==0 & ranks$comp_imp==0),]
cor.test(ranks2$pcr_rank,ranks2$comp_rank,method="spearman")

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

## 