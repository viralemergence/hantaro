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
pm_brts=readRDS("pm brts.rds")

## get net AUC
mean(c(sapply(pcr_brts,function(x) x$testAUC),sapply(comp_brts,function(x) x$testAUC)))
se(c(sapply(pcr_brts,function(x) x$testAUC),sapply(comp_brts,function(x) x$testAUC)))

## for cites
mean(sapply(pm_brts,function(x) x$testAUC))
se(sapply(pm_brts,function(x) x$testAUC))

## clean
rm(pm_brts)

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

## ranks for table S2
ts2=ranks
ts2$feature=ts2$var
ts2=ts2[c("feature","pcr_imp","comp_imp","pcr_rank","comp_rank")]

## Table S2
setwd("~/Desktop/hantaro/figs")
write.csv(ts2,"Table S2.csv")

## correlate
cor.test(ranks$pcr_rank,ranks$comp_rank,method="spearman")

## trim to non-zero
ranks2=ranks[-which(ranks$pcr_imp==0 & ranks$comp_imp==0),]

## rerank
ranks2=ranks2[order(ranks2$pcr_imp,decreasing=T),]
ranks2$pcr_rank=1:nrow(ranks2)
ranks2=ranks2[order(ranks2$comp_imp,decreasing=T),]
ranks2$comp_rank=1:nrow(ranks2)
cor.test(ranks2$pcr_rank,ranks2$comp_rank,method="spearman")

## set color
col='grey70'

## figure 2A
set.seed(1)
f2A=ggplot(adata,aes(response,auc))+
  geom_boxplot(width=0.5,alpha=0.25,colour=col,fill=col)+
  geom_jitter(width=0.1,colour=col,size=3,alpha=1)+
  scale_x_discrete(labels=c("infection","competence"))+
  guides(colour=F)+
  ylim(0.8,1)+
  theme_bw()+
  labs(x="response variable",
       y="BRT AUC")+
  theme(axis.text=element_text(size=10),
        axis.text.x=element_text(size=12),
        axis.title=element_text(size=12))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
  guides(colour=F,fill=F)

## figure 2B
f2B=ggplot(ranks2,aes(pcr_rank,comp_rank))+
  #geom_point()+
  geom_label(aes(label=var),size=2,fill=col,alpha=0.2)+
  #scale_y_reverse(limits=c(max(c(ranks$comp_rank,ranks$pcr_rank))+3,0))+
  #scale_x_reverse(limits=c(max(c(ranks$comp_rank,ranks$pcr_rank))+3,0))+
  scale_y_reverse(limits=c(max(c(ranks2$comp_rank,ranks2$pcr_rank))+4,0))+
  scale_x_reverse(limits=c(max(c(ranks2$comp_rank,ranks2$pcr_rank))+4,0))+
  #geom_abline(slope=1,linetype=2,size=0.5)+
  theme_bw()+
  labs(x="feature rank for PCR",
       y="feature rank for competence")+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))

## Figure 2
library(ggpubr)
setwd("~/Desktop/hantaro/figs")
png("Figure 2.png",width=7,height=3.5,units="in",res=300)
ggarrange(f2A,f2B,ncol=2,widths=c(1,1),
          labels=c("(A)","(B)"),
          label.x=c(0.21,0.18),
          label.y=0.97,
          font.label=list(face="plain",size=12))
dev.off()

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