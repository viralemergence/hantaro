## hantaro 04: processing rodent hantavirus BRT
## danbeck@ou.edu

## clean environment & plots
rm(list=ls()) 
graphics.off()

## libraries
library(tidyr)
library(ggplot2)
library(sciplot)
library(fastDummies)

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

## function for compiling across BRTs for a given predictor, all else equal
pdp_agg=function(mod,feature){
  
  ## just the plot function
  pdep=plot(mod$mod,feature,
            return.grid=T,
            n.trees=mod$best,
            plot=F,
            continuous.resolution=200,
            type="response")
  
  ## add seed
  pdep$seed=unique(mod$roc$seed)
  
  ## save predictor
  pdep$predictor=pdep[feature][,1]
  
  ## order
  pdep=pdep[order(pdep$predictor),]
  
  ## get rank
  pdep$rank=1:nrow(pdep)
  
  ## save yhat
  pdep$yhat=pdep$y
  
  ## return
  return(pdep)
  
}

## function to plot
pdp_plot=function(bmods,feature){
  
  ## pdp_agg
  agg=do.call(rbind,lapply(bmods,function(x) pdp_agg(x,feature)))
  
  ## get class of the feature
  cl=class(data[feature][,1])
  
  ## if else based on type
  if(cl%in%c("numeric","integer")){
    
    ## get element-wise means
    x=with(agg,tapply(predictor,rank,mean))
    y=with(agg,tapply(yhat,rank,mean))
    
    ## save as mean
    pmean=data.frame(predictor=x,yhat=y)
    
    ## get yrange
    yrange=range(agg$yhat,pmean$yhat,na.rm=T)
    
    ## get histogram
    hi=hist(data[feature][,1],breaks=30,plot=F)
    hi=with(hi,data.frame(breaks[1:(length(breaks)-1)],counts))
    names(hi)=c("mids","counts")
    
    ## ggplot it
    ggplot(agg,aes(predictor,yhat,group=seed))+
      
      ## add histogram
      geom_segment(data=hi,inherit.aes=F,
                   aes(x=mids,xend=mids,
                       y=yrange[1],yend=plotrix::rescale(counts,yrange)),
                   size=1,colour=col,alpha=0.25)+
      
      ## add lines
      geom_line(size=1,alpha=0.5,colour=col)+
      
      ## add mean
      geom_line(data=pmean,size=2,inherit.aes=F,
                aes(predictor,yhat))+
      
      ## theme
      theme_bw()+
      theme(axis.text=element_text(size=6),
            axis.title=element_text(size=7))+
      theme(axis.title.x=element_text(margin=margin(t=5,r=0,b=0,l=0)))+
      theme(axis.title.y=element_text(margin=margin(t=0,r=5,b=0,l=0)))+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
      labs(x=feature,y="marginal effect")+
      scale_y_continuous(labels=scales::number_format(accuracy=0.01))
    
    ## end numeric
  }else{ ## factor-based plot
    
    ## get element-wise means
    y=with(agg,tapply(yhat,predictor,mean))
    
    ## save as mean
    #pmean=data.frame(predictor=x,yhat=y)
    pmean=data.frame(y)
    names(pmean)="yhat"
    pmean$predictor=rownames(pmean)
    rownames(pmean)=NULL
    
    ## make temp data
    temp=data
    temp$predictor=temp[feature][,1]
    
    ## do nothing
    agg=agg
    pmean=pmean
    temp=temp
    
    ## get yrange
    yrange=range(agg$yhat,pmean$yhat,na.rm=T)
    
    ## fix temp to yrange
    temp$yhat=ifelse(temp$hPCR==1,max(yrange),min(yrange))
    
    ## ggplot with rug
    set.seed(1)
    ggplot(agg,aes(predictor,yhat,group=seed))+
      
      ## add individual BRTs
      geom_jitter(size=1,alpha=0.5,colour=col,width=0.1)+
      
      ## add mean
      geom_point(data=pmean,size=2,inherit.aes=F,shape=15,
                 aes(predictor,yhat))+
      
      ## add rug
      geom_rug(data=temp,inherit.aes=F,
               aes(predictor,yhat),
               sides="b",position="jitter",
               colour=col,alpha=0.25,
               na.rm=T)+
      
      ## theme
      theme_bw()+
      theme(axis.text=element_text(size=6),
            axis.title=element_text(size=7))+
      theme(axis.title.x=element_text(margin=margin(t=5,r=0,b=0,l=0)))+
      theme(axis.title.y=element_text(margin=margin(t=0,r=5,b=0,l=0)))+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
      labs(x=feature,y="marginal effect")+
      scale_y_continuous(limits=c(yrange[1]-0.01,yrange[2]+0.01),
                         labels=scales::number_format(accuracy=0.01))
    
  }
  
}

## load files
setwd("~/Desktop/hantaro/data/clean files")
data=read.csv('hantaro cleaned response and traits.csv')

## classify true negatives
data$type=ifelse(data$studies>0 & data$hPCR==0 & data$competence==0,"true negative","other")

## make binary columns for genus
dums=dummy_cols(data["gen"])

## unique
dums=dums[!duplicated(dums$gen),]

## ensure all factor
for(i in 1:ncol(dums)){
  
  ## column as factor
  dums[,i]=factor(dums[,i])
  
}

## merge
data=merge(data,dums,by="gen",all.x=T)
rm(dums)

## top PCR
ranks2=ranks2[order(ranks2$pcr_rank),]
p1=pdp_plot(pcr_brts,ranks2$var[1])
p2=pdp_plot(pcr_brts,ranks2$var[2])
p3=pdp_plot(pcr_brts,ranks2$var[3])
p4=pdp_plot(pcr_brts,ranks2$var[4])
p5=pdp_plot(pcr_brts,ranks2$var[5])
p6=pdp_plot(pcr_brts,ranks2$var[6])
p7=pdp_plot(pcr_brts,ranks2$var[7])
p8=pdp_plot(pcr_brts,ranks2$var[8])
p9=pdp_plot(pcr_brts,ranks2$var[9])
p10=pdp_plot(pcr_brts,ranks2$var[10])

## top competence
ranks2=ranks2[order(ranks2$comp_rank),]
c1=pdp_plot(comp_brts,ranks2$var[1])
c2=pdp_plot(comp_brts,ranks2$var[2])
c3=pdp_plot(comp_brts,ranks2$var[3])
c4=pdp_plot(comp_brts,ranks2$var[4])
c5=pdp_plot(comp_brts,ranks2$var[5])
c6=pdp_plot(comp_brts,ranks2$var[6])
c7=pdp_plot(comp_brts,ranks2$var[7])
c8=pdp_plot(comp_brts,ranks2$var[8])
c9=pdp_plot(comp_brts,ranks2$var[9])
c10=pdp_plot(comp_brts,ranks2$var[10])

## compile
library(patchwork)
setwd("~/Desktop/hantaro/figs")
png("Figure S1.png",width=4,height=10,units="in",res=300)
p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+
  c1+c2+c3+c4+c5+c6+c7+c8+c9+c10+plot_layout(nrow=10,ncol=2,byrow=F)
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