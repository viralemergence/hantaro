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
library(caper)
library(ape)
library(phylofactor)
library(treeio)
library(ggtree)
library(plotrix)
library(rstatix)
library(ggrepel)
library(ggpubr)

## load files
setwd("~/Desktop/hantaro/data/clean files")
pcr_brts=readRDS("pcr brts.rds")
comp_brts=readRDS("comp brts.rds")
pm_brts=readRDS("pm brts.rds")

## index non-missing
pcr_keep=which(!is.na(sapply(pcr_brts,function(x) x$testAUC)))
comp_keep=which(!is.na(sapply(comp_brts,function(x) x$testAUC)))

## all
keep=intersect(pcr_keep,comp_keep)

## trim
pcr_brts=pcr_brts[keep]
comp_brts=comp_brts[keep]

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
                 response=c(rep('infection',n),rep('competence',n)),
                 seed=c(sapply(pcr_brts,function(x) x$seed),
                        sapply(comp_brts,function(x) x$seed)))
rm(n)

## factor
adata$response=factor(adata$response,levels=c('infection','competence'))

## long to wide
adata2=spread(adata,response,auc)

## difference
adata2$diff=adata2$competence-adata2$infection

## paired t test
summary(lm(diff~1,data=adata2))

## formal test for reporting
t.test(adata2$infection,adata2$competence,paired=T)

## t-test
t.test(auc~response,data=adata,
       alternative='two.sided',
       var.equal=F,paired=T)  

## effect size
cohens_d(auc~response,data=adata,paired=T)

## make jitter position
adata$x=as.numeric(factor(adata$response))
set.seed(1)
adata$xj=jitter(adata$x,0.5)

## dataset of means
amean=data.frame(mean=tapply(adata$auc,adata$response,mean))
amean$se=tapply(adata$auc,adata$response,se)
amean$lower=amean$mean-(1.96*amean$se)
amean$upper=amean$mean+(1.96*amean$se)
amean$response=rownames(amean)

## x
amean$response=factor(amean$response,levels=levels(adata$response))
amean$x=as.numeric(factor(amean$response))

## plot with segments
set.seed(1)
f2A=ggplot(adata)+
  #geom_violin(aes(x=x,y=auc,group=x),trim=T,scale="count",width=0.5)+
  geom_boxplot(aes(x=x,y=auc,group=x),width=0.25,alpha=0.25)+
  geom_line(aes(x=xj,y=auc,group=seed),alpha=0.25)+
  geom_point(aes(x=xj,y=auc),size=1.5,alpha=1)+
  scale_x_continuous(breaks=c(1,2),
                     labels=levels(adata$response),
                     limits=c(0.5,2.5))+
  theme_bw()+
  labs(x="response variable",
       y="model performance (AUC)")+
  theme(axis.text=element_text(size=10),
        axis.text.x=element_text(size=12),
        axis.title=element_text(size=12))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
  guides(colour=F)#+
  
  ## add mean
  #geom_line(data=amean,aes(x=x,y=mean),size=1.5)+
  #geom_segment(data=amean,aes(x=x,xend=x,y=lower,yend=upper),size=1.5)+
  #geom_point(data=amean,aes(x=x,y=mean),size=3,shape=15)

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

## figure 2A
set.seed(1)
# f2A=ggplot(adata,aes(response,auc))+
#   geom_boxplot(width=0.5,alpha=0.25,colour=col,fill=col)+
#   geom_jitter(width=0.1,colour=col,size=3,alpha=1)+
#   scale_x_discrete(labels=c("infection","competence"))+
#   guides(colour=F)+
#   ylim(0.8,1)+
#   theme_bw()+
#   labs(x="response variable",
#        y="BRT AUC")+
#   theme(axis.text=element_text(size=10),
#         axis.text.x=element_text(size=12),
#         axis.title=element_text(size=12))+
#   theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
#   theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
#   theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
#   guides(colour=F,fill=F)

## identify features with high residuals
ranks2$resid=abs(resid(lm(comp_rank~pcr_rank,data=ranks2)))

## flag if resid>10
ranks2$select=ifelse(ranks2$resid>10,"yes","no")

## flag if consistently low or consistently high
n=7
ranks2$select=ifelse(ranks2$comp_rank<n & ranks2$pcr_rank<n,"yes",ranks2$select)
ranks2$select=ifelse(ranks2$comp_rank%in%tail(1:nrow(ranks2),n) & 
                       ranks2$pcr_rank%in%tail(1:nrow(ranks2),n),"yes",ranks2$select)

## flag if high or low ranks
# rnk=c(head(ranks2$comp_rank,n),tail(ranks2$comp_rank,n))
# ranks2$select=ifelse(ranks2$comp_rank%in%rnk,"yes",ranks2$select)

## just yes
rset=ranks2[ranks2$select=="yes",]

## rset as ranks2
rset=ranks2
rset$var=ifelse(rset$select=="yes",rset$var,"")

## figure 2B
set.seed(1)
f2B=ggplot(ranks2,aes(pcr_rank,comp_rank))+
  #geom_label(data=rset,aes(label=var),size=2,fill=col,alpha=0.2)+
  geom_text_repel(data=rset,aes(label=var),
                  size=2,
                  force=4,
                  #nudge_y=-2,
                  #nudge_x=1,
                  direction="both",
                  segment.size=0.5,
                  segment.color="grey")+
  geom_point()+
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
setwd("~/Desktop/hantaro/figs")
png("Figure 2.png",width=7,height=3.5,units="in",res=300)
ggarrange(f2A,f2B,ncol=2,widths=c(1,1),
          labels=c("(A)","(B)"),
          label.x=c(0.21,0.18),
          label.y=0.97,
          font.label=list(face="plain",size=12))
dev.off()

## pdp
detach("package:purrr", unload=TRUE)
library(pdp)
library(gbm)

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
                   size=1,colour="grey",alpha=0.25)+
      
      ## add lines
      geom_line(size=1,alpha=0.25,colour="grey")+
      
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
      geom_jitter(size=1,alpha=0.25,colour="grey",width=0.1)+
      
      ## add mean
      geom_point(data=pmean,size=2,inherit.aes=F,shape=15,
                 aes(predictor,yhat))+
      
      ## add rug
      geom_rug(data=temp,inherit.aes=F,
               aes(predictor,yhat),
               sides="b",position="jitter",
               colour="grey",alpha=0.25,
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

## clean
rm(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,
   c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,
   ranks,ranks2,ts2,f2A,f2B,adata)

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

## add study
apreds=merge(apreds,data[c("treename","studies")],by="treename")

## make positivity
apreds$positivity=ifelse(apreds$hPCR==1 & apreds$type=="PCR",1,
       ifelse(apreds$competence==1 & apreds$type=='competence',1,0))

## make type
apreds$cat=ifelse(apreds$studies==0,"unsampled",
                  ifelse(apreds$positivity==1,"positive","negative"))

## type
library(plyr)
apreds$type=factor(apreds$type,levels=c("PCR","competence"))
apreds$type2=revalue(apreds$type,c("PCR"="infection"))

## long to wide
apreds2=spread(apreds[c('treename','type','cpred')],type,cpred)
comp_apreds$comp=comp_apreds$competence

## merge
apreds2=merge(apreds2,comp_apreds[c("treename","hPCR","comp")],by="treename")

## fix names
names(apreds2)=c("treename","pred_pcr","pred_comp","PCR","competence")

## classify true negatives
data$type=ifelse(data$studies>0 & data$hPCR==0 & data$competence==0,"true negative","other")

## with data
apreds2=merge(apreds2,data[c("treename","type",'studies',"fam","gen")],by='treename')

## fix type
apreds2$cat=ifelse(apreds2$studies==0,"unsampled",
                   ifelse(apreds2$PCR==0 & apreds2$competence==0,"negative","positive"))

## fix cat
apreds2$cat=factor(apreds2$cat,c("positive",'negative','unsampled'))
apreds$cat=factor(apreds$cat,levels=levels(apreds2$cat))

## figure 3a
library(awtools)
cc=mpalette[2:4] 
cc=rev(cc)
f3A=ggplot(apreds,aes(cpred))+
  geom_density(aes(fill=cat,colour=cat),alpha=0.5)+
  facet_wrap(~type2,ncol=1,strip.position='top',scales="free_y")+
  theme_bw()+
  theme(legend.position="top")+
  labs(x=expression(paste("predicted probability (",italic(P),") of hosting")))+
  xlim(0,1)+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        legend.title=element_text(size=12),
        legend.text=element_text(size=11),
        strip.text=element_text(size=11),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(20,20,20,20))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
  scale_colour_manual(values=cc)+
  scale_fill_manual(values=cc)+
  guides(colour=guide_legend(title="(A) orthohantavirus positivity"),
         fill=guide_legend(title="(A) orthohantavirus positivity"))
f3A

## scatterplot
f3B=ggplot(apreds2,aes(pred_pcr,pred_comp))+
  geom_point(alpha=0.5,size=2,aes(colour=cat,fill=cat))+
  geom_smooth(method='gam',colour="grey")+
  labs(x=expression(paste(italic(P),' from infection models')),
       y=expression(paste(italic(P),' from competence models')))+
  theme_bw()+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        legend.title=element_text(size=12),
        legend.text=element_text(size=11),
        strip.text=element_text(size=11),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(20,20,20,20))+
  theme(legend.position="top")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
  scale_colour_manual(values=cc)+
  scale_fill_manual(values=cc)+
  guides(colour=guide_legend(title="(A) orthohantavirus positivity"),
         fill=guide_legend(title="(A) orthohantavirus positivity"))

## combine
setwd("~/Desktop/hantaro/figs")
png("Figure 3.png",width=6.5,height=4,units="in",res=300)
f3=ggarrange(f3A,f3B,common.legend=T)
f3
dev.off()

## save
preds=apreds2
preds$fam=NULL
preds$gen=NULL

## write file
setwd("~/Desktop/hantaro/data/clean files")
write.csv(preds,"hantaro predictions.csv")

## test correlation
cor(apreds2$pred_pcr,apreds2$pred_comp,method='spearman')
cor.test(apreds2$pred_pcr,apreds2$pred_comp,method='spearman')

## load phylogeny
setwd("~/Desktop/hantaro/data/clean files")
rtree=readRDS('rodent phylo trim.rds')

## setdiff
apreds2$tree=ifelse(apreds2$treename%in%setdiff(apreds2$treename,rtree$tip.label),'cut','keep')
table(apreds2$tree)

## trim
bdata=apreds2[-which(apreds2$tree=='cut'),]

## match
bdata=bdata[match(rtree$tip.label,bdata$treename),]

## save
bdata$label=bdata$treename
bdata$Species=bdata$treename

## merge
cdata=comparative.data(phy=rtree,data=bdata,names.col=treename,vcv=T,na.omit=F,warn.dropped=T)

## fix
cdata$data$tree=NULL

## lambda
pcr_lmod=pgls(pred_pcr~1,data=cdata,lambda="ML")
comp_lmod=pgls(pred_comp~1,data=cdata,lambda="ML")
summary(pcr_lmod)
summary(comp_lmod)
# moderate phylogenetic signal in predictions

## taxonomy
cdata$data$taxonomy=paste(cdata$data$fam,cdata$data$gen,cdata$data$Species,sep='; ')

## set taxonomy
taxonomy=data.frame(cdata$data$taxonomy)
names(taxonomy)="taxonomy"
taxonomy$Species=rownames(cdata$data)
taxonomy=taxonomy[c("Species","taxonomy")]
taxonomy$taxonomy=as.character(taxonomy$taxonomy)

## Holm rejection procedure
HolmProcedure <- function(pf,FWER=0.05){
  
  ## get split variable
  cs=names(coef(pf$models[[1]]))[-1]
  split=ifelse(length(cs)>1,cs[3],cs[1])
  
  ## obtain p values
  if (pf$models[[1]]$family$family%in%c('gaussian',"Gamma","quasipoisson")){
    pvals <- sapply(pf$models,FUN=function(fit) summary(fit)$coefficients[split,'Pr(>|t|)'])
  } else {
    pvals <- sapply(pf$models,FUN=function(fit) summary(fit)$coefficients[split,'Pr(>|z|)'])
  }
  D <- length(pf$tree$tip.label)
  
  ## this is the line for Holm's sequentially rejective cutoff
  keepers <- pvals<=(FWER/(2*D-3 - 2*(0:(pf$nfactors-1))))
  
  
  if (!all(keepers)){
    nfactors <- min(which(!keepers))-1
  } else {
    nfactors <- pf$nfactors
  }
  return(nfactors)
}

## get species in a clade
cladeget=function(pf,factor){
  spp=pf$tree$tip.label[pf$groups[[factor]][[1]]]
  return(spp)
}

## summarize pf object 
pfsum=function(pf){
  
  ## get formula
  chars=as.character(pf$frmla.phylo)[-1]
  
  ## response
  resp=chars[1]
  
  ## fix
  resp=ifelse(resp=='cbind(pos, neg)','prevalence',resp)
  
  ## holm
  hp=HolmProcedure(pf)
  
  ## save model
  model=chars[2]
  
  ## set key
  setkey(pf$Data,'Species')
  
  ## make data
  dat=data.frame(pf$Data)
  
  ## make clade columns in data
  for(i in 1:hp){
    
    dat[,paste0(resp,'_pf',i)]=ifelse(dat$Species%in%cladeget(pf,i),'factor','other')
    
  }
  
  ## make data frame to store taxa name, response, mean, and other
  results=data.frame(matrix(ncol=6, nrow = hp))
  colnames(results)=c('factor','taxa','tips','node',"clade",'other')
  
  ## set taxonomy
  taxonomy=dat[c('Species','taxonomy')]
  taxonomy$taxonomy=as.character(taxonomy$taxonomy)
  
  ## loop
  for(i in 1:hp){
    
    ## get taxa
    tx=pf.taxa(pf,taxonomy,factor=i)$group1
    
    ## get tail
    tx=sapply(strsplit(tx,'; '),function(x) tail(x,1))
    
    ## combine
    tx=paste(tx,collapse=', ')
    
    # save
    results[i,'factor']=i
    results[i,'taxa']=tx
    
    ## get node
    tips=cladeget(pf,i)
    node=ggtree::MRCA(pf$tree,tips)
    results[i,'tips']=length(tips)
    results[i,'node']=ifelse(is.null(node) & length(tips)==1,'species',
                             ifelse(is.null(node) & length(tips)!=1,NA,node))
    
    ## get means
    ms=(tapply(dat[,resp],dat[,paste0(resp,'_pf',i)],mean))
    
    ## add in
    results[i,'clade']=ms['factor']
    results[i,'other']=ms['other']
    
  }
  
  ## return
  return(list(set=dat,results=results))
}

## fix palette
AlberPalettes <- c("YlGnBu","Reds","BuPu", "PiYG")
AlberColours <- sapply(AlberPalettes, function(a) RColorBrewer::brewer.pal(5, a)[4])
afun=function(x){
  a=AlberColours[1:x]
  return(a)
}

## make low and high
pcols=afun(2)

## pcr predictions
set.seed(1)
pcrpred_pf=gpf(Data=cdata$data,tree=cdata$phy,
            frmla.phylo=pred_pcr~phylo,
            family=gaussian,algorithm='phylo',nfactors=5,min.group.size=5)

## comp predictions
set.seed(1)
comppred_pf=gpf(Data=cdata$data,tree=cdata$phy,
               frmla.phylo=pred_comp~phylo,
               family=gaussian,algorithm='phylo',nfactors=6,min.group.size=5)

## summarize
pcrpred_pf_results=pfsum(pcrpred_pf)$results
comppred_pf_results=pfsum(comppred_pf)$results

## add model
pcrpred_pf_results$model="infection"
comppred_pf_results$model="competence"

## bind
predpfs=rbind.data.frame(pcrpred_pf_results,comppred_pf_results)

## round
predpfs$clade=round(predpfs$clade,2)
predpfs$other=round(predpfs$other,2)

## write
setwd("~/Desktop/hantaro/figs")
write.csv(predpfs,"Table S3.csv")

## combine tree and data
dtree=treeio::full_join(as.treedata(cdata$phy),cdata$data,by="label")

## plot base tree
pbase=ggtree(dtree,layout="fan",branch.length="none",size=0.25)

## get tree data
tdata=pbase$data

## tips only
tdata=tdata[which(tdata$isTip==T),]

## set x max
xmax=max(tdata$x)+10

## make data frame
samp=data.frame(x=tdata$x,
                y=tdata$y,
                yend=tdata$y,
                xend_pcr=rescale(cdata$data$pred_pcr,c(max(tdata$x),xmax)),
                xend_comp=rescale(cdata$data$pred_comp,c(max(tdata$x),xmax)),
                pred_pcr=(cdata$data$pred_pcr),
                pred_comp=(cdata$data$pred_comp),
                treename=tdata$label)

## merge with cat
samp=merge(samp,apreds2[c("treename","cat")],by="treename",all.x=T)

## pcr
gg=pbase
for(i in 1:nrow(pcrpred_pf_results)){
  
  gg=gg+
    geom_hilight(node=pcrpred_pf_results$node[i],
                 alpha=ifelse(pcrpred_pf_results$tips[i]/Ntip(cdata$phy)<0.5,0.5,0.25),
                 fill="black")
}

## add preds
p1=gg+
  geom_segment(data=samp,aes(x=x,y=y,xend=xend_pcr,yend=yend,colour=cat),size=0.75)+
  scale_colour_manual(values=c(col,viridis(2,option="E",end=0.8)))+
  scale_fill_manual(values=c(col,viridis(2,option="E",end=0.8)))+
  guides(colour=F)

## competence
gg=pbase
for(i in 1:nrow(comppred_pf_results)){
  
  gg=gg+
    geom_hilight(node=comppred_pf_results$node[i],
                 alpha=ifelse(comppred_pf_results$tips[i]/Ntip(cdata$phy)<0.5,0.5,0.15),
                 fill="black")
}

## add preds
p2=gg+
  geom_segment(data=samp,aes(x=x,y=y,xend=xend_comp,yend=yend,colour=cat),size=0.75)+
  scale_colour_manual(values=c(col,viridis(2,option="E",end=0.8)))+
  scale_fill_manual(values=c(col,viridis(2,option="E",end=0.8)))+
  guides(colour=F)

## combine
f3C=p1+p2
f3C=ggarrange(p1,p2,
              labels=c("(B) infection predictions","(C) competence predictions"),
              label.x=c(-0.03,-0.1),
              label.y=0.1,
              font.label=list(face="plain",size=13))

## revise fig 3
setwd("~/Desktop/hantaro/figs")
png("Figure 3-2.png",width=7,height=7.25,units="in",res=300)
#f3+f3C+plot_layout(nrow=2,heights=c(1.25,1))
ggarrange(f3,f3C,nrow=2,heights=c(1.1,1))
#f3B+f3C+plot_layout(nrow=2,heights=c(1,1.5))
#f3B|(p1/p2)+plot_layout(widths=c(1.5,1))
dev.off()

# ## revise
# f3B+f3C+plot_layout(nrow=2)
