## hantaro 03: rodent hantavirus BRT
## danbeck@ou.edu

## clean environment & plots
rm(list=ls()) 
graphics.off()

## libraries
library(gbm)
library(fastDummies)
library(rsample)
library(ROCR)
library(sciplot)
library(ggplot2)
library(pdp)
library(PresenceAbsence)
library(tidyr)
library(viridis)
library(caper)
library(phylofactor)
library(ggtree)
library(treeio)

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

## mode function
mode.prop <- function(x) {
  ux <- unique(x[is.na(x)==FALSE])
  tab <- tabulate(match(na.omit(x), ux))
  max(tab)/length(x[is.na(x)==FALSE])
}

## assess variation across columns
vars=data.frame(apply(data,2,function(x) mode.prop(x)),
                apply(data,2,function(x) length(unique(x))))

## get names
vars$variables=rownames(vars)
names(vars)=c("var","uniq","column")
vars$var=round(vars$var,2)

## if homogenous (100%)
vars$keep=ifelse(vars$var<1,"keep","cut")
vars$keep=ifelse(vars$column%in%c('hPCR','competence',"fam"),'keep',vars$keep)
vars=vars[order(vars$keep),]

## trim
keeps=vars[-which(vars$keep=="cut"),]$column

## drop if no variation
data=data[keeps]
rm(keeps,vars)

## assess missing values
mval=data.frame(apply(data,2,function(x) length(x[!is.na(x)])/nrow(data)))

## get names
mval$variables=rownames(mval)
names(mval)=c("comp","column")

## if have at least 25% values, keep
mval$comp=round(mval$comp,2)
mval$keep=ifelse(mval$comp>=0.25,"keep","cut")
mval=mval[order(mval$keep),]
keeps=mval[-which(mval$keep=="cut"),]$column

## order
mval=mval[order(mval$comp),]

## drop if not well represented
data=data[keeps]
rm(mval,keeps)

## drop unnecessary columns
data$X=NULL
data$superres=NULL
data$Flag=NULL
data$Zdiv=NULL
data$References=NULL
data$WOS_HITS=NULL
data$traitname=NULL
data$Rodents=NULL
data$MSW05_Order=NULL
data$MSW05_Genus=NULL
data$MSW05_Species=NULL

## make simplified
set=data
set$treename=NULL
set$tip=NULL
set$X1.1_ActivityCycle=factor(set$X1.1_ActivityCycle)
set$X12.2_Terrestriality=factor(set$X12.2_Terrestriality)
set$X6.2_TrophicLevel=factor(set$X6.2_TrophicLevel)
set$gen=NULL
set$fam=NULL
set$type=NULL

## remove studies
set$studies=NULL

## covariates
ncol(set)-2

## coverage table s1
ts1=data.frame(apply(set,2,function(x) length(x[!is.na(x)])/nrow(set)))

## get names
ts1$variables=rownames(ts1)
names(ts1)=c("comp","column")
rownames(ts1)=NULL

## sort
ts1=ts1[order(ts1$column),]
ts1$comp=round(ts1$comp,2)

## trim disease data
ts1=ts1[!ts1$column%in%c("hPCR","competence"),]
ts1$feature=ts1$column
ts1$column=NULL
ts1$coverage=ts1$comp
ts1$comp=NULL

## write file
setwd("~/Desktop/hantaro/figs")
write.csv(ts1,"Table S1.csv")

## function to use different data partitions
brt_part=function(seed,response){
  
  ## make new data
  ndata=set
  
  ## correct response
  ndata$response=ndata[response][,1]
  
  ## remove raw
  ndata$hPCR=NULL
  ndata$competence=NULL
  
  ## fix cites if response
  if(response=="cites"){
    
    ## plus 1 for 0
    ndata$cites=ifelse(ndata$cites==0,1,ndata$cites)
    
  }else{
    
    ndata=ndata
    
  }
  
  ## use rsample to split
  set.seed(seed)
  split=initial_split(ndata,prop=0.8,strata="response")
  
  ## test and train
  dataTrain=training(split)
  dataTest=testing(split)
  
  ## yTest and yTrain
  yTrain=dataTrain$response
  yTest=dataTest$response
  
  ## dist
  dist=ifelse(response=="cites","poisson","bernoulli")
  
  ## n.trees
  nt=ifelse(response=="cites",10000,5000)
  
  ## BRT
  set.seed(1)
  gbmOut=gbm(response ~ . ,data=dataTrain,
             n.trees=nt,
             distribution=dist,
             shrinkage=0.001,
             interaction.depth=3,
             n.minobsinnode=4,
             cv.folds=5,class.stratify.cv=TRUE,
             bag.fraction=0.5,train.fraction=1,
             n.cores=1,
             verbose=F)

  ## performance
  par(mfrow=c(1,1),mar=c(4,4,1,1))
  best.iter=gbm.perf(gbmOut,method="cv")
  
  ## predict with test data
  preds=predict(gbmOut,dataTest,n.trees=best.iter,type="response")
  
  ## AUC on train
  auc_train=gbm.roc.area(yTrain,predict(gbmOut,dataTrain,n.trees=best.iter,type="response"))
  
  ## AUC on test
  auc_test=gbm.roc.area(yTest,predict(gbmOut,dataTest,n.trees=best.iter,type="response"))
  
  ## skip if poisson
  if(response=="cites"){
    
    perf=NA
    
  }else{
    
  ## ROC
  pr=prediction(preds,dataTest$response)
  perf=performance(pr,measure="tpr",x.measure="fpr")
  perf=data.frame(perf@x.values,perf@y.values)
  names(perf)=c("fpr","tpr")
  
  ## add seed
  perf$seed=seed
  }
  
  ## relative importance
  bars=summary(gbmOut,n.trees=best.iter,plotit=F)
  bars$rel.inf=round(bars$rel.inf,2)
  
  ## predict with cites
  preds=predict(gbmOut,data,n.trees=best.iter,type="response")
  pred_data=data[c("tip",'treename',"fam","gen","hPCR","competence")]
  pred_data$pred=preds
  pred_data$type=response
  
  ## predict with mean cites
  pdata=data
  pdata$cites=mean(pdata$cites)
  pred_data$cpred=predict(gbmOut,pdata,n.trees=best.iter,type="response")
  
  ## sort
  pred_data=pred_data[order(pred_data$pred,decreasing=T),]
  
  ## print
  print(paste("BRT ",seed," done; test AUC = ",auc_test,sep=""))
  
  ## save outputs
  return(list(mod=gbmOut,
              best=best.iter,
              trainAUC=auc_train,
              testAUC=auc_test,
              roc=perf,
              rinf=bars,
              predict=pred_data,
              traindata=dataTrain,
              testdata=dataTest))
}

## apply across 10 splits each
smax=10
pcr_brts=lapply(1:smax,function(x) brt_part(seed=x,response="hPCR"))
comp_brts=lapply(1:smax,function(x) brt_part(seed=x,response="competence"))

## write to files
setwd("~/Desktop/hantaro/data/clean files")
saveRDS(pcr_brts,"pcr brts.rds")
saveRDS(comp_brts,"comp brts.rds")

## run wos brts
pm_brts=lapply(1:smax,function(x) brt_part(seed=x,response="cites"))

## write
saveRDS(pm_brts,"pm brts.rds")

## average predictions: PCR
pcr_apreds=lapply(pcr_brts,function(x) x$predict)
pcr_apreds=do.call(rbind,pcr_apreds)

## aggregate
pcr_apreds=data.frame(aggregate(pred~treename,data=pcr_apreds,mean),
                      aggregate(cpred~treename,data=pcr_apreds,mean)['cpred'],
                      aggregate(hPCR~treename,data=pcr_apreds,prod)["hPCR"],
                      aggregate(competence~treename,data=pcr_apreds,prod)["competence"])

## type
pcr_apreds$type='PCR'

## average predictions: competence
comp_apreds=lapply(comp_brts,function(x) x$predict)
comp_apreds=do.call(rbind,comp_apreds)

## aggregate
comp_apreds=data.frame(aggregate(pred~treename,data=comp_apreds,mean),
                       aggregate(cpred~treename,data=comp_apreds,mean)['cpred'],
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

## extract AUC and type
adata=data.frame(auc=c(sapply(pcr_brts,function(x) x$testAUC),
                 sapply(comp_brts,function(x) x$testAUC)),
           response=c(rep('PCR',smax),rep('competence',smax)))

## factor
adata$response=factor(adata$response,levels=c('PCR','competence'))

## t test
t.test(auc~response,data=adata,
       alternative='two.sided',
       var.equal=F)                 

## make color
col=rgb(red=125,green=118,blue=170,max=255)

## visualize
p1=ggplot(adata,aes(response,auc))+
  geom_boxplot(width=0.5,alpha=0.25,colour=col,fill=col)+
  geom_jitter(width=0.1,colour=col,size=3,alpha=1)+
  scale_x_discrete(labels=c("infection","competence"))+
  guides(colour=F)+
  ylim(0.84,1)+
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

## check iterations of the tree
table(sapply(pcr_brts,function(x) x$best)<pcr_brts[[1]]$mod$n.trees)
table(sapply(comp_brts,function(x) x$best)<comp_brts[[1]]$mod$n.trees)

## aggregate ROCs
rocs=lapply(pcr_brts,function(x) x$roc)
pcr_rocs=do.call(rbind,rocs)

## for comp
rocs=lapply(comp_brts,function(x) x$roc)
comp_rocs=do.call(rbind,rocs)

## plot
ggplot(pcr_rocs,aes(fpr,tpr,group=seed))+
  theme_bw()+
  geom_line(size=1,alpha=0.5,colour="black")+
  geom_abline(intercept=0,slope=1,size=0.15)+
  
  ## theme
  labs(x="false positive rate",y="true positive rate")+
  theme(axis.title=element_text(size=10),
        axis.text=element_text(size=10))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=5,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=5,b=0,l=0)))
ggplot(comp_rocs,aes(fpr,tpr,group=seed))+
  theme_bw()+
  geom_line(size=1,alpha=0.5,colour="black")+
  geom_abline(intercept=0,slope=1,size=0.15)+
  
  ## theme
  labs(x="false positive rate",y="true positive rate")+
  theme(axis.title=element_text(size=10),
        axis.text=element_text(size=10))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=5,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=5,b=0,l=0)))

## variable importance
vinf=lapply(pcr_brts,function(x) x$rinf)
pcr_vinf=do.call(rbind,vinf)

## comp
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

## type
vdata_pcr$type="PCR"
vdata_comp$type="competence"

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

## correlate
cor.test(ranks$pcr_imp,ranks$comp_imp,method="spearman")

## drop zero
ranks=ranks[-which(ranks$pcr_imp<0.01 & ranks$comp_imp<0.01),]
cor.test(ranks$pcr_imp,ranks$comp_imp,method="spearman")

## plot
p2=ggplot(ranks,aes(pcr_rank,comp_rank))+
  #geom_point()+
  geom_label(aes(label=var),size=2,fill=col,alpha=0.2)+
  scale_y_reverse(limits=c(max(c(ranks$comp_rank,ranks$pcr_rank))+3,0))+
  scale_x_reverse(limits=c(max(c(ranks$comp_rank,ranks$pcr_rank))+3,0))+
  #geom_abline(slope=1,linetype=2,size=0.5)+
  theme_bw()+
  labs(x="feature rank for PCR",
       y="feature rank for competence")+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))

## combine
library(patchwork)
setwd("~/Desktop/rodent comp")
png("model compare.png",width=7,height=3.5,units="in",res=300)
p1+p2
dev.off()

## same for relative importance
library(scales)
ggplot(ranks,aes(pcr_imp,comp_imp))+
  #geom_point()+
  #geom_text_repel(aes(label=var),force=1)+
  geom_label(aes(label=var),size=4,fill="grey90",alpha=0.7)+
  #scale_x_sqrt()+
  #scale_y_sqrt()+
  scale_x_continuous(trans=modulus_trans(p=0.01,offset=1))+
  scale_y_continuous(trans=modulus_trans(p=0.01,offset=1))+
  geom_abline(intercept=0,slope=1,linetype=2)+
  theme_bw()+
  labs(x="relative importance for PCR BRTs",
       y="relative importance for competence BRTs")+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))

## function for compiling across BRTs for a given predictor, all else equal
pdp_agg=function(mod,feature){
  
  ## get plot
  # pdep=pdp::partial(object=mod$mod,
  #                   pred.var=feature,
  #                   n.trees=mod$best,plot=F,
  #                   train=mod$traindata)
  
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
pdp_plot=function(bmods,feature,feature2=NULL){
  
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
                   size=1,colour=viridis(1),alpha=0.25)+
      
      ## add lines
      geom_line(size=1,alpha=1,colour="grey60")+
      
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
      labs(x=feature2,y="marginal effect")+
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
      geom_jitter(size=1,alpha=1,colour="grey60",width=0.1)+
      
      ## add mean
      geom_point(data=pmean,size=2,inherit.aes=F,shape=15,
                 aes(predictor,yhat))+
      
      ## add rug
      geom_rug(data=temp,inherit.aes=F,
               aes(predictor,yhat),
               sides="b",position="jitter",
               colour=viridis(1),alpha=0.25,
               na.rm=T)+
      
      ## theme
      theme_bw()+
      theme(axis.text=element_text(size=6),
            axis.title=element_text(size=7))+
      theme(axis.title.x=element_text(margin=margin(t=5,r=0,b=0,l=0)))+
      theme(axis.title.y=element_text(margin=margin(t=0,r=5,b=0,l=0)))+
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
      labs(x=feature2,y="marginal effect")+
      scale_y_continuous(limits=c(yrange[1]-0.01,yrange[2]+0.01),
                         labels=scales::number_format(accuracy=0.01))
    
  }
  
}

## sort ranks
ranks=ranks[order(ranks$comp_rank,decreasing=F),]

## compare function
pcomp=function(x){

  library(patchwork)
  a=pdp_plot(pcr_brts,x,x)
  a=a+labs(title="PCR")
  b=pdp_plot(comp_brts,x,x)
  b=b+labs(title="competence")
  a+b
}

## run through
pcomp(as.character(ranks$var)[1])
pcomp(as.character(ranks$var)[2])
pcomp(as.character(ranks$var)[3])
pcomp(as.character(ranks$var)[4])
pcomp(as.character(ranks$var)[5])
pcomp(as.character(ranks$var)[6])
pcomp(as.character(ranks$var)[7])
pcomp(as.character(ranks$var)[8])
pcomp(as.character(ranks$var)[9])
pcomp(as.character(ranks$var)[10])
pcomp(as.character(ranks$var)[11])
pcomp(as.character(ranks$var)[12])
pcomp(as.character(ranks$var)[13])

## average predictions: PCR
pcr_apreds=lapply(pcr_brts,function(x) x$predict)
pcr_apreds=do.call(rbind,pcr_apreds)

## aggregate
pcr_apreds=data.frame(aggregate(pred~treename,data=pcr_apreds,mean),
                  aggregate(cpred~treename,data=pcr_apreds,mean)['cpred'],
                  aggregate(hPCR~treename,data=pcr_apreds,prod)["hPCR"],
                  aggregate(competence~treename,data=pcr_apreds,prod)["competence"])

## type
pcr_apreds$type='PCR'

## average predictions: competence
comp_apreds=lapply(comp_brts,function(x) x$predict)
comp_apreds=do.call(rbind,comp_apreds)

## aggregate
comp_apreds=data.frame(aggregate(pred~treename,data=comp_apreds,mean),
                       aggregate(cpred~treename,data=comp_apreds,mean)['cpred'],
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

## compare
png("model pred.png",width=7,height=3.5,units="in",res=300)
p3=ggplot(apreds2,aes(PCR,competence))+
  geom_point(alpha=0.5,colour="black")+
  geom_smooth(method='gam',colour=col)+
  labs(x='predictions from PCR',
       y='predictions from competence')+
  theme_bw()+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))
p3+p3
dev.off()
cor(apreds2$PCR,apreds2$competence,method='spearman')
cor.test(apreds2$PCR,apreds2$competence,method='spearman')

## vis
ggplot(apreds,aes(pred,cpred,colour=factor(hPCR)))+geom_point()

## density
ggplot(apreds,aes(cpred))+
  geom_density(aes(colour=factor(hPCR),fill=factor(hPCR)),alpha=0.5)+
  facet_wrap(~type)+
  xlim(0,1)+
  labs(x=expression(paste("mean predicted probability of infection")))+
  theme_bw()+
  geom_vline(xintercept=0.5,linetype=2,size=0.5)+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=9),
        legend.text=element_text(size=7),
        legend.title=element_text(size=8))+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(legend.key.size = unit(0.5,"cm"),
        legend.box.margin=margin(t=0,r=0,b=-2.5,l=0),
        legend.margin=margin(0,0,0,0))+
  scale_fill_viridis_d(end=0.6)+
  scale_colour_viridis_d(end=0.6)+
  theme(legend.position="top")

## shortlist cpred
nrow(apreds[apreds$cpred>0.5 & apreds$hPCR==0,])

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

## logit
cdata$data$logit_PCR=car::logit(cdata$data$PCR)
cdata$data$logit_comp=car::logit(cdata$data$competence)

## lambda
pcr_lmod=pgls(logit_PCR~1,data=cdata,lambda="ML")
comp_lmod=pgls(logit_comp~1,data=cdata,lambda="ML")
summary(pcr_lmod)
summary(comp_lmod)
## moderate phylogenetic signal in predictions

## relate prediction
lmod=pgls(competence~PCR,data=cdata,lambda="ML")
summary(lmod)
lmod=pgls(logit_comp~logit_PCR,data=cdata,lambda="ML")
summary(lmod)
ggplot(cdata$data,aes(logit_PCR,logit_comp))+geom_point()+
  geom_smooth(method='lm')+
  labs(x='PCR-based predictions',
       y='competence-based predictions')

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

## predictions
set.seed(1)
pred_pf=gpf(Data=cdata$data,tree=cdata$phy,
           frmla.phylo=competence~phylo,
           family=gaussian,algorithm='phylo',nfactors=5,min.group.size=5)

## summarize
pred_pf_results=pfsum(pred_pf)$results

## save tree
cdata$data$hPCR2=factor(cdata$data$hPCR)
dtree=treeio::full_join(as.treedata(cdata$phy),cdata$data,by="label")

## make base
base=ggtree(dtree,size=0.5,layout='circular',branch.length='none')

## plot
gg=base+
  geom_tippoint(shape=21,size=1.5,
                aes(fill=hPCR2))+
  scale_fill_manual(values=c('white','black'))

## add clades
plus=1
pplus=plus
for(i in 1:nrow(pred_pf_results)){
  
  gg=gg+
    geom_hilight(node=pred_pf_results$node[i],
                 alpha=0.25,
                 fill=ifelse(pred_pf_results$clade>
                               pred_pf_results$other,pcols[2],pcols[1])[i])+
    geom_cladelabel(node=pred_pf_results$node[i],
                    label=pred_pf_results$factor[i],
                    offset=plus,
                    offset.text=pplus)
}
