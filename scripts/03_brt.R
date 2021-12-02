## hantaro 03: rodent hantavirus BRT
## danbeck@ou.edu
<<<<<<< Updated upstream
## updated 12/01/2021
=======
## updated 11/30/2021
>>>>>>> Stashed changes

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
library(caret) 
library(InformationValue)
library(mgcv)

## load files
setwd("~/Desktop/hantaro/data/clean files")
data=read.csv('hantaro cleaned response and traits.csv')

## classify true negatives
data$type=ifelse(data$studies>0 & data$hPCR==0 & data$competence==0,"true negative","other")

## tabulate PCR and isolation
set=data
set$pcr=ifelse(set$hPCR==0,"PCR negative","PCR positive")
set$iso=ifelse(set$competence==0,"no isolation","isolation")
table(set$pcr,set$iso)

## which species is competent but no PCR record?
set[set$hPCR==0 & set$competence==1,"treename"]
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

## visualize distribution of NA
setwd("~/Desktop/hantaro/figs")
png("Figure S1.png",width=4,height=4,units="in",res=600)
ggplot(mval[!mval$column%in%c("gen","treename","studies","hPCR","competence","tip","fam"),],
       aes(comp))+
  geom_histogram(bins=50)+
  geom_vline(xintercept=0.25,linetype=2,size=0.5)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
  labs(y="frequency",
       x="trait coverage across muroid rodent species")+
  scale_x_continuous(labels = scales::percent)
dev.off()

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

## hyperparameter tuning ifelse
hok="ok"
if(hok!="ok"){
  
## hyperparameter grid
hgrid=expand.grid(n.trees=5000,
                  interaction.depth=c(2,3,4),
                  shrinkage=c(0.01,0.001,0.0005),
                  n.minobsinnode=4,
                  seed=seq(1,10,by=1))

## fix trees
hgrid$n.trees=ifelse(hgrid$shrinkage<0.001,hgrid$n.trees*3,hgrid$n.trees)

## trees, depth, shrink, min, prop
hgrid$id=with(hgrid,paste(n.trees,interaction.depth,shrinkage,n.minobsinnode))

## sort by id then seed
hgrid=hgrid[order(hgrid$id,hgrid$seed),]

## now add rows
hgrid$row=1:nrow(hgrid)

## factor id
hgrid$id2=factor(as.numeric(factor(hgrid$id)))

## function to assess each hyperpar combination
hfit=function(row,response){
  
  ## make new data
  ndata=set
  
  ## correct response
  ndata$response=ndata[response][,1]
  
  ## remove raw
  ndata$hPCR=NULL
  ndata$competence=NULL
  
  ## use rsample to split
  set.seed(hgrid$seed[row])
  split=initial_split(ndata,prop=0.7,strata="response")
  
  ## test and train
  dataTrain=training(split)
  dataTest=testing(split)
  
  ## yTest and yTrain
  yTrain=dataTrain$response
  yTest=dataTest$response
  
  ## BRT
  set.seed(1)
  gbmOut=gbm(response ~ . ,data=dataTrain,
             n.trees=hgrid$n.trees[row],
             distribution="bernoulli",
             shrinkage=hgrid$shrinkage[row],
             interaction.depth=hgrid$interaction.depth[row],
             n.minobsinnode=hgrid$n.minobsinnode[row],
             cv.folds=5,class.stratify.cv=TRUE,
             bag.fraction=0.5,train.fraction=1,
             n.cores=1,
             verbose=F)
  
  ## performance
  par(mfrow=c(1,1),mar=c(4,4,1,1))
  best.iter=gbm.perf(gbmOut,method="cv")
  
  ## predict with test data
  preds=predict(gbmOut,dataTest,n.trees=best.iter,type="response")
  
  ## known
  result=dataTest$response
  
  ## sensitiviy and specificity
  sen=InformationValue::sensitivity(result,preds)
  spec=InformationValue::specificity(result,preds)
  
  ## AUC on train
  auc_train=gbm.roc.area(yTrain,predict(gbmOut,dataTrain,n.trees=best.iter,type="response"))
  
  ## AUC on test
  auc_test=gbm.roc.area(yTest,predict(gbmOut,dataTest,n.trees=best.iter,type="response"))
  
  ## print
  print(paste("hpar row ",row," done; test AUC is ",auc_test,sep=""))
  
  ## save outputs
  return(list(best=best.iter,
              trainAUC=auc_train,
              testAUC=auc_test,
              spec=spec,
              sen=sen,
              wrow=row))
}

## run the function
hpars=lapply(1:nrow(hgrid),function(x) hfit(x,response="hPCR"))

## get results
hresults=data.frame(sapply(hpars,function(x) x$trainAUC),
                    sapply(hpars,function(x) x$testAUC),
                    sapply(hpars,function(x) x$spec),
                    sapply(hpars,function(x) x$sen),
                    sapply(hpars,function(x) x$wrow),
                    sapply(hpars,function(x) x$best))
names(hresults)=c("trainAUC","testAUC",
                  "spec","sen","row","best")

## combine and save
hsearch=merge(hresults,hgrid,by="row")

## save
hsearch$type="PCR"

## rerun for competence
hpars=lapply(1:nrow(hgrid),function(x) hfit(x,response="competence"))

## get results
hresults=data.frame(sapply(hpars,function(x) x$trainAUC),
                    sapply(hpars,function(x) x$testAUC),
                    sapply(hpars,function(x) x$spec),
                    sapply(hpars,function(x) x$sen),
                    sapply(hpars,function(x) x$wrow),
                    sapply(hpars,function(x) x$best))
names(hresults)=c("trainAUC","testAUC",
                  "spec","sen","row","best")

## combine and save
csearch=merge(hresults,hgrid,by="row")

## assign data type
csearch$type="competence"

## combine
search=rbind.data.frame(csearch,hsearch)
search$type=factor(search$type,levels=c("PCR","competence"))

## export
setwd("~/Desktop/hantaro/figs")
write.csv(search,"par tuning data summary.csv")

}else{
  
## load
setwd("~/Desktop/hantaro/figs")
search=read.csv("par tuning data summary.csv")
  
}

## factor parameters
search$shrinkage=factor(search$shrinkage)
lvl=rev(sort(unique(search$shrinkage)))
search$shrinkage=factor(search$shrinkage,levels=lvl); rm(lvl)

## factor other
search$interaction.depth=factor(search$interaction.depth)

## PCR beta regression for AUC
mod=gam(testAUC~interaction.depth*shrinkage,
        data=search[search$type=="PCR",],method="REML",family=betar)
anova(mod)

## virus isolation
mod=gam(testAUC~interaction.depth*shrinkage,
        data=search[search$type=="competence",],method="REML",family=betar)
anova(mod)

## PCR beta regression for sensitivity
mod=gam(sen~interaction.depth*shrinkage,
            data=search[search$type=="PCR",],method="REML",family=betar)
anova(mod)

## virus isolation
mod=gam(sen~interaction.depth*shrinkage,
        data=search[search$type=="competence",],method="REML",family=betar)
anova(mod)

## PCR beta regression for specificity
mod=gam(spec~interaction.depth*shrinkage,
        data=search[search$type=="PCR",],method="REML",family=betar)
anova(mod)

## virus isolation
mod=gam(spec~interaction.depth*shrinkage,
        data=search[search$type=="competence",],method="REML",family=betar)
anova(mod)

## fix type
search$type=plyr::revalue(search$type,
                    c("PCR"="RT-PCR",
                      "competence"="virus isolation"))

## recast from wide to long
search2=gather(search,measure,value,testAUC:sen)

## revalue and factor
search2$measure=plyr::revalue(search2$measure,
                              c("sen"="sensitivity",
                                "spec"="specificity",
                                "testAUC"="test AUC"))
search2$measure=factor(search2$measure,
                       levels=c("test AUC","sensitivity","specificity"))

## visualize
setwd("~/Desktop/hantaro/figs")
png("Figure S2.png",width=5,height=8,units="in",res=600)
set.seed(1)
ggplot(search2,aes(shrinkage,value,
                   colour=interaction.depth,fill=interaction.depth))+
  geom_boxplot(alpha=0.25)+
  geom_point(alpha=0.75,
             position = position_jitterdodge(dodge.width=0.75))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
  facet_grid(measure~type,scales="free_y",switch="y")+
  theme(strip.placement="outside",
        strip.background=element_blank())+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        strip.text=element_text(size=12))+
  theme(legend.position="top")+
  scale_color_brewer(palette="Pastel2")+
  scale_fill_brewer(palette="Pastel2")+
  guides(colour=guide_legend(title="interaction depth"),
         fill=guide_legend(title="interaction depth"))+
  labs(y=NULL,
       x="learning rate")+
  scale_y_continuous(n.breaks=4)
dev.off()

## clean
rm(search,search2,hok,mod)

## brt function to use different data partitions
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
  split=initial_split(ndata,prop=0.7,strata="response")
  
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
  
  ## known
  result=dataTest$response
  
  ## sensitiviy and specificity
  sen=InformationValue::sensitivity(result,preds)
  spec=InformationValue::specificity(result,preds)
  
  ## AUC on train
  auc_train=gbm.roc.area(yTrain,predict(gbmOut,dataTrain,n.trees=best.iter,type="response"))
  
  ## AUC on test
  auc_test=gbm.roc.area(yTest,predict(gbmOut,dataTest,n.trees=best.iter,type="response"))
  
  ## skip if poisson
  if(response=="cites"){
    
    perf=NA
    
  }else{
  
    ## inner loop if yTest is all 0
    if(var(yTest)==0){
      
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
              spec=spec,
              sen=sen,
              roc=perf,
              rinf=bars,
              predict=pred_data,
              traindata=dataTrain,
              testdata=dataTest,
              seed=seed))
}

## apply across 10 splits each
smax=101
pcr_brts=lapply(1:smax,function(x) brt_part(seed=x,response="hPCR"))
comp_brts=lapply(1:smax,function(x) brt_part(seed=x,response="competence"))

## write to files
setwd("~/Desktop/hantaro/data/clean files")
saveRDS(pcr_brts,"pcr brts.rds")
saveRDS(comp_brts,"comp brts.rds")

## run wos brts
pm_brts=lapply(1:(smax-1),function(x) brt_part(seed=x,response="cites"))

## write
saveRDS(pm_brts,"pm brts.rds")