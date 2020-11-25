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