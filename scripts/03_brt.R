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

## if mostly homogenous (97%)
vars$keep=ifelse(vars$var<0.97,"keep","cut")
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

## if missing 25% or more
mval$comp=round(mval$comp,2)
mval$keep=ifelse(mval$comp>=0.25,"keep","cut")
mval=mval[order(mval$keep),]
keeps=mval[-which(mval$keep=="cut"),]$column

## order
mval=mval[order(mval$comp),]

## drop if not well represented
data=data[keeps]
rm(mval,keeps)

## make binary columns for genus
dums=dummy_cols(data["gen"])

## unique
dums=dums[!duplicated(dums$gen),]

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

## make simplified
set=data
set$treename=NULL
set$tip=NULL
set$X1.1_ActivityCycle=factor(set$X1.1_ActivityCycle)
set$X12.2_Terrestriality=factor(set$X12.2_Terrestriality)
set$X6.2_TrophicLevel=factor(set$X6.2_TrophicLevel)
set$gen=NULL
set$fam=NULL

## remove studies
set$studies=NULL

## function to use different data partitions
brt_part=function(seed,response){
  
  ## make new data
  ndata=set
  
  ## correct response
  ndata$response=ndata[response][,1]
  
  ## remove raw
  ndata$hPCR=NULL
  ndata$competence=NULL
  
  ## use rsample to split
  set.seed(seed)
  split=initial_split(ndata,prop=0.8,strata="response")
  
  ## test and train
  dataTrain=training(split)
  dataTest=testing(split)
  
  ## yTest and yTrain
  yTrain=dataTrain$response
  yTest=dataTest$response
  
  ## BRT
  set.seed(1)
  gbmOut=gbm(response ~ . ,data=dataTrain,
             n.trees=5000,
             distribution="bernoulli",
             shrinkage=0.001,
             interaction.depth=4,
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
  
  ## ROC
  pr=prediction(preds,dataTest$response)
  perf=performance(pr,measure="tpr",x.measure="fpr")
  perf=data.frame(perf@x.values,perf@y.values)
  names(perf)=c("fpr","tpr")
  
  ## add seed
  perf$seed=seed
  
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

## apply across 10 splits
smax=10
pcr_brts=lapply(1:smax,function(x) brt_part(seed=x,response="hPCR"))
comp_brts=lapply(1:smax,function(x) brt_part(seed=x,response="competence"))

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

## visualize
ggplot(adata,aes(response,auc,colour=response))+
  geom_boxplot()+
  geom_jitter(width=0.1)+
  guides(colour=F)+
  theme_bw()

## summarize AUC on test
hist(sapply(brts,function(x) x$testAUC))
mean(sapply(brts,function(x) x$testAUC))
se(sapply(brts,function(x) x$testAUC))

## check iterations of the tree
max(sapply(brts,function(x) x$best))
hist(sapply(brts,function(x) x$best),xlim=c(0,brts[[1]]$mod$n.trees))
abline(v=brts[[1]]$mod$n.trees,col="red")
table(sapply(brts,function(x) x$best)<brts[[1]]$mod$n.trees)

## aggregate ROCs
rocs=lapply(brts,function(x) x$roc)
rocs=do.call(rbind,rocs)

## plot
ggplot(rocs,aes(fpr,tpr,group=seed))+
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
vinf=lapply(brts,function(x) x$rinf)
vinf=do.call(rbind,vinf)

## as percent
vinf$rel.inf=vinf$rel.inf/100

## aggregate mean
vdata=data.frame(aggregate(rel.inf~var,data=vinf,mean),
                 aggregate(rel.inf~var,data=vinf,se)["rel.inf"])
names(vdata)=c("var","rel.inf","rse")
vdata=vdata[order(vdata$rel.inf,decreasing=T),]

## make rmin and rmax
vdata$rmin=vdata$rel.inf-vdata$rse
vdata$rmax=vdata$rel.inf+vdata$rse

## minimize
vset=vdata[-which(round(vdata$rel.inf,3)<=0.01),]

## vtops
vtops=vset[order(vset$rel.inf,decreasing=T),]

## plot
ggplot(vset,aes(reorder(var,rel.inf),rel.inf))+
  geom_errorbar(aes(ymin=rmin,ymax=rmax),width=0)+
  #scale_y_reverse()+
  coord_flip()+
  geom_point()+
  scale_y_continuous(labels = scales::percent_format(accuracy=1))+
  guides(colour=guide_legend(title="",override.aes=list(linetype=0)))+
  
  ## theme
  theme_bw()+
  theme(legend.position="top")+
  theme(axis.text=element_text(size=9),
        #axis.text.x=element_text(angle=67.5,hjust=1),
        axis.title=element_text(size=10),
        legend.text=element_text(size=8),
        legend.key.size = unit(0.5,"cm"),
        legend.box.margin=margin(t=-2.5,r=0,b=-2.5,l=0),
        legend.margin=margin(0,0,0,0))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
  labs(x="",y="relative importance")

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

## plot the top features
pdp_plot(brts,as.character(vtops$var[1]),as.character(vtops$var[1]))
pdp_plot(brts,as.character(vtops$var[2]),as.character(vtops$var[2]))
pdp_plot(brts,as.character(vtops$var[3]),as.character(vtops$var[3]))
pdp_plot(brts,as.character(vtops$var[4]),as.character(vtops$var[4]))
pdp_plot(brts,as.character(vtops$var[5]),as.character(vtops$var[5]))
pdp_plot(brts,as.character(vtops$var[6]),as.character(vtops$var[6]))
pdp_plot(brts,as.character(vtops$var[7]),as.character(vtops$var[7]))
pdp_plot(brts,as.character(vtops$var[8]),as.character(vtops$var[8]))
pdp_plot(brts,as.character(vtops$var[9]),as.character(vtops$var[9]))
pdp_plot(brts,as.character(vtops$var[10]),as.character(vtops$var[10]))
pdp_plot(brts,as.character(vtops$var[11]),as.character(vtops$var[11]))
pdp_plot(brts,as.character(vtops$var[12]),as.character(vtops$var[12]))
pdp_plot(brts,as.character(vtops$var[13]),as.character(vtops$var[13]))

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

## compare
ggplot(apreds2,aes(PCR,competence))+geom_point()+
  geom_smooth(method='gam')+
  labs(x='PCR-based predictions',
       y='competence-based predictions')
cor(apreds2$PCR,apreds2$competence,method='spearman')

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
