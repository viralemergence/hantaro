## hantaro 02: rodent hantavirus phylof
## danbeck@ou.edu

## clean environment & plots
rm(list=ls()) 
graphics.off()

## libraries
library(ape)
library(caper)
library(phylofactor)
library(data.table)
library(treeio)
library(ggtree)

## load files
setwd("~/Desktop/Github/hantaro/data/clean files")
data=read.csv('hantaro cleaned response and traits.csv')
rtree=readRDS('rodent phylo trim.rds')

## setdiff
data$tree=ifelse(data$treename%in%setdiff(data$treename,rtree$tip),'cut','keep')

## trim
bdata=data[-which(data$tree=='cut'),]

## match
bdata=bdata[match(rtree$tip.label,bdata$treename),]

## save
bdata$label=bdata$treename
bdata$Species=bdata$treename

## merge
cdata=comparative.data(phy=rtree,data=bdata,names.col=treename,vcv=T,na.omit=F,warn.dropped=T)

## fix
cdata$data$tree=NULL

## proportion of each
nrow(data)
round(prop.table(table(data$hPCR)),4)*100
round(prop.table(table(data$competence)),4)*100

## phylogenetic signal in response
## D of 0 = Brownian model, D of 1 = random (no phylogenetic signal)
set.seed(1)
mod1=phylo.d(cdata,binvar=hPCR,permut=1000); mod1
set.seed(1)
mod2=phylo.d(cdata,binvar=competence,permut=1000); mod2

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

## PCR
set.seed(1)
pcr_pf=gpf(Data=cdata$data,tree=cdata$phy,
           frmla.phylo=hPCR~phylo,
           family=binomial,algorithm='phylo',nfactors=3,min.group.size=5)

## summarize
pcr_pf_results=pfsum(pcr_pf)$results

## competence
set.seed(1)
hc_pf=gpf(Data=cdata$data,tree=cdata$phy,
          frmla.phylo=competence~phylo,
          family=binomial,algorithm='phylo',nfactors=2,min.group.size=5)

## summarize
hc_pf_results=pfsum(hc_pf)$results

## save tree
cdata$data$infect=factor(cdata$data$hPCR)
cdata$data$comp=factor(cdata$data$competence)
dtree=treeio::full_join(as.treedata(cdata$phy),cdata$data,by="label")

## fix palette
AlberPalettes <- c("YlGnBu","Reds","BuPu", "PiYG")
AlberColours <- sapply(AlberPalettes, function(a) RColorBrewer::brewer.pal(5, a)[4])
afun=function(x){
  a=AlberColours[1:x]
  return(a)
}

## make low and high
pcols=afun(2)

## set x max
plus=1
pplus=plus+1

## fix taxa
pcr_pf_results$taxa[1]="subclade~of~italic(Peromyscus)"
pcr_pf_results$taxa[2]="italic(Oligoryzomys)"

## ggtree
gg=ggtree(dtree,size=0.25)+
  geom_tippoint(aes(colour=infect),shape=15)+
  scale_colour_manual(values=c("grey80","black"))+
  guides(colour=F)

## add clades
for(i in 1:nrow(pcr_pf_results)){
  
  gg=gg+
    geom_hilight(node=pcr_pf_results$node[i],
                 alpha=0.25,
                 fill=ifelse(pcr_pf_results$clade>
                               pcr_pf_results$other,pcols[2],pcols[1])[i])+
    geom_cladelabel(node=pcr_pf_results$node[i],
                    label=pcr_pf_results$taxa[i],
                    offset=pplus,
                    hjust=0.75,
                    offset.text=pplus*2,
                    parse=T,
                    angle=90)
}
pcr_gg=gg

## plot competence
comp_gg=ggtree(dtree,size=0.25)+
  geom_tippoint(aes(colour=comp),shape=15)+
  scale_colour_manual(values=c("grey80","black"))+
  guides(colour=F)

## print
library(ggpubr)
setwd("~/Desktop/Github/hantaro/figs")
png("Figure 1.png",width=6,height=6,units="in",res=300)
ggarrange(pcr_gg,comp_gg,ncol=2,widths=c(1.2,1),
          labels=c("(a) RT-PCR","(b) virus isolation"),
          label.x=c(-0.1,-0.2),
          font.label=list(face="plain",size=12))
dev.off()

# ## add in hantavirus studies for PCR
# set.seed(1)
# pcr_pf_study=gpf(Data=cdata$data,tree=cdata$phy,
#            frmla.phylo=hPCR~phylo,
#            weights=cdata$data$studies,
#            family=binomial,algorithm='phylo',nfactors=5,min.group.size=5)
# HolmProcedure(pcr_pf_study)
# 
# ## for competence
# set.seed(1)
# hc_pf_study=gpf(Data=cdata$data,tree=cdata$phy,
#                  frmla.phylo=competence~phylo,
#                  weights=cdata$data$studies,
#                  family=binomial,algorithm='phylo',nfactors=5,min.group.size=5)
# HolmProcedure(hc_pf_study)
# 
# ## log1p pubmed cites
# cdata$data$logcites=log1p(cdata$data$cites)
# 
# ## add in pubmed for PCR
# set.seed(1)
# pcr_pf_pm=gpf(Data=cdata$data,tree=cdata$phy,
#                  frmla.phylo=hPCR~phylo,
#                  weights=cdata$data$logcites,
#                  family=binomial,algorithm='phylo',nfactors=5,min.group.size=5)
# HolmProcedure(pcr_pf_pm)
# 
# ## for competence
# set.seed(1)
# hc_pf_pm=gpf(Data=cdata$data,tree=cdata$phy,
#                 frmla.phylo=competence~phylo,
#                 weights=cdata$data$logcites,
#                 family=binomial,algorithm='phylo',nfactors=3,min.group.size=5)
# HolmProcedure(hc_pf_pm)
# 
# ## model studies and citations themselves
# set.seed(1)
# study_pf=gpf(Data=cdata$data,tree=cdata$phy,
#              frmla.phylo=studies~phylo,
#              family=poisson,algorithm='phylo',nfactors=5,min.group.size=5)
# HolmProcedure(study_pf)
# pfsum(study_pf)$results
# 
# ## citations
# set.seed(1)
# pm_pf=gpf(Data=cdata$data,tree=cdata$phy,
#              frmla.phylo=cites~phylo,
#              family=poisson,algorithm='phylo',nfactors=5,min.group.size=5)
# HolmProcedure(pm_pf)
# pfsum(pm_pf)$results
