## rodent hantavirus phylofactor
## danbeck@iu.edu

## clean environment & plots
rm(list=ls()) 
graphics.off()

## packages
library(ape)
library(phylofactor)
library(tidyverse)
library(data.table)
library(ggtree)
library(plyr)
library(tidyr)
library(caper)

## load hanta data
setwd("~/Desktop/massinnombres")
data=read.csv('MullEdgelist.csv',header=T,na.strings='')

## remove NA
data$Rodent.Host=as.character(data$Rodent.Host)
data$Rodent.Host=ifelse(data$Rodent.Host=='NA',NA,data$Rodent.Host)
data=data[!is.na(data$Rodent.Host),]

## fix detection columns
data$GenomeS=ifelse(is.na(data$GenomeS),0,1)
data$GenomeM=ifelse(is.na(data$GenomeM),0,1)
data$GenomeL=ifelse(is.na(data$GenomeL),0,1)

## fix isolation
data$Virus.Isolation=ifelse(is.na(data$Virus.Isolation),0,1)

## load trait data
traits=read.csv('pan for PNAS Dryad.csv',header=T)

## get family columns
fams=traits[which(names(traits)=='Abrocomidae'):which(names(traits)=='Thryonomyidae')]
fams=as.matrix(fams)
fams=apply(fams,2,function(x) ifelse(x>0,1,0))
traits$family=factor(fams%*%(1:ncol(fams)),labels = colnames(fams))

## remove binary
traits=traits[-c(which(names(traits)=='Abrocomidae'):which(names(traits)=='Thryonomyidae'))]

## check length of species names
table(sapply(strsplit(as.character(traits$Rodents),'_'),length))

## reduce data hosts
x=strsplit(as.character(data$Rodent.Host),' ')
data$genus=sapply(x,function(x) x[1])
data$species=sapply(x,function(x) x[2])

## make host name
data$host_name=paste(data$genus,data$species)

## make PCR
data$PCR=rowSums(data[c('GenomeS','GenomeM','GenomeL')])
data$PCR=ifelse(data$PCR>0,1,0)

## aggregate
rdata=data.frame(aggregate(PCR~host_name,data=data,sum,na.rm=T),
                 aggregate(Virus.Isolation~host_name,data=data,sum,na.rm=T)['Virus.Isolation'])

## binary
rdata$PCR=ifelse(rdata$PCR>0,1,0)
rdata$Virus.Isolation=ifelse(rdata$Virus.Isolation>0,1,0)

## match standard name
rdata$traitname=gsub(' ','_',rdata$host_name)

## fix
rdata$traitname=revalue(rdata$traitname,
                        c('Bolomys_lasiurus'='Necromys_lasiurus',
                          'Bolomys_obscurus'='Necromys_obscurus',
                          'Loxodontomys_micopus'='Loxodontomys_micropus'))

## setdiff
setdiff(rdata$traitname,traits$Rodents)

## note true 0/1 daa
rdata$true0=1

## merge
traits$traitname=traits$Rodents
data=merge(rdata,traits,by='traitname',all=T)

## clean
rm(x,traits,rdata,fams)

## replace NAs
data$PCR=replace_na(data$PCR,0)
data$Virus.Isolation=replace_na(data$Virus.Isolation,0)
data$true0=replace_na(data$true0,0)

## clean
data$host_name=NULL
data$X=NULL

## load supertree
setwd("~/Desktop/cleanbats_betacov/raw data")
tree=readRDS('Full Supertree.rds')

## get tips
n1=tree$tip.label
n2=data$traitname
keep=n1[n1%in%n2]

## trim tree
rtree=keep.tip(tree,keep)

## assess which names are missing
data$missing=ifelse(data$traitname%in%rtree$tip.label,'keep','cut')

## trim
bdata=data[-which(data$missing=='cut'),]
bdata=bdata[match(rtree$tip.label,bdata$traitname),]

## merge
cdata=comparative.data(phy=rtree,data=bdata,names.col=traitname,vcv=T,na.omit=F,warn.dropped=T)

## add correct names
cdata$data$hostnames=rownames(cdata$data)
cdata$data$Species=rownames(cdata$data)

## taxonomy
cdata$data$taxonomy=paste(cdata$data$family,cdata$data$MSW05_Genus,cdata$data$hostnames,sep='; ')

## Holm rejection procedure
HolmProcedure <- function(pf,FWER=0.05){
  if (pf$models[[1]]$family$family%in%c('gaussian',"Gamma","quasipoisson")){
    pvals <- sapply(pf$models,FUN=function(fit) summary(fit)$coefficients['phyloS','Pr(>|t|)'])
  } else {
    pvals <- sapply(pf$models,FUN=function(fit) summary(fit)$coefficients['phyloS','Pr(>|z|)'])
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

## set taxonomy
taxonomy=data.frame(cdata$data$taxonomy)
names(taxonomy)="taxonomy"
taxonomy$Species=rownames(cdata$data)
taxonomy=taxonomy[c("Species","taxonomy")]
taxonomy$taxonomy=as.character(taxonomy$taxonomy)

## gpf for Virus.Isolation
set.seed(1)
vpf=gpf(Data=cdata$data,tree=cdata$phy,
       frmla.phylo=Virus.Isolation~phylo,
       family=binomial,algorithm='phylo',nfactors=5)
keep=HolmProcedure(vpf)

## gpf for PCR
set.seed(1)
ppf=gpf(Data=cdata$data,tree=cdata$phy,
        frmla.phylo=PCR~phylo,
        family=binomial,algorithm='phylo',nfactors=5)
keep=HolmProcedure(ppf)
pf.tree(ppf,factors=1:keep,size=0.1,layout="circular")$ggplot

## set key
setkey(ppf$Data,'Species')

## summarize
pf.taxa(ppf,taxonomy,factor=1)$group1

## get clade
cdata$data$pf1=ifelse(cdata$data$hostnames%in%cladeget(ppf,1),'clade','other')

## table
table(cdata$data$pf1,cdata$data$PCR)
prop.table(table(cdata$data$pf1,cdata$data$PCR))*100

## refit with only this clade
set.seed(1)
ppf_final=gpf(Data=cdata$data,tree=cdata$phy,
        frmla.phylo=PCR~phylo,
        family=binomial,algorithm='phylo',nfactors=keep)

## predict
cdata$data$pred=predict(ppf_final,type='response')
hist(cdata$data$pred)