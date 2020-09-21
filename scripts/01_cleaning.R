## rodent hantavirus cleaning
## danbeck@iu.edu

## clean environment & plots
rm(list=ls()) 
graphics.off()

## libraries
library(ape)
library(plyr)
library(caper)
library(phylofactor)
library(data.table)

## load in hantavirus data
setwd("~/Desktop/massinnombres/data")
hdata=read.csv("Muroid Taxonomy and Traits with Hantavirus Evidence.csv",header=T)

## synonym IUCN
hdata$IUCNname=hdata$Notes.regarding.IUCN.status

## load in traits
traits=read.csv("pan for PNAS Dryad.csv",header=T)
traits$X=NULL

## remove family
start='Abrocomidae'
end='Thryonomyidae'
trim=which(names(traits)==start):which(names(traits)==end)
traits=traits[,-trim]
rm(trim,start,end)

## remove duplicates
traits=traits[!duplicated(traits$Rodents),]

## load Upham phylogeny
setwd("~/Desktop/massinnombres/phylo")
tree=read.nexus('MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre')

## load in taxonomy
taxa=read.csv('taxonomy_mamPhy_5911species.csv',header=T)

## simplify to genus
gtaxa=taxa[!duplicated(taxa$gen),]
gtaxa=gtaxa[c('gen','fam','ord','clade','higher')]

## simplify datasets
hdata$PCR=hdata$X..papers.PCR.pos
hdata$PCRneg=hdata$X..papers.PCR.neg..whether.seropos.or.seroneg.
hdata$isolation=hdata$X..papers.isolated

## trim
hdata=hdata[c('Species_Name','PCR','PCRneg','isolation',
              'IUCNname')]

## geography
hdata$geography='New World'

## studies
hdata$studies=rowSums(hdata[c('PCR','PCRneg','isolation')],na.rm=T)

## hantavirus PCR
hdata$hPCR=ifelse(hdata$PCR>0,1,0)
hdata$hPCR=ifelse(is.na(hdata$hPCR),0,hdata$hPCR)

## virus isolation
hdata$competence=ifelse(hdata$isolation>0,1,0)
hdata$competence=ifelse(is.na(hdata$competence),0,hdata$competence)

## genus
hdata$gen=sapply(strsplit(as.character(hdata$Species_Name),'_'),function(x) x[1])

## combine
data=merge(hdata,gtaxa,by='gen')

## clean
rm(taxa,hdata,gtaxa)

## make simple names
data$tip=data$Species_Name
tree$tip=sapply(strsplit(tree$tip.label,'_'),function(x) paste(x[1],x[2],sep='_'))

## are all hdata in tree
data$intree=ifelse(data$tip%in%setdiff(data$tip,tree$tip),
                   'missing','upham')

## are all hdata in traits
data$intraits=ifelse(data$tip%in%setdiff(data$tip,
                                         traits$Rodents),
                     'missing','traits')

## just missing names
miss=data[c('tip','IUCNname','intree','intraits')]
miss=miss[miss$intree=='missing' | miss$intraits=='missing',]
miss=miss[order(miss$intree,miss$intraits),]

## export
setwd("~/Desktop/massinnombres/data")
write.csv(miss,'hantaro name mismatch.csv')
rm(miss)

## load in revised names
fix=read.csv('hantaro name mismatch_edit.csv',header=T)

## merge
fix=fix[c('tip','treename','traitname')]
data=merge(data,fix,by='tip',all.x=T)

## if blank, NA
data$traitname=ifelse(data$traitname=='',NA,as.character(data$traitname))
data$treename=ifelse(data$treename=='',NA,as.character(data$treename))

## if NA, tip
data$traitname=ifelse(is.na(data$traitname),as.character(data$tip),as.character(data$traitname))
data$treename=ifelse(is.na(data$treename),as.character(data$tip),as.character(data$treename))

## simplify
data=data[c('tip','studies','hPCR','competence','gen','fam','treename','traitname')]

## fix duplicates in phylogeny
set=data[c('studies','hPCR','competence','treename')]

## aggregate
set=aggregate(.~treename,set,sum)

## merge with meta
set2=merge(set,data[!duplicated(data$treename),c('tip','gen','fam','treename','traitname')],all.x=T,by='treename')

## simplify
data=set2
rm(set,set2,fix)

## fix binomial
data$hPCR=ifelse(data$hPCR>0,1,0)
data$competence=ifelse(data$competence>0,1,0)

## merge traits
traits$traitname=traits$Rodents
traits$trait=1
data=merge(data,traits,by='traitname',all.x=T)
rm(traits)

## setdiff
data$tree=ifelse(data$treename%in%setdiff(data$treename,tree$tip),'cut','keep')

## trim
bdata=data[-which(data$tree=='cut'),]

## fix tree
rtree=tree
rtree$tip.label=rtree$tip

## trim
rtree=keep.tip(rtree,bdata$treename)

## fix
rtree$tip=NULL

## fix
rtree=makeLabel(rtree)

## match
bdata=bdata[match(rtree$tip.label,bdata$treename),]

## save
bdata$label=bdata$treename
bdata$Species=bdata$treename

## merge
cdata=comparative.data(phy=rtree,data=bdata,names.col=treename,vcv=T,na.omit=F,warn.dropped=T)

## fix
cdata$data$tree=NULL

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
    
    dat[,paste0(resp,'_pf',i)]=ifelse(dat$genus%in%cladeget(pf,i),'factor','other')
    
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

## PCR
set.seed(1)
pcr_pf=gpf(Data=cdata$data,tree=cdata$phy,
        frmla.phylo=hPCR~phylo,
        family=binomial,algorithm='phylo',nfactors=5,min.group.size=5)

## summarize
pcr_pf_results=pfsum(pcr_pf)$results

## competence
set.seed(1)
hc_pf=gpf(Data=cdata$data,tree=cdata$phy,
           frmla.phylo=competence~phylo,
           family=binomial,algorithm='phylo',nfactors=5,min.group.size=1)

## summarize
hc_pf_results=pfsum(hc_pf)$results
