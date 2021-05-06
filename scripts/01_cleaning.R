## hantaro 01: rodent hantavirus cleaning
## danbeck@ou.edu

## clean environment & plots
rm(list=ls()) 
graphics.off()

## libraries
library(ape)
library(plyr)

## load in hantavirus data
setwd("~/Desktop/hantaro/data/raw")
hdata=read.csv("Muroid Taxonomy and Traits with Hantavirus Evidence.csv",header=T)

## fix missing PCR studies
hdata$X..papers.PCR.pos=ifelse(hdata$Species_Name=="Akodon_simulator",1,hdata$X..papers.PCR.pos)
hdata$X..papers.PCR.pos=ifelse(hdata$Species_Name=="Calomys_callosus",1,hdata$X..papers.PCR.pos)

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
setwd("~/Desktop/hantaro/phylo")
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
setwd("~/Desktop/hantaro/data/names")
write.csv(miss,'hantaro name mismatch.csv')
rm(miss)

## load in revised names with Than edits
fix=read.csv('hantaro name mismatch_edit_TM.csv',header=T)

## merge
fix=fix[c('tip','treename','traitname','proxy')]
data=merge(data,fix,by='tip',all.x=T)

## if blank, NA
data$traitname=ifelse(data$traitname=='',NA,as.character(data$traitname))
data$treename=ifelse(data$treename=='',NA,as.character(data$treename))

## if NA, tip
data$treename=ifelse(is.na(data$treename),as.character(data$tip),as.character(data$treename))

## trait name includes proxy
#data$traitname=ifelse(is.na(data$traitname),as.character(data$tip),as.character(data$traitname))
data$traitname=ifelse(data$intraits=='missing' & is.na(data$traitname),as.character(data$proxy),
                 ifelse(data$intraits=='missing' & !is.na(data$traitname),as.character(data$traitname),
                        as.character(data$tip)))

## simplify
data=data[c('tip','studies','hPCR','competence','gen','fam','treename','traitname')]
rm(fix)

## fix duplicates in phylogeny
set=data[c('studies','hPCR','competence','treename')]

## which treenames occur multiple times
unique(set$treename)[table(set$treename)>1]

## aggregate
set=aggregate(.~treename,set,sum)

## merge with meta
set2=merge(set,data[!duplicated(data$treename),c('tip','gen','fam','treename','traitname')],all.x=T,by='treename')

## simplify
data=set2
rm(set,set2)

## fix binomial
data$hPCR=ifelse(data$hPCR>0,1,0)
data$competence=ifelse(data$competence>0,1,0)

## merge traits
traits$traitname=traits$Rodents
traits$trait=1
data=merge(data,traits,by='traitname',all.x=T)
rm(traits)

## fix tree
rtree=tree
rtree$tip.label=rtree$tip

## trim
#rtree=keep.tip(rtree,data$treename)
rtree=keep.tip(rtree,rtree$tip.label[rtree$tip.label%in%data$treename])

## fix
rtree$tip=NULL

## fix
rtree=makeLabel(rtree)

## get ed
library(picante)
ed=evol.distinct(rtree,type='equal.splits')

## treename
ed$treename=ed$Species
ed$Species=NULL

## rename
ed$ed_equal=ed$w
ed$w=NULL

## merge into data
data=merge(data,ed,by='treename',all.x=T)
rm(ed)

## cleaning
data$trait=NULL

## pubmed citations
library(easyPubMed)

## function
counter=function(name){
  as.numeric(as.character(get_pubmed_ids(gsub('_','-',name))$Count))
}
citations=c()

## loop through
for(i in 1:length(data$treename)) {
  citations[i]=counter(data$treename[i])
  print(i)
}

## compile
cites=data.frame(treename=data$treename,
                 cites=citations)

## merge
data=merge(data,cites,by='treename')

## clean
rm(cites,citations,i,counter)

## export files
setwd("~/Desktop/hantaro/data/clean files")
write.csv(data,'hantaro cleaned response and traits.csv')
saveRDS(rtree,'rodent phylo trim.rds')