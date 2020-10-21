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
setwd("~/Github/hantaro/data/clean files")
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
#dums=dummy_cols(data["gen"])

## unique
#dums=dums[!duplicated(dums$gen),]

## merge
#data=merge(data,dums,by="gen",all.x=T)
#rm(dums)

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
#set$gen=NULL
set$fam=NULL

## remove studies
set$studies=NULL

################## Colin grabs the wheel and drives the car into a lake

library(dbarts)
library(embarcadero)
library(magrittr)
library(tidyverse)

set %<>% as_tibble()

#na.ratio <- unlist(map(set, ~sum(is.na(.))))/nrow(set) # A way of checking this using purrr

#set %<>% dplyr::select(!names(na.ratio[na.ratio>0.25])) # Drop the things again 

# That shouldn't work if the above stuff works correctly

#model0 <- bart(set[,3:ncol(set)], set$hPCR, keeptrees = TRUE)

#set <- na.omit(set)

fixer <- bart(x.train = set[,4:ncol(set)],
              y.train = set$hPCR,
              keeptrees = TRUE)

plot(fixer)

# X26.1_GR_Area_km2 X26.3_GR_MinLat_dd X26.4_GR_MidRangeLat_dd X26.6_GR_MinLong_dd X27.2_HuPopDen_Mean_n.km2 MammRich 

fixer.pred <- colMeans(dbarts:::predict.rbart(fixer, set, group.by = set$gen))
mean(fixer.pred)
mean(set$hPCR)

# This is so few variables! I think something ain't right here

model1 <- rbart_vi(hPCR ~ X26.1_GR_Area_km2 + X26.3_GR_MinLat_dd + X26.4_GR_MidRangeLat_dd + X26.6_GR_MinLong_dd + X27.2_HuPopDen_Mean_n.km2 + MammRich,
                   group.by = gen,
                   data = set,
                   offset = (mean(set$hPCR) - 0.50),
                   n.samples = 1000, n.burn = 100, n.chains = 1, 
                   n.threads = 1, n.trees = 200,
                   keepTrees = TRUE)
summary(model1)
model1.pred <- colMeans(dbarts:::predict.rbart(model1, set, group.by = set$gen))
mean(model1.pred)

