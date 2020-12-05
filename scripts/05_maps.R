
library(classInt)
library(tidyverse)
library(raster)
library(rgdal)
library(dismo)
library(XML)
library(maps)
library(sp)
library(velox)


setwd("~/Github/hantaro")

set.seed(12345)

# Read in the data

pred <- read_csv("./data/clean files/hantaro predictions.csv")

# 1. Threshold the results

library(PresenceAbsence)

###################### INTERMISSION: THRESHOLD IMPACTS

ts.p <- optimal.thresholds(data.frame(pred[,c('treename','PCR','pred_pcr')]),
                            threshold = 10001,
                            opt.methods = c(2,4,5,10),
                            req.sens = 0.90,
                            na.rm = TRUE)

cut.p <- function(x) {sum(pred$pred_pcr[pred$PCR==0] > x)}

sapply(unlist(ts.p[2]), cut.p)


ts.c <- optimal.thresholds(data.frame(pred[,c('treename','competence','pred_comp')]),
                           threshold = 10001,
                           opt.methods = c(2,4,5,10),
                           req.sens = 0.90,
                           na.rm = TRUE)

cut.c <- function(x) {sum(pred$pred_comp[pred$competence==0] > x)}

sapply(unlist(ts.c[2]), cut.c)




###################### MOVE FORWARD

t.pcr <- optimal.thresholds(data.frame(pred[,c('treename','PCR','pred_pcr')]),
                             threshold = 10001,
                             opt.methods = 10,
                             req.sens = 0.95,
                             na.rm = TRUE)

t.comp <- optimal.thresholds(data.frame(pred[,c('treename','competence','pred_comp')]),
                            threshold = 10001,
                            opt.methods = 10,
                            req.sens = 0.95,
                            na.rm = TRUE)

# Threshold the results to binary outputs

pred %>%
  mutate(bin_comp = (pred_comp > t.comp$pred_comp),
         bin_pcr = (pred_pcr) > t.pcr$pred_pcr) -> pred

# How many predicted undiscovered hosts by PCR?

table(pred$pred_pcr[pred$PCR==0] > t.pcr$pred_pcr)

# How many predicted undiscovered hosts by competence

table(pred$pred_comp[pred$competence==0] > t.comp$pred_comp)

# Do the predicted competence hosts overlap with PCR

pred %>% filter(competence==0, pred_comp > t.comp$pred_comp)

# Pull out the relevant lists 

pred %>% filter(competence==1) %>% pull(treename) %>% gsub("_"," ",.) -> known.comp

pred %>% filter(PCR==1) %>% pull(treename) %>% gsub("_"," ",.) -> known.pcr

pred %>% filter(bin_comp==1) %>% pull(treename) %>% gsub("_"," ",.) -> pred.comp

pred %>% filter(bin_pcr==1) %>% pull(treename) %>% gsub("_"," ",.) -> pred.pcr

sort(pred.pcr[!(pred.pcr %in% known.pcr)])

sort(pred.comp[!(pred.comp %in% known.comp)])

# 2. Let's make some maps?

library(fasterize)
library(rgdal)
library(raster)
library(sf)

iucn <- st_read(dsn = "C:/Users/cjcar/Dropbox/HowManyHelminths2019", layer='TERRESTRIAL_MAMMALS')

r <- disaggregate(getData("worldclim",var="alt",res=2.5)*0,2) # Make a blank raster

# Create four layers

iucn.1 <- iucn[iucn$binomial %in% known.comp,] 
iucn.2 <- iucn[iucn$binomial %in% known.pcr,] 
iucn.3 <- iucn[iucn$binomial %in% pred.comp,] 
iucn.4 <- iucn[iucn$binomial %in% pred.pcr,] 

map.knc <- (fasterize(iucn.1, r, fun="sum"))
map.knp <- (fasterize(iucn.2, r, fun="sum"))
map.prc <- (fasterize(iucn.3, r, fun="sum"))
map.prp <- (fasterize(iucn.4, r, fun="sum"))

fix <- function(x) {sum(x,r,na.rm=TRUE)+r} # This adds zeros for the continental area

map.knc <- fix(map.knc)
map.knp <- fix(map.knp)
map.prc <- fix(map.prc)
map.prp <- fix(map.prp)

raster::stack(map.knp, map.knc, map.prp, map.prc) %>% 
  crop(c(-170,-25,-90,90)) %>% 
  raster::trim() -> maps

names(maps) <- c('KnownPCR', 'KnownComp', 'PredPCR', 'PredComp')

# Generate the actual visualization

library(rasterVis)
library(RColorBrewer)

mycolors <- colorRampPalette(rev(brewer.pal(10,"Spectral")))(21)
mycolors[1] <- "#C0C0C0"

rasterVis::levelplot(maps,  
                     col.regions = mycolors,
                     #at = seq(0, 15, 1),
                     alpha = 0.5, 
                     scales=list(alternating=FALSE),
                     par.strip.text=list(cex=0),
                     xlab = NULL, ylab = NULL,
                     maxpixels = 5e6)

##########################################################################

# Time for bivariate 

foot <- raster('~/GitHub/hantaro/data/footprint/wildareas-v3-2009-human-footprint.tif')

memory.limit(20000000000)

foot <- projectRaster(foot, maps[[1]])

foot[foot>50] <- 0

foot <- foot + maps[[1]]*0

# Generate an "all reservoirs" layer

iucn.all <- iucn[iucn$binomial %in% c(known.pcr,
                                    pred.pcr),]

map.all <- (fasterize(iucn.all, r, fun="sum"))

fix <- function(x) {sum(x,r,na.rm=TRUE)+r}

map.all <- fix(map.all)

map.all <- map.all + foot*0

# Generate the bivariate square

colmat<-function(nquantiles=10, upperleft=rgb(0,150,235, maxColorValue=255), upperright=rgb(130,0,80, maxColorValue=255), bottomleft="grey", bottomright=rgb(255,230,15, maxColorValue=255), xlab="x label", ylab="y label"){
  my.data<-seq(0,1,.01)
  my.class<-classIntervals(my.data,n=nquantiles,style="quantile")
  my.pal.1<-findColours(my.class,c(upperleft,bottomleft))
  my.pal.2<-findColours(my.class,c(upperright, bottomright))
  col.matrix<-matrix(nrow = 101, ncol = 101, NA)
  for(i in 1:101){
    my.col<-c(paste(my.pal.1[i]),paste(my.pal.2[i]))
    col.matrix[102-i,]<-findColours(my.class,my.col)}
  plot(c(1,1),pch=19,col=my.pal.1, cex=0.5,xlim=c(0,1),ylim=c(0,1),frame.plot=F, xlab=xlab, ylab=ylab,cex.lab=1.3)
  for(i in 1:101){
    col.temp<-col.matrix[i-1,]
    points(my.data,rep((i-1)/100,101),pch=15,col=col.temp, cex=1)}
  seqs<-seq(0,100,(100/nquantiles))
  seqs[1]<-1
  col.matrix<-col.matrix[c(seqs), c(seqs)]}

col.matrix<-colmat(nquantiles=5, upperleft="#be64ac", upperright="#3b4994",
                   bottomleft="#e8e8e8", bottomright="#5ac8c8", 
                   ylab = "Human footprint", xlab = "Hantavirus hosts")

# Generate the bivariate map


bivariate.map<-function(rasterx, rastery, colormatrix=col.matrix, nquantiles=10){
  quanmean<-getValues(rasterx)
  temp<-data.frame(quanmean, quantile=rep(NA, length(quanmean)))
  brks<-with(temp, quantile(unique(temp),na.rm=TRUE, probs = c(seq(0,1,1/nquantiles))))
  r1<-within(temp, quantile <- cut(quanmean, breaks = brks, labels = 2:length(brks),include.lowest = TRUE))
  quantr<-data.frame(r1[,2]) 
  quanvar<-getValues(rastery)
  temp<-data.frame(quanvar, quantile=rep(NA, length(quanvar)))
  brks<-with(temp, quantile(unique(temp),na.rm=TRUE, probs = c(seq(0,1,1/nquantiles))))
  r2<-within(temp, quantile <- cut(quanvar, breaks = brks, labels = 2:length(brks),include.lowest = TRUE))
  quantr2<-data.frame(r2[,2])
  as.numeric.factor<-function(x) {as.numeric(levels(x))[x]}
  col.matrix2<-colormatrix
  cn<-unique(colormatrix)
  for(i in 1:length(col.matrix2)){
    ifelse(is.na(col.matrix2[i]),col.matrix2[i]<-1,col.matrix2[i]<-which(col.matrix2[i]==cn)[1])}
  cols<-numeric(length(quantr[,1]))
  for(i in 1:length(quantr[,1])){
    a<-as.numeric.factor(quantr[i,1])
    b<-as.numeric.factor(quantr2[i,1])
    cols[i]<-as.numeric(col.matrix2[b,a])}
  r<-rasterx
  r[1:length(r)]<-cols
  return(r)}

layers <- raster::stack(foot, maps$PredPCR)
layers <- embarcadero::bigstack(layers, 5)
bivmap<-bivariate.map(layers[[2]], layers[[1]], colormatrix = col.matrix, nquantiles = 5)
plot(bivmap, frame.plot = TRUE, axes = F, box = T, add = F, legend = F, col = as.vector(col.matrix), asp = 1)
map(interior = T, add = T)

layers <- raster::stack(foot, maps$PredComp)
layers <- embarcadero::bigstack(layers, 5)
bivmap<-bivariate.map(layers[[2]], layers[[1]], colormatrix = col.matrix, nquantiles = 5)
plot(bivmap, frame.plot = TRUE, axes = F, box = T, add = F, legend = F, col = as.vector(col.matrix), asp = 1)
map(interior = T, add = T)
