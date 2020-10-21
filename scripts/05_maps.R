
library(tidyverse)

setwd("~/Github/hantaro")

# Read in the data

pred <- read_csv("./data/clean files/hantaro predictions.csv")

# 1. Threshold the results

library(PresenceAbsence)

t.pcr <- optimal.thresholds(data.frame(pred[,c('treename','PCR','pred_pcr')]),
                             threshold = 10001,
                             opt.methods = 10,
                             req.sens = 0.9,
                             na.rm = TRUE)

t.comp <- optimal.thresholds(data.frame(pred[,c('treename','competence','pred_comp')]),
                            threshold = 10001,
                            opt.methods = 10,
                            req.sens = 0.9,
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
                     at = seq(0, 15, 1),
                     alpha = 0.5, 
                     scales=list(alternating=FALSE),
                     par.strip.text=list(cex=0),
                     xlab = NULL, ylab = NULL,
                     maxpixels = 5e6)