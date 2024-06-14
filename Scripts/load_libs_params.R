#PURPOSE:
# To load packages and set common parameters for data processing and mapping for BBRKC bycatch SDMS

#AUTHOR:
# Emily Ryznar - NOAA/NMFS RACE Shellfish Assessment Program (emily.ryznar@noaa.gov)

### LOAD PACKAGES  --------------------------------------------------------------
library(sf)
library(ggmap)
library(rgdal)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(ggpubr)
library(raster)
library(terra)
library(lubridate)
library(akgfmaps)
library(corrplot)
library(usdm)
library(dismo)
library(biomod2)
library(gam)
library(gbm) 
library(pROC)
library(ggrepel)
library(geosphere)
library(usdm)
library(viridis)

### SET SPATIAL DETAILS ---------------------------------------------------------
  crs.latlon <- "epsg:4326" #lat lon crs
  
  map.crs <- "EPSG:3338"
  
  in.crs = "+proj=longlat +datum=NAD83"

# LOAD SPATIAL LAYERS -----------------------------------------------------------
  region_layers <- akgfmaps::get_base_layers(select.region = "bs.south", set.crs="auto")
  
  survey_gdb <- "./Data/SAP_layers.gdb"
  
  readOGR(dsn=survey_gdb,layer="BristolBaySurveyStrata") %>%
    vect(crs = crs.latlon) %>%
    project(map.crs) -> BB_strata
  
  st_read("./Data/Closure areas/BLZ.shp") %>%
    vect() -> BLZ1 
  st_read("./Data/Closure areas/sepoct_closure.shp") %>%
    vect() -> janfebsepoct_closure
  st_read("./Data/Closure areas/sepoct_closure_wsub.shp") %>%
    vect() -> janfebsepoct_closure_wsub
  st_read("./Data/Closure areas/aprmay_closure.shp") %>%
    vect() -> aprmay_closure
  st_read("./Data/Closure areas/aprmay_closure_wsub.shp") %>%
    vect() -> aprmay_closure_wsub
  st_read("./Data/Closure areas/RKCSA.shp") %>%
    vect() -> RKCSA
  st_read("./Data/Closure areas/RKCSA_sub.shp") %>%
    vect() -> RKCSA_sub
  st_read("./Data/Closure areas/ns_trawl.shp") %>%
    vect() -> ns_trawl
  st_read("./Data/Closure areas/fivesxtn.shp") %>%
    vect() -> fivesxtn
  st_read("./Data/Closure areas/NBBTCA.shp") %>%
    vect() -> NBBTCA
  st_read("./Data/Closure areas/RKCSA_subarea.shp") %>%
    vect() -> RKCSA_subarea

# OTHER PARAMS ------------------------------------------------------------------
# Best training/testing iterations from model performance based on AUC, RMSE, and RHO
  lm_iter <- 9 #legal male
  im_iter <- 10 #immature male
  mf_iter <- 8 #mature female
  imf_iter <- 1 # immature female
  
# Thresholds
  thres <- read.csv("./Output/thresh.selection.csv")
  
  lm.thres <- thres %>%
    filter(category == "Legal male") %>%
    pull(max.kappa.thresh)
  
  im.thres <- thres %>%
    filter(category == "Immature male") %>%
    pull(max.kappa.thresh)
  
  mf.thres <- thres %>%
    filter(category == "Mature female") %>%
    pull(max.kappa.thresh)
  
  imf.thres <- thres %>%
    filter(category == "Immature female") %>%
    pull(max.kappa.thresh)
  