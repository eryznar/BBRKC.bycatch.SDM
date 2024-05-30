#PURPOSE:
# To compile data to date for BBRKC bycatch SDMs and process as needed

#AUTHOR:
# Emily Ryznar - NOAA/NMFS RACE Shellfish Assessment Program (emily.ryznar@noaa.gov)

### LOAD PROCESSING PARAMETERS --------------------------------------------------
source("./Scripts/load_libs_params.R")

years <- c(1995:2019, 2021:2023)

# PROCESS AND STACK SPATIAL COVARIATES ------------------------------------------

  # Read in raster layers
  ice<- rast("./Data/ice.tif") #Jan-Feb and Mar-Apr
  sst<- rast("./Data/sst.tif") #Jan-Dec in two month increments
  bt<- rast("./Data/bt95.23.tif") #May-Jun, Jul-Aug 
  lm_sap<- rast("./Data/lm_sap_CPUE.EBS.tif") %>%
    mask(BB_strata)
  im_sap<- rast("./Data/im_sap_CPUE.EBS.tif") %>%
    mask(BB_strata)
  mf_sap<- rast("./Data/mf_sap_CPUE.EBS.tif")%>%
    mask(BB_strata)
  imf_sap<- rast("./Data/imf_sap_CPUE.EBS.tif")%>%
    mask(BB_strata)
  sed <- rast("./Data/EBS_phi_1km.grd") #already has crs defined, same across all years
  depth <- rast("./Data/z_Layer1.tif") 
  yfs <- rast("./Data/yfs_cpue95.23.tif") #yellowfin sole cpue from summer survey
  rs <- rast("./Data/rs_cpue95.23.tif") #northern rock sole cpue from summer survey

  #Set input crs for rasters that don't already have them
  crs(sst) <- in.crs
  crs(ice) <- in.crs

  #Project rasters to the same crs
  sst2 <- terra::project(sst, map.crs)
  ice2 <- terra::project(ice, map.crs)
  sed2 <- terra::project(sed, map.crs)
  lm_sap2 <- terra::project(lm_sap, map.crs)
  im_sap2 <- terra::project(im_sap, map.crs)
  mf_sap2 <- terra::project(mf_sap, map.crs)
  imf_sap2 <- terra::project(imf_sap, map.crs)
  bt2 <- terra::project(bt, map.crs)
  depth2 <- terra::project(depth, map.crs)
  yfs2 <- terra::project(yfs, map.crs)
  rs2 <- terra::project(rs, map.crs)

  #Set extents and crop rasters that don't match the coarsest raster
  rast_ext <- c(-1500000, -170000, 539823, 1600000)
  
  sst3 <- crop(sst2, rast_ext)
  ice3 <- crop(ice2, rast_ext)
  bt3 <- crop(bt2, rast_ext)
  sed3 <- crop(sed2, rast_ext)
  lm_sap3 <- crop(lm_sap2, rast_ext)
  im_sap3 <- crop(im_sap2, rast_ext)
  mf_sap3 <- crop(mf_sap2, rast_ext)
  imf_sap3 <- crop(imf_sap2, rast_ext)
  depth3 <- crop(depth2, rast_ext)
  yfs3 <- crop(yfs2, rast_ext)
  rs3 <- crop(rs2, rast_ext)


  #Resample to match bt raster
  resample(sst3, bt3) -> sst4
  resample(ice3, bt3) -> ice4
  resample(depth3, bt3) -> depth4
  resample(lm_sap3, bt3) -> lm_sap4
  resample(im_sap3, bt3) -> im_sap4
  resample(mf_sap3, bt3) -> mf_sap4
  resample(imf_sap3, bt3) -> imf_sap4
  resample(sed3, bt3) -> sed4
  resample(yfs3, bt3) -> yfs4
  resample(rs3, bt3) -> rs4

  #Create raster stack of Environmental variables with the same crs, resolution, and extent, crop to Bristol Bay
  BB_extent <- c(-900000, -190000, 539823, 1050000)
  
  c(sst4, ice4, bt3, rep(sed4, length(years)), lm_sap4, rep(depth4,length(years)), yfs4, rs4) %>%
    terra::mask(BB_strata) %>%
    crop(BB_extent) -> lm_preds
  names(lm_preds[[253:280]]) <- paste("Sed", years)
  names(lm_preds[[308:335]]) <- paste("Depth", years)# set names with years for sediment and depth layers
  
  c(sst4, ice4, bt3, rep(sed4, length(years)), im_sap4, rep(depth4,length(years)), yfs4, rs4) %>%
    terra::mask(BB_strata) %>%
    crop(BB_extent) -> im_preds
  names(im_preds[[253:280]]) <- paste("Sed", years) # set names with years for sediment layers
  names(im_preds[[308:335]]) <- paste("Depth", years)# set names with years for sediment and depth layers
  
  c(sst4, ice4, bt3, rep(sed4, length(years)), mf_sap4, rep(depth4,length(years)), yfs4, rs4) %>%
    terra::mask(BB_strata) %>%
    crop(BB_extent) -> mf_preds
  names(mf_preds[[253:280]]) <- paste("Sed", years) # set names with years for sediment layers
  names(mf_preds[[308:335]]) <- paste("Depth", years)# set names with years for sediment and depth layers
  
  c(sst4, ice4, bt3, rep(sed4, length(years)), imf_sap4, rep(depth4,length(years)), yfs4, rs4) %>%
    terra::mask(BB_strata) %>%
    crop(BB_extent) -> imf_preds
  names(imf_preds[[253:280]]) <- paste("Sed", years) # set names with years for sediment layers
  names(imf_preds[[308:335]]) <- paste("Depth", years)# set names with years for sediment and depth layers
  
 # Load directed fishery target catch data
  # Legal males
  catch_lm <- read.csv("./Data/rkc_gfbycatch.csv")
  
  TAC_wide <- catch_lm %>% # Split TAC by rock sole and yellowfin sole
    dplyr::select(year, Target, TAC) %>%
    unique() %>%
    pivot_wider(., id_cols = year, names_from = Target, values_from = TAC) %>%
    rename("RS_TAC" = `Rock Sole - BSAI`, "YFS_TAC" = `Yellowfin Sole - BSAI`)
  
  catch_lm <- right_join(catch_lm, TAC_wide) %>%
    dplyr::select(!c(Target, TAC)) %>%
    filter(is.infinite(YFS_dfish) == FALSE, is.infinite(RS_dfish)==FALSE)
  
  # Immature males
  catch_im <- read.csv("./Data/rkc_gfbycatch_immale.csv") %>%
    right_join(., TAC_wide) %>%
    dplyr::select(!c(Target, TAC)) %>%
    filter(is.infinite(YFS_dfish) == FALSE, is.infinite(RS_dfish)==FALSE)
  
  # Mature females
  catch_mf <- read.csv("./Data/rkc_gfbycatch_mfem.csv") %>%
    right_join(., TAC_wide) %>%
    dplyr::select(!c(Target, TAC)) %>%
    filter(is.infinite(YFS_dfish) == FALSE, is.infinite(RS_dfish)==FALSE)
  
  # Immature females
  catch_imf <- read.csv("./Data/rkc_gfbycatch_imfem.csv") %>%
    right_join(., TAC_wide) %>%
    dplyr::select(!c(Target, TAC)) %>%
    filter(is.infinite(YFS_dfish) == FALSE, is.infinite(RS_dfish)==FALSE)
  
 # Write function to process target catch data into spatial rasters
  dfish_rasts <- function(catch, predict_yr, period, preds){
    # Specify month #s for filtering
    if(period == "Jan/Feb"){
      mths = 1:2
    } else if(period == "Apr/May"){
      mths = 4:5
    } else{
      mths = 9:10
    }
    
    data2 <- catch %>%
      filter(month %in% mths) %>%
      rename(y = lat, x = lon) %>%
      st_as_sf(coords = c("x", "y"), crs = in.crs) %>%
      st_transform(., map.crs)
    
    as.data.frame(data2) %>% dplyr::select(!c(X, geometry)) -> ll
    
    data3 <- cbind(st_coordinates(data2), ll) %>%
      rename(x = X, y = Y)
    
    for(ii in 1:length(predict_yr)){
      # Create rasters for yfs and rs directed fishing counts
      data3 %>%
        filter(year == predict_yr[ii]) %>%
        dplyr::select(x, y, YFS_dfish) -> yfs_dfish
      
      data3 %>%
        filter(year == predict_yr[ii]) %>%
        dplyr::select(x, y, RS_dfish) -> rs_dfish
      
      # Create empty raster of BB strata
      BB_rast <- rast(st_as_sf(BB_strata), res = 4971.196)
      ext(BB_rast) <- ext(preds)
      BB_rast2 <- raster(BB_rast)
      
      # Fill empty BB rast with yfs and rs directed fishing counts
      YFSdfish_rast <- rasterize(yfs_dfish[,c("x", "y")], raster(preds), yfs_dfish$YFS_dfish, sum, background = 0) %>%
        rast() %>%
        mask(BB_strata) 
      
      RSdfish_rast <- rasterize(rs_dfish[,c("x", "y")], raster(preds), rs_dfish$RS_dfish, sum, background = 0) %>%
        rast() %>%
        mask(BB_strata)
      
      c(YFSdfish_rast, RSdfish_rast) -> dfish_rast
      
      names(dfish_rast) = paste(period, c("YFS_dfish", "RS_dfish"), predict_yr[ii])
      
      rasts = c(dfish_rast, rasts)
      
    }
    return(rasts)
  }
  
 # Run function for each sex/size-maturity category (yellowfin sole and rock sole catch should be = for each)
  # Legal males
    predict_yr <- c(1997:2019, 2021:2023)
    preds <- lm_preds
    catch <- catch_lm
    rasts <- c()
    
    # Run function
    c("Jan/Feb", "Apr/May", "Sep/Oct") %>%
      purrr::map(~dfish_rasts(catch, predict_yr, .x, preds)) -> out
    
    # Join directed fishery rasters with rest of predictors, save
    lm_preds <- c(preds, rast(out))
    
    #writeRaster(lm_preds, "./Data/lm_preds_df.tif", overwrite = TRUE)
    
  # Immature males
    preds <- im_preds
    catch <- catch_im
    rasts <- c()
    
    c("Jan/Feb", "Apr/May", "Sep/Oct") %>%
      purrr::map(~dfish_rasts(catch, predict_yr, .x, preds)) -> out.im
    
    # Join directed fishery rasters with rest of predictors
    im_preds <- c(preds, rast(out.im))
    
    #writeRaster(im_preds, "./Data/im_preds_df.tif", overwrite = TRUE)
    
  # Mature females
    preds <- mf_preds
    catch <- catch_mf
    rasts <- c()
    
    c("Jan/Feb", "Apr/May", "Sep/Oct") %>%
      purrr::map(~dfish_rasts(catch, predict_yr, .x, preds)) -> out.mf
    
    # Join directed fishery rasters with rest of predictors
    mf_preds <- c(preds, rast(out.mf))
    
    #writeRaster(mf_preds, "./Data/mf_preds_df.tif", overwrite = TRUE)
    
  # Immature females
    preds <- imf_preds
    catch <- catch_imf
    rasts <- c()
    
    c("Jan/Feb", "Apr/May", "Sep/Oct") %>%
      purrr::map(~dfish_rasts(catch, predict_yr, .x, preds)) -> out.imf
    
    # Join directed fishery rasters with rest of predictors
    imf_preds <- c(preds, rast(out.imf))
    
    #writeRaster(imf_preds, "./Data/imf_preds_df.tif", overwrite = TRUE)
    
 