#PURPOSE:
# To spatially predict to BBRKC bycatch across trawlable area in Bristol Bay

#AUTHOR:
# Emily Ryznar - NOAA/NMFS RACE Shellfish Assessment Program (emily.ryznar@noaa.gov)

### LOAD PROCESSING PARAMETERS --------------------------------------------------
source("./Scripts/load_libs_params.R")

### WRITE FUNCTION TO GENERATURE PREDICTION RASTERS -----------------------------
predict_rast<- function(preds, model_b, model_p, train, test, period, predict_yr){
  
  for(ii in 1:length(predict_yr)){

    # Filter data (train or test) by prediction year and period
    data2 <- rbind(train, test) %>%
      filter(predict_year == predict_yr[ii], Period == period)
    
    # Specify past year based on present year and prediction period for subsetting
    pst <- ifelse(predict_yr[ii] == 2021, c(predict_yr[ii]-2, predict_yr[-length(predict_yr[ii])]), 
                  c(predict_yr[ii]-1, predict_yr[-length(predict_yr[ii])]))
    
    # Subset predictors into past and present predictors, specify labels based on year and prediction period
    if (period == "Sep/Oct"){
      preds <- subset(preds, grep("Jan/Feb |Apr/May ", names(preds), invert = TRUE))
      preds_past <- subset(preds, grep("Sep_Oct SST |Nov_Dec SST ", names(preds)))
      preds_pres <- subset(preds, grep("Sep_Oct SST |Nov_Dec SST ", names(preds), invert = TRUE))
      
      labs <- c("Sep_Oct_SST", "Nov_Dec_SST", "Jan_Feb_SST", "Mar_Apr_SST", "May_Jun_SST", "Jul_Aug_SST", "Jan_Feb_Ice",
                "Mar_Apr_Ice", "May_Aug_BT", "Sed", "SAP_Count", "Depth", "YFS_CPUE", "RS_CPUE", "YFS_dfish", "RS_dfish")
      
      preds2 <- c(terra::subset(preds_past, grep(pst, names(preds_past), value = T)),
                  terra::subset(preds_pres, grep(predict_yr[ii], names(preds_pres))))
      
    } else if (period == "Apr/May"){
      preds <- subset(preds, grep("Jan/Feb |Sep/Oct ", names(preds), invert = TRUE))
      
      preds_past <- subset(preds, grep("May_Jun SST |Jul_Aug SST |Sep_Oct SST |Nov_Dec SST |SAP Count |BT |YFS CPUE |RS CPUE ", 
                                       names(preds)))
      preds_pres <- subset(preds, grep("May_Jun SST |Jul_Aug SST |Sep_Oct SST |Nov_Dec SST |SAP Count |BT |YFS CPUE |RS CPUE ", 
                                       names(preds), invert = TRUE))
      
      labs <- c("May_Jun_SST", "Jul_Aug_SST", "Sep_Oct_SST", "Nov_Dec_SST", "May_Aug_BT", "SAP_Count", "YFS_CPUE", 
                "RS_CPUE", "Jan_Feb_SST", "Mar_Apr_SST", "Jan_Feb_Ice", "Mar_Apr_Ice", "Sed", "Depth", "YFS_dfish", "RS_dfish")
      
      preds2 <- c(terra::subset(preds_past, grep(pst, names(preds_past), value = T)),
                  terra::subset(preds_pres, grep(predict_yr[ii], names(preds_pres))))
      
    } else{
      preds <- subset(preds, grep("Apr/May |Sep/Oct ", names(preds), invert = TRUE))
      
      preds_past <- subset(preds, grep("Jan/Feb ", names(preds), invert = TRUE))
      
      preds_pres <- subset(preds, grep("Jan/Feb ", names(preds)))
      
      labs <- c("Jan_Feb_SST", "Mar_Apr_SST", "May_Jun_SST", "Jul_Aug_SST", "Sep_Oct_SST", "Nov_Dec_SST",
                "Jan_Feb_Ice", "Mar_Apr_Ice", "May_Aug_BT", "Sed", "SAP_Count", "Depth", "YFS_CPUE", "RS_CPUE", "YFS_dfish", "RS_dfish")
      
      preds2 <- c(terra::subset(preds_past, grep(pst, names(preds_past), value = T)),
                  terra::subset(preds_pres, grep(predict_yr[ii], names(preds_pres))))
    }
    
    # Set names of preds
    names(preds2) <- labs
    
    # specify closure areas by year
    if(predict_yr[ii] %in% c(1996, 2022:2023) & period %in% c("Jan/Feb", "Sep/Oct")){ # rkc df fishery closed
      mm <- janfebsepoct_closure_wsub
    } else if(predict_yr[ii] %in% c(1997, 2003) & period %in% c("Jan/Feb", "Sep/Oct")){
      mm <- BLZ1
    } else if(predict_yr[ii] %in% c(1996, 2022:2023) & period == "Apr/May"){ # rkc df fishery closed
      mm <- aprmay_closure_wsub
    } else if(predict_yr[ii] %in% c(1997:2021) & period == "Apr/May"){
      mm <- aprmay_closure
    } else{
      mm <- janfebsepoct_closure
    }
    
    # mask preds by closure areas
    preds3 <- preds2 %>%
      mask(mm) #Closure area by year

    
    if(period == "Jan/Feb"){
      mnths = 1:2
    } else if(period == "Apr/May"){
      mnths = 4:5
    } else{
      mnths = 9:10
    }

    #set up plotting features
    panel_extent <- data.frame(y = c(54, 60), #Bristol Bay #was 53.5
                               x = c(-167.5, -158)) %>%
      akgfmaps::transform_data_frame_crs(out.crs = map.crs) #crs was coldpool:::ebs_proj_crs
    
    
    # Specify TAC constant
    data2$ELV_SWP <- as.factor(data2$ELV_SWP)
    data2$Period <- as.factor(data2$Period)
    
    add <- data.frame(RS_TAC = NA, YFS_TAC = NA, ELV_SWP = NA, Period = NA)
    
    add$RS_TAC <- data2 %>%
      pull(RS_TAC) %>%
      unique()
    
    add$YFS_TAC <- data2 %>%
      pull(YFS_TAC) %>%
      unique()
    
    add$ELV_SWP <- data2 %>%
      pull(ELV_SWP) %>%
      unique()
    
    add$Period <- data2 %>%
      pull(Period) %>%
      unique
    
    
    # # Create rasters for yfs and rs directed fishing counts
    # data2 %>%
    #   dplyr::select(x, y, YFS_dfish) -> yfs_dfish
    # 
    # data2 %>%
    #   dplyr::select(x, y, RS_dfish) -> rs_dfish
    # 
    # # Create empty raster of BB strata
    # BB_rast <- rast(st_as_sf(BB_strata), res = 4971.196)
    # ext(BB_rast) <- ext(preds2)
    # BB_rast2 <- raster(BB_rast)
    # 
    # # Fill empty BB rast with yfs and rs directed fishing counts
    # YFSdfish_rast <- rasterize(yfs_dfish[,c("x", "y")], raster(preds2), yfs_dfish$YFS_dfish, sum, background = 0) %>%
    #   rast() %>%
    #   mask(BB_strata) %>%
    #   mask(mm)
    # 
    # RSdfish_rast <- rasterize(rs_dfish[,c("x", "y")], raster(preds2), rs_dfish$RS_dfish, sum, background = 0) %>%
    #   rast() %>%
    #   mask(BB_strata) %>%
    #   mask(mm)
    # 
    # # Combine yfs and rs rasters with rest of raster stack, set name
    # c(preds3, as.numeric(YFSdfish_rast), as.numeric(RSdfish_rast)) -> preds4 #changing TAC to factor here reducing prob of occurance but removes NAs
    # names(preds4)[15:16] <- c("YFS_dfish", "RS_dfish")
    # 
    # 
    # newdat <- cbind(crds(preds4), na.omit(as.data.frame(preds4))) %>%
    #             mutate(RS_TAC = add$RS_TAC, YFS_TAC = add$YFS_TAC, ELV_SWP = add$ELV_SWP, Period = add$Period)
    # 
    # Generate spatial model predictions for Bristol Bay management area extent using raster predictors
    #PA
    
    spatpred_b <- predict(preds3, # raster stack 
                          model_b, # fitted model
                          n.trees=model_b$gbm.call$best.trees, # see help
                          const = add,
                          #factors = ff,
                          na.rm = TRUE,
                          type="response", 
                          ext = panel_extent) %>%
      mask(BB_strata)  %>%
      mask(mm) %>%
      na.omit()
    
    
    #Abund
    spatpred_p <- predict(preds3, # raster stack 
                          model_p, # fitted model
                          n.trees=model_p$gbm.call$best.trees, # see help
                          const = add,
                          na.rm = TRUE,
                          #factors = ff,
                          type="response", 
                          ext = panel_extent) %>%
      mask(BB_strata)  %>%
      mask(mm)%>%
      na.omit()
    
    thres <- mean(as.data.frame(spatpred_b)$lyr1)
    
    spatpred_b_df <- cbind(crds(spatpred_b$lyr1), as.data.frame(spatpred_b$lyr1)) %>%
      dplyr::rename("log_count" = lyr1) %>%
      mutate(year = predict_yr[ii], period = period) %>%
      mutate(PA = ifelse(log_count > thres, 1, 0))
    
    spatpred_p_df <- cbind(crds(spatpred_p$lyr1), as.data.frame(spatpred_p$lyr1)) %>%
      dplyr::rename("log_count" = lyr1) %>%
      mutate(year = predict_yr[ii], period = period) 
    
    pred_df <- spatpred_p_df %>%
      mutate(count = log_count*spatpred_b_df$PA)
    
    pred_rast <- spatpred_b$lyr1*spatpred_p$lyr1
    names(pred_rast) <- paste0("Count.", predict_yr[ii], ".", period)
    
    spatpred_df <- rbind(spatpred_df, pred_df)
    spatpred_rast <- c(pred_rast)
    spatpred_b_rast <- c(spatpred_b)
    spatpred_p_rast <- c(spatpred_p)
    
  }
  return(list(spatpred_df = spatpred_df, spatpred_b_rast, spatpred_p_rast, spatpred_rast)) 
  #spatpred_rast = spatpred_rast, spatpred_p = spatpred_p_rast, 
  #spatpred_b = spatpred_b_rast))
}

predict_yr <- c(1997:2019, 2021:2023)

### RUN FUNCTION ----------------------------------------------------------------
 # LEGAL MALES
   # Set parameters
    model_b <- readRDS("./Models/BRT_legalmale_modelb_BEST.CPUE.rda")
    model_p <- readRDS("./Models/BRT_legalmale_modelp_BEST.CPUE.rda")
  
    preds <- rast("./Data/lm_preds_df.tif")
    spatpred_df <- data.frame()
    spatpred_rast <- list()
    train <- lm_train <- read.csv("./Data/legalmale_trainCPUE.csv") %>%
                filter(iter == lm_iter)
    test <- lm_test <- read.csv("./Data/legalmale_testCPUE.csv") %>%
                filter(iter == lm_iter)
    
  # Run function
    c("Jan/Feb", "Apr/May", "Sep/Oct") %>%
    purrr::map(~predict_rast(preds, model_b, model_p, train, test, .x, predict_yr)) -> out

   
  # By phase 
    rbind(out[[1]]$spatpred_df, out[[2]]$spatpred_df, out[[3]]$spatpred_df) %>%
      mutate(phase = case_when((year %in% 1997:2005) ~ "1997:2005",
                               (year %in% 2006:2014) ~ "2006:2014",
                               (year %in% 2015:2023) ~ "2015:2023")) %>%
      group_by(phase, x, y) %>%
      reframe(mean_count = mean(count)) -> df


    pp <- c("1997:2005", "2006:2014", "2015:2023")

  # Set up list to store maps
    spatpred_list <- list()
  
  # Calculate encounter percentiles
    for(ii in 1:length(pp)){
      spatpred_df_yr <-  df %>%
        filter(phase %in% pp[ii]) %>%
        group_by(x, y) %>%
        reframe(mean_count = mean(mean_count))
      
      # Find plotting breaks
      quantiles = c(.05, .25, .5, .75)
      quants<-sort(unique(c(0,quantiles,1)))
      
      threshold <- 0.0513
      
      sample <- stats::na.omit(spatpred_df_yr$mean_count)
      sample[sample <= threshold] <- NA
      perc.breaks <- stats::quantile(sample, probs = quants, na.rm = TRUE, names = FALSE)
      perc.breaks[1]<-0
      perc.breaks[length(perc.breaks)]<-Inf
      
      # Make raster again
      spatpred_df_yr %>%
        rast() %>%
        raster() -> spatpred_yr
      
      # Set crs
      crs(spatpred_yr) <- map.crs
      
      # Cut the prediction map by the breaks
      perc.map <- raster::cut(spatpred_yr, breaks = perc.breaks)
      
      # set up the factor maps
      perc.vals <- raster::getValues(perc.map)
      perc.vals[perc.vals == 1] <- NA
      
      # convert the raster to polygons
      percpoly0 <- stars::st_as_stars(perc.map)
      percpoly <- sf::st_as_sf(percpoly0,merge = TRUE)
      percpoly2 <- percpoly[percpoly$layer != 1, ]
      
      # we'll need a new outline
      perc.dummy.raster <- raster::raster(perc.map)
      perc.vals2 <- is.na(perc.vals) == F
      perc.dummy.raster <- raster::setValues(perc.dummy.raster, values = perc.vals2)
      
      percdummy0 <- stars::st_as_stars(perc.dummy.raster)
      percdummy <- sf::st_cast(sf::st_as_sf(percdummy0, merge = TRUE))
      percdummy2 <- sf::st_transform(percdummy, sf::st_crs(map.crs))
      
      # Dropping the smallest areas
      percdummy.poly <- sf::st_cast(percdummy2, "POLYGON")
      areas <- sf::st_area(percdummy.poly)
      
      outside <- order(areas, decreasing = T)[1]
      toosmall <- which(as.numeric(areas) < 10^8)
      
      perc.x <- percdummy2$layer[-c(outside, toosmall)]
      perc.y <- percdummy2$geometry[-c(outside, toosmall)]
      percdummy3 <- sf::st_sf(perc.x, perc.y) %>%
        vect() %>%
        crop(BB_strata)
      
      # Set up plot boundary
      plot.boundary.untrans <- data.frame(y = c(54.25, 59.25),
                                          x = c(-167.5, -158)) # plot boundary unprojected
      
      plot.boundary <- plot.boundary.untrans %>%
        sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
        sf::st_transform(crs = map.crs) %>%
        sf::st_coordinates() %>%
        as.data.frame() %>%
        dplyr::rename(x = X, y = Y) # plot boundary projected
      
      
      # Set up year label size
      size = 3.5
      lw = 1
      
      if(pp[ii] == "1997:2005"){
        labs = "1997-\n2005"
      }else if(pp[ii] == "2006:2014"){
        labs = "2006-\n2014"
      }else{
        labs = "2015-\n2023"
      }
      
    
      year_untrans <- data.frame(lab = labs, x = -158.3, y = 55.3)
      
      year_lab <- year_untrans %>%
        sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
        sf::st_transform(crs = map.crs) %>%
        cbind(st_coordinates(.)) %>%
        as.data.frame()
      
      
      # Map
      ggplot2::ggplot() +
        #ggplot2::geom_sf(data = survey.sf, fill = "grey95")+
        ggplot2::geom_sf(data = percpoly2, ggplot2::aes(fill = as.factor(layer)), col = NA) +
        ggplot2::geom_sf(data = st_as_sf(percdummy3),fill=NA, size = .3) +
        ggplot2::geom_sf(data = st_as_sf(ns_trawl),
                         fill = NA,
                         color = "darkturquoise",
                         linewidth = lw)+
        ggplot2::geom_sf(data = st_as_sf(NBBTCA),
                         fill = NA,
                         color = "darkturquoise",
                         linewidth = lw)+
        ggplot2::geom_sf(data = st_as_sf(fivesxtn),
                         fill = NA,
                         color = "darkturquoise",
                         linewidth = lw)+
        ggplot2::geom_sf(data = st_as_sf(RKCSA_sub),
                         fill = NA,
                         color = "darkturquoise",
                         linewidth = lw)+
        ggplot2::geom_sf(data = st_as_sf(RKCSA),
                         fill = NA,
                         color = "darkturquoise",
                         linewidth = lw)+
        ggplot2::geom_sf(data = st_as_sf(BB_strata),
                         fill = NA,
                         color = "black",
                         linewidth = 1.75)+
        ggplot2::geom_sf(data = region_layers$akland,
                         fill = "grey70",
                         color = "black")+
        
        geom_text(data = year_lab, aes(x=X, y=Y, label= lab), fontface = "bold", size=size) +
        #ggtitle(title)+
        
        coord_sf(xlim = plot.boundary$x,
                 ylim = plot.boundary$y)+
        viridis::scale_fill_viridis(option = "cividis", discrete = T, name = "Percentiles", labels = c("95%", "75%", "50%", "25%")) +
        ggplot2::theme_bw() +
        ggplot2:: theme(
          panel.border = ggplot2::element_rect(color = "black", fill = NA),
          panel.background = ggplot2::element_rect(fill = NA, color = "black"),
          legend.key = ggplot2::element_rect(fill = NA, color = "grey30"),
          #legend.position = legend.pos,
          panel.grid.major = element_blank(),
          axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 10),
          legend.text = ggplot2::element_text(size = 11), legend.title = ggplot2::element_text(size = 11),
          legend.position = "bottom", plot.title = element_text(size = 18),
          plot.background = ggplot2::element_rect(fill = "white", color = "white")) -> spatpred_list[[ii]]
      
    }
  
  # Arrange plots
  ggarrange(spatpred_list[[1]],
            spatpred_list[[2]],
            spatpred_list[[3]],
            nrow=1, ncol=3, common.legend = TRUE, legend = "bottom") -> phase.perc_rast

  ggsave(plot = phase.perc_rast, "./Figures/legalemalebycatch.phaserast.png", height=3, width=8.5, units="in")

 # IMMATURE MALES
   # Set parameters
    model_b <- readRDS("./Models/BRT_immaturemale_modelb_BEST.CPUE.rda")
    model_p <- readRDS("./Models/BRT_immaturemale_modelp_BEST.CPUE.rda")
    
    preds <- rast("./Data/im_preds_df.tif")
    spatpred_df <- data.frame()
    spatpred_rast <- list()
    train <- im_train <- read.csv("./Data/immaturemale_trainCPUE.csv") %>%
      filter(iter == im_iter)
    test <- im_test <- read.csv("./Data/immaturemale_testCPUE.csv") %>%
      filter(iter == im_iter)
  
    # Run function
    c("Jan/Feb", "Apr/May", "Sep/Oct") %>%
      purrr::map(~predict_rast(preds, model_b, model_p, train, test, .x, predict_yr)) -> out
    
    # By phase 
    rbind(out[[1]]$spatpred_df, out[[2]]$spatpred_df, out[[3]]$spatpred_df) %>%
      mutate(phase = case_when((year %in% 1997:2005) ~ "1997:2005",
                               (year %in% 2006:2014) ~ "2006:2014",
                               (year %in% 2015:2023) ~ "2015:2023")) %>%
      group_by(phase, x, y) %>%
      reframe(mean_count = mean(count)) -> df
    
    pp <- c("1997:2005", "2006:2014", "2015:2023")
  
  # Set up list to store maps
    spatpred_list <- list()
  
  # Calculate encounter percentiles
    for(ii in 1:length(pp)){
      spatpred_df_yr <-  df %>%
        filter(phase %in% pp[ii]) %>%
        group_by(x, y) %>%
        reframe(mean_count = mean(mean_count))
      
      # Find plotting breaks
      quantiles = c(.05, .25, .5, .75)
      quants<-sort(unique(c(0,quantiles,1)))
      
      threshold <- 0.0513
      
      sample <- stats::na.omit(spatpred_df_yr$mean_count)
      sample[sample <= threshold] <- NA
      perc.breaks <- stats::quantile(sample, probs = quants, na.rm = TRUE, names = FALSE)
      perc.breaks[1]<-0
      perc.breaks[length(perc.breaks)]<-Inf
      
      # Make raster again
      spatpred_df_yr %>%
        rast() %>%
        raster() -> spatpred_yr
      
      # Set crs
      crs(spatpred_yr) <- map.crs
      
      # Cut the prediction map by the breaks
      perc.map <- raster::cut(spatpred_yr, breaks = perc.breaks)
      
      # set up the factor maps
      perc.vals <- raster::getValues(perc.map)
      perc.vals[perc.vals == 1] <- NA
      
      # convert the raster to polygons
      percpoly0 <- stars::st_as_stars(perc.map)
      percpoly <- sf::st_as_sf(percpoly0,merge = TRUE)
      percpoly2 <- percpoly[percpoly$layer != 1, ]
      
      # we'll need a new outline
      perc.dummy.raster <- raster::raster(perc.map)
      perc.vals2 <- is.na(perc.vals) == F
      perc.dummy.raster <- raster::setValues(perc.dummy.raster, values = perc.vals2)
      
      percdummy0 <- stars::st_as_stars(perc.dummy.raster)
      percdummy <- sf::st_cast(sf::st_as_sf(percdummy0, merge = TRUE))
      percdummy2 <- sf::st_transform(percdummy, sf::st_crs(map.crs))
      
      # Dropping the smallest areas
      percdummy.poly <- sf::st_cast(percdummy2, "POLYGON")
      areas <- sf::st_area(percdummy.poly)
      
      outside <- order(areas, decreasing = T)[1]
      toosmall <- which(as.numeric(areas) < 10^8)
      
      perc.x <- percdummy2$layer[-c(outside, toosmall)]
      perc.y <- percdummy2$geometry[-c(outside, toosmall)]
      percdummy3 <- sf::st_sf(perc.x, perc.y) %>%
        vect() %>%
        crop(BB_strata)
      
      # Set up plot boundary
      plot.boundary.untrans <- data.frame(y = c(54.25, 59.25),
                                          x = c(-167.5, -158)) # plot boundary unprojected
      
      plot.boundary <- plot.boundary.untrans %>%
        sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
        sf::st_transform(crs = map.crs) %>%
        sf::st_coordinates() %>%
        as.data.frame() %>%
        dplyr::rename(x = X, y = Y) # plot boundary projected
      
      
      # Set up year label size
      size = 3.5
      lw = 1
      
      if(pp[ii] == "1997:2005"){
        labs = "1997-\n2005"
      }else if(pp[ii] == "2006:2014"){
        labs = "2006-\n2014"
      }else{
        labs = "2015-\n2023"
      }
      
      
      year_untrans <- data.frame(lab = labs, x = -158.3, y = 55.3)
      
      year_lab <- year_untrans %>%
        sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
        sf::st_transform(crs = map.crs) %>%
        cbind(st_coordinates(.)) %>%
        as.data.frame()
      
      
      # Map
      ggplot2::ggplot() +
        #ggplot2::geom_sf(data = survey.sf, fill = "grey95")+
        ggplot2::geom_sf(data = percpoly2, ggplot2::aes(fill = as.factor(layer)), col = NA) +
        ggplot2::geom_sf(data = st_as_sf(percdummy3),fill=NA, size = .3) +
        ggplot2::geom_sf(data = st_as_sf(ns_trawl),
                         fill = NA,
                         color = "darkturquoise",
                         linewidth = lw)+
        ggplot2::geom_sf(data = st_as_sf(NBBTCA),
                         fill = NA,
                         color = "darkturquoise",
                         linewidth = lw)+
        ggplot2::geom_sf(data = st_as_sf(fivesxtn),
                         fill = NA,
                         color = "darkturquoise",
                         linewidth = lw)+
        ggplot2::geom_sf(data = st_as_sf(RKCSA_sub),
                         fill = NA,
                         color = "darkturquoise",
                         linewidth = lw)+
        ggplot2::geom_sf(data = st_as_sf(RKCSA),
                         fill = NA,
                         color = "darkturquoise",
                         linewidth = lw)+
        ggplot2::geom_sf(data = st_as_sf(BB_strata),
                         fill = NA,
                         color = "black",
                         linewidth = 1.75)+
        ggplot2::geom_sf(data = region_layers$akland,
                         fill = "grey70",
                         color = "black")+
        
        geom_text(data = year_lab, aes(x=X, y=Y, label= lab), fontface = "bold", size=size) +
        #ggtitle(title)+
        
        coord_sf(xlim = plot.boundary$x,
                 ylim = plot.boundary$y)+
        viridis::scale_fill_viridis(option = "cividis", discrete = T, name = "Percentiles", labels = c("95%", "75%", "50%", "25%")) +
        ggplot2::theme_bw() +
        ggplot2:: theme(
          panel.border = ggplot2::element_rect(color = "black", fill = NA),
          panel.background = ggplot2::element_rect(fill = NA, color = "black"),
          legend.key = ggplot2::element_rect(fill = NA, color = "grey30"),
          #legend.position = legend.pos,
          panel.grid.major = element_blank(),
          axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 10),
          legend.text = ggplot2::element_text(size = 11), legend.title = ggplot2::element_text(size = 11),
          legend.position = "bottom", plot.title = element_text(size = 18),
          plot.background = ggplot2::element_rect(fill = "white", color = "white")) -> spatpred_list[[ii]]
      
    }
  
  # Arrange plots
    ggarrange(spatpred_list[[1]],
              spatpred_list[[2]],
              spatpred_list[[3]],
              nrow=1, ncol=3, common.legend = TRUE, legend = "bottom") -> phase.perc_rast
    
    ggsave(plot = phase.perc_rast, "./Figures/immaturemalebycatch.phaserast.png", height=3, width=8.5, units="in")
    
 # MATURE FEMALES
   # Set parameters
    model_b <- readRDS("./Models/BRT_maturefemale_modelb_BEST.CPUE.rda")
    model_p <- readRDS("./Models/BRT_maturefemale_modelp_BEST.CPUE.rda")
    
    preds <- rast("./Data/mf_preds_df.tif")
    spatpred_df <- data.frame()
    spatpred_rast <- list()
    train <- mf_train <- read.csv("./Data/maturefemale_trainCPUE.csv") %>%
      filter(iter == mf_iter)
    test <- mf_test <- read.csv("./Data/maturefemale_testCPUE.csv") %>%
      filter(iter == mf_iter)
    
  # Run function
    c("Jan/Feb", "Apr/May", "Sep/Oct") %>%
      purrr::map(~predict_rast(preds, model_b, model_p, train, test, .x, predict_yr)) -> out
    
  # By phase 
    rbind(out[[1]]$spatpred_df, out[[2]]$spatpred_df, out[[3]]$spatpred_df) %>%
      mutate(phase = case_when((year %in% 1997:2005) ~ "1997:2005",
                               (year %in% 2006:2014) ~ "2006:2014",
                               (year %in% 2015:2023) ~ "2015:2023")) %>%
      group_by(phase, x, y) %>%
      reframe(mean_count = mean(count)) -> df
    
    
    pp <- c("1997:2005", "2006:2014", "2015:2023")
    
  # Set up list to store maps
    spatpred_list <- list()
  
  # Calculate encounter percentiles
    for(ii in 1:length(pp)){
      spatpred_df_yr <-  df %>%
        filter(phase %in% pp[ii]) %>%
        group_by(x, y) %>%
        reframe(mean_count = mean(mean_count))
      
      # Find plotting breaks
      quantiles = c(.05, .25, .5, .75)
      quants<-sort(unique(c(0,quantiles,1)))
      
      threshold <- 0.0513
      
      sample <- stats::na.omit(spatpred_df_yr$mean_count)
      sample[sample <= threshold] <- NA
      perc.breaks <- stats::quantile(sample, probs = quants, na.rm = TRUE, names = FALSE)
      perc.breaks[1]<-0
      perc.breaks[length(perc.breaks)]<-Inf
      
      # Make raster again
      spatpred_df_yr %>%
        rast() %>%
        raster() -> spatpred_yr
      
      # Set crs
      crs(spatpred_yr) <- map.crs
      
      # Cut the prediction map by the breaks
      perc.map <- raster::cut(spatpred_yr, breaks = perc.breaks)
      
      # set up the factor maps
      perc.vals <- raster::getValues(perc.map)
      perc.vals[perc.vals == 1] <- NA
      
      # convert the raster to polygons
      percpoly0 <- stars::st_as_stars(perc.map)
      percpoly <- sf::st_as_sf(percpoly0,merge = TRUE)
      percpoly2 <- percpoly[percpoly$layer != 1, ]
      
      # we'll need a new outline
      perc.dummy.raster <- raster::raster(perc.map)
      perc.vals2 <- is.na(perc.vals) == F
      perc.dummy.raster <- raster::setValues(perc.dummy.raster, values = perc.vals2)
      
      percdummy0 <- stars::st_as_stars(perc.dummy.raster)
      percdummy <- sf::st_cast(sf::st_as_sf(percdummy0, merge = TRUE))
      percdummy2 <- sf::st_transform(percdummy, sf::st_crs(map.crs))
      
      # Dropping the smallest areas
      percdummy.poly <- sf::st_cast(percdummy2, "POLYGON")
      areas <- sf::st_area(percdummy.poly)
      
      outside <- order(areas, decreasing = T)[1]
      toosmall <- which(as.numeric(areas) < 10^8)
      
      perc.x <- percdummy2$layer[-c(outside, toosmall)]
      perc.y <- percdummy2$geometry[-c(outside, toosmall)]
      percdummy3 <- sf::st_sf(perc.x, perc.y) %>%
        vect() %>%
        crop(BB_strata)
      
      # Set up plot boundary
      plot.boundary.untrans <- data.frame(y = c(54.25, 59.25),
                                          x = c(-167.5, -158)) # plot boundary unprojected
      
      plot.boundary <- plot.boundary.untrans %>%
        sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
        sf::st_transform(crs = map.crs) %>%
        sf::st_coordinates() %>%
        as.data.frame() %>%
        dplyr::rename(x = X, y = Y) # plot boundary projected
      
      
      # Set up year label size
      size = 3.5
      lw = 1
      
      if(pp[ii] == "1997:2005"){
        labs = "1997-\n2005"
      }else if(pp[ii] == "2006:2014"){
        labs = "2006-\n2014"
      }else{
        labs = "2015-\n2023"
      }
      
      
      year_untrans <- data.frame(lab = labs, x = -158.3, y = 55.3)
      
      year_lab <- year_untrans %>%
        sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
        sf::st_transform(crs = map.crs) %>%
        cbind(st_coordinates(.)) %>%
        as.data.frame()
      
      
      # Map
      ggplot2::ggplot() +
        #ggplot2::geom_sf(data = survey.sf, fill = "grey95")+
        ggplot2::geom_sf(data = percpoly2, ggplot2::aes(fill = as.factor(layer)), col = NA) +
        ggplot2::geom_sf(data = st_as_sf(percdummy3),fill=NA, size = .3) +
        ggplot2::geom_sf(data = st_as_sf(ns_trawl),
                         fill = NA,
                         color = "darkturquoise",
                         linewidth = lw)+
        ggplot2::geom_sf(data = st_as_sf(NBBTCA),
                         fill = NA,
                         color = "darkturquoise",
                         linewidth = lw)+
        ggplot2::geom_sf(data = st_as_sf(fivesxtn),
                         fill = NA,
                         color = "darkturquoise",
                         linewidth = lw)+
        ggplot2::geom_sf(data = st_as_sf(RKCSA_sub),
                         fill = NA,
                         color = "darkturquoise",
                         linewidth = lw)+
        ggplot2::geom_sf(data = st_as_sf(RKCSA),
                         fill = NA,
                         color = "darkturquoise",
                         linewidth = lw)+
        ggplot2::geom_sf(data = st_as_sf(BB_strata),
                         fill = NA,
                         color = "black",
                         linewidth = 1.75)+
        ggplot2::geom_sf(data = region_layers$akland,
                         fill = "grey70",
                         color = "black")+
        
        geom_text(data = year_lab, aes(x=X, y=Y, label= lab), fontface = "bold", size=size) +
        #ggtitle(title)+
        
        coord_sf(xlim = plot.boundary$x,
                 ylim = plot.boundary$y)+
        viridis::scale_fill_viridis(option = "cividis", discrete = T, name = "Percentiles", labels = c("95%", "75%", "50%", "25%")) +
        ggplot2::theme_bw() +
        ggplot2:: theme(
          panel.border = ggplot2::element_rect(color = "black", fill = NA),
          panel.background = ggplot2::element_rect(fill = NA, color = "black"),
          legend.key = ggplot2::element_rect(fill = NA, color = "grey30"),
          #legend.position = legend.pos,
          panel.grid.major = element_blank(),
          axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 10),
          legend.text = ggplot2::element_text(size = 11), legend.title = ggplot2::element_text(size = 11),
          legend.position = "bottom", plot.title = element_text(size = 18),
          plot.background = ggplot2::element_rect(fill = "white", color = "white")) -> spatpred_list[[ii]]
      
    }
    
  # Arrange plots
    ggarrange(spatpred_list[[1]],
              spatpred_list[[2]],
              spatpred_list[[3]],
              nrow=1, ncol=3, common.legend = TRUE, legend = "bottom") -> phase.perc_rast
    
    ggsave(plot = phase.perc_rast, "./Figures/maturefemalebycatch.phaserast.png", height=3, width=8.5, units="in")
 
 # IMMATURE FEMALES
  # Set parameters
    model_b <- readRDS("./Models/BRT_immaturefemale_modelb_BEST.CPUE.rda")
    model_p <- readRDS("./Models/BRT_immaturefemale_modelp_BEST.CPUE.rda")
    
    preds <- rast("./Data/imf_preds_df.tif")
    spatpred_df <- data.frame()
    spatpred_rast <- list()
    train <- imf_train <- read.csv("./Data/immaturefemale_trainCPUE.csv") %>%
      filter(iter == imf_iter)
    test <- imf_test <- read.csv("./Data/immaturefemale_testCPUE.csv") %>%
      filter(iter == imf_iter)
    
  # Run function
    c("Jan/Feb", "Apr/May", "Sep/Oct") %>%
      purrr::map(~predict_rast(preds, model_b, model_p, train, test, .x, predict_yr)) -> out
    
  # By phase 
    rbind(out[[1]]$spatpred_df, out[[2]]$spatpred_df, out[[3]]$spatpred_df) %>%
      mutate(phase = case_when((year %in% 1997:2005) ~ "1997:2005",
                               (year %in% 2006:2014) ~ "2006:2014",
                               (year %in% 2015:2023) ~ "2015:2023")) %>%
      group_by(phase, x, y) %>%
      reframe(mean_count = mean(count)) -> df
    
    
    pp <- c("1997:2005", "2006:2014", "2015:2023")
    
  # Set up list to store maps
    spatpred_list <- list()
  
  # Calculate encounter percentiles
    for(ii in 1:length(pp)){
      spatpred_df_yr <-  df %>%
        filter(phase %in% pp[ii]) %>%
        group_by(x, y) %>%
        reframe(mean_count = mean(mean_count))
      
      # Find plotting breaks
      quantiles = c(.05, .25, .5, .75)
      quants<-sort(unique(c(0,quantiles,1)))
      
      threshold <- 0.0513
      
      sample <- stats::na.omit(spatpred_df_yr$mean_count)
      sample[sample <= threshold] <- NA
      perc.breaks <- stats::quantile(sample, probs = quants, na.rm = TRUE, names = FALSE)
      perc.breaks[1]<-0
      perc.breaks[length(perc.breaks)]<-Inf
      
      # Make raster again
      spatpred_df_yr %>%
        rast() %>%
        raster() -> spatpred_yr
      
      # Set crs
      crs(spatpred_yr) <- map.crs
      
      # Cut the prediction map by the breaks
      perc.map <- raster::cut(spatpred_yr, breaks = perc.breaks)
      
      # set up the factor maps
      perc.vals <- raster::getValues(perc.map)
      perc.vals[perc.vals == 1] <- NA
      
      # convert the raster to polygons
      percpoly0 <- stars::st_as_stars(perc.map)
      percpoly <- sf::st_as_sf(percpoly0,merge = TRUE)
      percpoly2 <- percpoly[percpoly$layer != 1, ]
      
      # we'll need a new outline
      perc.dummy.raster <- raster::raster(perc.map)
      perc.vals2 <- is.na(perc.vals) == F
      perc.dummy.raster <- raster::setValues(perc.dummy.raster, values = perc.vals2)
      
      percdummy0 <- stars::st_as_stars(perc.dummy.raster)
      percdummy <- sf::st_cast(sf::st_as_sf(percdummy0, merge = TRUE))
      percdummy2 <- sf::st_transform(percdummy, sf::st_crs(map.crs))
      
      # Dropping the smallest areas
      percdummy.poly <- sf::st_cast(percdummy2, "POLYGON")
      areas <- sf::st_area(percdummy.poly)
      
      outside <- order(areas, decreasing = T)[1]
      toosmall <- which(as.numeric(areas) < 10^8)
      
      perc.x <- percdummy2$layer[-c(outside, toosmall)]
      perc.y <- percdummy2$geometry[-c(outside, toosmall)]
      percdummy3 <- sf::st_sf(perc.x, perc.y) %>%
        vect() %>%
        crop(BB_strata)
      
      # Set up plot boundary
      plot.boundary.untrans <- data.frame(y = c(54.25, 59.25),
                                          x = c(-167.5, -158)) # plot boundary unprojected
      
      plot.boundary <- plot.boundary.untrans %>%
        sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
        sf::st_transform(crs = map.crs) %>%
        sf::st_coordinates() %>%
        as.data.frame() %>%
        dplyr::rename(x = X, y = Y) # plot boundary projected
      
      
      # Set up year label size
      size = 3.5
      lw = 1
      
      if(pp[ii] == "1997:2005"){
        labs = "1997-\n2005"
      }else if(pp[ii] == "2006:2014"){
        labs = "2006-\n2014"
      }else{
        labs = "2015-\n2023"
      }
      
      
      year_untrans <- data.frame(lab = labs, x = -158.3, y = 55.3)
      
      year_lab <- year_untrans %>%
        sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
        sf::st_transform(crs = map.crs) %>%
        cbind(st_coordinates(.)) %>%
        as.data.frame()
      
      
      # Map
      ggplot2::ggplot() +
        #ggplot2::geom_sf(data = survey.sf, fill = "grey95")+
        ggplot2::geom_sf(data = percpoly2, ggplot2::aes(fill = as.factor(layer)), col = NA) +
        ggplot2::geom_sf(data = st_as_sf(percdummy3),fill=NA, size = .3) +
        ggplot2::geom_sf(data = st_as_sf(ns_trawl),
                         fill = NA,
                         color = "darkturquoise",
                         linewidth = lw)+
        ggplot2::geom_sf(data = st_as_sf(NBBTCA),
                         fill = NA,
                         color = "darkturquoise",
                         linewidth = lw)+
        ggplot2::geom_sf(data = st_as_sf(fivesxtn),
                         fill = NA,
                         color = "darkturquoise",
                         linewidth = lw)+
        ggplot2::geom_sf(data = st_as_sf(RKCSA_sub),
                         fill = NA,
                         color = "darkturquoise",
                         linewidth = lw)+
        ggplot2::geom_sf(data = st_as_sf(RKCSA),
                         fill = NA,
                         color = "darkturquoise",
                         linewidth = lw)+
        ggplot2::geom_sf(data = st_as_sf(BB_strata),
                         fill = NA,
                         color = "black",
                         linewidth = 1.75)+
        ggplot2::geom_sf(data = region_layers$akland,
                         fill = "grey70",
                         color = "black")+
        
        geom_text(data = year_lab, aes(x=X, y=Y, label= lab), fontface = "bold", size=size) +
        #ggtitle(title)+
        
        coord_sf(xlim = plot.boundary$x,
                 ylim = plot.boundary$y)+
        viridis::scale_fill_viridis(option = "cividis", discrete = T, name = "Percentiles", labels = c("95%", "75%", "50%", "25%")) +
        ggplot2::theme_bw() +
        ggplot2:: theme(
          panel.border = ggplot2::element_rect(color = "black", fill = NA),
          panel.background = ggplot2::element_rect(fill = NA, color = "black"),
          legend.key = ggplot2::element_rect(fill = NA, color = "grey30"),
          #legend.position = legend.pos,
          panel.grid.major = element_blank(),
          axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 10),
          legend.text = ggplot2::element_text(size = 11), legend.title = ggplot2::element_text(size = 11),
          legend.position = "bottom", plot.title = element_text(size = 18),
          plot.background = ggplot2::element_rect(fill = "white", color = "white")) -> spatpred_list[[ii]]
      
    }
    
  # Arrange plots
    ggarrange(spatpred_list[[1]],
              spatpred_list[[2]],
              spatpred_list[[3]],
              nrow=1, ncol=3, common.legend = TRUE, legend = "bottom") -> phase.perc_rast
    
    ggsave(plot = phase.perc_rast, "./Figures/immaturefemalebycatch.phaserast.png", height=3, width=8.5, units="in")
    