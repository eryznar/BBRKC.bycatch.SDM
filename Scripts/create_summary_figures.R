# PURPOSE:
# To generate several summary plots for Bristol Bay red king crab bycatch SDMs

# AUTHOR: 
# Emily Ryznar - NOAA/NMFS RACE Shellfish Assessment Program (emily.ryznar@noaa.gov)

# LOAD PROCESSING PARAMETERS --------------------------------------------------
source("./Scripts/load_libs_params.R")

# LOAD DATA ----------------------------------------------------------------------------------------------
  # Raw data
    # Legal males
    catch_lm <- read.csv("./Data/rkc_gfbycatch.csv")
    
      TAC_wide <- catch_lm %>%
        dplyr::select(year, Target, TAC) %>%
        unique() %>%
        pivot_wider(., id_cols = year, names_from = Target, values_from = TAC) %>%
        dplyr::rename("RS_TAC" = `Rock Sole - BSAI`, "YFS_TAC" = `Yellowfin Sole - BSAI`)
    
    catch_lm <- right_join(catch_lm, TAC_wide) %>%
      dplyr::select(!c(TAC)) %>%
      dplyr::rename("total" = legal_total, "total_extrap" = legal_extrap, "total_extrap_rd" = legal_extrap_rd)
  
    # Immature males
    catch_im <- read.csv("./Data/rkc_gfbycatch_immale.csv") %>%
      right_join(., TAC_wide) %>%
      dplyr::select(!c(TAC)) %>%
      dplyr::rename("total" = immale_total, "total_extrap" = immale_extrap, "total_extrap_rd" = immale_extrap_rd)
    
    # Mature females
    catch_mf <- read.csv("./Data/rkc_gfbycatch_mfem.csv") %>%
      right_join(., TAC_wide) %>%
      dplyr::select(!c(TAC)) %>%
      dplyr::rename("total" = mfem_total, "total_extrap" = mfem_extrap, "total_extrap_rd" = mfem_extrap_rd)
  
    # Immature females
    catch_imf <- read.csv("./Data/rkc_gfbycatch_imfem.csv") %>%
      right_join(., TAC_wide) %>%
      dplyr::select(!c(TAC)) %>%
      dplyr::rename("total" = imfem_total, "total_extrap" = imfem_extrap, "total_extrap_rd" = imfem_extrap_rd)

  # Training/testing data 
    # Legal male
    lm_train <- read.csv("./Data/legalmale_trainCPUE.csv") %>%
      filter(iter == lm_iter)
    lm_test <- read.csv("./Data/legalmale_testCPUE.csv") %>%
      filter(iter == lm_iter)
    
    # Immature male
    im_train <- read.csv("./Data/immaturemale_trainCPUE.csv") %>%
      filter(iter == im_iter)
    im_test <- read.csv("./Data/immaturemale_testCPUE.csv") %>%
      filter(iter == im_iter)
    
    # Mature female
    mf_train <- read.csv("./Data/maturefemale_trainCPUE.csv") %>%
      filter(iter == mf_iter)
    mf_test <- read.csv("./Data/maturefemale_testCPUE.csv") %>%
      filter(iter == mf_iter)
    
    # Immature female
    imf_train <- read.csv("./Data/immaturefemale_trainCPUE.csv") %>%
      filter(iter == imf_iter)
    imf_test <- read.csv("./Data/immaturefemale_testCPUE.csv") %>%
      filter(iter == imf_iter)

  # Best models 
    # Legal male
    lm_modelb <- readRDS("./Models/BRT_legalmale_modelb_BEST.CPUE.rda")
    lm_modelp <- readRDS("./Models/BRT_legalmale_modelp_BEST.CPUE.rda")
    
    # Immature male
    im_modelb <- readRDS("./Models/BRT_immaturemale_modelb_BEST.CPUE.rda")
    im_modelp <- readRDS("./Models/BRT_immaturemale_modelp_BEST.CPUE.rda")
    
    # Mature female
    mf_modelb <- readRDS("./Models/BRT_maturefemale_modelb_BEST.CPUE.rda")
    mf_modelp <- readRDS("./Models/BRT_maturefemale_modelp_BEST.CPUE.rda")
    
    # Immature female
    imf_modelb <- readRDS("./Models/BRT_immaturefemale_modelb_BEST.CPUE.rda")
    imf_modelp <- readRDS("./Models/BRT_immaturefemale_modelp_BEST.CPUE.rda")

  
# SPECIFY FUNCTIONS ---------------------------------------------------------------------------------------
  # Generate prediction data frame using the testing data (to be used in sum_plots function below)
    # @param model_b: best binomial model
    # @param model_p: best poisson model
    # @param test: testing data (options = "lm_test", "im_test", "mf_test", "imf_test")
    # @return Returns a dataframe of spatial predictions at fishing locations in the testing data
    
    pred_df_test <- function(model_b, model_p, test){
    
    # Calculate predictions
    pred_b <- suppressWarnings(predict.gbm(model_b, # fitted model to predict
                          test, # data to predict to
                          n.trees=model_b$gbm.call$best.trees, # see help
                          type="response")) # predict probabilities
    
    pred_p <- suppressWarnings(predict.gbm(model_p, # fitted model to predict
                          test %>% filter(extrap_log_rd>0), # data to predict to
                          n.trees=model_p$gbm.call$best.trees, # see help
                          type="response")) # predict probabilities 
    
    #For delta model, multiply predicted output from model_b by model_p
    data.frame(predict_year = test$predict_year,
               period = test$Period,
               x = test$x,
               y = test$y,
               log_PA = pred_b) -> pred_b
    
    thres <- mean(pred_b$log_PA)
    
    pred_b %>%
      mutate(thres_PA = ifelse(log_PA > thres, 1, 0)) -> pred_b #using mean probability of predicted 
                                                            #species presence as threshold (via Liu et al. 2005)
    
    data.frame(predict_year = (test %>% filter(extrap_log_rd>0))$predict_year,
               period = (test %>% filter(extrap_log_rd>0))$Period,
               x = (test %>% filter(extrap_log_rd>0))$x,
               y = (test %>% filter(extrap_log_rd>0))$y,
               log_count = pred_p) -> pred_p
    
    right_join(pred_p, pred_b, by = c("predict_year", "period", "x", "y"), relationship = "many-to-many") %>%
      mutate(log_count = log_count * thres_PA) %>%
      replace_na(list(log_count = 0)) -> pred_df
  
    return(pred_df)
    
  }
  
  # Generate a several summary plots for model output and each sex-size/maturity category
    # @param model_b: best binomial model
    # @param model_p: best poisson model
    # @param pred_df: predictions generated using testing data generated by 'pred_df_test' function
    # @param train: training data (options = "lm_train", "im_train", "mf_train", "imf_train")
    # @param test: testing data (options = "lm_test", "im_test", "mf_test", "imf_test")
    # @param mat_sex: RKC sex-size/maturity category (options = "Legal_male", "Immature_male", "Mature female", "Immature female")
    # @param raw_dat: catch data (options = "catch_lm", "catch_im", "catch_mf", "catch_imf")
    # @param Period: prediction period over which to generate plots (options = "Jan/Feb", "Apr/May", "Sep/Oct", "All")
    # @return Returns several summary figures/metrics/dataframes, including weighted coordinates dataframe,
    # weighted coordinates COD plot, great circle distance between predicted and observed, great circle distance between CODs,
    # variable influence plots, timeseries of observed bycatch plot, top influential variables, AUC for occurrence performance,
    # rho and PDE for abundance performance
    
    sum_plots <- function(model_b, model_p, pred_df, train, test, mat_sex, raw_dat, period){
    
    if(period == "Jan/Feb"){
      pred_df <- filter(pred_df, period == "Jan/Feb")
      train <- filter(train, Period == "Jan/Feb")
      test <- filter(test, Period == "Jan/Feb")
    } else if(period == "Sep/Oct"){
      pred_df <- filter(pred_df, period == "Sep/Oct")
      train <- filter(train, Period == "Sep/Oct")
      test <- filter(test, Period == "Sep/Oct")
    } else if(period == "Apr/May"){
      pred_df <- filter(pred_df, period == "Apr/May")
      train <- filter(train, Period == "Apr/May")
      test <- filter(test, Period == "Apr/May")
    } else{
      pred_df <- pred_df
      train <- train
      test <- test
    }
    
    # Add A,B,C,D label
    if(mat_sex == "Legal male"){
      fig_lab = "A"
    } else if(mat_sex == "Immature male"){
      fig_lab = "B"
    } else if(mat_sex == "Mature female"){
      fig_lab = "C"
    } else{
      fig_lab = "D"
    }
    
    # Plot 1: Predicted vs. observed (test data) -----------------------------
    pred <- suppressWarnings(predict.gbm(model_p, # fitted model to predict
                        test %>% filter(extrap_log_rd>0), # data to predict to
                        n.trees=model_p$gbm.call$best.trees, # see help
                        type="response")) # predict probabilities
    
    obs <- test$extrap_log_rd[which(test$extrap_log_rd >0)]
    
    RMSE <- sqrt(mean((obs-pred)^2))
    
    d <- as.data.frame(cbind(obs, pred))
    
    suppressWarnings(cor.test(d$obs, d$pred, method = "spearman")) -> cor_out
    
    vals.ranked <- data.frame(cbind(rank(d$obs, ties.method = 'average'),
                                    rank(d$pred, ties.method = 'average')))
    
    colnames(vals.ranked) <- c('obs', 'preds')
    rho <- cov(vals.ranked) / (sd(vals.ranked$obs) * sd(vals.ranked$preds))
    
  
    # Calculate PDE
    PDE <- (model_p$self.statistics$null-model_p$self.statistics$resid)/model_p$self.statistics$null
    
    # Calculate AUC
    pred <- suppressWarnings(predict.gbm(model_b, # fitted model to predict
                        test, # data to predict to
                        n.trees=model_b$gbm.call$best.trees, # see help
                        type="response")) # predict probabilities
    
    obs <- test$PA 
    
    d <- cbind(obs, pred)
    pres <- d[d[,1]==1, 2]
    abs <- d[d[,1]==0, 2]
    e <- dismo::evaluate(p=pres, a=abs)
    AUC <- e@auc
    
    rocobj <- roc(test$PA, pred, AUC = TRUE)
    
    tt <- data.frame(FPR = 1-rocobj$specificities, TPR = rocobj$sensitivities)
      
    # Plot 2: Timeseries of bycatch by mat-sex ------------------------------
    # Calculate mean catch per haul per hour and CI
    ts_dat <- raw_dat %>%
      group_by(year) %>%
      reframe(N = n(),
              var = (var(total_extrap)/N),
              sd = sqrt(var),
              ci = 1.96 * sd,
              total = (sum(total_extrap)/1000),
              avg = mean(total_extrap),
              log_tot = log(total+1),
              log_avg = log(avg + 1),
              log_ci = log(ci+1))
    
    # Specify plotting breaks
    if(mat_sex == "Legal male"){
      breaks = c(2, 4, 6, 8, 10)
    } else if(mat_sex == "Immature male"){
      breaks = c(0.5, 1, 1.5, 2, 2.5)
    } else if(mat_sex == "Mature female"){
      breaks = c(2, 4, 6, 8, 10, 12)
    } else{
      breaks = c(0.5, 1, 1.5, 2, 2.5)
    }
    
    # Plot
    ggplot(ts_dat, aes(x = year, y = avg)) +
      geom_ribbon(ts_dat, mapping =aes(ymin = avg - ci, ymax = avg+ci), color = NA, fill = "#6051B0",
                  alpha = 0.25)+
      geom_line(color = "#6051B0", linewidth = 1) +
      theme_bw() +
      xlab("Year") +
      ylab("Catch per hour") +
      scale_x_continuous(breaks = seq(min(ts_dat$year), max(ts_dat$year), by = 2),
                         labels= seq(min(ts_dat$year), max(ts_dat$year), by = 2))+
      annotate(geom = "text", label = fig_lab, hjust = -0.2, vjust = 1.2,
               x = -Inf, y = Inf,
               fontface = "bold", color = "black", size = 8) +
      theme(axis.title = element_text(size = 20), axis.text.y = element_text(size = 15),
            axis.text.x=element_text(angle=40, hjust=1, size = 15))  -> bycatch_ts
  
    # Plot 3: Mean lat/lon for predicted vs. observed distribution -----------
    pred_df %>% # You could predict with the full dataset (training/testing), but better to use testing for out-of-sample ability
      mutate(grp = "Predicted") %>%
      dplyr::select(predict_year, period, x, y, log_count, grp) %>%
      dplyr::rename(year = predict_year) -> pred_df
    
    rbind(train, test) %>% # you want to include the full dataset for observed CODs, otherwise it may be biased by the training/testing split
      mutate(grp = "Observed") %>%
      dplyr::rename(log_count = extrap_log_rd, year = predict_year, period = Period) %>%
      dplyr::select(year,period, x, y, log_count, grp) -> obs_df
    
    rbind(pred_df, obs_df) %>%
      mutate(log_count = log_count) %>%
      #filter(log_count >0) %>%
      group_by(year, grp) %>%
      dplyr::summarise(mean_lat = mean(y),
                       mean_lon = mean(x),
                       sum_count = sum(log_count),
                       chk_lat = weighted.mean(y, log_count),
                       chk_lon = weighted.mean(x, log_count),
                       mean_weighted_lat = sum(y*log_count)/sum(log_count),
                       mean_weighted_lon = sum(x*log_count)/sum(log_count),
                       N = n()) %>%
      ungroup() %>%
      mutate(mean_weighted_lat = ifelse(is.na(mean_weighted_lat) == TRUE, mean_lat, mean_weighted_lat),
             mean_weighted_lon = ifelse(is.na(mean_weighted_lon) == TRUE, mean_lon, mean_weighted_lon)) %>%
      mutate(shrt = paste0("'", sprintf('%02d', year %% 100))) -> wtd_crds
    
    # Get region layers
    region_layers <- suppressWarnings(akgfmaps::get_base_layers(select.region = "bs.south", set.crs="auto"))
    
   
    # Calculate euclidean distance between obs and pred, set jitter to euc threshold
    obs.wtd <- wtd_crds %>% 
      filter(grp == "Observed") %>% 
      dplyr::select(year, mean_weighted_lat, mean_weighted_lon) %>%
      sf::st_as_sf(coords = c(x = "mean_weighted_lon", y = "mean_weighted_lat"), crs = sf::st_crs(map.crs)) %>%
      sf::st_transform(crs = "epsg:4326") %>%
      cbind(st_coordinates(.)) %>%
      as.data.frame()
    
    pred.wtd <- wtd_crds %>% 
      filter(grp == "Predicted") %>% 
      dplyr::select(year, mean_weighted_lat, mean_weighted_lon) %>%
      sf::st_as_sf(coords = c(x = "mean_weighted_lon", y = "mean_weighted_lat"), crs = sf::st_crs(map.crs)) %>%
      sf::st_transform(crs = "epsg:4326") %>%
      cbind(st_coordinates(.)) %>%
      as.data.frame()
    

    # Calculate Great Circle distance
    spDists(x= as.matrix(obs.wtd[,-c(1,4)]), y = as.matrix(pred.wtd[,-c(1,4)]), longlat = TRUE, diagonal = TRUE) -> sp_df
    
    data.frame(year = unique(wtd_crds$year), dist = as.numeric(sp_df),
               total = ts_dat[ts_dat$year!=2020,]$total) -> sp.dist
    
    # Specify plotting years that need to be jittered based on GC distance < 5
    jit.yr <- filter(sp.dist, dist <5) %>%
                pull(year)
    
    
    # Set breaks and panel extents
    breaks.x = c(-165, -164, -163, -162)
    labels.x = paste0(c(165, 164, 163, 162), "°W")
    
    breaks.y = c(58, 57, 56)
    labels.y = paste0(c(58, 57, 56), "°N")
    
    panel_extent <- data.frame(y = c(55.4, 58.4), #Bristol Bay # was 54
                               x = c(-165.2, -161.9)) %>%
      akgfmaps::transform_data_frame_crs(out.crs = coldpool:::ebs_proj_crs)
    
    
    # Plot
    ggplot() +
      geom_sf(data = st_as_sf(BB_strata), fill=NA, size=0.8, color="black")+
      scale_color_manual(values = c("#6051B0", "#999900"), name = element_blank())+
      ggplot2::geom_sf(data = region_layers$akland, 
                       fill = "grey70", 
                       color = "black") +
      ggplot2::coord_sf(xlim = panel_extent$x, 
                        ylim = panel_extent$y) + 
      geom_text(data = wtd_crds %>% filter(!(year %in% jit.yr)), mapping = aes(mean_weighted_lon, mean_weighted_lat, label = shrt, color = grp), 
                      fontface = "bold", size = 4)+
      geom_text_repel(data = wtd_crds %>% filter((year %in% jit.yr)), mapping = aes(mean_weighted_lon, mean_weighted_lat, label = shrt, color = grp), 
                fontface = "bold", size = 4, force = 0.25, segment.color = NA, max.overlaps = 40)+
      #geom_text(data = labs, aes(x=x, y=y, label= lab), fontface = "bold", size = c(5, 5) , color = c("#6051B0", "#999900"))+
      ggplot2::scale_x_continuous(name = "Longitude", 
                                  breaks = breaks.x, 
                                  labels = labels.x) + 
      ggplot2::scale_y_continuous(name = "Latitude", 
                                  breaks = breaks.y, labels = labels.y) +
      annotate(geom = "text", label = fig_lab, hjust = -0.2, vjust = 1.2,
               x = -Inf, y = Inf,
               fontface = "bold", color = "black", size = 8) +
      annotate(geom = "text", label = "Observed", hjust = -0.2, vjust = 36,
               x = -Inf, y = Inf,
               fontface = "bold", color = "#6051B0", size = 5) +
      annotate(geom = "text", label = "Predicted", hjust = -0.2, vjust = 37.5,
               x = -Inf, y = Inf,
               fontface = "bold", color = "#999900", size = 5) +
      #ggtitle(mat_sex)+
      theme_bw() +
      theme(legend.position = "none", legend.text = element_text(size = 8), axis.title = element_text(size = 24), 
            axis.text = element_text(size = 20), plot.title = element_text(size = 24), legend.background = element_blank(),  
            legend.key = element_rect(fill = NA), aspect.ratio = 1) -> wtd_crds_plot
    
   
    # Plot 3: Covariate importance -------------------------------------------
    # Get covariate summaries
    summary(model_b) -> sum_b
    summary(model_p) -> sum_p
    
    # Bind covariate summaries, calculate mean and SD relative influence by covariate
    rbind(sum_b, sum_p) %>%
      group_by(var) %>%
      dplyr::summarise(sum_inf = mean(rel.inf),
                       se = (sd(rel.inf))/sqrt(2)) -> inf_sum
    
    # Rank covariates by mean relative influence, select the top 4
    top <- inf_sum[order(inf_sum$sum_inf, decreasing = TRUE), ] %>%
      slice(1:4) %>%  # Order data descending and select top 4 vals
      mutate(var = case_when((var == "SAP_Count") ~ "BBRKC abundance",
                             (var == "YFS_dfish") ~ "Yellowfin sole fishery cpue",
                             (var == "RS_dfish") ~ "Rock sole fishery cpue",
                             (var == "YFS_TAC") ~ "Yellowfin sole quota",
                             (var == "RS_TAC") ~ "Rock sole quota",
                             (var == "Jul_Aug_SST") ~ "Jul/Aug SST",
                             (var == "May_Jun_SST") ~ "May/Jun SST",
                             (var == "Sep_Oct_SST") ~ "Sep/Oct SST",
                             (var == "RS_CPUE") ~ "Rock sole summer survey CPUE",
                             (var == "Depth") ~ "Depth"))
                             
    # Specify plotting parameters by sex-size/maturity                      
    if(mat_sex == "Legal male"){
      xlabs = c("Yellowfin sole\nfishery CPUE", "BBRKC summer\nsurvey CPUE", "Rock sole\nfishery CPUE",
                  "Rock sole summer\nsurvey CPUE")
      val1 = 0
      val2 = 2.5
      
    } else if(mat_sex == "Immature male"){
      xlabs = c("Yellowfin sole\nfishery CPUE", "Rock sole\nfishery CPUE", "BBRKC summer\nsurvey CPUE", 
                "Rock sole summer\nsurvey CPUE")
      val1 = 1.5
      val2 = 1.25
      
    } else if(mat_sex == "Mature female"){
      xlabs = c("Yellowfin sole\nfishery CPUE", "Rock sole\nfishery CPUE", "BBRKC summer\nsurvey CPUE",
                "Depth")
      val1 = 1.5
      val2 = 2
      
      
    } else{
      xlabs = c("Yellowfin sole\nfishery CPUE", "Rock sole\nfishery CPUE", "BBRKC summer\nsurvey CPUE",
                "Sep/Oct SST")
      val1 = 0
      val2 = 3
    }
    
    # Plot
    ggplot(top, aes(reorder(var, -sum_inf), sum_inf)) +
      geom_bar(stat = "identity", fill = c("#312271", "#48398F", "#6051B0", "#7A6FBA")) +
      theme_bw() +
      xlab(element_blank()) +
      geom_errorbar(aes(ymin=sum_inf - se, ymax=sum_inf+se),size=0.8, width=0.3)+
      ylab("Relative influence (%)") +
      scale_x_discrete(labels = function(x) stringr::str_wrap(xlabs, width = 16))+
      annotate(geom = "text", label = fig_lab, hjust = -0.2, vjust = 1.2,
               x = -Inf, y = Inf,
               fontface = "bold", color = "black", size = 8) +
      theme(axis.title = element_text(size = 22), axis.text.y = element_text(size = 20),
            axis.text.x=element_text(size = 18)) -> var_inf_plot
            #axis.text.x=element_text(angle=40, hjust=1, size = 13)) -> var_inf_plot
    
    
    # Plot 4: GC distance between obs and predicted bycatch CODs -------------------
    gc.obs <- vector()
    
    for(ii in 1:(nrow(obs.wtd)-1)){
      spDists(x = as.matrix(obs.wtd[ii, -c(1,4)]), y = as.matrix(obs.wtd[ii+1, -c(1, 4)]), 
              longlat = TRUE, diagonal = TRUE) -> gc.obs[ii]
    }
    
    data.frame(year = obs.wtd$year[-1], gc.obs.dist = gc.obs) -> gc.obs.dist

    sp.dist %>%
      reframe(mean.COD.dist = mean(dist),
              se.COD.dist = (sd(dist))/(sqrt(nrow(sp.dist))),
              max.COD.dist = max(dist)) -> gc.COD.dist
    
    # Plot 5: Response curves -----------------------------------------------------------
      # Extract variable importance information for top 6 model_b vars
      summary(model_b) %>%
        slice(1:4) -> top_b
      
      top_b %>%
        dplyr::select(var) %>%
        pull() %>%
        map_df(~plot(model_b, .x, return.grid = TRUE)) %>%
        reshape::melt(id.vars = "y", variable_name = "var") %>%
        filter(is.na(value) == "FALSE") %>% 
        right_join(., top_b) %>%
        dplyr::rename(Fitted = y, Value = value) %>%
        mutate(var = case_when((var == "SAP_Count") ~ "BBRKC summer survey CPUE",
                               (var == "YFS_dfish") ~ "Yellowfin sole fishery CPUE",
                               (var == "RS_dfish") ~ "Rock sole fishery CPUE",
                               (var == "YFS_TAC") ~ "Yellowfin sole quota",
                               (var == "RS_TAC") ~ "Rock sole quota",
                               (var == "Jul_Aug_SST") ~ "Jul/Aug SST",
                               (var == "May_Jun_SST") ~ "May/Jun SST",
                               (var == "RS_CPUE") ~ "Rock sole summer survey CPUE",
                               (var == "Jan_Feb_Ice") ~ "Jan/Feb ice",
                               (var == "Sed") ~ "Sediment",
                               (var == "May_Aug_BT") ~ "May-Aug bottom temperature",
                               (var == "Depth") ~ "Depth",
                               (var == "Sep_Oct_SST") ~ "Sep/Oct SST")) -> top_b_df
      
      # Specify plotting labels
      unique(top_b_df$var) -> top_b_names
      
      data.frame(names = top_b_names, rel.inf = top_b$rel.inf) %>%
        mutate(lab = paste0(names, " (", round(rel.inf, 1), "%)")) -> b_labs
   
      # Extract variable importance information for top 6 model_p vars
      summary(model_p) %>%
        slice(1:4) -> top_p
      
      top_p %>%
        dplyr::select(var) %>%
        pull() %>%
        map_df(~plot(model_p, .x, return.grid = TRUE)) %>%
        reshape::melt(id.vars = "y", variable_name = "var") %>%
        filter(is.na(value) == "FALSE") %>% 
        right_join(., top_p) %>%
        dplyr::rename(Fitted = y, Value = value) %>%
        mutate(var = case_when((var == "SAP_Count") ~ "BBRKC summer survey CPUE",
                               (var == "YFS_dfish") ~ "Yellowfin sole fishery CPUE",
                               (var == "RS_dfish") ~ "Rock sole fishery CPUE",
                               (var == "YFS_TAC") ~ "Yellowfin sole quota",
                               (var == "RS_TAC") ~ "Rock sole quota",
                               (var == "Jul_Aug_SST") ~ "Jul/Aug SST",
                               (var == "May_Jun_SST") ~ "May/Jun SST",
                               (var == "RS_CPUE") ~ "Rock sole summer survey CPUE",
                               (var == "YFS_CPUE") ~ "Yellowfin sole summer survey CPUE",
                               (var == "Jan_Feb_Ice") ~ "Jan/Feb ice",
                               (var == "Sed") ~ "Sediment",
                               (var == "May_Aug_BT") ~ "May-Aug bottom temperature",
                               (var == "Depth") ~ "Depth",
                               (var == "Sep_Oct_SST") ~ "Sep/Oct SST"))  -> top_p_df
      
      # Specify plotting labels
      unique(top_p_df$var) -> top_p_names
      
      data.frame(names = top_p_names, rel.inf = top_p$rel.inf) %>%
        mutate(lab = paste0(names, " (", round(rel.inf, 1), "%)")) -> p_labs
      
      # Plot
      ggplot(top_b_df, aes(x = Value, y = Fitted)) +
        geom_line(linewidth = 1, color = "#6051B0") +
        facet_wrap(~factor(var, levels = top_b_names, labels = b_labs$lab), scales = "free_x", ncol = 4, nrow = 1)+
        ggtitle(fig_lab)+
        theme_bw() +
        theme(panel.spacing.x = unit(4, "mm"), strip.text = element_text(size = 8)) -> b_response
    
      
      ggplot(top_p_df, aes(x = Value, y = Fitted)) +
        geom_line(linewidth = 1, color = "#999900") +
        facet_wrap(~factor(var, levels = top_p_names, labels = p_labs$lab), scales = "free_x", ncol = 4, nrow = 1)+
        theme_bw() +
        theme(panel.spacing.x = unit(4, "mm"), strip.text = element_text(size = 8)) -> p_response
   
      return(list(wtd_crds_df = wtd_crds, wtd_crds_plot = wtd_crds_plot, sp.dist = sp.dist, 
                  gc.obs.dist = gc.obs.dist,gc.COD.dist = gc.COD.dist, var_inf_plot = var_inf_plot, 
                  bycatch_ts = bycatch_ts, top = top,AUC= AUC, 
                  rho = cor_out$estimate, PDE = PDE))   
    }

  # Generate dot plot of sampling distribution
    # @param resp_data: catch data (options = "catch_lm", "catch_im", "catch_mf", "catch_imf")
    # @param predict_yr: years over which to map sampling distribution
    # @return Returns sampling distribution plot and summary table
    
    MakeDotPlot <- function(resp_data, predict_yr){
    
    presence_df <- data.frame()
    absence_df <- data.frame()
    highdensity_df <- data.frame()
    
    
    for (ii in 1:length(predict_yr)){
      # Filter response data by season and plotting years
      resp_data %>%
        filter(year == predict_yr[ii]) -> resp_data2
      
      # Specify high density quantile
      hd <- stats::quantile(resp_data2 %>% filter(total_extrap >0) %>% dplyr::pull(total_extrap), .9)
      
      # Specify high density, presence, and absence dfs
      presence = resp_data2[resp_data2[, "total_extrap"] > 0, ]
      absence = resp_data2[resp_data2[, "total_extrap"] == 0, ]
      highdensity = resp_data2[resp_data2[, "total_extrap"] >= hd, ]
      
      
      rbind(presence_df, presence) -> presence_df
      rbind(absence_df, absence) -> absence_df
      rbind(highdensity_df, highdensity) -> highdensity_df
      
    }
    
    # Transform data frame crs', make into sf objects
    presence_df %>%
      sf::st_as_sf(coords = c("lon", "lat"), crs = in.crs) %>%
      sf::st_transform(sf::st_crs(map.crs)) %>%
      vect(.) %>%
      mask(., BB_strata)%>%
      sf::st_as_sf() -> pres_df
    
    absence_df %>%
      sf::st_as_sf(coords = c("lon", "lat"), crs = in.crs) %>%
      sf::st_transform(sf::st_crs(map.crs)) %>%
      vect() %>%
      mask(., BB_strata)%>%
      sf::st_as_sf()-> abs_df
    
    highdensity_df %>%
      sf::st_as_sf(coords = c("lon", "lat"), crs = in.crs) %>%
      sf::st_transform(sf::st_crs(map.crs)) %>%
      vect() %>%
      mask(., BB_strata)%>%
      sf::st_as_sf()-> hd_df
    
    # Set up plot boundary
    plot.boundary.untrans <- data.frame(y = c(54.25, 59.25), 
                                        x = c(-167.5, -158)) # plot boundary unprojected
    
    plot.boundary <- plot.boundary.untrans %>%
      sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
      sf::st_transform(crs = map.crs) %>%
      sf::st_coordinates() %>%
      as.data.frame() %>%
      dplyr::rename(x = X, y = Y) # plot boundary projected
    
    # Specify dot features
    abs.name = "absent"
    pres.name = "present"
    hd.name = "top 10%"
    abs.col = "#FDE333"
    pres.col = "#009B95"
    hd.col = "#4B0055" 
    abs.shape = 16
    pres.shape = 1
    hd.shape = 16
    abs.size = 2
    pres.size = 2
    hd.size = 2
    pres.fac <- 2
    abs.fac <- 1
    hd.fac <- 3
    
    # Now go through and set up the dot locations and add them to legend
    if (is.data.frame(abs_df)) {
      leg.name <- abs.name
      leg.col <- abs.col
      leg.shape <- abs.shape
      leg.size <- abs.size
      abs.fac <- 1
    } else {
      abs.fac <- 0
    }
    
    if (is.data.frame(pres_df)) {
      leg.name <- c(leg.name, pres.name)
      leg.col <- c(leg.col, pres.col)
      leg.shape <- c(leg.shape, pres.shape)
      leg.size <- c(leg.size, pres.size)
      pres.fac <- abs.fac + 1
    }
    
    if (is.data.frame(hd_df)) {
      
      leg.name <- c(leg.name, hd.name)
      leg.col <- c(leg.col, hd.col)
      leg.shape <- c(leg.shape, hd.shape)
      leg.size <- c(leg.size, hd.size)
      hd.fac <- pres.fac + 1
    }
    
    rbind(abs_df %>% mutate(type = "absent"),
          pres_df %>% mutate(type = "present"),
          hd_df %>% mutate(type = "hd")) -> PA_df
    
    # Map
    ggplot2::ggplot() +
      ggplot2::geom_sf(data = PA_df %>% filter(type == "absent"), 
                       alpha = .15, size = abs.size, shape = abs.shape, ggplot2::aes(color = factor(abs.fac)))+
      
      ggplot2::geom_sf(data = PA_df %>% filter(type == "present"), 
                       size = pres.size, ggplot2::aes(color = factor(pres.fac)), shape = pres.shape, stroke = .8)+
      ggplot2::geom_sf(data = PA_df %>% filter(type == "hd"), size = hd.size, 
                       shape = hd.shape, ggplot2::aes(color = factor(hd.fac)))+
      ggplot2::geom_sf(data = st_as_sf(BB_strata),
                       fill = NA,
                       color = "black",
                       linewidth = 2)+
      ggplot2::geom_sf(data = region_layers$akland, 
                       fill = "grey70", 
                       color = "black")+
      labs(title = "Bycatch Sampling Distribution")+
      coord_sf(xlim = plot.boundary$x,
               ylim = plot.boundary$y)+
      ggplot2::scale_color_manual(name = NULL, values = leg.col, labels = leg.name) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(shape = leg.shape, size = leg.size)))+
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(color = "black", fill = NA),
        panel.background = ggplot2::element_rect(fill = NA, color = "black"),
        legend.key = ggplot2::element_rect(fill = NA, color = "grey30"),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        plot.title = element_text(size = 24),
        plot.subtitle = element_text(size = 14),
        axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 20),
        legend.text = ggplot2::element_text(size = 20), legend.title = ggplot2::element_text(size = 12),
        plot.background = ggplot2::element_rect(fill = "white", color = "white")) -> sum_dotplot
  
    
    # Summary table
    right_join(
      PA_df %>%
        group_by(year, type) %>%
        reframe(n = n()),
      PA_df %>%
        group_by(year) %>%
        reframe(Total = n()), by = "year") -> sum_table
    
    return(list(sum_dotplot = sum_dotplot, sum_table = sum_table))
  }

# RUN FUNCTIONS ---------------------------------------------------------------------------------------------------------------------
  predict_yr <- c(1997:2019, 2021:2023) 

 # LEGAL MALES 
  model_b <- lm_modelb
  model_p <- lm_modelp
  test <- lm_test
  train <- lm_train
  mat_sex <- "Legal male"
  raw_dat <- catch_lm
  period <- "All"
  
  # Generate prediction df from models using test data
  pred_df_test(lm_modelb, lm_modelp, lm_test) -> pred_df

  # Generate summary plots
  sum_plots(lm_modelb, lm_modelp, pred_df, lm_train, lm_test, "Legal male", catch_lm, "All") -> lm_plot_out

 # IMMATURE MALES 
  model_b <- im_modelb
  model_p <- im_modelp
  test <- im_test
  train <- im_train
  mat_sex <- "Immature male"
  raw_dat <- catch_im
  period <- "All"
  
  # Generate prediction df from models using test data
  pred_df_test(im_modelb, im_modelp, im_test) -> pred_df
  
  # Generate summary plots
  sum_plots(im_modelb, im_modelp, pred_df, im_train, im_test, "Immature male", catch_im, "All") -> im_plot_out

 # MATURE FEMALES
  model_b <- mf_modelb
  model_p <- mf_modelp
  test <- mf_test
  train <- mf_train
  mat_sex <- "Mature female"
  raw_dat <- catch_mf
  period <- "All"
  
  # Generate prediction df from models using test data
  pred_df_test(mf_modelb, mf_modelp, mf_test) -> pred_df
    
  # Generate summary plots
  sum_plots(mf_modelb, mf_modelp, pred_df, mf_train, mf_test, "Mature female", catch_mf, "All") -> mf_plot_out
    
 # IMMATURE FEMALES
  model_b <- imf_modelb
  model_p <- imf_modelp
  test <- imf_test
  train <- imf_train
  mat_sex <- "Immature female"
  raw_dat <- catch_imf
  period <- "All"
  
  # Generate prediction df from models using test data
  pred_df_test(imf_modelb, imf_modelp, imf_test) -> pred_df
    
  # Generate summary plots
  sum_plots(imf_modelb, imf_modelp, pred_df, imf_train, 
            imf_test, "Immature female", catch_imf, "All") -> imf_plot_out

# WEIGHTED LAT TIMESERIES WITH SURVEY ------------------------------------------------
    # Bind bycatch coordinates from sum_plots above
    rbind(lm_plot_out$wtd_crds_df %>% mutate(MAT_SEX = "Legal male"),
          im_plot_out$wtd_crds_df %>% mutate(MAT_SEX = "Immature male"),
          mf_plot_out$wtd_crds_df %>% mutate(MAT_SEX = "Mature female"),
          imf_plot_out$wtd_crds_df %>% mutate(MAT_SEX = "Immature female")) %>%
      sf::st_as_sf(coords = c(x = "mean_weighted_lon", y = "mean_weighted_lat"), crs = sf::st_crs(map.crs)) %>%
      sf::st_transform(crs = "epsg:4326") %>%
      cbind(st_coordinates(.)) %>%
      as.data.frame() -> wtd_crds_df
  
    # Read in survey CODs
    read.csv("./Data/BBRKC.SAP.wtdcrds.csv") %>%
      rename(Y = mean_weighted_lat, year = AKFIN_SURVEY_YEAR) %>%
      mutate(grp = "Survey",
             MAT_SEX = case_when((MAT_SEX == "Immature Female")~"Immature female",
                                 (MAT_SEX == "Mature Female") ~ "Mature female",
                                 (MAT_SEX == "Immature Male") ~ "Immature male",
                                 (MAT_SEX == "Legal Male") ~ "Legal male")) %>%
      filter(year >1996) -> surv.wtdcrds
  
    # Bind survey and bycatch CODs
    rbind(surv.wtdcrds %>% dplyr::select(year, MAT_SEX, Y, grp),
          wtd_crds_df %>% dplyr::select(year, MAT_SEX, Y, grp)) %>%
      rbind(data.frame(year = rep(2020, 12), MAT_SEX = rep(c("Immature female", "Immature male",
                                                             "Legal male", "Mature female"),3),
                       Y = rep(NA, 12), grp = rep(c("Survey", "Observed", "Predicted"), 4))) -> wtd.df
  
    # Plot ts
    ggplot() +
      geom_line(wtd.df, mapping = aes(x = year, y = Y, color = grp, linetype = grp), linewidth = 1.25)+
      theme_bw()+
      scale_color_manual(values = c("#6051B0", "#6051B0", "#999900"),
                         labels = c("Observed bycatch", "Predicted bycatch", "NMFS survey"),
                         name = "none")+
      scale_linetype_manual(values = c("solid", "twodash", "solid"),
                            labels = c("Observed bycatch", "Predicted bycatch", "NMFS survey"),
                            name = "none")+
      scale_x_continuous(breaks = seq(min(wtd_crds_df$year), max(wtd_crds_df$year), by = 2),
                         labels= seq(min(wtd_crds_df$year), max(wtd_crds_df$year), by = 2))+
      ylab("Center of distribution (latitude)")+
      xlab("Year")+
      ylim(c(55.5, 58.8))+
      facet_wrap(~ factor(MAT_SEX, levels = c("Legal male", "Immature male",
                                              "Mature female", "Immature female")))+
      guides(color= guide_legend(nrow = 1), linetype = guide_legend(nrow = 3))+
      #ggtitle(paste(mat_sex, "bycatch"))+
      #scale_y_continuous(breaks = breaks)+
      theme(axis.title = element_text(size = 17), axis.text.y = element_text(size = 17),
            axis.text.x=element_text(angle=40, hjust=1, size = 17),
            legend.direction = "horizontal",
            strip.text = element_text(size = 17),
            legend.position = "top",
            legend.background = element_rect(fill = NA, color = NA),
            legend.key = element_rect(fill = NA),
            legend.title = element_blank(),
            legend.text = element_text(size = 12)) -> wtd.ts
    
# CORRELATIONS BETWEEN OBSERVED-PREDICTED, OBSERVED-SURVEY, CORRECTING FOR AUTOCORRELATIONS VIA PYPER AND PETERMAN (1998) ----------------------------
    # Filter crds by mat sex
    wtd.df %>%
      dplyr::filter(MAT_SEX == "Legal male") %>%
      na.omit() -> lm.crds
    
    wtd.df %>%
      dplyr::filter(MAT_SEX == "Immature male") %>%
      na.omit() -> im.crds
    
    wtd.df %>%
      dplyr::filter(MAT_SEX == "Mature female") %>%
      na.omit() -> mf.crds
    
    wtd.df %>%
      dplyr::filter(MAT_SEX == "Immature female") %>%
      na.omit() -> imf.crds
  
    #Pyper and Peterman code written by Franz Mueter, last updated 2000
    # Required function for cor.test.PP:
    N.effective <- function(mat) {
      # written by Franz Mueter   Last modified: 23 October 2000
      # function to compute effective sample size for pairwise correlations among 
      # autocorrelated variables. 
      # Based on Pyper & Peterman (1998). CJFAS 55:2127-2140.  Eq. (1)
      # where summation was done over j = 1, ..., N/5
      # 
      # mat a matrix of variables, one variable per column
      #
      # function to compute simple estimates of autocorrelation up to lag N/5: 
      # (Eq. 7 in Pyper & Peterman)
      ar.fun <- function(x, max.lag = ceiling(sum(!is.na(x))/5)) {
        res <- rep(NA, max.lag)
        n <- length(x)
        for(i in 1.:max.lag) {
          x.bar <- mean(x, na.rm = T)
          res[i] <- ((n/(n - i)) * sum((x[1:(n - i)] - x.bar) * (x[(i + 1):n] - x.bar), na.rm
                                       = T))/sum((x - x.bar)^2, na.rm = T)
        }
        res
      }
      AR <- apply(mat, 2., ar.fun)
      k <- ncol(mat)
      if(is.matrix(AR)) {
        AR1 <- vector("list", k)
        for(i in 1:k) AR1[[i]] <- AR[, i]
        AR <- AR1  
      }
      N <- t(!is.na(mat)) %*% (!is.na(mat))
      N.lags <- ceiling(N/5.)
      N.eff <- matrix(0., k, k)
      # constrain effective N to smaller than or equal to actual N:
      for(i in 1.:k) {
        for(j in 1.:i) {
          lags <- 1.:N.lags[i, j]
          Nij <- N[i, j]
          N.eff[i, j] <- round((1./Nij + 2./Nij * sum((Nij - lags)/Nij * AR[[i]][lags] * AR[[
            j]][lags]))^(-1.))
        }
      }
      j <- N.eff > N
      N.eff[j] <- N[j]
      N.eff + t(N.eff) - diag(diag(N.eff))
    }
    
    # Function to test for significant correlation between two time series
    # in the presence of autocorrelation in one or both of the series:
    cor.test.PP <- function(x, y) {
      # Function to test for significant correlations between x and y, which may be autocorrelated, 
      # using modified Chelton method after Pyper & Peterman (1998)
      # Eqn. 3 with N*-2 degrees of freedom
      N.eff <- N.effective(cbind(x, y))[1, 2]
      r <- cor(x, y, use="pair")
      fun <- function(alpha, N, r) {
        t2 <- qt(1 - alpha/2, N - 2)^2
        sqrt(t2/(t2 + N - 2)) - abs(r)
      }
      p.value <- uniroot(fun, c(1e-015, 0.9999), N = N.eff, r = r)$root
      cat("Two-sided test\n\n")
      c(correlation = r, P.value = p.value, N.eff= N.eff)
    }
    
    # Pivot dfs wider
    pivot_wider(lm.crds, values_from = Y, names_from = grp) -> lm.crds2
    pivot_wider(im.crds, values_from = Y, names_from = grp) -> im.crds2
    pivot_wider(mf.crds, values_from = Y, names_from = grp) -> mf.crds2
    pivot_wider(imf.crds, values_from = Y, names_from = grp) -> imf.crds2
  
    # Run functions
      # Observed vs Predicted
      cor.test.PP(lm.crds2$Observed, lm.crds2$Predicted) -> lm.corr.obs_pred
      cor.test.PP(im.crds2$Observed, im.crds2$Predicted) -> im.corr.obs_pred
      cor.test.PP(mf.crds2$Observed, mf.crds2$Predicted) -> mf.corr.obs_pred
      cor.test.PP(imf.crds2$Observed, imf.crds2$Predicted) -> imf.corr.obs_pred
      
      # Observed vs. Survey
      cor.test.PP(lm.crds2$Observed, lm.crds2$Survey) -> lm.corr.obs_surv
      cor.test.PP(im.crds2$Observed, im.crds2$Survey) -> im.corr.obs_surv
      cor.test.PP(mf.crds2$Observed, mf.crds2$Survey) -> mf.corr.obs_surv
      cor.test.PP(imf.crds2$Observed, imf.crds2$Survey) -> imf.corr.obs_surv
  
# SAMPLE SIZE BY PERIOD/MATSEX TABLE --------------------------------------------
  rbind(rbind(lm_train, lm_test) %>%
          mutate(MAT_SEX = "Legal male"),
        rbind(im_train, im_test) %>%
          mutate(MAT_SEX = "Immature male"),
        rbind(mf_train, mf_test) %>%
          mutate(MAT_SEX = "Mature female"),
        rbind(imf_train, imf_test) %>%
          mutate(MAT_SEX = "Immature female")) %>%
        group_by(Period, MAT_SEX) %>%
        reframe(N_pres = sum(PA == 1),
                N_abs = sum(PA == 0)) -> sum_table

# SAMPLING DISTRIBUTION ----------------------------------------------------------
 # Create empty dataframes to store presences, absences, and high density points
  presence_df <- data.frame()
  absence_df <- data.frame()
  highdensity_df <- data.frame()
  
 # Bind catch data across sex-size/maturity categories, filter by prediction periods
  resp_data <- rbind(catch_lm %>% dplyr::select(year, month, lat, lon, total_extrap), 
                     catch_im %>% dplyr::select(year, month, lat, lon, total_extrap), 
                     catch_mf %>% dplyr::select(year, month, lat, lon, total_extrap), 
                     catch_imf %>% dplyr::select(year, month, lat, lon, total_extrap)) %>%
    filter(month %in% c(1:2, 4:5, 9:10)) %>%
    group_by(year, month, lat, lon) %>%
    reframe(total_extrap = sum(total_extrap))

 # Loop through prediction years, fill presence, absence, and high density dfs
  for (ii in 1:length(predict_yr)){
    # Filter response data by season and plotting years
    resp_data %>%
      filter(year == predict_yr[ii]) -> resp_data2
    
    # Specify high density quantile
    hd <- stats::quantile(resp_data2 %>% filter(total_extrap >0) %>% dplyr::pull(total_extrap), .9)
    
    # Specify high density, presence, and absence dfs
    presence = resp_data2[resp_data2[, "total_extrap"] > 0, ]
    absence = resp_data2[resp_data2[, "total_extrap"] == 0, ]
    highdensity = resp_data2[resp_data2[, "total_extrap"] >= hd, ]
    
    
    rbind(presence_df, presence) -> presence_df
    rbind(absence_df, absence) -> absence_df
    rbind(highdensity_df, highdensity) -> highdensity_df
    
  }
  
 # Transform each data frame crs', make into sf objects, mask to Bristol Bay strata
  presence_df %>%
    sf::st_as_sf(coords = c("lon", "lat"), crs = in.crs) %>%
    sf::st_transform(sf::st_crs(map.crs)) %>%
    vect(.) %>%
    mask(., BB_strata)%>%
    sf::st_as_sf() -> pres_df
  
  absence_df %>%
    sf::st_as_sf(coords = c("lon", "lat"), crs = in.crs) %>%
    sf::st_transform(sf::st_crs(map.crs)) %>%
    vect() %>%
    mask(., BB_strata)%>%
    sf::st_as_sf()-> abs_df
  
  highdensity_df %>%
    sf::st_as_sf(coords = c("lon", "lat"), crs = in.crs) %>%
    sf::st_transform(sf::st_crs(map.crs)) %>%
    vect() %>%
    mask(., BB_strata)%>%
    sf::st_as_sf()-> hd_df
  
 # Set up plot boundary
  plot.boundary.untrans <- data.frame(y = c(54.25, 59.25), 
                                      x = c(-167.5, -158)) # plot boundary unprojected
  
  plot.boundary <- plot.boundary.untrans %>%
    sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
    sf::st_transform(crs = map.crs) %>%
    sf::st_coordinates() %>%
    as.data.frame() %>%
    dplyr::rename(x = X, y = Y) # plot boundary projected
  
 # Specify dot features
    abs.name = "absent"
    pres.name = "present"
    hd.name = "top 10%"
    abs.col = "#FDE333"
    pres.col = "#009B95"
    hd.col = "#4B0055" 
    abs.shape = 16
    pres.shape = 1
    hd.shape = 16
    abs.size = 2
    pres.size = 2
    hd.size = 2
    pres.fac <- 2
    abs.fac <- 1
    hd.fac <- 3
  
 # Set up the dot locations and add them to legend
  if (is.data.frame(abs_df)) {
    leg.name <- abs.name
    leg.col <- abs.col
    leg.shape <- abs.shape
    leg.size <- abs.size
    abs.fac <- 1
  } else {
    abs.fac <- 0
  }
  
  if (is.data.frame(pres_df)) {
    leg.name <- c(leg.name, pres.name)
    leg.col <- c(leg.col, pres.col)
    leg.shape <- c(leg.shape, pres.shape)
    leg.size <- c(leg.size, pres.size)
    pres.fac <- abs.fac + 1
  }
  
  if (is.data.frame(hd_df)) {
    
    leg.name <- c(leg.name, hd.name)
    leg.col <- c(leg.col, hd.col)
    leg.shape <- c(leg.shape, hd.shape)
    leg.size <- c(leg.size, hd.size)
    hd.fac <- pres.fac + 1
  }
  
  rbind(abs_df %>% mutate(type = "absent"),
        pres_df %>% mutate(type = "present"),
        hd_df %>% mutate(type = "hd")) -> PA_df
  
 # Map
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = PA_df %>% filter(type == "absent"), 
                     alpha = .15, size = abs.size, shape = abs.shape, ggplot2::aes(color = factor(abs.fac)))+
    
    ggplot2::geom_sf(data = PA_df %>% filter(type == "present"), 
                     size = pres.size, ggplot2::aes(color = factor(pres.fac)), shape = pres.shape, stroke = .8)+
    ggplot2::geom_sf(data = PA_df %>% filter(type == "hd"), size = hd.size, 
                     shape = hd.shape, ggplot2::aes(color = factor(hd.fac)))+
        # ggplot2::geom_sf(data = st_as_sf(area512), 
    #                  fill = NA, 
    #                  color = "purple",
    #                  linewidth = 1.5)+
    ggplot2::geom_sf(data = st_as_sf(BB_strata),
                     fill = NA,
                     color = "black",
                     linewidth = 2)+
    ggplot2::geom_sf(data = region_layers$akland, 
                     fill = "grey70", 
                     color = "black")+
    # ggplot2::geom_sf(data = st_as_sf(RKCSA),
    #                  fill = NA,
    #                  color = "red",
    #                  linewidth = 1.5)+
    labs(title = "Bycatch Sampling Distribution")+
    # ggplot2::geom_sf(data = st_as_sf(RKCSA_sub),
    #                  fill = NA,
    #                  color = "red",
    #                  linewidth = 1.5)+
    coord_sf(xlim = plot.boundary$x,
             ylim = plot.boundary$y)+
    ggplot2::scale_color_manual(name = NULL, values = leg.col, labels = leg.name) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(shape = leg.shape, size = leg.size)))+
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA),
      panel.background = ggplot2::element_rect(fill = NA, color = "black"),
      legend.key = ggplot2::element_rect(fill = NA, color = "grey30"),
      legend.position = "bottom",
      panel.grid.major = element_blank(),
      plot.title = element_text(size = 24),
      plot.subtitle = element_text(size = 14),
      axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 20),
      legend.text = ggplot2::element_text(size = 20), legend.title = ggplot2::element_text(size = 12),
      plot.background = ggplot2::element_rect(fill = "white", color = "white")) -> sum_dotplot
 

# SAMPLING DISTRIBUTION FOR MATURE FEMALES --------------------------------------
  presence_df <- data.frame()
  absence_df <- data.frame()
  highdensity_df <- data.frame()
  
  resp_data <- rbind(catch_mf %>% dplyr::select(year, month, lat, lon, total_extrap)) %>% 
    #filter(month %in% c(1:2, 4:5, 9:10)) %>%
    group_by(year, month, lat, lon) %>%
    reframe(total_extrap = sum(total_extrap))
  
  for (ii in 1:length(predict_yr)){
    # Filter response data by season and plotting years
    resp_data %>%
      filter(year == predict_yr[ii]) -> resp_data2
    
    # Specify high density quantile
    hd <- stats::quantile(resp_data2 %>% filter(total_extrap >0) %>% dplyr::pull(total_extrap), .9)
    
    # Specify high density, presence, and absence dfs
    presence = resp_data2[resp_data2[, "total_extrap"] > 0, ]
    absence = resp_data2[resp_data2[, "total_extrap"] == 0, ]
    highdensity = resp_data2[resp_data2[, "total_extrap"] >= hd, ]
    
    
    rbind(presence_df, presence) -> presence_df
    rbind(absence_df, absence) -> absence_df
    rbind(highdensity_df, highdensity) -> highdensity_df
    
  }
  
  # Transform data frame crs', make into sf objects
  presence_df %>%
    sf::st_as_sf(coords = c("lon", "lat"), crs = in.crs) %>%
    sf::st_transform(sf::st_crs(map.crs)) %>%
    vect(.) %>%
    mask(., BB_strata)%>%
    sf::st_as_sf() -> pres_df
  
  absence_df %>%
    sf::st_as_sf(coords = c("lon", "lat"), crs = in.crs) %>%
    sf::st_transform(sf::st_crs(map.crs)) %>%
    vect() %>%
    mask(., BB_strata)%>%
    sf::st_as_sf()-> abs_df
  
  highdensity_df %>%
    sf::st_as_sf(coords = c("lon", "lat"), crs = in.crs) %>%
    sf::st_transform(sf::st_crs(map.crs)) %>%
    vect() %>%
    mask(., BB_strata)%>%
    sf::st_as_sf()-> hd_df
  
  # Set up plot boundary
  plot.boundary.untrans <- data.frame(y = c(54.25, 59.25), 
                                      x = c(-167.5, -158)) # plot boundary unprojected
  
  plot.boundary <- plot.boundary.untrans %>%
    sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
    sf::st_transform(crs = map.crs) %>%
    sf::st_coordinates() %>%
    as.data.frame() %>%
    dplyr::rename(x = X, y = Y) # plot boundary projected
  
  # Specify dot features
  abs.name = "absent"
  pres.name = "present"
  hd.name = "top 10%"
  abs.col = "#FDE333"
  pres.col = "#009B95"
  hd.col = "#4B0055" 
  abs.shape = 16
  pres.shape = 1
  hd.shape = 16
  abs.size = 2
  pres.size = 2
  hd.size = 2
  pres.fac <- 2
  abs.fac <- 1
  hd.fac <- 3
  
  # Now go through and set up the dot locations and add them to legend
  if (is.data.frame(abs_df)) {
    leg.name <- abs.name
    leg.col <- abs.col
    leg.shape <- abs.shape
    leg.size <- abs.size
    abs.fac <- 1
  } else {
    abs.fac <- 0
  }
  
  if (is.data.frame(pres_df)) {
    leg.name <- c(leg.name, pres.name)
    leg.col <- c(leg.col, pres.col)
    leg.shape <- c(leg.shape, pres.shape)
    leg.size <- c(leg.size, pres.size)
    pres.fac <- abs.fac + 1
  }
  
  if (is.data.frame(hd_df)) {
    
    leg.name <- c(leg.name, hd.name)
    leg.col <- c(leg.col, hd.col)
    leg.shape <- c(leg.shape, hd.shape)
    leg.size <- c(leg.size, hd.size)
    hd.fac <- pres.fac + 1
  }
  
  rbind(abs_df %>% mutate(type = "absent"),
        pres_df %>% mutate(type = "present"),
        hd_df %>% mutate(type = "hd")) -> PA_df
  
  # Map
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = PA_df %>% filter(type == "absent"), 
                     alpha = .15, size = abs.size, shape = abs.shape, ggplot2::aes(color = factor(abs.fac)))+
    
    ggplot2::geom_sf(data = PA_df %>% filter(type == "present"), 
                     size = pres.size, ggplot2::aes(color = factor(pres.fac)), shape = pres.shape, stroke = .8)+
    ggplot2::geom_sf(data = PA_df %>% filter(type == "hd"), size = hd.size, 
                     shape = hd.shape, ggplot2::aes(color = factor(hd.fac)))+
    # ggplot2::geom_sf(data = st_as_sf(area512), 
    #                  fill = NA, 
    #                  color = "purple",
    #                  linewidth = 1.5)+
    ggplot2::geom_sf(data = st_as_sf(BB_strata),
                     fill = NA,
                     color = "black",
                     linewidth = 2)+
    ggplot2::geom_sf(data = region_layers$akland, 
                     fill = "grey70", 
                     color = "black")+
    # ggplot2::geom_sf(data = st_as_sf(RKCSA),
    #                  fill = NA,
    #                  color = "red",
    #                  linewidth = 1.5)+
    labs(title = "Bycatch Sampling Distribution")+
    # ggplot2::geom_sf(data = st_as_sf(RKCSA_sub),
    #                  fill = NA,
    #                  color = "red",
    #                  linewidth = 1.5)+
    coord_sf(xlim = plot.boundary$x,
             ylim = plot.boundary$y)+
    ggplot2::scale_color_manual(name = NULL, values = leg.col, labels = leg.name) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(shape = leg.shape, size = leg.size)))+
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA),
      panel.background = ggplot2::element_rect(fill = NA, color = "black"),
      legend.key = ggplot2::element_rect(fill = NA, color = "grey30"),
      legend.position = "bottom",
      panel.grid.major = element_blank(),
      plot.title = element_text(size = 24),
      plot.subtitle = element_text(size = 14),
      axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 20),
      legend.text = ggplot2::element_text(size = 20), legend.title = ggplot2::element_text(size = 12),
      plot.background = ggplot2::element_rect(fill = "white", color = "white")) -> sum_dotplot
  
 
  