# PURPOSE:
# To fit boosted regression trees to model Bristol Bay red king crab bycatch occurrence and abundance

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

# SPECIFY FUNCTIONS ---------------------------------------------------------------------------------------
  # Fit boosted regression tree models to random iterative splits of 80/20 train/test of catch data
    # @param train: training data (options = "lm_train", "im_train", "mf_train", "imf_train")
    # @param test: testing data (options = "lm_test", "im_test", "mf_test", "imf_test")
    # @param matsex: red king crab sex-size/maturity category (options = "legalmale", "immaturemale", "maturefemale", "immaturefemale")
    # @param iteration: random training/testing data split iteration (options = 1:10)
    # @return Returns data frame of model performance diagnostics for binomial (bycatch occurrence) and Poisson (bycatch abundance) models
  
    model_iter <- function(train, test, matsex, iteration){

        # Filter by iteration
        train2 <- train %>% 
          filter(iter %in% iteration, is.infinite(YFS_dfish) == "FALSE", is.infinite(RS_dfish) == FALSE) %>%
          mutate(ELV_SWP = as.factor(ELV_SWP),
                 Period = as.factor(Period))
        
        test2 <- test %>% 
          filter(iter %in% iteration, is.infinite(YFS_dfish) == "FALSE", is.infinite(RS_dfish) == FALSE)%>%
          mutate(ELV_SWP = as.factor(ELV_SWP),
                 Period = as.factor(Period))
        
        # Specify non-predictor variables
        non_preds<- c("x", "y", "predict_year", "extrap_log_rd", "PA", "iter")
        
        train2 <- train2[,-1]
        test2 <- test2[,-1]
        
        
        # Fit binomial model
        model_b <- gbm.step(data = train2, 
                            gbm.x = which(!colnames(train2) %in% non_preds), # columns of predictors
                            gbm.y = which(colnames(train2) == "PA"), # column of response
                            family = "bernoulli", # for PA data
                            tree.complexity = 5, # model interactions, etc?
                            learning.rate = 0.005, # influence of each tree 
                            bag.fraction = 0.5) 
        #saveRDS(model_b, paste0("./Data/model_b_CPUE", iteration, "_", matsex, ".rda"))
        
        ntreesb <- model_b$n.trees
        
        model_p <- gbm.step(data = train2 %>% filter(extrap_log_rd >0), 
                            gbm.x = which(!colnames(train2) %in% non_preds), # columns of predictors 
                            gbm.y = which(colnames(train2) == "extrap_log_rd"), # column of response
                            family = "poisson", # for count data
                            tree.complexity = 5, # model interactions, etc?
                            learning.rate = 0.005, # influence of each tree 
                            bag.fraction = 0.5)
        #saveRDS(model_p, paste0("./Data/model_p_CPUE", iteration, "_", matsex, ".rda"))
        
        ntreesp <- model_p$n.trees
        
        # Calculate AUC for binomial
        pred <- suppressWarnings(predict.gbm(model_b, # fitted model to predict
                            test2, # data to predict to
                            n.trees=model_b$gbm.call$best.trees, # see help
                            type="response")) # predict probabilities
        
        obs <- test2$PA 
        
        d <- cbind(obs, pred)
        pres <- d[d[,1]==1, 2]
        abs <- d[d[,1]==0, 2]
        e <- dismo::evaluate(p=pres, a=abs)
        AUC <- e@auc
        
        
        # Calculate RMSE and spearman's for positive model
        pred <- suppressWarnings(predict.gbm(model_p, # fitted model to predict
                            test2 %>% filter(extrap_log_rd>0), # data to predict to
                            n.trees=model_p$gbm.call$best.trees, # see help
                            type="response")) # predict probabilities
        
        obs <- test2$extrap_log_rd[which(test2$extrap_log_rd >0)]
        
        RMSE <- sqrt(mean((obs-pred)^2))
        
        d <- as.data.frame(cbind(obs, pred))
        
        ggplot(d, aes(x = obs, y = pred))+
          geom_point()+
          geom_smooth(method = "lm") + 
          theme_classic()+
          xlab("Observed vals")+
          ylab("Predicted vals") -> pred_obs_plot
        
        cor.test(d$obs, d$pred, method = "spearman") -> cor_out
        
        rho <- cor_out$estimate
        p <- cor_out$p.value
        
        # Bind metrics into a data frame
        eval_df <- data.frame(ntreesb = ntreesb, ntreesp = ntreesp, AUC = AUC, RMSE = RMSE, rho = rho, p = p, 
                              iter = iteration)
        
        return(list(eval_df))
    }

# RUN FUNCTION --------------------------------------------------------------------------------------------
  #Specify # of iterations
  iteration <- 1:10
  
  # LEGAL MALES
    matsex = "legalmale"
    train = lm_train
    test = lm_test

    # Fit and save models for each training/testing split iteration
    iteration %>%
      map_df(~model_iter(train, test, matsex, .x)) -> iter_out

    # Save model diagnostic data frame
    eval_df <- data.frame(ntreesb = iter_out$ntreesb, ntreesp = iter_out$ntreesp, AUC = iter_out$AUC,
                          RMSE = iter_out$RMSE, rho = iter_out$rho, p = iter_out$p, iter = iter_out$iter)

  # IMMATURE MALES 
    matsex = "immaturemale"
    train = im_train
    test = im_test

    # Fit and save models for each training/testing split iteration
      iteration %>%
        map_df(~model_iter(train, test, matsex, .x)) -> iter_out
      
    # Save model diagnostic data frame
    lm_eval_df <- data.frame(ntreesb = iter_out$ntreesb, ntreesp = iter_out$ntreesp, AUC = iter_out$AUC,
                            RMSE = iter_out$RMSE, rho = iter_out$rho, p = iter_out$p, iter = iter_out$iter)

    #saveRDS(lm_eval_df, "./Output/Bycatch_legalmale_iterCPUE.rda")
    
  # MATURE FEMALES
    matsex = "maturefemale"
    train = mf_train
    test = mf_test
    
    # Fit and save models for each training/testing split iteration
    iteration %>%
      map_df(~model_iter(train, test, matsex, .x)) -> iter_out
    
    # Save model diagnostic data frame
    im_eval_df <- data.frame(ntreesb = iter_out$ntreesb, ntreesp = iter_out$ntreesp, AUC = iter_out$AUC,
                          RMSE = iter_out$RMSE, rho = iter_out$rho, p = iter_out$p, iter = iter_out$iter)
    
    #saveRDS(im_eval_df, "./Output/Bycatch_immaturemale_iterCPUE.rda")
    
  # MATURE FEMALES
    matsex = "immaturefemale"
    train = imf_train
    test = imf_test
    
    # Fit and save models for each training/testing split iteration
    iteration %>%
      map_df(~model_iter(train, test, matsex, .x)) -> iter_out
    
    # Save model diagnostic data frame
    mf_eval_df <- data.frame(ntreesb = iter_out$ntreesb, ntreesp = iter_out$ntreesp, AUC = iter_out$AUC,
                          RMSE = iter_out$RMSE, rho = iter_out$rho, p = iter_out$p, iter = iter_out$iter)
  
    #saveRDS(mf_eval_df, "./Output/Bycatch_maturefemale_iterCPUE.rda")

  # IMMATURE FEMALES
    matsex = "maturefemale"
    train = mf_train
    test = mf_test
    
    # Fit and save models for each training/testing split iteration
    iteration %>%
      map_df(~model_iter(train, test, matsex, .x)) -> iter_out
    
    # Save model diagnostic data frame
    imf_eval_df <- data.frame(ntreesb = iter_out$ntreesb, ntreesp = iter_out$ntreesp, AUC = iter_out$AUC,
                          RMSE = iter_out$RMSE, rho = iter_out$rho, p = iter_out$p, iter = iter_out$iter)
    
      #saveRDS(imf_eval_df, "./Output/Bycatch_immaturefemale_iterCPUE.rda")
    
# LOAD and SAVE BEST MODELS ---------------------------------------------------------------
  # LEGAL MALE
    lm_eval_df <- readRDS("./Output/Bycatch_legalmale_iterCPUE.rda")
    
    lm_iter = 9 # Best training/testing iteration based on AUC, RMSE, and RHO

    # # Identify and save best models (iteration # 9) (THESE CAN BE SAVED/ACCESSED IN MODEL_ITER FUNCTION ABOVE)
    # lm_modelb_BEST <- readRDS("./Models/model_b_CPUE9_legalmale.rda") %>%
    #   saveRDS("./Models/BRT_legalmale_modelb_BEST.CPUE.rda")
    # 
    # lm_modelp_BEST <- readRDS("./Models/model_p_CPUE9_legalmale.rda") %>%
    #   saveRDS("./Models/BRT_legalmale_modelp_BEST.CPUE.rda")
    
  # IMMATURE MALE
    im_eval_df <- readRDS("./Output/Bycatch_immaturemale_iterCPUE.rda")
    
    im_iter = 10 # Best training/testing iteration based on AUC, RMSE, and RHO
    
    # # Identify and save best models (iteration # 10)
    # im_modelb_BEST <- readRDS("./Models/model_b_CPUE10_immaturemale.rda") %>%
    #   saveRDS("./Models/BRT_immaturemale_modelb_BEST.CPUE.rda")
    # 
    # im_modelp_BEST <- readRDS("./Models/model_p_CPUE10_immaturemale.rda") %>%
    #   saveRDS("./Models/BRT_immaturemale_modelp_BEST.CPUE.rda")

  # MATURE FEMALE
    mf_eval_df <- readRDS("./Output/Bycatch_maturefemale_iterCPUE.rda")
    
    mf_iter = 8 # Best training/testing iteration based on AUC, RMSE, and RHO
    
    # # Identify and save best models (iteration # 8)
    # mf_modelb_BEST <- readRDS("./Models/model_b_CPUE8_maturefemale.rda") %>%
    #   saveRDS("./Models/BRT_maturefemale_modelb_BEST.CPUE.rda")
    # 
    # mf_modelp_BEST <- readRDS("./Models/model_p_CPUE8_maturefemale.rda") %>%
    #   saveRDS("./Models/BRT_maturefemale_modelp_BEST.CPUE.rda")
  
  # IMMATURE FEMALE
    imf_eval_df <- readRDS("./Output/Bycatch_immaturefemale_iterCPUE.rda")
    
    imf_iter = 1 # Best training/testing iteration based on AUC, RMSE, and RHO
  
    # # Identify and save best models (iteration # 1)
    # imf_modelb_BEST <- readRDS("./Models/model_b_CPUE1_immaturefemale.rda") %>%
    #   saveRDS("./Models/BRT_immaturefemale_modelb_BEST.CPUE.rda")
    # 
    # imf_modelp_BEST <- readRDS("./Models/model_p_CPUE1_immaturefemale.rda") %>%
    #   saveRDS("./Models/BRT_immaturefemale_modelp_BEST.CPUE.rda")
    

