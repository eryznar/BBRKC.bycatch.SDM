


thresh.select <- function(model_b, test){
  
  # Calculate testing data prevalance
  test.prev <- nrow(test %>% filter(PA == 1))/nrow(test)
  
  # Make prediction using test data
  pred <- suppressWarnings(predict.gbm(model_b, # fitted model to predict
                                       test, # data to predict to
                                       n.trees=model_b$gbm.call$best.trees, # see help
                                       type="response")) # predict probabilities
  
  # Calculate AUC
  obs <- test$PA 
  
  d <- cbind(obs, pred)
  pres <- d[d[,1]==1, 2]
  abs <- d[d[,1]==0, 2]
  e <- dismo::evaluate(p=pres, a=abs)
  AUC <- e@auc
  
  # Extract unique predicted probabilities
  thresholds <- unique(pred)
  
  # Sort the thresholds in descending order
  thresholds <- sort(thresholds, decreasing = TRUE)
  
  # Create empty vectors to store diagnostics
  tpr <- numeric(length(thresholds))
  fpr <- numeric(length(thresholds))
  gmean <- numeric(length(thresholds))
  kappa <- numeric(length(thresholds))
  pred.prev <- numeric(length(thresholds))
  
  
  # Calculate diagnostics for each threshold
  for (i in seq_along(thresholds)) {
    
    threshold <- thresholds[i]
    
    
    pred2 <- cbind(prediction = pred,
                        test) %>%
      mutate(pred_state = factor(ifelse(prediction >= threshold,
                                        1, 0),
                                 levels = c(1, 0)),
             obs_state = factor(PA, levels = c(1, 0)))
    
       # Create a confusion matrix
      cm <- table(Predicted = pred2$pred_state, 
                  Actual = pred2$obs_state)
      
      # Calculate TPR and FPR
      if(sum(dim(cm))==4){
        tpr[i] <- cm["1", "1"] / (cm["1", "1"] + cm["0", "1"])
        fpr[i] <- cm["1", "0"] / (cm["1", "0"] + cm["0", "0"])
        #gmean[i] <- sqrt(tpr[i]*(1-fpr[i]))
        gmean[i] <- tpr[i]-fpr[i]
        kappa[i] = cohen.kappa(x = cbind(pred2$pred_state, pred2$obs_state))$kappa
        pred.prev[i] = nrow(pred2 %>% filter(pred_state==1))/nrow(pred2)
      }
    }
    # Create data frame of diagnostics
    data_frame(thresh = thresholds, tpr = tpr, fpr = fpr, 
               gmean = gmean, kappa = kappa, pred.prev = pred.prev) -> diag
    
    # Select threshold with max kappa
    diag %>%
      na.omit() %>%
      filter(kappa == max(kappa)) %>%
      dplyr::select(thresh, pred.prev) -> tt
    
    return(list(diag = diag, t.select = data.frame(AUC = AUC, max.kappa.thresh = tt$thresh,
                                        test.prev = test.prev, pred.prev = tt$pred.prev)))
}

# Legal males
model_b <- readRDS("./Models/BRT_legalmale_modelb_BEST.CPUE.rda")
test <- read.csv("./Data/lm_test.csv") 

thresh.select(model_b, test) -> lm.out

# Immature males
model_b <- readRDS("./Models/BRT_immaturemale_modelb_BEST.CPUE.rda")
test <- read.csv("./Data/im_test.csv") 

thresh.select(model_b, test) -> im.out

# Mature females
model_b <- readRDS("./Models/BRT_maturefemale_modelb_BEST.CPUE.rda")
test <- read.csv("./Data/mf_test.csv") 

thresh.select(model_b, test) -> mf.out

# Immature females
model_b <- readRDS("./Models/BRT_immaturefemale_modelb_BEST.CPUE.rda")
test <- read.csv("./Data/imf_test.csv") 

thresh.select(model_b, test) -> imf.out


rbind(lm.out$t.select %>% mutate(category = "Legal male"),
      im.out$t.select %>% mutate(category = "Immature male"),
      mf.out$t.select %>% mutate(category = "Mature female"),
      imf.out$t.select %>% mutate(category = "Immature female")) -> thresh.out


rbind(lm.out$diag %>% mutate(category = "Legal male"),
      im.out$diag %>% mutate(category = "Immature male"),
      mf.out$diag %>% mutate(category = "Mature female"),
      imf.out$diag %>% mutate(category = "Immature female")) -> diagnostic.out

write.csv(thresh.out, "./Output/thresh.selection.csv")
write.csv(diagnostic.out, "./Output/model.diag.csv")


read.csv("./Output/model.diag.csv") -> tt
read.csv("./Output/thresh.selection.csv") -> kk

tt
head(tt)

kk %>%
  dplyr::select(test.prev, category) -> hh

right_join(hh, tt, by = "category") %>%
  dplyr::select(!X) %>%
  mutate(prev.diff = abs(test.prev - pred.prev)) -> gg

gg %>%
  group_by(category) %>%
  filter(prev.diff == min(prev.diff)) -> ff
