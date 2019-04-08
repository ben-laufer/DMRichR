#' methylLearn
#' @description Performs random forest machine learning with repeated k-fold cross validation on significant DMRs to identify the most important DMRs
#' @param bsseq Smoothed bsseq object
#' @param regions Genomic ranges object
#' @param names Ordered sample names
#' @param groups Ordered test covariate information for each sample
#' @param k.fold number of folds in the cross-validation
#' @param repeats number of repeats in the cross-validation
#' @param ... Additional arguments passed onto randomForest::randomForest()
#' @return ?
#' @references \url{https://cran.rstudio.com/web/packages/randomForestExplainer/vignettes/randomForestExplainer.html}
#' @references \url{https://stats.stackexchange.com/questions/82162/cohens-kappa-in-plain-english}
#' @references \url{https://explained.ai/rf-importance/}
#' @references \url{https://machinelearningmastery.com/feature-selection-with-the-caret-r-package/}
#' @references \url{http://ml-tutorials.kyrcha.info/rfe.html}
#' @import bsseq
#' @import tidyverse
#' @import caret
#' @import randomForest
#' @import randomForestExplainer
#' @import e1071
#' @export methylLearn
methylLearn <- function(data = getMeth(BSseq = bs.filtered.bsseq, regions = sigRegions, type = "smooth", what = "perRegion") %>% as.matrix() %>% t(),
                        groups = bs.filtered.bsseq %>% pData() %>% as_tibble() %>% pull(!!testCovariate),
                        k.fold = 10,
                        repeats = 3,
                        ...){

# install(c("randomForest", "randomForestExplainer", "caret"))

# stopifnot(suppressMessages(sapply(c("randomForest", "randomForestExplainer", "caret", "e1071"), require, character.only = TRUE)))
# setwd("/Users/blaufer/Box Sync/DMRseq/DS")
# load("bismark.RData")
# load("bsseq.RData")
# load("DMRs.RData")
# load("GO.RData")
# load("Blocks.RData")


  # Tidy --------------------------------------------------------------------

  # Create a column to replace names for variable importantance
  colnames(data) <- cbind(as.data.frame(seqnames(sigRegions)), ranges(sigRegions)) %>%
    as_tibble() %>%
    dplyr::select(value,start,end) %>%
    tidyr::unite("bed", c("value","start","end"), sep = ".") %>%
    as.matrix()
  
  # groupMatrix <- groups %>% as.character() %>% as.matrix()
  # colnames(groupMatrix) <- "groups"
  # data <- cbind(groupMatrix, data)
  # rm(groups)
  
  # Random forest: Standard approach ----------------------------------------
  
  # Random forest cross-valdidation for feature selection  
  feature <- randomForest::rfcv(trainx = data,
                                trainy = groups,
                                cv.fold = k.fold,
                                ...)
  # Implement random forest
  forest <- randomForest(groups ~ .,
                         data = data,
                         importance = T,
                         ...)

  # Extract variable importance measure
  variableImportance <- importance(forest,
                                   type = 1,
                                   scale = F) %>%
    as.data.frame() %>% 
    tibble::rownames_to_column() %>% 
    dplyr::as_tibble() %>%
    dplyr::arrange(dplyr::desc(MeanDecreaseAccuracy))
  
  # Html report
  explain_forest(forest,
                 interactions = F, # True creates an error
                 data = data) 

  # Caret approach: RFE and RF ----------------------------------------------

  fitRfModel <- function(data){
    set.seed(5)
    
    # Remove highly correlated DMRs 
    
    # cor.idx <- data %>%
    #   cor() %>%
    #   findCorrelation(cutoff = 0.75,
    #                   exact = T,
    #                   names = T)
    # 
    # data <- data %>%
    #   tibble::as_tibble() %>%
    #   dplyr::select(-cor.idx) %>%
    #   as.matrix()
    
    # Recursive feature elimination
    # Automatically select a subset of the most predictive features, where a Random Forest algorithm is used on each iteration to evaluate the model. 
    
    rfeMe <- data %>%
      rfe(.,
          groups,
          rfeControl = rfeControl(functions = rfFuncs,
                                  method = "repeatedcv",
                                  number = k.fold,
                                  repeats = repeats)
      )
    
    # Summarize
    print(rfeMe)
    
    # List top predictors
    predictors(rfeMe)
    
    # Plot
    plot(rfeMe, type=c("g", "o"))
    
    # Fit rf on selected features
    model <- train(as.formula(paste("groups", paste(rfeMe$optVariables, collapse = " + "), sep = " ~ ")), # groups ~ .
                   data = data,
                   method = "rf",
                   tuneGrid = expand.grid(.mtry = length(rfeMe$optVariables)),
                   trControl = trainControl(method = "repeatedcv", 
                                            number = k.fold,
                                            repeats = repeats,
                                            verboseIter = TRUE,
                                            returnResamp = "final", #all
                                            savePredictions = "final", #all
                                            classProbs = TRUE) 
                   )
    return(model)
  }
  
  model <- fitRfModel(data)

  # Confusion matrix of cross-validated results
  confusionMatrix <- confusionMatrix.train(model,
                                           norm = "none")
  
  # Feature selection using variable importance
  varImpList <- varImp(object = model)
  
}
