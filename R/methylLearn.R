#' methylLearn
#' @description Performs random forest machine learning with k-fold cross validation on significant DMRs
#' @param bsseq Smoothed bsseq object
#' @param regions Genomic ranges object
#' @param names Ordered sample names
#' @param groups Ordered test covariate information for each sample
#' @param k.fold number of folds in the cross-validation
#' @param ... Additional arguments passed onto randomForest::randomForest()
#' @return ?
#' @references \url{https://cran.rstudio.com/web/packages/randomForestExplainer/vignettes/randomForestExplainer.html}
#' @references \url{https://stats.stackexchange.com/questions/82162/cohens-kappa-in-plain-english}
#' @references \url{https://explained.ai/rf-importance/}
#' @import bsseq
#' @import tidyverse
#' @import caret
#' @import randomForest
#' @import randomForestExplainer
#' @import e1071
#' @export methylLearn
methylLearn <- function(data = getMeth(BSseq = bs.filtered.bsseq, regions = sigRegions, type = "smooth", what = "perRegion") %>% as.matrix() %>% t(),
                        groups = bs.filtered.bsseq %>% pData() %>% as_tibble() %>% pull(!!testCovariate),
                        k.fold = 5,
                        ...){

# install(c("randomForest", "randomForestExplainer", "caret"))

# stopifnot(suppressMessages(sapply(c("randomForest", "randomForestExplainer", "caret", "e1071"), require, character.only = TRUE)))
# setwd("/Users/blaufer/Box Sync/DMRseq/DS")
# load("bismark.RData")
# load("bsseq.RData")
# load("DMRs.RData")
# load("GO.RData")
# load("Blocks.RData")

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
  
  # 5 fold cross validation random forest model -------------------------------------------------
  
  fitRfModel <- function(data){
    set.seed(5)
    model <- train(groups ~ .,
                   data = data,
                   method = "rf",
                   trControl = trainControl(method = "cv", 
                                            number = 5, 
                                            verboseIter = TRUE,
                                            returnResamp = "final",#all
                                            savePredictions = "final", #"all"
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


  # Recursive feature elimination  ------------------------------------------

  # Automatically selecting a subset of the most predictive features, where a Random Forest algorithm is used on each iteration to evaluate the model. 
  # https://machinelearningmastery.com/feature-selection-with-the-caret-r-package/
  
  rfeMe <- function(data){
    set.seed(5)
    results <- rfe(data,
                   groups,
                   rfeControl = rfeControl(functions = rfFuncs,
                                           method = "cv",
                                           number = 5)
    )
    return(results)
  }
  
  results <- rfeMe(data)
  
  # Summarize
  results
  
  # List top predictors
  predictors(results)
  
  # Plot
  plot(results, type=c("g", "o"))
  
}
