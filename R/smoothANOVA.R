#' smoothANOVA
#' @description Perform an ANOVA on smoothed methylation values averaged globally or across chromosomes
#' @param smoothAvg Tibble of average smoothed methylation values and covariates from getGlobal() or getChrom()
#' @return A list summariazing the ANOVA(s), where a clean statistical structure means no random effects, multiple testing corrections, or type III ANOVAs are needed
#' @references \url{https://cran.r-project.org/web/packages/broom/vignettes/broom_and_dplyr.html}
#' @import lsmeans
#' @import tidyverse
#' @import broom
#' @import openxlsx
#' @export smoothANOVA
smoothANOVA <- function(smoothAvg = smoothAvg){
  message("Peforming ANOVA")
  
  # Global tidier 
  globalPipe <- . %>% 
    tidy %>% 
    list("Anova" = .,
         "input" = smoothAvg) %>% 
    return()
  
  # Chrom tidier
  chromTidy <- function(){
    tidyFit <- models %>%
      dplyr::select(chromosome, tidyFit) %>%
      unnest
    tidyAnova <- models %>%
      dplyr::select(chromosome, tidyAnova) %>%
      unnest() %>%
      dplyr::select(chromosome, term, p.value)  %>%
      spread(key = term, value = p.value)  %>%
      dplyr::select(-Residuals)
    pairWise <- models %>%
      dplyr::select(chromosome, pairWise) %>%
      unnest  %>%
      mutate(fdr = p.adjust(p.value, method = 'fdr'))
    return(list("pairWise" = pairWise,
                "Anova p-values" = tidyAnova,
                "lm" = tidyFit,
                "input" = smoothAvg))
  }
  
  # getGlobal()  
  if(!("adjustCovariate" %in% names(smoothAvg)) &
     !("matchCovariate" %in% names(smoothAvg)) &
     !("chromosome" %in% colnames(smoothAvg))){
    aov(CpG_Avg ~ testCovariate, data = smoothAvg) %>%
      globalPipe
    
  }else if("adjustCovariate" %in% names(smoothAvg) &
           (!("matchCovariate" %in% names(smoothAvg)) | (length(levels(smoothAvg$matchCovariate))) <= 1) &
           !("chromosome" %in% colnames(smoothAvg))){
      aov(CpG_Avg ~ testCovariate + adjustCovariate, data = smoothAvg) %>%
      globalPipe
    
  }else if(!("adjustCovariate" %in% names(smoothAvg)) &
           ("matchCovariate" %in% names(smoothAvg) | (!(length(levels(smoothAvg$matchCovariate))) <= 1)) &
           !("chromosome" %in% colnames(smoothAvg))){
    aov(CpG_Avg ~ testCovariate + matchCovariate, data = smoothAvg) %>%
      globalPipe
    
  }else if("adjustCovariate" %in% names(smoothAvg) &
           ("matchCovariate" %in% names(smoothAvg) & (!(length(levels(smoothAvg$matchCovariate))) <= 1)) &
           !("chromosome" %in% colnames(smoothAvg))){
    aov(CpG_Avg ~ testCovariate + matchCovariate + adjustCovariate, data = smoothAvg) %>%
      globalPipe
  
  # getChrom()
  }else if(!("adjustCovariate" %in% names(smoothAvg)) &
           !("matchCovariate" %in% names(smoothAvg)) &
           ("chromosome" %in% colnames(smoothAvg))){
    models <- smoothAvg %>%
      nest(-chromosome) %>%
      mutate(
        fit = map(data, ~ lm(CpG_Avg ~ testCovariate, data = .x)),
        tidyFit = map(fit, tidy),
        anova = map(data, ~ aov(CpG_Avg ~ testCovariate , data = .x)),
        tidyAnova = map(anova, tidy),
        pairWise = map(data, ~ lm(CpG_Avg ~ testCovariate, data = .x) %>%
                         ref.grid() %>%
                         lsmeans(~testCovariate) %>%
                         pairs() %>%
                         summary())
      )
    chromTidy()
    
  }else if("adjustCovariate" %in% names(smoothAvg) &
           (!("matchCovariate" %in% names(smoothAvg)) | (length(levels(smoothAvg$matchCovariate))) <= 1) &
           ("chromosome" %in% colnames(smoothAvg))){
    models <- smoothAvg %>%
      nest(-chromosome) %>%
      mutate(
        fit = map(data, ~ lm(CpG_Avg ~ testCovariate + adjustCovariate, data = .x)),
        tidyFit = map(fit, tidy),
        anova = map(data, ~ aov(CpG_Avg ~ testCovariate + adjustCovariate, data = .x)),
        tidyAnova = map(anova, tidy),
        pairWise = map(data, ~ lm(CpG_Avg ~ testCovariate + adjustCovariate, data = .x) %>%
                         ref.grid() %>%
                         lsmeans(~testCovariate) %>%
                         pairs() %>%
                         summary())
      )
    chromTidy()
    
  }else if(!("adjustCovariate" %in% names(smoothAvg)) &
           ("matchCovariate" %in% names(smoothAvg) | (!(length(levels(smoothAvg$matchCovariate))) <= 1)) &
           ("chromosome" %in% colnames(smoothAvg))){
    models <- smoothAvg %>%
      nest(-chromosome) %>%
      mutate(
        fit = map(data, ~ lm(CpG_Avg ~ testCovariate + matchCovariate, data = .x)),
        tidyFit = map(fit, tidy),
        anova = map(data, ~ aov(CpG_Avg ~ testCovariate + matchCovariate, data = .x)),
        tidyAnova = map(anova, tidy),
        pairWise = map(data, ~ lm(CpG_Avg ~ testCovariate + matchCovariate, data = .x) %>%
                         ref.grid() %>%
                         lsmeans(~testCovariate) %>%
                         pairs() %>%
                         summary())
      )
    chromTidy()
    
  }else if("adjustCovariate" %in% names(smoothAvg) &
           ("matchCovariate" %in% names(smoothAvg) & (!(length(levels(smoothAvg$matchCovariate))) <= 1)) &
           ("chromosome" %in% colnames(smoothAvg))){
    models <- smoothAvg %>%
      nest(-chromosome) %>%
      mutate(
        fit = map(data, ~ lm(CpG_Avg ~ testCovariate + adjustCovariate + matchCovariate, data = .x)),
        tidyFit = map(fit, tidy),
        anova = map(data, ~ aov(CpG_Avg ~ testCovariate + adjustCovariate + matchCovariate, data = .x)),
        tidyAnova = map(anova, tidy),
        pairWise = map(data, ~ lm(CpG_Avg ~ testCovariate + adjustCovariate + matchCovariate, data = .x) %>%
                         ref.grid() %>%
                         lsmeans(~testCovariate) %>%
                         pairs() %>%
                         summary())
      )
    chromTidy()
  }
}
