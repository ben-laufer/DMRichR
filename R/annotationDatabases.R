#' annotationDatabases
#' @description Assigns Bioconductor annotation databases (BSgenome, TxDb, org.db).
#' @param genome Character string of genome symbol (i.e. "hg38").
#' @param ensembl A logical indicating whether Ensembl annotations should be used. 
#' @return BSgenome, TxDb, org.db for genome of interest are loaded and
#'  assigned to the global environment.
#' @import BiocManager
#' @importFrom dplyr case_when
#' @importFrom glue glue
#' @import AnnotationHub
#' @import ensembldb
#' @export annotationDatabases
#' 
annotationDatabases <- function(genome = genome,
                                ensembl = FALSE){
  packages <- dplyr::case_when(genome == "hg38" ~ c("BSgenome.Hsapiens.UCSC.hg38",
                                                    "TxDb.Hsapiens.UCSC.hg38.knownGene",
                                                    "org.Hs.eg.db"),
                               genome == "hg19" ~ c("BSgenome.Hsapiens.UCSC.hg19",
                                                    "TxDb.Hsapiens.UCSC.hg19.knownGene",
                                                    "org.Hs.eg.db"),
                               genome == "mm10" ~ c("BSgenome.Mmusculus.UCSC.mm10",
                                                    "TxDb.Mmusculus.UCSC.mm10.knownGene",
                                                    "org.Mm.eg.db"),
                               genome == "mm9" ~ c("BSgenome.Mmusculus.UCSC.mm9",
                                                   "TxDb.Mmusculus.UCSC.mm9.knownGene",
                                                   "org.Mm.eg.db"),
                               genome == "rheMac10" ~ c("BSgenome.Mmulatta.UCSC.rheMac10",
                                                        "TxDb.Mmulatta.UCSC.rheMac10.refGene",
                                                        "org.Mmu.eg.db"),
                               genome == "rheMac8" ~ c("BSgenome.Mmulatta.UCSC.rheMac8",
                                                       "TxDb.Mmulatta.UCSC.rheMac8.refGene",
                                                       "org.Mmu.eg.db"),
                               genome == "rn6" ~ c("BSgenome.Rnorvegicus.UCSC.rn6",
                                                   "TxDb.Rnorvegicus.UCSC.rn6.refGene",
                                                   "org.Rn.eg.db"),
                               genome == "danRer11" ~ c("BSgenome.Drerio.UCSC.danRer11",
                                                        "TxDb.Drerio.UCSC.danRer11.refGene",
                                                        "org.Dr.eg.db"),
                               genome == "galGal6" ~ c("BSgenome.Ggallus.UCSC.galGal6",
                                                       "TxDb.Ggallus.UCSC.galGal6.refGene",
                                                       "org.Gg.eg.db"),
                               genome == "bosTau9" ~ c("BSgenome.Btaurus.UCSC.bosTau9",
                                                       "TxDb.Btaurus.UCSC.bosTau9.refGene",
                                                       "org.Bt.eg.db"),
                               genome == "panTro6" ~ c("BSgenome.Ptroglodytes.UCSC.panTro6",
                                                       "TxDb.Ptroglodytes.UCSC.panTro6.refGene",
                                                       "org.Pt.eg.db"),
                               genome == "dm6" ~ c("BSgenome.Dmelanogaster.UCSC.dm6",
                                                   "TxDb.Dmelanogaster.UCSC.dm6.ensGene",
                                                   "org.Dm.eg.db"),
                               genome == "susScr11" ~ c("BSgenome.Sscrofa.UCSC.susScr11",
                                                        "TxDb.Sscrofa.UCSC.susScr11.refGene",
                                                        "org.Ss.eg.db"),
                               genome == "canFam3" ~ c("BSgenome.Cfamiliaris.UCSC.canFam3",
                                                       "TxDb.Cfamiliaris.UCSC.canFam3.refGene",
                                                       "org.Cf.eg.db"),
                               # TAIR10 is an "annotation release" based on TAIR9.
                               genome == "TAIR10" ~ c("BSgenome.Athaliana.TAIR.TAIR9",
                                                      "TxDb.Athaliana.BioMart.plantsmart28",
                                                      "org.At.tair.db"),
                               genome == "TAIR9" ~ c("BSgenome.Athaliana.TAIR.TAIR9",
                                                     "TxDb.Athaliana.BioMart.plantsmart28",
                                                     "org.At.tair.db")
  )
  
  new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){
    glue::glue("Installing {new.packages}")
    suppressMessages(BiocManager::install(new.packages, ask = FALSE, quiet = TRUE))
    cat("Done", "\n")
  }
  print(glue::glue("Loading {packages}"))
  stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))
  
  if(genome == "hg38"){
    assign("goi", BSgenome.Hsapiens.UCSC.hg38, envir = .GlobalEnv)
    assign("TxDb", TxDb.Hsapiens.UCSC.hg38.knownGene, envir = .GlobalEnv)
    assign("annoDb", "org.Hs.eg.db", envir = .GlobalEnv)
  }else if(genome == "hg19"){
    assign("goi", BSgenome.Hsapiens.UCSC.hg19, envir = .GlobalEnv) 
    assign("TxDb", TxDb.Hsapiens.UCSC.hg19.knownGene, envir = .GlobalEnv)
    assign("annoDb", "org.Hs.eg.db", envir = .GlobalEnv)
  }else if(genome == "mm10"){
    assign("goi", BSgenome.Mmusculus.UCSC.mm10, envir = .GlobalEnv)
    assign("TxDb", TxDb.Mmusculus.UCSC.mm10.knownGene, envir = .GlobalEnv)
    assign("annoDb", "org.Mm.eg.db", envir = .GlobalEnv)
  }else if(genome == "mm9"){
    assign("goi", BSgenome.Mmusculus.UCSC.mm9, envir = .GlobalEnv)
    assign("TxDb", TxDb.Mmusculus.UCSC.mm9.knownGene, envir = .GlobalEnv)
    assign("annoDb", "org.Mm.eg.db", envir = .GlobalEnv)
  }else if(genome == "rheMac10"){
    assign("goi", BSgenome.Mmulatta.UCSC.rheMac10, envir = .GlobalEnv)
    assign("TxDb", TxDb.Mmulatta.UCSC.rheMac10.refGene, envir = .GlobalEnv)
    assign("annoDb", "org.Mmu.eg.db", envir = .GlobalEnv)
  }else if(genome == "rheMac8"){
    assign("goi", BSgenome.Mmulatta.UCSC.rheMac8, envir = .GlobalEnv)
    assign("TxDb", TxDb.Mmulatta.UCSC.rheMac8.refGene, envir = .GlobalEnv)
    assign("annoDb", "org.Mmu.eg.db", envir = .GlobalEnv)
  }else if(genome == "rn6"){
    assign("goi", BSgenome.Rnorvegicus.UCSC.rn6, envir = .GlobalEnv)
    assign("TxDb", TxDb.Rnorvegicus.UCSC.rn6.refGene, envir = .GlobalEnv)
    assign("annoDb", "org.Rn.eg.db", envir = .GlobalEnv)
  }else if(genome == "danRer11"){
    assign("goi", BSgenome.Drerio.UCSC.danRer11, envir = .GlobalEnv)
    assign("TxDb", TxDb.Drerio.UCSC.danRer11.refGene, envir = .GlobalEnv)
    assign("annoDb", "org.Dr.eg.db", envir = .GlobalEnv)
  }else if(genome == "galGal6"){
    assign("goi", BSgenome.Ggallus.UCSC.galGal6, envir = .GlobalEnv)
    assign("TxDb", TxDb.Ggallus.UCSC.galGal6.refGene, envir = .GlobalEnv)
    assign("annoDb", "org.Gg.eg.db", envir = .GlobalEnv)
  }else if(genome == "bosTau9"){
    assign("goi", BSgenome.Btaurus.UCSC.bosTau9, envir = .GlobalEnv)
    assign("TxDb", TxDb.Btaurus.UCSC.bosTau9.refGene, envir = .GlobalEnv)
    assign("annoDb", "org.Bt.eg.db", envir = .GlobalEnv)
  }else if(genome == "panTro6"){
    assign("goi", BSgenome.Ptroglodytes.UCSC.panTro6, envir = .GlobalEnv)
    assign("TxDb", TxDb.Ptroglodytes.UCSC.panTro6.refGene, envir = .GlobalEnv)
    assign("annoDb", "org.Pt.eg.db", envir = .GlobalEnv)
  }else if(genome == "dm6"){
    assign("goi", BSgenome.Dmelanogaster.UCSC.dm6, envir = .GlobalEnv)
    assign("TxDb", TxDb.Dmelanogaster.UCSC.dm6.ensGene, envir = .GlobalEnv)
    assign("annoDb", "org.Dm.eg.db", envir = .GlobalEnv)
  }else if(genome == "susScr11"){
    assign("goi", BSgenome.Sscrofa.UCSC.susScr11, envir = .GlobalEnv)
    assign("TxDb", TxDb.Sscrofa.UCSC.susScr11.refGene, envir = .GlobalEnv)
    assign("annoDb", "org.Ss.eg.db", envir = .GlobalEnv)
  }else if(genome == "canFam3"){
    assign("goi", BSgenome.Cfamiliaris.UCSC.canFam3, envir = .GlobalEnv)
    assign("TxDb", TxDb.Cfamiliaris.UCSC.canFam3.refGene, envir = .GlobalEnv)
    assign("annoDb", "org.Cf.eg.db", envir = .GlobalEnv)
  }else if(genome %in% c("TAIR9", "TAIR10")){
    assign("goi", BSgenome.Athaliana.TAIR.TAIR9, envir = .GlobalEnv)
    assign("TxDb", TxDb.Athaliana.BioMart.plantsmart28, envir = .GlobalEnv)
    assign("annoDb", "org.At.tair.db", envir = .GlobalEnv)
  }else{
    stop(glue::glue("{genome} is not supported, please choose either hg38, hg19, mm10, mm9, \\
    rheMac10, rheMac8, rn6, danRer11, galGal6, bosTau9, panTro6, dm6, susScr11, canFam3, TAIR10, \\
    or TAIR9 [Case Sensitive]"))
  }
  
  if(ensembl == TRUE){
    if(genome %in% c("hg38", "mm10", "rheMac10", "rheMac8", "rn6", "danRer11", "galGal6",
                     "bosTau9", "dm6", "susScr11", "canFam3")){
      
      print(glue::glue("EnsemblDb annotations for {genome} will be used."))
      # AnnotationHub::display(AnnotationHub::query(AnnotationHub::AnnotationHub(), c("EnsDb", "Pan troglodytes")))
      
      ahCode <- dplyr::case_when(genome == "hg38" ~ "AH83216",
                                 genome == "mm10" ~ "AH83247",
                                 genome == "rheMac10" ~ "AH83244",
                                 genome == "rheMac8" ~ "AH73903",
                                 genome == "rn6" ~ "AH83318",
                                 genome == "danRer11" ~ "AH83189",
                                 genome == "galGal6" ~ "AH83209",
                                 genome == "bosTau9" ~ "AH83145",
                                 genome == "dm6" ~ "AH83185",
                                 genome == "susScr11" ~ "AH83340",
                                 genome == "canFam3" ~ "AH78741")
      
      print(glue::glue("Your AnnotationHub code is {ahCode}."))
      TxDb <- AnnotationHub::AnnotationHub()[[ahCode]]
      
    }else if(genome %in% c("hg19", "mm9", "panTro6", "TRAIR9", "TAIR10")){
      
      print(glue::glue("EnsemblDB annotations for {genome} are not available, \\
                   the default Bioconductor TxDb will be used."))
    }
  }
}
