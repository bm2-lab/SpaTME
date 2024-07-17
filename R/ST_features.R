setClass("STFeature",
         slots = list(
           Position = "data.frame",
           CellAnno = "data.frame",
           CellCodis = "data.frame",
           ImmInfil = "data.frame",
           GsetSig = "list",
           SpeGenes = "list",
           LRscore = "data.frame",
           PhenoScore = "data.frame",
           Annotation = "data.frame",
           Others = "ANY"))

#' Create STFeature object
#' @param st Seurat object for ST data.
#' @param stf STFeature object. If available, overwrite it with specified
#'   features.
#' @param assay Name of \code{assay} in \code{st}.
#' @param slot Name of \code{slot} in \code{st}.
#' @param cell.anno A data.frame containing annotation of cell types. Rows are
#'   spots and columns are cell types. The values can be abundances or
#'   proportions of cell types.
#' @param norm Logical value. If TRUE, the sum of abundances of each spot will
#'   normalized to 1.
#' @param min.prop Cell types with normalized abundance below this value will be
#'   set to 0. Default is 0.
#' @param init.fea Features to be initialized. It should be selected from
#'   \code{c("Position", "CellCodis", "SpeGenes")}(Default).
#' @param verbose Logical value. If TRUE, print messages.
#' @return A STFeature object.
#' @export
#' @importFrom Seurat GetAssayData
#' @examples
CreateStfObj <- function(st,
                         stf = NULL,
                         assay = "SCT",
                         slot = "data",
                         cell.anno = NULL,
                         norm = T,
                         min.prop = 0,
                         init.fea = c("Position", "CellAnno", "CellCodis"),
                         verbose = T){
  if (is.null(stf)) {
    stf <- new("STFeature")
  }
  if ("Position" %in% init.fea) {
    if (! all(c("x", "y") %in% colnames(st@meta.data))) {
      stop("x,y information are missing in meta.data, please run STcoordCheck !")
    }
    stf@Position <- st@meta.data[, c("x", "y")]
  }
  if ("CellAnno" %in% init.fea) {
    if (is.null(cell.anno)) {
      stop("Cell annotation should be provided with argument 'cell.anno'!")
    }
    if (length(setdiff(rownames(cell.anno), colnames(st)))) {
      if (verbose) {
        warning("Some spots/cells in cell.anno not present in st data.")
      }
      cell.anno <- cell.anno[intersect(rownames(cell.anno), colnames(st)), , drop = F]
    }
    new.data <- data.frame(row.names = colnames(st))
    new.data[, colnames(cell.anno)] <- cell.anno[rownames(new.data), , drop = F]
    norm.data <- new.data/rowSums(new.data)
    norm.data[norm.data <= min.prop] <- 0
    # re-norm
    norm.data <- norm.data/rowSums(norm.data)
    if (norm) {
      stf@CellAnno <- data.frame(norm.data, check.names = F)
    } else {
      stf@CellAnno <- data.frame(new.data, check.names = F)
    }
  }
  if ("CellCodis" %in% init.fea) {
    stf@CellCodis <- CalCellCodis(stf@CellAnno, sort = T)
  }
  return(stf)
}

#' Calculation of Shannon's Entropy
#'
#' @param x a numeric probability vector.
#' @param base a character string specifying the logarithm unit. The default
#'   value is 2.
#'
#' @return a numeric value representing Shannon's Entropy
CalEntropy <- function(x, base = 2){
  # Check the the validity of input probability distributions
  if (anyNA(x)){
    stop("The input vector includes NA values.")
  }
  if (max(x) >1 | min(x) <0) {
    stop("The probability values are not between: [0,1].")
  }
  if (sum(x) > 1.001) {
    stop("The probability distribution does not sum to 1.")
  }
  entropy <- -sum(x[x != 0]*log(x[x != 0], base = base))
  return(entropy)
}

#' Measuring the degree of immune infiltration using two measures: immune
#' enrichment and immune diversity
#' @param data A matrix or data.frame containing proportions of
#'   immune cells. Rows are spots and columns are cell types.
#' @param base base a character string specifying the logarithm unit. The
#'   default value is 2.
#' @param min.prop Immune cells with a proportion less than this value will be set
#'   to 0.
#' @return A data.frame containing two columns: Imm.diversity and
#'   Imm.infiltration
#' @export
#' @examples
CalImmInfiltration <- function(data,
                               base = 2,
                               min.prop = 0.05){
  if (anyNA(data)){
    stop("The input data includes NA values.")
  }
  data[data < min.prop] <- 0
  # immune enrichment
  Imm.enrichment <- rowSums(data)
  # immune diversity
  prob <- data/Imm.enrichment
  prob[Imm.enrichment == 0,] <- 0
  Imm.diversity <- apply(prob, 1, CalEntropy, base = base)
  Imm.infiltration <- cbind(Imm.enrichment, Imm.diversity)
  return(data.frame(Imm.infiltration))
}

#' Calculation of cell co-distribution scores.
#' @param data A matrix or data.frame containing abundances or proportions of
#'   cell types. Rows are spots and columns are cell types.
#' @param sort Logical value indicating whether to sort the names of cells.
#' @return A data.frame of co-distribution scores for all cell-cell pairs.
#' @export
#' @examples
CalCellCodis <- function(data, sort = T){
  data <- Norm01(data)
  CellCodis <- combn(x = colnames(data),
                     m = 2,
                     FUN = function(x){
                       data[, x[1]] * data[, x[2]]
                     })
  colnames(CellCodis) <- combn(x = colnames(data),
                               m = 2,
                               FUN = function(x, s = sort){
                                 if (s){
                                   paste(sort(c(x[1], x[2])), collapse = "_")
                                 } else{
                                   paste(c(x[1], x[2]), collapse = "_")
                                 }
                               },s = sort)
  rownames(CellCodis) <- rownames(data)
  CellCodis <- data.frame(CellCodis, check.names = F)
  return(CellCodis)
}

#' Normalizing numeric data between 0 and 1.
#' @param data A numeric vector, matrix or data.frame.
#' @param constant.v Constent value when there is no varience in \code{data}.
#'   Default is 0.
#' @export
#' @return Normalized data.
Norm01 <- function(data, constant.v = 0){
  if (! is.null(dim(data))) {
    res <- apply(data, 2, function(x){
      if (length(unique(x)) ==1) {
        x <- rep(constant.v, length(x))
      } else {
        (x-min(x[!is.na(x)]))/(max(x[!is.na(x)])-min(x[!is.na(x)]))
      }
    })
    res <- as.matrix(res)
  } else {
    if (length(unique(data)) ==1) {
      res <- rep(constant.v, length(data))
    } else {
      res <- (data-min(data[!is.na(data)]))/(max(data[!is.na(data)])-min(data[!is.na(data)]))
    }
  }
  return(res)
}

#' Enrichment scores of gene sets.
#'
#' @param expr Expression data. Rows are genes and Columns are spots.
#' @param geneset A list of vectors of features; each entry should be named with
#'   feature.
#' @param method Method for enrichment analysis.  Available options are:
#' \itemize{
#'   \item{"AddModuleScore"}: Calculate the average expression levels of each
#'   gene set for each spot, subtracted by the aggregated expression of control
#'   feature sets. This method is conducted using \pkg{Seurat}.
#'   \item{"UCell"}: Calculate module enrichment scores using \pkg{UCell}.
#'   \item{"AUCell"}: Calculates the 'AUC' for each gene-set in each spot uising
#'   functions in \pkg{AUCell}.
#'   \item{"gsva"}: Estimates enrichment scores
#'   using \pkg{GSVA}.
#' }
#' @param gsva.method Method employed in the estimation of \pkg{GSVA}. default
#'   is \code{ssgsva}.
#' @param scale Logical value. If TRUE, scaling and centering the scores for all spots.
#' @param verbose Logical value. If TRUE, print messages.
#' @importFrom SeuratObject PackageCheck
#' @return A data.frame containing the enrichment scores.
#' @export
GsetScore <- function(expr,
                      geneset,
                      method = "AddModuleScore",
                      gsva.method = "ssgsea",
                      scale = F,
                      verbose = T){
  # method can be ssGSEA, AUCell or AddModuleScore

  if (method == "AddModuleScore") {
    empty.sets <- lapply(geneset, function(x){
      length(intersect(x, rownames(expr))) > 0
    })
    geneset2 <- geneset[unlist(empty.sets)]
    temp <- CreateSeuratObject(counts = expr)
    n1 <- ncol(temp@meta.data)
    temp <- suppressWarnings(AddModuleScore(temp,
                                            features = geneset2,
                                            name = names(geneset2),
                                            verbose = verbose))
    n2 <- ncol(temp@meta.data)
    score <- temp@meta.data[, (n1+1) : n2, drop = F]
    colnames(score) <- names(geneset2)
    res <- matrix(data = NA,
                  ncol = length(geneset),
                  nrow = nrow(score))
    rownames(res) <- rownames(score)
    colnames(res) <- names(geneset)
    res <- data.frame(res, check.names = F)
    res[, colnames(score)] <- score
    score <- res
  }
  if (method == "AUCell" ) {
    if (!PackageCheck('AUCell', error = FALSE)) {
      stop("Please install AUCell!")
    }
    cells_rankings <- AUCell::AUCell_buildRankings(exprMat = expr,
                                           plotStats=FALSE,
                                           verbose = verbose)
    cells_AUC <- AUCell::AUCell_calcAUC(geneSets = geneset,
                                rankings = cells_rankings,
                                verbose = verbose)
    score <- AUCell::getAUC(cells_AUC)
    score <- data.frame(t(score), check.names = F)
  }
  if (method == "gsva") {
    if (!PackageCheck('GSVA', error = FALSE)) {
      stop("Please install GSVA!")
    }
    score <- GSVA::gsva(expr = expr,
                  gset.idx.list = geneset,
                  method = gsva.method,
                  verbose = verbose)
    score <- data.frame(t(score), check.names = F)
  }
  if (method == "UCell") {
    if (!PackageCheck('UCell', error = FALSE)) {
      stop("Please install UCell!")
    }
    temp <- CreateSeuratObject(counts = expr)
    n1 <- ncol(temp@meta.data)
    temp <- suppressWarnings(UCell::AddModuleScore_UCell(temp,
                                                         features = geneset,
                                                         name = NULL))
    n2 <- ncol(temp@meta.data)
    score <- temp@meta.data[, (n1+1) : n2, drop = F]
  }
  if (scale){
    score <- apply(score, 2, scale)
    score <- data.frame(score, check.names = F)
  }
  return(score)
}

#' Gene set types in MSigDB.
#' @return A character vector represents types of gene sets in "MSigDB" database.
#' @export
#'
#' @examples
#' msigdb_types()
msigdb_types <- function(){
  msigdb_type <- c(HALLMARK = "H",
                   KEGG = "C2",
                   REATOME = "C2",
                   MOTIF = "C3",
                   BIOCARTA = "C2",
                   "GO:BP" = "C5",
                   "GO:CC" = "C5",
                   "GO:MF" = "C5",
                   HPO = "C5",
                   "TFT:GTRD" = "C3",
                   "TFT:TFT_Legacy" = "C3")
  return(msigdb_type)
}
#' Getting gene sets and accessing enrichment scores.
#' @description SpaTME enables exploration of signatures or gene sets collected
#'   from literatures and databases.
#' @param expr Expression data. Rows are genes and Columns are spots.
#' @param stf STFeature object to store the enrichment scores.
#' @param category Category of gene sets. it should be one of "CuratedSig" and
#'   "MSigDB".
#' @param types Types of gene sets. If \code{category} is "CuratedSig", it
#'   should be one of "Immune" and "Tumor". If \code{category} is "MSigDB", the
#'   types of gene sets can be checked with msigdb_types function.
#' @param subtype Subtypes of gene sets. For "Immune" type, it should be one of
#'   "Immunoresponse", "PancancerCellAtlas", and "TLS; for "Tumor" type, it
#'   should be one of "CancerMPs", "RecurTumorModules", "CancerSEA", and "Fges".
#' @param species Species name. Default is Mus musculus.
#' @param method Method for enrichment analysis, passed on to \code{GsetScore}
#' @param scale Logical value. If TRUE, scaling and centering the scores for all
#'   spots.
#' @param verbose Logical value. If TRUE, print messages.
#' @return If \code{stf} is NULL, return a data with enrichment scores; else,
#'   return stf with \code{GsetSig} slot added.
#' @export
#' @importFrom msigdbr msigdbr
#' @examples
GetGsetSigScore <- function(expr,
                            stf = NULL,
                            category = "CuratedSig",
                            types = c("Immune", "Tumor"),
                            subtype = NULL,
                            species = "Homo sapiens",
                            method = "AddModuleScore",
                            scale = F,
                            verbose = T){
  if (! category %in% c("CuratedSig", "MSigDB")) {
    stop("Category must be CuratedSig or MSigDB!")
  }
  if (!is.null(stf)) {
    gset_score <- stf@GsetSig[[category]]
  } else {
    gset_score <- list()
  }
  if (category == "CuratedSig") {
    data(CuratedSig.lt)
    if (! all(types %in% c("Immune", "Tumor"))) {
      stop("types for CuratedSig are only Immune and Tumor")
    }
    for (fea in types){
      if (verbose) {
        message(paste0("Calculating ",fea , "-related signatures"))
      }
      subtypes <- subtype %||% names(CuratedSig.lt[[fea]])
      temp <- setdiff(subtypes, names(CuratedSig.lt[[fea]]))
      if (length(temp)) {
        warning("subtypes ", paste(temp, collapse = ", "), "are not available for ", fea)
        subtypes <- intersect(subtypes, names(CuratedSig.lt[[fea]]))
      }
      for (sub in subtypes){
        gset <- CuratedSig.lt[[fea]][[sub]]
        gset_score[[fea]][[sub]] <- GsetScore(expr = expr,
                                                  geneset = gset,
                                                  method = method,
                                                  scale = scale,
                                                  verbose = verbose)
      }
    }
  }
  if (category == "MSigDB") {
    msigdb_type <- msigdb_types()
    if (! all(types %in% names(msigdb_type))) {
      stop("Please provide correct subtypes for MSigDB! Check it using msigdb_types()")
    }
    for(type in types){
      if (type %in% c("BIOCARTA", "KEGG", "REACTOME")) {
        type <- paste0("CP:", type)
      }
      if (verbose) {
        message(paste0("Calculating ",type , " geneset"))
      }
      if(type == "HALLMARK"){
        db <- msigdbr(species = species,
                      category = msigdb_type[type])
      } else{
        db <- msigdbr(species = species,
                      category = msigdb_type[type],
                      subcategory = types)
      }
      db <- split(x = db$gene_symbol,
                  f = db$gs_name)
      gset_score[[type]]  <- GsetScore(expr = expr,
                                       geneset = db,
                                       method = method,
                                       scale = scale,
                                       verbose = verbose)

    }
  }
  if (! is.null(stf)) {
    stf@GsetSig[[category]] <- gset_score
    return(stf)
  } else {
    return(gset_score)
  }
}

