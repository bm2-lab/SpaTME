
#' Coordinate check for ST data
#'
#' @param st Seurat object for ST data. The coordinate information will be
#'   checked using the \code{metadata}, \code{images}, or spot names.
#' @param reset Logical indicating whether the existing 'x' and 'y' columns in
#'   \code{metadata} should be reset.
#' @param platform Spatial transcriptomic platform. Specify 'Visium' for hex
#'   lattice geometry or 'ST' for square lattice geometry.
#' @param hexagon.correct Logical value, indicating whether to correct the
#'   coordinate position so that it conforms to the hexagon. Used only when
#'   \code{platform == "visium"}.
#' @param hexagon.direct Direction of hexagon. Available options are:
#' \itemize{
#'  \item{"horizontal"} : spots to left and right, two above, two below.
#'  \item{"vertical"} : spots to top and bottom , two left and two right.
#'  }
#' @param verbose Logical value. If TRUE, print messages.
#'
#' @return Seurat object of ST data with checked coordinates in \code{metadata} named
#'   as x' and 'y'.
#' @export
STcoordCheck <- function(st,
                         reset = F,
                         platform = "Visium",
                         hexagon.correct = F,
                         hexagon.direct = c("horizontal", "vertical"),
                         verbose = T){
  # hexagon.direct:
  # horizontal, spots to left and right, two above, two below
  # vertical, spots to top and bottom , two left and two right
  if (reset) {
    st@meta.data <- st@meta.data[, setdiff(colnames(st@meta.data), c("x", "y"))]
  }
  if (all(c("x", "y") %in% colnames(st@meta.data))){
    if (verbose){
      message("Coordinate is ok!")
    }
  }else {
    if (length(st@images) !=0){
      pos <- st@images[[1]]@coordinates[, c("row", "col")]
      if (verbose){
        message("Get coordinate information from @image")
      }
    }else{
      # check if the spot names contain the coordinate information
      if (grepl("x", colnames(st)[1])){
        pos <- sapply(colnames(st),function(x)strsplit(x, "x", fixed = T))
        pos <- data.frame(pos, check.names = F)
        pos <- t(pos)
        if (verbose){
          message("Get coordinate information from spot names")
        }
      }else {
        stop("Please provide coordinate information(x, y) in meta.data !")
      }
    }
    colnames(pos) <- c("x", "y")
    pos <- apply(pos, 2, as.numeric)
    rownames(pos) <- colnames(st)
    st <- AddMetaData(st,
                      metadata = data.frame(pos))
  }
  if (platform == "Visium"){
    pos <- st@meta.data[, c("x", "y")]
    if (hexagon.correct) {
      library(dplyr)
      if (verbose){
        message("Converting coordinate position to hexagon arrangement!")
      }
      if (hexagon.direct == "horizontal"){
        st@meta.data[, "y"] <- st@meta.data[, "y"]*(sqrt(3)) %>%
          round(2)
      } else if (hexagon.direct == "vertical"){
        st@meta.data[, "x"] <- st@meta.data[, "x"]*(sqrt(3)) %>%
          round(2)
      } else {
        stop("Please provide hexagon.direct value for Visium : horizontal or vertical")
      }
    } else if ( ! all(pos[, "x"] %% 1 ==0 & pos[, "y"] %% 1 ==0)) {
      warning("Coordinates do not fit hexagon shape ! It is recommended to run with hexagon.direct = True")
    }
  }
  print("Completing coordinate check!")
  return(st)
}

#' Title Preprocessing of ST data with SpaTME
#'
#' @param se A Seurat object of ST data to preprocess.
#' @param assay name of assay to use.
#' @param mt.qc Spots with percentage of reads that map to the mitochondrial
#'   genome more than this value will be filtered. If NULL, no filtering.
#' @param norm.SCT Logical value. If TRUE, normalize data with SCTransform.
#'   If FALSE, running NormalizeData, FindVariableFeatures, and ScaleData
#'   Sequentially.
#' @param variable.features.n Number of features to select as top variable features.
#' @param dims Dimensions of reduction to use.
#' @param cluster.resolution Value of the resolution parameter.
#' @param verbose Logical value. If TRUE, print messages.
#' @importFrom Seurat PercentageFeatureSet SCTransform NormalizeData
#'   FindVariableFeatures ScaleData RunPCA FindNeighbors FindClusters RunTSNE
#'   RunUMAP
#' @return A Seurat object with preprocessing results (like PCA, clustering, and UMAP).
#' @export
SePreprocess <- function(se,
                         assay = "RNA",
                         mt.qc = 30,
                         norm.SCT = TRUE,
                         variable.features.n = 3000,
                         dims = 1:20,
                         cluster.resolution = 0.8,
                         verbose = F
){
  if (mt.qc){
    se[["percent.mt"]] = PercentageFeatureSet(se, assay = assay, pattern = "^MT-")
    se <- subset(x = se,
                 subset = percent.mt < mt.qc)
  }
  # Normalization
  if (norm.SCT){
    se <- SCTransform(object = se,
                      assay = assay,
                      variable.features.n = variable.features.n,
                      verbose = verbose)
  }else{
    se <- NormalizeData(object = se,
                      normalization.method = norm.method,
                      verbose = verbose)
    se <- FindVariableFeatures(object = se,
                               nfeatures = variable.features.n,
                               verbose = verbose)
    se <- ScaleData(object = se,
                    verbose = verbose)
  }
  # dimensional reduction and clustering
  se <- RunPCA(object = se,
               verbose = verbose)
  se <- FindNeighbors(object = se,
                      dims = dims,
                      verbose = verbose)
  se <- FindClusters(object = se,
                     resolution = cluster.resolution,
                     verbose = verbose)
  se <- RunTSNE(object = se,
                dims = dims)
  se <- RunUMAP(object = se,
                dims = dims,
                verbose = verbose)
}

