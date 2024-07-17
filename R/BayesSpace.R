#' Title Domain annotation with BayesSpace
#'
#' @param st Seurat object for ST data.
#' @param assay Name of assay.
#' @param slot Name of slot.
#' @param platfrom Spatial transcriptomic platform. Specify 'Visium' for hex
#'   lattice geometry or 'ST' for square lattice geometry.
#' @param num_clu Number of clusters. If NULL, it will be selected
#'   automatically.
#' @param qs A vector of continuous integers. Given the range to select the
#'   number of cluster. Only used when \code{num.clu = NULL}. Default 3-7.
#' @param show_qplot Logical value. If TURE, plots the pseudo-log-likelihood
#' to help choose num_clu.
#' @param verbose Logical value. If TURE, print messages.
#' @param seed Set the seed.
#' @param ... Other arguments passed on to \emph{spatialCluster}.
#' @return Returns a modified \code{sce} with cluster assignments stored in
#'   \code{colData} under the name \code{spatial.cluster}.
#' @importFrom grDevices pdf dev.off
#' @importFrom SeuratObject PackageCheck
#' @export
BayesCluster <- function(st,
                         assay = "Spatial",
                         slot = "count",
                         platfrom = "Visium",
                         num_clu = NULL,
                         qs = seq(3, 7),
                         show_qplot = T,
                         verbose = T,
                         seed = 123,
                         ...
){
  # check package.
  if (!PackageCheck('BayesSpace', error = FALSE)) {
    stop("Please install BayesSpace!")
  }
  if (!PackageCheck('SingleCellExperiment', error = FALSE)) {
    stop("Please install SingleCellExperiment!")
  }
  coldata <- st@meta.data
  if ( ! all(c("row", "col") %in% colnames(coldata))) {
    if (all(c("x", "y") %in% colnames(coldata))) {
      colnames(coldata)[colnames(coldata) == "x"] = "row"
      colnames(coldata)[colnames(coldata) == "y"] = "col"
    } else {
      stop("Spot coordinates must be provided as columns named row(x) and col(y) in colData")
    }
  }
  if (verbose) {
    print("Starting clustering with BayesSpace")
  }
  data <- GetAssayData(object = st,
                       assay = assay,
                       slot = slot)
  data <- as.matrix(data)
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = data),
                              colData = coldata)
  set.seed(seed)
  sce <- BayesSpace::spatialPreprocess(sce, platform = platfrom)
  if (is.null(num_clu)){
    if (verbose) {
      message("Number of cluster is not provided, qTune will used for selecting best q value")
    }
    set.seed(seed)
    sce <- BayesSpace::qTune(sce = sce,
                             qs = qs,
                             platform = platfrom)
    if (show_qplot) {
      print(BayesSpace::qPlot(sce))
    }
    # select best q
    logliks <- attr(sce, "q.logliks")
    logliks$loglik <- Norm01(-(logliks$loglik))
    cross <- sapply(1:(nrow(logliks) - 1) , function(i) {
      logliks[i, "loglik"] - logliks[i + 1, "loglik"]
      })
    if (any(cross < 0)) {
      names(cross) <- logliks$q[1:(nrow(logliks) - 1)]
      cross <- cross[cross < 0]
      num_clu <- min(as.numeric(names(cross)))
    } else {
      cross <- sapply(1:(length(cross) - 1), function(i) {
        1 - cross[i+1]/cross[i]
      })
      cross <- unlist(cross)
      num_clu <- logliks$q[which.max(cross) + 1]
    }
    if (verbose) {
      print(paste0("The recommended cluster number is ", num_clu))
    }
  }
  set.seed(seed)
  sce <- BayesSpace::spatialCluster(sce = sce,
                                    q = num_clu,
                                    platform = platfrom,
                                    ...)
  return(sce)
  }
