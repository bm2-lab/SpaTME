#' WGCNA analysis
#'
#' @param data Expression data. Rows are selected genes, and columns are spots.
#' @param outdir Directory to save output files.
#' @param power.seq A vector of soft thresholding powers for which the scale
#'   free topology fit indices are to be calculated.
#' @param power soft-thresholding power for network construction. If NULL, it
#'   will be selected from \code{power.seq}.
#' @param verbose integer level of verbosity. Zero means silent, higher values
#'   make the output progressively more and more verbose. Default is 3.
#' @param ... Other arguments passed to \code{blockwiseModules} function in \pkg{WGCNA}.
#' @return A list containing multiply components.
#' @importFrom grDevices pdf dev.off
#' @importFrom SeuratObject PackageCheck
#' @export
#' @examples
BulidCGnet <- function(data,
                       outdir = NULL,
                       power.seq = NULL,
                       power = NULL,
                       verbose = 3,
                       ...){
  # check package.
  if (!PackageCheck('WGCNA', error = FALSE)) {
    stop("Please install WGCNA!")
  }
  if (is.null(outdir)) {
    message("The output directory is recommended to save the result files!")
  }
  # check expression
  gsg = WGCNA::goodSamplesGenes(data, verbose = verbose)
  if (!gsg$allOK) {
    if (verbose) {
      if (sum(!gsg$goodGenes)>0)
        printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
      if (sum(!gsg$goodSamples)>0)
        printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
    }
    data = data[gsg$goodSamples, gsg$goodGenes]
  }
  # select soft threshold
  if (is.null(power.seq)) {
    power.seq <- c(c(1:10), seq(from = 12, to=20, by=2))
  }
  message("Selecting soft threshold...")
  sft <- WGCNA::pickSoftThreshold(data,
                                  powerVector = power.seq,
                                  verbose = verbose)
  # Plot the results:
  # Scale-free topology fit index as a function of the soft-thresholding power
  if (is.null(power)) {
    power <- sft$powerEstimate
    if (is.na(power)) {
      warning("No soft-thresholding power selected!")
    } else {
      message(paste0("soft-thresholding power is selected: ", power))
    }
  }
  if ( ! is.null(outdir)){
    pdf(paste0(outdir, "/pickSoftThreshold.pdf"))
    cex1 = 0.9
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
         main = paste("Scale independence"))
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels=power,cex=cex1,col="#c32f27")
    # this line corresponds to using an R^2 cut-off of h
    abline(h=0.90,col="#c32f27")
    # Mean connectivity as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab="Soft Threshold (power)",ylab="Mean Connectivity",
         main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=power, cex=cex1,col="#c32f27")
    dev.off()
  }
  # One-step network construction and module detection
  message("Starting module detection...")
  cor <- WGCNA::cor
  net = WGCNA::blockwiseModules(datExpr = data,
                                power = power,
                                TOMType = "unsigned",
                                numericLabels = TRUE,
                                verbose = verbose,
                                ...)
  if ( ! is.null(outdir)){
    saveRDS(net, paste0(outdir, "/net.rds"))
  }
  moduleColors <- LabelMapcolor(labels = as.character(net$colors),
                                assgin.col = c("0" = "#e9ecef"))
  plot <- WGCNA::plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                             "Module colors",
                             dendroLabels = FALSE, hang = 0.03,
                             addGuide = TRUE, guideHang = 0.05)
  if ( ! is.null(outdir)){
    # plot
    pdf(paste0(outdir, "/Cluster_Dendrogram.pdf"))
    print(plot)
    dev.off()
  } else {
    print(plot)
  }
  cor<-stats::cor
  return(net)
}
