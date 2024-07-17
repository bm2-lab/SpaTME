# Internal function to convert count to TPM.
Count2TPM <- function(counts,geneLen) {
  common.genes <- intersect(rownames(counts),rownames(geneLen))
  x <- counts[common.genes,]/geneLen[common.genes,]
  return(t(t(x)*1e6/colSums(x)))
}

#' Preparing input data for NMF analysis.
#'
#' @param data Expression data. Rows are genes and columns are samples.
#' @param totpm Logical value. If TURE, convert the counts to TPM values.
#' @param genelen A data.frame containing the gene length information. Only used when \code{totpm = TURE}.
#' @param dolog Logical value. If TURE, conduct log2 transformation.
#' @param pt_gene_exp Only keep genes that are detected in a minimum fraction (pt_gene_exp) of samples.
#' @param base_gene_exp Only genes with expression above this value are considered to be expressed.
#' @param pseudogene.filt Logical value indicating whether to filt pseudogene.
#'
#' @return A data.frame with pre-processed expression profiles.
#' @export
#' @examples
NMF_bulk_input <- function(data,
                           totpm = F,
                           genelen = NULL,
                           dolog = F,
                           pt_gene_exp = 0.1,
                           base_gene_exp = 0,
                           pseudogene.filt = F){
  data <- data[! is.na(rowSums(data)), ]
  if(totpm){
    print("starting convert counts to TPM")
    data <- Count2TPM(data, geneLen = genelen)
  }
  if(dolog){
    data <- log2(data + 1)
  }
  g.pt <- rowSums(data > base_gene_exp)/ncol(data) > pt_gene_exp
  data <- data[g.pt,]
  if(pseudogene.filt){
    data(pseudogenes)
    data <- data[! rownames(data) %in% pseudogenes,]
  }
  return(data)
}

#' NMF analysis to identify phenotype-related factors
#'
#' @param expr Expression data. Rows are phenotype-related genes and columns are
#'   samples.
#' @param rank A single numeric value or a numeric vector. Number of factors.
#' @param nrun Number of runs to perform. Default is 30.
#' @param seed a single numeric that is used to seed the random number
#'   generator, before generating a random starting point. Default is 123456.
#' @param min_cophenetic Minimum cophenetic coefficient. If the rank is a numeric
#'   vector, the cophenetic coefficient is used to measure the stability of the
#'   clusters obtained. The selected number of ranks should be above this value.
#'   Defualt is 0.95.
#' @param name_prefix A character specified the prefix name of factors.
#' @param verbose Logical value. If TURE, print messages.
#' @param return.all Logical value
#'
#' @return If the parameter return.all is TURE and the rank is a numeric vector,
#'   return a list containing all the NMFfit objects of each rank and their
#'   cophenetic coefficients; else, return an object of class NMFfit.
#' @export
#' @import NMF
#' @examples
RunNMFtest <- function(expr,
                     rank,
                     nrun = 30,
                     seed = 123456,
                     min_cophenetic = 0.95,
                     name_prefix = "Factor",
                     return.all = F,
                     verbose = T){
  # check gene or samples with all 0
  which.r <- which(apply(expr, 1, sd) == 0)
  if (length(which.r) != 0) {
    expr <- expr[-c(which.r),]
    warning(paste0("Features: ", paste0(rownames(expr)[which.r], collapse = "/"), " with no variance are removed!"))
  }
  expr <- t(scale(t(expr)))
  expr[expr < 0] <- 0
  which.c <- which(apply(expr, 2, sd) == 0)
  if (length(which.c) != 0) {
    expr <- expr[, -c(which.c)]
    warning(paste0("Samples: ", paste0(colnames(expr)[which.c], collapse = "/"), " with no variance are removed!"))
  }
  if (verbose) {
    message("starting nmf analysis...")
    print(Sys.time())
  }
  nmf <- nmf(x = expr,
             rank = rank,
             nrun = nrun,
             seed = seed,
             .options = "v")
  if (length(rank) > 1) {
    for (i in 1:length(nmf$fit)) {
      name <- paste0(name_prefix, "_", 1:ncol(nmf$fit[[i]]@fit@W))
      colnames(nmf$fit[[i]]@fit@W) <- name
      rownames(nmf$fit[[i]]@fit@H) <- name
    }
  } else {
    name <- paste0(name_prefix, "_", 1:ncol(nmf@fit@W))
    colnames(nmf@fit@W) <- name
    rownames(nmf@fit@H) <- name
  }
  if ( length(rank) > 1) {
    if (verbose) {
      message("Best rank will be seleted from ", rank)
    }
    # rank selection
    cophs <- data.frame(n_clusters = names(nmf$consensus),
                        Cophenetic = nmf$measures$cophenetic)
    best_rank <- biggest_drop(x = cophs,
                              min_cophenetic = min_cophenetic)
    if (verbose) {
      message("Selected rank is :", best_rank)
    }
    if (return.all) {
      res <- list(nmf = nmf,
                  cophs = cophs,
                  best_rank = best_rank)
    } else {
      res <- nmf$fit[[best_rank]]
    }
  } else {
    res <- nmf
  }
  return(res)
}

biggest_drop <- function(x, min_cophenetic = 0.95)
{
  cross = c(sapply(1:(nrow(x) - 1) , function(i) round(x[i,"Cophenetic"], 2) - round(x[i + 1, "Cophenetic"], 2)), NA)
  names(cross) = x$n_clusters
  cross = cross[x[, "Cophenetic"] >= min_cophenetic]
  cross <- na.omit(cross)
  if(length(cross) ==0){
    best.clu <- NA
    print("There is no best cluster number can be seleted with the given cutoff of min_cophenetic !")
  }
  else{
    max.drop = max(cross)
    best.clu <- max(names(cross[cross == max.drop]))
  }
  return(best.clu)
}

#' Confirming phenotype-related NMF factors
#' @param nmfobj An object of class NMFfit.
#' @param phenotype A numeric or character vector representing the annotation of spots or
#' phenotype of samples.
#' @param ... Other arguments passed to \code{PhenoAssoFeatures}.
#'
#' @return A character vector represent association between phenotype and factors.
#' @export
#'
#' @examples
HPhenoAsso <- function(nmfobj, phenotype, method, ...){
  train_H <- nmfobj@fit@H
  asso_res <- PhenoAssoFeatures(data = train_H,
                                phenotype = phenotype,
                                method = method,
                                ...)
  if(method == "cox"){
    p.sig <- asso_res[, "coef.p"]
    test.sig <- asso_res[, "cox.test.p"]
    W_type <- ifelse(p.sig <= 0.05 & test.sig >= 0.1,
                     ifelse(p.sig >= 0.01, "*",
                            ifelse(p.sig >= 0.001, "**", "***")), "Nonsig")
    W_type <- paste0(ifelse(asso_res[, "coef"] <0, "Pos", "Neg"), "_", W_type)
  }
  if (method == "auc") {
    W_type <- paste0("AUC_", asso_res[, "myAUC"])}
  if(method == "cor"){
    W_type <- paste0("corr", round(asso_res[, "corr"], 2))
  }
  if (method == "glm"){
      p.sig <- asso_res[, "p.value"]
      W_type <- ifelse(p.sig <= 0.05,
                       ifelse(p.sig >= 0.01, "*",
                              ifelse(p.sig >= 0.001, "**", "***")), "Nonsig")
      W_type <- paste0(ifelse(asso_res[, "coef"] >0, "Pos", "Neg"), W_type)
  }
  if (method == "wilcox") {
    p.sig <- asso_res[, "p.value"]
    W_type <- ifelse(p.sig <= 0.05,
                     ifelse(p.sig >= 0.01, "*",
                            ifelse(p.sig >= 0.001, "**", "***")), "Nonsig")
    W_type <- paste0(ifelse(asso_res[, "avg.diff"] >0, "Pos", "Neg"), W_type)
  }
  names(W_type) <- rownames(asso_res)
  return(W_type)
}
#
#' NMF Multiplicative Updates for Kullback-Leibler Divergence
#'
#' @param v target matrix.
#' @param w current basis matrix.
#' @param h current coefficient matrix.
#' @param nbterms number of fixed basis terms.
#' @param ncterms number of fixed coefficient terms.
#' @param copy logical that indicates if the update should be made on the original
#' matrix directly (\code{FALSE}) or on a copy (\code{TRUE} - default).
#' @return a matrix of the same dimension as the input matrix to update.
std.divergence.update.h <- function(v, w, h, nbterms=0L, ncterms=0L, copy=TRUE)
{
  .Call("divergence_update_H", v, w, h, nbterms, ncterms, copy, PACKAGE='NMF')
}
# Internal function to predict H matrix in new data.
#' @references Luca, B.A., et al. Atlas of clinically distinct cell states and
#'   ecosystems across human solid tumors. Cell 184, 5482-5496 e5428 (2021).
#' @import NMF
NMFpredict <- function(W, new_data)
{
  trainig_gene_set = unique(rownames(W))
  new_data = new_data[match(trainig_gene_set, rownames(new_data)),]
  rownames(new_data) = trainig_gene_set
  new_data[is.na(new_data)] = 0
  to_predict = as.matrix(new_data)
  # check samples with non-expressed genes
  final.sam <- colnames(to_predict)
  out <- colnames(to_predict)[colSums(to_predict)==0]
  to_predict <- to_predict[, setdiff(colnames(to_predict), out)]
  my_method <- function (i, v, x, copy = FALSE, eps = .Machine$double.eps, ...)
  {
    w <- .basis(x)
    h <- .coef(x)
    nb <- nbterms(x)
    nc <- ncterms(x)
    h <- std.divergence.update.h(v, w, h, nbterms = nb, ncterms = nc, copy = copy)
    if (i%%10 == 0) {
      h <- pmax.inplace(h, eps, icterms(x))
    }
    if (copy) {
      .coef(x) <- h
    }
    return(x)
  }
  ws = W
  ws <- ws[apply(to_predict, 1, function(x) var(x) > 0),]
  to_predict = to_predict[apply(to_predict, 1, function(x) var(x) > 0),]
  ws = as.matrix(ws)

  dummy = rnmf(ncol(W), to_predict)

  my.seeding.method <- function(model, target){
    basis(model) <- ws
    coef(model) <- dummy@H
    # return updated object
    return(model)
  }
  nmf_method <- NMFStrategy('my-method', 'brunet', Update = my_method, objective = 'KL', Stop='connectivity')
  new_nmf = nmf(to_predict, ncol(W), nrun = 1, method = nmf_method, seed = my.seeding.method, .opt='P1')
  if (length(out) >0) {
    H <- new_nmf@fit@H
    H = as.data.frame(H)
    H[, out] <- 0
    H <- as.matrix(H[, final.sam])
    new_nmf@fit@H <- H
  }
  new_nmf
}
#' Prediction of factors in ST data
#'
#' @param st a seurat object for ST data.
#' @param assay Name of \code{assay} in \code{st}.
#' @param slot Name of \code{slot} in \code{st}.
#' @param W The W matrix in NMF analysis.
#' @param plotfile File name to save this plot. Only supporting pdf files.
#' @param numCol Number of columns in the plot grid.
#' @param verbose Logical value. If TURE, print messages.
#' @param ... Other arguments passed to \code{SpotVisualize}.
#'
#' @return A NMFfit object with predicted new H matrix.
#' @export
#' @importFrom cowplot plot_grid
#' @importFrom grDevices pdf dev.off
#' @importFrom Seurat GetAssayData
#' @examples
PredNMFinST <- function(st,
                        W,
                        assay = "SCT",
                        slot = "data",
                        plotfile = NULL,
                        numCol = 1,
                        verbose = T,
                        ...){
  if (verbose) {
    message("Predict expression of factors in spatial data...")
  }
  newdata <- GetAssayData(object = st,
                          assay = assay,
                          slot = slot)
  newdata <- as.matrix(newdata)
  newdata <- newdata[rowSums(newdata) >0,]
  newdata <- t(scale(t(newdata)))
  newdata[newdata < 0] = 0
  #
  nmf_pred <- NMFpredict(W, newdata)
  if (! is.null(colnames(W))) {
    name <- colnames(W)
  } else {
    name <- paste0("Factor_", 1:ncol(W))
  }
  colnames(nmf_pred@fit@W) <- name
  rownames(nmf_pred@fit@H) <- name
  # plot
  H = nmf_pred@fit@H
  plot.list = c()
  for(i in rownames(H)){
    plot.list[[i]] <- SpotVisualize(st,
                       meta = scale(H[i,])[, 1],
                       return = T,
                       title = i,
                       ...)
    plot = plot_grid(plotlist = plot.list, ncol = numCol)
  }
  if (! is.null(plotfile)) {
    pdf(plotfile)
    print(plot)
    dev.off()
  } else {
    print(plot)
  }
  return(nmf_pred)
}

#' Meta-genes of each factors
#'
#' @param ref_W A data.frame contains the W matrix obtained from NMF analysis.
#' @param top_num Number of top genes.
#'
#' @return A vector of factors named with meta-genes.
#' @export
#'
#' @examples
FactorMetagenes <- function(ref_W,
                            top_num = 200){
  # sd_n :1.96sd: 0.01 pvalue; 1.65sd:0.1 pvalue
  if(is.null(colnames(ref_W))){
    colnames(ref_W) <- paste0("Factor_", 1:ncol(ref_W))
  }
  mg <- apply(ref_W, 2, function(x){
    temp <- order(x, decreasing = T)[1: top_num]
    return(c(rownames(ref_W)[temp]))
  })
  mg_vt <- reshape2::melt(mg)
  mg_vt <- setNames(mg_vt$Var2, mg_vt$value)
  return(mg_vt)
}

#' Enrichment analysis for meta-genes of factors
#'
#' @param mg_vt The vector of meta-genes.
#' @param fun One of "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or
#'   "enrichPathway".
#' @param plot A logical value indicating whether to plot.
#' @param show_num The number of terms to be displayed.
#' @param simplify A logical value indicating whether to simplify output from
#'   compareCluster by removing redundancy of enriched GO terms.
#' @param savefile File name to save this plot (pdf).
#' @param wrap_width Positive integer giving target line width in characters. A
#'   width less than or equal to 1 will put each word on its own line.
#' @param p.width The width of pdf file.
#' @param p.height The height of pdf file.
#' @param ... Other arguments passed on to \code{compareCluster} in \pkg{clusterProfiler}.
#' @importFrom ggplot2 ggplot scale_y_discrete theme
#' @return A compareClusterResult object.
#' @export
#'
#' @examples
FactorEnrichAnalysis <- function(mg_vt,
                                 fun = "enrichGO",
                                 plot = T,
                                 show_num = 5,
                                 simplify = T,
                                 savefile = NULL,
                                 wrap_width = 50,
                                 p.width = 8, p.height = 10,
                                 ...){
  # mg_vt; group vector named by metagene
  metagenes <- data.frame(factor = mg_vt,
                          genes = names(mg_vt))
  # enrichment analysis
  # symbol convert to entreid
  gid <- clusterProfiler::bitr(unique(metagenes$genes), 'SYMBOL', 'ENTREZID', OrgDb='org.Hs.eg.db',)
  metagenes <- full_join(metagenes, gid, by=c('genes' = 'SYMBOL'))
  metagenes <- metagenes[! is.na(metagenes$ENTREZID), ]
  x = clusterProfiler::compareCluster(ENTREZID ~ factor, data = metagenes, fun = fun, ...)
  if (simplify & fun == "enrichGO") {
    x <- x %>% clusterProfiler::simplify()
  }
  if(plot){
    p <- enrichplot::dotplot(x, showCategory = show_num) +
      theme(axis.text.x = element_text(angle=45, hjust=1),
            axis.text = element_text(size = 18)) +
      scale_y_discrete(labels=function(x) stringr::str_wrap(x, width= wrap_width))
    print(p)
  }
  if(! is.null(savefile)){
    pdf(savefile, width = p.width, height = p.height)
    print(p)
    dev.off()
  }
  return(x)
}

