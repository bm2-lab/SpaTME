#' Calculation of effective expression of ligand and receptor proteins.
#'
#' @param Expr Expression data. Rows are genes and columns are spots.
#' @param LR_input A data.frame containing information of ligand and receptor
#'   proteins.
#' @param st_pos A data.frame with coordinate information in two columns named
#'   'x' and 'y'.
#' @param do.dist Logical value. If TURE, the expression of secreted proteins
#'   will re-calculated using neighboring spots based on distance.
#' @param r.dist A numeric value. Distance to define the neighboring spots.
#'   Default is 4.
#' @param long.dist.method Methods to integrate expression from spot and it's
#'   neighbors. It should be one of 'weighted' and 'mean'.
#' @param adjust A numeric value controlling the weights of neighboring spots
#'   with difference distance. Default is 2.
#' @param na.rm Logical value. If TURE, the proteins with no gene in \code{Expr}
#'   will be removed.
#' @param verbose Logical value. If TURE, print messages.
#' @return A data.frame of expression data. Rows are ligand and receptor
#'   proteins, and columns are spots.
#' @importFrom foreach foreach %do%
#' @export

BuildLRExprssion <- function(Expr,
                             LR_input,
                             st_pos = NULL,
                             do.dist = T,
                             r.dist = 4,
                             long.dist.method = "weighted",
                             adjust = 2,
                             na.rm = T,
                             verbose = T){
  # do.dist : whether to calculate weighted expression of secreted proteins based on distance.
  if (do.dist) {
    if (is.null(st_pos)) {
      stop("Position must be provided!")
    }
  }
  if (verbose) {
    message("Calculating the ligand/receptor protein expression value.")
    pb <- txtProgressBar(min = 0,
                         max = nrow(LR_input),
                         style = 3)
  }
  na_lr <- c()
  LR_expression <- foreach(i =  1:nrow(LR_input), .combine = "rbind") %do% {
    comm.gene <- intersect(x = rownames(Expr),
                           y = unlist(LR_input[i, c("gene_1", "gene_2", "gene_3")]))
    if (length(comm.gene) > 0) {
      LR.expr <- Expr[comm.gene , , drop = F]
      LR.expr <- apply(LR.expr, 2, min)
    } else {
      LR.expr <- NA
      na_lr <- c(na_lr, rownames(LR_input)[i])
    }
    if (verbose) {
      setTxtProgressBar(pb = pb, value = i)
    }
    return(LR.expr)
  }
  if (verbose) {
    close(con = pb)
  }
  if (length(na_lr) != 0 ) {
    if (verbose) {
      message("There are ", nrow(LR_input) - length(na_lr), " L/R proteins with detected genes")
    }
  }
  rownames(LR_expression) <- rownames(LR_input)
  if (na.rm) {
    LR_expression <- LR_expression[setdiff(rownames(LR_expression), na_lr),]
  }
  #
  if (do.dist) {
    temp_LR_expr <- LR_expression[LR_input[rownames(LR_expression), "protein_type"] != "only_menbrane", ]
    if (verbose) {
      message("Calculating weighted expression of ", nrow(temp_LR_expr), "secreted proteins according neighbors ...")
    }
    nei.pairs <- FindNeispots(st_pos = st_pos,
                              r.dist = r.dist,
                              return.list = F)# calculating Euclidean Distance
    dist <- nei.pairs[, "dist"]
    names(dist) <- nei.pairs[, "spot2"]
    reg.lt <- split(x = as.factor(dist),
                    f = nei.pairs[, "spot1"])
    if (long.dist.method == "weighted") {
      d <- unique(unlist(reg.lt))
      d <- as.numeric(paste(d))
      weights <- KernelWeight(d = d,
                              adjust = adjust)
      weights <- weights[levels(reg.lt[[1]])]
    }
    temp_LR_expr <- UnitEnrich(data = t(temp_LR_expr),
                               d.lt = reg.lt,
                               method = long.dist.method,
                               weight = weights,
                               verbose = verbose
    )
    temp_LR_expr <- round(temp_LR_expr, digits = 3)
    temp_LR_expr <- t(temp_LR_expr)
    LR_expression[rownames(temp_LR_expr),] <- temp_LR_expr
  }
  LR_expression <- data.frame(LR_expression, check.names = F)
  return(LR_expression)
}

#' Calculation of ligand-receptor interaction (LRI) scores.
#'
#' @param interaction_input A data.frame containing information of
#'   ligand-receptor interactions.
#' @param LR_expression A data.frame containing expression data of ligand and
#'   receptor proteins.
#' @param na.rm Logical value. If TURE, LRIs with with all NA values will be
#'   removed.
#' @param p The number of CPU cores used for parallel computing.
#' @param verbose Logical value. If TURE, print messages.
#'
#' @return A data.frame of LRI scores. Rows are LRIs, and columns are spots.
#' @export
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom doParallel stopImplicitCluster
#' @importFrom parallel makeCluster stopCluster
#' @importFrom SeuratObject PackageCheck
CalLRIScore <- function(interaction_input,
                       LR_expression,
                       na.rm = T,
                       p = NULL,
                       verbose = T){
  if (! is.null(p)){
    if ( !PackageCheck("doSNOW")){
      stop("Please install package: doSNOW or set verbose to FALSE")
    } else {
      suppressMessages({cl = makeCluster(p)
      doSNOW::registerDoSNOW(cl)})
      `%myinfix%` <- `%dopar%`
    }
  } else {
    `%myinfix%` <- `%do%`
  }
  opts <- NULL
  if (verbose) {
    message("Starting to calculate scores for ", nrow(LR_expression)," LR interaction...")
    pb <- txtProgressBar(min = 0,
                         max = nrow(LR_expression),
                         style = 3)
    if ( !is.null(p)){
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
    }
  }
  #
  na_lr1 <- 0
  system.time({
  LR_score <- foreach(i = 1:nrow(interaction_input),
                           .combine = "rbind",
                           .options.snow = opts) %myinfix% {
                             lr_a <- interaction_input[i, "partner_a"]
                             lr_b <- interaction_input[i, "partner_b"]
                             if (all(c(lr_a, lr_b) %in% rownames(LR_expression))) {
                               score <- LR_expression[lr_a, ] * LR_expression[lr_b, ]
                             } else{
                               score <- NA
                               na_lr1 <- na_lr1 + 1
                             }
                             if (verbose & is.null(p)) {
                               setTxtProgressBar(pb = pb, value = i)
                             }
                             return(score)
  }
  })
  if (verbose) {
    close(con = pb)
  }
  rownames(LR_score) <- rownames(interaction_input)
  colnames(LR_score) <- colnames(LR_expression)

  if (na.rm) {
    LR_score <- LR_score[rowSums(! is.na(LR_score)) > 0,]
  }
  if (! is.null(p)){
    stopImplicitCluster()
    stopCluster(cl)
  }
  return(data.frame(LR_score, check.names = F))
}

