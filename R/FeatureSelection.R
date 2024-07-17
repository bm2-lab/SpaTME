###############################################################
# Internal function
###############################################################
# Cox test
#' @importFrom survival coxph Surv cox.zph
#' @importFrom progressr progressor
Coxtest <- function(data,
                    phenotype,
                    set_progress,
                    my_sapply){
  # phenotype, Should be a two-column matrix with columns named 'time' and 'status'.
  # The latter is a binary variable, with '1' indicating event (e.g.recurrence of cancer or death), and '0' indicating right censored.
  cox.genes <- set_progress({
    p <- progressor(along = 1:nrow(data))
    my_sapply(X = 1:nrow(data),
              FUN = function(i) {
                p(sprintf("num = %g", i))
                x <- unlist(data[i, ])
                cox = suppressWarnings(coxph(Surv(phenotype$time, phenotype$status) ~ x))
                coef = summary(cox)$coefficients[c(1,5)]
                c.index = cox$concordance["concordance"]
                #
                cox.test = cox.zph(cox)
                cox.test.p = cox.test$table["GLOBAL","p"]
                return(c(coef, c.index, cox.test.p))
              }
              )
  })
  cox.genes <- t(cox.genes)
  rownames(cox.genes) <- rownames(data)
  colnames(cox.genes) <- c("coef", "coef.p", "c.index", "cox.test.p")
  cox.genes <- data.frame(cox.genes)
  return(cox.genes)
}

# Regression models
#' @importFrom progressr progressor
Glmtest <- function(data,
                    phenotype,
                    family,
                    set_progress,
                    my_sapply){
  # phenotype:
  # Continuous dependent variable. Should be a quantitative vector for family = gaussian.
  # Binary group indicator vector. Should be either a 0-1 encoded vector or a factor with two levels for family = binomial.
  y <- phenotype
  asso.fea <- set_progress({
    p <- progressor(along = 1:nrow(data))
    my_sapply(X = 1:nrow(data),
              FUN = function(i){
                p(sprintf("num = %g", i))
                x <- unlist(data[i, ])
                model <- glm(formula = y ~ x, family = family)
                coef <- summary(model)
                coef <- coef$coefficients["x", ]
              })
  })
  asso.fea <- t(asso.fea)
  rownames(asso.fea) <- rownames(data)
  colnames(asso.fea) <- c("coef",  "Std.Error", "z.value",  "p.value")
  asso.fea <- data.frame(asso.fea, check.names = F)
  return(asso.fea)
}

# Wilcoxon Rank Sum test
#' @importFrom progressr progressor
WilcoxTest <- function(
    data,
    phenotype,
    set_progress,
    my_sapply
    ) {
  if (class(phenotype) != "factor") {
    phenotype <- factor(phenotype)
  }
  group <- levels(phenotype)
  res <- set_progress({
    p <- progressor(along = 1:nrow(data))
    my_sapply(X = 1:nrow(x = data),
              FUN = function(i) {
                p(sprintf("num = %g", i))
                y <- unlist(data[i, ])
                pval <- wilcox.test(y ~ phenotype,)$p.value
                avg.diff <- ExpMean(
                  x = as.numeric(x = data[i, phenotype == group[2]])
                ) - ExpMean(
                  x = as.numeric(x = data[i, phenotype == group[1]])
                )
                return(c(pval, avg.diff))
    })
  })
  res <- t(res)
  colnames(res) <- c("p.value", "avg.diff")
  return(data.frame(res, row.names = rownames(x = data), check.names = F))
}

ExpMean <- function(x, ...) {
  if (inherits(x = x, what = 'AnyMatrix')) {
    return(apply(X = x, FUN = function(i) {log(x = mean(x = exp(x = i) - 1) + 1)}, MARGIN = 1))
  } else {
    return(log(x = mean(x = exp(x = x) - 1) + 1))
  }
}

#  AUC metric
#' @importFrom progressr progressor
AUCMarkerTest <- function(data1,
                          data2,
                          set_progress,
                          my_sapply,
                          verbose = T) {
  myAUC <- set_progress({
    p <- progressor(along = 1:nrow(data1))
    my_sapply(
    X = 1:nrow(data1),
    FUN = function(i) {
      p(sprintf("num = %g", i))
      res <- DifferentialAUC(
        x = as.numeric(x = data1[i, ]),
        y = as.numeric(x = data2[i, ])
      )
      return(res)
    }
  )
  })
  myAUC <- unlist(myAUC)
  myAUC[is.na(x = myAUC)] <- 0
  avg_diff <- unlist(x = lapply(
    X = 1:nrow(data1),
    FUN = function(x) {
      return(
        ExpMean(
          x = as.numeric(x = data2[x, ])
        ) - ExpMean(
          x = as.numeric(x = data1[x, ])
        )
      )
    }
  ))
  toRet <- data.frame(cbind(myAUC, avg_diff), row.names = row.names(data1))
  toRet <- toRet[rev(x = order(toRet$myAUC)), ]
  return(toRet)
}
# AUC metric
#' @importFrom ROCR prediction performance
DifferentialAUC <- function(x, y) {
  prediction.use <- prediction(
    predictions = c(x, y),
    labels = c(rep(x = 0, length(x)), rep(x = 1, length(y))),
    label.ordering = 0:1
  )
  perf.use <- performance(prediction.obj = prediction.use, measure = "auc")
  auc.use <- round(x = perf.use@y.values[[1]], digits = 3)
  return(auc.use)
}

MarkerTest <- function(data,
                       phenotype,
                       auc.power = 0,
                       set_progress,
                       my_sapply,
                       verbose = T
) {
  if (class(phenotype) != "factor") {
    phenotype <- factor(phenotype)
  }
  group <- levels(phenotype)
  to.return <- AUCMarkerTest(
    data1 = data[, phenotype == group[1], drop = FALSE],
    data2 = data[, phenotype == group[2], drop = FALSE],
    set_progress,
    my_sapply,
    verbose = verbose
  )
  to.return$power <- (to.return$myAUC - 0.5) * 2
  to.return <- to.return[abs(to.return$power) >= auc.power, ]
  return(to.return)
}

# Differential expression using DESeq2
#' @importFrom SeuratObject PackageCheck
DESeq2DETest <- function(
    data,
    phenotype,
    contrast,
    verbose = TRUE
) {
  if (!PackageCheck('DESeq2', error = FALSE)) {
    stop("Please install DESeq2!")
  }
  group.info <- data.frame(row.names = colnames(data),
                           group = factor(phenotype))
  group.info$wellKey <- rownames(x = group.info)
  dds1 <- DESeq2::DESeqDataSetFromMatrix(
    countData = data,
    colData = group.info,
    design = ~ group
  )
  dds1 <- DESeq2::estimateSizeFactors(object = dds1)
  dds1 <- DESeq2::estimateDispersions(object = dds1, fitType = "local")
  dds1 <- DESeq2::nbinomWaldTest(object = dds1)
  res <- DESeq2::results(
    object = dds1,
    contrast = contrast,
    alpha = 0.05
  )
  colnames(res)[which(colnames(res)=="pvalue")] <- "p.value"
  return(res)
}
# correlation
CorrTest <- function(data,
                     phenotype,
                     cor.method){
  res <- apply(data, 1, function(x){
    corr <- cor.test(x = x,
             y = phenotype,
             method = cor.method)
    return(c(corr$estimate, corr$p.value))
  })
  res <- t(res)
  colnames(res) <- c("corr", "p.value")
  return(data.frame(res, check.names = F))
}
#
#' Identification of phenotype or annotation associated features.
#' @description
#' To explore the association between features and specific phenotypes of
#' patients or annotations of spots, this function utilizes various methods that
#' can be selected based on the type of the target variable.
#' @param data A data matrix or data frame; col are samples and row are features.
#' @param phenotype A numeric or character vector representing the annotation of spots or
#' phenotype of samples.
#' @param method Method to use. Available options are:
#' \itemize{
#'  \item{"cox"} ：Fits a Cox proportional hazards regression model. This is used
#'  when the phenotype is survival time of patients.
#'  \item{"wilcox"} : Identifies differential features between two
#'  groups of spots or samples using a Wilcoxon Rank Sum test.
#'  \item{"glm"} Fits generalized linear models.
#'  \item{"auc"} A metric obtained from ROC analysis.
#'  \item{"deseq2"} Identifies differential features between two groups using
#'  \pkg{DESeq2}. Only integer values are supported.
#'  \item{"cor"} Test for association/correlation between paired samples/spots. It is used for continuous variables.
#'  }
#' @param family a character string naming a family function when \code{method} is "glm".
#' @param cor.method Method used for correlation coefficient. One of "pearson",
#'   "kendall", or "spearman", can be abbreviated.
#' @param deseq2_contrast Argument specifies the groups for comparison. It is
#'   passed to \code{results} function in \pkg{DESeq2}.
#' @param p.adj Logical value indicating whether to adjust the p values.
#' @param adj.method Method to adjust p values. It should be one of
#'   \code{c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr",
#'   "none")}. Default is "BH".
#' @param p.cut A cutoff of p values. It is ignored if \code{p.adj = TRUE}. Default is 0.05.
#' @param p.adj.cut A cutoff of adjusted p values. It is ignored if \code{p.adj = FALSE}. Default is 0.05.
#' @param auc.power A cutoff of absolute value of auc power. It is used when
#'   \code{method = "auc"}. The auc power is defined as (AUC-0.5) * 2. Default is 0.4.
#' @param cox.test.cut A cutoff of p values in proportional hazards assumption
#'   test. Features with p value above this value are kept. Default is 0.05.
#' @param logFC A cutoff of log(fold change). Default is 0.5
#' @param verbose Logical value. If TRUE, print messages.
#' @details
#' Except for DESeq2 method, this calculation is performed separately for each feature in the given data.
#' @return A data.frame with features as rows, and associated statistics as
#'   columns (p-values, AUC power, etc.), depending on the methods.
#' @export
#' @importFrom foreach foreach %do%
#' @importFrom progressr with_progress without_progress
#' @examples
PhenoAssoFeatures <- function(data,
                              phenotype,
                              method,
                              family = NULL,
                              cor.method = "pearson",
                              setNaN = 1,
                              p.adj = F,
                              deseq2_contrast = NULL,
                              adj.method = "BH",
                              p.cut = 0.05,
                              p.adj.cut = 0.05,
                              auc.power = 0.4,
                              cox.test.cut = 0.05,
                              logFC = 0.5,
                              verbose = T){
  if ( ! method %in% c("cox", "glm", "auc", "wilcox", "deseq2", "cor")) {
    stop("method must be one of cox, glm, auc, wilcox, deseq2 and cor !")
  }
  my_sapply <- sapply
  set_progress = ifelse(test = verbose,
                        yes = with_progress,
                        no = without_progress)
  if (verbose) {
    message("Starting feature selection with method ", method, " ", family, "...")
  }
  label.lt <- list(phenotype = phenotype)
  if (class(phenotype) %in% c("character", "factor")) {
    if (length(unique(phenotype)) > 2) {
      label.names <- sort(unique(phenotype))
      label.names <- as.character(label.names)
      message(paste0("Feature selection will be perfromed for each group :", paste(label.names, collapse = ",")))
      label.lt <- sapply(label.names, USE.NAMES = T, function(x){
        label <- ifelse(phenotype == x, 1, 0)
        return(list(label))
        })
    }
  }
  res.lt <- list()
  for (n in names(label.lt)) {
    phenotype <- label.lt[[n]]
    asso_res <- switch (method,
                        "cox" = Coxtest(data = data,
                                        phenotype = phenotype,
                                        set_progress = set_progress,
                                        my_sapply = my_sapply),
                        "glm" = Glmtest(data = data,
                                        phenotype = phenotype,
                                        family = family,
                                        set_progress = set_progress,
                                        my_sapply = my_sapply),
                        "auc" = MarkerTest(data = data,
                                           phenotype = phenotype,
                                           auc.power = auc.power,
                                           set_progress = set_progress,
                                           my_sapply = my_sapply,
                                           verbose = verbose),
                        "wilcox" = WilcoxTest(data = data,
                                              phenotype = phenotype,
                                              set_progress = set_progress,
                                              my_sapply = my_sapply),
                        "deseq2" = DESeq2DETest(data = data,
                                                phenotype = phenotype,
                                                contrast = deseq2_contrast,
                                                verbose = verbose),
                        "cor" = CorrTest(data = data,
                                         phenotype = phenotype,
                                         cor.method = cor.method)
    )
    # filter with p value
    if(method == "cox"){
      if (p.adj) {
        asso_res[, "p.adj"] <- p.adjust(asso_res[, "coef.p"],
                                        method = adj.method)
        asso_res <- subset(asso_res,
                           subset = p.adj <= p.adj.cut & cox.test.p > cox.test.cut)
      } else {
        asso_res <- subset(asso_res,
                           subset = coef.p <= p.cut & cox.test.p > cox.test.cut)
      }

    } else{
      if (method %in% c("glm", "wilcox", "deseq2", "cor")){
        asso_res[is.nan(asso_res[,"p.value"]), "p.value"] <- setNaN
        if (p.adj) {
          asso_res[, "p.adj"] <- p.adjust(asso_res[, "p.value"],
                                          method = adj.method)
          asso_res <- subset(asso_res, subset = p.adj <= p.adj.cut)
        } else {
          asso_res <- subset(asso_res, subset = p.value <= p.cut)
        }
      }
    }
    # foldchange
    if ("log2FoldChange" %in% colnames(asso_res)){
      asso_res <- asso_res[abs(asso_res[, "log2FoldChange"]) >= logFC,]
    }
    asso_res <- cbind(features = rownames(asso_res), asso_res)
    res.lt[[n]] <- asso_res
  }
  if (length(res.lt) == 1) {
    res <- res.lt[[1]]
  } else {
    res <- foreach(n = names(res.lt), .combine = "rbind") %do% {
      if (nrow(res.lt[[n]]) != 0){
        res.lt[[n]] <- cbind(res.lt[[n]], group = n)
      }
    }
  }
  return(res)
}
