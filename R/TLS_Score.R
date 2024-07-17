#' Calculation of TLS scores. TLS features stored in
#' @param data A data.frame should contain three TLS features, including
#'   co-local scores of "Plasma/B.cells_T.cells" and enrichment scores of two
#'   signatures: "LC.50sig", "imprint.65sig".
#' @param st_pos A data.frame with coordinate information in two columns.
#' @param cluster Cluster labels of spots.
#' @param ... Arguments passed to other methods.
#'
#' @return Lists containing distinct forms of TLS features and the final TLS scores.
#' @export
#' @seealso [RegFeaEnrich]
#' @examples
CalTLSfea <- function(data,
                      st_pos,
                      cluster,
                      ...) {
  library(dplyr)
  library(foreach)
  unit.fea <- RegFeaEnrich(data = data,
                           st_pos = st_pos,
                           cluster = cluster,
                           ...)
  rownames(unit.fea) <- rownames(data)
  # Calculating the random values
  rep <- ceiling(1000/nrow(st_pos))
  random.fea <- c()
  num <- 0
  for(i in 1:rep){
    set.seed(i)
    pos <- st_pos[sample(1:nrow(st_pos), nrow(st_pos), replace = F), ]
    temp.fea <- RegFeaEnrich(data = data,
                             st_pos = pos,
                             cluster = NULL,
                             ...)
    if (i == rep) {
      random.fea <- rbind(random.fea, temp.fea[1:(1000 - num),])
    } else {
      random.fea <- rbind(random.fea, temp.fea)
    }
    num <- nrow(random.fea)
  }
  # significance
  signi.fea <- foreach(reg = 1: nrow(unit.fea), .combine = rbind, .final = data.frame) %do% {
    p <- foreach(cell = 1:ncol(unit.fea), .combine = "c") %do%{
      sum(random.fea[, cell] >= unit.fea[reg, cell])/nrow(random.fea)
    }
    p <- (-log10(p+0.001))
    #ifelse(p>=2,1,0)
  }
  signi.fea <- apply(signi.fea, 2, function(x){
    if(any(x !=0)){
      res <- Norm01(x)
    }else {
      res <- x
    }
    return(res)
  })
  norm.t <- apply(unit.fea, 2, function(x){
    if(any(x !=0)){
      res <- Norm01(x)
    }else {
      res <- x
    }
    return(res)
  })
  signi.fea <- signi.fea*norm.t
  rownames(signi.fea) <- rownames(unit.fea)
  colnames(signi.fea) <- colnames(unit.fea)
  TLS.score <- apply(signi.fea, 1, function(x)exp(mean(log(x))))
  return(list(orig.fea = data,
              unit.fea = data.frame(unit.fea, check.names = F),
              random.fea = data.frame(random.fea, check.names = F),
              signi.fea = data.frame(signi.fea, check.names = F),
              TLS.score = TLS.score))
}
