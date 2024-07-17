# Internal function to find neighbors for each spot within given clusters
#' @param st_pos A data.frame with coordinate information in two columns.
#' @param cluster Cluster labels of spots.
#' @param r.dist A numeric value. Distance to define the neighboring spots. Default
#'   is 2.
#' @param return.list Logical value. If TURE, return the neighbors of each spot
#'   using list.
#' @return List or matrix.
FindNeispots <- function(st_pos,
                         cluster = NULL,
                         r.dist = 2,
                         return.list = T){
  dist.mtrx <- as.matrix(dist(st_pos))
  dist.mtrx <- round(dist.mtrx, digits = 2)
  spot.nei <- foreach(i = rownames(dist.mtrx), .combine = "rbind") %do% {
    d <- dist.mtrx[i, ]
    data.frame(spot1 = i,
          spot2 = rownames(dist.mtrx)[which(d <= r.dist)],
          dist = d[d <= r.dist])
  }
  rownames(spot.nei) <- NULL
  if ( ! is.null(cluster)) {
    if (is.null(names(cluster))) {
      names(cluster) <- rownames(st_pos)
    }
    spot.nei <- spot.nei[cluster[spot.nei$spot1] == cluster[spot.nei$spot2],]
  }
  if ( return.list) {
    spot.nei <- split(x = spot.nei$spot2,
                      f = spot.nei$spot1)
    }
  return(spot.nei)
}

# Internal function to get weights based on distance.
# @param d A numeric vector of distance.
# @param adjust Adjustment of weights.
# @return A vector of weights named with distance (d).
KernelWeight <- function(d,
                         adjust = 1){
  if (length(d) == 1) {
    weight = 1
  } else {
    sd <- adjust * max(d)/3  # d = 3En = 3sd; sd != 0
    weight <- dnorm(d, sd = sd)
    names(weight) <- d
    # weight sum keep to 1
    weight <- weight/sum(weight)
  }
  return(weight)
}

#' Calculating enriched value in unit region.
#' @description For each spot, the unit is defined as the central spot and its
#'   neighboring spots up to n-layers. Features within this unit could be
#'   integrated.
#' @param data Data containing features of spots.
#' @param d.lt A lists of spots grouped by distance.
#' @param method Methods to aggregate values from neighbor layers.
#' @param weight Weights for each layer.
#' @param layer.method Methods to aggregate values for each neighbor layer.
#' @param verbose Logical value. If TURE, print messages.
#' @return Integrated value of feature in this unit.
#' @name unit
#' @rdname unit
#' @importFrom foreach foreach %do%
UnitEnrich <- function(data,
                       d.lt,
                       method = "mean",
                       weight = NULL,
                       layer.method = "mean",
                       verbose = F
){
  data <- data.frame(data, check.names = F)
  if (verbose) {
    pb <- txtProgressBar(min = 0,
                         max = length(d.lt),
                         style = 3)
    i = 0
  }
  if (method == "mean") {
    es <- foreach(d = d.lt, .combine = "rbind") %do% {
      if (verbose) {
        i = i + 1
        setTxtProgressBar(pb = pb, value = i)
      }
      colMeans(data[names(d),, drop = F])
    }
  }
  if (method == "weighted") {
    if (is.null(weight)) {
      stop("weight param must be provided for weighted method !")
    }
    fun <- switch(layer.method,
                  "mean" = mean,
                  "sum" = sum)
    es <- foreach(d = d.lt, .combine = "rbind") %do% {
      score <- apply(data[names(d), , drop = F], 2, function(x){
        d_mean <- tapply(x, d, FUN = fun)
        d_mean[is.na(d_mean)] <- 0
        d_mean %*% weight
      })
      if (verbose) {
        i = i + 1
        setTxtProgressBar(pb = pb, value = i)
      }
      return(score)
    }
  }
  if (verbose) {
    close(con = pb)
  }
  rownames(es) <- names(d.lt)
  return(es)
}

#' Integration of features for all defined regions in sample.
#' @param data Data containing features of spots.
#' @param st_pos A data.frame with coordinate information in two columns.
#' @param cluster Clusters of spots.
#' @param r.dist A numeric value. Distance to define the neighboring spots.
#' @param method Methods to aggregate values from neighbor layers.
#' @param layer.method Methods to aggregate values for each neighbor layer.
#' @param adjust Adjustment of weights.
#' @param verbose Logical value. If TURE, print messages.
#' @return Integrated values of features for each unit.
#' @export
#' @rdname unit
RegFeaEnrich <- function(data,
                         st_pos,
                         cluster = NULL,
                         r.dist = 2,
                         method = "weighted",
                         layer.method = "mean",
                         adjust = 2,
                         verbose = T){
  data <- data.frame(data, check.names = F)
  rownames(data) <- rownames(st_pos)
  nei.pairs <- FindNeispots(st_pos = st_pos,
                            r.dist = r.dist,
                            cluster = cluster,
                            return.list = F)# calculating Euclidean Distance
  dist <- nei.pairs[, "dist"]
  names(dist) <- nei.pairs[, "spot2"]
  d.lt <- split(x = as.factor(dist),
                f = nei.pairs[, "spot1"])
  #
  if (method == "weighted") {
    uniq.d <- unique(as.numeric(dist))
    weight <- KernelWeight(d = uniq.d,
                            adjust = adjust)
  }
  res <- UnitEnrich(data = data,
                     d.lt = d.lt,
                     method = method,
                     weight = weight,
                    layer.method = layer.method,
                    verbose = verbose)
  rownames(res) <- rownames(data)
  return(res)
}

