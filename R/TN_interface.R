# Internal function to calculate I score.
CalIscoreUnit <- function(x,
                         w.mtr,
                         verbose = T){
  n <- length(x)
  if (n != ncol(w.mtr)) {
    stop("objects of different length")
  }
  if (n < 2){
    I <- 0
  } else{
    m <- x - mean(x)
    W <- sum(w.mtr)
    #Z <- matrix(m,nrow = n)%*%matrix(m,ncol = n)
    Z <- matrix(NA, nrow = n, ncol = n)
    for(i in 1:n){
      for(j in 1:n){
        Z[i, j] <- (m[i] + m[j])^2
      }
    }
    I <- (1/W)*(sum((w.mtr*Z)))
  }
  return(I)
}

# Internal function to filter and abstract spatially connected graph.
#' @importFrom igraph graph_from_data_frame components

FiltConnectGraph <- function(st_pos,
                             nei.dist,
                             spot_class = NULL,
                             min.spots = 8,
                             min.class.pt = 0,
                             max.class.pt = 0){
  if (nrow(st_pos) < min.spots) {
    keep_g <- NULL
  } else {
    d <- FindNeispots(st_pos,
                      r.dist = nei.dist,
                      return.list = F)
    g <- graph_from_data_frame(d = d,
                               vertices = data.frame(name = rownames(st_pos)))
    comp <- components(g)$membership
    comp <- split(names(comp), comp)
    keep_g <- comp[sapply(comp, length) >= min.spots]
    if (! is.null(spot_class)){
      spot_class <- as.factor(spot_class)
      class_pt <- lapply(keep_g, function(x){
        table(spot_class[x])/length(x)})
      keep_g <- keep_g[lapply(class_pt, max) <= max.class.pt]
      keep_g <- keep_g[lapply(class_pt, min) >= min.class.pt]
    }
  }
  return(keep_g)
}

#' Tumor-normal interface (TNI) scores
#' @param ES Tumor cell abundances or other features representing the enrichment
#'   of tumor cells.
#' @param st_pos A data.frame with coordinate information in two columns.
#' @param cluster Cluster labels of spots.
#' @param r.dist A numeric value. Distance to define the neighboring spots.
#'   Default is 2.
#' @param norm01 Logical value indicating whether normalizing the TNI scores
#'   within smaple to range 0-1.
#' @return A vector of TNI scores.
#' @export
#' @examples
TNIscore <- function(ES,
                       st_pos,
                       cluster = NULL,
                       r.dist = 2,
                       norm01 = F){
  names(ES) <- rownames(st_pos)
  reg.lt <- FindNeispots(st_pos,
                         r.dist = r.dist,
                         cluster = NULL,
                         return.list = T)
  reg.lt <- reg.lt[names(ES)]
  I.score <- foreach(i = reg.lt, .combine = "c") %do%{
    D <- 1/as.matrix(dist(st_pos[i,]))
    diag(D) <- 0
    unit_I <- CalIscoreUnit(x = ES[i],
                              w.mtr = D)
    return(unit_I)
  }
  I.score <- log2(1 + I.score)
  if (! is.null(cluster)){
    cluster <- as.character(cluster)
    names(cluster) <- rownames(st_pos)
    range <- foreach(i = reg.lt, .combine = "c") %do%{
      clu.es <- tapply(ES[i], cluster[i], mean)
      (max(clu.es) - min(clu.es))/exp(2*min(clu.es))
    }
  } else{
    range <- foreach(i = reg.lt, .combine = "c") %do%{
      es <- sort(ES[i], decreasing = T)
      if (length(es) > 2){
        up.es <- mean(es[1:2])
        low.es <- mean(es[length(es)-1:length(es)])
        d <- (up.es - low.es)/exp(2*min(low.es))
      } else{
        d <- 0
      }
      return(d)
    }
  }
  edge.score <- apply(cbind(I.score, range), 1, function(x)exp(mean(log(x))))
  if (norm01) {
    edge.score <- Norm01(edge.score)
  }
  names(edge.score) <- names(ES)
  return(edge.score)
}

#' Identification of tumor-normal interface (TNI) regions
#' @param x A vector of TNI scores.
#' @param st_pos A data.frame with coordinate information in two columns.
#' @param r.dist A numeric value. Distance to define the neighboring spots.
#'   Default is 2.
#' @param maxval A numeric value. Spots with TNI scores above this value are
#'   defiend as TNI spots.
#' @param minval A numeric value. Spots with TNI scores below this value are
#'   defiend as non-TNI spots.
#' @param boun.nei.n A numeric value. The number of neighboring spots belong to
#'   TNI.
#' @param candi.step An integer indicating the number of steps to add
#'   \code{candi} spots when counting neighboring TNI spots. The higher the
#'   value, the more contiguous the TNI area obtained. Default is 1.
#' @param max.iter.num The maximum number of times to update the spot labels.
#' @param connect.dist Distance to define the spatially connected TNI regions.
#' @param min.spots A numeric value. Connected TNI regions with the number of
#'   spots that fall below this value will be filtered out.

#' @param verbose Logical value. If TURE, print messages.
#'
#' @return A vector of labels named with spots.
#' @export
#' @examples
DefineTNIregion <- function(x,
                     st_pos,
                     r.dist = 2,
                     maxval = 0.08,
                     minval = 0.03,
                     boun.nei.n = 2,
                     candi.step = 1,
                     max.iter.num = 20,
                     min.spots = 10,
                     connect.dist = 2,
                     verbose = T){
  if (is.null(names(x))) {
    names(x) <- rownames(st_pos)
  }
  if (verbose) {
    print("Identifing TNI spots from edge scores ...")
  }
  reg.lt <- FindNeispots(st_pos = st_pos,
                         r.dist = r.dist,
                         return.list = T)
  spotype <- ifelse(x > maxval,
                    "TNI",
                    ifelse(x < minval, "nTNI", "candi"))
  spotype <- factor(spotype, levels = c("TNI", "nTNI", "candi"))
  for (it in 1:max.iter.num){
    type.n.raw <- table(spotype)
    for (candi in names(spotype[spotype == "candi"])){
      nei.spots <- candi
      for (k in 1:candi.step) {
        nei <- unlist(reg.lt[nei.spots])
        nei.spots <- append(nei.spots,
                            nei[spotype[nei] != "nTNI"])
      }
      nei.spots <- unique(unlist(reg.lt[nei.spots]))
      boun.n <- sum(spotype[nei.spots] == "TNI")
      if (boun.n >= boun.nei.n) {
        spotype[candi] = "TNI"
      }
    }
    type.n.new <- table(spotype)
    if (all(type.n.raw == type.n.new)){
      break
    } else {
      type.n.raw = type.n.new
    }
  }
  spotype[spotype == "candi"] = "nTNI"
  if (verbose) {
    print("removing unconnected TNI regions ...")
  }
  # filter connected regions
  edgespots <- names(spotype)[spotype == "TNI"]
  if (length(edgespots) == 0) {
    warning("There is no spot defined as TNI!")
    final.edgespots = NULL
  } else {
    final.edgespots <- FiltConnectGraph(st_pos = st_pos[edgespots,],
                                        nei.dist = connect.dist,
                                        min.spots = min.spots)
    final.edgespots <- unname(unlist(final.edgespots))
    # plot
    if(verbose){
      print("Completing the identificaton of TNI spots !")
    }
  }
  alltype <- rep("nTNI", nrow(st_pos))
  names(alltype) <- rownames(st_pos)
  alltype[final.edgespots] <- "TNI"
  return(alltype)
}

# Internal function to perform unsupervised clustering.
#' @param x Numeric matrix of data.
#' @param clu_num The number of cluster.
#' @param method The method used for clustering. One of
#' "kmeans", "hclust" and "pam".
#' @return A vector of clusters.
#' @importFrom cluster pam
UnsupCluster <- function(x,
                         clu_num,
                         method){
  # cluster method is one of kmeans, hclust and pam
  if (! method %in% c("kmeans", "hclust", "pam")){
    stop("Clustering method must be one of kmeans, hclust and pam")
  }
  if (method == "kmeans"){
      clu <- kmeans(x = x,
                    centers = clu_num,
                    iter.max = 20)
      groups <- clu$cluster
  }
  if (method == "hclust" ){
    d <- dist(x, method = "euclidean")
    fit <- hclust(d,
                  method ="ward.D")
    groups <- cutree(fit,
                     k = clu_num)
  }
  if (method == "pam"){
    pk <- pam(x,
              k = clu_num,
              metric = "euclidean")
    groups <- pk$clustering
  }
  return(groups)
  }


#' Grouping TNI spots into distinct types.
#' @param TNI_pos Position of TNI spots. A data.frame with coordinate
#'   information in two columns.
#' @param cluster A vector of cluster labels of TNI spots.
#' @param seed Seed for clustering spots.
#' @param type_n Number of types to cluster.
#' @param cluster_method Clustering method, one of "kmeans", "hclust", and
#'   "pam".
#' @param nei.dist A numeric value. Distance to define the neighboring spots.
#'   Default is 4.
#' @param connect.dist A numeric value. Distance to define the spatially
#'   connected TNI regions within each type. Default is 4.
#' @param min.type.spots A numeric value. Types with number of spots less than this value will be filtered.
#' @param disconnect.spots A numeric value. For each TNI type, connected TNI regions with the number of
#'   spots that fall below this value will be filtered out. Default is 3.
#' @param verbose Logical value. If TURE, print messages.
#' @importFrom fpc pamk
#' @return A vector of TNI types named with spots.
#' @export
#' @examples
GroupTNItypes <- function(TNI_pos,
                          cluster,
                          seed = 123,
                          type_n = NULL,
                          cluster_method = "pam",
                          nei.dist = 4,
                          connect.dist = 4,
                          min.type.spots = 0,
                          disconnect.spots = 3,
                          verbose = T){
  if (verbose) {
    message("Separating tumor boundary to different type according it's label composition of spot-neighbors")
  }
  names(cluster) <- rownames(TNI_pos)
  spot_label <- as.factor(cluster)
  # boundary consist of diff-clusters were separated
  nei.spots <- FindNeispots(TNI_pos,
                            r.dist = nei.dist,
                            return.list = T)
  clu.mt <- foreach(i = names(nei.spots), .combine = "rbind") %do% {
    te <- table(spot_label[nei.spots[[i]]])
    #min.c <- floor(length(nei.spots[[i]])/5)
    #te[te <= min.c] <- 0
    te <- te/sum(te)
    return(te)
  }
  rownames(clu.mt) <- names(nei.spots)
  if ( is.null(type_n)){
    set.seed(seed = seed)
    pamk.clu <- pamk(clu.mt)
    type_n <- pamk.clu$nc
    if (verbose) {
      message("pamk: best number of clusters is ", pamk.clu$nc)
    }
  }
  set.seed(seed = seed)
  type <- UnsupCluster(x = clu.mt,
                       clu_num = type_n,
                       method = cluster_method)
  # remove type with spots less than cutoff
  if (verbose) {
    message("Filtering short boundary type and remove discontinuous spots")
  }
  lt <- split(names(type), type)
  type <- lapply(lt, function(x){
    FiltConnectGraph(st_pos = TNI_pos[x, ],
                     nei.dist = connect.dist,
                     min.spots = disconnect.spots)
  })
  type <- type[lapply(type, length) >= min.type.spots]
  # names
  final.type <- foreach(n = 1:length(type),
                        .combine = "c") %do% {
                          x <- unlist(type[[n]])
                          setNames(object = rep(paste0("type_", n), length(x)),
                                   nm = x)
                          }
  return(final.type)
}
#' Categories of TNI regions
#'
#' @param type A vector of type labels of TNI spots.
#' @param cluster A vector of cluster labels of TNI spots.
#' @param ES Tumor cell abundances or other features representing the enrichment
#'   of tumor cells.
#'
#' @return A vector of TNI categories named with spots.
#' @export
#'
#' @examples
TNIClass <- function(type, cluster, ES){
  types <- sort(unique(type))
  types <- types[types != "others"]
  for (t in types) {
    t.boun.spots <- type == t
    t.es.mean <- mean(ES[t.boun.spots])
    clu.mean <- tapply(ES[t.boun.spots], cluster[t.boun.spots], mean)
    boun.class <- ifelse(clu.mean[cluster[t.boun.spots]] > t.es.mean, "T_boun", "N_boun")
    type[t.boun.spots] <- boun.class
  }
  return(type)
}

