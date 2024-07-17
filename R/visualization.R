default.colors <- c("#336699", "#679436", "#f4a261", "#70a288","#b95b13","#264653", "#a3d5ff", "#60afff","#f9c74f",
                    "#3cb0cd","#fcefb4","#c9a227","#057a55","#0077b6", "#2a9d8f","#87153b", "#797d62","#c32f27","#d199f1","#6247aa",
                    "#023047","#6c809a","#004346", "#4b1d3f","#ffa69e","#00838f","#dbb68f","#5e6472","#82a0bc", "#ced4da",  "#a3b18a")
#' @title Colors of labels.
#' @param labels A vector of labels to set colors. It should be a character or
#'   factor vector.
#' @param na.col Color for NA values.
#' @param assgin.col A vector of labels assigned colors by users. Colors and
#'   labels included in this vector will not be assigned again.
#' @param colors Characters vector of colors assigned to labels. If NULL, the
#'   default colors are used.
#' @param seed An integer to set the seed.
#' @param random Logical indicating whether colors are assigned to labels
#'   sequentially or randomly.
#' @return A vector of colors named with labels
#' @export
#'
#' @examples
#' x <- c("A", "B", "C", "D")
#' cols <- LabelMapcolor(x, assgin.col = c(A = "red"))
LabelMapcolor <- function(labels,
                          na.col = "grey",
                          assgin.col = NULL,
                          colors = NULL,
                          seed = NULL,
                          random = F) {
  if (! (is.character(labels) | is.factor(labels))) {
    stop("labels need to be character or factor vector")
  }
  if (is.factor(labels)) {
    labels <- as.character(labels)
  }
  if (is.null(colors)) {
    colors <- default.colors
  }
  if (any(duplicated(colors))) {
    warning("There are duplicated colors!")
  }
  key <- unique(labels)
  key <- key[! is.na(key)]
  key <- key[! key %in% names(assgin.col)]
  colors <- setdiff(colors, assgin.col)
  if (length(key) > length(colors)) {
    warning("The number of labels are more than that of colors, and some colors will be repeated!")
  }
  if (length(key) >=1) {
    if (random) {
      if (! is.null(seed)){
        set.seed(seed = seed)
      }
      mapcol <- sample(colors, length(key))
      names(mapcol) <- key
    } else {
      map <- suppressWarnings(cbind(key, colors))
      mapcol <- map[1:length(key), 2]
      names(mapcol) <- map[1:length(key), 1]
    }
    mapcol <- c(mapcol, assgin.col)
  } else {
    mapcol <- assgin.col
  }
  # add na and assigned colors
  fmapcol <- mapcol[labels]
  fmapcol[is.na(labels)] <- na.col
  return(fmapcol)
}

#' @title Visualization of ST spots with features
#' @param st Seurat object for ST data. If available, the coordinate information
#'   should be included in the \code{metadata} as separate columns labeled 'x'
#'   and 'y'.
#' @param pos A data.frame with two columns named as 'x' and 'y', representing
#'   the coordinate information; One of st and pos must be provided.
#' @param meta A vector of the feature to visualize. If st is available, it can
#'   be one of the feature name in the metadata.
#' @param feature A character of gene name in \code{data} slot of st. Provide
#'   either \code{meta} or \code{feature}, not both.
#' @param size Size of the spots.
#' @param savefile File name to save this plot. Only supporting pdf files.
#' @param return Logical indicating if the plot should be print or returned.
#' @param cha.col Colors of meta if the meta is a character vector. By default,
#'   ggplot2 assigns colors.
#' @param title Title of plot.
#' @param legend.name Name of the legend.
#' @param p.width The width of pdf file. The default value is 10.
#' @param p.height The height of pdf file. The default value is 10.
#' @param size.factor A numeric value, scale the size of spots based on the
#'   p.width and p.height. It will be ignored if size is set. The default value
#'   is 0.02.
#' @param f.color Fill colors of spots.
#' @param scale_y_reverse Logical indicating whether reversing the y position.
#' @param scale_x_reverse Logical indicating whether reversing the x position.
#' @param limits Scaling limits. Any points outside these limits will not be
#'   plotted.
#' @param na.col Color for NA values.
#' @param ... Other arguments passed on to \code{theme} in \pkg{ggplot}.
#'
#' @return A ggplot object if \code{return = TRUE}.
#' @export
#' @examples
#' @importFrom ggplot2 ggplot aes geom_point theme_void coord_fixed labs theme
#'   scale_y_reverse scale_x_reverse scale_colour_gradientn scale_color_manual
#'   element_text
#' @importFrom grDevices pdf dev.off
SpotVisualize <- function(st = NULL , pos = NULL, meta = NULL, feature = NULL, size = 2.5,
                          savefile = NULL, return = F, cha.col = NULL,
                          title = "NULL", legend.name = NULL,
                          p.width = 10, p.height = 10, size.factor = 0.02,
                          f.color = c("#0077b6","lightyellow","#c32f27"),
                          scale_y_reverse = F, scale_x_reverse = F,
                          limits = NULL, na.col = "#ced4da", ...){
  if (is.null(pos)) {
    if ( ! is.null(st)) {
      pos <- st@meta.data[, c("x", "y")]
    }
    else {
      stop("Spatial object (st) or pos must be provided with columns named as x and y!")
    }
  }
  if (! is.null(meta)) {
    if (length(meta) == 1){
      if ( ! is.null(st)) {
        label = st@meta.data[, meta]
      }
    } else {
      label = meta
    }
  }
  if (! is.null(feature)) {
    if (is.null(st)) {
      stop(paste0("st object must provided for feature ", feature))
    } else {
      label = st@assays$SCT@data[feature,]
    }
  }
  if(is.null(size)){
    max_num_spots <- max(c(max(pos[, "x"])-min(pos[, "x"])),
                         c(max(pos[, "y"])-min(pos[, "y"])))
    min_len <- min(p.width, p.height)
    size <- min_len/(size.factor*max_num_spots)
    size <- round(size, 1)
  }
  # if col is meta
  gp = ggplot(data = data.frame(pos, label), aes(x = x, y = y)) +
    theme_void() +
    coord_fixed() +
    labs(title = title, fill = legend.name, color = legend.name) +
    theme(plot.title = element_text(hjust = 0.5, size = 30),
          legend.text = element_text(size = 18),
          legend.title = element_text(size=20)) +
    theme(...)
  if (scale_y_reverse) {
    gp = gp + scale_y_reverse()
  }
  if (scale_x_reverse) {
    gp = gp + scale_x_reverse()
  }
  p <- gp + geom_point(aes(colour = label), size = size)
  if (class(label) == "numeric") {
    if (is.null(limits)) {
      limits <- c(min(label[!is.na(label)]), max(label[!is.na(label)]))
    }
      p <- p + scale_colour_gradientn(colors = f.color, limits = limits)
  }
  if (class(label) == "character") {
    if (! is.null(cha.col)) {
      p <- p + scale_color_manual(values = cha.col, na.value = na.col)
    }
  }
  if (! is.null(savefile)) {
    pdf(savefile, height = p.height, width = p.width)
    print(p)
    dev.off()
  } else {
    if (return) {
      return(p)
    } else {
      print(p)
    }
  }
}
#
#' Visualization of spots with tumor cell abundances and tumor-normal interface
#' (TNI) regions
#' @param abun Numeric vector of tumor cell abundances.
#' @param label labels of spots.
#' @param l_nshow labels represent nonTNI. Spots with these labels will not be
#'   highlighted.
#' @param st Seurat object for ST data. If available, the coordinate information
#'   should be included in the \code{metadata} as separate columns labeled 'x' and 'y'.
#' @param pos A data.frame with two columns named as 'x' and 'y', representing
#'   the coordinate information; One of st and pos must be provided.
#' @param size Size of the spots.
#' @param shape A numerical value to set the shape of the TNI spots. Default: 21.
#' @param stroke A numerical value to set the stroke thickness.
#' @param title Title of plot.
#' @param legend.name Name of the legend.
#' @param savefile File name to save this plot. Only supporting pdf files.
#' @param return Logical indicating if the plot should be print or returned.
#' @param scale_y_reverse Logical indicating whether reversing the y position.
#' @param scale_x_reverse Logical indicating whether reversing the x position.
#' @param f.color Fill colors of spots.
#' @param limits Scaling limits. Any points outside these limits will not be
#'   plotted.
#' @param p.width The width of pdf file. The default value is 10.
#' @param p.height The height of pdf file. The default value is 10.
#' @param line_col The color of the border of TNI spots.
#' @importFrom ggplot2 ggplot aes geom_point theme_void coord_fixed labs theme
#'   scale_y_reverse scale_x_reverse scale_colour_gradientn scale_color_manual
#'   element_text
#' @importFrom grDevices pdf dev.off
#' @export
#' @examples
AbunTNIPlot <- function(abun, label,
                         l_nshow,
                         st = NULL,
                         pos = NULL,size = 2.5,
                         shape = 21,
                         stroke = 1.2,
                         title = NULL,
                         legend.name = NULL,
                         savefile = NULL, return = F,
                         scale_y_reverse = F,
                         scale_x_reverse = F,
                         p.width = 10, p.height = 10,
                         f.color = c("#0077b6","lightyellow","#c32f27"),
                         limits = NULL,
                         line_col = "black"){
  if (is.null(pos)) {
    if ( ! is.null(st)) {
      pos <- st@meta.data[, c("x", "y")]
    }
    else {
      stop("Spatial object (st) or pos must be provided !")
    }
  }
  if (is.null(limits)) {
    limits <- c(0, max(abun))
  }
  stroke = rep(stroke, length(label))
  stroke[label %in% l_nshow] = 0
  gp = ggplot(data = data.frame(pos, abun, label, stroke), aes(x = x, y = y)) +
    theme_void() +
    labs(title = title, fill = legend.name) +
    coord_fixed() +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.text = element_text(size = 16),
          legend.title = element_text(size= 18))
  if (scale_y_reverse) {
    gp = gp + scale_y_reverse()
  }
  if (scale_x_reverse) {
    gp = gp + scale_x_reverse()
  }
  mid.p <- ifelse(is.null(limits), NA, max(limits)/2)
  p <- gp + geom_point(aes(color = label, fill = abun, stroke = stroke),
                           size = size, shape = shape) +
    scale_fill_gradientn(colors = f.color, limits = limits) +
    scale_color_manual(values = line_col)
  if (! is.null(savefile)) {
    pdf(savefile, height = p.height, width = p.width)
    print(p)
    dev.off()
  } else {
    if (return) {
      return(p)
    } else {
      print(p)
    }
  }
}
# visualization of cell proportion
#' Spatial plot of cell distribution
#'
#' @param decon_mtrx A numeric matrix or data.frame containing cell abundances or
#'   percents for each spot. Rows are spots and columns are cell types.
#' @param st_pos A data.frame with coordinate information in two columns named 'x'
#'   and 'y'.
#' @param tarCells Cell types selected to display. If NULL, all cell types in
#'   decon_mtrx are used.
#' @param savefile File name to save this plot. Only supporting pdf files.
#' @param separate Logical value. If TRUE, the default, each cell types will
#'   plot separately. If FALSE, the distribution of all selected cell types are
#'   showed with pie charts for each spot.
#' @param f.color Fill colors of spots.
#' @param e.color Edge color of pie.
#' @param pie_scale Amount to scale the pie size.
#' @param pie_color A character vector. Colors of cell types in pie.
#' @param size Size of spots when separate is TRUE.
#' @param size.auto Logical indicating if adjusting the size automatically based
#'   on sample size. It is useful when drawing multiple samples.
#' @param size.factor A numeric value, scale the size of spots when \code{size.auto = TRUE}.
#' @param numCol The number of columns in one pdf file if \code{separate = TRUE}.
#' @param p.width The width of pdf file. The default value is 10.
#' @param p.height The height of pdf file. The default value is 10.
#' @param scale_y_reverse Logical indicating whether reversing the y position.
#' @param scale_x_reverse Logical indicating whether reversing the x position.
#'
#' @export
#' @importFrom cowplot plot_grid
#' @importFrom scatterpie geom_scatterpie
#' @importFrom ggplot2 ggplot aes geom_point theme_void coord_fixed labs theme
#'   scale_y_reverse scale_x_reverse scale_colour_gradientn scale_color_manual
#'   element_text
#'   scale_fill_manual
#' @importFrom grDevices pdf dev.off
#' @examples
PlotCellSpot = function(decon_mtrx,
                        st_pos,
                        tarCells = NULL,
                        savefile = NULL,
                        separate = T,
                        f.color = c("#0077b6","lightyellow","#c32f27"),
                        e.color = "black",
                        pie_scale = 0.7,
                        pie_color = NULL,
                        size = 5,
                        size.auto = F,
                        size.factor = 0.36,
                        numCol = 4,
                        p.width = 10,
                        p.height = 10,
                        scale_y_reverse = F,
                        scale_x_reverse = F){
  gp = ggplot() +
    coord_fixed() +
    theme_void()  +
    theme(plot.title = element_text(hjust = 0.5, size = 24),
          legend.text = element_text(size = 18),
          legend.title = element_text(size=20))
  if(scale_y_reverse){
    gp = gp + scale_y_reverse()
  }
  if(scale_x_reverse){
    gp = gp + scale_x_reverse()
  }
  if(size.auto){
    max_num_spots <- max(c(max(st_pos[, "x"])-min(st_pos[, "x"])),
                  c(max(st_pos[, "y"])-min(st_pos[, "y"])))
    #numRow <- ceiling(length(tarCells)/numCol)
    size <- size.factor * log2(max_num_spots)
    size <- round(size, 2)
  }
  plotSingle = function(cell, title, type = "pie"){
    p = gp + labs(title = title)
    if(type == "pie"){
      others <- setdiff(colnames(decon_mtrx), cell)
      if (length(others)){
        data <- data.frame(decon_mtrx[, cell, drop = F],
                           others = rowSums(decon_mtrx[, others, drop = F]))
      } else {
        data <- decon_mtrx
      }
      p = p + geom_scatterpie(data = data.frame(st_pos, data, check.names = F),
                              aes(x = x,y = y),
                              cols = colnames(data),
                              pie_scale = pie_scale,
                              color = e.color
                              )
      if ( ! is.null(pie_color)) {
        p <- p + scale_fill_manual(values = pie_color)
      }
    }
    if(type == "point"){
      p = p +
        geom_point(data = data.frame(st_pos, percent = decon_mtrx[, cell]),
                   aes(x = x, y = y, colour = percent), size = size)
    }
    if(length(f.color)>1){
      p <- p + scale_colour_gradientn(values = seq(0,1,0.25),
                                     colours = f.color)
    }
    return(p)
  }
  if (is.null(tarCells)) {
    tarCells <- colnames(decon_mtrx)
  }
  if(separate){
    plot.list = c()
    for(cell in tarCells){
      if(cell %in% colnames(decon_mtrx)){
        plot.list[[cell]] = plotSingle(cell,
                                       title = cell,
                                       type = "point")
      }
    }
    plot = plot_grid(plotlist = plot.list, ncol = numCol)
  } else {
    plot = plotSingle(tarCells, title = "percent for cells", type = "pie")
  }
  if(! is.null(savefile)){
    pdf(file = savefile, width = p.width, height = p.height)
    print(plot)
    dev.off()
  }
  else{return(plot)}
}
