#' Construct a sample variogram
#'
#' @param df A data frame containing the columns "col", "row", and the effect to be plotted.
#' @param effect The name of the effect to be plotted.
#' @param plot When TRUE (default), the sample variogram is displayed graphically.
#'   Otherwise, a data frame is returned.
#' @param min.np The minimum number of pairs that semi-variances will be displayed for.
#'   By default, \code{min.np = 30}.
#'
#' @return Graphic of the sample variogram, where the x- and y- axes display the row and
#'   column displacements and the z-axis displays the semi-variance (variogram ordinates).
#'   When \code{plot = FALSE}, a data frame is returned with the column and row displacements
#'   as well as the estimated semi-variances.
#'
#' @examples
#' # Plot the sample variogram for the spatial error component in the data frame error_df.
#'
#' sample.vario(
#'   error_df,
#'   effect = "e_spat",
#'   plot = TRUE,
#' )
#' @export
#'
sample.vario <- function(df,
                         effect,
                         plot = TRUE,
                         min.np = 30){

  colnames(df) <- tolower(colnames(df))
  if (any(!c("col", "row") %in% colnames(df))) {
    stop("'df' must contain columns 'col' and 'row', and the effect to be plotted.")
  }

  sample.df <- data.frame(col = df[["col"]],
                          row = df[["row"]],
                          effect = df[[effect]])
  sample.df <- sample.df[order(sample.df$col, sample.df$row),]

  col.dis <- abs(outer(as.numeric(trimws(sample.df$col)), as.numeric(trimws(sample.df$col)), FUN = "-"))
  row.dis <- abs(outer(as.numeric(trimws(sample.df$row)), as.numeric(trimws(sample.df$row)), FUN = "-"))
  var.mat <- outer(sample.df$effect, sample.df$effect, FUN = "-")^2/2
  sample.df <- data.frame(col.dis = col.dis[upper.tri(col.dis, diag = T)],
                          row.dis = row.dis[upper.tri(row.dis, diag = T)],
                          semi.var = var.mat[upper.tri(var.mat, diag = T)])
  sample.df <- sample.df[order(sample.df$col.dis, sample.df$row.dis),]

  sample.df <- data.frame(col.dis = rep(unique(sample.df$col.dis), each = length(unique(sample.df$row.dis))),
                          row.dis = unique(sample.df$row.dis),
                          semi.var = c(with(sample.df, tapply(semi.var, list(row.dis,col.dis), function(x) mean(x, na.rm=T)))),
                          np = c(with(sample.df, tapply(semi.var, list(row.dis,col.dis), function(x) length(x[!is.na(x)])))))

  if(!plot){sample.df$col.dis <- factor(sample.df$col.dis)
  sample.df$row.dis <- factor(sample.df$row.dis)
  print(sample.df)}

  if(plot){
    require(lattice)
    lattice.options(layout.heights=list(bottom.padding=list(x=-1), top.padding=list(x=-1.5)),
                    layout.widths=list(left.padding=list(x=-1.25), right.padding=list(x=-3)))
    par(mar=c(5.1,4.1,4.1,2.1))
    p1 <- wireframe(semi.var ~ row.dis * col.dis, data = sample.df[sample.df$np >= min.np,], drape = T, colorkey = F, zoom = 0.97, cuts = 30,
                    screen = list(z = 30, x = -60, y = 0), aspect = c(1, 0.66),
                    scales = list(distance = c(1.2, 1.2, 0.5), arrows = F, cex = 0.7, col = "black"),
                    zlab = list(label = paste("Semi-variance"), cex = 0.9, rot = 90, just = c(0.5,-2.25)),
                    xlab = list(label = paste("Row displacement"), cex = 0.9, rot = 19, just = c(0.5,-0.75)),
                    ylab = list(label = paste("Column displacement"), cex = 0.9, rot = -49, just = c(0.5,-0.75)),
                    par.settings = list(axis.line = list(col = "transparent"), clip = list(panel = "off")))
    print(p1)
  }
}

