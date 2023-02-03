#' Construct a theoretical variogram
#'
#' @param n_cols The total number of columns.
#' @param n_rows The total number of rows.
#' @param var_R The total error variance. By default, \code{var_R = 1}.
#' @param prop_spatial The proportion of spatial error variance to total error
#'   variance (spatial + random). By default, \code{prop_spatial = 0.5}.
#' @param col_cor The column autocorrelation value.
#' @param row_cor The row autocorrelation value.
#' @param plot When TRUE (default), the theoretical variogram is displayed graphically.
#'   Otherwise, a data frame is returned.
#'
#' @return Graphic of the theoretical variogram, where the x- and y- axes display the row and
#'   column displacements and the z-axis displays the semi-variance (variogram ordinates).
#'   When \code{plot = FALSE}, a data frame is returned with the column and row displacements
#'   as well as the theoretical semi-variances.
#'
#' @examples
#' # Plot a theoretical variogram for a field with 10 columns and 20 rows,
#' # using column and row autocorrelations of 0.4 and 0.8.
#'
#' theoretical_variogram(
#'   n_cols = 10,
#'   n_rows = 20,
#'   var_R = 1,
#'   prop_spatial = 0.5,
#'   col_cor = 0.4,
#'   row_cor = 0.8,
#'   plot = TRUE
#' )
#' @export
#'
theoretical_variogram <- function(n_cols,
                                  n_rows,
                                  var_R = 1,
                                  prop_spatial = 0.5,
                                  col_cor,
                                  row_cor,
                                  plot = TRUE){

prop_rand <- 1 - prop_spatial
col.displacement <- rep(0:(n_cols-1), each = n_rows)
row.displacement <- rep(0:(n_rows-1), times = n_cols)
theoretical.df <- data.frame(col.dis = col.displacement,
                             row.dis = row.displacement,
                             semi.var = var_R * (prop_rand + prop_spatial * (1 - col_cor^(col.displacement) * row_cor^(row.displacement))))
theoretical.df$semi.var[1] <- 0
theoretical.df$col.dis <- as.numeric(trimws(theoretical.df$col.dis))
theoretical.df$row.dis <- as.numeric(trimws(theoretical.df$row.dis))

if(!plot){theoretical.df$col.dis <- factor(theoretical.df$col.dis)
          theoretical.df$row.dis <- factor(theoretical.df$row.dis)
          print(theoretical.df)}

if(plot){
    require(lattice)
    lattice.options(layout.heights=list(bottom.padding=list(x=-1), top.padding=list(x=-1.5)),
                    layout.widths=list(left.padding=list(x=-1.25), right.padding=list(x=-3)))
    par(mar=c(5.1,4.1,4.1,2.1))
    p1 <- wireframe(semi.var ~ row.dis * col.dis, data = theoretical.df, drape = T, colorkey = F, zoom = 0.97, cuts = 30,
                    screen = list(z = 30, x = -60, y = 0), aspect = c(1, 0.66),
                    scales = list(distance = c(1.2, 1.2, 0.5), arrows = F, cex = 0.7, col = "black"),
                    zlab = list(label = paste("Semi-variance"), cex = 0.9, rot = 90, just = c(0.5,-2.25)),
                    xlab = list(label = paste("Row displacement"), cex = 0.9, rot = 19, just = c(0.5,-0.75)),
                    ylab = list(label = paste("Column displacement"), cex = 0.9, rot = -49, just = c(0.5,-0.75)),
                    par.settings = list(axis.line = list(col = "transparent"), clip = list(panel = "off")))
    print(p1)
}
}

