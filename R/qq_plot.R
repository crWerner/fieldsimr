#' Construct a qqplot
#'
#' @param df A data frame containing the columns "col", "row", and the effect to be plotted.
#' @param effect The name of the effect to be plotted.
#' @param labels When TRUE (default), the column and row labels are inserted onto the qqplot.
#'   Otherwise, data points without labels are plotted. In this case, the data frame does not require
#'   the columns "col" and "row", just the effect to be plotted.
#' @param plot When TRUE (default), the qqplot is displayed graphically.
#'   Otherwise, a data frame is returned.
#'
#' @return Graphic of the qqplot, where the x- and y- axes display the theoretical and
#'   sample quantiles. When \code{plot = FALSE}, a data frame is returned with the column and row displacements
#'   as well as the estimated semi-variances.
#'
#' @examples
#' # Plot the qqplot for the spatial error component in the data frame error_df.
#'
#' qq_plot(
#'   error_df,
#'   effect = "e_spat",
#'   labels = TRUE,
#'   plot = TRUE,
#' )
#' @export
#'
#'

qq_plot <- function(df,
                    effect,
                    labels = TRUE,
                    plot = TRUE){

  require(ggplot2)
  if(!labels){qq.df <- data.frame(resid = df[[effect]])
  p1 <- ggplot(qq.df, aes(sample = resid)) + stat_qq()
  qq.df <- data.frame(sample = ggplot_build(p1)$data[[1]]$sample,
                      theoretical = ggplot_build(p1)$data[[1]]$theoretical)
  if(!plot){print(qq.df)}
  if(plot){p1 <- ggplot(data = qq.df, aes(x = theoretical, y = sample)) +
    stat_qq_line(data = qq.df, aes(sample = sample), colour = "steelblue", linewidth = 0.75, inherit.aes = F) +
    geom_point(size = 2) + labs(y = "Sample quantiles", x = "Theoretical quantiles") +
    theme(title = element_text(size = 10),
          axis.title.x = element_text(size = 12),
          axis.text.x = element_text(size = 9),
          axis.title.y = element_text(size = 12),
          axis.text.y = element_text(size = 9))
  print(p1)}
  }

  if(labels){
    colnames(df) <- tolower(colnames(df))
    if (any(!c("col", "row") %in% colnames(df))) {
      stop("'df' must contain columns 'col' and 'row' in order to plot labels")
    }
    qq.df <- data.frame(col = df[["col"]],
                        row = df[["row"]],
                        resid = df[[effect]])
    qq.df$col <- factor(as.numeric(trimws(qq.df$col)))
    qq.df$row <- factor(as.numeric(trimws(qq.df$row)))
    p1 <- ggplot(qq.df, aes(sample = resid)) + stat_qq()
    qq.df <- data.frame(col = qq.df$col[order(qq.df$resid)],
                        row = qq.df$row[order(qq.df$resid)],
                        sample = ggplot_build(p1)$data[[1]]$sample,
                        theoretical = ggplot_build(p1)$data[[1]]$theoretical)
    qq.df <- qq.df[order(qq.df$col, qq.df$row),]
    rownames(qq.df) <- NULL
    if(!plot){print(qq.df)}
    if(plot){qq.df$name <- paste0(qq.df$col, ":", qq.df$row)
    p1 <- ggplot(data = qq.df, aes(x = theoretical, y = sample, label = name)) +
      stat_qq_line(data = qq.df, aes(sample = sample), colour = "steelblue", linewidth = 0.75, inherit.aes = F) +
      geom_text(size = 4) + labs(y = "Sample quantiles", x = "Theoretical quantiles") +
      ggtitle(label = "Residuals indexed as Col:Row") +
      theme(title = element_text(size = 10),
            axis.title.x = element_text(size = 12),
            axis.text.x = element_text(size = 9),
            axis.title.y = element_text(size = 12),
            axis.text.y = element_text(size = 9))
    print(p1)}
  }
}
