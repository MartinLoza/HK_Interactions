
#' TextSize
#'
#' @param size Text size
#'
#' @return A theme with the selected text size.
#' @export
TextSize <- function(size = 10 ){
  tmp_theme <- theme(text = element_text(size = size))
  return(tmp_theme)
}

#' LegendPosition
#'
#' @param position Legend position
#'
#' @return Theme with the legened in the selected position. Available positions are :"none", "right", "left", "top", and "bottom".
#' @export
LegendPosition <- function(position = "right"){
  tmp_theme <- theme(legend.position = position)
}

#
Y.Axis <- function(rotation = 0, h_adjust = 0 ){
  t <- theme(axis.text.x = element_text(angle = rotation,hjust = h_adjust ))
  return(t)
}

# NoAxis
NoAxis <- function(keep_lines = FALSE){
  t <- theme(axis.title = element_blank(),
             axis.text = element_blank(),
             axis.ticks = element_blank())
  if(keep_lines == FALSE){
    t <- t + theme(axis.line = element_blank())
  }
  return(t)
}


## My histogram 
MyHistogram <- function(.data = NULL, x = NULL, fill = NULL, n_bins = 80, ...){
  p <- .data %>% ggplot(mapping = aes_string(x = x, fill = fill)) + geom_histogram(bins = n_bins, ...) + theme_bw()
  return(p)
}

#MyBoxPlot
MyBoxPlot <- function(.data = NULL, x = NULL, y = NULL, fill = NULL, color = NULL, ...){
  p <- .data %>% ggplot(mapping = aes_string(x = x, y = y, fill = fill, color = color)) + geom_boxplot( ...) + theme_bw()
  return(p)
}

LegendDotSize <- function(size = 5){
  return(guides(colour = guide_legend(override.aes = list(size=size))))
}

## Function name: MyViolinPlot
## input:     .data,   data frame containing Distance and celltype information
##            
## output:    Violin plot.
MyViolinPlot <- function(.data = NULL, x = NULL, y = NULL, color_violin = NULL, fill_violin = NULL,
                         boxplot = FALSE, box_width = 0.2, jitter = FALSE, adjust = 1, ...){
  
  if(is.null(ylim)){
    ylim <- c(min(.data[[y]], na.rm = TRUE), max(.data[[y]], na.rm = TRUE))
  }
  
  p <- ggplot(data = .data, aes_string(x = x, y = y, color = color_violin, fill = fill_violin, ...)) +
    geom_violin(trim = TRUE, adjust = adjust, ...) +
    theme(axis.text.x = element_text(angle = 50, hjust = 1))
  
    p <- p + theme_bw() 
  
  return(p)
}

