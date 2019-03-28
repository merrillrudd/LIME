#' Theme for ggplot2 
#'
#' \code{mytheme} Preferred theme for ggplot2
#' @author Darcy Webber

#' @param base_size default text size for figure, default=14
#' @param base_family font type for figure, default=""
#' 
#' @return theme for use with figures from ggplot2
#' @export
mytheme <- function (base_size = 14, base_family = "") 
{
    theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(axis.title.x = element_text(margin = margin(10,0,0,0)),
          #axis.title.x = element_text(vjust = -1.5),
          #axis.title.y = element_text(margin = margin(0,20,0,0)),
          #axis.title.y = element_text(vjust = -0.1),
          axis.text = element_text(size = rel(0.8)),
          axis.ticks = element_line(colour = "black"), 
          legend.key = element_rect(colour = "grey80"),
          panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour = "grey50"),
          panel.grid.major = element_line(colour = "grey90", size = 0.2),
          panel.grid.minor = element_line(colour = "grey98", size = 0.5),
          strip.background = element_rect(fill = "grey80", colour = "grey50", size = 0.2))
}