#' Print letter on a figure
#'
#' \code{print.letter} Print a letter to label a panel on a figure
#'
#' @author T.A. Branch
#' @param label label to display on figure panel
#' @param xy x and y in terms of percentages of total axes, for plotting text
#' @param ... other arguments to 'text'
#' @importFrom graphics par text
#' @return text on a figure panel
#' @export

print_letter <- function(label="(a)",xy=c(0.1,0.925),...) { 
            tmp <- par("usr") 
            text.x <- tmp[1]+xy[1]*diff(tmp[1:2]) #x position, diff=difference 
            text.y <- tmp[3]+xy[2]*diff(tmp[3:4]) #y position 
            text(x=text.x, y=text.y, labels=label,...) 
}