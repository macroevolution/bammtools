#############################################################
#
#	transparentColor <- function(...)
#
#	Internal function allows for defining named colors with opacity
#	alpha = opacity
#

##' @title Define colors with transparency
##'
##' @description Converts a named color and opacity and returns the proper RGB
##'     code.
##'
##' @param namedColor A color name.
##' @param alpha A transparency value between 0 and 1, where 0 is fully
##'     transparent.
##'
##' @details This function is used internally by
##'     \code{\link{plotRateThroughTime}}.
##'
##' @return Returns the transparent color in RGB format.
##'
##' @author Pascal Title
##' @keywords manip
##' @export
transparentColor <- function(namedColor, alpha = 0.8) {
	res <- col2rgb(namedColor) / 255;
	return(rgb(red = res[1,], green = res[2,], blue = res[3,], alpha = alpha));
}
