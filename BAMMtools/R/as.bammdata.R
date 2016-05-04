
as.bammdata <- function(x, ...) {
	if ("bammdata" %in% class(x)) {
		return(x);
	}
	UseMethod("as.bammdata");
}
##' @export
as.bammdata.credibleshiftset <- function(x, ...) {
	obj <- x;
	class(obj) <- "bammdata";
	return(obj);
}
