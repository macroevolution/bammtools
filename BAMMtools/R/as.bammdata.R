
as.bammdata <- function(x, ...) {
	if (length(class(x)) == 1 && inherits(x, "bammdata")) {
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
