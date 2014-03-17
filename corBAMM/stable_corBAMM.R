corBAMM <- function(value = 1, ephy, form = ~1, fixed = TRUE) 
{
	if (class(ephy) != "bammdata")
	    stop("arg 'ephy' is not of class 'bammdata'");
	attr(value, "formula") <- form;
	attr(value, "fixed") <- fixed;
	attr(value, "bammdata") <- ephy;
	class(value) <- c("corBAMM", "corStruct");
	return(value);
}

Initialize.corBAMM <- function(object, data, ...)
{
    form <- formula(object);
    if (!is.null(getGroupsFormula(form))) {
        attr(object, "groups") <- getGroups(object, form, data = data);
        attr(object, "Dim") <- Dim(object, attr(object, "groups"));
    } 
    else { 
        attr(object, "Dim") <- Dim(object, as.factor(rep(1, nrow(data))));
    }
    attr(object, "covariate") <- getCovariate(object, data = data);
    attr(object, "Eb") <- getCohortMatrix(attr(object, "bammdata"));
    
    phy <- as.phylo.bammdata(attr(object, "bammdata"));
    if (is.null(data)) data <- parent.frame();
    if (nrow(data) != length(phy$tip.label))
        stop("number of observations and number of tips in the tree are not equal.");
    if (is.null(rownames(data))) {
        warning("No rownames supplied in data frame, data taken to be in the same order as in tree");
        attr(object, "index") <- 1:dim(data)[1];
    } 
    else {
        index <- match(rownames(data), phy$tip.label);
        if (any(is.na(index))) {
            warning("Rownames in data frame do not match tree tip names; data taken to be in the same order as in tree");
            attr(object, "index") <- 1:dim(data)[1];
        } 
        else {
            attr(object, "index") <- index;
        }
    }
    return(object);	
}

corMatrix.corBAMM <- function(object, covariate = getCovariate(object), corr = TRUE, ...)
{
    if (!("corBAMM" %in% class(object)))
        stop("arg 'object' not of class 'corBAMM'");
    if (!any(attr(object, "index")))
        stop("arg 'object' has not been initialized");
    if (is.null(Eb <- attr(object, "Eb")))
        Eb <- getCohortMatrix(attr(object, "bammdata"));
    index <- attr(object, "index");
    return(Eb[index, index]);
}

coef.corBAMM <- function(object, unconstrained = TRUE, ...)
{
    if (!("corBAMM" %in% class(object)))
        stop('object is not of class "corBAMM"')
    return(numeric(0));
}
