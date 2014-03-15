corBAMM <- function(value = 1, ephy, form = ~1, fixed = FALSE) 
{
	if (class(ephy) != "bammdata")
	    stop("arg 'ephy' is not of class 'bammdata'");
	if (value < 0 || value > 1) 
	    stop("kappa must be between 0 and 1");
	attr(value, "formula") <- form;
	attr(value, "fixed") <- fixed;
	attr(value, "bammdata") <- ephy;
	class(value) <- c("corBAMM", "corStruct");
	return(value);
}

Initialize.corBAMM <- function(object, data, ...)
{
    ## The same as in Initialize corStruct:
    form <- formula(object);
    ## Obtaining the group information, if any
    if (!is.null(getGroupsFormula(form))) {
        attr(object, "groups") <- getGroups(object, form, data = data);
        attr(object, "Dim") <- Dim(object, attr(object, "groups"));
    } 
    else { # no groups
        attr(object, "Dim") <- Dim(object, as.factor(rep(1, nrow(data))));
    }
    ## Obtaining the covariate(s)
    attr(object, "covariate") <- getCovariate(object, data = data);
    
    ## The same as in Initialize corPhyl but with addition of as.phylo.bammdata
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
    if (!corr)
        stop("arg 'corr' must be 'TRUE' for class 'corBAMM'");
    ephy <- attr(object, "bammdata");
    phy <- as.phylo.bammdata(ephy);
    vm <- vcv.phylo(phy, corr = TRUE);
    vb <- getCohortMatrix(ephy);
    kap <- object[1];
    vm <- kap*vm + (1-kap)*vb;
    index <- attr(object, "index");
    return(vm[index, index]);
}

coef.corBAMM <- function(object, unconstrained = TRUE, ...)
{
    if (unconstrained) {
        if (attr(object, "fixed")) { 
            return(numeric(0));
        }
        else 
            return(object[1]);
    }
    aux <- object[1];
    names(aux) <- "kappa";
    return(aux);
}
