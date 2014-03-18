## A file to put relatively stable working versions for corBAMM stuff

##########
# Implements corBAMM error structure
##########
corBAMM <- function(value = 1, ephy, form = ~1) 
{
	if (class(ephy) != "bammdata")
	    stop("arg 'ephy' is not of class 'bammdata'");
	attr(value, "formula") <- form;
	attr(value, "fixed") <- TRUE;
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

##########
# Implements the linear mixture error structure:
# kappa*(corBrownian) + (1-kappa)*corBAMM
# kappa is estimated during model fitting
##########

corBAMM <- function(value = 0.5, ephy, form = ~1, fixed = FALSE)
{
    if (class(ephy) != "bammdata")
	    stop("arg 'ephy' is not of class 'bammdata'");
    attr(value, "formula") <- form;
	attr(value, "fixed") <- fixed;
	attr(value, "bammdata") <- ephy;
    class(value) <- c("corBAMM", "corStruct");
    return (value);    
}

Initialize.corBAMM <- function(object, data, ...)
{
    form <- formula(object);
    if (is.null(data)) data <- parent.frame();
    if (!is.null(getGroupsFormula(form))) {
        attr(object, "groups") <- getGroups(object, form, data = data);
        attr(object, "Dim") <- Dim(object, attr(object, "groups"));
    } 
    else {
        attr(object, "Dim") <- Dim(object, as.factor(rep(1, nrow(data))));
    }

    bammdata <- attr(object, "bammdata");
    phy <- as.phylo.bammdata(bammdata);
    attr(object, "Eb") <- getCohortMatrix(bammdata);
    attr(object, "Em") <- vcv.phylo(phy, corr = TRUE);
  
    if (nrow(data) != length(phy$tip.label))
        stop("number of observations and number of tips in the tree are not equal");
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
    
    return (object);   
}

corMatrix.corBAMM <- function(object, covariate = getCovariate(object), corr = TRUE, ...) 
{
    if (!("corBAMM" %in% class(object)))
        stop("arg 'object' not of class 'corBAMM'");
    if (!any(attr(object, "index")))
        stop("arg 'object' has not been initialized");
    if (is.null(Eb <- attr(object, "Eb")))
        Eb <- getCohortMatrix(attr(object, "bammdata"));
    if (is.null(Em <- attr(object, "Em")))
        Em <- vcv.phylo(as.phylo.bammdata(attr(object, "bammdata")), corr = TRUE);
    index <- attr(object, "index");
    if (object[1] < 1e-7) 
        object[1] <- 0.0
    else if (object[1] > 1)
        object[1] <- 1.0;
    kappa <- object[1];
    return(kappa*Em[index,index] + (1-kappa)*Eb[index, index]);
}
