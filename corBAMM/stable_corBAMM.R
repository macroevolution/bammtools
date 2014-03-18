## A file to put relatively stable working versions for corBAMM stuff

##########
# Implements the linear mixture error structure
#           psi*lambda*Em + (1-psi)*Eb
# \params
#    psi    :: mixture coefficient
#    lambda :: Pagel's lambda
# \constants
#    Em     :: phylogenetic correlation matrix
#    Eb     :: BAMM correlation matrix
# 
# Parameters are supplied via the 'value' argument. If value = numeric(0),
# the default is to estimate psi (initial value = 0.5) and to treat lambda
# as a fixed parameter equal to 1.
#
# if fixed = TRUE parameters are treated as fixed and not estimated. 
#
# Examples:
#    value = numeric(0),   fixed = FALSE        psi = 0.5(estimated), lambda = 1(fixed)   default
#    value = c(1,0),       fixed = TRUE         psi = 1(fixed), lambda = 0(fixed)
#    value = c(0,1),       fixed = TRUE         psi = 0(fixed), lambda = 1(fixed)
#    value = c(0.5, 0.75)  fixed = FALSE        psi = 0.5(estimated), lambda = 0.75(estimated)
##########

corBAMM <- function(value = numeric(0), ephy, form = ~1, fixed = FALSE)
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
    val <- as.vector(object);
    if (length(val) > 0) {
        if (val[1] < 0 || val[1] > 1)
            stop("psi must be between 0 and 1");
        if (val[1] == 0)
            val[1] <- .Machine$double.eps;
        val[1] <- log(val[1]) - log(1-val[1]);
        if (length(val) > 1) {
            if (val[2] < 0 || val[2] > 1)
                stop("lambda must be between 0 and 1");
            if (val[2] == 0)
                val[2] <- .Machine$double.eps;
            val[2] <- log(val[2]) - log(1-val[2]);
        }
    }
    else {
        val <- log(1);
    }
    oldAttr <- attributes(object);
    object <- val;
    attributes(object) <- oldAttr;
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
    psi <- 1/(1+exp(-object[1]));
    if (length(as.vector(object)) == 1) {
        Em <- psi*Em[index,index] + (1-psi)*Eb[index, index];
    }
    else {
        lambda <- 1/(1+exp(-object[2]));
        Em <- psi*lambda*Em[index,index] + (1-psi)*Eb[index, index];
    }
    diag(Em) <- 1;
    return (Em);
}

# This is basically coef
coef.corBAMM <- function (object, unconstrained = TRUE, ...) 
{
    lx <- function(x) 1/(1+exp(-x));
    if (unconstrained) {
        if (attr(object, "fixed")) 
            return(numeric(0))
        else return(lx(as.vector(object)));
    }
    aux <- as.vector(object);
    if (length(aux) > 1)
        names(aux) <- c("psi", "lambda")
    else
        names(aux) <- "psi";
    return (lx(aux));
}









