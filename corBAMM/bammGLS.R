

logit <- function(p) return(log(p) - log(1-p));
invlogit <- function(x) return(1 / (exp(-x) + 1));

# fixed: vector specifying fixed parameters.

BAMMgls <- function (formula, data, Vb, Vp, lambda = 1, psi = 1, fixed=NULL, param.CI = 0.95, control = list(fnscale = -1)) 
{
	
#	formula <- y ~ x
#	data <- dff;
#	Vp <- vv;
#	Vb <- vb;
#	control = list(fnscale = -1)	
#    lambda <- 1;
#    psi <- 0;
#    fixed <- "psi"	


	# rescale cov matrices:    
    if (max(Vb) > 0){
      Vb <- Vb / max(Vb);  	
    }

    Vp <- Vp / max(Vp);
   	
	Dfun <- function(Cmat) {
        iCmat <- solve(Cmat, tol = .Machine$double.eps)
        svdCmat <- La.svd(iCmat)
        D <- svdCmat$u %*% diag(sqrt(svdCmat$d)) %*% t(svdCmat$v)
        return(t(D))
    }
    
    #if (!inherits(data, "comparative.data")) 
    #   stop("data is not a 'comparative' data object.")
    
    #dname <- deparse(substitute(data))
    call <- match.call()
    miss <- model.frame(formula, data, na.action = na.pass)
    miss.na <- apply(miss, 1, function(X) (any(is.na(X))))
    if (any(miss.na)) {
        miss.names <- data$phy$tip.label[miss.na]
        data <- data[-which(miss.na), ]
    }
    m <- model.frame(formula, data)
    y <- m[, 1]
    x <- model.matrix(formula, m)
    k <- ncol(x)
    namey <- names(m)[1]
    xVar <- apply(x, 2, var)[-1]
    badCols <- xVar < .Machine$double.eps
    if (any(badCols)) 
        stop("Model matrix contains columns with zero variance: ", 
            paste(names(xVar)[badCols], collapse = ", "))
    
#### Compute cohort matrix and phylo vcv matrix here if not specified:    
        
    nm <- rownames(data)
    n <- nrow(data)
    if (!is.null(param.CI)) {
        if (!is.numeric(param.CI) || param.CI <= 0 || param.CI > 
            1) 
            stop("param.CI is not a number between 0 and 1.")
    }
     
   # parVals <- list(lambda = lambda, psi = psi)
 
	if (length(fixed) < 2) {
    	
    	## Test here to make sure valid initial lambda and psi are specified.
    	
 		parVals <- c(lambda, psi)
        names(parVals) <- c("lambda", "psi")
	
		optimPar <- parVals[setdiff(names(parVals), fixed)];
        fixedPar <- parVals[fixed]

        
        bounds.low <- rep(0, length(optimPar));
        bounds.up <- rep(1, length(optimPar));
        
        
        #BAMMgls.likelihood(optimPar, fixedPar, y, x, Vb, Vp, optim.output=T)
        
        optim.param.vals <- optim(optimPar, fn = BAMMgls.likelihood, 
            method = "L-BFGS-B", control = control, Vb = Vb, Vp=Vp, y = y, x = x, fixedPar = fixedPar, 
            optim.output = TRUE, logit.params=F, lower=bounds.low, upper=bounds.up)
        if (optim.param.vals$convergence != "0") {
            stop("Problem with optim:", optim.param.vals$convergence, 
                optim.param.vals$message)
        }
        
        #pars.transformed <- invlogit(optim.param.vals$par);
        
        fixedPar <- c(optim.param.vals$par, fixedPar)
        fixedPar <- fixedPar[c("lambda", "psi")]
        
    } else {
        fixedPar <- c(lambda, psi);
        names(fixedPar) <- c("lambda", "psi")
        
    }
    
    ll <- BAMMgls.likelihood(optimPar = NULL, fixedPar = fixedPar, 
        y, x, Vb = Vb, Vp = Vp, optim.output = FALSE)
        
    log.lik <- ll$ll
    Vt <- BAMMgls.mxTransform(Vp = Vp, Vb = Vb, fixedPar)
    aic <- -2 * log.lik + 2 * k
    aicc <- -2 * log.lik + 2 * k + ((2 * k * (k + 1))/(n - k - 
        1))
    coeffs <- ll$mu
    names(coeffs) <- colnames(x)
    varNames <- names(m)
    pred <- x %*% ll$mu
    res <- y - pred
    
    
    ## Check: What is Dfun doing???
    D <- Dfun(Vt)
    pres <- D %*% res
    fm <- list(coef = coeffs, aic = aic, log.lik = log.lik)
    
    #######
    ####### checked through here...
    
    RMS <- ll$s2
    RSSQ <- ll$s2 * (n - k)
    xdummy <- matrix(rep(1, length(y)))
    nullMod <- BAMMgls.likelihood(optimPar = NULL, fixedPar = fixedPar, 
        y, xdummy, Vb = Vb, Vp = Vp, optim.output = FALSE)
    NMS <- nullMod$s2
    NSSQ <- nullMod$s2 * (n - 1)
    errMat <- t(x) %*% solve(Vt) %*% x
    errMat <- solve(errMat) * RMS[1]
    sterr <- diag(errMat)
    sterr <- sqrt(sterr)
    
    
    
    RET <- list(model = fm, formula = formula, call = call, RMS = RMS, 
        NMS = NMS, NSSQ = NSSQ[1], RSSQ = RSSQ[1], aic = aic, 
        aicc = aicc, n = n, logLik = log.lik, k = k, sterr = sterr, fitted = pred, 
        residuals = res, phyres = pres, x = x, data = data, varNames = varNames, 
        y = y, param = fixedPar, namey = namey, 
        Vt = Vt)
    class(RET) <- "BAMMgls"
    if (any(miss.na)) {
        RET$na.action <- structure(which(miss.na), class = "omit", 
            .Names = miss.names)
    }
    
    # # confidence intervals not yet incorporated.
    #
    #
    # if (!is.null(param.CI) && any(mlVals)) {
        # param.CI.list <- list(kappa = NULL, lambda = NULL, delta = NULL)
        # mlNames <- names(mlVals)[which(mlVals)]
        # for (param in mlNames) {
            # param.CI.list[[param]] <- pgls.confint(RET, param, 
                # param.CI)
        # }
        # RET$param.CI <- param.CI.list
    # }
    return(RET)
}



BAMMgls.likelihood <- function (optimPar, fixedPar, y, x, Vb, Vp, optim.output = TRUE, names.optim = NULL, logit.params=FALSE) 
{
    get.coeffs <- function(Y, iV, X) {
        xVix <- crossprod(X, iV %*% X)
        xViy <- crossprod(X, iV %*% Y)
        mu <- solve(xVix, tol = .Machine$double.eps) %*% xViy
        return(mu)
    }
    est.var <- function(y, iV, x, mu) {
        e <- y - x %*% mu
        s2 <- crossprod(e, iV %*% e)
        n <- length(y)
        k <- length(x[1, ])
        return(s2/(n - k))
    }
    if (!is.null(names.optim)) 
        names(optimPar) <- names.optim
    
    if (logit.params){
    	optimPar <- invlogit(optimPar)
    }
    
    allPar <- c(optimPar, fixedPar)
 
	V <- BAMMgls.mxTransform(Vp = Vp, Vb = Vb, allPar)
 
 	iV <- solve(V, tol = .Machine$double.eps)
    mu <- get.coeffs(y, iV, x)
    s2 <- est.var(y, iV, x, mu)
    n <- nrow(x)
    k <- ncol(x)
    logDetV <- determinant(V, logarithm = TRUE)$modulus[1]

    ll <- -n/2 * log(2 * pi) - n/2 * log((n - k) * s2/n) - logDetV/2 - 
        n/2

    if (optim.output) 
        return(ll)
    else return(list(ll = ll, mu = mu, s2 = s2))
}


BAMMgls.mxTransform <- function (Vp, Vb, fixedPar) 
{
 
	d <- diag(Vp);
	Vp <- Vp*fixedPar['lambda'];
	diag(Vp) <- d;
	
	Vp[Vb >= 0.95] <- 0;
	Vb[Vb < 0.95] <- 0;
	V <- Vp + Vb;
	
	#V <- (1 - fixedPar['psi'])*Vp + fixedPar['psi'] * Vb;
    attr(V, "mxTransform") <- fixedPar
    return(V)
}


summary.BAMMgls <- function (object, ...) 
{
    ans <- list(call = object$call)
    class(ans) <- "summary.BAMMgls"
    p <- object$k
    n <- object$n
    rdf <- n - p
    ans$df <- c(p, rdf)
    r <- object$phyres
    rss <- object$RSSQ
    resvar <- rss/rdf
    ans$sigma <- sqrt(resvar)
    ans$residuals <- r
    cf <- object$model$coef
    se <- object$sterr
    t <- cf/se
    coef <- cbind(cf, se, t, 2 * (1 - pt(abs(t), rdf)))
    colnames(coef) <- c("Estimate", "Std. Error", "t value", 
        "Pr(>|t|)")
    ans$coefficients <- coef
    ans$param <- object$param
    ans$mlVals <- object$mlVals
   #if (!is.null(object$param.CI)) 
   #    ans$param.CI <- object$param.CI
    if (!is.null(object$na.action)) 
        ans$na.action <- object$na.action
    ans$fstatistic <- c(value = ((object$NSSQ - object$RSSQ)/object$RMS)/(object$k - 
        1), numdf = p - 1, dendf = rdf)
    ans$r.squared <- (object$NSSQ - object$RSSQ)/object$NSSQ
    ans$adj.r.squared <- (object$NMS - object$RMS)/object$NMS
    ans$logLik <- as.numeric(object$logLik);
    return(ans)
}




