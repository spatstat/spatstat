# hackglmm.R
#  $Revision: 1.4 $ $Date: 2016/02/11 00:48:01 $

hackglmmPQL <- 
function (fixed, random, family, data, correlation, weights,
    control, niter = 10, verbose = TRUE, subset, ..., reltol=1e-3)
{
    if (is.character(family))
        family <- get(family)
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    m <- mcall <- Call <- match.call()
    nm <- names(m)[-1]
    keep <- is.element(nm, c("weights", "data", "subset", "na.action"))
    for (i in nm[!keep]) m[[i]] <- NULL
    allvars <- if (is.list(random))
        allvars <- c(all.vars(fixed), names(random), unlist(lapply(random,
            function(x) all.vars(formula(x)))))
    else c(all.vars(fixed), all.vars(random))
    Terms <- if (missing(data))
        terms(fixed)
    else terms(fixed, data = data)
    off <- attr(Terms, "offset")
    if (length(off <- attr(Terms, "offset")))
        allvars <- c(allvars, as.character(attr(Terms, "variables"))[off +
            1])
    Call$fixed <- eval(fixed)
    Call$random <- eval(random)
    m$formula <- as.formula(paste("~", paste(allvars, collapse = "+")))
    environment(m$formula) <- environment(fixed)
    m$drop.unused.levels <- TRUE
    m[[1]] <- as.name("model.frame")
    mf <- eval.parent(m)
    off <- model.offset(mf)
    if (is.null(off))
        off <- 0
    w <- model.weights(mf)
    if (is.null(w))
        w <- rep(1, nrow(mf))
    wts <- mf$wts <- w
    if(missing(subset)) 
      fit0 <- glm(formula = fixed, family = family, data = mf,
                  weights = wts, ...)
    else {
    # hack to get around peculiar problem with `subset' argument
      glmmsubset <- eval(expression(subset), data)
      if(length(glmmsubset) != nrow(mf)) {
        if(sum(glmmsubset) != nrow(mf))
          stop("Internal error: subset vector is wrong length")
        message("(Fixing subset index..)")
        glmmsubset <- glmmsubset[glmmsubset]
      }
      mf$glmmsubset <- glmmsubset
      fit0 <- glm(formula = fixed, family = family, data = mf,
                  weights = wts, subset=glmmsubset, ...)
    } 
    w <- fit0$prior.weights
    eta <- fit0$linear.predictor
    zz <- eta + fit0$residuals - off
    wz <- fit0$weights
    fam <- family
    nm <- names(mcall)[-1]
    keep <- is.element(nm, c("fixed", "random", "data", "subset",
        "na.action", "control"))
    for (i in nm[!keep]) mcall[[i]] <- NULL
    fixed[[2]] <- quote(zz)
    mcall[["fixed"]] <- fixed
    mcall[[1]] <- as.name("lme")
    mcall$random <- random
    mcall$method <- "ML"
    if (!missing(correlation))
        mcall$correlation <- correlation
    mcall$weights <- quote(varFixed(~invwt))
    mf$zz <- zz
    mf$invwt <- 1/wz
    mcall$data <- mf
    for (i in 1:niter) {
        if (verbose)
            cat("iteration", i, "\n")
        fit <- eval(mcall)
        etaold <- eta
        eta <- fitted(fit) + off
        if (sum((eta - etaold)^2) < (reltol^2) * sum(eta^2))
            break
        mu <- fam$linkinv(eta)
        mu.eta.val <- fam$mu.eta(eta)
        mf$zz <- eta + (fit0$y - mu)/mu.eta.val - off
        wz <- w * mu.eta.val^2/fam$variance(mu)
        mf$invwt <- 1/wz
        mcall$data <- mf
    }
    attributes(fit$logLik) <- NULL
    fit$call <- Call
    fit$family <- family
    fit$logLik <- as.numeric(NA)
    oldClass(fit) <- c("glmmPQL", oldClass(fit))
    fit
}

