#' @title Evaluate Trait Predictability via Cross Validation
#' @description The cv function evaluates trait predictability based on eight GS methods via k-fold cross validation.
#' The trait predictability is defined as the squared Pearson correlation coefficient between the observed and the predicted trait values.
#' @param fix a design matrix of the fixed effects.
#' @param gena a matrix (n x m) of additive genotypes for the training population.
#' @param gend a matrix (n x m) of domiance genotypes for the training population. Default is NULL.
#' @param y a vector(n x 1) of the phenotypic values.
#' @param nfold the number of folds. Default is 5.
#' @param nTimes the number of independent replicates for the cross-validation. Default is 1.
#' @param CPU the number of CPU.
#' @param method eight GS methods including "GBLUP", "BayesB", "RKHS", "PLS", "LASSO", "EN", "XGBOOST", "RF".
#' Users may select one of these methods or all of them simultaneously with "ALL". Default is "GBLUP".
#' @param seed the random number. Default is 123.
#' @return Trait predictability
#' @examples
#' \donttest{
#' ## load example data from hypred package
#' data(hybrid_phe)
#' data(input_geno)
#'
#' ## convert original genotype
#' inbred_gen <- convertgen(input_geno, type = "hmp2")
#'
#' ##additive model infer the additive and dominance genotypes of hybrids
#' gena <- infergen(inbred_gen, hybrid_phe)$add
#' gend <- infergen(inbred_gen, hybrid_phe)$dom
#'
#' ##additive model
#' R2<-cv(fix=NULL,gena,gend=NULL,y=hybrid_phe[,3],method ="GBLUP",nfold=5,nTimes=1,seed=123,CPU=1)
#'
#' ##additive-dominance model
#' R2<-cv(fix=NULL,gena,gend,y=hybrid_phe[,3],method ="GBLUP",nfold=5,nTimes=1,seed=123,CPU=1)}
#' @export
cv <- function(fix = NULL, gena, gend = NULL, y, method = "GBLUP", nfold = 5, nTimes = 1, seed = 123, CPU = 1) {
    cv.GBLUP <- function(fix, gena, gend, y, nfold, seed, CPU) {
        print("Predict by GBLUP...")
        set.seed(seed + time)
        foldid <- sample(rep(1:nfold, ceiling(n/nfold))[1:n])
        if (is.null(gend)) {
            kk <- list(kin(gena))
        } else {
            kk <- list(kin(gena), kin(gend))
        }
        n <- length(y)
        y <- as.matrix(y)
        if (is.null(fix)) {
            fix <- matrix(1, n, 1)
        } else {
            fix <- as.matrix(fix)
        }
        yobs <- NULL
        yhat <- NULL
        fold <- NULL
        id <- NULL
        k11 <- list()
        k21 <- list()
        g <- length(kk)
        cl.cores <- detectCores()
        if (cl.cores <= 2 || CPU <= 1) {
            cl.cores <- 1
        } else if (cl.cores > 2) {
            if (cl.cores > 10) {
                cl.cores <- 10
            } else {
                cl.cores <- detectCores() - 1
            }
        }
        cl <- makeCluster(cl.cores)
        registerDoParallel(cl)
        res <- foreach(k = 1:nfold, .multicombine = TRUE, .combine = "rbind", .packages = c("predhy")) %dopar% {
            i1 <- which(foldid != k)
            i2 <- which(foldid == k)
            x1 <- fix[i1, , drop = F]
            y1 <- y[i1, , drop = F]
            x2 <- fix[i2, , drop = F]
            y2 <- y[i2, , drop = F]
            for (i in 1:g) {
                k11[[i]] <- kk[[i]][i1, i1]
                k21[[i]] <- kk[[i]][i2, i1]
            }
            parm <- mixed(fix = x1, y = y1, kk = k11)
            G21 <- 0
            for (i in 1:g) {
                G21 <- G21 + k21[[i]] * parm$var[i]
            }
            v_i <- parm$v_i
            beta <- parm$beta
            y3 <- x2 %*% beta + G21 %*% v_i %*% (y1 - x1 %*% beta)
            yobs <- y2
            yhat <- y3
            data.frame(yobs, yhat)
        }
        stopCluster(cl)
        R2 <- cor(res$yobs, res$yhat)^2
        print("Predict by GBLUP...ended.")
        return(R2)
    }
    cv.bayesb <- function(fix, gena, gend, y, nfold, seed, CPU) {
        print("Predict by BayesB...")
        set.seed(seed + time)
        foldid <- sample(rep(1:nfold, ceiling(n/nfold))[1:n])
        if (!requireNamespace("BGLR", quietly = TRUE)) {
            stop("BGLR needed for this function to work. Please install it.", call. = FALSE)
        }
        # library(BGLR)
        if (is.null(gend)) {
            x <- gena
        } else {
            x <- cbind(gena, gend)
        }

        if (is.null(fix)) {
            fix <- matrix(1, n, 1)
        } else {
            fix <- as.matrix(fix)
        }
        yhat <- NULL
        yobs <- NULL
        cl.cores <- detectCores()
        if (cl.cores <= 2 || CPU <= 1) {
            cl.cores <- 1
        } else if (cl.cores > 2) {
            if (cl.cores > 10) {
                cl.cores <- 10
            } else {
                cl.cores <- detectCores() - 1
            }
        }
        cl <- makeCluster(cl.cores)
        registerDoParallel(cl)
        res <- foreach(k = 1:nfold, .multicombine = TRUE, .combine = "rbind", .packages = c("BGLR")) %dopar% {
            yNa <- y
            whichNa <- which(foldid == k)
            yNa[whichNa] <- NA
            eta <- list(list(X = x, model = "BayesB"))
            fm <- BGLR(y = yNa, ETA = eta, verbose = F)
            yhat <- fm$yHat[whichNa]
            yobs <- y[whichNa]
            data.frame(yobs, yhat)
        }
        stopCluster(cl)
        pr2 <- cor(res$yobs, res$yhat)^2
        print("Predict by BayesB...ended.")
        return(pr2)
    }
    cv.rkhsmk <- function(fix, gena, gend, y, nfold, seed, CPU) {
        print("Predict by RKHS...")
        set.seed(seed + time)
        foldid <- sample(rep(1:nfold, ceiling(n/nfold))[1:n])
        # library(BGLR)
        if (is.null(gend)) {
            x <- gena
        } else {
            x <- cbind(gena, gend)
        }
        if (is.null(fix)) {
            fix <- matrix(1, n, 1)
        } else {
            fix <- as.matrix(fix)
        }
        X <- scale(x, center = TRUE, scale = TRUE)
        D <- (as.matrix(dist(X, method = "euclidean"))^2)/ncol(X)
        h <- 0.5 * c(1/5, 1, 5)
        ETA <- list(list(K = exp(-h[1] * D), model = "RKHS"), list(K = exp(-h[2] * D), model = "RKHS"), list(K = exp(-h[3] *
            D), model = "RKHS"), list(X = fix, model = "FIXED"))
        yhat <- NULL
        yobs <- NULL
        cl.cores <- detectCores()
        if (cl.cores <= 2) {
            cl.cores <- 1
        } else if (cl.cores > 2 || CPU <= 1) {
            if (cl.cores > 10) {
                cl.cores <- 10
            } else {
                cl.cores <- detectCores() - 1
            }
        }
        cl <- makeCluster(cl.cores)
        registerDoParallel(cl)
        res <- foreach(k = 1:nfold, .multicombine = TRUE, .combine = "rbind", .packages = c("BGLR")) %dopar% {
            yNa <- y
            whichNa <- which(foldid == k)
            yNa[whichNa] <- NA
            fm <- BGLR(y = yNa, ETA = ETA, verbose = F)
            yhat <- fm$yHat[whichNa]
            yobs <- y[whichNa]
            data.frame(yobs, yhat)
        }
        stopCluster(cl)
        pr2 <- cor(res$yobs, res$yhat)^2
        print("Predict by RKHS...ended.")
        return(pr2)
    }
    cv.pls <- function(fix, gena, gend, y, nfold, seed, CPU) {
        print("Predict by PLS...")
        set.seed(seed + time)
        foldid <- sample(rep(1:nfold, ceiling(n/nfold))[1:n])
        if (!requireNamespace("pls", quietly = TRUE)) {
            stop("pls needed for this function to work. Please install it.", call. = FALSE)
        }
        if (is.null(gend)) {
            x <- gena
        } else {
            x <- cbind(gena, gend)
        }
        z <- cbind(fix, x)
        # library(pls)
        yp <- NULL
        yo <- NULL
        cl.cores <- detectCores()
        if (cl.cores <= 2 || CPU <= 1) {
            cl.cores <- 1
        } else if (cl.cores > 2) {
            if (cl.cores > 10) {
                cl.cores <- 10
            } else {
                cl.cores <- detectCores() - 1
            }
        }
        cl <- makeCluster(cl.cores)
        registerDoParallel(cl)
        res <- foreach(k = 1:nfold, .multicombine = TRUE, .combine = "rbind", .packages = c("pls")) %dopar% {
            # for (k in 1:nfold) {
            id1 <- which(foldid != k)
            id2 <- which(foldid == k)
            x1 <- z[id1, ]
            x2 <- z[id2, ]
            y1 <- y[id1]
            y2 <- y[id2]
            pls.fit <- plsr(y1 ~ x1, ncomp = 5, validation = "CV")
            # plot(RMSEP(pls.fit), legendpos = 'topright')
            nn <- as.numeric(which.min(tt <- RMSEP(pls.fit)$val[1, , ][-1]))
            ## predicting
            yhat <- predict(pls.fit, newdata = x2, ncomp = nn)
            yp <- as.numeric(yhat)
            yo <- y2
            data.frame(yp, yo)
        }
        stopCluster(cl)
        pr2 <- cor(res$yo, res$yp)^2
        print("Predict by PLS...ended.")
        return(pr2)
    }
    cv.lasso <- function(fix, gena, gend, y, nfold, seed, CPU) {
        print("Predict by LASSO...")
        set.seed(seed + time)
        foldid <- sample(rep(1:nfold, ceiling(n/nfold))[1:n])
        if (!requireNamespace("glmnet", quietly = TRUE)) {
            stop("glmnet needed for this function to work. Please install it.", call. = FALSE)
        }
        # library(glmnet)
        if (is.null(gend)) {
            x <- gena
        } else {
            x <- cbind(gena, gend)
        }
        z <- cbind(fix, x)
        yyhat <- NULL
        yyobs <- NULL
        cl.cores <- detectCores()
        if (cl.cores <= 2 || CPU <= 1) {
            cl.cores <- 1
        } else if (cl.cores > 2) {
            if (cl.cores > 10) {
                cl.cores <- 10
            } else {
                cl.cores <- detectCores() - 1
            }
        }
        cl <- makeCluster(cl.cores)
        registerDoParallel(cl)
        res <- foreach(k = 1:nfold, .multicombine = TRUE, .combine = "rbind", .packages = c("glmnet")) %dopar% {
            # for (k in 1:nfold) {
            iindex <- which(foldid == k)
            xx <- z[-iindex, ]
            yy <- y[-iindex]
            fit0 <- cv.glmnet(x = xx, y = yy)
            lambda <- fit0$lambda.min
            ffit <- glmnet(x = xx, y = yy, lambda = lambda)
            tmp <- predict(ffit, newx = z[iindex, ])
            yyhat <- as.numeric(tmp)
            yyobs <- as.numeric(as.matrix(y[iindex]))
            data.frame(yyhat, yyobs)
        }
        stopCluster(cl)
        pr2 <- cor(res$yyobs, res$yyhat)^2
        print("Predict by LASSO...ended.")
        return(pr2)
    }
    cv.elas <- function(fix, gena, gend, y, nfold, seed, CPU) {
        print("Predict by EN...")
        # library(glmnet)
        set.seed(seed + time)
        foldid <- sample(rep(1:nfold, ceiling(n/nfold))[1:n])
        if (is.null(gend)) {
            x <- gena
        } else {
            x <- cbind(gena, gend)
        }
        z <- cbind(fix, x)
        yyhat <- NULL
        yyobs <- NULL
        cl.cores <- detectCores()
        if (cl.cores <= 2 || CPU <= 1) {
            cl.cores <- 1
        } else if (cl.cores > 2) {
            if (cl.cores > 10) {
                cl.cores <- 10
            } else {
                cl.cores <- detectCores() - 1
            }
        }
        cl <- makeCluster(cl.cores)
        registerDoParallel(cl)
        res <- foreach(k = 1:nfold, .multicombine = TRUE, .combine = "rbind", .packages = c("glmnet")) %dopar% {
            iindex <- which(foldid == k)
            xx <- z[-iindex, ]
            yy <- y[-iindex]
            fit0 <- cv.glmnet(x = xx, y = yy, alpha = 0.5)
            lambda <- fit0$lambda.min
            ffit <- glmnet(x = xx, y = yy, lambda = lambda, alpha = 0.5)
            tmp <- predict(ffit, newx = z[iindex, ])
            yyhat <- as.numeric(tmp)
            yyobs <- as.numeric(as.matrix(y[iindex]))
            data.frame(yyhat, yyobs)
        }
        stopCluster(cl)
        pr2 <- cor(res$yyobs, res$yyhat)^2
        print("Predict by EN...ended.")
        return(pr2)
    }
    cv.xgboost <- function(fix, gena, gend, y, nfold, seed, CPU) {
        print("Predict by XGBOOST...")
        set.seed(seed + time)
        foldid <- sample(rep(1:nfold, ceiling(n/nfold))[1:n])
        if (!requireNamespace("xgboost", quietly = TRUE)) {
            stop("xgboost needed for this function to work. Please install it.", call. = FALSE)
        }
        # library(xgboost)
        if (is.null(gend)) {
            x <- gena
        } else {
            x <- cbind(gena, gend)
        }
        z <- cbind(fix, x)
        yp <- NULL
        yo <- NULL
        cl.cores <- detectCores()
        if (cl.cores <= 2 || CPU <= 1) {
            cl.cores <- 1
        } else if (cl.cores > 2) {
            if (cl.cores > 10) {
                cl.cores <- 10
            } else {
                cl.cores <- detectCores() - 1
            }
        }
        cl <- makeCluster(cl.cores)
        registerDoParallel(cl)
        res <- foreach(k = 1:nfold, .multicombine = TRUE, .combine = "rbind", .packages = c("xgboost")) %dopar% {
            id1 <- which(foldid != k)
            id2 <- which(foldid == k)
            x1 <- z[id1, ]
            x2 <- z[id2, ]
            y1 <- y[id1]
            y2 <- y[id2]
            xg <- xgboost(data = x1, label = y1, nrounds = 1000, eta = 0.07)
            ## predicting
            yhat <- predict(xg, x2)
            yp <- yhat
            yo <- y2
            data.frame(yp, yo)
        }
        stopCluster(cl)
        pr2 <- cor(res$yo, res$yp)^2
        print("Predict by XGBOOST...ended.")
        return(pr2)
    }
    cv.rf <- function(fix, gena, gend, y, nfold, seed, CPU) {
        print("Predict by RF...")
        set.seed(seed + time)
        foldid <- sample(rep(1:nfold, ceiling(n/nfold))[1:n])
        if (!requireNamespace("randomForest", quietly = TRUE)) {
            stop("randomForest needed for this function to work. Please install it.", call. = FALSE)
        }
        # library(randomForest)
        if (is.null(gend)) {
            x <- gena
        } else {
            x <- cbind(gena, gend)
        }
        z <- cbind(fix, x)
        yobs <- NULL
        yhat <- NULL
        cl.cores <- detectCores()
        if (cl.cores <= 2 || CPU <= 1) {
            cl.cores <- 1
        } else if (cl.cores > 2) {
            if (cl.cores > 10) {
                cl.cores <- 10
            } else {
                cl.cores <- detectCores() - 1
            }
        }
        cl <- makeCluster(cl.cores)
        registerDoParallel(cl)
        res <- foreach(k = 1:nfold, .multicombine = TRUE, .combine = "rbind", .packages = c("randomForest")) %dopar%
            {
                id1 <- which(foldid != k)
                id2 <- which(foldid == k)
                x1 <- z[id1, ]
                x2 <- z[id2, ]
                y1 <- y[id1]
                y2 <- y[id2]
                fit <- randomForest(x = x1, y = y1, ntree = 500)
                ypred <- predict(fit, x2)
                yobs <- y2
                yhat <- ypred
                data.frame(yhat, yobs)
            }
        stopCluster(cl)
        pr2 <- cor(res$yobs, res$yhat)^2
        print("Predict by RF...ended.")
        return(pr2)
    }
    cvres <- NULL
    cvres <- list(cvres)
    k <- NULL
    if (is.null(gend)) {
        print("additive model")
    } else {
        print("additive-dominance model")
    }
    for (time in 1:nTimes) {
        cat(time)
        n <- length(y)
        if ((method == "GBLUP") | (method == "BayesB") | (method == "RKHS") | (method == "PLS") | (method == "LASSO") |
            (method == "EN") | (method == "XGBOOST") | (method == "RF") | (method == "ALL")) {
            if (method == "GBLUP") {
                cvres[[time]] <- cv.GBLUP(fix, gena, gend, y, nfold, seed, CPU)
            }
            if (method == "BayesB") {
                cvres[[time]] <- cv.bayesb(fix, gena, gend, y, nfold, seed, CPU)
            }
            if (method == "RKHS") {
                cvres[[time]] <- cv.rkhsmk(fix, gena, gend, y, nfold, seed, CPU)
            }
            if (method == "PLS") {
                cvres[[time]] <- cv.pls(fix, gena, gend, y, nfold, seed, CPU)
            }
            if (method == "LASSO") {
                cvres[[time]] <- cv.lasso(fix, gena, gend, y, nfold, seed, CPU)
            }
            if (method == "EN") {
                cvres[[time]] <- cv.elas(fix, gena, gend, y, nfold, seed, CPU)
            }
            if (method == "XGBOOST") {
                cvres[[time]] <- cv.xgboost(fix, gena, gend, y, nfold, seed, CPU)
            }
            if (method == "RF") {
                cvres[[time]] <- cv.rf(fix, gena, gend, y, nfold, seed, CPU)
            }
            if (method == "ALL") {
                print("Predict by ALL...")
				GBLUP <- cv.GBLUP(fix, gena, gend, y, nfold, seed, CPU)
				BayesB <- cv.bayesb(fix, gena, gend, y, nfold, seed, CPU)
                RKHS <- cv.rkhsmk(fix, gena, gend, y, nfold, seed, CPU)
                PLS <- cv.pls(fix, gena, gend, y, nfold, seed, CPU)
                LASSO <- cv.lasso(fix, gena, gend, y, nfold, seed, CPU)
                EN <- cv.elas(fix, gena, gend, y, nfold, seed, CPU)
                XGBOOST <- cv.xgboost(fix, gena, gend, y, nfold, seed, CPU)
                RF <- cv.rf(fix, gena, gend, y, nfold, seed, CPU)
                cvres[[time]] <- cbind(GBLUP, BayesB, RKHS, PLS, LASSO, EN, XGBOOST, RF)
                print("Predict by ALL...ended.")
            }
        } else {
            stop("Please choose a predict method: GBLUP, BayesB, RKHS, PLS, LASSO, EN, XGBOOST, RF.")
        }  # End of all predict methods
    }
    res <- do.call(rbind, cvres)
    return(res)
}
