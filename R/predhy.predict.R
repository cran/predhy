#' @title Predict the Performance of Hybrids
#' @description Predict all potential crosses of a given set of parents using a subset of crosses as the training sample.
#' @param inbred_gen a matrix for genotypes of parental lines in numeric format, coded as 1, 0 and -1. The row.names of inbred_gen must be provied. It can be obtained from the original genotype using  \code{\link{convertgen}} function.
#' @param hybrid_phe a data frame with three columns. The first column and the second column are the names of male and female parents of the corresponding hybrids, respectively; the third column is the phenotypic values of hybrids.
#' The names of male and female parents must match the rownames of inbred_gen. Missing (NA) values are not allowed.
#' @param method eight GS methods including "GBLUP", "BayesB", "RKHS", "PLS", "LASSO", "EN", "XGBOOST", "RF".
#' Users may select one of these methods. Default is "GBLUP".
#' @param model the prediction model. There are two options: model = "A" for the additive model, model = "AD" for the additive-dominance model. Default is model = "A".
#' @param select the selection of hybrids based on the prediction results. There are three options: select = "all", which selects all potential crosses. select = "top", which selects the top n crosses. select = "bottom", which selects the bottom n crosses. The n is determined by the param number.
#' @param number the number of selected top or bottom hybrids, only when select = "top" or select = "bottom".
#' @return a data frame of prediction results with two columns. The first column denotes the names of male and female parents of the predicted hybrids, and the second column denotes the phenotypic values of the predicted hybrids.
#' @examples
#' \donttest{
#' ## load example data from hypred package
#' data(hybrid_phe)
#' data(input_geno)
#' inbred_gen <- convertgen(input_geno, type = "hmp2")
#'
#' ## infer the additive and dominance genotypes of hybrids
#' gena <- infergen(inbred_gen, hybrid_phe)$add
#' gend <- infergen(inbred_gen, hybrid_phe)$dom
#'
#' pred<-predhy.predict(inbred_gen,hybrid_phe,method="LASSO",model="A",select="top",number="100")
#' pred<-predhy.predict(inbred_gen,hybrid_phe,method="LASSO",model="AD",select="all")
#'  }
#' @export
predhy.predict <- function(inbred_gen, hybrid_phe, method = "GBLUP", model = "A", select = "top", number = "100") {
    gena <- infergen(inbred_gen, hybrid_phe)$add
    gend <- infergen(inbred_gen, hybrid_phe)$dom
    y <- hybrid_phe[, 3]
    if ((method == "PLS") | (method == "XGBOOST") | (method == "BayesB") | (method ==
        "LASSO") | (method == "GBLUP") | (method == "RKHS") | (method == "RF") |
        (method == "EN")) {
        if (method == "PLS") {
            predict.pls <- function(inbred_gen, hybrid_phe, model = NULL) {
                if (!requireNamespace("pls", quietly = TRUE)) {
                  stop("pls needed for this function to work. Please install it.",
                    call. = FALSE)
                }
                # library(pls)
                predparent_gen <- t(inbred_gen)
                pred_name <- colnames(predparent_gen)
                predparent_gen <- as.matrix(predparent_gen)
                phe_name <- NULL
                ynew <- NULL
                if (model == "AD") {
                  x <- cbind(gena, gend)
                  pls.fit <- plsr(y ~ x, ncomp = 5, validation = "CV")
                  nn <- as.numeric(which.min(tt <- RMSEP(pls.fit)$val[1, , ][-1]))
                  for (i in 1:(ncol(predparent_gen) - 1)) {
                    ha1 <- t((predparent_gen[,i] + predparent_gen[,-(1:i)])/2)
                    ha2 <- t((predparent_gen[,i] - predparent_gen[,-(1:i)])/2)
                    X <- cbind(ha1, ha2)
                    yhat <- predict(pls.fit, newdata = X, ncomp = nn)
                    ynew <- c(ynew, yhat)
                    test_name <- paste(pred_name[i], pred_name[-(1:i)], sep = "/")
                    phe_name <- c(phe_name, test_name)
                  }
                  ynew <- as.matrix(ynew)
                  row.names(ynew) <- phe_name
                  return(ynew)
                }
                if (model == "A") {
                  x <- gena
                  pls.fit <- plsr(y ~ x, ncomp = 5, validation = "CV")
                  nn <- as.numeric(which.min(tt <- RMSEP(pls.fit)$val[1, , ][-1]))
                  for (i in 1:(ncol(predparent_gen) - 1)) {
                    ha1 <- t((predparent_gen[,i] + predparent_gen[,-(1:i)])/2)
                    X <- ha1
                    yhat <- predict(pls.fit, newdata = X, ncomp = nn)
                    ynew <- c(ynew, yhat)
                    test_name <- paste(pred_name[i], pred_name[-(1:i)], sep = "/")
                    phe_name <- c(phe_name, test_name)
                  }
                  ynew <- as.matrix(ynew)
                  row.names(ynew) <- phe_name
                  return(ynew)
                }
            }
            print("Predict by PLS...")
            if (model == "A") {
                print("additive model")
                predict_pls <- predict.pls(inbred_gen, hybrid_phe, model = "A")
                Results <- predict_pls
            } else {
                print("additive-dominance model")
                predict_pls <- predict.pls(inbred_gen, hybrid_phe, model = "AD")
                Results <- predict_pls
            }
            print("Predict by PLS...ended.")
        }
        if (method == "XGBOOST") {
            predict.xgboost <- function(inbred_gen, hybrid_phe,model = NULL) {
                if (!requireNamespace("xgboost", quietly = TRUE)) {
                  stop("xgboost needed for this function to work. Please install it.",
                    call. = FALSE)
                }
                # library(xgboost)
                predparent_gen <- t(inbred_gen)
                pred_name <- colnames(predparent_gen)
                predparent_gen <- as.matrix(predparent_gen)
                phe_name <- NULL
                ynew <- NULL
                if (model == "AD") {
                  x <- cbind(gena, gend)
                  xg <- xgboost(data = x, label = y, nrounds = 1000, eta = 0.07)
                  for (i in 1:(ncol(predparent_gen) - 1)) {
                    ha1 <- t((predparent_gen[, i] + predparent_gen[, -(1:i)])/2)
                    ha2 <- t((predparent_gen[, i] - predparent_gen[, -(1:i)])/2)
                    X <- cbind(ha1, ha2)
                    yhat <- predict(xg, X)
                    ynew <- c(ynew, yhat)
                    test_name <- paste(pred_name[i], pred_name[-(1:i)], sep = "/")
                    phe_name <- c(phe_name, test_name)
                  }
                  ynew <- as.matrix(ynew)
                  row.names(ynew) <- phe_name
                  return(ynew)
                }
                if (model == "A") {
                  x <- gena
                  xg <- xgboost(data = x, label = y, nrounds = 1000, eta = 0.07)
                  for (i in 1:(ncol(predparent_gen) - 1)) {
                    ha1 <-  t((predparent_gen[, i] + predparent_gen[, -(1:i)])/2)
                    X <- ha1
                    yhat <- predict(xg, X)
                    ynew <- c(ynew, yhat)
                    test_name <- paste(pred_name[i], pred_name[-(1:i)], sep = "/")
                    phe_name <- c(phe_name, test_name)
                  }
                  ynew <- as.matrix(ynew)
                  row.names(ynew) <- phe_name
                  return(ynew)
                }
            }
            print("Predict by XGBOOST...")
            if (model == "A") {
                print("additive model")
                predict_xgboost <- predict.xgboost(inbred_gen, hybrid_phe,model = "A")
                Results <- predict_xgboost
            } else {
                print("additive-dominance model")
                predict_xgboost <- predict.xgboost(inbred_gen, hybrid_phe,model = "AD")
                Results <- predict_xgboost
            }
            print("Predict by XGBOOST...ended")
        }
        if (method == "BayesB") {
            predict.bayesb <- function(inbred_gen, hybrid_phe,model = NULL) {
                if (!requireNamespace("BGLR", quietly = TRUE)) {
                  stop("BGLR needed for this function to work. Please install it.",
                    call. = FALSE)
                }
                # library(BGLR)
                predparent_gen <- t(inbred_gen)
                pred_name <- colnames(predparent_gen)
                predparent_gen <- as.matrix(predparent_gen)
                phe_name <- NULL
                ynew <- NULL
                if (model == "AD") {
                  for (i in 1:(ncol(predparent_gen) - 1)) {
                    ha1 <-  t((predparent_gen[, i] + predparent_gen[, -(1:i)])/2)
                    X1 <- rbind(gena, ha1)
                    ha2 <-  t((predparent_gen[, i] - predparent_gen[, -(1:i)])/2)
                    X2 <- rbind(gend, ha2)
                    X <- cbind(X1, X2)
                    yNa <- as.matrix(rep(NA, nrow(ha1)))
                    yNA <- c(y, yNa)
                    yNA <- as.matrix(yNA)
                    eta <- list(list(X = X, model = "BayesB"))
                    fm <- BGLR(y = yNA, ETA = eta, verbose = F)
                    yhat <- fm$yHat
                    y <- as.matrix(y)
                    yhat1 <- yhat[-c(1:nrow(y))]
                    ynew <- c(ynew, yhat1)
                    test_name <- paste(pred_name[i], pred_name[-(1:i)], sep = "/")
                    phe_name <- c(phe_name, test_name)
                  }
                  ynew <- as.matrix(ynew)
                  row.names(ynew) <- phe_name
                  return(ynew)
                }
                if (model == "A") {
                  for (i in 1:(ncol(predparent_gen) - 1)) {
                    ha1 <-  t((predparent_gen[, i] + predparent_gen[, -(1:i)])/2)
                    X <- rbind(gena, ha1)
                    yNa <- as.matrix(rep(NA, nrow(ha1)))
                    yNA <- c(y, yNa)
                    yNA <- as.matrix(yNA)
                    eta <- list(list(X = X, model = "BayesB"))
                    fm <- BGLR(y = yNA, ETA = eta, verbose = F)
                    yhat <- fm$yHat
                    y <- as.matrix(y)
                    yhat1 <- yhat[-c(1:nrow(y))]
                    ynew <- c(ynew, yhat1)
                    test_name <- paste(pred_name[i], pred_name[-(1:i)], sep = "/")
                    phe_name <- c(phe_name, test_name)
                  }
                  ynew <- as.matrix(ynew)
                  row.names(ynew) <- phe_name
                  return(ynew)
                }
            }
            print("Predict by BayesB...")
            if (model == "A") {
                print("additive model")
                predict_bayesb <- predict.bayesb(inbred_gen, hybrid_phe, model = "A")
                Results <- predict_bayesb
            } else {
                print("additive-dominance model")
                predict_bayesb <- predict.bayesb(inbred_gen, hybrid_phe,model = "AD")
                Results <- predict_bayesb
            }
            print("Predict by BayesB ...ended")
        }
        if (method == "LASSO") {
            predict.lasso <- function(inbred_gen, hybrid_phe,model = NULL) {
                if (!requireNamespace("glmnet", quietly = TRUE)) {
                  stop("glmnet needed for this function to work. Please install it.",
                    call. = FALSE)
                }
                # library(glmnet)
                predparent_gen <- t(inbred_gen)
                pred_name <- colnames(predparent_gen)
                predparent_gen <- as.matrix(predparent_gen)
                phe_name <- NULL
                ynew <- NULL
                if (model == "AD") {
                  x <- cbind(gena, gend)
                  fit0 <- cv.glmnet(x, y)
                  lambda <- fit0$lambda.min
                  ffit <- glmnet(x, y, lambda = lambda)
                  for (i in 1:(ncol(predparent_gen) - 1)) {
                    ha1 <- t((predparent_gen[, i] + predparent_gen[, -(1:i)])/2)
                    ha2 <- t((predparent_gen[, i] - predparent_gen[, -(1:i)])/2)
                    X <- cbind(ha1, ha2)
                    yhat <- predict(ffit, newx = X)
                    ynew <- c(ynew, yhat)
                    test_name <- paste(pred_name[i], pred_name[-(1:i)], sep = "/")
                    phe_name <- c(phe_name, test_name)
                  }
                  ynew <- as.matrix(ynew)
                  row.names(ynew) <- phe_name
                  return(ynew)
                }
                if (model == "A") {
                  x <- gena
                  fit0 <- cv.glmnet(x, y)
                  lambda <- fit0$lambda.min
                  ffit <- glmnet(x, y, lambda = lambda)
                  for (i in 1:(ncol(predparent_gen) - 1)) {
                    ha1 <- t((predparent_gen[, i] + predparent_gen[, -(1:i)])/2)
                    X <- ha1
                    yhat <- predict(ffit, newx = X)
                    ynew <- c(ynew, yhat)
                    test_name <- paste(pred_name[i], pred_name[-(1:i)], sep = "/")
                    phe_name <- c(phe_name, test_name)
                  }
                  ynew <- as.matrix(ynew)
                  row.names(ynew) <- phe_name
                  return(ynew)
                }
            }
            print("Predict by LASSO ...")
            if (model == "A") {
                print("additive model")
                predict_lasso <- predict.lasso(inbred_gen, hybrid_phe, model = "A")
                Results <- predict_lasso
            } else {
                print("additive-dominance model")
                predict_lasso <- predict.lasso(inbred_gen, hybrid_phe,model = "AD")
                Results <- predict_lasso
            }
            print("Predict by LASSO ...ended")
        }
        if (method == "RF") {
            predict.rf <- function(inbred_gen, hybrid_phe,model = NULL) {
                if (!requireNamespace("randomForest", quietly = TRUE)) {
                  stop("randomForest needed for this function to work. Please install it.",
                    call. = FALSE)
                }
                # library(randomForest)
                predparent_gen <- t(inbred_gen)
                pred_name <- colnames(predparent_gen)
                predparent_gen <- as.matrix(predparent_gen)
                ynew <- NULL
                phe_name <- NULL
                if (model == "AD") {
                  x <- cbind(gena, gend)
                  fit <- randomForest(x = x, y = y, ntree = 500)
                  for (i in 1:(ncol(predparent_gen) - 1)) {
                    ha1 <- t((predparent_gen[, i] + predparent_gen[, -(1:i)])/2)
                    ha2 <- t((predparent_gen[, i] - predparent_gen[, -(1:i)])/2)
                    X <- cbind(ha1, ha2)
                    yhat <- predict(fit, X)
                    ynew <- c(ynew, yhat)
                    test_name <- paste(pred_name[i], pred_name[-(1:i)], sep = "/")
                    phe_name <- c(phe_name, test_name)
                  }
                  ynew <- as.matrix(ynew)
                  row.names(ynew) <- phe_name
                  return(ynew)
                }
                if (model == "A") {
                  x <- gena
                  fit <- randomForest(x = x, y = y, ntree = 500)
                  for (i in 1:(ncol(predparent_gen) - 1)) {
                    ha1 <- t((predparent_gen[, i] + predparent_gen[, -(1:i)])/2)
                    X <- ha1
                    yhat <- predict(fit, X)
                    ynew <- c(ynew, yhat)
                    test_name <- paste(pred_name[i], pred_name[-(1:i)], sep = "/")
                    phe_name <- c(phe_name, test_name)
                  }
                  ynew <- as.matrix(ynew)
                  row.names(ynew) <- phe_name
                  return(ynew)
                }
            }
            print("Predict by RF ...")
            if (model == "A") {
                print("additive model")
                predict_rf <- predict.rf(inbred_gen, hybrid_phe,model = "A")
                Results <- predict_rf
            } else {
                print("additive-dominance model")
                predict_rf <- predict.rf(inbred_gen, hybrid_phe,model = "AD")
                Results <- predict_rf
            }
            print("Predict by RF ...ended")
        }
        if (method == "EN") {
            predict.EN <- function(inbred_gen, hybrid_phe,model = NULL) {
                # library(glmnet)
                predparent_gen <- t(inbred_gen)
                pred_name <- colnames(predparent_gen)
                predparent_gen <- as.matrix(predparent_gen)
                ynew <- NULL
                phe_name <- NULL
                if (model == "AD") {
                  x <- cbind(gena, gend)
                  fit0 <- cv.glmnet(x = x, y = y, alpha = 0.5)
                  lambda <- fit0$lambda.min
                  ffit <- glmnet(x = x, y = y, lambda = lambda, alpha = 0.5)
                  for (i in 1:(ncol(predparent_gen) - 1)) {
                    ha1 <-  t((predparent_gen[, i] + predparent_gen[, -(1:i)])/2)
                    ha2 <-  t((predparent_gen[, i] - predparent_gen[, -(1:i)])/2)
                    X <- cbind(ha1, ha2)
                    yhat <- predict(ffit, newx = X)
                    ynew <- c(ynew, yhat)
                    test_name <- paste(pred_name[i], pred_name[-(1:i)], sep = "/")
                    phe_name <- c(phe_name, test_name)
                  }
                  ynew <- as.matrix(ynew)
                  row.names(ynew) <- phe_name
                  return(ynew)
                }
                if (model == "A") {
                  x <- gena
                  fit0 <- cv.glmnet(x = x, y = y, alpha = 0.5)
                  lambda <- fit0$lambda.min
                  ffit <- glmnet(x = x, y = y, lambda = lambda, alpha = 0.5)
                  for (i in 1:(ncol(predparent_gen) - 1)) {
                    ha1 <-  t((predparent_gen[, i] + predparent_gen[, -(1:i)])/2)
                    X <- ha1
                    yhat <- predict(ffit, newx = X)
                    ynew <- c(ynew, yhat)
                    test_name <- paste(pred_name[i], pred_name[-(1:i)], sep = "/")
                    phe_name <- c(phe_name, test_name)
                  }
                  ynew <- as.matrix(ynew)
                  row.names(ynew) <- phe_name
                  return(ynew)
                }
            }
            print("Predict by EN ...")
            if (model == "A") {
                print("additive model")
                predict_EN <- predict.EN(inbred_gen, hybrid_phe,model = "A")
                Results <- predict_EN
            } else {
                print("additive-dominance model")
                predict_EN <- predict.EN(inbred_gen, hybrid_phe,model = "AD")
                Results <- predict_EN
            }
            print("Predict by EN ...ended")
        }
        if (method == "GBLUP") {
            predict.GBLUP <- function(fix = NULL, fixnew = NULL, inbred_gen, hybrid_phe, model = NULL) {
                n <- nrow(hybrid_phe)
                if (model == "AD") {
                  gena <- infergen(inbred_gen, hybrid_phe)[[1]]
                  ka <- kin(gena)
                  gend <- infergen(inbred_gen, hybrid_phe)[[2]]
                  kd <- kin(gend)
                  parm <- mixed(fix = fix, y = y, kk = list(ka, kd))
                  v_i <- parm$v_i
                  beta <- parm$beta
                  va <- parm$var[1]
                  vd <- parm$var[2]
                  ve <- parm$ve
                  ka21 <- NULL
                  kd21 <- NULL
                  phe_name <- NULL
                  predparent_gen <- t(inbred_gen)
                  pred_name <- colnames(predparent_gen)
                  predparent_gen <- as.matrix(predparent_gen)
                  for (i in 1:(ncol(predparent_gen) - 1)) {
                    ha1 <-  t((predparent_gen[, i] + predparent_gen[, -(1:i)])/2)
                    ka2 <- tcrossprod(ha1, gena)/ncol(gena)
                    ka21 <- rbind(ka21, ka2)
                    hd1 <- abs(t((predparent_gen[, i] - predparent_gen[, -(1:i)])/2))
                    kd2 <- tcrossprod(hd1, gend)/ncol(gend)
                    kd21 <- rbind(kd21, kd2)
                    test_name <- paste(pred_name[i], pred_name[-(1:i)], sep = "/")
                    phe_name <- c(phe_name, test_name)
                  }
                  if (is.null(fix)) {
                    fix <- matrix(1, n, 1)
                  } else {
                    fix <- as.matrix(fix)
                  }
                  n1 <- nrow(kd21)
                  if (is.null(fixnew)) {
                    fixnew <- matrix(1, n1, 1)
                  } else {
                    fixnew <- as.matrix(fixnew)
                  }
                  G21 <- ka21 * va + kd21 * vd
                  pred_phe <- fixnew %*% beta + G21 %*% v_i %*% (y - fix %*% beta)
                  row.names(pred_phe) <- phe_name
                  return(pred_phe)
                }
                if (model == "A") {
                  gena <- infergen(inbred_gen, hybrid_phe)[[1]]
                  ka <- kin(gena)
                  parm <- mixed(fix = fix, y = y, kk = list(ka))
                  v_i <- parm$v_i
                  beta <- parm$beta
                  va <- parm$var[1]
                  ve <- parm$ve
                  ka21 <- NULL
                  phe_name <- NULL
                  predparent_gen <- t(inbred_gen)
                  pred_name <- colnames(predparent_gen)
                  for (i in 1:(ncol(predparent_gen) - 1)) {
                    ha1 <-  t((predparent_gen[, i] + predparent_gen[, -(1:i)])/2)
                    ka2 <- tcrossprod(ha1, gena)/ncol(gena)
                    ka21 <- rbind(ka21, ka2)
                    test_name <- paste(pred_name[i], pred_name[-(1:i)], sep = "/")
                    phe_name <- c(phe_name, test_name)
                  }
                  if (is.null(fix)) {
                    fix <- matrix(1, n, 1)
                  } else {
                    fix <- as.matrix(fix)
                  }
                  n1 <- nrow(ka21)
                  if (is.null(fixnew)) {
                    fixnew <- matrix(1, n1, 1)
                  } else {
                    fixnew <- as.matrix(fixnew)
                  }
                  G21 <- ka21 * va
                  pred_phe <- fixnew %*% beta + G21 %*% v_i %*% (y - fix %*% beta)
                  row.names(pred_phe) <- phe_name
                  return(pred_phe)
                }
            }
            print("Predict by GBLUP ...")
            if (model == "A") {
                print("additive model")
                predict_GBLUP <- predict.GBLUP(fix = NULL, fixnew = NULL, inbred_gen,
                  hybrid_phe, model = "A")
                Results <- predict_GBLUP
            } else {
                print("additive-dominance model")
                predict_GBLUP <- predict.GBLUP(fix = NULL, fixnew = NULL, inbred_gen,
                  hybrid_phe, model = "AD")
                Results <- predict_GBLUP
            }
            print("Predict by GBLUP ...ended")
        }
        if (method == "RKHS") {
            predict.rkhsmk <- function(inbred_gen,hybrid_phe,model = NULL) {
                # library(BGLR)
                predparent_gen <- t(inbred_gen)
                pred_name <- colnames(predparent_gen)
                predparent_gen <- as.matrix(predparent_gen)
                ynew <- NULL
                phe_name <- NULL
                if (model == "AD") {
                  for (i in 1:(ncol(predparent_gen) - 1)) {
                    ha1 <-  t((predparent_gen[, i] + predparent_gen[, -(1:i)])/2)
                    X1 <- rbind(gena, ha1)
                    ha2 <-  t((predparent_gen[, i] - predparent_gen[, -(1:i)])/2)
                    X2 <- rbind(gend, ha2)
                    X <- cbind(X1, X2)
                    yNa <- as.matrix(rep(NA, nrow(ha1)))
                    yNA <- c(y, yNa)
                    yNA <- as.matrix(yNA)
                    M <- scale(X, center = TRUE, scale = TRUE)
                    D <- (as.matrix(dist(M, method = "euclidean"))^2)/ncol(X)
                    h <- 0.5 * c(1/5, 1, 5)
                    ETA <- list(list(K = exp(-h[1] * D), model = "RKHS"), list(K = exp(-h[2] *
                      D), model = "RKHS"), list(K = exp(-h[3] * D), model = "RKHS"))
                    fm <- BGLR(y = yNA, ETA = ETA, verbose = F)
                    yhat <- fm$yHat
                    y <- as.matrix(y)
                    yhat1 <- yhat[-c(1:nrow(y))]
                    ynew <- c(ynew, yhat1)
                    test_name <- paste(pred_name[i], pred_name[-(1:i)], sep = "/")
                    phe_name <- c(phe_name, test_name)
                  }
                  ynew <- as.matrix(ynew)
                  row.names(ynew) <- phe_name
                  return(ynew)
                }
                if (model == "A") {
                  for (i in 1:(ncol(predparent_gen) - 1)) {
                    ha1 <-  t((predparent_gen[, i] + predparent_gen[, -(1:i)])/2)
                    X1 <- rbind(gena, ha1)
                    X <- X1
                    yNa <- as.matrix(rep(NA, nrow(ha1)))
                    yNA <- c(y, yNa)
                    yNA <- as.matrix(yNA)
                    M <- scale(X, center = TRUE, scale = TRUE)
                    D <- (as.matrix(dist(M, method = "euclidean"))^2)/ncol(X)
                    h <- 0.5 * c(1/5, 1, 5)
                    ETA <- list(list(K = exp(-h[1] * D), model = "RKHS"), list(K = exp(-h[2] *
                      D), model = "RKHS"), list(K = exp(-h[3] * D), model = "RKHS"))
                    fm <- BGLR(y = yNA, ETA = ETA, verbose = F)
                    yhat <- fm$yHat
                    y <- as.matrix(y)
                    yhat1 <- yhat[-c(1:nrow(y))]
                    ynew <- c(ynew, yhat1)
                    test_name <- paste(pred_name[i], pred_name[-(1:i)], sep = "/")
                    phe_name <- c(phe_name, test_name)
                  }
                  ynew <- as.matrix(ynew)
                  row.names(ynew) <- phe_name
                  return(ynew)
                }
            }
            print("Predict by RKHS ...")
            if (model == "A") {
                print("additive model")
                predict_rkhsmk <- predict.rkhsmk(inbred_gen, hybrid_phe,model = "A")
                Results <- predict_rkhsmk
            } else {
                print("additive-dominance model")
                predict_rkhsmk <- predict.rkhsmk(inbred_gen, hybrid_phe, model = "AD")
                Results <- predict_rkhsmk
            }
            print("Predict by RKHS ...ended")
        }
        if (select == "all") {
            Results_select <- as.data.frame(Results)
            colnames(Results_select) <- paste("all_", nrow(Results), sep = "")
        } else if (select == "top") {
            Results_select <- as.data.frame(sort(Results[, 1], decreasing = T)[c(1:number)])
            names(Results_select) <- paste("top_", number, sep = "")
        } else if (select == "bottom") {
            Results_select <- as.data.frame(sort(Results[, 1], decreasing = F)[c(1:number)])
            colnames(Results_select) <- paste("bottom_", number, sep = "")
        }
        return(Results_select)
    } else {
        stop("Please choose a predict method: GBLUP, BayesB, RKHS, PLS, LASSO, EN, XGBOOST, RF.")
    }  # End of all predict methods
    return(Results)
}

