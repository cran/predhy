#' @title Predict the Performance of Hybrids
#' @description Predict all potential crosses of a given set of parents using a subset of crosses as the training sample.
#' @param inbred_gen a matrix for genotypes of parental lines in numeric format, coded as 1, 0 and -1. The row.names of inbred_gen must be provied. It can be obtained from the original genotype using  \code{\link{convertgen}} function.
#' @param hybrid_phe a data frame with three columns. The first column and the second column are the names of male and female parents of the corresponding hybrids, respectively; the third column is the phenotypic values of hybrids.
#' The names of male and female parents must match the rownames of inbred_gen. Missing (NA) values are not allowed.
#' @param inbred_phe a matrix (n x 2) of inbred_phe phenotypic.Default is NULL.
#' @param male_name a vector of the names of male parents.
#' @param female_name a vector of the names of female parents.
#' @param method eight GS methods including "GBLUP", "BayesB", "RKHS", "PLS", "LASSO", "EN", "XGBoost", "LightGBM".
#' Users may select P of these methods. Default is "GBLUP".
#' @param model the prediction model. There are two options: model = "A" for the additive model, model = "AD" for the additive-dominance model,model = "A-P" for the additive-phenotypic model,model = "AD-P" for the additive-dominance-phenotypic model. Default is model = "A".
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
#' pred<-predhy.predict_NCII(inbred_gen,hybrid_phe,method="LASSO",model="A")
#' pred<-predhy.predict_NCII(inbred_gen,hybrid_phe,method="LASSO",model = "AD",select="all")
#'  }
#' @export
predhy.predict_NCII <- function (inbred_gen, hybrid_phe, inbred_phe=NULL,male_name=hybrid_phe[,1],female_name=hybrid_phe[,2],
                             method = "GBLUP", model = "A",select = "top", number = "100"){
  gena <- infergen(inbred_gen, hybrid_phe)$add
  gend <- infergen(inbred_gen, hybrid_phe)$dom
  y <- hybrid_phe[, 3]
  m <- as.character(male_name)
  f <- as.character(female_name)
  inbred_m <- inbred_gen[row.names(inbred_gen) %in% m ,]
  inbred_f <- inbred_gen[row.names(inbred_gen) %in% f ,]
  predparent_gen_m <- as.matrix(t(inbred_m))
  predparent_gen_f <- as.matrix(t(inbred_f))
  if (is.null(inbred_phe)) {
    print(message("no phenotypic values of inbreds"))
  }else {
  inbred_phe <- as.matrix(inbred_phe)
  pena <- infergen(inbred_phe, hybrid_phe)$add
  pend <- infergen(inbred_phe, hybrid_phe)$dom
  inbredphe<-cbind(pena,pend)
  inbred_phe_m<-inbred_phe[row.names(inbred_phe) %in% m ,]
  inbred_phe_f<-inbred_phe[row.names(inbred_phe) %in% f ,]
  predparent_phe_m <- as.matrix(t(inbred_phe_m))
  predparent_phe_f <- as.matrix(t(inbred_phe_f))}
  pred_name_m <- colnames(predparent_gen_m)
  pred_name_f <- colnames(predparent_gen_f)
  if ((method == "PLS") | (method == "XGBoost") | (method == "BayesB") | (method ==
        "LASSO") | (method == "GBLUP") | (method == "RKHS") | (method == "LightGBM") |
        (method == "EN")) {
    if(method=='PLS'){
      predict.pls <- function(model = NULL) {
        if (!requireNamespace("pls", quietly = TRUE)) {
          stop("pls needed for this function to work. Please install it.", 
               call. = FALSE)
        }
        phe_name <- NULL
        ynew <- NULL
        if (model == "AD") {
          x <- cbind(gena, gend)
          pls.fit <- plsr(y ~ x, ncomp = 5, validation = "CV")
          nn <- as.numeric(which.min(tt <- RMSEP(pls.fit)$val[1, , ][-1]))
          for (i in 1:(ncol(predparent_gen_m))) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
            ha2 <- t(abs((predparent_gen_m[, i] - predparent_gen_f)/2))
            X <- cbind(ha1, ha2)
            yhat <- predict(pls.fit, newdata = X, ncomp = nn)
            ynew <- c(ynew, yhat)
            test_name <- paste(pred_name_m[i], pred_name_f,sep = "/")
            phe_name <- c(phe_name, test_name)
          }
          ynew <- as.matrix(ynew)
          row.names(ynew) <- phe_name
          return(ynew)
        }
		if (model == "AD-P") {
          x <- cbind(gena, gend,inbredphe)
          pls.fit <- plsr(y ~ x, ncomp = 5, validation = "CV")
          nn <- as.numeric(which.min(tt <- RMSEP(pls.fit)$val[1, , ][-1]))
          for (i in 1:(ncol(predparent_gen_m))) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
            ha2 <- t(abs((predparent_gen_m[, i] - predparent_gen_f)/2))
			p1<-t((predparent_phe_m[,i] + predparent_phe_f)/2)
			p2<-t(abs((predparent_phe_m[,i] - predparent_phe_f)/2))
            X <- cbind(ha1,ha2,p1,p2)
            yhat <- predict(pls.fit, newdata = X, ncomp = nn)
            ynew <- c(ynew, yhat)
            test_name <- paste(pred_name_m[i], pred_name_f,sep = "/")
            phe_name <- c(phe_name, test_name)
          }
          ynew <- as.matrix(ynew)
          row.names(ynew) <- phe_name
          return(ynew)
        }
        if (model == "A") {
          x <- gena
          pls.fit <- plsr(y ~ x, ncomp = 5, validation = "CV")
          nn <- as.numeric(which.min(tt <- RMSEP(pls.fit)$val[1,, ][-1]))
          for (i in 1:(ncol(predparent_gen_m))) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
            X <- ha1
            yhat <- predict(pls.fit, newdata = X, ncomp = nn)
            ynew <- c(ynew, yhat)
            test_name <- paste(pred_name_m[i], pred_name_f,sep = "/")
            phe_name <- c(phe_name, test_name)
          }
          ynew <- as.matrix(ynew)
          row.names(ynew) <- phe_name
          return(ynew)
        }
		if (model == "A-P") {
          x <- cbind(gena,inbredphe)
          pls.fit <- plsr(y ~ x, ncomp = 5, validation = "CV")
          nn <- as.numeric(which.min(tt <- RMSEP(pls.fit)$val[1,, ][-1]))
          for (i in 1:(ncol(predparent_gen_m))) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
			p1<-t((predparent_phe_m[,i] + predparent_phe_f)/2)
			p2<-t(abs((predparent_phe_m[,i] - predparent_phe_f)/2))
            X <- cbind(ha1,p1,p2)
            yhat <- predict(pls.fit, newdata = X, ncomp = nn)
            ynew <- c(ynew, yhat)
            test_name <- paste(pred_name_m[i], pred_name_f,sep = "/")
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
        predict_pls <- predict.pls(model = "A")
        Results <- predict_pls
      }
	  else if (model == "A-P") {
        print("additive-phenotypic model")
        predict_pls <- predict.pls(model = "A-P")
        Results <- predict_pls
      }
	  else if (model == "AD-P") {
        print("additive-dominance-phenotypic model")
        predict_pls <- predict.pls(model = "AD-P")
        Results <- predict_pls
      }
	  else {
        print("additive-dominance model")
        predict_pls <- predict.pls(model = "AD")
        Results <- predict_pls
      }
      print("Predict by PLS...ended.")}
    if (method == "BayesB") {
      predict.bayesb <- function(model = NULL) {
        if (!requireNamespace("BGLR", quietly = TRUE)) {
          stop("BGLR needed for this function to work. Please install it.", 
               call. = FALSE)
        }
        phe_name <- NULL
        ynew <- NULL
        if (model == "AD") {
          for (i in 1:(ncol(predparent_gen_m))) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
            X1 <- rbind(gena, ha1)
            ha2 <- t(abs((predparent_gen_m[, i] - predparent_gen_f)/2))
            X2 <- rbind(gend, ha2)
            yNa <- as.matrix(rep(NA, nrow(ha1)))
            yNA <- c(y, yNa)
            yNA <- as.matrix(yNA)
            eta <- list(list(X = X1, model = "BayesB"),list(X = X2, model = "BayesB"))
            fm <- BGLR(y = yNA, ETA = eta, verbose = F)
            yhat <- fm$yHat
            y <- as.matrix(y)
            yhat1 <- yhat[-c(1:nrow(y))]
            ynew <- c(ynew, yhat1)
            test_name <- paste(pred_name_m[i], pred_name_f,sep = "/")
            phe_name <- c(phe_name, test_name)
          }
          ynew <- as.matrix(ynew)
          row.names(ynew) <- phe_name
          return(ynew)
        }
		if (model == "AD-P") {
		for (i in 1:(ncol(predparent_gen_m))) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
            X1 <- rbind(gena, ha1)
			ha2 <- t(abs((predparent_gen_m[, i] - predparent_gen_f)/2))
            X2 <- rbind(gend, ha2)
			p1<-t((predparent_phe_m[,i] + predparent_phe_f)/2)
			p2<-t(abs((predparent_phe_m[,i] - predparent_phe_f)/2))
			P<-cbind(p1,p2)
			P1<-rbind(inbredphe,P)
            yNa <- as.matrix(rep(NA, nrow(ha1)))
            yNA <- c(y, yNa)
            yNA <- as.matrix(yNA)
            eta <- list(list(X = X1, model = "BayesB"),list(X = X2, model = "BayesB"),list(X = P1, model = "BayesB"))
            fm <- BGLR(y = yNA, ETA = eta, verbose = F)
            yhat <- fm$yHat
            y <- as.matrix(y)
            yhat1 <- yhat[-c(1:nrow(y))]
            ynew <- c(ynew, yhat1)
            test_name <- paste(pred_name_m[i], pred_name_f,sep = "/")
            phe_name <- c(phe_name, test_name)
          }
          ynew <- as.matrix(ynew)
          row.names(ynew) <- phe_name
          return(ynew)
        }
		
		if (model == "A-P") {
		for (i in 1:(ncol(predparent_gen_m))) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
            X <- rbind(gena, ha1)
			p1<-t((predparent_phe_m[,i] + predparent_phe_f)/2)
			p2<-t(abs((predparent_phe_m[,i] - predparent_phe_f)/2))
			P<-cbind(p1,p2)
			P1<-rbind(inbredphe,P)
            yNa <- as.matrix(rep(NA, nrow(ha1)))
            yNA <- c(y, yNa)
            yNA <- as.matrix(yNA)
            eta <- list(list(X = X, model = "BayesB"),list(X = P1, model = "BayesB"))
            fm <- BGLR(y = yNA, ETA = eta, verbose = F)
            yhat <- fm$yHat
            y <- as.matrix(y)
            yhat1 <- yhat[-c(1:nrow(y))]
            ynew <- c(ynew, yhat1)
            test_name <- paste(pred_name_m[i], pred_name_f,sep = "/")
            phe_name <- c(phe_name, test_name)
          }
          ynew <- as.matrix(ynew)
          row.names(ynew) <- phe_name
          return(ynew)
        }
        if (model == "A") {
          for (i in 1:(ncol(predparent_gen_m))) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
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
            test_name <- paste(pred_name_m[i], pred_name_f,sep = "/")
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
        predict_bayesb <- predict.bayesb(model = "A")
        Results <- predict_bayesb
      }
	  else if (model == "A-P") {
        print("additive-phenotypic model")
        predict_bayesb <- predict.bayesb(model = "A-P")
        Results <- predict_bayesb
      }
	  else if (model == "AD-P") {
        print("additive-dominance-phenotypic model")
        predict_bayesb <- predict.bayesb(model = "AD-P")
        Results <- predict_bayesb
      }
      else {
        print("additive-dominance model")
        predict_bayesb <- predict.bayesb(model = "AD")
        Results <- predict_bayesb
      }
      print("Predict by BayesB ...ended")
    }
    if (method == "LASSO") {
      predict.lasso <- function(model = NULL) {
        if (!requireNamespace("glmnet", quietly = TRUE)) {
          stop("glmnet needed for this function to work. Please install it.", 
               call. = FALSE)
        }
        phe_name <- NULL
        ynew <- NULL
        if (model == "AD") {
          x <- cbind(gena, gend)
          fit0 <- cv.glmnet(x, y)
          lambda <- fit0$lambda.min
          ffit <- glmnet(x, y, lambda = lambda)
          for (i in 1:ncol(predparent_gen_m)) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
            ha2 <- t(abs((predparent_gen_m[, i] - predparent_gen_f)/2))
            X <- cbind(ha1, ha2)
            yhat <- predict(ffit, newx = X)
            ynew <- c(ynew, yhat)
            test_name <- paste(pred_name_m[i], pred_name_f,sep = "/")
            phe_name <- c(phe_name, test_name)
          }
          ynew <- as.matrix(ynew)
          row.names(ynew) <- phe_name
          return(ynew)
        }
		if (model == "AD-P") {
          x <- cbind(gena, gend,inbredphe)
          fit0 <- cv.glmnet(x, y)
          lambda <- fit0$lambda.min
          ffit <- glmnet(x, y, lambda = lambda)
          for (i in 1:ncol(predparent_gen_m)) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
            ha2 <- t(abs((predparent_gen_m[, i] - predparent_gen_f)/2))
			p1<-t((predparent_phe_m[,i] + predparent_phe_f)/2)
			p2<-t(abs((predparent_phe_m[,i] - predparent_phe_f)/2))
            X <- cbind(ha1, ha2,p1,p2)
            yhat <- predict(ffit, newx = X)
            ynew <- c(ynew, yhat)
            test_name <- paste(pred_name_m[i], pred_name_f,sep = "/")
            phe_name <- c(phe_name, test_name)
          }
          ynew <- as.matrix(ynew)
          row.names(ynew) <- phe_name
          return(ynew)
        }
        if (model == "A-P") {
          x <- cbind(gena,inbredphe)
          fit0 <- cv.glmnet(x, y)
          lambda <- fit0$lambda.min
          ffit <- glmnet(x, y, lambda = lambda)
          for (i in 1:ncol(predparent_gen_m)) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
			p1<-t((predparent_phe_m[,i] + predparent_phe_f)/2)
			p2<-t(abs((predparent_phe_m[,i] - predparent_phe_f)/2))
            X <- cbind(ha1,p1,p2)
            yhat <- predict(ffit, newx = X)
            ynew <- c(ynew, yhat)
            test_name <- paste(pred_name_m[i], pred_name_f, sep = "/")
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
          for (i in 1:ncol(predparent_gen_m)) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
            X <- ha1
            yhat <- predict(ffit, newx = X)
            ynew <- c(ynew, yhat)
            test_name <- paste(pred_name_m[i], pred_name_f, sep = "/")
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
        predict_lasso <- predict.lasso(model = "A")
        Results <- predict_lasso
      }
	  else if (model == "A-P") {
        print("additive-phenotypic model")
	    predict_lasso <- predict.lasso(model = "A-P")
        Results <- predict_lasso
      }
	  else if (model == "AD-P") {
        print("additive-dominance-phenotypic model")
	    predict_lasso <- predict.lasso(model = "AD-P")
        Results <- predict_lasso
      }
      else {
        print("additive-dominance model")
        predict_lasso <- predict.lasso(model = "AD")
        Results <- predict_lasso
      }
      print("Predict by LASSO ...ended")
    }
	
    if(method=='GBLUP'){
      predict.GBLUP <- function(fix = NULL, fixnew = NULL,model=NULL) {
        n <- nrow(hybrid_phe)
        if (model == "AD") {
          ka <- kin(gena)
          kd <- kin(gend)
          parm <- mixed(fix = fix, y = y, kk = list(ka,kd))
          v_i <- parm$v_i
          beta <- parm$beta
          va <- parm$var[1]
          vd <- parm$var[2]
          ve <- parm$ve
          ka21 <- NULL
          kd21 <- NULL
          phe_name <- NULL
          for (i in 1:(ncol(predparent_gen_m))) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
            ka2 <- tcrossprod(ha1, gena)/ncol(gena)
            ka21 <- rbind(ka21, ka2)
            hd1 <- t(abs((predparent_gen_m[, i] - predparent_gen_f)/2))
            kd2 <- tcrossprod(hd1, gend)/ncol(gend)
            kd21 <- rbind(kd21, kd2)
            test_name <- paste(pred_name_m[i], pred_name_f,sep = "/")
            phe_name <- c(phe_name, test_name)
          }
          if (is.null(fix)) {
            fix <- matrix(1, n, 1)
          }else {
            fix <- as.matrix(fix)
          }
          n1 <- nrow(kd21)
          if (is.null(fixnew)) {
            fixnew <- matrix(1, n1, 1)
          }else {
            fixnew <- as.matrix(fixnew)
          }
          
          G21 <- ka21 * va + kd21 * vd
          pred_phe <- fixnew %*% beta + G21 %*% v_i %*% 
            (y - fix %*% beta)
          
          row.names(pred_phe) <- phe_name
          return(pred_phe)
        }
        if (model == "AD-P") {
          ka <- kin(gena)
          kd <- kin(gend)
		  kp<-kin(inbredphe)
          parm <- mixed(fix = fix, y = y, kk = list(ka,kd,kp))
          v_i <- parm$v_i
          beta <- parm$beta
          va <- parm$var[1]
          vd <- parm$var[2]
		  vp<-parm$var[3]
          ve <- parm$ve
          ka21 <- NULL
          kd21 <- NULL
		  kp2<-NULL
          phe_name <- NULL
          for (i in 1:(ncol(predparent_gen_m))) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
            ka2 <- tcrossprod(ha1, gena)/ncol(gena)
            ka21 <- rbind(ka21, ka2)
            hd1 <- t(abs((predparent_gen_m[, i] - predparent_gen_f)/2))
            kd2 <- tcrossprod(hd1, gend)/ncol(gend)
            kd21 <- rbind(kd21, kd2)
			p1<-t((predparent_phe_m[,i] + predparent_phe_f)/2)
			p2<-t(abs((predparent_phe_m[,i] - predparent_phe_f)/2))
			P1<-cbind(p1,p2)
			kp1 <- tcrossprod(P1, inbredphe)/ncol(inbredphe)
            kp2 <- rbind(kp2, kp1)
            test_name <- paste(pred_name_m[i], pred_name_f,sep = "/")
            phe_name <- c(phe_name, test_name)
          }
          if (is.null(fix)) {
            fix <- matrix(1, n, 1)
          }else {
            fix <- as.matrix(fix)
          }
          n1 <- nrow(kd21)
          if (is.null(fixnew)) {
            fixnew <- matrix(1, n1, 1)
          }else {
            fixnew <- as.matrix(fixnew)
          }
          
          G21 <- ka21 * va + kd21 * vd+ kp2 * vp
          pred_phe <- fixnew %*% beta + G21 %*% v_i %*% 
            (y - fix %*% beta)
          
          row.names(pred_phe) <- phe_name
          return(pred_phe)
        }
		if (model == "A-P") {
          ka <- kin(gena)
		  kp<-kin(inbredphe)
          parm <- mixed(fix = fix, y = y, kk = list(ka,kp))
          v_i <- parm$v_i
          beta <- parm$beta
          va <- parm$var[1]
		  vp<-parm$var[2]
          ve <- parm$ve
          ka21 <- NULL
          kp2 <- NULL
          phe_name <- NULL
          for (i in 1:(ncol(predparent_gen_m))) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
            ka2 <- tcrossprod(ha1, gena)/ncol(gena)
            ka21 <- rbind(ka21, ka2)
			p1<-t((predparent_phe_m[,i] + predparent_phe_f)/2)
			p2<-t(abs((predparent_phe_m[,i] - predparent_phe_f)/2))
			P1<-cbind(p1,p2)
			kp1 <- tcrossprod(P1, inbredphe)/ncol(inbredphe)
            kp2 <- rbind(kp2, kp1)
            test_name <- paste(pred_name_m[i], pred_name_f,sep = "/")
            phe_name <- c(phe_name, test_name)
          }
          if (is.null(fix)) {
            fix <- matrix(1, n, 1)
          }else {
            fix <- as.matrix(fix)
          }
          n1 <- nrow(kd21)
          if (is.null(fixnew)) {
            fixnew <- matrix(1, n1, 1)
          }else {
            fixnew <- as.matrix(fixnew)
          }
          
          G21 <- ka21 * va + kp2 * vp
          pred_phe <- fixnew %*% beta + G21 %*% v_i %*% 
            (y - fix %*% beta)
          
          row.names(pred_phe) <- phe_name
          return(pred_phe)
        }
        if (model == "A") {
          ka <- kin(gena)
          parm <- mixed( y = y, kk = list(ka))
          v_i <- parm$v_i
          beta <- parm$beta
          va <- parm$var[1]
          ve <- parm$ve
          ka21 <- NULL
          phe_name <- NULL
          for(i in 1:ncol(predparent_gen_m)){
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
            ka2 <- tcrossprod(ha1, gena)/ncol(gena)
            ka21 <- rbind(ka21, ka2)
            test_name <- paste(pred_name_m[i], pred_name_f,sep = "/")
            phe_name <- c(phe_name, test_name)}
          
          if (is.null(fix)) {
            fix <- matrix(1, n, 1)
          }else {
            fix <- as.matrix(fix)}
          
          n1 <- nrow(ka21)
          if (is.null(fixnew)) {
            fixnew <- matrix(1, n1, 1)
          }else{
            fixnew <- as.matrix(fixnew)}
          
          G21 <- ka21 * va
          pred_phe <- fixnew %*% beta + G21 %*% v_i %*% 
            (y - fix %*% beta)
          row.names(pred_phe) <- phe_name
          return(pred_phe)
        }
      }
      print("Predict by GBLUP ...")
      if (model == "A") {
        print("additive model")
        predict_GBLUP <- predict.GBLUP(fix = NULL, fixnew = NULL, model = "A")
        Results <- predict_GBLUP
      }
	  else if(model == "A-P") {
        print("additive-phenotypic model")
	  predict_GBLUP <- predict.GBLUP(fix = NULL, fixnew = NULL, model = "A-P")
        Results <- predict_GBLUP
      }
	  else if(model == "AD-P") {
        print("additive-dominance-phenotypic model")
	  predict_GBLUP <- predict.GBLUP(fix = NULL, fixnew = NULL, model = "AD-P")
        Results <- predict_GBLUP
      }
	  else {
        print("additive-dominance model")
        predict_GBLUP <- predict.GBLUP(fix = NULL, fixnew = NULL, model = "AD")
        Results <- predict_GBLUP
      }
      print("Predict by GBLUP ...ended")
      
    }
	if (method == "LightGBM") {
            predict.lightgbm <- function(model = NULL) {
                if (!requireNamespace("lightgbm", quietly = TRUE)) {
                  stop("lightgbm needed for this function to work. Please install it.",
                    call. = FALSE)
                }
                # library(lightgbm)
        ynew <- NULL
        phe_name <- NULL
        if (model == "AD") {
          x <- cbind(gena, gend)
		  colnames(x)<-NULL
				  dtrain <- lgb.Dataset(data=x,label=y)
                     params = list(boosting_type = 'gbdt',objective="regression",
                                   learning_rate=0.03 )
                     fit <- lgb.train(
                       params = params,   
                       data=dtrain,   
                       nrounds = 500 ,
                       verbose = -1) 
          for (i in 1:ncol(predparent_gen_m)) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
            ha2 <- t(abs((predparent_gen_m[, i] - predparent_gen_f)/2))
            X <- cbind(ha1, ha2)
			colnames(X)<-NULL
            yhat <- predict(fit, X)
            ynew <- c(ynew, yhat)
            test_name <- paste(pred_name_m[i], pred_name_f,sep = "/")
            phe_name <- c(phe_name, test_name)
          }
          ynew <- as.matrix(ynew)
          row.names(ynew) <- phe_name
          return(ynew)
        }
		 if (model == "AD-P") {
          x <- cbind(gena,gend,inbredphe)
		  colnames(x)<-NULL
				  dtrain <- lgb.Dataset(data=x,label=y)
                     params = list(boosting_type = 'gbdt',objective="regression",
                                   learning_rate=0.03 )
                     fit <- lgb.train(
                       params = params,   
                       data=dtrain,   
                       nrounds = 500 ,
                       verbose = -1) 
          for (i in 1:ncol(predparent_gen_m)) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
            ha2 <- t(abs((predparent_gen_m[, i] - predparent_gen_f)/2))
			p1<-t((predparent_phe_m[,i] + predparent_phe_f)/2)
			p2<-t(abs((predparent_phe_m[,i] - predparent_phe_f)/2))
            X <- cbind(ha1, ha2,p1,p2)
			colnames(X)<-NULL
            yhat <- predict(fit, X)
            ynew <- c(ynew, yhat)
            test_name <- paste(pred_name_m[i], pred_name_f,sep = "/")
            phe_name <- c(phe_name, test_name)
          }
          ynew <- as.matrix(ynew)
          row.names(ynew) <- phe_name
          return(ynew)
        }
        if (model == "A") {
          x <- gena
		  colnames(x)<-NULL
		  dtrain <- lgb.Dataset(data=x,label=y)
                     params = list(boosting_type = 'gbdt',objective="regression",
                                   learning_rate=0.03 )
                     fit <- lgb.train(
                       params = params,   
                       data=dtrain,   
                       nrounds = 500 ,
                       verbose = -1) 
          for (i in 1:ncol(predparent_gen_m)) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
            X <- ha1
			colnames(X)<-NULL
            yhat <- predict(fit, X)
            ynew <- c(ynew, yhat)
            test_name <- paste(pred_name_m[i], pred_name_f,sep = "/")
            phe_name <- c(phe_name, test_name)
          }
          ynew <- as.matrix(ynew)
          row.names(ynew) <- phe_name
          return(ynew)
        }
		if (model == "A-P") {
          x <- cbind(gena,inbredphe)
		  colnames(x)<-NULL
		  dtrain <- lgb.Dataset(data=x,label=y)
                     params = list(boosting_type = 'gbdt',objective="regression",
                                   learning_rate=0.03 )
                     fit <- lgb.train(
                       params = params,   
                       data=dtrain,   
                       nrounds = 500 ,
                       verbose = -1) 
          for (i in 1:ncol(predparent_gen_m)) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
			p1<-t((predparent_phe_m[,i] + predparent_phe_f)/2)
			p2<-t(abs((predparent_phe_m[,i] - predparent_phe_f)/2))
            X <- cbind(ha1,p1,p2)
			colnames(X)<-NULL
            yhat <- predict(fit, X)
            ynew <- c(ynew, yhat)
            test_name <- paste(pred_name_m[i], pred_name_f,sep = "/")
            phe_name <- c(phe_name, test_name)
          }
          ynew <- as.matrix(ynew)
          row.names(ynew) <- phe_name
          return(ynew)
        }
      }
	  print("Predict by LightGBM ...")
            if (model == "A") {
                print("additive model")
                predict_lightgbm <- predict.lightgbm(model = "A")
                Results <- predict_lightgbm
            }
            else if (model == "A-P") {
                print("additive-phenotypic model")
                predict_lightgbm <- predict.lightgbm(model = "A-P")
                Results <- predict_lightgbm
            }
            else if (model == "AD-P") {
                print("additive-dominance-phenotypic model")
                predict_lightgbm <- predict.lightgbm(model = "AD-P")
                Results <- predict_lightgbm
            }
			else {
                print("additive-dominance model")
                predict_rf <- predict.lightgbm(model = "AD")
                Results <- predict_lightgbm
            }
            print("Predict by LightGBM ...ended")
        }
    if (method == "EN") {
      predict.EN <- function(model = NULL) {
        ynew <- NULL
        phe_name <- NULL
        if (model == "AD") {
          x <- cbind(gena, gend)
          fit0 <- cv.glmnet(x = x, y = y, alpha = 0.5)
          lambda <- fit0$lambda.min
          ffit <- glmnet(x = x, y = y, lambda = lambda, 
                         alpha = 0.5)
          for (i in 1:ncol(predparent_gen_m)) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
            ha2 <- t(abs((predparent_gen_m[, i] - predparent_gen_f)/2))
            X <- cbind(ha1, ha2)
            yhat <- predict(ffit, newx = X)
            ynew <- c(ynew, yhat)
            test_name <- paste(pred_name_m[i], pred_name_f, sep = "/")
            phe_name <- c(phe_name, test_name)
          }
          ynew <- as.matrix(ynew)
          row.names(ynew) <- phe_name
          return(ynew)
        }
		if (model == "AD-P") {
          x <- cbind(gena, gend,inbredphe)
          fit0 <- cv.glmnet(x = x, y = y, alpha = 0.5)
          lambda <- fit0$lambda.min
          ffit <- glmnet(x = x, y = y, lambda = lambda, 
                         alpha = 0.5)
          for (i in 1:ncol(predparent_gen_m)) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
            ha2 <- t(abs((predparent_gen_m[, i] - predparent_gen_f)/2))
			p1<-t((predparent_phe_m[,i] + predparent_phe_f)/2)
			p2<-t(abs((predparent_phe_m[,i] - predparent_phe_f)/2))
            X <- cbind(ha1,ha2,p1,p2)
            yhat <- predict(ffit, newx = X)
            ynew <- c(ynew, yhat)
            test_name <- paste(pred_name_m[i], pred_name_f, sep = "/")
            phe_name <- c(phe_name, test_name)
          }
          ynew <- as.matrix(ynew)
          row.names(ynew) <- phe_name
          return(ynew)
        }
		if (model == "A-P") {
          x <- cbind(gena,inbredphe)
          fit0 <- cv.glmnet(x = x, y = y, alpha = 0.5)
          lambda <- fit0$lambda.min
          ffit <- glmnet(x = x, y = y, lambda = lambda, 
                         alpha = 0.5)
          for (i in 1:ncol(predparent_gen_m)) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
			p1<-t((predparent_phe_m[,i] + predparent_phe_f)/2)
			p2<-t(abs((predparent_phe_m[,i] - predparent_phe_f)/2))
            X <- cbind(ha1,p1,p2)
            yhat <- predict(ffit, newx = X)
            ynew <- c(ynew, yhat)
            test_name <- paste(pred_name_m[i], pred_name_f, sep = "/")
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
          ffit <- glmnet(x = x, y = y, lambda = lambda, 
                         alpha = 0.5)
          for (i in 1:ncol(predparent_gen_m)) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
            X <- ha1
            yhat <- predict(ffit, newx = X)
            ynew <- c(ynew, yhat)
            test_name <- paste(pred_name_m[i], pred_name_f, sep = "/")
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
        predict_EN <- predict.EN(model = "A")
        Results <- predict_EN
      }
	  else if (model == "A-P") {
       print("additive-phenotypic model")
	   predict_EN <- predict.EN(model = "A-P")
        Results <- predict_EN
      }
	  else if (model == "AD-P") {
       print("additive-dominance-phenotypic model")
	   predict_EN <- predict.EN(model = "AD-P")
        Results <- predict_EN
      } 
      else {
        print("additive-dominance model")
        predict_EN <- predict.EN(model = "AD")
        Results <- predict_EN
      }
      print("Predict by EN ...ended")
    }
	
    if (method == "RKHS") {
      predict.rkhsmk <- function(model = NULL) {
        ynew <- NULL
        phe_name <- NULL
        if (model == "AD") {
          for (i in 1:ncol(predparent_gen_m)) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
            X1 <- rbind(gena, ha1)
            ha2 <- t(abs((predparent_gen_m[, i] - predparent_gen_f)/2))
            X2 <- rbind(gend, ha2)
            yNa <- as.matrix(rep(NA, nrow(ha1)))
            yNA <- c(y, yNa)
            yNA <- as.matrix(yNA)
			M <- scale(X1, center = TRUE, scale = TRUE)
            D <- (as.matrix(dist(M, method = "euclidean"))^2)/ncol(X1)
            M1 <- scale(X2, center = TRUE, scale = TRUE)
            D1 <- (as.matrix(dist(M1, method = "euclidean"))^2)/ncol(X2)
            h <- 0.5 * c(1/5, 1, 5)
            ETA <- list(list(K = exp(-h[1] * D), model = "RKHS"), 
            list(K = exp(-h[2] * D), model = "RKHS"), list(K = exp(-h[3] * D), model = "RKHS"),list(K = exp(-h[1] * D1), model = "RKHS"), 
            list(K = exp(-h[2] * D1), model = "RKHS"), list(K = exp(-h[3] * D1), model = "RKHS"))
            fm <- BGLR(y = yNA, ETA = ETA, verbose = F)
            yhat <- fm$yHat
            y <- as.matrix(y)
            yhat1 <- yhat[-c(1:nrow(y))]
            ynew <- c(ynew, yhat1)
            test_name <- paste(pred_name_m[i], pred_name_f, sep = "/")
            phe_name <- c(phe_name, test_name)
          }
          ynew <- as.matrix(ynew)
          row.names(ynew) <- phe_name
          return(ynew)
        }
        if (model == "AD-P") {
          for (i in 1:ncol(predparent_gen_m)) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
            X1 <- rbind(gena, ha1)
            ha2 <- t(abs((predparent_gen_m[, i] - predparent_gen_f)/2))
            X2 <- rbind(gend, ha2)
			p1<-t((predparent_phe_m[,i] + predparent_phe_f)/2)
			p2<-t(abs((predparent_phe_m[,i] - predparent_phe_f)/2))
			P<-cbind(p1,p2)
			P1<-rbind(inbredphe,P)
            yNa <- as.matrix(rep(NA, nrow(ha1)))
            yNA <- c(y, yNa)
            yNA <- as.matrix(yNA)
			M <- scale(X1, center = TRUE, scale = TRUE)
            D <- (as.matrix(dist(M, method = "euclidean"))^2)/ncol(X1)
		    M1 <- scale(X2, center = TRUE, scale = TRUE)
            D1 <- (as.matrix(dist(M1, method = "euclidean"))^2)/ncol(X2)
            h <- 0.5 * c(1/5, 1, 5)
            ETA <- list(list(K = exp(-h[1] * D), model = "RKHS"), 
            list(K = exp(-h[2] * D), model = "RKHS"), list(K = exp(-h[3] * D), model = "RKHS"),list(K = exp(-h[1] * D1), model = "RKHS"), 
            list(K = exp(-h[2] * D1), model = "RKHS"), list(K = exp(-h[3] * D1), model = "RKHS"), list(X = P1, model = "BayesB"))
            fm <- BGLR(y = yNA, ETA = ETA, verbose = F)
            yhat <- fm$yHat
            y <- as.matrix(y)
            yhat1 <- yhat[-c(1:nrow(y))]
            ynew <- c(ynew, yhat1)
            test_name <- paste(pred_name_m[i], pred_name_f, sep = "/")
            phe_name <- c(phe_name, test_name)
          }
          ynew <- as.matrix(ynew)
          row.names(ynew) <- phe_name
          return(ynew)
        }
		if (model == "A-P") {
          for (i in 1:ncol(predparent_gen_m)) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
            X1 <- rbind(gena, ha1)
            X <- X1
			p1<-t((predparent_phe_m[,i] + predparent_phe_f)/2)
			p2<-t(abs((predparent_phe_m[,i] - predparent_phe_f)/2))
			P<-cbind(p1,p2)
			P1<-rbind(inbredphe,P)
            yNa <- as.matrix(rep(NA, nrow(ha1)))
            yNA <- c(y, yNa)
            yNA <- as.matrix(yNA)
			M <- scale(X1, center = TRUE, scale = TRUE)
            D <- (as.matrix(dist(M, method = "euclidean"))^2)/ncol(X)
            h <- 0.5 * c(1/5, 1, 5)
            ETA <- list(list(K = exp(-h[1] * D), model = "RKHS"), 
            list(K = exp(-h[2] * D), model = "RKHS"), list(K = exp(-h[3] * D), model = "RKHS"), list(X = P1, model = "BayesB"))
            fm <- BGLR(y = yNA, ETA = ETA, verbose = F)
            yhat <- fm$yHat
            y <- as.matrix(y)
            yhat1 <- yhat[-c(1:nrow(y))]
            ynew <- c(ynew, yhat1)
            test_name <- paste(pred_name_m[i], pred_name_f, sep = "/")
            phe_name <- c(phe_name, test_name)
          }
          ynew <- as.matrix(ynew)
          row.names(ynew) <- phe_name
          return(ynew)
        }
        if (model == "A") {
          for (i in 1:ncol(predparent_gen_m)) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
            X1 <- rbind(gena, ha1)
            X <- X1
            yNa <- as.matrix(rep(NA, nrow(ha1)))
            yNA <- c(y, yNa)
            yNA <- as.matrix(yNA)
            M <- scale(X, center = TRUE, scale = TRUE)
            D <- (as.matrix(dist(M, method = "euclidean"))^2)/ncol(X)
            h <- 0.5 * c(1/5, 1, 5)
            ETA <- list(list(K = exp(-h[1] * D), model = "RKHS"), 
                        list(K = exp(-h[2] * D), model = "RKHS"), 
                        list(K = exp(-h[3] * D), model = "RKHS"))
            fm <- BGLR(y = yNA, ETA = ETA, verbose = F)
            yhat <- fm$yHat
            y <- as.matrix(y)
            yhat1 <- yhat[-c(1:nrow(y))]
            ynew <- c(ynew, yhat1)
            test_name <- paste(pred_name_m[i], pred_name_f, sep = "/")
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
        predict_rkhsmk <- predict.rkhsmk(model = "A")
        Results <- predict_rkhsmk
      }
	  else if (model == "A-P") {
         print("additive-phenotypic model")
		predict_rkhsmk <- predict.rkhsmk(model = "A")
        Results <- predict_rkhsmk
      } 
	  else if (model == "AD-P") {
         print("additive-dominance-phenotypic model")
		predict_rkhsmk <- predict.rkhsmk(model = "A")
        Results <- predict_rkhsmk
      } 
      else {
        print("additive-dominance model")
        predict_rkhsmk <- predict.rkhsmk(model = "AD")
        Results <- predict_rkhsmk
      }
      print("Predict by RKHS ...ended")
    }
    if (method == "XGBoost") {
      predict.xgboost <- function(model = NULL) {
        if (!requireNamespace("xgboost", quietly = TRUE)) {
          stop("xgboost needed for this function to work. Please install it.", 
               call. = FALSE)
        }
        phe_name <- NULL
        ynew <- NULL
        if (model == "AD") {
          x <- cbind(gena, gend)
          xg <- xgboost(data = x, label = y,
                                   colsample_bytree=0.9,
                                   eta=0.02,
                                   min_child_weight=11,
                                   nrounds=1150,
                                   subsample=0.8,
                                   set.seed(123),verbose = FALSE)
          for (i in 1:ncol(predparent_gen_m)) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
            ha2 <- t(abs((predparent_gen_m[, i] - predparent_gen_f)/2))
            X <- cbind(ha1, ha2)
            yhat <- predict(xg, X)
            ynew <- c(ynew, yhat)
            test_name <- paste(pred_name_m[i], pred_name_f,sep = "/")
            phe_name <- c(phe_name, test_name)
          }
          ynew <- as.matrix(ynew)
          row.names(ynew) <- phe_name
          return(ynew)
        }
		if (model == "AD-P") {
          x <- cbind(gena, gend,inbredphe)
          xg <- xgboost(data = x, label = y,
                                   colsample_bytree=0.9,
                                   eta=0.02,
                                   min_child_weight=11,
                                   nrounds=1150,
                                   subsample=0.8,
                                   set.seed(123),verbose = FALSE)
          for (i in 1:ncol(predparent_gen_m)) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
            ha2 <- t(abs((predparent_gen_m[, i] - predparent_gen_f)/2))
			p1<-t((predparent_phe_m[,i] + predparent_phe_f)/2)
			p2<-t(abs((predparent_phe_m[,i] - predparent_phe_f)/2))
            X <- cbind(ha1, ha2,p1,p2)
            yhat <- predict(xg, X)
            ynew <- c(ynew, yhat)
            test_name <- paste(pred_name_m[i], pred_name_f,sep = "/")
            phe_name <- c(phe_name, test_name)
          }
          ynew <- as.matrix(ynew)
          row.names(ynew) <- phe_name
          return(ynew)
        }
		if (model == "A-P") {
          x <- cbind(gena,inbredphe)
          xg <- xgboost(data = x, label = y,
                                   colsample_bytree=0.9,
                                   eta=0.02,
                                   min_child_weight=11,
                                   nrounds=1150,
                                   subsample=0.8,
                                   set.seed(123),verbose = FALSE)
          for (i in 1:ncol(predparent_gen_m)) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
			p1<-t((predparent_phe_m[,i] + predparent_phe_f)/2)
			p2<-t(abs((predparent_phe_m[,i] - predparent_phe_f)/2))
            X <- cbind(ha1,p1,p2)
            yhat <- predict(xg, X)
            ynew <- c(ynew, yhat)
            test_name <- paste(pred_name_m[i], pred_name_f,sep = "/")
            phe_name <- c(phe_name, test_name)
          }
          ynew <- as.matrix(ynew)
          row.names(ynew) <- phe_name
          return(ynew)
        }
        if (model == "A") {
          x <- gena
          xg <- xgboost(data = x, label = y, nrounds = 1000,eta = 0.07)
          for (i in 1:(ncol(predparent_gen_m))) {
            ha1 <- t((predparent_gen_m[, i] + predparent_gen_f)/2)
            X <- ha1
            yhat <- predict(xg, X)
            ynew <- c(ynew, yhat)
            test_name <- paste(pred_name_m[i], pred_name_f,sep = "/")
            phe_name <- c(phe_name, test_name)
          }
          ynew <- as.matrix(ynew)
          row.names(ynew) <- phe_name
          return(ynew)
        }
      }
      print("Predict by XGBoost...")
      if (model == "A") {
        print("additive model")
        predict_xgboost <- predict.xgboost(model = "A")
        Results <- predict_xgboost
      }
	  else if (model == "A-P") {
        print("additive-phenotypic model")
	    predict_xgboost <- predict.xgboost(model = "A")
        Results <- predict_xgboost
      }
	  else if (model == "AD-P") {
        print("additive-dominance-phenotypic model")
	    predict_xgboost <- predict.xgboost(model = "A")
        Results <- predict_xgboost
      }
      else {
        print("additive-dominance model")
        predict_xgboost <- predict.xgboost(model = "AD")
        Results <- predict_xgboost
      }
      print("Predict by XGBoost...ended")
    }
    if (select == "all") {
      Results_select <- as.data.frame(Results)
      colnames(Results_select) <- paste("all_", nrow(Results), 
                                        sep = "")
    }
    else if (select == "top") {
      Results_select <- as.data.frame(sort(Results[, 1], 
                                           decreasing = T)[c(1:number)])
      names(Results_select) <- paste("top_", number, sep = "")
    }
    else if (select == "bottom") {
      Results_select <- as.data.frame(sort(Results[, 1], 
                                           decreasing = F)[c(1:number)])
      colnames(Results_select) <- paste("bottom_", number, 
                                        sep = "")
    }
    return(Results_select)
  }
  else {
    stop("Please choose a predict method: GBLUP, BayesB, RKHS, PLS, LASSO, EN, XGBoost, LightGBM.")
  }
  return(Results)
}
