###############################################################################
# the following functions are to implement the parameter estimation of the    #
# AOPR and ridge regression.                                                  #
###############################################################################


#' @description This function is used for cross-validation and not directly used\r
#'  by users
#' @param nfold The number of fold in cross validation.
#' @param N The sample size
CV_testlist <- function(nfold = 5,N){
  x <- sample(N,N)
  testlist <- list()
  for (i in 1:nfold) {
    ntest <- round(length(x)/(nfold-i+1))
    testlist[[i]] <- x[1:ntest]
    x <- x[-(1:ntest)]
  }
  testlist
}

#' @description This function selects the appropriate penalized coefficient by prefixed 
#' coefficient, iterative estimation, cross-validation, generalized cross-validation and 
#' OLS estimation, and implements the parameter estimator of ridge regression.
#' @author 
#' @param formula An object of class "formula".
#' @param data A data frame. 
#' @param lambda A non-negative numeric referring to the penalized coefficient.
#' @param lambdatype A character indicating the method of selecting lambda.\r
#' lambdatype must be one of 'prefixed','gcv','iterative','ols','cv'.
#' @param lambdas A vector of positive numeric representing the candidate lambda.
#' @param maxiter The maximal iterative numbers when lambdatype is "iterative".
#' @param k The number of fold for cross validation.
#' @param seed The random seed for the random sample in cross-validation.
ridge <- function(formula,data,lambda = 0,lambdatype,maxiter = 100,
                  lambda_sep = 1e-10,lambdas = seq(0,1,0.1),k = 10,seed=95){
  if(missing(lambda) & lambdatype == "prefixed") 
    stop("When lambdatype is 'prefixed', lambda can not be missed") 
  if(!missing(lambda) & !lambdatype == "prefixed") 
    stop("When lambda is not missed, lambdatype must be 'prefixed'") 
  if(!lambdatype %in% c("prefixed","gcv","iterative","ols","cv"))
    stop("lambdatype must be one of 'prefixed','gcv','iterative','ols','cv'")

  dat0 <- as.matrix(model.frame(formula,data)) 
  sigma <- apply(dat0, 2, sd) 
  miu <- apply(dat0, 2, mean) 
  sxy <- sigma[1]/sigma[-1]   
  sddat <- scale(dat0)        
  y <- sddat[,1];  x <- sddat[,-1]; nx <- dim(x)[2]; n <- dim(x)[1]
  txx <- t(x)%*%x
  residue <- lm(y~x)
  residue <- var(residue$residuals)*(n-1)/(n-nx)
  
  eigens <- eigen(txx)
  eigenvector <- eigens[[2]]
  eigenvalue <- eigens[[1]]
  
  gcvs <- cvs <- NULL
  if(lambdatype == "gcv"){
    lambda_iterate <- beta_iterate <- NULL
    for (i in lambdas) {
      sx <- solve(txx+diag(rep(i,nx))) %*% t(x)
      pk <- sum(diag(x%*%sx)) 
      RSSi <- sum((y-x%*%sx%*%y)^2)
      gcvi <- RSSi/(n*(1-pk/n)^2)
      gcvs <- c(gcvs,gcvi)
    }
    lambda <- lambdas[which.min(gcvs)]
  } 
  else if(lambdatype == "cv"){
    lambda_iterate <- beta_iterate <- NULL
    if(!missing(seed)) set.seed(seed)
    folds <- CV_testlist(nfold = k,N = n)
    y_cv <- dat0[,1][unlist(folds)]
    for (i in lambdas) {
      predict_values <- NULL
      for(j in 1:k){
        fold_test <- matrix(dat0[folds[[j]],],ncol = nx + 1)      
        fold_train <- dat0[-folds[[j]],]  

        tr_sigma <- apply(fold_train, 2, sd) 
        tr_miu <- apply(fold_train, 2, mean) 
        tr_sxy <- tr_sigma[1]/tr_sigma[-1] 
        tr_sddat <- scale(fold_train) 
        tr_y <- tr_sddat[,1];  tr_x <- tr_sddat[,-1]
        tr_nx <- dim(tr_x)[2]; tr_n <- dim(tr_x)[1]
        tr_txx <- t(tr_x)%*%tr_x
        
        beta <- solve(tr_txx+diag(rep(i,tr_nx))) %*% t(tr_x) %*%tr_y 
        beta <- beta[,1] * tr_sxy 
        intercept <- tr_miu[1]-sum(tr_miu[-1]*beta) 
        names(intercept) <- "intercept"
        predict_value <- fold_test[,-1] %*% beta + intercept 
        predict_values <- c(predict_values,predict_value)
      }
      cvi <- sqrt(mean((predict_values - y_cv)^2))
      cvs <- c(cvs,cvi)
    }
    lambda <- lambdas[which.min(cvs)]
  }
  else if(lambdatype == "prefixed"){
    lambda <- lambda
    lambda_iterate <- beta_iterate <- NULL
  } 
  else {
    maxiter <- ifelse(lambdatype=="iterative",maxiter,1)
    lambdai <- 0;lambda_iterate <- 0;beta_iterate <- NULL
    for (i in 1:maxiter) {
      betai <- solve(txx+diag(rep(lambdai,nx))) %*% t(x) %*%y
      alpha <- abs(t(eigenvector)%*%betai[,1])
      maxalpha <- max(alpha)
      ki <- residue/maxalpha^2
      beta_iterate <- rbind(beta_iterate,betai[,1]*sxy) 
      if(abs(lambdai-ki) < lambda_sep)  break() 
      else lambdai <- ki
      if(lambdai > 1e+10) break()
      lambda_iterate <- c(lambda_iterate,lambdai)
    }
    lambda <- lambdai
  } 
  beta <- solve(txx+diag(nx)*lambda) %*% t(x) %*%y
  beta <- beta[,1] * sxy
  intercept <- miu[1]-sum(miu[-1]*beta)
  names(intercept) <- "intercept"
  fitvalue <- dat0[,-1] %*% beta + intercept
  R2 <- 1-var(fitvalue[,1]-dat0[,1])/sigma[1]^2
  M <- diag(sxy)%*%solve(txx+diag(rep(lambda,nx)))%*%txx
  cov <- residue*M%*%solve(txx)%*%t(M)
  list(beta = beta,intercept = intercept,fitvalue = fitvalue[,1],cov = cov,R2 = R2,
       lambda = lambda,cvs = cvs,gcvs = gcvs,lambda_iterate = lambda_iterate,
       beta_iterate=beta_iterate)
}

###for example using ridge -------
# n <- 20
# x1 <- runif(n,-1,1)
# x2 <- x1 + rnorm(n)
# x3 <- x1 + rnorm(n)
# x4 <- x1 + rnorm(n)
# y <- 0.1 * x1 + 0.2*x2 + 0.3*x3 + 0.4*x4 + rnorm(n,0,4)
# data0 <- data.frame(x1 = x1,x2 = x2,x3 = x3,x4 = x4,y = y)
# cor(data0)
# formula0 <- y ~ x1 + x2 + x3 + x4
# res <- ridge(formula = formula0,data = data0,lambdatype = "cv",lambdas = 0:20,k = 10,seed = 5)
### end ridge example ------------



#' @description This function selects the appropriate penalized coefficient by prefixed 
#' coefficient, iterative estimation, cross-validation, generalized cross-validation and OLS estimation, 
#' and implements the parameter estimator of average ordinary least squares centered penalized regression.
#' @author 
#' @param formula An object of class "formula".
#' @param data A data frame. 
#' @param lambda A non-negative numeric referring to the penalized coefficient.
#' @param lambdatype A character indicating the method of selecting lambda.\r
#' lambdatype must be one of 'prefixed','gcv','iterative','ols','cv'.
#' @param lambdas A vector of positive numeric representing the candidate lambda.
#' @param maxiter The maximal iterative numbers when lambdatype is "iterative".
#' @param k The number of fold for cross validation.
#' @param seed The random seed for the random sample in cross-validation.
aopr <- function(formula,data,lambda = 0,lambdatype,maxiter = 100,
                 lambda_sep = 1e-10,lambdas = seq(0,1,0.1),k = 10,seed=95){
  
  if(missing(lambda) & lambdatype == "prefixed") 
    stop("When lambdatype is 'prefixed', lambda can not be missed") 
  if(!missing(lambda) & !lambdatype == "prefixed") 
    stop("When lambda is missed, lambdatype must be 'prefixed'") 
  if(!lambdatype %in% c("prefixed","gcv","iterative","ols","cv"))
    stop("lambdatype must be one of 'prefixed','gcv','iterative','ols','cv'")
  
  dat0 <- as.matrix(model.frame(formula,data)) 
  sigma <- apply(dat0, 2, sd) 
  miu <- apply(dat0, 2, mean) 
  sxy <- sigma[1]/sigma[-1] 
  sddat <- scale(dat0) 
  y <- sddat[,1];  x <- sddat[,-1]; nx <- dim(x)[2]; n <- dim(x)[1]
  txx <- t(x)%*%x
  residue <- lm(y~x)
  residue <- var(residue$residuals)*(n-1)/(n-nx)
  
  eigens <- eigen(txx)
  eigenvector <- eigens[[2]]
  eigenvalue <- eigens[[1]]
  
  sumx <- apply(x, 1, sum)
  sop <- sum(sumx * y)/sum(sumx^2)
  centre <- rep(sop,nx)
  d <- rep(1,nx)
  h <- t(eigenvector)%*%d%*%solve(t(d)%*%txx%*%d)%*%t(d)%*%txx%*%eigenvector
  
  gcvs <- cvs <- NULL
  if(lambdatype == "gcv"){
    lambda_iterate <- beta_iterate <- NULL
    for (i in lambdas) {
      sx <- solve(txx+diag(rep(i,nx))) %*% (t(x)+i*d%*%solve(t(d)%*%txx%*%d)%*%t(d)%*%t(x))
      pk <- sum(diag(x%*%sx))  
      RSSi <- sum((y-x%*%sx%*%y)^2)
      gcvi <- RSSi/(n*(1-pk/n)^2)
      gcvs <- c(gcvs,gcvi)
    }
    lambda <- diag(rep(lambdas[which.min(gcvs)],nx)) 
  } else if(lambdatype == "cv"){
    lambda_iterate <- beta_iterate <- NULL
    if(!missing(seed)) set.seed(seed)
    folds <- CV_testlist(nfold = k,N = n) 
    y_cv <- dat0[,1][unlist(folds)]
    for (i in lambdas) {
      predict_values <- NULL
      for(j in 1:k){
        fold_test <- matrix(dat0[folds[[j]],],ncol = nx + 1)      
        fold_train <- dat0[-folds[[j]],]  
        
        tr_sigma <- apply(fold_train, 2, sd) 
        tr_miu <- apply(fold_train, 2, mean) 
        tr_sxy <- tr_sigma[1]/tr_sigma[-1] 
        tr_sddat <- scale(fold_train) 
        tr_y <- tr_sddat[,1];  tr_x <- tr_sddat[,-1]; tr_nx <- dim(tr_x)[2]; tr_n <- dim(tr_x)[1]
        tr_txx <- t(tr_x)%*%tr_x
        
        tr_sumx <- apply(tr_x, 1, sum)
        tr_sop <- sum(tr_sumx * tr_y)/sum(tr_sumx^2)
        tr_centre <- rep(tr_sop,tr_nx)
        
        lambda <- diag(rep(i,tr_nx))
        beta <- solve(tr_txx+lambda) %*% (t(tr_x) %*%tr_y+lambda%*%tr_centre);
        beta <- beta[,1] * tr_sxy
        intercept <- tr_miu[1]-sum(tr_miu[-1]*beta)
        names(intercept) <- "intercept"
        predict_value <- fold_test[,-1] %*% beta + intercept
        predict_values <- rbind(predict_values,predict_value)
      }
      cvi <- sqrt(mean((predict_values - y_cv)^2))
      cvs <- c(cvs,cvi)
    }
    lambda <- diag(rep(lambdas[which.min(cvs)],nx))
  } else if(lambdatype == "prefixed"){
    lambda_iterate <- beta_iterate <- NULL
    lambda <- diag(rep(lambda,nx))
  } else {
    maxiter <- ifelse(lambdatype=="iterative",maxiter,1)
    sigmabi <- residue/eigenvalue 
    lambdai <- rep(0,nx);lambda_iterate <- lambdai;beta_iterate <- NULL
    pcentre <- (t(eigenvector)%*%centre)[,1]
    for (i in 1:maxiter) {
      betai <- (solve(diag(eigenvalue)+diag(lambdai)) %*% (t(x%*%eigenvector)%*%y+lambdai*pcentre))[,1]
      bi <- betai 
      f <- apply(h, 1, function(x) sum(x^2*sigmabi)+sum(bi*x)^2)
      hb <- (h%*%bi)[,1]
      s <- bi^2 + f - 2*bi*hb
      ki <- ((1-diag(h))*eigenvalue*residue)/(eigenvalue*s-diag(h)*residue)
      ki[] <- max(min(ki),0)
      beta_iterate <- rbind(beta_iterate,betai)
      if(sum(abs(lambdai- ki)) < lambda_sep*nx |sum(abs(ki)) > 1e+6 )  break()
      else lambdai <- ki
      lambda_iterate <- rbind(lambda_iterate,lambdai)
    }
    lambda <- eigenvector%*%diag(ki)%*%t(eigenvector)
  } 
  
  beta <- solve(txx+lambda) %*% (t(x) %*%y+lambda%*%centre)
  beta <- beta[,1] * sxy
  intercept <- miu[1]-sum(miu[-1]*beta)
  names(intercept) <- "intercept"
  fitvalue <- dat0[,-1] %*% beta + intercept
  R2 <- 1-var(fitvalue[,1]-dat0[,1])/sigma[1]^2
  M <- diag(sxy)%*%solve(txx+lambda)%*%(txx+lambda%*%d%*%solve(t(d)%*%txx%*%d)%*%t(d)%*%txx)
  cov <- residue*M%*%solve(txx)%*%t(M)
  list(beta = beta,intercept = intercept,fitvalue = fitvalue[,1],cov = cov,R2 = R2,
       lambda = lambda,cvs = cvs,gcvs = gcvs,lambda_iterate = lambda_iterate,
       beta_iterate=beta_iterate)
}

# for example using aopr -------
# n <- 20
# x1 <- runif(n,-1,1)
# x2 <- x1 + rnorm(n)
# x3 <- x1 + rnorm(n)
# x4 <- x1 + rnorm(n)
# y <- 0.1 * x1 + 0.2*x2 + 0.3*x3 + 0.4*x4 + rnorm(n,0,4)
# data0 <- data.frame(x1 = x1,x2 = x2,x3 = x3,x4 = x4,y = y)
# cor(data0)
# formula0 <- y ~ x1 + x2 + x3 + x4
# res <- aopr(formula = formula0,data = data0,lambdatype = "cv",lambdas = 0:20,k = 10,seed = 5)
### end aopr example ------------