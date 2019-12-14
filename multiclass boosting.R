multi.adaboost <- function(formula,data,boos=TRUE,mfinal=100,coeflearn="Zhu",priors,cost.matrix,maxd=1){  
  if (!(as.character(coeflearn) %in% c("Freund","Breiman","Zhu"))) {
    stop("coeflearn must be 'Freund','Breiman' or 'Zhu' ")
  }
  formula <- as.formula(formula)
  vardep <- data[,as.character(formula[[2]])]
  n <- length(data[,1])
  nclases <- nlevels(vardep)
  pesos <- rep(1/n,n)
  guardarpesos <- array(0,c(n,mfinal))
  w <- rep(1/n,n)
  data <- cbind(pesos,data)
  arboles <- list()
  pond <- rep(0,mfinal)
  pred <- data.frame(rep(0,n))
  arboles[[1]] <- rpart(formula,data=data[,-1],control=rpart.control(minsplit=1,
                                                                            cp=-1,maxdepth=30))
  nvar <- dim(varImp(arboles[[1]],surrogates=FALSE,competes=FALSE))[1]
  imp <- array(0,c(mfinal,nvar))
  for (m in 1:mfinal) {
    if (boos == TRUE) {
      k <- 1
      while (k == 1) {
        boostrap <- sample(1:n,replace=TRUE,prob=pesos)
        fit <- rpart(formula,data=data[boostrap,-1],parms=list(prior=priors,loss=cost.matrix),
                     control=rpart.control(maxdepth=maxd))
        k <- length(fit$frame$var)
      }
      flearn <- predict(fit,newdata=data[,-1],type="class")
      ind <- as.numeric(vardep != flearn)
      err <- sum(pesos * ind)
    }
    if (boos == FALSE) {
      w <<- pesos
      fit <- rpart(formula=formula,data=data[,-1],
                   weights=w,parms=list(prior=priors,loss=cost.matrix),control=rpart.control(maxdepth=maxd))
      flearn <- predict(fit,data=data[,-1],type="class")
      ind <- as.numeric(vardep != flearn)
      err <- sum(pesos * ind)
    }
    c <- log((1-err)/err)
    if (coeflearn == "Breiman") {
      c <- (1/2) * c
    }
    if (coeflearn == "Zhu") {
      c <- c+log(nclases-1)
    }
    guardarpesos[,m] <- pesos
    pesos <- pesos * exp(c * ind)
    pesos <- pesos/sum(pesos)
    maxerror <- 0.5
    eac <- 0.001
    if (coeflearn == "Zhu") {
      maxerror <- 1-1/nclases
    }
    if (err >= maxerror) {
      pesos <- rep(1/n,n)
      maxerror <- maxerror-eac
      c <- log((1-maxerror)/maxerror)
      if (coeflearn == "Breiman") {
        c <- (1/2) * c
      }
      if (coeflearn == "Zhu") {
        c <- c+log(nclases-1)
      }
    }
    if (err == 0) {
      pesos <- rep(1/n,n)
      c <- log((1-eac)/eac)
      if (coeflearn == "Breiman") {
        c <- (1/2) * c
      }
      if (coeflearn == "Zhu") {
        c <- c+log(nclases-1)
      }
    }
    arboles[[m]] <- fit
    pond[m] <- c
    if (m == 1) {
      pred <- flearn
    }
    else {
      pred <- data.frame(pred,flearn)
    }
    if (length(fit$frame$var) > 1) {
      k <- varImp(fit,surrogates=FALSE,competes=FALSE)
      imp[m,] <- k[sort(row.names(k)),]
    }
    else {
      imp[m,] <- rep(0,nvar)
    }
  }
  classfinal <- array(0,c(n,nlevels(vardep)))
  for (i in 1:nlevels(vardep)) {
    classfinal[,i] <- matrix(as.numeric(pred == levels(vardep)[i]),
                              nrow=n) %*% as.vector(pond)
  }
  predclass <- rep("O",n)
  predclass[] <- apply(classfinal,1,FUN=which.max)
  imppond <- as.vector(as.vector(pond) %*% imp)
  imppond <- imppond/sum(imppond) * 100
  names(imppond) <- sort(row.names(k))
  votosporc <- classfinal/apply(classfinal,1,sum)
  ans <- list(formula=formula,trees=arboles,weights=pond,
              votes=classfinal,prob=votosporc,class=predclass,
              importance=imppond)
  attr(ans,"vardep.summary") <- summary(vardep,maxsum=700)
  mf <- model.frame(formula=formula,data=data[,-1])
  terms <- attr(mf,"terms")
  ans$terms <- terms
  ans$call <- match.call()
  class(ans) <- "boosting"
  ans
}

