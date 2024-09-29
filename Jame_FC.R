smoothing.cluts <- function(data=stop("No data supplied"), K = 5, q = 6, h = 2, use.covariates=FALSE,stand.cov=TRUE, 
                            index.cov=NULL, random=TRUE, B=NULL, svd=TRUE, lambda=1.4e-4, EM.maxit=50, EMstep.tol=1e-6, Mstep.maxit=20,Mstep.tol=1e-4, EMplot=TRUE,trace=FALSE,n.cores=NULL){
  ## lambda=0.000140625
  ## function to fit functional clustering model to data...
  ## data is a list of at least five objects (vectors) named as 'x', 'time', 'timeindex', 'curve', 'grid':
  ##      i) suppose we observe N independent subjects different (e.g, years), each consisting of n_i observations (time points for ith subject). 'x' is a vector of length (\sum_i n_i) with the first n_1 elements representing the observations of the first subject, followed by n_2 observations of the second subject, etc;
  ##      ii) 'time' is a (\sum_i n_i) vector of time points (usually time points are scaled to [0,1] for each subject),
  ##      iii) 'timeindex' is a (\sum_i n_i) vector of time indices from T possible from 'grid'. So each observation has a corresponding location (time index) within [0,1] uniquely specified time points.
  ##      iv)  'curve' is a (\sum_i n_i) vector of integers from 1,..., N, specifying the curve (subject number) for each observation in 'x'
  ##      v)   'grid' is a T vector of all unique time points (values within [0,1] interval) for all N years, needed for estimation of the B-spline coefficients in fda::eval.basis().  'timeindex' and 'grid' together give the timepoint for each subject (curve).
  ##      vi) if supplied, 'covariates' is an \eqn{N x r} matrix (or data frame) of scalar covariates (finite-dimensional covariates).
  ## h is a positive integer, parameter vector dimension in low-dimensionality representation of the curves (spline coefficients), h should be <= K, number of clusters. (default: h=2)
  ## q is number of splines (default: q=6); 
  ## K is a number of clusters. (default: K=3);
  ## random - TRUE/FALSE,  if TRUE the cluster belongings is given by uniform distribution, otherwise k-means is used to initialize cluster belongings;
  ## B is an N x q matrix of spline coefficients,  the spline approximation of the yearly curves based on p number of splines. Usually the matrix is supplied. If you supply it, check so that it has correct dimension, correct number of splines. If B=NULL (default), the coefficients are estimated using fda:: create.bspline.basis;
  ## svd is TRUE/FALSE, whether SVD decomposition should be used for the matrix of spline coefficients;
  ## use.covariates is TRUE/FALSE, whether covariates should be used when modelling;
  ## stand.cov tells whether covariates should be standardized (centering and dividing by the standard deviation)
  ## index.cov a vector of indices, indicating the selected covariates to be used when modelling. If NULL (default) all present covariates are used
  ## lambda is a positive real number, smoothing parameter value to be used when estimating B-spline coefficients
  ## trace is logical indicating whether to print the current values of sig² and sig²_covariates at each iteration of M step.
  ## n.cores is a number of cores to be used for parallel computing
  
  if (!is.list(data)) stop("data should be supplied as a list")
  if (length(data) < 3) stop("not enough data supplied; data should be a list of at least three elements")
  if ("x" %in% names(data) ==0) stop("no 'x' object in the data supplied")
  if ("time" %in% names(data) ==0) stop("no 'time' object in the data supplied")
  if ("curve" %in% names(data) ==0) stop("no 'curve' object in the data supplied")
  if ("grid" %in% names(data) ==0) 
    # stop("no 'grid' object in the data supplied")
    data$grid <- sort(unique(data$time))
  if ("timeindex" %in% names(data) ==0){ #stop("no 'timeindex' object in the data supplied")
    fwhich <- function(x,grid) which(x==grid)
    data$timeindex <- sapply(data$time,fwhich, grid=data$grid)   
  }
  
  object <- estimate.smoothing(data=data, K = K, q = q, h = h, random=random, B=B,  svd=svd, use.covariates=use.covariates,
                               stand.cov=stand.cov,index.cov=index.cov,lambda=lambda, EM.maxit=EM.maxit, EMstep.tol=EMstep.tol, Mstep.maxit=Mstep.maxit, Mstep.tol=Mstep.tol, EMplot=EMplot,trace=trace,n.cores=n.cores) 
  
  object$pars <- list()
  if (object$initials$moc) { ## model with covariates...
    object$pars$lambda.zero <- object$parameters$lambda.zero
    object$pars$Delta <- object$parameters$tag.GammaK
    object$pars$alpha <- object$parameters$alpha
    object$pars$Lambda <- object$parameters$Lambda
    object$pars$vk <- object$parameters$tag.estimates
    colnames(object$pars$vk)<- colnames(object$data$covariates) ## selected covariates if index.cov != Null
    rownames(object$pars$vk)<- paste("cluster:",1:nrow(object$pars$vk),sep="")
    object$pars$probs <- object$parameters$tag.pi
    object$pars$sig2 <- object$parameters$sigma[1]
    object$pars$sig2x <- object$parameters$sigma[2]
  } else{ ## model without covariates...
    object$pars$lambda.zero <- object$parameters$lambda.zero
    object$pars$Gamma <- object$parameters$GammaK
    object$pars$alpha <- object$parameters$alpha
    object$pars$Lambda <- object$parameters$Lambda
    object$pars$probs <- object$parameters$pi 
    object$pars$sig2 <- object$parameters$sigma
  }
  
  parameters <- NULL ## binding the variable locally to the function, so the 'R CMD check' has nothing to complain about 
  object <- within(object, rm(parameters))
  
  object$data <- data ## to include data$grid and data$timeindex (when not supplied) to be used with plot.mocca()
  object$nobs <- object$initials$N ## needed for BIC() to work
  
  class(object)<-"mocca"
  object
} ## end mocca


###############################################
## EM-algorithm to estimate model parameters... 
###############################################

estimate.smoothing <- function(data, K = 5, q = 6, h = 2, random=TRUE, B=NULL, svd=TRUE,use.covariates=FALSE,stand.cov=TRUE,
                               index.cov=NULL, lambda=1.4e-4, EM.maxit=50, EMstep.tol=1e-8, Mstep.maxit=10,Mstep.tol=1e-4, EMplot=TRUE,
                               trace=TRUE,n.cores=NULL){
  ## function to fit mocca model to data via EM algorithm 
  ## to find maximum log likelihood model parameter estimates...
  # if (is.null(n.cores))
  #   nc <- parallel::detectCores()
  # else nc <- n.cores
  # cl <- parallel:: makeCluster(nc-1)
  # invisible(gc())
  # doParallel:: registerDoParallel(cl)
  
  ## get starting values of model parameters...
  initfit.mod <- initial.mocca(data=data, K=K, q=q,h=h, random=random, B=B, svd=svd, use.covariates=use.covariates,
                               stand.cov=stand.cov,index.cov=index.cov, lambda=lambda)
  data <- initfit.mod$data.aug
  design <- initfit.mod$design
  initials <- initfit.mod$initials
  
  ## The first "E-step" with the initial values obtained from initial.mocca()...
  vars0 <- initfit.mod$vars
  parameters0 <- initfit.mod$parameters
  score.hist <- NULL ## keeping loglik and variances at each iteration for further tracking if needed
  ll.old <- NULL
  ## The first "M-step" maximization step...
  parameters1 <- Mstep.mocca(parameters=parameters0, vars=vars0, data=data, initials=initials, design=design,
                             Mstep.maxit=Mstep.maxit, Mstep.tol=Mstep.tol,trace=trace)
  ## get current estimate of sigma...
  sigma <- variances.mocca(parameters1,vars0,data,design)
  parameters1$sigma <- sigma
  ## score.hist <- rbind(score.hist,c(sigma,loglik.EMmocca(parameters1,data,design)))
  ll.new <- loglik.EMmocca(data,parameters1,design)
  score.hist <- rbind(score.hist,c(sigma,ll=ll.new))
  dimnames(score.hist)[[2]][dim(score.hist)[2]]<-'ll'
  ## dimnames(score.hist)[[2]][1]<-'variance splines'
  ## dimnames(score.hist)[[2]][2]<-'variance covariates'
  ## dimnames(score.hist)[[2]][3]<-'variance covariates un-adjusted'
  ## dimnames(score.hist)[[2]][4]<-'trace of covariates'
  
  vars1 <- Estep.mocca(parameters1,vars0,data,design,initials) 
  ## sigma.old <- sigma[1]+sigma[2]
  ## sigma.ny <- 1
  ll.old <- ll.new+ll.new*.5
  iter <- 0    
  
  if(EMplot){
    oldpar <- par(no.readonly = TRUE)    
    on.exit(par(oldpar))   
    par(mfrow=c(3,3))
    cluster.mean <- matrix(0,K,dim(design$FullS)[1])
    for(k in 1:K) cluster.mean[k,] <- design$FullS %*% (parameters1$lambda.zero + parameters1$Lambda %*% parameters1$alpha[k,  ])
    
    ylim <- c(min(cluster.mean)-.1*min(cluster.mean), max(cluster.mean)+.1*max(cluster.mean))   
    
    if (initials$moc){ ## model with covariates...
      ll <- round(score.hist[iter+1,'ll']) ## current loglik
      sig2 <- round(score.hist[iter+1,1],3) ## current sig^2
      sig2x <- round(score.hist[iter+1,2],3) ## current sig^2 of covariates, sig2_x
      plot(data$grid,data$grid,type='n',ylab="Cluster Mean", ylim=ylim,
           xlab=bquote("ll" ==.(ll) ~ ", " ~ sigma^2 ==.(sig2) ~ "and" ~ sigma[x]^2 ==.(sig2x)) )
      legend('topright',legend=round(parameters1$tag.pi,5),lty=1,col=c(1:K)+1)
      ## plot(data$grid,data$grid,type='n',ylab="Cluster Mean",xlab=paste('ll =',round(score.hist[iter+1,'ll']),'\n sigmas =',round(score.hist[iter+1,1],3), 'and',round(score.hist[iter+1,2],3)),ylim=c(-45,55))
    } else{ ## model without covariates...
      ## plot(data$grid,data$grid,type='n',ylab="Cluster Mean", xlab=expression(\sigma^2"=" ) ylim=c(-45,55)) ##xlab=paste('ll =',round(score.hist[iter+1,'ll']),'\n sigmas =',round(score.hist[iter+1,1],3)), ylim=c(-45,55))
      ll <- round(score.hist[iter+1,'ll']) ## current loglik
      sig2 <- round(score.hist[iter+1,1],3) ## current sig^2
      plot(data$grid,data$grid,type='n',ylab="Cluster Mean", ylim=ylim,
           xlab=bquote("ll" ==.(ll) ~ "and" ~ sigma^2 ==.(sig2)) )
      legend('topright',legend=round(parameters1$pi,5),lty=1,col=c(1:K)+1)
    }
    title(paste(Sys.time(),'with iteration no:',iter))
    grid()
    
    ##legend('bottomleft',legend=paste('Rel diff =',round(abs(sigma.old - sigma.ny)/sigma.ny,5)))
    legend('bottomleft',legend=paste('loglik rel diff =',round(abs(ll.old - ll.new)/(.1 + abs(ll.new)),5)))
    lines(data$grid,design$FullS%*%parameters1$lambda.zero,lty=2,lwd=3);for(k in 1:K) lines(data$grid,cluster.mean[k,], col = k+1, lwd = 2)
  } ## end if (EMplot)
  
  
  nu <- Sys.time()
  upprepa <- 0
  
  old.prop <- parameters1$tag.pi
  ## while((abs(sigma.old - sigma.ny)/sigma.ny > EMstep.tol) & (iter <= EM.maxit)){ 
  while(abs(ll.new - ll.old)/(.1 + abs(ll.new)) > EMstep.tol & (iter <= EM.maxit)){ 
    ## checking for loglik convergence instead...
    iter <- iter+1
    ll.old <- ll.new
    parameters1 <- Mstep.mocca(parameters=parameters1,vars=vars1,data=data,initials=initials,design=design, Mstep.maxit=Mstep.maxit,  Mstep.tol=Mstep.tol,trace=trace)
    sigma <- variances.mocca(parameters1,vars1,data,design)
    parameters1$sigma <- sigma 
    ## sigma.ny <- sigma[1]+sigma[2]
    ll.new <- loglik.EMmocca(data,parameters1,design)
    ## need to check if new step increases log likelihood ??? (check not done yet)...
    if (iter > 1 & ll.new < ll.old) #print(paste("with iter=", iter, " decrease in loglik"))
      warning(paste("with iter=", iter, " decrease in loglik"))
    
    score.hist <- rbind(score.hist,c(sigma,ll=ll.new))
    vars1 <- Estep.mocca(parameters1,vars1,data,design,initials)
    
    if(EMplot){
      cluster.mean <- matrix(0,K,dim(design$FullS)[1])
      for(k in 1:K) cluster.mean[k,] <- design$FullS %*% (parameters1$lambda.zero + parameters1$Lambda %*% parameters1$alpha[k,  ])
      ylim <- c(min(cluster.mean)-.1*min(cluster.mean), max(cluster.mean)+.1*max(cluster.mean)) 
      ##plot(data$grid,data$grid,type='n',ylab="Cluster Mean",xlab=paste('ll =',round(score.hist[iter+1,'ll']),'\n sigmas =',round(score.hist[iter+1,1],3), 'and',round(score.hist[iter+1,2],3)),ylim=c(-45,55))
      if (length(sigma)==2){ ## model with covariates...
        ll <- round(score.hist[iter+1,'ll'],3) ## current loglik
        sig2 <- round(score.hist[iter+1,1],3) ## current sig^2
        sig2x <- round(score.hist[iter+1,2],3) ## current sig^2 of covariates, sig2_x
        plot(data$grid,data$grid,type='n',ylab="Cluster Mean", ylim=ylim,
             xlab=bquote("ll" ==.(ll) ~ ", " ~ sigma^2 ==.(sig2) ~ "and" ~ sigma[x]^2 ==.(sig2x)) )
        new.prop <- parameters1$tag.pi
      } else{ ## model without covariates...
        ll <- round(score.hist[iter+1,'ll'],3) ## current loglik
        sig2 <- round(score.hist[iter+1,1],3) ## current sig^2
        plot(data$grid,data$grid,type='n',ylab="Cluster Mean", ylim=ylim,
             xlab=bquote("ll" ==.(ll) ~ "and" ~ sigma^2 ==.(sig2)) )
        new.prop <- parameters1$pi
      }
      
      title(paste(format(round(Sys.time()-nu,2)),'with iteration no:',iter))
      grid()
      legend('topright',legend=round(new.prop,5),lty=1,col=c(1:K)+1)
      legend('bottomleft',legend=paste('loglik rel diff =',round(abs(ll.old - ll.new)/(.1 + abs(ll.new)),3)))
      legend('topleft',legend=paste('Num of changing curves =',round(sum(initials$N*abs(old.prop-new.prop)),2) ))
      lines(data$grid,design$FullS%*%parameters1$lambda.zero,lty=2,lwd=3);for(k in 1:K) lines(data$grid,cluster.mean[k,], col = k+1, lwd = 2)
      old.prop <- new.prop 
    } ## end if (EMplot)
    
    max.ll<-max(score.hist[,'ll'])
    if(score.hist[iter,'ll']>=max.ll){
      upprepa<-upprepa+1
      ## svar.ll<-list(par=parameters1,var=vars1,data=data,design=design,score.hist=score.hist)
      if (initials$moc){ ## model with covariates...
        if(upprepa==5) {
          klass<-apply(vars1$tag.piigivej,1,order)[K,] ## [ant.kl,]
          piigivej <- matrix(0, initials$N, K)
          piigivej[col(piigivej)==klass]<-1
          vars1$tag.piigivej<-piigivej
          parameters1$lambda.zero<-round(parameters1$lambda.zero)
          parameters1$Lambda<-round(parameters1$Lambda)
          parameters1$alpha<-round(parameters1$alpha)
          parameters1$tag.GammaK<-round(parameters1$tag.GammaK,1)
          vars1 <- Estep.mocca(parameters1,vars1,data,design,initials)
          upprepa<-0
        } else vars1$tag.piigivej<-round(vars1$tag.piigivej,1)
      } else { ## model without covariates...
        if(upprepa==5) {
          klass<-apply(vars1$piigivej,1,order)[K,] ## [ant.kl,]
          piigivej <- matrix(0, initials$N, K)
          piigivej[col(piigivej)==klass]<-1
          vars1$piigivej<-piigivej
          parameters1$lambda.zero<-round(parameters1$lambda.zero)
          parameters1$Lambda<-round(parameters1$Lambda)
          parameters1$alpha<-round(parameters1$alpha)
          parameters1$GammaK<-round(parameters1$GammaK,1)
          vars1 <- Estep.mocca(parameters1,vars1,data,design,initials)
          upprepa<-0
        } else vars1$piigivej<-round(vars1$piigivej,1)
      }
    }
    nu<-Sys.time()
  } ## end while loop
  
  # if (!is.null(cl)){
  #   parallel::stopCluster(cl) 
  #   invisible(gc())
  #   rm(cl)
  #   invisible(gc())
  # }
  
  dimnames(score.hist)[[1]]<-paste('iter no',1:dim(score.hist)[1])
  if (abs(ll.new - ll.old)/(.1 + abs(ll.new)) < EMstep.tol) {
    conv <- 0 ## full convergence
  } else if (iter > EM.maxit) conv <- 1 ## maximum number of iterations reached
  
  
  list(loglik=ll.new, sig2=sigma, conv=conv, iter=iter, parameters=parameters1,vars=vars1,data=data,design=design, score.hist=score.hist, initials=initials) 
} ## end estimate.mocca


######################################################
## starting step to initialize required parameters...
######################################################

initial.mocca <- function(data, K = 3, q = 6, h = 2, random=FALSE, B=NULL, svd=TRUE,use.covariates=FALSE,stand.cov=TRUE,
                          index.cov=NULL, lambda=1.4e-4){
  ## Initial E-step of the EM-algorithm.
  ## the function to initialize all necessary parameters: all parameters except Gamma_k. It can start in two ways. One is by running k-means several times on the initial penalized spline approximation of the data using the pre-specified penalty, lambda, and p cubic spline functions with G cluster centroids and keep the best fit. The second is to uniformly assign a cluster belongings for each object. Both methods gives the first guesses of the posterior probabilities, \hat{\pi}_{k|i}, for each object which in turn gives the start values \hat{\pi}_k, k=1,...,G, from the mean value of \hat{\pi}_{k|i}, i=1,...,N. This k-mean fit is then used to estimate the initial parameters, \pi^[0], \sigma^{2}_(0),\lambda_0^[0],\Lambda^[0],\alpha}_k^[0], and also \pi_{k|i}, \gamma_i, \gamma_i\gamma_i^T, cov(gamma).
  ## data is a list of five objects (vectors) named as 'x', 'time', 'timeindex', 'curve', 'grid':
  ##      i) 'suppose we observe N independent subjects different (e.g, years), each consisting of n_i observations (time points for ith subject). 'x' is a vector of length \sum_i n_i with the first n_1 elements representing the observations of the first subject, followed by n_2 observations of the second subject, etc;
  ##      ii) 'time' is a (\sum_i n_i) vector of time points (usually time points are scaled to [0,1] for each subject),
  ##      iii) 'timeindex' is a (\sum_i n_i)vector of time indices from T possible from 'grid'. So each observation has a corresponding location (time index) within [0,1] uniquely specified time points.
  ##      iv)  'curve' is a (\sum_i n_i) vector of integers from 1,..., N, specifying the curve (subject number) for each observation in 'x'
  ##      v)   'grid' is a T vector of all unique time points (values within [0,1] interval) for all N years, needed for estimation of the B-spline coefficients in fda::eval.basis();
  ##      vi) if supplied, 'covariates' is an \eqn{N x r} matrix (or data frame) of scalar covariates (finite-dimensional covariates).
  ## h is a positive integer, parameter vector dimension in low-dimensionality representation of the curves (spline coefficients), h should be <= K, number of clusters. (default: h=2)
  ## q is number of splines (default: q=6); 
  ## K is a number of clusters. (default: K=3);
  ## random - TRUE/FALSE,  if TRUE the cluster belongings is given by uniform distribution, otherwise k-means is used to initialize cluster belongings;
  ## B is an N x p matrix of spline coefficients,  the spline approximation of the yearly curves based on p number of splines. Usually the matrix is supplied. If you supply it, check so that it has correct dimension, correct number of splines. If B=NULL (default), the coefficients are estimated using fda:: create.bspline.basis;
  ## svd is TRUE/FALSE, whether SVD decomposition should be used for the matrix of spline coefficients;
  ## use.covariates is TRUE/FALSE, whether covariates should be used when modelling;
  ## lambda is a positive real number, smoothing parameter value to be used when estimating B-spline coefficients
  
  S <- FullS <- NULL ## initialize matrix of basis functions evaluated at supplied points/grid
  ## get spline basis matrix...
  FullS <- fda:: eval.basis(data$grid,fda::create.bspline.basis(breaks=seq(0,1,len=q-2),norder=4))
  FullS.bmat <- FullS
  if(svd) {
    svd.FullS <- svd(FullS)
    d.inv <- rep(0,q)  # initial vector of inverse singular values
    good <- svd.FullS$d >= max(svd.FullS$d)*.Machine$double.eps^.5
    d.inv[good] <- 1/svd.FullS$d[good]
    DVt <- svd.FullS$d*t(svd.FullS$v)
    inv.DVt <- t(d.inv*t(svd.FullS$v))
    ## DV.t <- diag(svd.FullS$d) %*% t(svd.FullS$v)
    ## inv.DV.t<- solve(DV.t)
    FullS <- svd.FullS$u
  }
  dimFullS<-dim(FullS) ## T x q
  N <- length(table(data$curve)) ## number of subjects/years
  
  ## Here, if covariates are introduced. Need to augment the FullS matrix so that it 
  ## contains the diagonal matrix of ones to account for the covariates. 
  ## Also the covariates should be added to the 'x' vector (vector of curve observations)...
  
  ## first, check if covariates should be used in the model...
  moc <- FALSE ## signaling if model is with covariates (TRUE) or without (FALSE)
  if (!use.covariates) 
    covariates <- NULL
  else if (is.null(data[["covariates"]])) { ## then check if covariates are supplied...
    covariates <- NULL
    warning("no covariates supplied, clustering model without covariates is used instead")
  } else{ 
    covariates <- data$covariates
    moc <- TRUE
    if (!is.null(index.cov)) { ## subset of covariates should be used
      ind <- index.cov
      if (length(ind) > dim(covariates)[2]) 
        stop("Number of selected covariates is large than the number of covariates supplied")
      covariates <- as.matrix(covariates[,ind]) 
      rownames(covariates) <- rownames(covariates)
      colnames(covariates) <- colnames(data$covariates)[ind]
    }
  }
  
  ## if (is.null(data[["covariates"]])) 
  ##       covariates <- NULL
  ##   else covariates <- data$covariates
  
  if (!is.null(data$years.names))
    years.names <- data$years.names
  else if(!is.null(data$covariates))
    years.names <- dimnames(data$covariates)[[1]] ## years names are needed for trendplot() function (given as dimnames for 'varve' data, returns NULL for other data sets)
  else years.names <- NULL
  
  if(!is.null(covariates)){
    dimcov<-dim(covariates)
    if(is.null(dimcov)) {covariates<-matrix(covariates,ncol=1); dimcov<-dim(covariates)}
    if (dimcov[1]!=N) stop("Number of rows in the supplied matrix of covariates should be equal to N, number of curves/subjects")
    covariates <- as.matrix(covariates)
    
    if (stand.cov) ## covariates should be standardized by subtracting mean and dividing by std
      covariates <- scale(covariates)
    
    zero.vector.long<-matrix(0,nrow=dimFullS[1],ncol=dimcov[2])
    # zero.vector.short<-rep(0,len=dimFullS[2]*dimcov[2])
    # zero.vector.short<-matrix(zero.vector.short,nrow=dimcov[2])
    zero.vector.short<-matrix(0,nrow=dimcov[2],ncol=dimFullS[2])
    diagmatr<-diag(1,dimcov[2])
    newFullS<-cbind(FullS,zero.vector.long)
    # small.matr<-cbind(zero.vector.short,diagmatr)
    # newFullS<-rbind(newFullS,small.matr)# tag.FullS
    newFullS<-rbind(newFullS,cbind(zero.vector.short,diagmatr))# matrix S: covariance.FullS
    data.and.covariates<-NULL
    new.timeindex<-NULL
    new.curve<-NULL
    temp<-data$timeindex[data$curve==1]
    add.values<-temp[length(temp)]+1
    add.values<-add.values:c(add.values+dimcov[2]-1)
    for(i in 1:N) {
      data.and.covariates<-c(data.and.covariates,c(data$x[data$curve==i],as.numeric(covariates[i,])))  # covariance.x
      new.timeindex<-c(new.timeindex,c(data$timeindex[data$curve==i],add.values))# covariance.timeindex
      temp<-data$curve[data$curve==i]
      new.curve<-c(new.curve,c(temp,rep(temp[1],len=dimcov[2])))#covariance.curve
    }
    Snew <- newFullS[new.timeindex,]
  } else{ ## no covariates included...
    data.and.covariates<-new.timeindex<-new.curve<-Snew<-NULL
    dimcov<-c(0,0)
  }
  
  S <- FullS[data$timeindex,] ## without covariates
  StS.inv <- solve(t(S) %*% S)
  ## qrR <- qr.R(qr(S))
  ## StS.inv <-  solve(t(qrR)%*% qrR)
  
  
  ## get initial estimates of B-spline coefficients if they are not provided...
  if(is.null(B)){
    B <- matrix(0,N,q)
    ## lambda<-0.000140625
    bsp.length <- q-2
    bsp1 <- fda::create.bspline.basis(norder=4, breaks=seq(0,1,len=bsp.length))
    fdPar2 <- fda::fdPar(bsp1, Lfdobj=2, lambda=lambda)
    ## for (i in 1:N) B[i,] <- t(fda::smooth.basis(data$time[data$curve==i], data$x[data$curve==i], fdPar2)$fd$coefs)
    get.beta0 <- function(ind,data,fdPar2){
      fda::smooth.basis(data$time[data$curve==ind], data$x[data$curve==ind], fdPar2)$fd$coefs
    }
    B <- t(sapply(1:N,get.beta0,data,fdPar2)) ## consider parallelizing here...!!!
    ## Spara<<-B
  }
  
  ## Use k-means to get initial cluster memberships from B...
  if(K > 1){
    if(!random){
      if(use.covariates) klass <- kmeans(x=cbind(B,covariates), centers=K, iter.max=100, nstart=500)$cluster ## , algorithm="Lloyd"
      else klass <- kmeans(x=B, centers=K, iter.max=100, nstart=500)$cluster ## kmeans() from stats r package
    } else 
      klass <- sample(K, N, replace =TRUE) ## cluster belongings via uniform distribution
  }  else klass <- rep(1, N)
  
  ## initialize estimates for the posterior probs of cluster membership...
  piigivej <- matrix(0, N, K) 
  piigivej[col(piigivej) == klass] <- 1 ## N x K matrix of 0s and 1s
  ## calculate coefficients for cluster means...
  classmean <- matrix(0,K,q)
  for (k in 1:K)
    classmean[k,] <- apply(B[klass==k,],2,mean)
  ## initialize lambda.zero, Lambda and alpha as defined in the paper...
  lambda.zero <- apply(classmean, 2, mean) ## q-vector of cluster mean spline coefficients
  lambda.zero <- as.vector(lambda.zero)
  Lambda <- as.matrix(svd(scale(classmean, scale = FALSE))$v[,1:h]) ##  q x h transformation matrix for coefs dimensionallity reduction 
  alpha <- scale(classmean, scale = FALSE)%*%Lambda ## K x h matrix, where ith row is alpha_k, an h-vector of low-dimentional representation of coeffs (curves) for kth cluster
  pii <- apply(piigivej,2,mean) ## a K-vector of initial posterior probs of cluster membership 
  ## calculate estimates for gamma and gamma.gammaT...
  gamma <- t(t(B) - lambda.zero - (Lambda %*% t(alpha[klass, ]))) ## N x q matrix, each row is a q-vector of gammai
  gprodK <- array(0,c(N,K,q,q)) ## pki*gamma.gammaT is a q x q matrix for each observ (out of N) and each cluster, where gamma is a q-vector
  for(i in 1:N){
    for(j in 1:K) gprodK[i,j,,] <- pii[j]*gamma[i,] %*% t(gamma[i,]) ## pki*gamma.gammaT
  }
  gamma <- array(gamma[, rep(1:sum(q), rep(K, sum(q)))], c(N, K, sum(q))) ## N x K x q array , initialization of gammai, a q-vector for each observations, cluster specific. copies of gamma from above as initialization
  
  ## below 'tag' stands for 'covariates'...
  gcovK <- array(0, c(N,K,q,q)) ## initial array of covariances of gamma set as array of 0.1's, cluster specific # 0s
  if(!is.null(covariates)){
    tag.classmean <- matrix(0,K,dimcov[2]) ## cluster means for covariates
    for (k in 1:K)
      tag.classmean[k,] <-  apply(as.matrix(covariates[klass==k,]),2,mean)
    tag.lambda.zero <- c(lambda.zero,apply(tag.classmean, 2, mean)) ## lambda.zero and covariates together
    tag.lambda.zero <- as.vector(tag.lambda.zero)
    tag.Lambda <- as.matrix(svd(scale(cbind(classmean,tag.classmean), scale =FALSE))$v[, 1:h]) ## with covariates
    tag.alpha <- scale(cbind(classmean,tag.classmean), scale = FALSE)%*% tag.Lambda ## with covariates
    tag.gamma <- t(t(cbind(B,covariates)) - tag.lambda.zero - (tag.Lambda %*% t(tag.alpha[klass,  ]))) ## gamma with covariates
    tag.gprodK <- array(0,c(N,K,q+dimcov[2],q+dimcov[2])) ## array of gamma.gammaT*pi with covariates
    for(i in 1:N){
      for(j in 1:K) tag.gprodK[i,j,,] <- tag.gamma[i,] %*% t(tag.gamma[i,])*pii[j] ## gamma.gammaT*pi
    }
    tag.gamma <- array(tag.gamma[, rep(1:sum(q+dimcov[2]), rep(K, sum(q+dimcov[2])))], c(N, K, sum(q+dimcov[2]))) ## array of gamma with coavriates
    tag.gcovK<- array(0, c(N,K,q+dimcov[2],q+dimcov[2])) ## initial array of covariances of gamma, set as .1's # 0s
  } else{
    tag.gamma<-tag.gprodK<-tag.gcovK<-tag.piigivej<-tag.S<-newFullS<-NULL
  }
  
  if(svd) {  
    lambda.zero <- DVt %*% as.matrix(lambda.zero,ncol=1)  
    Lambda <- DVt %*% Lambda
  }
  initfit.mod <- list(
    design = list(FullS = FullS, S = S, StS.inv=StS.inv, tag.S=Snew,  tag.FullS = newFullS, FullS.bmat = FullS.bmat), 
    parameters = list(lambda.zero = lambda.zero, Lambda = Lambda, alpha = alpha), 
    vars =  list(gamma=gamma, tag.gamma=tag.gamma, gprodK = gprodK, tag.gprodK=tag.gprodK, gcovK = gcovK,  
                 tag.gcovK=tag.gcovK, piigivej = piigivej, tag.piigivej=piigivej),
    data.aug =  list(x=data$x, time=data$time, timeindex=data$timeindex, curve=data$curve, tag.x=data.and.covariates,
                     tag.curve=new.curve, tag.timeindex=new.timeindex, covariates=covariates, covariates.original=data$covariates, grid = data$grid, B=B), 
    initials = list(q=q,N=N, Q=q+dimcov[2], random=random,h=h, K=K,r =dimcov[2],moc=moc,years.names=years.names) 
  )
  ## above 'parameters' is a list of estimated objects to be extended in M-step, whereas
  ## 'vars' is a list of estimated objects coming from E-step 
  initfit.mod
} ## end initial.mocca


################################
## M step of the EM algorithm... 
################################

lamda0.mocca <- function(N,K,q,x,S,StS.inv,curve,Lambda,alpha,tag.gamma,tag.piigivej){ 
  ## function to calculate current estimate of lambda0 withing M-step...
  i <- NULL ## binding the variable locally to the function, so the 'R CMD check' has nothing to complain about 
  y <- foreach::foreach(i = 1:N) %do% {
    xi <- x[curve==i]
    ni<-length(xi)
    SLag <- S[curve==i,]%*%(Lambda%*%t(alpha)+t(tag.gamma[i,,1:q]))
    summa<-matrix(0,nrow=ni,ncol=1)
    for(k in 1:K) summa<-summa+matrix(SLag[,k]*tag.piigivej[i,k])
    Si.dev<-matrix(t(S[curve==i,])%*%(matrix(xi)-summa))
    return(Si.dev)
  }
  totsum<-apply(simplify2array(y)[,1,],1,sum)
  StS.inv%*%totsum
} ## end lamda0.mocca



Mstep.mocca <- function(parameters, vars, data, initials, design, Mstep.maxit=10,Mstep.tol=1e-4,trace=TRUE) {
  ## This function implements the M step of the EM algorithm.
  ## 'parameters' is a list of objects from  initial.mocca();
  ## 'vars' is a list of objects from  initial.mocca();
  ## 'data' is a list of supplied data from  initial.mocca();
  ## 'initials' is a list of objects from initial.mocca();
  ## 'design' is a list of objects from  initial.mocca();
  ## 'Mstep.maxit' is a positive scalar which gives the maximum number of iterations for an inner loop of the parameter estimation;
  ## 'Mstep.tol' is the tolerance to use within iterative procedure to estimate model parameters;
  ##  'trace' is logical indicating whether to print the current values of sig² and sig²_covariates at each iteration.
  ## 1) for the model with functional data but without covariates, the estimates of six parameter objects are updated here: 
  ##    pi and Gamma explicitely, while lambda.0, alpha and Lambda via sequential iterative procedure, and  finally
  ##    the estimate of sigma^2 is updated.
  ## 2) for the model with functional data and covariates, additionally, estimates of the covariates means for each cluster, v1,...,vG, should be updated; as well as the common variance of the covariances. The estimation of the covariance matrix of the covariances (to account for random deviation within cluster) is within estimation of Gamma which is augmented with this covariance matrix. in overall, there are eight objects to be updated.
  ## the routine returns a list of the estimated parameters.
  
  ## first, set up variable names...
  alpha <- parameters$alpha
  Lambda <- parameters$Lambda  ## q x h  matrix
  gcovK<-vars$gcovK ## cov of gamma cluster specific
  gprodK<-vars$gprodK ## pki*gamma.gammaT
  curve <- data$curve
  x <- data$x
  S <- design$S
  StS.inv <- design$StS.inv
  q <- initials$q ## spline dimension, colnum of S
  Q <- initials$Q ## q plus number of covariates
  n <- length(curve) ## total number of observations (N*ni)
  gamma<-vars$gamma[,,1:q]
  N <- initials$N ## number of subjects/years
  K <- dim(alpha)[1] ## number of clusters 
  h <- dim(alpha)[2] ## dimension of alpha, h<=K,  for spline dimensionallity reduction 
  
  
  if(initials$moc){ ## if the model with scalar covariates... 
    ## additional setting up variable names...
    tag.gcovK <- vars$tag.gcovK;
    tag.gprodK <- vars$tag.gprodK
    covariates <- data$covariates
    tag.piigivej <- vars$tag.piigivej; 
    tag.S <- design$tag.S
    ## gamma<-vars$tag.gamma[,,1:q]
    tag.gamma <- vars$tag.gamma;
    # N <- dim(tag.gamma)[1]; 
    
    ## get updated estimate of the parameter pi=(pi_1,...,pi_K), probabilities of cluster membership, via
    ## explicit formula (15) in the paper... 
    parameters$tag.pi <- apply(tag.piigivej, 2, mean) 
    
    ## get updated estimate of (rank p ??) estimate of Gamma...
    parameters$GammaK <- array(0,c(K,q,q))
    parameters$tag.GammaK <- array(0,c(K,Q,Q))   
    ## indm is the Identity matrix for the dimension of the splines repeated num of obs times.
    tag.indm <- matrix(rep(c(rep(c(1, rep(0, Q - 1)), N), 0), Q)[1:(N*Q^2)], N * Q, Q)
    tag.piigivej.q <- tag.piigivej[sort(rep(1:N,Q)),]
    for(k in 1:K){
      tag.indm.q <- tag.indm
      for(j in 1:Q) tag.indm.q[,j]<-tag.indm.q[,j]*tag.piigivej.q[,k]
      ## tag.gprod <- NULL
      ## for(j in 1:N) tag.gprod<-rbind(tag.gprod,tag.gprodK[j,k,,])
      # tag.gprod <- foreach::foreach(j =1:N, .combine = rbind) %dopar% tag.gprodK[j,k,,]
      tag.gprod_t <- tag.gprodK[1,k,,]
      if(N > 1){
        for(j in 2:N){
          tag.gprod_t <- rbind(tag.gprod_t, tag.gprodK[j,k,,])
        }
      }
      tag.gprod <- tag.gprod_t
      ind <- is.na(tag.gprod)
      if (sum(ind)!=0) 
        tag.gprod[ind] <- 0
      ind <- is.na(tag.indm.q)
      if (sum(ind)!=0) tag.indm.q[ind] <- 0
      if (is.na(parameters$tag.pi[k])) parameters$tag.pi[k] <- 1e-5
      tag.gsvd <- svd(t(tag.gprod) %*% tag.indm.q/(N*parameters$tag.pi[k]))
      good<- tag.gsvd$d >= max(tag.gsvd$d)*.Machine$double.eps^.5
      tag.gsvd$d[!good]<-0
      ## parameters$tag.GammaK[k,,] <- tag.gsvd$u %*% diag(tag.gsvd$d) %*% t(tag.gsvd$u) ??? CHECK expression!!!!!!!!!!!!!!!!!!!!!!!
      parameters$tag.GammaK[k,,] <- tag.gsvd$u %*% (tag.gsvd$d*t(tag.gsvd$u)) ## explicit estimate of Delta.k(augmented Gamma.k)
    }
    
    ## get estimates for the covariates, vk, (modelled as xi=vk+deltai +ei), vk is an estimated covariate cluster mean... 
    tag.estimates<-NULL
    for(k in 1:K) {
      ans<-c(apply(vars$tag.piigivej[,k]*(as.matrix(covariates[,])-vars$tag.gamma[,k,(q+1):Q]),2,sum)/c(N*parameters$tag.pi)[k]) # vk      
      tag.estimates<-rbind(tag.estimates,ans)
    }
    parameters$tag.estimates<-tag.estimates ## explicit estimate of vk
    
    ## next is a loop that iteratively calculates the updated estimates of lambda.0, alpha and then Lambda and stops untill convergence
    ## of estimates...
    # max.lambda.zero.old <- 2;  max.lambda.zero <- 1  
    # max.alpha.old <- 2;  max.alpha <- 1
    # max.Lambda.old <- 2;  max.Lambda <- 1
    para.old <- c(2,2,2)
    para <- c(1,1,1)
    
    loopit <-1
    while((max(abs(para - c(para.old))) > Mstep.tol*max(abs(para + c(para.old)))/2) & (loopit < Mstep.maxit)){
      ## while((abs(sigma.old[1] - sigma[1])/sigma[1] > Mstep.tol) & (loopit < Mstep.maxit)){
      ## change to check for lambdazero, alpha, Lambda convergence...??
      ## if (max(abs(beta - c(betaold))) > control$steptol.fit*max(abs(beta + c(betaold)))/2) 
      para.old <- para
      ## get estimate of lambda0...
      lambda.zero <- lamda0.mocca(N,K,q,x,S,StS.inv,curve,Lambda,alpha,tag.gamma,tag.piigivej)  ## q-vector    
      ## calculate estimate of alpha...
      xcent <- x - S %*% lambda.zero
      S.Lam <- S %*% Lambda
      for(i in 1:K) {
        S.Lam.pi <- S.Lam * tag.piigivej[curve, i]
        LtStSL <- t(S.Lam.pi)%*% S.Lam
        if(sum(tag.piigivej[, i]) > 1e-4){
          tgamma.i <- t(tag.gamma[curve,i,1:q])
          diagonalen <- foreach(j = 1:N, .combine = c) %do% diag(S[curve==j,]%*%tgamma.i[,curve==j])
          alpha[i,  ] <- solve(LtStSL) %*% t(S.Lam.pi) %*% (xcent - diagonalen) ## K x h matrix
        } else print("Warning: empty cluster")
      }
      ## calculate Lambda given alpha. This is done by iterating
      ## through each column of Lambda holding the other columns fixed...
      for(m in 1:h) {
        pi.alphasq <- apply(t(tag.piigivej) * (alpha^2)[, m], 2,sum)[curve]
        pi.alpha <- apply(t(tag.piigivej) * alpha[, m], 2, sum)[curve]
        tS.Lambda <- t(S.Lam)
        if(h != 1) {
          temp <- NULL
          for(i in 1:K) {
            temp <- cbind(temp, as.vector(rep(1, h - 1) %*% 
                                            matrix((tS.Lambda * alpha[i,  ])[ - m,  ], h - 1,dim(S)[1])) * alpha[i, m])
          }
          otherLambda <- (temp * tag.piigivej[curve,  ])%*%rep(1, K)
        } else otherLambda <- 0
        gamma.pi.alpha <- apply(tag.gamma[,,1:q] * as.vector(tag.piigivej) * rep(alpha[, m], rep(N, K)), c(1, 3), sum)[curve,  ]
        Lambda[, m] <- solve(t(S * pi.alphasq) %*% S) %*% t(S) %*% (xcent * pi.alpha - otherLambda - (S *gamma.pi.alpha) %*% rep(1, sum(q))) ## q x h matrix
      }
      ## calculate max values of each estimated para's...
      para <- c(max(lambda.zero), max(alpha),max(Lambda))
      
      loopit <- loopit + 1
      if (trace){ ## printing current values of max estimates...
        para0<- round(para, 3)
        too <- setNames(c("lambda0","alpha","Lambda"), para0)
        paste(too, names(too), sep = "=", collapse = "; ")
        print(too)
      }
    } ## end while loop
  } else{ ## model without scalar covariates...
    ## gamma<-vars$gamma[,,1:q]
    piigivej <- vars$piigivej
    ## get updated estimate of the probabilities of cluster membership...
    parameters$pi <- apply(piigivej, 2, mean);
    
    ## get updated estimate of (rank p ??) estimate of Gamma...
    parameters$GammaK <- array(0,c(K,q,q))
    ## indm is the Identity matrix for the dimension of the splines repeted num of obs times...
    indm <- matrix(rep(c(rep(c(1, rep(0, q - 1)), N), 0), q)[1:(N*q^2)], N * q, q)
    piigivej.q<-piigivej[sort(rep(1:N,q)),]
    for(k in 1:K){
      indm.q <- indm
      for(j in 1:q) indm.q[,j] <- indm.q[,j]*piigivej.q[,k]
      gprod <- foreach::foreach(j =1:N, .combine = rbind) %do% gprodK[j,k,,]
      gsvd <- svd(t(gprod) %*% indm.q/(N*parameters$pi[k]))
      good <- gsvd$d >= max(gsvd$d)*.Machine$double.eps^.5
      gsvd$d[!good]<-0
      ## parameters$GammaK[k,,] <- gsvd$u %*% diag(gsvd$d) %*% t(gsvd$u) 
      parameters$GammaK[k,,] <- gsvd$u %*% (gsvd$d*t(gsvd$u))  ## equation (18), explicit estimate of Gamma.k. CHECK HERE, uduT??? or udVt? !!!!!
    }
    
    ## loop that iteratively calculates the updated estimates of lambda.0, alpha and then Lambda and stops untill convergence
    ## of estimates...
    para.old <- c(2,2,2)
    para <- c(1,1,1)
    loopit<-1
    while((max(abs(para - c(para.old))) > Mstep.tol*max(abs(para + c(para.old)))/2) & (loopit < Mstep.maxit)){
      ## while((abs(sigma.old[1] - sigma[1])/(.1 + abs(sigma[1])) > Mstep.tol) & (loopit < Mstep.maxit)){
      ## while((abs(sigma.old[1] - sigma[1])/sigma[1] > Mstep.tol) & (loopit < Mstep.maxit)){
      para.old <- para
      ## get estimate of lambda0...
      lambda.zero <- lamda0.mocca(N,K,q,x,S,StS.inv,curve,Lambda,alpha,gamma,piigivej)
      
      ## calculate alpha...
      xcent <- x - S %*% lambda.zero
      S.Lam <- S %*% Lambda
      for(i in 1:K) {
        S.Lam.pi <- S.Lam * piigivej[curve, i]
        LtStSL <- t(S.Lam.pi)%*% S.Lam
        if(sum(piigivej[, i]) > 1e-4){
          tgamma.i <- t(gamma[curve,i,1:q])
          diagonalen <- foreach::foreach(j = 1:N, .combine = c) %do% diag(S[curve==j,]%*%tgamma.i[,curve==j])
          alpha[i,  ] <- solve(LtStSL) %*% t(S.Lam.pi) %*% (xcent - diagonalen)
        } else print("Warning: empty cluster")
      }
      ## calculate Lambda given alpha. This is done by iterating
      ## through each column of Lambda holding the other columns fixed...
      for(m in 1:h) {
        pi.alphasq <- apply(t(piigivej) * (alpha^2)[, m], 2,sum)[curve]
        pi.alpha <- apply(t(piigivej) * alpha[, m], 2, sum)[curve]
        tS.Lambda <- t(S.Lam)
        if(h!= 1){
          temp <- NULL
          for(i in 1:K) {
            temp <- cbind(temp, as.vector(rep(1, h - 1) %*% 
                                            matrix((tS.Lambda * alpha[i,  ])[ - m,  ], h - 1,dim(S)[1])) * alpha[i, m])
          }
          otherLambda <- (temp * piigivej[curve,  ])%*%rep(1, K)
        } else otherLambda <- 0
        gamma.pi.alpha <- apply(gamma[,,1:q] * as.vector(piigivej) * rep(alpha[, m], rep(N, K)), c(1, 3), sum)[curve,  ]
        Lambda[, m] <- solve(t(S * pi.alphasq) %*% S) %*% t(S) %*% (xcent * pi.alpha - otherLambda - (S *gamma.pi.alpha) %*% rep(1, sum(q)))
      }
      ## calculate max values of each estimated para's...
      para <- c(max(lambda.zero), max(alpha),max(Lambda))
      loopit <- loopit + 1
    } ## end while loop
  } ## end else, model without covariates
  ##  stopCluster(cl)
  parameters$lambda.zero <- as.vector(lambda.zero);
  parameters$alpha <- alpha;
  parameters$Lambda <- Lambda;
  parameters
} ## end Mstep.mocca




################################
## E step of the EM algorithm... 
################################

Estep.mocca <- function(parameters, variables, data, design, initials) { ## variables=vars
  ## This function performs the E step of the EM algorithm by
  ## calculating the expected values of gamma and gamma %*% t(gamma) and piigivej.
  ## given the current parameter estimates.
  ## It does it in two ways depending on whether there are any covariates or not. 
  ## 'parameters' is a list of several objects from M-step, current parameter values to be estimated
  ## 'variables' is a list of eight objects, current variables values
  ## 'data' is a list of ten objects, original list of data plus additional objects from initial step
  ## 'design' is a list of objects, design matrices and their variants
  ## Above 'parameters' is a list of estimated objects coming from M-step, whereas  
  ## 'vars' is a list of estimated objects coming from the previous E-step. 
  
  ## below 'tag' stands for 'covariates'...
  if(initials$moc){ ## if there are covariates included...
    ## the variance of gamma is allowed to vary between clusters and covariance matrix is included.
    N <- initials$N #dim(variables$tag.gamma)[1]
    K <- initials$K # dim(variables$tag.gamma)[2]
    Q <- initials$Q # dim(variables$tag.gamma)[3] ## q+r (spline dimension plus number of covariates)
    r <- initials$r ## number of scalar covariates
    tag.GammaK <- parameters$tag.GammaK # dim K x Q x Q, 
    tag.estimates<- parameters$tag.estimates
    # variables$tag.gammaK  dim N x K x K x Q
    # variables$tag.gprodK   dim N x K x q x Q
    # variables$tag.gcovK    dim N x K x q x Q
    
    ## Start with calculating the expected values of the parameters...
    Lambda.alpha <- as.vector(parameters$lambda.zero) + parameters$Lambda %*% t(parameters$alpha)
    tag.Lambda.alpha <- rbind(Lambda.alpha,t(tag.estimates))
    variables$tag.gprodK <- array(0,c(N,K,Q,Q))
    variables$tag.gcovK <- array(0,c(N,K,Q,Q))
    tag.S <- design$tag.S
    curve <- data$curve; tag.curve <- data$tag.curve
    ## Calculate expected value of gamma...
    j <- NULL ## binding the variable locally to the function, so the 'R CMD check' has nothing to complain about 
    temp <- foreach::foreach(j = 1:N, .packages = "mvtnorm") %do% { #1
      ##  for(j in 1:2) { #1
      tag.Sj <- tag.S[tag.curve == j,]
      tag.nj <- sum(tag.curve == j)
      len.sigma<-length(parameters$sigma)
      tag.varm <- rep(parameters$sigma[1],tag.nj)
      ind <- (tag.nj-r+1):tag.nj
      tag.varm[ind] <- parameters$sigma[2]
      ## for(i in (tag.nj-r+1):tag.nj)  tag.varm[i] <- parameters$sigma[2]
      d.tag.invar <- rep(0,tag.nj)
      good <- tag.varm >= max(tag.varm)*.Machine$double.eps^.5
      d.tag.invar[good] <- 1/tag.varm[good]
      ## tag.invvar <- diag(d.tag.invar)
      Ccov.GammaK <- tag.GammaK
      ## for(k in 1:K)
      ##    Ccov.GammaK0[k,,] <- tag.GammaK[k,,]-tag.GammaK[k,,]%*%t(tag.Sj)%*%tag.invvar%*%
      ##         solve(diag(tag.nj)+tag.Sj%*%tag.GammaK[k,,]%*%t(tag.Sj)%*%tag.invvar)%*%tag.Sj%*%tag.GammaK[k,,] ## p.21
      SGamma <- list()
      for(k in 1:K){
        SGamma[[k]] <- tag.Sj%*%tag.GammaK[k,,]
        terminv <- solve(diag(tag.nj)+SGamma[[k]]%*%t(d.tag.invar*tag.Sj))
        Ccov.GammaK[k,,] <- tag.GammaK[k,,]- t(SGamma[[k]])%*%(d.tag.invar*terminv)%*%SGamma[[k]]
      }
      
      tag.centx <- data$tag.x[data$tag.curve == j] - tag.Sj %*% tag.Lambda.alpha
      for(k in 1:K)
        variables$tag.gamma[j,k, ] <- t(Ccov.GammaK[k,,]%*%t(tag.Sj)%*%(d.tag.invar*tag.centx[,k]))
      
      ## Calculate pi i given j...
      tag.d <- NULL
      for(k in 1:K){ #3
        ## tag.covx <- tag.Sj%*%tag.GammaK[k,,] %*% t(tag.Sj) + diag(tag.varm)
        tag.covx <- SGamma[[k]]  %*% t(tag.Sj) + diag(tag.varm) 
        tag.d <- c(tag.d,dmvnorm(tag.centx[,k],rep(0,length(tag.centx[,k])),tag.covx) * parameters$tag.pi[k])
      } #3
      variables$tag.piigivej[j, ] <- tag.d/sum(tag.d)
      test<-tag.d/sum(tag.d)
      ## Calculate expected value of gamma %*% t(gamma)...
      for(k in 1:K) { #5
        variables$tag.gprodK[j,k,,] <- t(matrix(variables$tag.gamma[j,,],K,Q))%*%(matrix(variables$tag.gamma[j,,],K,Q)*
                                                                                    variables$tag.piigivej[j,]) + Ccov.GammaK[k,,]
        variables$tag.gcovK[j,k,,] <- Ccov.GammaK[k,,]
      } #5
      temp <- list(tag.gamma=variables$tag.gamma[j,,],tag.piigivej=variables$tag.piigivej[j,],
                   tag.gprodK=variables$tag.gprodK[j,,,],tag.gcovK=variables$tag.gcovK[j,,,])
    } #1 (end foreach:: )
    test<-list()
    test$tag.gamma<-simplify2array(lapply(temp,'[[',1))
    test$tag.gamma<-aperm(test$tag.gamma, c(3,1,2))
    test$tag.piigivej<-simplify2array(lapply(temp,'[[',2))
    test$tag.piigivej<-t(test$tag.piigivej)
    ind <- is.na(test$tag.piigivej)
    if (sum(ind)!=0) 
      test$tag.piigivej[ind] <- 0
    test$tag.gprodK<-simplify2array(lapply(temp,'[[',3))
    test$tag.gprodK<-aperm(test$tag.gprodK, c(4,1,2,3))
    test$tag.gcovK<-simplify2array(lapply(temp,'[[',4))
    test$tag.gcovK<-aperm(test$tag.gcovK, c(4,1,2,3))
    test
  } ####################################
  else{ ## no covariates included ...
    N <- dim(variables$gamma)[1]
    K <- dim(variables$gamma)[2]
    q <- dim(variables$gamma)[3] 
    GammaK <- parameters$GammaK # dim K x q x q
    ## variables$gammaK  dim N x K x K x q
    ## variables$gprodK   dim N x K x q x q
    ## variables$gcovK    dim N x K x q x q
    
    ## Start with calculating the expected values of the parameters...
    Lambda.alpha <- as.vector(parameters$lambda.zero) + parameters$Lambda %*% t(parameters$alpha)
    variables$gprodK <- array(0,c(N,K,q,q))
    variables$gcovK <- array(0,c(N,K,q,q))
    S <- design$S
    curve <- data$curve
    ## Calculate expected value of gamma...
    temp <- foreach::foreach(j = 1:N, .packages = "mvtnorm") %do% { #1
      ##  for(j in 1:2) { #1
      Sj <- S[curve == j,]
      nj <- sum(curve == j)
      varm <-  rep(parameters$sigma[1],nj) ## diag(parameters$sigma[1],nj)
      ## invvar <- diag(1/diag(varm))
      d.invvar <- rep(0,nj)
      good <- varm >= max(varm)*.Machine$double.eps^.5
      d.invvar[good] <- 1/varm[good]
      CGammaK <- GammaK
      SGamma <- list()
      for(k in 1:K){
        ## CGammaK[k,,] <- GammaK[k,,]-GammaK[k,,]%*%t(Sj)%*%invvar%*%
        ##       solve(diag(nj)+Sj%*%GammaK[k,,]%*%t(Sj)%*%invvar)%*%Sj%*%GammaK[k,,]
        SGamma[[k]] <- Sj%*%GammaK[k,,]
        terminv <- solve(diag(nj)+SGamma[[k]]%*%t(d.invvar*Sj))
        CGammaK[k,,] <- GammaK[k,,]- t(SGamma[[k]])%*%(d.invvar*terminv)%*%SGamma[[k]]
      }
      centx <- data$x[data$curve == j] - Sj %*% Lambda.alpha
      for(k in 1:K)
        variables$gamma[j,k, ] <- t(CGammaK[k,,]%*%t(Sj)%*%(d.invvar*centx[,k]))
      ## Calculate pi i given j...
      d <- NULL
      for(k in 1:K){ #3
        covx <- SGamma[[k]] %*% t(Sj) + diag(varm)
        d <- c(d,dmvnorm(centx[,k],rep(0,length(centx[,k])),covx) * parameters$pi[k])
      } #3
      variables$piigivej[j,  ] <- d/sum(d)
      test<-d/sum(d)
      ## Calculate expected value of gamma %*% t(gamma).
      for(k in 1:K) { #5
        variables$gprodK[j,k,,] <- t(matrix(variables$gamma[j,,],K,q))%*%(matrix(variables$gamma[j,,],K,q)*
                                                                            variables$piigivej[j,]) + CGammaK[k,,]
        variables$gcovK[j,k,,] <- CGammaK[k,,]
      } #5
      temp <- list(gamma=variables$gamma[j,,],piigivej=variables$piigivej[j,],
                   gprodK=variables$gprodK[j,,,],gcovK=variables$gcovK[j,,,])
    } #1
    test<-list()
    test$gamma<-simplify2array(lapply(temp,'[[',1))
    test$gamma<-aperm(test$gamma, c(3,1,2))
    test$piigivej<-simplify2array(lapply(temp,'[[',2))
    test$piigivej<-t(test$piigivej)
    test$gprodK<-simplify2array(lapply(temp,'[[',3))
    test$gprodK<-aperm(test$gprodK, c(4,1,2,3))
    test$gcovK<-simplify2array(lapply(temp,'[[',4))
    test$gcovK<-aperm(test$gcovK, c(4,1,2,3))
    test
  } ## end 'else' no covariates
} ## end Estep.mocca


#########################################################
## function to calculate log likelihood of the model...
#########################################################

loglik.EMmocca <- function(data,parameters,design){
  N<-length(unique(data$curve))
  Lambda.alpha <- as.vector(parameters$lambda.zero) + parameters$Lambda %*% t(parameters$alpha)
  
  ## First, if we have a single Gamma-matrix...
  if(is.null(parameters$GammaK)) {
    K<-length(parameters$pi)
    parameters$GammaK<-array(0,c(K,length(parameters$lambda.zero),length(parameters$lambda.zero)))
    for(i in 1:K) parameters$GammaK[i,,]<-parameters$Gamma
  }
  
  ## Next, if we have separate Gamma-matrices and covariates...
  j <- NULL ## binding the variable locally to the function, so the 'R CMD check' has nothing to complain about 
  if(!is.null(design$tag.S)) {
    ## print('hello')
    K <- length(parameters$tag.pi)
    tag.S<-design$tag.S
    antcov<-dim(data$covariates)[2]
    tag.GammaK <- parameters$tag.GammaK
    tag.Lambda.alpha <- rbind(Lambda.alpha,t(parameters$tag.estimates))
    summa <- foreach::foreach(j = 1:N, .combine = cbind, .packages = "mvtnorm") %do% {
      tag.Sj <- tag.S[data$tag.curve == j,  ]
      nj <- sum(data$tag.curve == j)
      varm <- diag(parameters$sigma[1], nj)
      nj<-nj+1
      for(sista in 1:antcov){
        nj<-nj-1
        varm[nj,nj]<-parameters$sigma[2] ## matrix R in the paper
      }
      centx <- data$tag.x[data$tag.curve == j] - tag.Sj %*% tag.Lambda.alpha
      d <- NULL
      for(k in 1:K){
        covx <- tag.Sj %*% tag.GammaK[k,,] %*% t(tag.Sj) + varm ## covariance matrix of the model with covariates (ui|zi)
        d <- c(d,dmvnorm(centx[,k],rep(0,length(centx[,k])),covx) * parameters$tag.pi[k])
      }
      d
    } ## end foreach
  } else{
    K <- length(parameters$pi)
    S<-design$S
    GammaK <- parameters$GammaK
    ## j <- NULL ## binding the variable locally to the function, so the 'R CMD check' has nothing to complain about 
    summa <- foreach:: foreach(j = 1:N, .combine = cbind, .packages = "mvtnorm") %do% {
      Sj <- S[data$curve == j,  ]
      nj <- sum(data$curve == j)
      varm <- diag(parameters$sigma, nj)
      centx <- data$x[data$curve == j] - Sj %*% Lambda.alpha
      d <- NULL
      for(k in 1:K){
        covx <- Sj %*% GammaK[k,,] %*% t(Sj) + varm ## covariance matrix of the model without covariates
        d <- c(d,dmvnorm(centx[,k],rep(0,length(centx[,k])),covx) * parameters$pi[k])
      }
      d
    } ## end foreach
  } ## end else
  summa<-apply(summa,2,sum)
  ind <- summa ==0
  if (sum(ind)!=0)
    summa <- summa+.0001*min(summa[!ind])
  summa<- log(summa)  
  summa<-sum(summa)
  summa
} ## end loglik.EMmocca





##############################################
## function to obtain residual variance estimates......
##############################################

fdvariance.mocca <- function(parameters, vars, data, design){
  ## function to only obtain a residual variance estimate for the functional data, sig^2, needed in M-step... 
  lambda.zero <- parameters$lambda.zero
  Lambda <- parameters$Lambda 
  alpha <- parameters$alpha
  S <- design$S
  curve <- data$curve
  q <- dim(S)[2]
  x <- data$x
  
  if(!is.null(design$tag.S)) { ## model with covariates...
    tag.gamma<- vars$tag.gamma
    tag.estimates <- parameters$tag.estimates
    tag.S<-design$tag.S
    tag.piigivej <- vars$tag.piigivej
    CovMatr <- data$covariates ## matrix of covariates
    Q <- dim(tag.S)[2]
    K <- dim(tag.piigivej)[2] ## number of clusters
    N <- dim(tag.piigivej)[1]
    gcov <- vars$tag.gcovK ## array of covariances of gamma, q x q matrices. Dimension of the array is c(N,K,q+r,q+r) 
    num.of.cov <- dim(tag.estimates)[2]
  } else{ ## model without covariates...
    tag.gamma<- vars$gamma
    tag.piigivej <- vars$piigivej
    K <- dim(vars$piigivej)[2]
    N <- dim(vars$piigivej)[1]
    gcov <- vars$gcovK ## array of covariances of gamma, q x q matrices. Dimension of the array is c(N,K,q,q)
  }
  
  if(length(dim(gcov))<=3){ ## not sure if needed...?????????????
    covariance<-array(0,c(N,q,q))
    k<-0
    for(i in seq(q,q*N,by=q)){
      k<-k+1
      covariance[k,,]<-gcov[,(i-(q-1)):i]
    }
    contr <- foreach::foreach(i = 1:N, .combine = c) %do% {
      residual <- matrix(x[curve==i],ncol=K,nrow=length(x[curve==i]))-(S[curve==i,]%*%
                                                                         (matrix(lambda.zero,ncol=K,nrow=q)+Lambda%*%t(alpha[,])+t(gamma[i,,])))
      ## trace <- 0
      trace <- sum(diag(S[curve==i,]%*%covariance[i,,]%*%t(S[curve==i,])))
      sum(tag.piigivej[i,]*( apply(residual^2,2,sum) + K*trace ) )
      
    }
  } else{
    gamma <- tag.gamma[,,1:q]
    GcovK <- list()
    ##   K <- dim(gcov)[2] ## needed to re-set from the previous K? K  - number of clusters?
    for(k in 1:K){
      covariance <- array(0,c(N,q,q))
      for(i in 1:N)
        covariance[i,,] <- gcov[i,k,1:q,1:q]
      GcovK[[k]]<-covariance
    }
    contr <- foreach::foreach(i = 1:N, .combine = c) %do% {
      residual <- matrix(x[curve==i],ncol=K,nrow=length(x[curve==i]))-(S[curve==i,]%*%
                                                                         (matrix(lambda.zero,ncol=K,nrow=q)+Lambda%*%t(alpha[,])+t(gamma[i,,])))
      trace <- 0
      for(k in 1:K)
        trace <- trace +sum(diag(S[curve==i,]%*%GcovK[[k]][i,,]%*%t(S[curve==i,])))*tag.piigivej[i,k]
      sum(tag.piigivej[i,]*apply(residual^2,2,sum)+trace)
    } ## estimate of sigma^2 
  } ## end else
  
  
  sigma <- sum(contr)/length(x)
  names(sigma) <- c("variance splines")
  sigma
} ## end fdvariance.mocca




variances.mocca <- function(parameters, vars, data, design){
  ## function to obtain variance estimates:
  ## 1) for model without covariates, the function returns a scalar, an estimate of sigma^2;
  ## 2) for model with covariates, the function returns two elements: 
  ##     estimates of 'variance splines' and 'variance covariates'
  
  lambda.zero <- parameters$lambda.zero
  Lambda <- parameters$Lambda 
  alpha <- parameters$alpha
  S <- design$S
  curve <- data$curve
  q <- dim(S)[2]
  x <- data$x
  
  if(!is.null(design$tag.S)) { ## model with covariates...
    tag.gamma<- vars$tag.gamma
    tag.estimates <- parameters$tag.estimates
    tag.S<-design$tag.S
    tag.piigivej <- vars$tag.piigivej
    CovMatr <- data$covariates ## matrix of covariates
    Q <- dim(tag.S)[2]
    K <- dim(tag.piigivej)[2] ## number of clusters
    N <- dim(tag.piigivej)[1]
    gcov <- vars$tag.gcovK ## array of covariances of gamma, q x q matrices. Dimension of the array is c(N,K,q+r,q+r) 
    num.of.cov <- dim(tag.estimates)[2]
  } else{ ## model without covariates...
    tag.gamma<- vars$gamma
    tag.piigivej <- vars$piigivej
    K <- dim(vars$piigivej)[2]
    N <- dim(vars$piigivej)[1]
    gcov <- vars$gcovK ## array of covariances of gamma, q x q matrices. Dimension of the array is c(N,K,q,q)
  }
  
  if(length(dim(gcov))<=3){ ## not sure if needed...?????????????
    covariance<-array(0,c(N,q,q))
    k<-0
    for(i in seq(q,q*N,by=q)){
      k<-k+1
      covariance[k,,]<-gcov[,(i-(q-1)):i]
    }
    contr <- foreach::foreach(i = 1:N, .combine = c) %do% {
      residual <- matrix(x[curve==i],ncol=K,nrow=length(x[curve==i]))-(S[curve==i,]%*%
                                                                         (matrix(lambda.zero,ncol=K,nrow=q)+Lambda%*%t(alpha[,])+t(gamma[i,,])))
      trace <- 0
      trace <- sum(diag(S[curve==i,]%*%covariance[i,,]%*%t(S[curve==i,])))
      sum(tag.piigivej[i,]*( apply(residual^2,2,sum) + K*trace ) )
    }
  } else{
    gamma <- tag.gamma[,,1:q]
    GcovK <- list()
    ##   K <- dim(gcov)[2] ## needed to re-set from the previous K? K  - number of clusters?
    for(k in 1:K){
      covariance <- array(0,c(N,q,q))
      for(i in 1:N)
        covariance[i,,] <- gcov[i,k,1:q,1:q]
      GcovK[[k]]<-covariance
    }
    contr <- foreach::foreach(i = 1:N, .combine = c) %do% {
      residual <- matrix(x[curve==i],ncol=K,nrow=length(x[curve==i]))-(S[curve==i,]%*%
                                                                         (matrix(lambda.zero,ncol=K,nrow=q)+Lambda%*%t(alpha[,])+t(gamma[i,,])))
      trace <- 0
      for(k in 1:K)
        trace <- trace +sum(diag(S[curve==i,]%*%GcovK[[k]][i,,]%*%t(S[curve==i,])))*tag.piigivej[i,k]
      sum(tag.piigivej[i,]*apply(residual^2,2,sum)+trace)
    } ## estimate of sigma^2 
  } ## end else
  
  if(!is.null(design$tag.S)){ ## model with covariates...
    temp<-gcov[,,(q+1):Q,(q+1):Q]
    tag.trace<-0
    ## for(k in 1:K){
    ##    temp.bind <- temp[,k,1,1]
    ##    if (num.of.cov>1){
    ##       for (jj in 2:num.of.cov) temp.bind <- cbind(temp.bind,temp[,k,jj,jj])
    ##    }
    ##    tag.trace <- tag.trace+apply(temp.bind,1,sum)*tag.piigivej[,k]
    ## }
    ## tag.trace<-sum(tag.trace)/(N*num.of.cov)
    
    if(num.of.cov==1){
      for(k in 1:K) tag.trace <- tag.trace+temp[,k]*tag.piigivej[,k]
      tag.trace <- sum(tag.trace)/N
    } else if(num.of.cov>1){
      for(k in 1:K){
        tag.trace <- tag.trace+apply(apply(temp[,k,1:num.of.cov,1:num.of.cov],c(1),diag),2,sum)*tag.piigivej[,k]
      }
      tag.trace <- sum(tag.trace)/(N*num.of.cov)
    }
    
    varians<-list()
    ## varians1<-list()
    for(k in 1:K){
      varians[[k]]<-(CovMatr[1:N,]-t(matrix(tag.estimates[k,],dim(CovMatr)[2],N))- vars$tag.gamma[1:N,k,c(q+1):Q])^2*vars$tag.piigivej[1:N,k] ## N x num.of.cov matrix 
      ##  varians1[[k]]<-(CovMatr[1:N,]-t(matrix(tag.estimates[k,],dim(CovMatr)[2],N)))^2*vars$tag.piigivej[1:N,k]
    }
    temp<-varians[[1]]
    ## dimcov<-dim(temp)
    ## temp1<-varians1[[1]]
    for(k in 2:K) temp <- temp+varians[[k]]  ##;temp1<-temp1+varians1[[k]]
    #temp <- apply(temp,2,sum);temp1 <- apply(temp1,2,sum)
    temp <- apply(temp,1,sum);
    ## temp1 <- apply(temp1,1,sum)
    varians<- sum(temp)/(N*num.of.cov)
    ## varians<- sum(temp)/(N*K) ## sum(temp)/(N*dimcov[2]);
    ## varians1<-sum(temp1)/(N*dimcov[2])
  } ## end if
  
  contr <- sum(contr)/length(x)
  if(!is.null(design$tag.S)){ ## if model with covariates, return  'sigma^2 splines', 'sigma^2 covariates'... 
    sigmas <- c(contr,varians+tag.trace)
    names(sigmas) <- c("variance splines","variance covariates")
  } else {
    sigmas <- contr
    names(sigmas) <- c("variance splines")
  }
  sigmas
} ## end variances.mocca





#########################
## loading functions...
##########################

print.fdaMocca.version <- function()
{ library(help=fdaMocca)$info[[1]] -> version
  version <- version[pmatch("Version",version)]
  um <- strsplit(version," ")[[1]]
  version <- um[nchar(um)>0][2]
  hello <- paste("This is fdaMocca ",version,".",sep="")
  packageStartupMessage(hello)
}


.onAttach <- function(...) { 
  print.fdaMocca.version()
}