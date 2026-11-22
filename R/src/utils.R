s.fpca <- function(eta, S, fpca.opt){
  sList <- list()
  etaList <- list()
  S.argvals <- 1:S
  
  for (nt in 1:nrow(eta)){
    ind <- which(!is.na(eta[nt,]))
    etaList[[nt]] <- eta[nt,ind]
    sList[[nt]] <- S.argvals[ind]
  }
  
  fpca <- suppressWarnings(FPCA(etaList, sList, optns = fpca.opt)) #fpca$cumFVE is cumulative variance
  psi <- fpca$phi 
  xiEst <- fpca$xiEst 
  #tau <- fpca$lambda
  pc.s <- ncol(xiEst)
  s.cov <- array(0, dim=c(pc.s, S, S))
  for(j in 1:pc.s){
    s.cov[j,,] <- psi[, j] %*% t(psi[,j]) #fpca$fittedCov
  }
  return(list(s.cov=s.cov, pc.s=pc.s, xiEst=xiEst, psi=psi))
}


b.lm <- function(size, Bt, pc.s, xiEst){
  b <- array(0, dim=c(size,ncol(Bt),pc.s))
  
  bspl.fit <- function(i){
    fit <- lm(bspl[i,] ~ 0 + Bt)
    b <- coef(fit)
    return(b)
  }
  
  T <- nrow(Bt)
  for (j in 1:pc.s){
    bspl <- matrix(xiEst[,j], nrow=size, ncol=T, byrow=TRUE)
    bspl.result <- lapply(1:size, function(i) bspl.fit(i))
    b[,,j] <- bspl.result %>% bind_rows() %>% as.matrix()
  }
  
  #store marginal covariance of T 
  t.cov <- array(0, dim=c(pc.s, T, T))
  ker <- array(0, dim=c(pc.s, ncol(Bt), ncol(Bt)))
  for(j in 1:pc.s){
    for (i in 1:size){
      ker[j,,] <- ker[j,,] + kronecker(t(b[i,,j]), b[i,,j])
    }
    t.cov[j,,] <- (1/size) * Bt %*% ker[j,,] %*% t(Bt)
  }
  return(list(t.cov=t.cov, ker=ker))
}


## Use bootstrap for covariance estimation
fmm2d.boot <- function(b, data, betaHat, n, formula, para,
                       S, robust, smoother, knots,
                       bandwidth.control, parallel,
                       pcb, scb, silence, var.bin){
  i <- 0
  if (b %% 10 == 0) cat(b, "\n")
  check_legal <- function(data.boot, var.bin){
    legal <- rep(NA, length(var.bin))
    for(i in 1:length(var.bin)) legal[i] <- all(data.boot[,var.bin[i]] == 0) | all(data.boot[,var.bin[i]] == 1)
    return(all(legal == FALSE))
  }
  
  while(TRUE){
    i <- i + 1
    sample.ind <- sample(unique(data$ID), size = n, replace = TRUE) #bootstrap id with replacement
    row.ind <- NULL
    if(para == "nonparametric"){
      for (id in 1:length(sample.ind)){
        row.ind <- c(row.ind, which(data$ID == sample.ind[id]))
      }
      data.boot <- data[row.ind,] #b_th dataset
    }else{
      stop("Other bootstrap methods are not avaliable so far.")
    }
    
    if(check_legal(data.boot, var.bin) | is.null(var.bin)){
      #if all binary covariates are valid or no binary covariate
      fmm2d.boot.result <- try(gfmm2d_gaussian(formula = formula, data = data.boot,
                                      S = S, robust = robust, smoother = smoother, knots = knots,
                                      bandwidth.control = bandwidth.control, parallel = parallel,
                                      pcb = FALSE, scb = FALSE, silence = TRUE), silent = TRUE)
      
      if ('try-error' %in% class(fmm2d.boot.result)){
        #if fmm2d encounters error
        next
      }else{
        break
      }
    }else{
      next
    }
    if(i == 100) stop("Please check the binary covariates if they contains both categories.")
    gc()
  }
  return(fmm2d.boot.result$betaHat)
}
