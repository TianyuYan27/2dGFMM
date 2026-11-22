#' Fit a generalized 2D functional mixed model for Gaussian responses
#'
#' @param formula two-sided formula object in lm() format, except that the response is a matrix
#' @param data dataframe containing variables in formula and ID column;
#' dataframe contains covariates (vec or data.frame) and response (data.frame)
#' @param S number of visits
#' @param robust whether M-estimation for robust regression is used
#' @param smoother sandwich smoother (sandwich) or tensor product smoother (te)
#' @param knots=c(Sknots,Tknots) the number of knots for S and T
#' @param analytic whether to use analytical inference, only for Gaussian response
#' @param bootstrap.opt bootstrapping options
#' @param bandwidth.control correlation correction and adjustment for confidence bands bandwidth control 
#' @param pcb whether to obtain covariance estimates and pointiwse confidence bands (PCB), default is TRUE
#' @param scb whether to obtain simultaneous confidence bands (SCB), default is FALSE
#' @param silence whether to suppress logging information during program running, default is FALSE
#' @param parallel whether to run parallel computing (TRUE only for Linux/Mac users)
#' @export
#' @return a list containing estimated beta(s,t), covariance matrix of beta(s,t), qn (zero for PCB and non-zero for SCB),
#' vm (intermediate results of bootstrap for p-values)

gfmm2d_gaussian <- function(formula,
                            data,
                            S,
                            robust = TRUE,
                            smoother = "sandwich",
                            knots = NULL,
                            analytic = FALSE,
                            bootstrap.opt = list(B = 50, M = 10000, method = "nonparametric"),
                            bandwidth.control = list(d = 2, adjust = NULL, correct = 1),
                            pcb = TRUE,
                            scb = FALSE,
                            parallel = FALSE,
                            silence = FALSE)
{
  
  ##dplyr  organize lapply results
  #parallel mclapply
  #mgcv tensor product smoothing
  #refund sandwich smoothing
  #fdapace PACE
  #robustbase robust estimator 
  #zoo interpolation
  
  list.of.packages <- c("dplyr", "parallel", "mgcv", "refund", "fdapace",
                        "Matrix", "mvtnorm", "MASS", "robustbase", "zoo")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  sapply(list.of.packages, require, character=TRUE)
  
  model_formula <- as.character(formula)
  stopifnot(model_formula[1] == "~" & length(model_formula) == 3)
  vars <- unlist(strsplit(model_formula[3], split = ".\\+."))
  
  n <- nrow(data[,model_formula[2]])/S
  T <- ncol(data[,model_formula[2]])
  T.argvals <- 1:T
  S.argvals <- 1:S
  Y <- lapply(1:n, function(i) data[,model_formula[2]][((i-1)*S+1):(i*S),])
  var.mat <- NULL
  var.bin <- NULL
  is.binary.all <- function(x) {
    length(levels(x)) == 2
  }
  for(v in vars){
    if(is.array(data[,v])) var.mat <- c(var.mat, v) #determine which covariate is bivariate function
    if(is.binary.all(data[,v])) var.bin  <- c(var.bin, v) #determine which covariate is binary 
  } 
  # whether we use analytic inference or bootstrap
  t_start <- proc.time()
  t_boot <- NULL 
  
  ##############Presmoothing
  # Ysm <- list()
  # for (i in 1:n) {
  #  ori  <- as.vector(t(Y[[i]]))
  #  sm <- lowess(ori, f=.08)$y  #stats::filter(ori, sides=2, filter = (rep(1, 10)/10)) ## smooth accelerometer dat  #lowess(ori, f=.08)$y
  #  sm[is.na(sm)] <- ori[is.na(sm)]
  #  Ysm[[i]] <- t(array(sm, dim=c(T,S)))
  # }
  # Y <- Ysm
  # 
  ##########################################################################################
  ## Step 1 Bivariate pointwise estimate
  ##########################################################################################
  if (!silence) cat("Step 1: Bivariate pointwise estimate \n")
  #function "bim" fit bivariate model at t and s
  bim <- function(s,t){
    #print(paste0("s=", s, " t=", t))
    Yst <- sapply(Y, function(x) x[s,t])
    Xst <- data[seq(s, (n-1)*S + s, by=S), vars]
    
    if(length(vars) == 1){ #if there is only one covariate
      if(!is.null(var.mat)){ #if the covariate is bivariate function
        Xst <- Xst[,t]
        Xst <- array(Xst, dim=c(n, 1)) 
        colnames(Xst) <- vars
      }else{
        Xst <- array(Xst, dim=c(n, 1)) 
        colnames(Xst) <- vars
      }
    } else{ #if more than one covariate
      if(!is.null(var.mat)){ #if at least one covariates have bivariate functions
        Xst.sub <- NULL
        for(p.mat in 1:length(var.mat)){
          Xst.sub <- cbind(Xst.sub, Xst[,var.mat[p.mat]][,t])
        }
        #Xst.sub <- as.matrix(Xst[,var.mat])[,seq(t, (length(var.mat)-1)*T + t, by=T)] #subset covariates which are bivariate functions
        colnames(Xst.sub) <- var.mat
        var.others <- names(Xst)[!names(Xst) %in% var.mat]
        Xst <- cbind(Xst[,var.others,drop=FALSE], Xst.sub)
      }
    }
    dat <- data.frame(cbind(Xst, Yst))
    vars <- colnames(model.matrix(as.formula(paste(model_formula[1], model_formula[3])), data = dat))[-1]
    
    if(robust){
      fit_bi <- suppressMessages(rlm(formula = as.formula(paste0("Yst ~ ", model_formula[3])), data=dat))
    }else{
      fit_bi <- suppressMessages(lm(formula = as.formula(paste0("Yst ~ ", model_formula[3])), data=dat))
    }
    betaTilde <- data.frame(t(coef(fit_bi)))
    designmat <- model.matrix(fit_bi, na.action = na.pass)
    Covst <- summary(fit_bi)$sigma^2
    colnames(betaTilde) <- c('intercept', vars)
    return(list(betaTilde = betaTilde, designmat = designmat, cov = Covst))
  }
  
  #massive bivariate model
  massm <- list() 
  if(parallel == TRUE){
    for(i in S.argvals) massm[[i]] <- mclapply(T.argvals, bim, s=i, mc.cores = detectCores() - 1)
  }else{
    for(i in S.argvals) massm[[i]] <- lapply(T.argvals, bim, s=i)
  }
  
  if(analytic){
    betaTilde <- lapply(S.argvals, function(i) lapply(massm[[i]], '[[', 1) %>% bind_rows()) %>% bind_rows()
    designmat <- lapply(S.argvals, function(i) lapply(massm[[i]], '[[', 2)) 
    RTilde <- lapply(S.argvals, function(i) lapply(massm[[i]], '[[', 3) %>% unlist()) %>% do.call("cbind",.)
  }else{
    betaTilde <- lapply(S.argvals, function(i) lapply(massm[[i]], '[[', 1) %>% bind_rows()) %>% bind_rows()
  }
  
  ##########################################################################################
  ## Step 2 Smoothing
  ##########################################################################################
  if (!silence) cat("Step 2: Bivariate smoothing \n")
  betaHat <- array(0, dim=c(T*S, ncol(betaTilde))) 
  RHat <- array(0, dim=c(T,S)) #R is the transponse of R from the paper
  lambda <- array(0, dim=c(2, ncol(betaTilde)))
  S2kS1 <- S2 <- S1 <- list()
  
  ### smooth raw covariance estimates (estimates to measurement error)
  if(is.null(knots)) knots <- c(min(round(S/1.5), 35), min(round(T/4), 35))
  Sknots <- knots[1]
  Tknots <- knots[2]

  #obtain S1 and S2 from Sandwich smoother
  sL <- refund:::pspline.setting(x=S.argvals, knots=refund:::select_knots(S.argvals,Sknots), p=3, m=2, periodicity=FALSE, weight=NULL)
  tL <- refund:::pspline.setting(x=T.argvals, knots=refund:::select_knots(T.argvals,Tknots), p=3, m=2, periodicity=FALSE, weight=NULL)

  B1 <- tL$B
  Bt1 <- Matrix(t(as.matrix(B1))) #B1^T
  P1 <- tL$P #D1^T %*% D1

  B2 <- sL$B
  Bt2 <- Matrix(t(as.matrix(B2))) #B2^T
  P2 <- sL$P #D1^T %*% D1

  #Sandwich/tensor product smoother
  for (p in 1:ncol(betaTilde)){
    sm <- array(betaTilde[[p]], dim=c(T,S))
    if (any(is.na(sm))) sm <- na.approx(sm, na.rm = F, rule =2)
    y <- fbps(sm, knots = c(Tknots,Sknots))
    if(smoother == "te"){
      te.sm <- data.frame(y = as.vector(sm), t=rep(1:T, S), s = rep(1:S, each=T)) 
      betaHat[,p] <- gam(y ~ te(t, s, k = c(Tknots,Sknots)), data = te.sm)$fitted.values
    }else if(smoother == "sandwich"){
      betaHat[,p] <- y$Yhat %>% as.vector()
    }else{ 
      stop("Smoother must be either sandwich smoother or tensor product smooths!")
    }
    lambda[,p] <- y$lambda
    S1[[p]] <- B1 %*% solve(Bt1 %*% B1 + lambda[1,p]*P1) %*% Bt1
    S2[[p]] <- B2 %*% solve(Bt2 %*% B2 + lambda[2,p]*P2) %*% Bt2
    S2kS1[[p]] <- kronecker(S2[[p]], S1[[p]]) #S2 \otimes S1
  }
  betaHat <- lapply(1:ncol(betaHat), function(p) array(betaHat[,p], dim =  c(T, S))) #betaHat is the transpose of betaHat from paper
  names(betaHat) <-  names(betaTilde)
  colnames(lambda) <- names(betaTilde)
  
  if(analytic){
    RHat <- bandwidth.control$correct * fbps(RTilde, knots = c(Tknots,Sknots))$Yhat
    RHat[which(RHat < 0)] <- 0
    rm(betaTilde, RTilde)
  }else{
    rm(betaTilde)
  }
  t_fit <- proc.time()
  
  ##########################################################################################
  ## Step 3 Covariance Estimation
  ##########################################################################################
  
  if(pcb){
    
    #fpca.opt A list of options control parameters specified by list, dataType: dense and regular (Dense); Very densely and regularly 
    #observed data: empirical mean and Densely recorded but irregular design, or contaminated with error: pre-smoothing for individual 
    #curves (DenseWithMV); Sparse random design (Sparse)
    
    fpca.opt <- list(dataType = 'Dense', methodSelectK = 'FVE')
    B <- bootstrap.opt$B #bootstrap times
    M <- bootstrap.opt$M #wild bootstrap residuals
    para <- bootstrap.opt$method
    
    ##########################################################################################
    ## Step 3.1 Centering
    ##########################################################################################
    if (!silence) cat("Step 3: Computing variance of bivariate estimates \n")
    if(analytic)
    {
      ## Use analytic solution for covariance estimation
      ## Get eta to do FPCA
      get_eta <- function(s,t){
        Yst <- sapply(Y, function(x) x[s,t])
        betaHat.st <- sapply(1:(length(vars)+1), function(p) betaHat[[p]][t,s]) %>% array(dim = c(length(vars)+1, 1))
        Ytilde <- designmat[[s]][[t]] %*% betaHat.st
        eta <- Yst - Ytilde
        return(eta = eta)
      }
      #massive eta
      massEta <- list()
      if(parallel == TRUE){
        for(i in S.argvals) massEta[[i]] <- mclapply(T.argvals, get_eta, s=i, mc.cores = detectCores() - 1)
      }else{
        for(i in S.argvals) massEta[[i]] <- lapply(T.argvals, get_eta, s=i)
      }
      
      eta <- NULL
      for(i in 1:n) {
        eta.i <- lapply(S.argvals, function(s) lapply(T.argvals, function(t) massEta[[s]][[t]][i])) %>% unlist() #ith sample
        eta <- rbind(eta, array(eta.i, dim=c(T,S)))
      } 
      
      ##########################################################################################
      ## Step 3.2 Computing marginal covariance w.r.t. S and eigenfunctions psi
      ##########################################################################################
      
      eigenS <- s.fpca(eta, S, fpca.opt)
      s.cov <- eigenS$s.cov
      pc.s <- eigenS$pc.s
      xiEst <- eigenS$xiEst
      
      ##########################################################################################
      ## Step 3.3 Computing B-spline smoothing xi(t) and marginal covariance w.r.t. T 
      ##########################################################################################
      
      Bt <- as.matrix(B1)
      bsplineT <- b.lm(n, Bt, pc.s, xiEst)
      t.cov <- bsplineT$t.cov
      
      ##########################################################################################
      ## Step 3.4 Computing four-dim covariance operator w.r.t. beta(s,t) 
      ##########################################################################################
      
      covBeta <- function(s,u,t,v){ #set s,u,t,v
        if(s==u & t==v){
          cEta <- RHat[t,s] + Reduce("+", lapply(1:pc.s,function(j) s.cov[j,s,u]*t.cov[j,t,v]))
        }else{
          cEta <- Reduce("+", lapply(1:pc.s,function(j) s.cov[j,s,u]*t.cov[j,t,v]))
        }
        xst1 <- designmat[[s]][[t]]
        xuv2 <- designmat[[u]][[v]]
        cBeta <- cEta * solve(t(xst1) %*% xst1) %*% t(xst1) %*% xuv2 %*% solve(t(xuv2) %*% xuv2) 
        cBeta
      }
      
      ## raw estimate of Cov(beta(s,t), beta(s,t))
      cov.beta.ts.tilde <-  array(0, dim = c(S*T, S*T, length(betaHat)))
      if(parallel == TRUE){
        ## Bootstrap helper function
        covBeta.parallel <- function(s.par, u.par){
          cov.beta <- array(0, dim = c(T, T, length(betaHat)))
          for(t in 1:T){
            for(v in t:T){
              tmp <- covBeta(s.par,u.par,t,v)
              for(p in 1:length(betaHat)){
                cov.beta[t, v, p] <- 
                  cov.beta[v, t, p] <- tmp[p,p]
              }
            }
          }
          return(cov.beta)
        }
        #parallel computing S domain
        tmp <- list()
        for(s in S.argvals) {
          tmp[[s]] <- mclapply(s:S, covBeta.parallel, s.par=s, mc.cores = detectCores() - 1)
          for(up in 1:length(tmp[[s]])){
            for(p in 1:length(betaHat)){
              # put each element to the right position
              u <- s - 1 + up
              cov.beta.ts.tilde[((s-1)*T+1):(s*T), ((u-1)*T+1):(u*T), p] <-
                cov.beta.ts.tilde[((u-1)*T+1):(u*T), ((s-1)*T+1):(s*T), p] <- tmp[[s]][[up]][,,p]  
            }
          }
        }
      }else{
        for(s in 1:S){
          for(u in s:S){
            for(t in 1:T){
              for(v in t:T){
                tmp <- covBeta(s,u,t,v)
                for(p in 1:length(betaHat)){
                  cov.beta.ts.tilde[(u-1)*T+v, (s-1)*T+t, p] <- 
                    cov.beta.ts.tilde[(u-1)*T+t, (s-1)*T+v, p] <- 
                    cov.beta.ts.tilde[(s-1)*T+v, (u-1)*T+t, p] <- 
                    cov.beta.ts.tilde[(s-1)*T+t, (u-1)*T+v, p] <- tmp[p,p]
                }
              }
            }
          }
        }
      }
      
    }else{
      
      # Bootstrap refits stay serial. The outer mclapply already uses every core;
      # a parallel gfmm2d inside each replicate forks another full set of processes.
      if (!silence) cat(paste("bootstrapping... \n"))
      
      if(parallel == TRUE){
        betaHat.boot <- mclapply(1:B, function(b) fmm2d.boot(b, data, betaHat, n, formula, para,
                                                             S, robust, smoother, knots,
                                                             bandwidth.control, FALSE,
                                                             pcb, scb, silence, var.bin), mc.cores = detectCores() - 1)
      }else{
        betaHat.boot <- lapply(1:B, function(b) fmm2d.boot(b, data, betaHat, n, formula, para,
                                                           S, robust, smoother, knots,
                                                           bandwidth.control, FALSE,
                                                           pcb, scb, silence, var.bin))
      }
      if (!silence) cat(paste("bootstrapping finished \n"))
      t_boot <- proc.time()
      
      gamma <- 1.2

      ## obtain bootstrap variance
      cov.beta.ts.tilde <- array(0, dim = c(T*S, T*S, length(betaHat)))
      for (p in 1:length(betaHat)){
        betaHat.boot.p <- lapply(1:B, function(b) betaHat.boot[[b]][[p]])
        cov.beta.ts.tilde[,,p] <- gamma*cov(do.call(rbind, lapply(betaHat.boot.p, function(x) c(x)))) ## account for within-subject correlation
      }
    }
    
    ## refined smoothing estimates covariance 
    cov.beta.ts.hat <- array(0, dim = c(T*S, T*S, length(betaHat)))
    for(p in 1:length(betaHat)){
      edcomp <- eigen(cov.beta.ts.tilde[,,p])
      eigen.positive <- which(edcomp$values > 0)
      if(length(eigen.positive) == S*T){
        ker.trimmed <- cov.beta.ts.tilde[,,p]
      }else{
        ker.trimmed <- edcomp$vectors[,eigen.positive] %*% diag(edcomp$values[eigen.positive]) %*% t(edcomp$vectors[,eigen.positive])
      }
      cov.beta.ts.hat[,,p] <- matrix(S2kS1[[p]] %*% ker.trimmed %*% t(S2kS1[[p]])) ## account for within-subject correlation
    }
    rm(ker.trimmed, cov.beta.ts.tilde)
    
    
    qn <- rep(0, length = length(betaHat))
    vm <- array(0, dim = c(length(betaHat), M))
    
    ##########################################################################################
    ## Step 3.4(optional) Simultaneous Confidence Bands 
    ##########################################################################################
    if(scb)
    {   
      if (!silence) cat("Step 3 (optional): Preparing simultaneous confidence bands \n")
      
      if(!exists("betaHat.boot")){
        #if no boostrapping betaHat before
        if (!silence) cat(paste("bootstrapping... \n"))
        
        if(parallel == TRUE){
          betaHat.boot <- mclapply(1:B, function(b) fmm2d.boot(b, data, betaHat, n, formula, para,
                                                               S, robust, smoother, knots,
                                                               bandwidth.control, FALSE,
                                                               pcb, scb, silence, var.bin), mc.cores = detectCores() - 1)
        }else{
          betaHat.boot <- lapply(1:B, function(b) fmm2d.boot(b, data, betaHat, n, formula, para,
                                                             S, robust, smoother, knots,
                                                             bandwidth.control, FALSE,
                                                             pcb, scb, silence, var.bin))
        }
        if (!silence) cat(paste("bootstrapping finished \n"))
      }
      
      #computing marginal decomposition on FPCA for beta.boot
      psi.boot <- list()
      ker.boot <- list()
      betaHat.boot.avg <- list()
      Bt <- as.matrix(B1)
      for (p in 1:length(qn)){
        betaHat.boot.p <- lapply(1:B, function(b) betaHat.boot[[b]][[p]])
        betaHat.boot.avg[[p]] <- Reduce("+", betaHat.boot.p)/B
        betaHat.boot.p.demean <- lapply(1:B, function(b) betaHat.boot.p[[b]] - betaHat.boot.avg[[p]]) #demean each betaHat.boot.ps
        eigenS.boot <- s.fpca(eta = do.call("rbind", betaHat.boot.p.demean), S = S, fpca.opt = fpca.opt)
        psi.boot[[p]] <- eigenS.boot$psi
        bsplineT.boot <- b.lm(size = B, Bt = Bt, pc.s = eigenS.boot$pc.s, xiEst = eigenS.boot$xiEst) 
        ker.boot[[p]] <- bsplineT.boot$ker
      }
      rm(betaHat.boot.p, betaHat.boot.p.demean)
      
      for(p in 1:length(qn)){
        if(analytic){
          sp <- sqrt(diag(cov.beta.ts.hat[,,p]))
        }else{
          sp <- sqrt(diag(cov.beta.ts.hat[,,p])/gamma)
        }
        qm <- lapply(1:M, function(m) array(0, dim=c(T,S)))
        JB <- dim(ker.boot[[p]])[1]
        for(j in 1:JB){
          #since extra variance is brought by B-spline regression,
          #we need to remove it, which reduces the varaince
          ifelse(is.null(bandwidth.control$adjust), Sig <- 1/(B*2.4)*ker.boot[[p]][j,,],
                 Sig <- 1/(B*bandwidth.control$adjust)*ker.boot[[p]][j,,])
          d <- bandwidth.control$d #length of removing tail 
          Sig.r.t <- Sig[-c(1:d,(ncol(Bt)-d+1):ncol(Bt)), -c(1:d,(ncol(Bt)-d+1):ncol(Bt))] #remove tail problems
          Bt.r <- Bt[,-c(1:d,(ncol(Bt)-d+1):ncol(Bt))] #remove tail problems
          #further smoothing covariance matrix to avoid large umat
          Sig.r <- fbps(Sig.r.t, knots = round(ncol(Bt.r)/1.5))$Yhat
          if(isSymmetric(Sig.r)){
            umat <- rmvnorm(M, mean = rep(0, ncol(Bt.r)), sigma = Sig.r)
          }else{
            umat <- rmvnorm(M, mean = rep(0, ncol(Bt.r)), sigma = Sig.r.t)
          }
          Btu <- apply(umat, 1, function(x) Bt.r %*% x)  
          qm_j <- lapply(1:M, function(m) kronecker(Btu[,m], t(psi.boot[[p]][,j])))
          qm <- mapply("+", qm_j, qm, SIMPLIFY = FALSE)
          betaHat.boot.p <- lapply(qm, function(q) q + betaHat.boot.avg[[p]]) 
        }
        qmdvar <- lapply(1:M, function(m) abs(betaHat.boot.p[[m]] - betaHat[[p]])/array(sp, dim=c(T,S)))
        qmstar <- lapply(1:M, function(m) max(qmdvar[[m]]))
        vm[p,] <- unlist(qmstar)
        qn[p] <- quantile(vm[p,], 0.95) 
      }
    }
    ## results for pcb/scb
    out <- list(betaHat = betaHat, betaHat.cov = cov.beta.ts.hat, qn = qn, vm = vm,
                time = gfmm2d_elapsed(t_start, t_fit, t_boot))
    class(out) <- "gfmm2d"
    return(out)
  }else{
    ## results for no pcb
    out <- list(betaHat = betaHat, time = gfmm2d_elapsed(t_start, t_fit, t_boot))
    class(out) <- "gfmm2d"
    return(out)
  }
}

## Wall-clock seconds. fit is pointwise estimation plus smoothing.
## bootstrap is the B subject-level refits. bands is the covariance
## assembly and the simultaneous-band Monte Carlo, which reuse those refits.
gfmm2d_elapsed <- function(t_start, t_fit, t_boot) {
  t_end <- proc.time()
  fit <- unname((t_fit - t_start)["elapsed"])
  if (is.null(t_boot)) {
    bootstrap <- 0
    bands <- unname((t_end - t_fit)["elapsed"])
  } else {
    bootstrap <- unname((t_boot - t_fit)["elapsed"])
    bands <- unname((t_end - t_boot)["elapsed"])
  }
  c(fit = fit, bootstrap = bootstrap, bands = bands,
    total = unname((t_end - t_start)["elapsed"]))
}
