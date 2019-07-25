"gibbSamplerGP" <- 
  function(priors = list(beta = list(m0 = 0, s0 = 1, dist = "gauss", constrained = FALSE),
                         thetaGP = list(k1 = list(a1 = 20, b1=10, a2=10, b2 = 500),
                                        k2=list(a1 = 20, b1=10, a2=10, b2 = 500,a3 = 20, b3=10, a4=10, b4 = 500)),
                         gamma = list(m0 = 0, s0 = 1, dist = "gauss"),
                         sigma2 = list(a0=2, b0=1, dist = "gamma"),
                         tau = list(a0 = 2, b0= 1, dist = "gamma"),
                         phi = list(inf = 0.1, sup = 1, dist = "unif"),
                         pi = list(alpha0 = matrix(1,ncol=1), dist = "dirichlet"),
                         rho = list(inf=0.1, sup = 3, dist = "unif")),
           y, basis, ID,coords,cov,SpTcov,Bds = NULL,
           N.run = 10000, nBatch = nBatch,
           debug = FALSE, print.res = FALSE,
           parallel = parallel, nCluster = nCluster, kernelList = NULL,
           ...){
    
    
    cov = cbind(cov, rep(1, nrow(y)))
    nBasis = 2
    nSites = length(ID)
    nCov = ncol(cov)
    nObs = nrow(y)
    nComp = ncol(priors$pi$alpha0)
    nKernels = length(priors$thetaGP)
    nParGP = length(unlist(priors$thetaGP)) / 2
    
    if(is.null(kernelList)){
      kernelList = list(k_longterm, k_seasonal)
      attr(kernelList, "name") <- c("long term", "seasonal")
      attr(kernelList, "parameters") <- list(c("q1","q2"), c("q3", "q4", "qs", "f"))
      attr(kernelList, "type") <- c("temporal", "temporal")
    }
    
    if (nKernels != length(kernelList)) {
      stop("Discrepancy between GP kernels and the number of hyperparameters")
    }
    
    if(is.null(SpTcov)){SpTcov = matrix(1, ncol = nComp, nrow = nObs)}
    
    SpTcov = as.matrix(SpTcov)
    
    nSpTCov = ncol(SpTcov)
    
    betaHist = matrix(NA, nrow = N.run, ncol = nSites * nBasis * nCov)
    gammaHist = matrix(NA, nrow = N.run, ncol = nSpTCov)
    thetaGPHist = matrix(NA, nrow = N.run, ncol = nParGP)
    sigma2Hist = matrix(NA, nrow = N.run, ncol = 1)
    tauHist = matrix(NA, nrow = N.run, ncol = 1)
    phiHist = matrix(NA, nrow = N.run, ncol = 1)
    rhoHist = matrix(NA, nrow = N.run, ncol = 1)
    piHist = array(NA, dim = c(N.run, nComp, nSites))
    zHist = array(NA, dim = c(N.run, nComp, nSites))
    logPostDist = matrix(NA, nrow = N.run,ncol = 1)
    logLike = matrix(NA, nrow = N.run,ncol = 1)
    dpriorsHist = matrix(NA, nrow = N.run, ncol = 8)
    yHatHist = matrix(NA, nrow = N.run, ncol = nObs)
    
    betaH = array(NA, dim = c(dim(betaHist), nBatch))
    gammaH = array(NA, dim = c(dim(gammaHist), nBatch))
    thetaGPH = array(NA, dim = c(dim(thetaGPHist),nBatch))
    sigma2H = array(NA, dim = c(dim(sigma2Hist), nBatch))
    tauH = array(NA, dim = c(dim(tauHist), nBatch))
    phiH = array(NA, dim = c(dim(phiHist), nBatch))
    rhoH = array(NA, dim = c(dim(rhoHist), nBatch))
    piH = array(NA, dim = c(dim(piHist), nBatch))
    zH = array(0, dim = c(dim(zHist), nBatch))
    zH = as.simple_sparse_array(zH)
    logPostDistH= array(NA, dim = c(N.run,1, nBatch))
    logLikeH= array(NA, dim = c(N.run,1, nBatch))
    yHatH = array(NA, dim = c(dim(yHatHist), nBatch))
    dpriorsH = array(NA, dim = c(dim(dpriorsHist), nBatch))
    
    Xtilde = simple_triplet_zero_matrix(nrow = nObs, ncol = nSites*nCov*nBasis)
    
    DistMat = as.matrix(dist(coords[match(coords$site, ID) ,2:3]))
    DistMatBeta = kronecker(DistMat,diag(nBasis))
    
    Mu = y$obs
    
    idxID = NULL
    for(si in 1:nSites){
      idxID = c(idxID, list(which(y$ID == ID[si])))
    }
    
    tic.clearlog()
    for(nb in 1:nBatch){
      tic()
      # Initiate RV
      
      # thetaGP0  = c(rgamma(1,priors$thetaGP$k1$a1,priors$thetaGP$k1$b1),
      #              rgamma(1,priors$thetaGP$k1$a2,priors$thetaGP$k1$b2),
      #              rgamma(1,priors$thetaGP$k2$a1,priors$thetaGP$k2$b1),
      #              rgamma(1,priors$thetaGP$k2$a2,priors$thetaGP$k2$b2),
      #              rgamma(1,priors$thetaGP$k2$a3,priors$thetaGP$k2$b3),
      #              rgamma(1,priors$thetaGP$k2$a4,priors$thetaGP$k2$b4))
      thetaGP0 = rgamma(nParGP, unlist(priors$thetaGP)[(1:nParGP)*2-1], unlist(priors$thetaGP)[(1:nParGP)*2])
      
      
      tau0 = rigamma(1,priors$tau$a0, priors$tau$b0)
      phi0 = runif(1,priors$phi$inf, priors$phi$sup)
      rho0 = runif(1,priors$rho$inf, priors$rho$sup)
      pi0 = matrix(round(t(apply(priors$pi$alpha0, 1, bayesm::rdirichlet)),5), ncol = nComp, byrow=F)
      z0 = matrix(apply(pi0, 1, function(x) rmultinom(1,1,x)), ncol=nComp, byrow=T)
      
      comp = apply(matrix(z0, nrow=nComp), 2 , which.max)
      
      #browser()
      betaPostMu =rep(priors$beta$m0,nSites * nBasis * nCov)
      betaPostSd = kronecker(diag(nBasis),tau0*exp(- phi0 * DistMat)* outer(comp, comp, FUN = function(x,y) as.numeric(x==y)))
      if(priors$beta$constrained){
        
        # Using Rue2005 alg.
        Ab = kronecker(t(z0), diag(nBasis))
        idxNC = which(rowSums(Ab)==0)
        if(length(idxNC)>0){
          Ab = matrix(Ab[-idxNC,], ncol = nSites*nBasis)
        }
        eb = matrix(c(rowSums(Ab)),ncol=1, nrow=nrow(Ab))
        eb[(1:(nrow(eb)/2))*2-1] <- 0
        L = t(chol(betaPostSd))
        zi = rmvnorm(1,rep(0, nrow(betaPostSd)), diag(1, nrow(betaPostSd)))
        v = solve(t(L),t(zi))
        xi = betaPostMu + v
        Vnk = (betaPostSd) %*% t(Ab)
        Wkk = Ab %*% Vnk
        Ukn = solve(Wkk) %*% t(Vnk)
        ci = Ab %*% xi - eb
        beta0 = xi - t(Ukn) %*% ci
        muBeta0 = betaPostMu - (betaPostSd) %*% t(Ab) %*% solve(Ab %*% (betaPostSd) %*% t(Ab)) %*% (Ab %*% betaPostMu - eb)
        #browser()
      }else{
        beta0 = rmvnorm(1, betaPostMu, betaPostSd)
      }
      
      gamma0 = rnorm(nSpTCov,  priors$gamma$m0, priors$gamma$s0)
      sigma20 = rigamma(1,priors$sigma2$a0, priors$sigma2$b0)
      
      
      betaHist[1,] = beta0
      gammaHist[1,] = gamma0
      sigma2Hist[1] = sigma20
      tauHist[1] = tau0
      rhoHist[1] = rho0
      phiHist[1] = phi0
      piHist[1,,] = t(pi0)
      zHist[1,,] = t(z0)
      thetaGPHist[1,] = thetaGP0
      wll_old = -Inf
      
      #if(nSpTCov>0){MuG = Mu - SpTcov %*% gammaHist[1,]}else{MuG = Mu}
      
      if(nComp>1){SpTcov = z0[match(y$ID, ID),]}
      
      if(nSpTCov>0){
        MuA = MuB = Mu - SpTcov %*% gammaHist[1,]
        MuG = Mu - crossprod_simple_triplet_matrix(t(Xtilde), betaHist[1,])
      }else{
        MuA = MuB = Mu
      }
      
      numLog = NULL
      for(t in 2:N.run){
        #print(t)
        
        # Updating alpha_j
        zCurr = matrix(zHist[t-1,,], nrow=nComp)
        
        idxC = ((si-1)*nBasis+1):(si*nBasis)
        
        MuCorr = MuA#diag(betaHist[t-1,match(y$ID, levels(y$ID))*nBasis]^(-1)) %*% (MuA - betaHist[t-1,match(y$ID, levels(y$ID))*nBasis-1])
        
        thetaGP_Prop = c(rgamma(nParGP, unlist(priors$thetaGP)[(1:nParGP)*2-1],
                                unlist(priors$thetaGP)[(1:nParGP)*2]),
                         sigma2Hist[t-1])
        
        
        
        wll = 0
        for(i in 1:nlevels(y$ID)){
          xList = list(matrix(y$date[y$ID == levels(y$ID)[i]], ncol=1),matrix(y$date[y$ID == levels(y$ID)[i]], ncol=1))
          wll  = wll + wrap_ll(theta = (thetaGP_Prop), y = MuCorr[y$ID == levels(y$ID)[i]], kernelList = kernelList, xList = xList)
        }
        
        ratio = exp(-wll_old +wll)
        
        thetaGPHist[t,] = thetaGPHist[t-1,]
        
        if(ratio > runif(1)){
          thetaGPHist[t,] = thetaGP_Prop[-length(thetaGP_Prop)]
          wll_old = wll
        }
        
        param = parToList(c(thetaGPHist[t,],sigma2Hist[t-1]), kernelList)
        
        ft = NULL
        for(j in 1:nComp){
          
          idxZ = which(y$ID %in% ID[which(apply(zCurr,2,which.max) == j)])
          if(length(idxZ)>0){
            
            xList = list(matrix(y$date[idxZ], ncol=1),matrix(y$date[idxZ], ncol=1))
            xPred = list(matrix(y$date, ncol=1),matrix(y$date, ncol=1))
            tt = GPpred(xd = xPred, x = xList,y = c(MuCorr[idxZ]), param = param, kernel = kernelList)
            ft = c(ft, list(tt$mp))
          }else{
            xList = list(matrix(y$date, ncol=1),matrix(y$date, ncol=1))
            xPred = list(matrix(y$date, ncol=1),matrix(y$date, ncol=1))
            tt = GPpred(xd = xPred, x = xList,y = c(MuCorr), param = param, kernel = kernelList)
            ft = c(ft, list(tt$mp))
          }
        }
        #browser()
        # Sample/Create latent temporal functionals
        
        for(si in 1:nSites){
          fTime = ft[[which.max(zCurr[,si])]]
          xt = cbind(1,fTime * cov)
          idxC = ((si-1)*nBasis+1):(si*nBasis)
          Xtilde[which(y$ID == ID[si]),idxC] = xt[which(y$ID == ID[si]),]
        }
        
        tXX = crossprod_simple_triplet_matrix(Xtilde)
        
        # Updating beta_j
        
        comp = apply(matrix(zHist[t-1,,], nrow=nComp), 2 , which.max)
        
        #SigmaBeta = kronecker(diag(nBasis),tauHist[t-1]*exp(- phiHist[t-1] * DistMat)* outer(comp, comp, FUN = function(x,y) as.numeric(x==y)))
        iSigmaBeta = kronecker(diag(nBasis),1/tauHist[t-1]*solve(exp(- phiHist[t-1] * DistMat)* outer(comp, comp, FUN = function(x,y) as.numeric(x==y))))
        
        iSigmaObs = tXX
        
        betaPostSd = solve(iSigmaBeta + iSigmaObs)
        ibetaPostSd = iSigmaBeta + iSigmaObs
        
        betaPostMu = betaPostSd %*% (crossprod_simple_triplet_matrix(Xtilde, MuB) + 
                                       iSigmaBeta %*% matrix(priors$beta$m0, ncol =1, nrow = nrow(iSigmaBeta)) )
        
        if(priors$beta$constrained){
          
          # Using Rue2005 alg.
          #Ab = kronecker(zCurr, matrix(c(0,1),ncol=2, nrow=1))
          Ab = kronecker(zCurr, diag(nBasis))
          idxNC = which(rowSums(Ab)==0)
          if(length(idxNC)>0){
            Ab = matrix(Ab[-idxNC,], ncol = nSites*nBasis)
          }
          eb = matrix(c(rowSums(Ab)),ncol=1, nrow=nrow(Ab))
          eb[(1:(nrow(eb)/2))*2-1] <- 0
          L = t(chol(ibetaPostSd))
          zi = rmvnorm(1,rep(0, nrow(betaPostSd)), diag(1, nrow(betaPostSd)))
          v = solve(t(L),t(zi))
          xi = betaPostMu + v
          Vnk = (betaPostSd) %*% t(Ab)
          Wkk = Ab %*% Vnk
          Ukn = solve(Wkk) %*% t(Vnk)
          ci = Ab %*% xi - eb
          betaHist[t,] = xi - t(Ukn) %*% ci
          
          # New cond dist parameters 
          
          muBetaC = betaPostMu - betaPostSd %*% t(Ab) %*% solve(Ab %*% betaPostSd %*% t(Ab)) %*% (Ab %*% betaPostMu - eb)
          sigBetaC = betaPostSd - betaPostSd %*% t(Ab) %*% solve(Ab %*% betaPostSd %*% t(Ab)) %*% Ab %*% betaPostSd
          EVD_beta = eigen(sigBetaC/2 + t(sigBetaC)/2)
          iEV_beta = 1 / EVD_beta$values
          iEV_beta[EVD_beta$values <= 0] = 0
          SigBM = EVD_beta$vectors %*% diag(iEV_beta) %*% t(EVD_beta$vectors)
          
        }else{
          betaHist[t,] = rmvnorm(1, betaPostMu, betaPostSd)
        }
        
        # Updating gamma_j
        
        
        Yp = crossprod_simple_triplet_matrix(t(Xtilde), betaHist[t,])
        
        if(nSpTCov>0){
          
          tWW = crossprod(SpTcov, SpTcov)
          iSigmaGamma = 1/priors$gamma$s0 * diag(ncol(gammaHist))
          
          iSigmaObsG = tWW
          
          gammaPostSd =solve(iSigmaGamma + iSigmaObsG)
          gammaPostMu = gammaPostSd %*% (crossprod(SpTcov, MuG) + 
                                           iSigmaGamma  %*% matrix(priors$gamma$m0, ncol =1, nrow = nrow(iSigmaGamma)))
          
          gammaHist[t,] = rmvnorm(1, gammaPostMu, gammaPostSd)
          
        }
        
        # Updating sigma
        
        if(nSpTCov>0){Ypp = Yp + SpTcov %*% gammaHist[t-1,]}else{Ypp = Yp}
        
        sigma2PostA = priors$sigma2$a0 + nObs / 2
        sigma2PostB = priors$sigma2$b0 + sum((Mu - Ypp)^2, na.rm=T) / 2
        
        #if(t==10){browser()}
        
        sigma2Hist[t] = rigamma(1, sigma2PostA, sigma2PostB)
        
        # Updating yHat
        # 
        yHatHist[t,] = Ypp
        
        # Updating tau
        
        remSite = which(comp %in% which(rowSums(zCurr) <= 1))
        #idxBeta = ((remSite-1)*nBasis+1):(remSite*nBasis)
        
        betaMu = matrix(priors$beta$m0, ncol=1, nrow = length(betaPostMu)) #betaPostMu
        SigmaBT = kronecker(diag(nBasis),exp(- phiHist[t-1] * DistMat))
        iSigmaBT = kronecker(diag(nBasis),solve(exp(- phiHist[t-1] * DistMat)))
        nConst = 0
        
        tauPostA = priors$tau$a0 + ((nSites-length(remSite))*nBasis - nConst) / 2
        if(!priors$beta$constrained){tauPostB = priors$tau$b0 + (t(as.matrix(betaHist[t,]) - betaMu) %*% iSigmaBT %*% as.matrix(betaHist[t,]- betaMu))/2}
        
        if(priors$beta$constrained){
          betaMu = muBeta0
          nConst = nrow(Ab)
          sigBetaC = SigmaBT - SigmaBT %*% t(Ab) %*% solve(Ab %*% SigmaBT %*% t(Ab)) %*% Ab %*% SigmaBT
          EVD_beta = eigen(sigBetaC/2 + t(sigBetaC/2))
          iEV_beta = 1 / EVD_beta$values
          iEV_beta[EVD_beta$values <= 1e-4] = 0
          SigBM = EVD_beta$vectors %*% diag(iEV_beta) %*% t(EVD_beta$vectors)
          iSigmaBT = SigBM
          tauPostB = priors$tau$b0 + (t(as.matrix(betaHist[t,]) - betaMu) %*% iSigmaBT %*% as.matrix(betaHist[t,]- betaMu))/2
        }
        
        # Only the clusters with more than one site provide information about tau and phi
        
        
        if(tauPostB < 0){browser()}
        tauHist[t] = rigamma(1,tauPostA, tauPostB)
        
        # Updating phi (MH Reject)
        
        phiProp = runif(1, priors$phi$inf, priors$phi$sup)
        
        #SigmaBetaN = kronecker(diag(nBasis),tauHist[t-1]*exp(- phiProp * DistMat))
        iSigmaBetaN = kronecker(diag(nBasis),1/tauHist[t-1]*solve(exp(- phiProp * DistMat)* outer(comp, comp, FUN = function(x,y) as.numeric(x==y))))
        #iSigmaBeta = kronecker(diag(nBasis),1/tauHist[t-1]*solve(exp(- phiHist[t-1] * DistMat)* outer(comp, comp, FUN = function(x,y) as.numeric(x==y))))
        betaPostSdN = solve(iSigmaBetaN + iSigmaObs)
        ibetaPostSdN = iSigmaBetaN + iSigmaObs
        
        # betaPostMuN = betaPostSdN %*% (sigma2Hist[t-1]^(-1) * crossprod_simple_triplet_matrix(Xtilde, MuG) + 
        #                                iSigmaBetaN %*% matrix(priors$beta$m0, ncol =1, nrow = nrow(iSigmaBetaN)) )
        if(priors$beta$constrained){
          sigBetaCN = ibetaPostSdN - ibetaPostSdN %*% t(Ab) %*% solve(Ab %*% ibetaPostSdN %*% t(Ab)) %*% Ab %*% ibetaPostSdN
          sigBetaCN = sigBetaCN/2 + t(sigBetaCN)/2
          EVD_betaN = eigen(sigBetaCN)
          iEV_betaN = 1 / EVD_betaN$values
          iEV_betaN[EVD_betaN$values <= 0] = 0
          SigBMN = EVD_betaN$vectors %*% diag(iEV_betaN) %*% t(EVD_betaN$vectors)
        }
        
        
        
        numLog = -1/2 * log(det(betaPostSd)) +  -1/2*t(betaHist[t,] - betaMu) %*% solve(betaPostSd) %*% (betaHist[t,] - betaMu)
        if(priors$beta$constrained){
          numLog   =  -(length(betaHist[t,])-1)/2 * log(2*pi) -
            1/2*t(betaHist[t,]- betaMu)%*% SigBM %*%  (betaHist[t,]- betaMu) -
            1/2*sum(log(EVD_beta$values[EVD_beta$values > 0]))
        }
        
        denLog = -1/2 * log(det(betaPostSdN)) +  -1/2*t(betaHist[t,] - betaMu) %*% ibetaPostSdN %*% (betaHist[t,] - betaMu)
        if(priors$beta$constrained){
          denLog   =  -(length(betaHist[t,])-1)/2 * log(2*pi) -
            1/2*t(betaHist[t,]- betaMu)%*% SigBMN %*%  (betaHist[t,]- betaMu) -
            1/2*sum(log(EVD_betaN$values[EVD_betaN$values > 0]))
        }
        
        ratio = exp(-numLog + denLog)
        
        phiHist[t] = phiHist[t-1]
        
        if(ratio > runif(1)){
          phiHist[t] = phiProp
          numLog = denLog
        }
        
        # Updating rho
        
        if(!is.null(Bds)){
          zz.temp = apply(zHist[t-1,,],2,which.max)
          out = MurraySchemePotts(beta.init = rhoHist[t-1], Coords = coords ,Nb = Bds ,n.colors = nComp, 
                                  col.obs = zz.temp, N.run = 100, range.beta = unlist(priors$rho[c('inf','sup')]))
          
          rhoHist[t] = mean(out$beta.sample, na.rm=T)
        }else{
          rhoHist[t] = runif(1,priors$rho$inf, priors$rho$sup)
        }
        # Updating z
        if(nSpTCov>0){MuG = Mu - SpTcov %*% gammaHist[t,]}else{MuG = Mu}
        
        if(!is.null(Bds)){vois = getCanStatF(Bds, zz.temp, nComp)}else{
          vois = matrix(0, nSites, nComp)
        }
        
        #if(t==2){browser()}
        
        #YpComp = crossprod_simple_triplet_matrix(t(Btilde),alphaHist[t,,])
        piCurr = matrix(piHist[t-1,,], nrow = nComp)
        if(is.null(Bds)){
          zObs = NULL
          logProb = matrix(NA, nSites, nComp)
          for(i in 1:nSites){
            idxID = which(y$ID == ID[i])
            zObsProb = NULL
            for(j in 1:nComp){
              #browser()
              YpComp = ft[[j]][idxID]# * betaHist[t,i*2] + betaHist[t,i*2-1]
              logProb[i,j] = sum(dnorm(Mu[idxID], YpComp, sigma2Hist[t], log=T), na.rm=T)* sqrt(2*pi*sigma2Hist[t]^2) + rhoHist[t] * vois[i,j]
            }
          }
        
        piCurr = matrix(piHist[t-1,,], nrow = nComp)
        
        logProb = logProb + t(log(piCurr))
        
        mlogProb = matrix(apply(logProb ,1, function(x) x - max(x, na.rm=T)), ncol=nComp, byrow=T)
        
        fProb = matrix(apply(mlogProb ,1, function(x) exp(x)/ sum(exp(x))), ncol=nComp, byrow=T)
        
        zHist[t,,] <- apply(fProb, 1, function(x) rmultinom(1,1,x))
      }else{
        zOLD <- zHist[t-1,,]
        piCurr = matrix(piHist[t-1,,], nrow = nComp)
        # Checquerboard updating
        logProb = matrix(NA, nSites, nComp)
        for(b in 1:length(Blocks)){
          for(i in Blocks[[b]]){
            idxID = which(y$ID == ID[i])
            zObsProb = NULL
            for(j in 1:nComp){
              YpComp = ft[[j]][idxID] * betaHist[t,i*2] + betaHist[t,i*2-1]
              logProb[i,j] = sum(dnorm(Mu[idxID], YpComp[idxID,j], sigma2Hist[t], log=T), na.rm=T) + rhoHist[t] * vois[i,j]
            }
            logProb[i,] = logProb[i,] + t(log(piCurr))[i,]
            mlogProb = matrix(apply(logProb,1, function(x) x - max(x, na.rm=T)), ncol=nComp, byrow=T)
            fProb = matrix(apply(mlogProb ,1, function(x) exp(x)/ sum(exp(x))), ncol=nComp, byrow=T)
            zHist[t,,i] <- rmultinom(1,1,fProb[i,])
            zOLD[,i] = zHist[t,,i]
          }
          zzO = apply(zOLD,2,which.max)
          vois = getCanStatF(Bds, zzO, nComp)
        }
      }
        # Updating pi
        
        tempAlpha=NULL
        for(i in 1:nSites){
          #browser()
          alpha = priors$pi$alpha0[i,] + zHist[t,,i]*length(which(y$ID == ID[i]))#zHist[t,,i]
          piHist[t,,i] = bayesm::rdirichlet(c(as.matrix(alpha)))
          tempAlpha = rbind(tempAlpha, alpha)
        }
        
        if(nSpTCov>0){
          MuA = MuB = Mu - SpTcov %*% gammaHist[t,]
          MuG = Mu - crossprod_simple_triplet_matrix(t(Xtilde), betaHist[t,])
        }else{
          MuA = MuB = Mu
          MuG = Mu - crossprod_simple_triplet_matrix(t(Xtilde), betaHist[t,])
        }
        
        # Likelihood and priors
        
        # Likelihood
        LogL = sum(dnorm(Ypp, Mu, sigma2Hist[t], log = T))
        
        # p(beta)
        betaPostMu = matrix(rep(priors$beta$m0, length(betaPostMu)), ncol=1)
        betaPostSd = diag(rep(priors$beta$s0, length(betaPostMu)))
        if(priors$beta$constrained){
          # Using Rue2005
          Ab = kronecker(zCurr, diag(nBasis))
          idxNC = which(rowSums(Ab)==0)
          if(length(idxNC)>0){
            Ab = matrix(Ab[-idxNC,], ncol = nSites*nBasis)
          }
          eb = matrix(c(rowSums(Ab)),ncol=1, nrow=nrow(Ab))
          eb[(1:(nrow(eb)/2))*2-1] <- 0
          
          # New cond dist parameters 
          muBetaC = betaPostMu - solve(betaPostSd) %*% t(Ab) %*% solve(Ab %*% solve(betaPostSd) %*% t(Ab)) %*% (Ab %*% betaPostMu - eb)
          sigBetaC = solve(betaPostSd) - solve(betaPostSd) %*% t(Ab) %*% solve(Ab %*% solve(betaPostSd) %*% t(Ab)) %*% Ab %*% solve(betaPostSd)
          EVD_beta = eigen(sigBetaC/2 + t(sigBetaC)/2)
          iEV_beta = 1 / EVD_beta$values
          iEV_beta[EVD_beta$values <= 0] = 0
          SigBM = EVD_beta$vectors %*% diag(iEV_beta) %*% t(EVD_beta$vectors)
          
          prBeta   =  -(length(betaHist[t,])-1)/2 * log(2*pi) -
            1/2*t(betaHist[t,]- muBetaC)%*% SigBM %*%  (betaHist[t,]- muBetaC) -
            1/2*sum(log(EVD_beta$values[EVD_beta$values > 1e-300]))
          
        }else{
          prBeta   = sum(mvtnorm::dmvnorm(betaHist[t,], betaPostMu, betaPostSd, log=T))
        }
        
        
        # p(gamma)
        prGamma = sum(dnorm(gammaHist[t,], priors$gamma$m0, priors$gamma$s0, log=T))
        
        # p(sigma), p(tau), p(phi), p(rho)
        prSigma2 = sum(dgamma(sigma2Hist[t], priors$sigma2$a0, priors$sigma2$b0, log=T))
        prTau = sum(dgamma(tauHist[t] ,priors$tau$a0, priors$tau$b0, log=T))
        prPhi = sum(log(dunif(phiHist[t], priors$phi$inf, priors$phi$sup)))
        prRho = sum(log(dunif(rhoHist[t], priors$rho$inf, priors$rho$sup)))
        #           
        
        # p(z | pi, rho), p(pi)
        tempAlpha[tempAlpha==0] = 1e-6
        prPi = prZ = 0
        for(i in 1:nSites){
          isZero = which(round(piHist[t,,i],5) > 0)
          prPi = prPi + 
            sum(log(
              mixtools::ddirichlet(matrix(round(piHist[t,isZero,i],6), ncol = length(isZero)) / sum(round(piHist[t,isZero,i],6)),
                                   alpha = tempAlpha[i,isZero] )
            ))
          prZ = prZ + log(fProb[i,which.max(zHist[t,,i])])
        }
        
        
        dpriorsHist[t, ] = cbind(prBeta, prGamma, prSigma2, prTau, prPhi, prZ, prRho, prPi)
        
        logPostDist[t] = LogL + prBeta + prGamma + prSigma2 + prTau + prPhi+ prPi + prZ + prRho
        logLike[t] = LogL
        
        if((round(t/2500) == (t/2500)) & print.res){print(t)}
      }
      
      betaH[,,nb] = betaHist
      gammaH[,,nb] = gammaHist 
      thetaGPH[,,nb] = thetaGPHist 
      sigma2H[,,nb] = sigma2Hist 
      tauH[,,nb] = tauHist
      phiH[,,nb] = phiHist
      rhoH[,,nb] = rhoHist
      piH[,,,nb] = piHist
      zH[,,,nb] = zHist
      logPostDistH[,,nb] = logPostDist
      logLikeH[,,nb] = logLike
      dpriorsH[,,nb] = dpriorsHist
      yHatH[,,nb] = yHatHist
      
      toc(log = TRUE, quiet = TRUE)
    }
    
    
    dimnames(betaH) <- list(1:N.run, 
                            paste(rep(ID, each =  nBasis), basis$tempPeriod),
                            1:nBatch)
    
    dimnames(dpriorsH) <- list(1:N.run, 
                               c("beta", "gamma", "sigma2", "tau", "phi", "z", "rho", "pi"),
                               1:nBatch)
    
    dimnames(gammaH) <- list(1:N.run, 
                             colnames(SpTcov),
                             1:nBatch)

    
    dimnames(piH) <- list(1:N.run, 
                          1:nComp,
                          ID,
                          1:nBatch)
    
    log.lst <- tic.log(format = FALSE)
    tic.clearlog()
    timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
    
    RET = list(theta = list(betaH = betaH, gammaH = gammaH, thetaGPH = thetaGPH,
                            sigma2H = sigma2H, tauH = tauH, phiH = phiH,rhoH = rhoH,
                            piH = piH, zH = zH),
               DistM = list(DM = DistMat, DMB = DistMatBeta),
               y = y, basis = basis, ID = ID ,coords = coords,cov = cov,sptcov = SpTcov,yHat = yHatH,
               priors = priors, kernelList = kernelList,
               execTime = timings, logPostDist = list(Post = logPostDistH, ll = logLikeH, dPriors = dpriorsH))
    
    return(RET)
    
  }