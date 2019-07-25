"gibbSamplerG" <- 
  function(priors = list(beta = list(m0 = 0, s0 = 1, dist = "gauss"),
                         alpha = list(m0 = 0, s0 = 1, dist = "gauss", constrained = FALSE),
                         gamma = list(m0 = 0, s0 = 1, dist = "gauss"),
                         sigma2 = list(a0=2, b0=1, dist = "gamma"),
                         tau = list(a0 = 2, b0= 1, dist = "gamma"),
                         phi = list(inf = 0.1, sup = 1, dist = "unif"),
                         pi = list(alpha0 = matrix(1,ncol=1), dist = "dirichlet"),
                         rho = list(inf=0.1, sup = 3, dist = "unif")),
           y, basis, ID,coords,cov,SpTcov,Bds = NULL,Blocks = NULL,
           N.run = 10000, nBatch = nBatch,
           debug = FALSE, print.res = FALSE,
           parallel = parallel, nCluster = nCluster, ...){
    
    
    siteRE = FALSE
    
    cov = cbind(cov, matrix(1, nrow = nrow(y), ncol = 1))
    
    if(basis$type == "penalized"){
      basisSplines = do.call(cbind,lapply(basis$splines, function(x) x$X))
    }else{
      basisSplines = do.call(cbind,basis$splines)
      basis2dSplines = do.call(cbind, basis$d2splines)
    }
    
    
    nBasis = length(basis$idx)
    nSites = length(ID)
    nSplines = ncol(basisSplines)
    nCov = ncol(cov)
    nObs = nrow(y)
    nComp = ncol(priors$pi$alpha0)
    
    if(is.null(SpTcov)){SpTcov = matrix(1, ncol = nCluster, nrow = nObs)}
    #if(is.null(SpTcov)){SpTcov = matrix(1, ncol = 1, nrow = nObs)}
    
    
    SpTcov = as.matrix(SpTcov)
    
    #idx.r = which(apply(SpTcov,2,function(x) diff(range(x)))==0)
    #if(length(idx.r) > 0){SpTcov = SpTcov[,-idx.r]}
    
    nSpTCov = ncol(SpTcov)
    
    alphaHist = array(NA, dim = c(N.run, nSplines,nComp))
    betaHist = matrix(NA, nrow = N.run, ncol = nSites * (nBasis) * nCov)
    betaIHist = matrix(NA, nrow = N.run, ncol = nSites)
    gammaHist = matrix(NA, nrow = N.run, ncol = nSpTCov)
    sigma2Hist = matrix(NA, nrow = N.run, ncol = 1)
    tauHist = matrix(NA, nrow = N.run, ncol = 1)
    phiHist = matrix(NA, nrow = N.run, ncol = 1)
    rhoHist = matrix(NA, nrow = N.run, ncol = 1)
    piHist = array(NA, dim = c(N.run, nComp, nSites))
    zHist = array(NA, dim = c(N.run, nComp, nSites))
    logPostDist = matrix(NA, nrow = N.run,ncol = 1)
    dpriorsHist = matrix(NA, nrow = N.run, ncol = 9)
    logLike = matrix(NA, nrow = N.run,ncol = 1)
    yHatHist = matrix(NA, nrow = N.run, ncol = nObs)
    
    alphaH = array(NA, dim = c(dim(alphaHist), nBatch))
    betaH = array(NA, dim = c(dim(betaHist), nBatch))
    betaIH = array(NA, dim = c(dim(betaIHist), nBatch))
    gammaH = array(NA, dim = c(dim(gammaHist), nBatch))
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
    
    #Xtilde = matrix(0,nrow = nObs, ncol = nSites*nCov*nBasis)
    Xtilde = simple_triplet_zero_matrix(nrow = nObs, ncol = nSites*nCov*nBasis)
    betaTilde = simple_triplet_zero_matrix(nrow = nrow(y), ncol = ncol(alphaHist))#matrix(0, nrow = nrow(y), ncol = ncol(alphaHist))
    
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
        tau0 = rigamma(1,priors$tau$a0, priors$tau$b0)
        phi0 = runif(1,priors$phi$inf, priors$phi$sup)
        rho0 = runif(1,priors$rho$inf, priors$rho$sup)
        if(debug){browser()}
        pi0 = matrix(round(t(apply(priors$pi$alpha0, 1, bayesm::rdirichlet)),5), ncol = nComp, byrow=F)
        z0 = matrix(apply(pi0*100, 1, function(x) rmultinom(1,1,x)), ncol=nComp, byrow=T)
        
        comp = apply(matrix(z0, nrow=nComp), 2 , which.max)
        
        if(priors$alpha$constrained){
          # Using Rue2005 alg.
          alphaPostMu = rep(priors$alpha$m0, nSplines)
          alphaPostSd = diag(priors$alpha$s0, nSplines)
          nConst = nBasis
          if(nBasis == 3){nConst = 2}
          A = matrix(0,nrow = nConst, ncol = length(alphaPostMu))
          for(cA in 1:(nConst)){
            #A[cA,basis$idx[[cA]]] = colMeans(basisSplines)[basis$idx[[cA]]]
            A[cA,basis$idx[[cA]]] = unlist(basis$isplines)[basis$idx[[cA]]]
          }
          
         # browser()
          
          e = matrix(c(0,0),ncol=1, nrow=(nConst))
          L = t(chol(solve(alphaPostSd)))
          zi = rmvnorm(1,rep(0, nrow(alphaPostSd)), diag(1, nrow(alphaPostSd)))
          v = solve(t(L),t(zi))
          xi = alphaPostMu + v
          Vnk = solve(alphaPostSd) %*% t(A)
          Wkk = A %*% Vnk
          Ukn = solve(Wkk) %*% t(Vnk)
          ci = A %*% xi - e
          alpha0 = rep(xi - t(Ukn) %*% ci, nComp)
        }else{
          alpha0 = rnorm(nSplines*nComp, priors$alpha$m0, priors$alpha$s0)
        }
        betaPostMu =rep(priors$beta$m0,nSites * (nBasis) * nCov)
        betaPostSd = kronecker(diag(nBasis),tau0*exp(- phi0 * DistMat)* outer(comp, comp, FUN = function(x,y) as.numeric(x==y)))
       #browser()
        if(priors$beta$constrained){
          # Using Rue2005 alg.
          Ab = kronecker(t(z0), diag(nBasis))
          idxNC = which(rowSums(Ab)==0)
          if(length(idxNC)>0){
            Ab = matrix(Ab[-idxNC,], ncol = nSites*nBasis)
          }
          eb = matrix(rowSums(Ab),ncol=1, nrow=nrow(Ab))
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
        betaI0 = rnorm(nSites, 0, 1)
        gamma0 = rnorm(nSpTCov,  priors$gamma$m0, priors$gamma$s0)
        sigma20 = rigamma(1,priors$sigma2$a0, priors$sigma2$b0)
        
        alphaHist[1,,] = alpha0
        betaHist[1,] = beta0
        betaIHist[1,] = betaI0
        gammaHist[1,] = gamma0
        sigma2Hist[1] = sigma20
        tauHist[1] = tau0
        rhoHist[1] = rho0
        phiHist[1] = phi0
        piHist[1,,] = t(pi0)
        zHist[1,,] = t(z0)
        
        betaS = matrix(betaHist[1,], ncol=nSites, byrow=F)
        
        for(nT in 1:(nBasis)){
          betaTilde[,basis$idx[[nT]]] = t(betaS[rep(nT, length(basis$idx[[nT]])),match(y$ID ,ID)])
        }
        
        Btilde = basisSplines * betaTilde
        
        if(nComp>1){SpTcov = z0[match(y$ID, ID),]}
        
        if(nSpTCov>0){
          MuA = MuB = Mu - SpTcov %*% gammaHist[1,]
          MuG = Mu - crossprod_simple_triplet_matrix(t(Xtilde), betaHist[1,])
        }else{
          MuA = MuB = Mu
          }
        
        Iter = matrix(diag(nSites)[match(y$ID, ID),], ncol = nSites)
        
        Yp0 = crossprod_simple_triplet_matrix(t(Xtilde), betaHist[1,]) + SpTcov %*% gammaHist[1,]
        Ypp0 = Yp0 + Iter %*% betaIHist[1,]
        
        numLog = NULL
        for(t in 2:N.run){
          
          
          #print(t)
          # 
          # Updating alpha_j
          
          zCurr = matrix(zHist[t-1,,], nrow=nComp)
          for(j in 1:nComp){
            idxZ = which(y$ID %in% ID[which(apply(zCurr,2,which.max) == j)])
            if(length(idxZ)>0){
              #BtildeZ = Btilde[idxZ,]
              BtildeZ = basisSplines[idxZ,]
              #tBB = crossprod_simple_triplet_matrix(BtildeZ)
              tBB = t(BtildeZ) %*% BtildeZ
              
              r = outer(basis$idx[[1]], basis$idx[[1]], FUN = function(X,Y){length(basis$idx[[1]])-abs(X-Y)})
              rr = outer(basis$idx[[1]], basis$idx[[1]], FUN = function(X,Y){abs(X-Y)})
              DistAlpha1 = matrix(mapply(min,r,rr), ncol = length(basis$idx[[1]]))^2
              if(nBasis==2){DistAlpha2 = outer(basis$idx[[2]], basis$idx[[2]], FUN = function(X,Y){abs(X-Y)})^2
              SigmaAlpha = rbind(cbind(exp(-DistAlpha1),matrix(0, ncol = length(basis$idx[[2]]), nrow = length(basis$idx[[1]]))), 
                          cbind(matrix(0, ncol = length(basis$idx[[1]]), nrow = length(basis$idx[[2]])), exp(-DistAlpha2)))}else{
                            SigmaAlpha =exp(-DistAlpha1)
                          }
              #browser()
              
              #iSigmaAlpha = 1/priors$alpha$s0 * diag(ncol(alphaHist))
              iSigmaAlpha = solve(priors$alpha$s0 *SigmaAlpha)
              
              #browser()
              
              #smPar = rgamma(1, 1e-8, 1)
              
              #print(t)
              
              #iSigmaAlpha = smPar*(t(basis2dSplines) %*% basis2dSplines)
              
              iSigmaObsB =  tBB
              #iSigmaObsB = tBB
              #print(t)
              #browser()
              alphaPostSd = sigma2Hist[t-1]*solve(iSigmaAlpha + iSigmaObsB)
              ialphaPostSd = sigma2Hist[t-1]^-1*(iSigmaAlpha + iSigmaObsB)
              #if(t > 450){browser()}
              #alphaPostMu = alphaPostSd %*% ( crossprod_simple_triplet_matrix(BtildeZ, MuA[idxZ]) + 
              #                                 iSigmaAlpha %*% matrix(priors$alpha$m0, ncol =1, nrow = nrow(iSigmaAlpha)))
              #  alphaPostMu = alphaPostSd %*% ( sigma2Hist[t-1]^(-1) * crossprod_simple_triplet_matrix(BtildeZ, MuA[idxZ]) + 
              #                                     iSigmaAlpha %*% matrix(priors$alpha$m0, ncol =1, nrow = nrow(iSigmaAlpha)))
               alphaPostMu = alphaPostSd %*% ( sigma2Hist[t-1]^(-1) *(t(BtildeZ) %*% MuA[idxZ]) +
                                                 iSigmaAlpha %*% matrix(priors$alpha$m0, ncol =1, nrow = nrow(iSigmaAlpha)))

            }else{
              alphaPostSd = diag(priors$alpha$s0, nSplines)
              ialphaPostSd = diag(priors$alpha$s0^-1, nSplines)
              alphaPostMu = matrix(priors$alpha$m0, nrow = nSplines, ncol=1)
            }
            if(priors$alpha$constrained){
              # Using Rue2005 alg.
              
              nConst = nBasis
              if(nBasis == 3){nConst = 2}
              
              A = matrix(0,nrow = nConst, ncol = length(alphaPostMu))
              for(cA in 1:(nConst)){
                #A[cA,basis$idx[[cA]]] = colMeans(basisSplines)[basis$idx[[cA]]]
                A[cA,basis$idx[[cA]]] = unlist(basis$isplines)[basis$idx[[cA]]]
                #A[2*cA,basis$idx[[cA]]] = 1#colMeans(basisSplines)[basis$idx[[cA]]]
              }
              #browser()
              e = matrix(c(0,0),ncol=1, nrow=(nConst))
              L = t(chol(ialphaPostSd))
              zi = rmvnorm(1,rep(0, nrow(alphaPostSd)), diag(1, nrow(alphaPostSd)))
              v = solve(t(L),t(zi))
              xi = alphaPostMu + v
              Vnk = ialphaPostSd %*% t(A)
              Wkk = A %*% Vnk
              Ukn = solve(Wkk) %*% t(Vnk)
              ci = A %*% xi - e
              alphaHist[t,,j] = xi - t(Ukn) %*% ci
              
              
              muAlphaC = alphaPostMu - ialphaPostSd %*% t(A) %*% solve(A %*% ialphaPostSd %*% t(A)) %*% (A %*% alphaPostMu - e)
              sigAlphaC = ialphaPostSd - ialphaPostSd %*% t(A) %*% solve(A %*% ialphaPostSd %*% t(A)) %*% A %*% ialphaPostSd
              #if(t > 250){browser()}
              #EVD_alpha = eigen(sigAlphaC)
              #iEV_alpha = 1 / EVD_alpha$values
              #iEV_alpha[EVD_alpha$values < 1e-300] = 0
              #SigAM = EVD_alpha$vectors %*% diag(iEV_alpha) %*% t(EVD_alpha$vectors)
              
            }else{
              alphaHist[t,,j] = rmvnorm(1, alphaPostMu, alphaPostSd)
            }
          }
          
          
          # Sample/Create latent temporal functionals
          for(si in 1:nSites){
            idxC = ((si-1)*nBasis+1):(si*nBasis)
            xt = NULL
            for(jj in 1:(nBasis)){
              fTime = as.matrix(basisSplines[,basis$idx[[jj]]]) %*% alphaHist[t,basis$idx[[jj]], which.max(zHist[t-1,,si])]
              xt = cbind(xt,fTime * cov)
            }
            #xt = cbind(xt, gammaHist[t-1,which.max(zHist[t-1,,si])])
            Xtilde[which(y$ID == ID[si]),idxC] = xt[which(y$ID == ID[si]),]
          }
          
          tXX = crossprod_simple_triplet_matrix(Xtilde)
          
          # Updating beta_j
          
          #browser()
          
          comp = apply(matrix(zHist[t-1,,], nrow=nComp), 2 , which.max)
          
          SigmaBeta = kronecker(diag(nBasis),tauHist[t-1]*exp(- phiHist[t-1] * DistMat)* outer(comp, comp, FUN = function(x,y) as.numeric(x==y)))
          iSigmaBeta = kronecker(diag(nBasis),1/tauHist[t-1]*solve(exp(- phiHist[t-1] * DistMat)* outer(comp, comp, FUN = function(x,y) as.numeric(x==y))))
          
          iSigmaObs = tXX
          #iSigmaObs = tXX
          
          betaPostSd = sigma2Hist[t-1]*solve(iSigmaBeta + iSigmaObs)
          ibetaPostSd = (iSigmaBeta + iSigmaObs)
          
          betaPostMu = betaPostSd %*% (sigma2Hist[t-1]^(-1) * crossprod_simple_triplet_matrix(Xtilde, MuB) + 
                                         iSigmaBeta %*% matrix(priors$beta$m0, ncol =1, nrow = nrow(iSigmaBeta)) )
          #betaPostMu = betaPostSd %*% (crossprod_simple_triplet_matrix(Xtilde, MuB) + 
          #                               iSigmaBeta %*% matrix(priors$beta$m0, ncol =1, nrow = nrow(iSigmaBeta)) )
          
          if(priors$beta$constrained){
            
            # Using Rue2005 alg.
            #Ab = kronecker(zCurr, matrix(c(0,1),ncol=2, nrow=1))
            Ab = kronecker(zCurr, diag(nBasis))
            idxNC = which(rowSums(Ab)==0)
            if(length(idxNC)>0){
              Ab = matrix(Ab[-idxNC,], ncol = nSites*nBasis)
            }
            eb = matrix(rowSums(Ab),ncol=1, nrow=nrow(Ab))
            L = t(chol(ibetaPostSd))
            zi = rmvnorm(1,rep(0, nrow(betaPostSd)), diag(1, nrow(betaPostSd)))
            v = solve(t(L),t(zi))
            xi = betaPostMu + v
            Vnk = (betaPostSd) %*% t(Ab)
            Wbkk = Ab %*% Vnk
            Ukn = solve(Wbkk) %*% t(Vnk)
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
          
          Yp = crossprod_simple_triplet_matrix(t(Xtilde), betaHist[t,])
          
          
          # Updating gamma_j
          
          betaS = matrix(betaHist[t,], ncol=nSites, byrow=F)
          # 
          for(nT in 1:(nBasis)){
            betaTilde[,basis$idx[[nT]]] = t(betaS[rep(nT, length(basis$idx[[nT]])),match(y$ID ,ID)])
          }
          # 
          Btilde = basisSplines * betaTilde
          
          #browser()
          if(t==2){browser()}
          if(nSpTCov>0){
            
            if(nComp>1){SpTcov = t(zCurr[,match(y$ID, ID)])}
            
            tWW = crossprod(SpTcov, SpTcov)
            iSigmaGamma = 1/priors$gamma$s0 * diag(ncol(gammaHist))
            
            iSigmaObsG = tWW
            #iSigmaObsG = tWW
            
            gammaPostSd = sigma2Hist[t-1]* solve(iSigmaGamma + iSigmaObsG)
            gammaPostMu = gammaPostSd %*% (sigma2Hist[t-1]^-1*crossprod(SpTcov, (MuG - Yp)) + 
                                             iSigmaGamma  %*% matrix(priors$gamma$m0, ncol =1, nrow = nrow(iSigmaGamma)))
            
            gammaHist[t,] = rmvnorm(1, gammaPostMu, gammaPostSd)
            
          }
          
          Yp = Yp + SpTcov %*% gammaHist[t-1,]
          
          # Site random effect
          
          Iter = matrix(diag(nSites)[match(y$ID, ID),], ncol = nSites)
          
          tII = crossprod(Iter, Iter)
          iSigmaBetaI = diag(nSites)
          
          iSigmaObsI = tII
          
          betaIPostSd = sigma2Hist[t-1]*solve(iSigmaBetaI + iSigmaObsI)
          ibetaIPostSd = iSigmaBetaI + iSigmaObsI
          
          betaIPostMu = betaIPostSd %*% (sigma2Hist[t-1]^-1*crossprod(Iter, (Mu - Yp)) + 
                                           iSigmaBetaI  %*% matrix(rep(0, nSites), ncol =1, nrow = nrow(iSigmaBetaI)))
          

          betaIHist[t,] = rmvnorm(1, betaIPostMu, betaIPostSd)
          
          # Updating sigma
          
          if(nComp > 0){Ypp = Yp + Iter %*% betaIHist[t,]
          }else{Ypp  = Yp}
          
          sigma2PostA = priors$sigma2$a0 + nObs / 2
          sigma2PostB = priors$sigma2$b0 + sum((Mu - Ypp)^2, na.rm=T) / 2 #+ nObs * priors$beta$s0 / (nObs + priors$beta$s0) * mean((Mu - Ypp0)^2)/2
          
          sigma2Hist[t] = rigamma(1, sigma2PostA, sigma2PostB)
         # if(sigma2Hist[t] > 1){browser()}
          
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
          if(!priors$beta$constrained){
            tauPostB = priors$tau$b0 + 
              (t(as.matrix(betaHist[t,]) - betaMu) %*% iSigmaBT %*% as.matrix(betaHist[t,]- betaMu))/2
          }
          
         # browser()
          
          if(priors$beta$constrained){
            #browser()
            
            idxNN = which(rowSums(zCurr)>0)
            
            Ab = kronecker(zCurr[idxNN,], diag(nBasis))
            if(length(idxNN)==1){
             Ab = t(Ab) 
            }
            if(nComp == 1){
              Ab = kronecker(zCurr, diag(nBasis))
            }
            betaMu = muBeta0
            nConst = nrow(Ab)
            sigBetaC = SigmaBT - SigmaBT %*% t(Ab) %*% solve((Ab) %*% SigmaBT %*% t(Ab)) %*% (Ab) %*% SigmaBT
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
          
          SigmaBetaN = kronecker(diag(nBasis),tauHist[t-1]*exp(- phiProp * DistMat))
          iSigmaBetaN = kronecker(diag(nBasis),1/tauHist[t-1]*solve(exp(- phiProp * DistMat)* outer(comp, comp, FUN = function(x,y) as.numeric(x==y))))
          #iSigmaBeta = kronecker(diag(nBasis),1/tauHist[t-1]*solve(exp(- phiHist[t-1] * DistMat)* outer(comp, comp, FUN = function(x,y) as.numeric(x==y))))
          betaPostSdN = solve(iSigmaBetaN + iSigmaObs)
          ibetaPostSdN = iSigmaBetaN + iSigmaObs
          
          # betaPostMuN = betaPostSdN %*% (sigma2Hist[t-1]^(-1) * crossprod_simple_triplet_matrix(Xtilde, MuG) + 
          #                                iSigmaBetaN %*% matrix(priors$beta$m0, ncol =1, nrow = nrow(iSigmaBetaN)) )
          if(priors$beta$constrained){
            sigBetaCN = ibetaPostSdN - ibetaPostSdN %*% t(Ab) %*% solve((Ab) %*% ibetaPostSdN %*% t(Ab)) %*% (Ab) %*% ibetaPostSdN
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

          if(!is.null(Bds)){vois = getCanStatF(Bds, zz.temp, nComp)}else{
            vois = matrix(0, nSites, nComp)
          }
          
          YpComp = crossprod_simple_triplet_matrix(t(Btilde),alphaHist[t,,])
          
          #browser()
          
          if(nSpTCov == nComp){
            YpComp = YpComp + matrix(gammaHist[t-1,], ncol = nSpTCov, nrow = nrow(YpComp), byrow = T)
          }else{
            YpComp = YpComp + SpTcov %*% gammaHist[t-1,]
          }
          
          #YpComp = YpComp #+ matrix(Iter %*% betaIHist[t,], ncol = nComp, nrow = nrow(YpComp), byrow = F)
          
          piCurr = matrix(piHist[t-1,,], nrow = nComp)
          
   
          if(is.null(Bds)){
          zObs = NULL
          logProb = matrix(NA, nSites, nComp)
          for(i in 1:nSites){
            idxID = which(y$ID == ID[i])
            zObsProb = NULL
            for(j in 1:nComp){
              logProb[i,j] = sum(dnorm(Mu[idxID], YpComp[idxID,j], sigma2Hist[t], log=T), na.rm=T) + rhoHist[t] * vois[i,j]
            }
          }
          
          piCurr = matrix(piHist[t-1,,], nrow = nComp)
          
          #if(min(rowSums(zCurr))==0){browser()}
          
          #if(t == 400){browser()}
          
          logProb = logProb + t(log(piCurr))
          
          mlogProb = matrix(apply(logProb ,1, function(x) x - max(x, na.rm=T)), ncol=nComp, byrow=T)
          
          fProb = matrix(apply(mlogProb ,1, function(x) exp(x)/ sum(exp(x))), ncol=nComp, byrow=T)
          
          z.Num = apply(fProb,1, function(x) sample(1:nComp,1,prob = x))
          
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
            alpha = priors$pi$alpha0[i,] +  zHist[t,,i]*length(which(y$ID == ID[i])) #zObs[i,]#
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
          
          # p(alpha)
          prAlpha = 0
          for(j in 1:nComp){
            alphaPostMu = matrix(rep(priors$alpha$m0, length(alphaPostMu)), ncol=1)
            alphaPostSd = diag(rep(priors$alpha$s0, length(alphaPostMu)))
            if(priors$alpha$constrained){
              # Using Rue2005
              nConst = 1 #nBasis
              A = matrix(0,nrow = nConst, ncol = length(alphaPostMu))
              for(cA in 1:(nConst)){
                A[cA,basis$idx[[cA]]] = colMeans(basisSplines)[basis$idx[[cA]]]
              }
              e = matrix(0,ncol=1, nrow=nConst)
              
              muAlphaC = alphaPostMu - solve(alphaPostSd) %*% t(A) %*% solve(A %*% solve(alphaPostSd) %*% t(A)) %*% (A %*% alphaPostMu - e)
              sigAlphaC = solve(alphaPostSd) - solve(alphaPostSd) %*% t(A) %*% solve(A %*% solve(alphaPostSd) %*% t(A)) %*% A %*% solve(alphaPostSd)
              
              EVD = eigen(sigAlphaC)
              iEV = 1 / EVD$values
              iEV[EVD$values < 1e-300] = 0
              SigM = EVD$vectors %*% diag(iEV) %*% t(EVD$vectors)
              
              prAlphaT   =  -(length(alphaHist[t,,j])-1)/2 * log(2*pi) -
                1/2*t(alphaHist[t,,j]- muAlphaC)%*% SigM %*%  (alphaHist[t,,j]- muAlphaC) -
                1/2*sum(log(EVD$values[EVD$values > 1e-300]))
              
            }else{
              prAlphaT   = sum(mvtnorm::dmvnorm(alphaHist[t,,j], alphaPostMu, alphaPostSd, log=T))
            }
            prAlpha = prAlpha + prAlphaT
          }
          
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
            eb = matrix(rowSums(Ab),ncol=1, nrow=nrow(Ab))
            #browser()
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
          prSigma2 = sum(dgamma(sigma2Hist[t], priors$sigma2$b0, priors$sigma2$a0, log=T))
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
          
          logPostDist[t] = LogL + prBeta + prGamma + prAlpha + 
            prSigma2 + prTau + prPhi+ prPi + prZ + prRho + prPi
          logLike[t] = LogL
          
          dpriorsHist[t, ] = cbind(prBeta, prGamma, prAlpha, prSigma2, prTau, prPhi, prZ, prRho, prPi)
          
          if((round(t/250) == (t/250)) & print.res){print(t)}
        }
        
        alphaH[,,,nb] = alphaHist
        betaH[,,nb] = betaHist
        betaIH[,,nb] = betaIHist
        gammaH[,,nb] = gammaHist 
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
      
      dimnames(alphaH) <- list(1:N.run, 
                               unlist(lapply(basis$idx, function(x) 1:length(x))),
                              1:nComp,
                              1:nBatch)
      
      dimnames(betaH) <- list(1:N.run, 
                               paste(rep(ID, each =  nBasis), basis$tempPeriod),
                               1:nBatch)
      
      dimnames(betaIH) <- list(1:N.run, 
                               ID,
                               1:nBatch)
      
      dimnames(dpriorsH) <- list(1:N.run, 
                                 c("beta", "gamma", "alpha", "sigma2", "tau", "phi", "z", "rho", "pi"),
                                 1:nBatch)
      
      dimnames(gammaH) <- list(1:N.run, 
                                 colnames(SpTcov),
                                 1:nBatch)
      
      # dimnames(zH) <- list(1:N.run, 
      #                      1:nComp,
      #                       ID,
      #                      1:nBatch)
      
      dimnames(piH) <- list(1:N.run, 
                           1:nComp,
                           ID,
                           1:nBatch)
      
      log.lst <- tic.log(format = FALSE)
      tic.clearlog()
      timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
    
    RET = list(theta = list(alphaH = alphaH,betaH = betaH, betaIH = betaIH, gammaH = gammaH, 
                            sigma2H = sigma2H, tauH = tauH, phiH = phiH,rhoH = rhoH,
                            piH = piH, zH = zH),
               DistM = list(DM = DistMat, DMB = DistMatBeta),
               y = y, basis = basis, ID = ID ,coords = coords,cov = cov,sptcov = SpTcov, yHat = yHatH,
               priors = priors, execTime = timings, logPostDist = list(Post = logPostDistH, ll = logLikeH, dPriors = dpriorsH))
    
    return(RET)
    
  }