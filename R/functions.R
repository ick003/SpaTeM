"SPTMData" <- 
  function(df.obs, tempBasis = "re", tempPeriod = "%m", nSplines = 3){
  
  SiteID = levels(df.obs$ID)
  
  t = df.obs$date
  
  BB = NULL
  idxT = NULL
  j=0
  maxIdx = 0
  for(i in tempPeriod){
    j = j +1 
    if(tempBasis == "bs"){
      if(j == 1){
        tTemp = as.Date(paste("2014-",format(t,"%m-%d"), sep = ""))
        BBT = bs(tTemp, df = nSplines[j],
                 Boundary.knots = c(min(tTemp)-diff(range(tTemp))/(1*nSplines[j]),max(tTemp)+diff(range(tTemp))/(1*nSplines[j])))
      }
      if(j == 2){
        BBT = bs(t, df = nSplines[j], 
                 Boundary.knots = c(min(t)-diff(range(t))/(1*nSplines[j]),max(t)+diff(range(t))/(1*nSplines[j])))
      }
      if(j == 3){
        BBT = matrix(rep(1, length(t)), ncol=1)
      }
      }
    if(tempBasis == "re"){
      BBT = sapply(BbaseTr, function(x) as.numeric(as.numeric(format(t,i)) == x))
    }
    
    BB = c(BB,list(BBT))
    idxt = 1:ncol(BBT) + maxIdx #(ncol(BB) - ncol(BBT))
    maxIdx = max(idxt)
    idxT = c(idxT,list(idxt))
}
  
  sptm.data = list(obs.data = df.obs, basis = BB, list.idx = idxT, tempBasis = tempBasis, tempPeriod = tempPeriod)
  class(sptm.data) <- 'sptm'
  return(sptm.data)
}
"SPTModel" <- 
  function(df.sptm, df.lu = NULL, coordinates, cov = NULL, SpTcov = NULL){
  
  
  sptm.model = list(data = df.sptm$obs.data,
                    basis = list(splines = df.sptm$basis, idx = df.sptm$list.idx, tempBasis = df.sptm$tempBasis, tempPeriod = df.sptm$tempPeriod),
                    lu = df.lu,
                    coord = coordinates,
                    covariates = cov,
                    SpTcov = SpTcov)
  
  class(sptm.model) <- 'stpmod'
  return(sptm.model)
  
}

"estimGibbs" <- 
  function(df.sptmod, model = 'noMixture', tempRE = "corr",
         priors = list(beta = list(m0 = 0, s0 = 1, dist = "gauss"),
                       alpha = list(m0 = 0, s0 = 1, dist = "gauss", constrained = FALSE),
                       thetaGP = list(k1 = list(a1 = 20, b1=10, a2=10, b2 = 500),k2=list(NULL)),
                       gamma = list(m0 = 0, s0 = 1, dist = "gauss"),
                       sigma2 = list(a0=2, b0=1, dist = "gamma"),
                       tau = list(a0 = 2, b0= 1, dist = "gamma"),
                       phi = list(inf = 0.1, sup = 1, dist = "unif"),
                       tauT = list(a0 = 2, b0= 1, dist = "gamma"),
                       phiT = list(inf = 0.1, sup = 1, dist = "unif"),
                       pi = list(alpha0 = matrix(1,nrow=1,ncol=1), dist = "dirichlet"),
                       rho = list(inf=0.1, sup=3, dist="unif")),
         mixture = TRUE, N.run = 10000,
         print.res = FALSE, debug=FALSE, meanFunc = NULL, Xcov = NULL,
         nBatch = 1, parallel = F, nCluster = NULL, kernelList = NULL){
  

  p.comp = df.sptmod$lu
  
  basis = df.sptmod$basis
  
  y = df.sptmod$data
  
  coords =  df.sptmod$coord
  
  cov = df.sptmod$covariates
  
  SpTcov = df.sptmod$SpTcov
  
  if(model == 'noMixture'){
    
    if(tempRE == 'notcorr'){
      res = gibbSamplerG(priors, y, basis = basis, ID = levels(y$ID), coords = coords, cov = cov,SpTcov = SpTcov,Bds = NULL,
                         N.run = N.run, debug, print.res = print.res, nBatch = nBatch,
                         parallel = parallel, nCluster = nCluster)
    }
    
    if(tempRE == 'corr'){
      res = gibbSamplerS2(priors, y, basis = basis, ID = levels(y$ID), coords = coords, cov = cov,
                          N.run = N.run, debug, print.res = print.res,
                          parallel = parallel, nCluster = nCluster)
    }
    
    
    RET = list(GibbsOut = res, model = model, tempRE = tempRE)
    
  }
  
  if(model == 'simpleMixture'){
    if(tempRE == 'notcorr'){
      res = gibbSamplerG(priors, y, basis = basis, ID = levels(y$ID), coords = coords, cov = cov,SpTcov = SpTcov,Bds = NULL,
                         N.run = N.run, debug, print.res = print.res, nBatch = nBatch,
                         parallel = parallel, nCluster = nCluster)
    }
    if(tempRE == 'gp'){
      res = gibbSamplerGP(priors, y, basis = basis, ID = levels(y$ID), coords = coords, cov = cov,SpTcov = SpTcov,Bds = NULL,
                         N.run = N.run, debug, print.res = print.res, nBatch = nBatch,
                         parallel = parallel, nCluster = nCluster, kernelList = kernelList)
    }
    if(tempRE == 'corr'){
      res = gibbSamplerM2(priors, y, basis = basis, ID = levels(y$ID), coords = coords, cov = cov,
                          N.run = N.run, debug, print.res = print.res,
                          parallel = parallel, nCluster = nCluster)
    }
    RET = list(GibbsOut = res, model = model, tempRE = tempRE)
    
  }
  
  if(model == 'spatialMixture'){
    
    if(debug){browser()}
    
    dd = deldir(df.sptmod$coord[,2], df.sptmod$coord[,3],suppressMsge=TRUE)
    tile.dd = tile.list(dd)
    
    Bds = comPts(tile.dd)
    
    res = gibbSamplerG(priors, y, basis = basis, ID = levels(y$ID), coords = coords, cov = cov,SpTcov = SpTcov,Bds= Bds,
                       N.run = N.run, debug, print.res = print.res, nBatch = nBatch,
                       parallel = parallel, nCluster = nCluster)
    
    # res = gibbSamplerMSpat(priors, y, basis = basis,ID = levels(y$ID), coords = coords,
    #                        cov = cov,Bds=Bds, N.run = N.run, nBatch,
    #                debug, print.res = print.res,
    #                parallel = parallel, nCluster = nCluster)
    # 
    RET = list(GibbsOut = res, model = model, tempRE = tempRE, Voronoi = list(dd, tile.dd, Bds))
  }
  
  class(RET) <- 'sptmRes'
  return(RET)
}

"gibbSamplerG" <- 
  function(priors = list(beta = list(m0 = 0, s0 = 1, dist = "gauss"),
                         alpha = list(m0 = 0, s0 = 1, dist = "gauss", constrained = FALSE),
                         gamma = list(m0 = 0, s0 = 1, dist = "gauss"),
                         sigma2 = list(a0=2, b0=1, dist = "gamma"),
                         tau = list(a0 = 2, b0= 1, dist = "gamma"),
                         phi = list(inf = 0.1, sup = 1, dist = "unif"),
                         pi = list(alpha0 = matrix(1,ncol=1), dist = "dirichlet"),
                         rho = list(inf=0.1, sup = 3, dist = "unif")),
           y, basis, ID,coords,cov,SpTcov,Bds = NULL,
           N.run = 10000, nBatch = nBatch,
           debug = FALSE, print.res = FALSE,
           parallel = parallel, nCluster = nCluster, ...){
    
    
    siteRE = FALSE
    
    cov = cbind(cov, matrix(1, nrow = nrow(y), ncol = 1))
    basisSplines = do.call(cbind,basis$splines)
    
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
    
    if(parallel){
      
      #cl <- makeCluster(nCluster) # create a cluster with nCores cores
      #registerDoParallel(cl)
      #registerDoSNOW(cl)
      #registerDoMC(nCluster)
      parallelRes = foreach(nb = 1:nBatch,
                            .combine = "c",
                            .packages = c("mvtnorm", "slam", "splines", "MASS","LearnBayes", "tictoc")) %dopar% {
                              
                              tic.clearlog()
                              tic()
                              if(priors$alpha$constrained){
                                # Using Rue2005 alg.
                                alphaPostMu = rep(priors$alpha$m0, nSplines)
                                alphaPostSd = diag(priors$alpha$s0, nSplines)
                                nConst = nBasis
                                if(nBasis == 2){nConst = 1}
                                A = matrix(0,nrow = nConst, ncol = length(alphaPostMu))
                                for(cA in 1:(nConst)){
                                  A[cA,basis$idx[[cA]]] = colMeans(basisSplines)[basis$idx[[cA]]]
                                }
                                e = matrix(0,ncol=1, nrow=(nConst))
                                L = t(chol(solve(alphaPostSd)))
                                zi = rmvnorm(1,rep(0, nrow(alphaPostSd)), diag(1, nrow(alphaPostSd)))
                                v = solve(t(L),t(zi))
                                xi = alphaPostMu + v
                                Vnk = solve(alphaPostSd) %*% t(A)
                                Wkk = A %*% Vnk
                                Ukn = solve(Wkk) %*% t(Vnk)
                                ci = A %*% xi - e
                                alpha0 = xi - t(Ukn) %*% ci
                              }else{
                                alpha0 = rnorm(nSplines, priors$alpha$m0, priors$alpha$s0)
                              }
                              betaPostMu =rep(priors$beta$m0,nSites * nBasis * nCov)
                              betaPostSd = diag(priors$beta$s0,nSites * nBasis * nCov)
                              if(priors$beta$constrained){
                                # Using Rue2005 alg.
                                A = matrix(0,nrow = 1, ncol = length(betaPostMu))
                                for(cA in 2:nBasis){
                                  A[cA-1,(0:(nSites-1))*nBasis + cA] = 1
                                }
                                e = matrix(1,ncol=1, nrow=1)
                                #e[2] = 1
                                L = t(chol(solve(betaPostSd)))
                                zi = rmvnorm(1,rep(0, nrow(betaPostSd)), diag(1, nrow(betaPostSd)))
                                v = solve(t(L),t(zi))
                                xi = betaPostMu + v
                                Vnk = solve(betaPostSd) %*% t(A)
                                Wkk = A %*% Vnk
                                Ukn = solve(Wkk) %*% t(Vnk)
                                ci = A %*% xi - e
                                beta0 = xi - t(Ukn) %*% ci
                              }else{
                                beta0 = rmvnorm(1, betaPostMu, betaPostSd)
                              }
                              gamma0 = rnorm(nSpTCov,  priors$gamma$m0, priors$gamma$s0)
                              sigma20 = rigamma(1,priors$sigma2$a0, priors$sigma2$b0)
                              tau0 = rigamma(1,priors$tau$a0, priors$tau$b0)
                              phi0 = runif(1,priors$phi$inf, priors$phi$sup)
                              
                              alphaHist[1,] = alpha0
                              betaHist[1,] = beta0
                              gammaHist[1,] = gamma0
                              sigma2Hist[1] = sigma20
                              tauHist[1] = tau0
                              phiHist[1] = phi0
                              
                              if(nSpTCov>0){MuG = Mu - SpTcov %*% gammaHist[1,]}else{MuG = Mu}
                              
                              numLog = NULL
                              for(t in 2:N.run){
                                
                                # Updating alpha_j
                                
                                betaS = matrix(betaHist[t-1,], ncol=nSites, byrow=F)
                                
                                for(nT in 1:nBasis){
                                  betaTilde[,basis$idx[[nT]]] = t(betaS[rep(nT, length(basis$idx[[nT]])),match(y$ID ,ID)])
                                }
                                
                                Btilde = basisSplines * betaTilde
                                
                                tBB = crossprod_simple_triplet_matrix(Btilde)
                                itBB = solve(tBB)
                                
                                SigmaAlpha = 1 / priors$alpha$s0 * diag(ncol(alphaHist))
                                iSigmaAlpha = priors$alpha$s0 * diag(ncol(alphaHist))
                                
                                SigmaObsB = sigma2Hist[t-1]* itBB
                                iSigmaObsB = 1/sigma2Hist[t-1]* tBB
                                
                                alphaPostSd = solve(iSigmaAlpha + iSigmaObsB)
                                alphaPostMu = alphaPostSd %*% (iSigmaObsB %*% (itBB %*% crossprod_simple_triplet_matrix(Btilde, MuG)))
                                
                                if(priors$alpha$constrained){
                                  # Using Rue2005 alg.
                                  nConst = 1#nBasis
                                  if(nBasis == 3){nConst = 1}
                                  A = matrix(0,nrow = nConst, ncol = length(alphaPostMu))
                                  for(cA in 1:(nConst)){
                                    A[cA,basis$idx[[cA]]] = colMeans(basisSplines)[basis$idx[[cA]]]
                                  }
                                  e = matrix(0,ncol=1, nrow=nConst)
                                  #e[2] = mean(MuG)
                                  L = t(chol(solve(alphaPostSd)))
                                  zi = rmvnorm(1,rep(0, nrow(alphaPostSd)), diag(1, nrow(alphaPostSd)))
                                  v = solve(t(L),t(zi))
                                  xi = alphaPostMu + v
                                  Vnk = solve(alphaPostSd) %*% t(A)
                                  Wkk = A %*% Vnk
                                  Ukn = solve(Wkk) %*% t(Vnk)
                                  ci = A %*% xi - e
                                  alphaHist[t,] = xi - t(Ukn) %*% ci
                                }else{
                                  alphaHist[t,] = rmvnorm(1, alphaPostMu, alphaPostSd)
                                }
                                # Update temporal latent functionals
                                
                                xt = NULL
                                for(jj in 1:nBasis){
                                  fTime = as.matrix(basisSplines[,basis$idx[[jj]]]) %*% alphaHist[t,basis$idx[[jj]]]
                                  xt = cbind(xt,fTime * cov)
                                }
                                
                                for(si in 1:nSites){
                                  idxC = ((si-1)*nBasis+1):(si*nBasis)
                                  Xtilde[idxID[[si]],idxC] = xt[idxID[[si]],]
                                }
                                
                                tXX = crossprod_simple_triplet_matrix(Xtilde)
                                
                                itXX = solve(tXX)
                                
                                # Updating beta_j
                                
                                SigmaBeta = kronecker(diag(nBasis),tauHist[t-1]*exp(- phiHist[t-1] * DistMat))
                                iSigmaBeta = kronecker(diag(nBasis),chol2inv(chol(tauHist[t-1]*exp(- phiHist[t-1] * DistMat))))
                                
                                SigmaObs = sigma2Hist[t-1] * itXX
                                iSigmaObs = 1/sigma2Hist[t-1]* tXX
                                
                                betaPostSd = solve(iSigmaBeta + iSigmaObs)
                                
                                betaPostMu = betaPostSd %*% (iSigmaObs %*% (itXX %*% crossprod_simple_triplet_matrix(Xtilde, MuG)))
                                
                                if(priors$beta$constrained){
                                  # Using Rue2005 alg.
                                  A = matrix(0,nrow = 1, ncol = length(betaPostMu))
                                  for(cA in 2:nBasis){
                                    A[cA-1,(0:(nSites-1))*nBasis + cA] = 1
                                  }
                                  e = matrix(1,ncol=1, nrow=1)
                                  #e[2] = 1
                                  L = t(chol(solve(betaPostSd)))
                                  zi = rmvnorm(1,rep(0, nrow(betaPostSd)), diag(1, nrow(betaPostSd)))
                                  v = solve(t(L),t(zi))
                                  xi = betaPostMu + v
                                  Vnk = solve(betaPostSd) %*% t(A)
                                  Wkk = A %*% Vnk
                                  Ukn = solve(Wkk) %*% t(Vnk)
                                  ci = A %*% xi - e
                                  betaHist[t,] = xi - t(Ukn) %*% ci
                                }else{
                                  betaHist[t,] = rmvnorm(1, betaPostMu, betaPostSd)
                                }
                                
                                # Updating gamma_j
                                
                                #Yp = matprod_simple_triplet_matrix(Btilde,alphaHist[t,])
                                Yp = crossprod_simple_triplet_matrix(t(Xtilde), betaHist[t,])
                                if(nSpTCov>0){
                                  
                                  tWW = crossprod(SpTcov, SpTcov)
                                  itWW = chol2inv(chol(tWW))
                                  
                                  SigmaGamma = sigma2Hist[t-1] * itWW
                                  iSigmaGamma = priors$gamma$s0 * diag(ncol(gammaHist))
                                  
                                  SigmaObsG = sigma2Hist[t-1]* itWW
                                  iSigmaObsG = 1/sigma2Hist[t-1]* tWW
                                  
                                  gammaPostSd =solve(iSigmaGamma + iSigmaObsG)
                                  gammaPostMu = gammaPostSd %*% (iSigmaObsG %*% (itWW %*% crossprod(SpTcov, (Mu - Yp))))
                                  
                                  gammaHist[t,] = rmvnorm(1, gammaPostMu, gammaPostSd)
                                  
                                }
                                
                                if(nSpTCov>0){MuG = Mu - SpTcov %*% gammaHist[t,]}else{MuG = Mu}
                                
                                # Updating sigma
                                
                                if(nSpTCov>0){Ypp = Yp + SpTcov %*% gammaHist[t,]}else{Ypp = Yp}
                                
                                sigma2PostA = priors$sigma2$a0 + nObs / 2
                                sigma2PostB = priors$sigma2$b0 + sum((Mu - Ypp)^2, na.rm=T) / 2
                                
                                sigma2Hist[t] = rigamma(1, sigma2PostA, sigma2PostB)
                                
                                # Updating tau
                                
                                iSigmaBT = kronecker(diag(nBasis),solve(exp(- phiHist[t-1] * DistMat)))
                                iSigmaBT = (iSigmaBT + t(iSigmaBT)) / 2
                                
                                tauPostA = priors$tau$a0 + nSites*nBasis / 2
                                tauPostB = priors$tau$b0 + t(as.matrix(betaHist[t,])) %*% iSigmaBT %*% as.matrix(betaHist[t,])/2
                                
                                tauHist[t] = rigamma(1,tauPostA, tauPostB)
                                
                                # Updating phi (MH Reject)
                                
                                phiProp = runif(1, priors$phi$inf, priors$phi$sup)
                                
                                SigmaBetaN = kronecker(diag(nBasis),tauHist[t-1]*exp(- phiProp * DistMat))
                                betaPostSdN = solve( solve(SigmaBetaN) + iSigmaObs)
                                betaPostMuN = betaPostSdN %*% (iSigmaObs %*% (itXX %*% crossprod_simple_triplet_matrix(Xtilde, MuG)))
                                
                                numLog = -1/2 * log(det(betaPostSd)) +  -1/2*t(betaHist[t,] - betaPostMu) %*% solve(betaPostSd) %*% (betaHist[t,] - betaPostMu)
                                
                                denLog = -1/2 * log(det(betaPostSdN)) +  -1/2*t(betaHist[t,] - betaPostMuN) %*% solve(betaPostSdN) %*% (betaHist[t,] - betaPostMuN)
                                
                                ratio = exp(-numLog + denLog)
                                
                                phiHist[t] = phiHist[t-1]
                                
                                if(ratio > runif(1)){
                                  phiHist[t] = phiProp
                                  numLog = denLog
                                }
                                # Likelihood and priors
                                LogL = sum(dnorm(Ypp, Mu, sigma2Hist[t], log = T))
                                alphaPostMu = matrix(rep(priors$alpha$m0, length(alphaPostMu)), ncol=1)
                                alphaPostSd = diag(rep(priors$alpha$s0, length(alphaPostMu)))
                                if(priors$alpha$constrained){
                                  # Using Rue2005 
                                  nConst = 1#nBasis
                                  if(nBasis == 3){nConst = 1}
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
                                  
                                  prAlpha   =  -(length(alphaHist[t,])-1)/2 * log(2*pi) -
                                    1/2*t(alphaHist[t,]- muAlphaC)%*% SigM %*%  (alphaHist[t,]- muAlphaC) -
                                    1/2*sum(log(EVD$values[EVD$values > 1e-300]))
                                  
                                }else{
                                  prAlpha   = sum(mvtnorm::dmvnorm(alphaHist[t,], alphaPostMu, alphaPostSd, log=T))
                                }
                                betaPostMu = matrix(rep(priors$beta$m0, length(betaPostMu)), ncol=1)
                                betaPostSd = diag(rep(priors$beta$s0, length(betaPostMu)))
                                if(priors$beta$constrained){
                                  # Using Rue2005 
                                  A = matrix(0,nrow = 1, ncol = length(betaPostMu))
                                  for(cA in 2:nBasis){
                                    A[cA-1,(0:(nSites-1))*nBasis + cA] = 1
                                  }
                                  e = matrix(1,ncol=1, nrow=1)
                                  
                                  muBetaC = betaPostMu - solve(betaPostSd) %*% t(A) %*% solve(A %*% solve(betaPostSd) %*% t(A)) %*% (A %*% betaPostMu - e)
                                  sigBetaC = solve(betaPostSd) - solve(betaPostSd) %*% t(A) %*% solve(A %*% solve(betaPostSd) %*% t(A)) %*% A %*% solve(betaPostSd)
                                  
                                  EVD = eigen(sigBetaC)
                                  iEV = 1 / EVD$values
                                  iEV[EVD$values < 1e-300] = 0
                                  SigM = EVD$vectors %*% diag(iEV) %*% t(EVD$vectors)
                                  
                                  prBeta   =  -(length(betaHist[t,])-1)/2 * log(2*pi) -
                                    1/2*t(betaHist[t,]- muBetaC)%*% SigM %*%  (betaHist[t,]- muBetaC) -
                                    1/2*sum(log(EVD$values[EVD$values > 1e-300]))
                                  
                                }else{
                                  prBeta   = sum(mvtnorm::dmvnorm(betaHist[t,], betaPostMu, betaPostSd, log=T))
                                }
                                prGamma = sum(dnorm(gammaHist[t,], priors$gamma$m0, priors$gamma$s0, log=T))
                                prSigma2 = sum(dgamma(sigma2Hist[t], priors$sigma2$a0, priors$sigma2$b0, log=T))
                                prTau = sum(dgamma(tauHist[t] ,priors$tau$a0, priors$tau$b0, log=T))
                                prPhi = sum(dunif(phiHist[t], priors$phi$inf, priors$phi$sup))
                                
                                logPostDist[t] = LogL + prAlpha + prBeta + prGamma + prSigma2 + prTau + prPhi
                                
                                
                                if((round(t/100) == (t/100)) & print.res){print(t)}
                              }
                              
                              toc(log = TRUE, quiet = TRUE)
                              log.lst <- tic.log(format = FALSE)
                              tic.clearlog()
                              timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
                              
                              list(alphaHist,
                                   betaHist,
                                   gammaHist, 
                                   sigma2Hist, 
                                   tauHist,
                                   phiHist,
                                   logPostDist,
                                   timings)
                              
                              
                            }
      stopCluster(cl)
      
      timings = NULL
      for(nb in 1:nBatch){
        alphaH[,,nb] = parallelRes[[1+8*(nb-1)]]
        betaH[,,nb] = parallelRes[[2+8*(nb-1)]]
        gammaH[,,nb] = parallelRes[[3+8*(nb-1)]]
        sigma2H[,,nb] = parallelRes[[4+8*(nb-1)]]
        tauH[,,nb] = parallelRes[[5+8*(nb-1)]]
        phiH[,,nb] = parallelRes[[6+8*(nb-1)]]
        logPostDistH[,,nb] = parallelRes[[7+8*(nb-1)]]
        timings = c(timings, parallelRes[[8+8*(nb-1)]])
      }
      
      
    }else{
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
            A[cA,basis$idx[[cA]]] = colMeans(basisSplines)[basis$idx[[cA]]]
          }
          e = matrix(0,ncol=1, nrow=(nConst))
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
        
        if(nSpTCov>0){MuG = Mu - SpTcov %*% gammaHist[1,]}else{MuG = Mu}
        
        Iter = diag(nSites)[match(y$ID, ID),]

        if(nComp>1){SpTcov = z0[match(y$ID, ID),]}
        Yp0 = crossprod_simple_triplet_matrix(t(Xtilde), betaHist[1,]) + SpTcov %*% gammaHist[1,]
        Ypp0 = Yp0 + Iter %*% betaIHist[1,]
        
        numLog = NULL
        for(t in 2:N.run){

          #browser()
          # print(t)
          # 
          # Updating alpha_j
          
          zCurr = matrix(zHist[t-1,,], nrow=nComp)
          for(j in 1:nComp){
            idxZ = which(y$ID %in% ID[which(apply(zCurr,2,which.max) == j)])
            if(length(idxZ)>0){
              BtildeZ = Btilde[idxZ,]
              tBB = crossprod_simple_triplet_matrix(BtildeZ)
             
              iSigmaAlpha = 1/priors$alpha$s0 * diag(ncol(alphaHist))
          
              #iSigmaObsB = 1/sigma2Hist[t-1]* tBB
              iSigmaObsB = tBB
          
              alphaPostSd = solve(iSigmaAlpha + iSigmaObsB)
              alphaPostMu = alphaPostSd %*% ( crossprod_simple_triplet_matrix(BtildeZ, MuG[idxZ]) + 
                                                iSigmaAlpha %*% matrix(priors$alpha$m0, ncol =1, nrow = nrow(iSigmaAlpha)))
              # alphaPostMu = alphaPostSd %*% ( sigma2Hist[t-1]^(-1) * crossprod_simple_triplet_matrix(BtildeZ, MuG[idxZ]) + 
              #                                   iSigmaAlpha %*% matrix(priors$alpha$m0, ncol =1, nrow = nrow(iSigmaAlpha)))
              
            }else{
              alphaPostSd = diag(priors$alpha$s0, nSplines)
              alphaPostMu = matrix(priors$alpha$m0, nrow = nSplines, ncol=1)
            }
          if(priors$alpha$constrained){
            # Using Rue2005 alg.
            
            nConst = nBasis
            if(nBasis == 3){nConst = 2}
            
            A = matrix(0,nrow = nConst, ncol = length(alphaPostMu))
            for(cA in 1:(nConst)){
              A[cA,basis$idx[[cA]]] = colMeans(basisSplines)[basis$idx[[cA]]]
            }
            e = matrix(0,ncol=1, nrow=(nConst))
            L = t(chol(solve(alphaPostSd)))
            zi = rmvnorm(1,rep(0, nrow(alphaPostSd)), diag(1, nrow(alphaPostSd)))
            v = solve(t(L),t(zi))
            xi = alphaPostMu + v
            Vnk = solve(alphaPostSd) %*% t(A)
            Wkk = A %*% Vnk
            Ukn = solve(Wkk) %*% t(Vnk)
            ci = A %*% xi - e
            alphaHist[t,,j] = xi - t(Ukn) %*% ci
            
            
            muAlphaC = alphaPostMu - solve(alphaPostSd) %*% t(A) %*% solve(A %*% solve(alphaPostSd) %*% t(A)) %*% (A %*% alphaPostMu - e)
            sigAlphaC = solve(alphaPostSd) - solve(alphaPostSd) %*% t(A) %*% solve(A %*% solve(alphaPostSd) %*% t(A)) %*% A %*% solve(alphaPostSd)
            
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
          
          #iSigmaObs = 1/sigma2Hist[t-1]* tXX
          iSigmaObs = tXX
          
          betaPostSd = solve(iSigmaBeta + iSigmaObs)
          ibetaPostSd = iSigmaBeta + iSigmaObs
          
          #betaPostMu = betaPostSd %*% (sigma2Hist[t-1]^(-1) * crossprod_simple_triplet_matrix(Xtilde, MuG) + 
          #                               iSigmaBeta %*% matrix(priors$beta$m0, ncol =1, nrow = nrow(iSigmaBeta)) )
          betaPostMu = betaPostSd %*% (crossprod_simple_triplet_matrix(Xtilde, MuG) + 
                                         iSigmaBeta %*% matrix(priors$beta$m0, ncol =1, nrow = nrow(iSigmaBeta)) )
          
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
          
          if(nSpTCov>0){
            
            if(nComp>1){SpTcov = t(zCurr[,match(y$ID, ID)])}
            
            tWW = crossprod(SpTcov, SpTcov)
            iSigmaGamma = 1/priors$gamma$s0 * diag(ncol(gammaHist))
            
            #iSigmaObsG = 1/sigma2Hist[t-1]* tWW
            iSigmaObsG = tWW
            
            gammaPostSd =solve(iSigmaGamma + iSigmaObsG)
            # gammaPostMu = gammaPostSd %*% (1/sigma2Hist[t-1]* crossprod(SpTcov, (Mu - Yp)) + 
            #                                  iSigmaGamma  %*% matrix(priors$gamma$m0, ncol =1, nrow = nrow(iSigmaGamma)))
            gammaPostMu = gammaPostSd %*% (crossprod(SpTcov, (Mu - Yp)) + 
                                             iSigmaGamma  %*% matrix(priors$gamma$m0, ncol =1, nrow = nrow(iSigmaGamma)))
            
            gammaHist[t,] = rmvnorm(1, gammaPostMu, gammaPostSd)
            
          }
          
          Yp = crossprod_simple_triplet_matrix(t(Xtilde), betaHist[t,]) + SpTcov %*% gammaHist[t,]
          
          # Site random effect
          
          Iter = diag(nSites)[match(y$ID, ID),]
           
          tII = crossprod(Iter, Iter)
          iSigmaBetaI = diag(nSites)
           
          iSigmaObsI = tII
           
          betaIPostSd =solve(iSigmaBetaI + iSigmaObsI)
          ibetaIPostSd = iSigmaBetaI + iSigmaObsI
        
          betaIPostMu = betaIPostSd %*% (crossprod(Iter, (Mu - Yp)) + 
                                           iSigmaBetaI  %*% matrix(rep(0, nSites), ncol =1, nrow = nrow(iSigmaBetaI)))
          
          Ab = kronecker(zCurr, diag(1))
          idxNC = which(rowSums(Ab)==0)
          if(length(idxNC)>0){
            Ab = matrix(Ab[-idxNC,], ncol = nSites)
          }
          eb = matrix(0,ncol=1, nrow=nrow(Ab))
          L = t(chol(ibetaIPostSd))
          zi = rmvnorm(1,rep(0, nrow(betaIPostSd)), diag(1, nrow(betaIPostSd)))
          v = solve(t(L),t(zi))
          xi = betaIPostMu + v
          Vnk = (betaIPostSd) %*% t(Ab)
          Wkk = Ab %*% Vnk
          Ukn = solve(Wkk) %*% t(Vnk)
          ci = Ab %*% xi - eb
          betaIHist[t,] = xi - t(Ukn) %*% ci

          # betaIHist[t,] = rmvnorm(1, betaIPostMu, betaIPostSd)
          
          # Updating sigma
          
          if(nSpTCov>0){Ypp = Yp + Iter %*% betaIHist[t,]}else{Ypp = Yp}
          #if(nSpTCov>0){Ypp = Yp}else{Ypp = Yp}
          
          sigma2PostA = priors$sigma2$a0 + nObs / 2
          sigma2PostB = priors$sigma2$b0 + sum((Mu - Ypp)^2, na.rm=T) / 2 #+ nObs * priors$beta$s0 / (nObs + priors$beta$s0) * mean((Mu - Ypp0)^2)/2
          
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
          if(!priors$beta$constrained){
            tauPostB = priors$tau$b0 + 
              (t(as.matrix(betaHist[t,]) - betaMu) %*% iSigmaBT %*% as.matrix(betaHist[t,]- betaMu))/2
            }

          if(priors$beta$constrained){
            #browser()
            
            idxNN = which(rowSums(zCurr)>0)
            
            Ab = kronecker(zCurr[idxNN,], diag(nBasis))
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
          MuG = Mu
          
          if(!is.null(Bds)){vois = getCanStatF(Bds, zz.temp, nComp)}else{
            vois = matrix(0, nSites, nComp)
          }
          
          YpComp = crossprod_simple_triplet_matrix(t(Btilde),alphaHist[t,,])
          if(nSpTCov == nComp){
            YpComp = YpComp + matrix(gammaHist[t,], ncol = nSpTCov, nrow = nrow(YpComp), byrow = T)
          }
          
          piCurr = matrix(piHist[t-1,,], nrow = nComp)
          zObs = NULL
          logProb = matrix(NA, nSites, nComp)
          for(i in 1:nSites){
            idxID = which(y$ID == ID[i])
            zObsProb = NULL
            for(j in 1:nComp){
              logProb[i,j] = sum(dnorm(MuG[idxID], YpComp[idxID,j], sigma2Hist[t], log=T), na.rm=T)* sqrt(2*pi*sigma2Hist[t]) + rhoHist[t] * vois[i,j]
              #if(sum(zCurr[j,])==0){logProb[i,j] = 0}
              zObsProb = cbind(zObsProb, (dnorm(MuG[idxID], YpComp[idxID,j], sigma2Hist[t], log=T)* sqrt(2*pi*sigma2Hist[t]) + t(log(piCurr))[i,j]))
              }
            zObs = rbind(zObs, rowSums(matrix(apply(exp(zObsProb), 1, function(x) rmultinom(1,1,x)), nrow = nComp)))
          }
          
          piCurr = matrix(piHist[t-1,,], nrow = nComp)
          
          #if(min(rowSums(zCurr))==0){browser()}
          
          logProb = logProb + t(log(piCurr))
          
          mlogProb = matrix(apply(logProb ,1, function(x) x - max(x, na.rm=T)), ncol=nComp, byrow=T)
          
          fProb = matrix(apply(mlogProb ,1, function(x) exp(x)/ sum(exp(x))), ncol=nComp, byrow=T)
          
          z.Num = apply(fProb,1, function(x) sample(1:nComp,1,prob = x))
          
          zHist[t,,] <- apply(fProb, 1, function(x) rmultinom(1,1,x))
          
          
          # Updating pi
          tempAlpha=NULL
          for(i in 1:nSites){
            #browser()
            alpha = priors$pi$alpha0[i,] + zObs[i,]#zHist[t,,i]
            piHist[t,,i] = bayesm::rdirichlet(c(as.matrix(alpha)))
            tempAlpha = rbind(tempAlpha, alpha)
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
            Ab = kronecker(zCurr, diag(2))
            idxNC = which(rowSums(Ab)==0)
            if(length(idxNC)>0){
              Ab = matrix(Ab[-idxNC,], ncol = nSites*nBasis)
            }
            eb = matrix(rowSums(Ab),ncol=1, nrow=nrow(Ab))
            
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
            prPi = prPi + sum(log(mixtools::ddirichlet(matrix(round(piHist[t,isZero,i],5), ncol = length(isZero)), alpha = tempAlpha[i,isZero])))
            prZ = prZ + log(fProb[i,which.max(zHist[t,,i])])
          }
          
          logPostDist[t] = LogL + prBeta + prGamma + prAlpha + 
            prSigma2 + prTau + prPhi+ prPi + prZ + prRho + prPi
          logLike[t] = LogL
          
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
        yHatH[,,nb] = yHatHist
        toc(log = TRUE, quiet = TRUE)
      }
      log.lst <- tic.log(format = FALSE)
      tic.clearlog()
      timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
    }
    
    
    
    RET = list(theta = list(alphaH = alphaH,betaH = betaH, betaIH = betaIH, gammaH = gammaH, 
                            sigma2H = sigma2H, tauH = tauH, phiH = phiH,rhoH = rhoH,
                            piH = piH, zH = zH),
               DistM = list(DM = DistMat, DMB = DistMatBeta),
               y = y, basis = basis, ID = ID ,coords = coords,cov = cov,sptcov = SpTcov, yHat = yHatH,
               priors = priors, execTime = timings, logPostDist = list(Post = logPostDistH, ll = logLikeH))
    
    return(RET)
    
  }


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
    
    if(is.null(SpTcov)){SpTcov = matrix(1, ncol = 1, nrow = nObs)}
    
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
        
        if(nSpTCov>0){MuG = Mu - SpTcov %*% gammaHist[1,]}else{MuG = Mu}
        
        numLog = NULL
        for(t in 2:N.run){
          #print(t)

          # Updating alpha_j
          zCurr = matrix(zHist[t-1,,], nrow=nComp)
          
          idxC = ((si-1)*nBasis+1):(si*nBasis)
          
          MuCorr = diag(betaHist[t-1,match(y$ID, levels(y$ID))*nBasis]^(-1)) %*% (MuG - betaHist[t-1,match(y$ID, levels(y$ID))*nBasis-1])

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
              #idxft = c(idxft, list(y$ID[idxZ]))
            }else{
              xList = list(matrix(y$date, ncol=1),matrix(y$date, ncol=1))
              xPred = list(matrix(y$date, ncol=1),matrix(y$date, ncol=1))
              tt = GPpred(xd = xPred, x = xList,y = c(MuCorr), param = param, kernel = kernelList)
              ft = c(ft, list(tt$mp))
            #dxft = c(idxft, list(NA))
            }
          }
          
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
          
          betaPostMu = betaPostSd %*% (crossprod_simple_triplet_matrix(Xtilde, MuG) + 
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
            gammaPostMu = gammaPostSd %*% (crossprod(SpTcov, (Mu - Yp)) + 
                                             iSigmaGamma  %*% matrix(priors$gamma$m0, ncol =1, nrow = nrow(iSigmaGamma)))
            
            gammaHist[t,] = rmvnorm(1, gammaPostMu, gammaPostSd)
            
          }
          
          # Updating sigma
          
          if(nSpTCov>0){Ypp = Yp + SpTcov %*% gammaHist[t,]}else{Ypp = Yp}
          
          sigma2PostA = priors$sigma2$a0 + nObs / 2
          sigma2PostB = priors$sigma2$b0 + sum((Mu - Ypp)^2, na.rm=T) / 2
          
          #if(t==10){browser()}
          
          sigma2Hist[t] = rigamma(1, sigma2PostA, sigma2PostB)
          
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
          zObs = NULL
          logProb = matrix(NA, nSites, nComp)
          for(i in 1:nSites){
            idxID = which(y$ID == ID[i])
            zObsProb = NULL
            for(j in 1:nComp){
              YpComp = ft[[j]][idxID] * betaHist[t,i*2] + betaHist[t,i*2-1]
              logProb[i,j] = sum(dnorm(MuG[idxID], YpComp, sigma2Hist[t], log=T), na.rm=T)* sqrt(2*pi*sigma2Hist[t]^2) + rhoHist[t] * vois[i,j]
              zObsProb = cbind(zObsProb, (dnorm(MuG[idxID], YpComp, sigma2Hist[t], log=T)* sqrt(2*pi*sigma2Hist[t]^2) + t(log(piCurr))[i,j]))
            }
            zObs = rbind(zObs, rowSums(matrix(apply((zObsProb), 1, function(x) rmultinom(1,1,exp(x-max(x)))), nrow = nComp)))
          }
          
          piCurr = matrix(piHist[t-1,,], nrow = nComp)
          
          logProb = logProb + t(log(piCurr))
          
          mlogProb = matrix(apply(logProb ,1, function(x) x - max(x, na.rm=T)), ncol=nComp, byrow=T)
          
          fProb = matrix(apply(mlogProb ,1, function(x) exp(x)/ sum(exp(x))), ncol=nComp, byrow=T)
          
          #z.Num = apply(fProb,1, function(x) sample(1:nComp,1,prob = x))
          
          zHist[t,,] <- apply(fProb, 1, function(x) rmultinom(1,1,x))
          
          
          # Updating pi
          tempAlpha=NULL
          for(i in 1:nSites){
            #browser()
            alpha = priors$pi$alpha0[i,] + zObs[i,]#zHist[t,,i]
            piHist[t,,i] = bayesm::rdirichlet(c(as.matrix(alpha)))
            tempAlpha = rbind(tempAlpha, alpha)
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
            prPi = prPi + ddirichlet(matrix(round(piHist[t,isZero,i],5), ncol = length(isZero)), alpha = tempAlpha[i,isZero], log=T, sum=T)
            prZ = prZ + log(fProb[i,which.max(zHist[t,,i])])
          }
          
          logPostDist[t] = LogL + prBeta + prGamma + prSigma2 + prTau + prPhi+ prPi + prZ + prRho + prPi
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
        toc(log = TRUE, quiet = TRUE)
      }
      log.lst <- tic.log(format = FALSE)
      tic.clearlog()
      timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))

    RET = list(theta = list(betaH = betaH, gammaH = gammaH, thetaGPH = thetaGPH,
                            sigma2H = sigma2H, tauH = tauH, phiH = phiH,rhoH = rhoH,
                            piH = piH, zH = zH),
               DistM = list(DM = DistMat, DMB = DistMatBeta),
               y = y, basis = basis, ID = ID ,coords = coords,cov = cov,sptcov = SpTcov,
               priors = priors, kernelList = kernelList,
               execTime = timings, logPostDist = list(Post = logPostDistH, ll = logLikeH))
    
    return(RET)
    
  }

"predictSPTM" <- 
function(SPTMresobj, transform = "none", posterior = F, newdata = NULL, keepRun = NULL){
  
  if(transform == "none"){
    fT = function(x){x}
  }
  if(transform == "log"){
    fT = function(x){exp(x)}
  }
  
    
    alphaH = SPTMresobj$GibbsOut$theta$alphaH
    betaH = SPTMresobj$GibbsOut$theta$betaH
    betaIH = SPTMresobj$GibbsOut$theta$betaIH
    sigma2H = SPTMresobj$GibbsOut$theta$sigma2H
    gammaH = SPTMresobj$GibbsOut$theta$gammaH
    tauH = SPTMresobj$GibbsOut$theta$tauH
    phiH = SPTMresobj$GibbsOut$theta$phiH
    tauTH = SPTMresobj$GibbsOut$theta$tauTH
    phiTH = SPTMresobj$GibbsOut$theta$phiTH
    piH = SPTMresobj$GibbsOut$theta$piH
    zH = SPTMresobj$GibbsOut$theta$zH
    yHat = SPTMresobj$GibbsOut$yHat
      
    y = SPTMresobj$GibbsOut$y
    basis = SPTMresobj$GibbsOut$basis
    basisSplines = do.call(cbind,basis$splines)
    
    ID = SPTMresobj$GibbsOut$ID
    cov = SPTMresobj$GibbsOut$cov
    coords = SPTMresobj$GibbsOut$coords[match(ID, SPTMresobj$GibbsOut$coords$site),]
    SpTcov = SPTMresobj$GibbsOut$sptcov
    
    tempBasis = SPTMresobj$GibbsOut$basis$tempBasis
    tempPeriod = SPTMresobj$GibbsOut$basis$tempPeriod
    
    nBasis = length(basis$idx)
    nSites = length(ID)
    nSplines = ncol(basisSplines)
    nCov = ncol(cov)
    nObs = nrow(y)
    Nrun = nrow(betaH)
    nHyper = 3 + 2 * as.numeric(SPTMresobj$tempRE == "corr")
    nComp = SPTMresobj$GibbsOut$priors$pi$nComp
    nSpTCov = ncol(SpTcov)
    if(is.null(keepRun)){
      keepRun = round(Nrun/2):Nrun
    }
    Xtilde = simple_triplet_zero_matrix(nrow = nObs, ncol = nSites*nCov*nBasis)
    Iter = diag(nSites)[match(y$ID, ID),]
    
    if(is.null(newdata)){
      
      newBetaH = NULL
     
      
      if(posterior){
        
        YpredH = colMeans(yHat[keepRun,,1])
        YpredHN = matrix(NA, nrow = length(keepRun), ncol = nObs)

        for(tt in 1:length(keepRun)){
          
          t = keepRun[tt]
          
          for(si in 1:nSites){
            idxC = ((si-1)*nBasis+1):(si*nBasis)
            xt = NULL
            for(jj in 1:(nBasis)){
              if(dim(alphaH)[3]>1){fTime = as.matrix(basisSplines[,basis$idx[[jj]]]) %*% alphaH[t,basis$idx[[jj]], apply(zH[t,,si,1], 1, which.max),1]}
              if(dim(alphaH)[3]==1){
                fTime = as.matrix(basisSplines[,basis$idx[[jj]]]) %*% alphaH[t,basis$idx[[jj]],1,1]
              }
              xt = cbind(xt,fTime * cov)
            }
            Xtilde[which(y$ID == ID[si]),idxC] = xt[which(y$ID == ID[si]),]
          }
          zCurr = matrix(zH[t,,,1], nrow = dim(zH)[2])
          if(dim(alphaH)[3]>1){SpTcov = t(zCurr[,match(y$ID,ID)])}
          
          Yp = crossprod_simple_triplet_matrix(t(Xtilde), betaH[t,,1]) + SpTcov %*% gammaH[t,,1]
          yHat = Yp + Iter %*% betaIH[t,,1]
          YpredHN[tt,] = yHat
          
        }
        
        YpredHN = YpredHN + rnorm(length(keepRun)*nObs, 0, mean(sigma2H[keepRun,,1]))
        
        YpredH025 = apply(YpredHN,2,quantile, 0.025)
        YpredH975 = apply(YpredHN,2,quantile, 0.975)
        
      }else{
        
        YpredH = colMeans(yHat[keepRun,,1])
        YpredH025 = apply(yHat[keepRun,,1] + rnorm(length(keepRun)*nObs, 0, mean(sigma2H[keepRun,,1])),2,quantile, 0.025)
        YpredH975 = apply(yHat[keepRun,,1] + rnorm(length(keepRun)*nObs, 0, mean(sigma2H[keepRun,,1])),2,quantile, 0.975)
        YpredHN =  colMeans(yHat[keepRun,,1]+ rnorm(length(keepRun)*nObs, 0, mean(sigma2H[keepRun,,1])))
      }
      
      
    }else{
      
      newLoc = newdata$loc
      newTime = newdata$time
      newCov = newdata$cov
      newID = newdata$ID
      
      nNewSites = nrow(newLoc)
      
      newZ = NULL
      for(nID in unique(newID)){
        if(nID %in% SPTMresobj$GibbsOut$ID){
          newZ = c(newZ, apply(apply(zH[keepRun,,,1], 2:3, mean),2,which.max)[match(nID, SPTMresobj$GibbsOut$ID)])
        }
      }
      YppredH = matrix(NA, nrow = length(keepRun), ncol = length(newTime))
      YppredHN = matrix(NA, nrow = length(keepRun), ncol = length(newTime))
      
      Xtilde = simple_triplet_zero_matrix(nrow = length(newTime), ncol = nNewSites*nCov*nBasis)
      newBB = NULL
      j=0
      for(i in tempPeriod){
        j = j +1 
        if(tempBasis == "bs"){
          if(j == 1){
            tTemp = as.Date(paste("2014-",format(newTime,"%m-%d"), sep = ""))
            nBBT = predict(basis$splines[[j]], tTemp)
          }
          if(j == 2){
            nBBT = predict(basis$splines[[j]], newTime)
            }
          if(j == 3){
            nBBT = matrix(rep(1, length(newTime)), ncol=1)
          }
        }
        if(tempBasis == "re"){
          nBBT = sapply(BbaseTr, function(x) as.numeric(as.numeric(format(t,i)) == x))
        }
        
        newBB = c(newBB,list(nBBT))
      }
      
      newsplines = do.call(cbind,newBB)
      
      Btilde = matrix(0, nrow = length(newTime), ncol = ncol(alphaH))
      
      BtildeO = matrix(0, nrow = nObs, ncol = ncol(alphaH))
      
      d11 = as.matrix(dist(coords[,2:3]))
      d12 = as.matrix(dist(rbind(coords[,2:3],newLoc)))[1:nSites,(nSites+1):(nSites + nNewSites)]
      d22 = as.matrix(dist(newLoc))
      
      newBetaH = matrix(NA,length(keepRun), ncol =  nBasis*nNewSites)
      
      colnames(newBetaH) <- rep(rownames(newLoc), each = nBasis)
 
      for(tt in 1:length(keepRun)){
        
        t = keepRun[tt]
        #if(t == 10){browser()}
        
        Sig11 = tauH[t,,1] * exp(-phiH[t,,1] * d11)
        Sig22 = tauH[t,,1] * exp(-phiH[t,,1] * d22)
        Sig12 = tauH[t,,1] * exp(-phiH[t,,1] * d12)
        SigP = Sig22 - t(Sig12) %*% solve(Sig11) %*% Sig12
        
        SigP = sapply(SigP, function(x) max(x,0))
        
        betaS = matrix(betaH[t,,1], ncol=nSites, byrow=F)
        muP = t(Sig12) %*% solve(Sig11) %*% t(betaS)
        
        SigB = kronecker(diag(nBasis), SigP)
        
        nBS = t(matrix(rmvnorm(1,c(muP),SigB), nrow = nNewSites))
        
        newBetaH[tt,] = nBS
        
        newBetaS = matrix(newBetaH[tt,], ncol=nNewSites, byrow=F, dimnames = list(paste0("beta-",1:nBasis), rownames(newLoc)))

        for(si in 1:nNewSites){
          for(nT in 1:nBasis){
            Btilde[which(newID == levels(newID)[si]),basis$idx[[nT]]] = newsplines[which(newID == levels(newID)[si]),basis$idx[[nT]]]*
              (newCov *
                 newBetaS[nT,which(colnames(newBetaS)==levels(newID)[si])]
               )
          }
        }

        YppredH[tt,] = Btilde %*% alphaH[t,,newZ,1] + gammaH[t,newZ,1] + betaIH[t,match(newID, SPTMresobj$GibbsOut$ID),1]
        
        for(si in 1:nNewSites){
          idxC = ((si-1)*nBasis+1):(si*nBasis)
          xt = NULL
          for(jj in 1:(nBasis)){
            if(dim(alphaH)[3]>1){fTime = as.matrix(newsplines[,basis$idx[[jj]]]) %*% alphaH[t,basis$idx[[jj]], apply(zH[t,,match(newID[si], SPTMresobj$GibbsOut$ID),1], 1, which.max),1]}
            if(dim(alphaH)[3]==1){
              fTime = as.matrix(newsplines[,basis$idx[[jj]]]) %*% alphaH[t,basis$idx[[jj]],1,1]
            }
            xt = cbind(xt,fTime * newCov)
          }
          Xtilde[which(newID == newID[si]),idxC] = xt[which(newID == newID[si]),]
        }
        
        Yp = crossprod_simple_triplet_matrix(t(Xtilde), newBetaH[tt,]) + gammaH[t,newZ,1]
        yHat = Yp + betaIH[t,match(newID, SPTMresobj$GibbsOut$ID),1]
        YppredHN[tt,] = yHat + rnorm(length(newTime),0,sigma2H[t,,1])
       
        #browser()
         
      }
      
      YpredH = colMeans(YppredH)
      YpredH025 = apply(YppredHN,2,quantile, 0.025)
      YpredH975 = apply(YppredHN,2,quantile, 0.975) 
      YpredHN = colMeans(YppredHN)
      
  }
  
  return(list(Yp = YpredH, YpwN = YpredHN, newB = newBetaH, CI025 = YpredH025, CI975 = YpredH975))
  
}

"plotSPTM" <- 
function(SPTMresobj, type = "prediction", color = NULL, cex.axis=NULL, n_x_lab = 7,
         siteDisp = NULL, siteSpTcov = NULL, keepRun = NULL){
  
    
    alphaH = SPTMresobj$GibbsOut$theta$alphaH
    thetaGPH = SPTMresobj$GibbsOut$theta$thetaGPH
    betaH = SPTMresobj$GibbsOut$theta$betaH
    gammaH = SPTMresobj$GibbsOut$theta$gammaH
    sigma2H = SPTMresobj$GibbsOut$theta$sigma2H
    tauH = SPTMresobj$GibbsOut$theta$tauH
    phiH = SPTMresobj$GibbsOut$theta$phiH
    tauTH = SPTMresobj$GibbsOut$theta$tauTH
    phiTH = SPTMresobj$GibbsOut$theta$phiTH
    piH = SPTMresobj$GibbsOut$theta$piH
    zH = SPTMresobj$GibbsOut$theta$zH
    yHat = SPTMresobj$GibbsOut$yHat
    
    kernelList = SPTMresobj$GibbsOut$kernelList
    
    if(SPTMresobj$model %in% "spatialMixture"){
      rhoH = SPTMresobj$GibbsOut$theta$rhoH
    }
    
    y = SPTMresobj$GibbsOut$y
    basis = SPTMresobj$GibbsOut$basis
    basisSplines = do.call(cbind,basis$splines)
    ID = SPTMresobj$GibbsOut$ID
    cov = SPTMresobj$GibbsOut$cov
    coords = SPTMresobj$GibbsOut$coords[match(SPTMresobj$GibbsOut$coords$site, ID),]
    Cells = SPTMresobj$Voronoi
    
    
    tempBasis = SPTMresobj$GibbsOut$basis$tempBasis
    tempPeriod = SPTMresobj$GibbsOut$basis$tempPeriod
    
    nBasis = length(basis$idx)
    nSites = length(ID)
    nSplines = ncol(basisSplines)
    nCov = ncol(cov)
    nObs = nrow(y)
    Nrun = nrow(betaH)
    nSpTcov = ncol(gammaH)
    nHyper = 3 + 2 * as.numeric(SPTMresobj$tempRE == "corr") + 6*as.numeric(SPTMresobj$tempRE == "gp")
    if(SPTMresobj$model %in% "spatialMixture"){
      nHyper = nHyper + 1
    }
    nComp = ncol(SPTMresobj$GibbsOut$priors$pi$alpha0)
    nBatch = dim(betaH)[3]
    
    if(is.null(keepRun)){
      keepRun = 1:Nrun
    }
    
    if(type == "residuals"){
      
      Yp = colMeans(yHat[-1,,1])
      
      residuals.df = y
      residuals.df$prediction = Yp
      residuals.df$residuals = Yp - y$obs
      
      zID = data.frame(ID = ID, Component = apply(zH[Nrun,,,1],3,which.max))
      
      residuals.df = merge(residuals.df, zID, by = "ID")
      
      par(mfrow =c(2,2),mar = c(5,4,4,1))
      
      x = residuals.df$date
      plot(x, rep(0, length(x)), type="n", main = "Residuals over time",ylab="residuals",xlab = "time",
           ylim = c(-max(abs(df$obs.data$obs), na.rm=T),max(abs(df$obs.data$obs), na.rm=T)))
      points(x, residuals.df$residuals)
      
      plot(0, type="n", xlim = range(y$obs, na.rm=T),ylab="predicted",xlab = "observed",
           ylim  = range(y$obs, na.rm=T), main = "Predicted vs Observed")
      abline(a=0,b=1)
      points(residuals.df$obs, residuals.df$prediction)
      
      hist(residuals.df$residuals, main = "Histogram of residuals")

      if(nlevels(residuals.df$Component)>1){
        boxplot(residuals ~ Component, data = residuals.df,main="Residuals per cluster")
      }else{
        boxplot(residuals ~ ID, data = residuals.df,main="Residuals per site")
      }
      
    }
    
    if(type == "tempBasis"){
      
      if(floor(sqrt(nBasis)) == sqrt(nBasis)){
        dim.plot = c(sqrt(nBasis), sqrt(nBasis))
      }else{
        dim.plot = c(floor(sqrt(nBasis)+1), floor(sqrt(nBasis)+1))
        reduction = floor((prod(dim.plot) - nBasis) / dim.plot[1])
        dim.plot = dim.plot - c(reduction, 0)
      }
      par(mfrow=c(nComp,nBasis), mar=c(2,2,1,1))
      
      if(is.null(color)){color = brewer.pal(max(nComp,3), "Set2")}
      
      if(SPTMresobj$tempRE == "gp"){
        nB = 1
        Intersect = mean(SPTMresobj$GibbsOut$theta$gammaH[,,nB])
        yTest = solve(diag(colMeans(SPTMresobj$GibbsOut$theta$betaH[,match(y$ID, levels(y$ID))*2,nB]))) %*% (y$obs - Intersect - colMeans(SPTMresobj$GibbsOut$theta$betaH[,match(y$ID, levels(y$ID))*2-1,nB]))
        
        par(mfrow=c(nBasis, nComp))
        Gr = apply(t(apply(SPTMresobj$GibbsOut$theta$piH[keepRun,,,nB],c(2,3),mean)),1,which.max)
        
        nK = length(attr(kernelList, "name"))
        
        param = parToList((c(colMeans(SPTMresobj$GibbsOut$theta$thetaGPH[keepRun,,nB]),
                             mean(SPTMresobj$GibbsOut$theta$sigma2H[keepRun,,nB]))), kernelList )
        
        for(ke in 1:nK){
        
          kernelListTMP = list(kernelList[[ke]])
          attr(kernelListTMP, "name") <- attr(kernelList, "name")[ke]
          attr(kernelListTMP, "parameters") <- attr(kernelList, "parameters")[ke]
          attr(kernelListTMP, "type") <- attr(kernelList, "type")[ke]
          
          paramTMP = parToList((c(unlist(param[[ke]]),
                              mean(SPTMresobj$GibbsOut$theta$sigma2H[keepRun,,nB]))), kernelListTMP )
          
        ttSea = NULL
        for(gr in sort(unique(Gr))){
          idxS = which(Gr == gr)
          yT = yTest[y$ID %in% levels(y$ID)[idxS]]
          
          dateNum = as.numeric(sort(unique(y$date)))
          xList = list(matrix(y$date[y$ID %in% levels(y$ID)[idxS]], ncol=1),
                       matrix(y$date[y$ID %in% levels(y$ID)[idxS]], ncol=1))
          
          xPred = list(matrix(seq(min(dateNum),max(dateNum),length.out = 100), ncol=1),
                       matrix(seq(min(dateNum),max(dateNum),length.out = 100), ncol=1))
          
          
          tt = GPpred(xd = xPred, x = xList,y = yT, param = paramTMP, kernel = kernelListTMP)
          ttP = GPpred(xd = xList, x = xList,y = yT, param = paramTMP, kernel = kernelListTMP)
          
          datePred = seq(min(y$date),max(y$date),length.out = 100)
          plot(y$date[y$ID %in% levels(y$ID)[idxS[1]]], y$obs[y$ID %in% levels(y$ID)[idxS[1]]], ylim =c(-4,3), 
               main = attr(kernelListTMP, "name"))
          points(datePred, 
                 tt$mp*colMeans(SPTMresobj$GibbsOut$theta$betaH[keepRun,,nB])[2*idxS[1]]+colMeans(SPTMresobj$GibbsOut$theta$betaH[keepRun,,nB])[2*idxS[1]-1] + Intersect, type="l", col="red")
          for(i in idxS){
            points(y$date[y$ID %in% levels(y$ID)[c(i)]], y$obs[y$ID %in% levels(y$ID)[c(i)]])
            points(datePred , 
                   tt$mp*colMeans(SPTMresobj$GibbsOut$theta$betaH[keepRun,,nB])[2*i]+colMeans(SPTMresobj$GibbsOut$theta$betaH[keepRun,,nB])[2*i-1] + 
                     Intersect, type="l", col="red")
          }
          points(datePred , tt$mp +Intersect, col="blue", type="l", lwd=2)
          ttSea = c(ttSea, list(ttP))
        }
        
        }
      }else{
      if(SPTMresobj$GibbsOut$basis$tempBasis == "re"){
        
        for(i in 1:nBasis){
          
          xBasis= format(SPTMresobj$GibbsOut$y$date,SPTMresobj$GibbsOut$basis$tempPeriod[i])
          xB = as.numeric(unique(xBasis))
          xBs = sort(xB)
          idxxB = match(xB, xBs)
          
          plot(xBs,rep(2*max(SPTMresobj$GibbsOut$theta$alphaH[,SPTMresobj$GibbsOut$basis$idx[[i]],]), length(xBs)),
               ylim =  quantile(SPTMresobj$GibbsOut$theta$alphaH[,SPTMresobj$GibbsOut$basis$idx[[i]],],c(0.001,0.999)),
               main = SPTMresobj$GibbsOut$basis$tempPeriod[i], ylab="")
          abline(h = 0, col="grey")
          for(nb in 1:nBatch){
            yBasis = SPTMresobj$GibbsOut$theta$alphaH[,SPTMresobj$GibbsOut$basis$idx[[i]],nb]
            points(xBs, colMeans(yBasis)[idxxB], col = color[i])
            points(xBs, colMeans(yBasis)[idxxB], col = color[i],type = "l", lty = 2)
            arrows(xBs, apply(yBasis,2,function(x) quantile(x,0.05))[idxxB],
                   xBs, apply(yBasis,2,function(x) quantile(x,0.95))[idxxB],
                   col = color[i], code = 0)
          }
        }
        
      }else{
        for(k in 1:nComp){
        
        for(i in 1:nBasis){
          t0 = y$date
          if(basis$tempPeriod[i] == "%m"){
            t0 = as.Date(paste("2014-",format(t0,"%m-%d"), sep = ""))
          }

          xPred = seq(min(t0), max(t0), "week")
          xBasis = predict(basis$splines[[i]],xPred)
          
          yPred = xBasis %*% t(alphaH[keepRun,basis$idx[[i]],k,1])
          
          plot(xPred,rowMeans(yPred),
               ylim =  quantile(yPred, c(0.01,0.99)),
               main = basis$tempPeriod[i],lty=1, ylab="", type="l", lwd=2,col = color[i])
          abline(h = 0, col="grey")
          polygon(c(xPred, rev(xPred)), c(apply(yPred,1,quantile,0.05), rev(apply(yPred,1,quantile,0.95))),
                  border=NA, col = adjustcolor(color[i], alpha.f = 0.15))
          
          if(nBatch > 1){
            for(nb in 2:nBatch){
                yPred = xBasis %*% t(alphaH[keepRun,basis$idx[[i]],k,nb])
                points(xPred,rowMeans(yPred),
                       ylab="", type="l", lwd=2,col = color[i], lty=nb)
                polygon(c(xPred, rev(xPred)), c(apply(yPred,1,quantile,0.05), rev(apply(yPred,1,quantile,0.95))),
                        border=NA, col = adjustcolor(color[i], alpha.f = 0.15))
            }
          }
          }
          

        }
      }
      }
    }
    
    if(type == "spatialBasis"){
      
      if(floor(sqrt(nBasis)) == sqrt(nBasis)){
        dim.plot = c(sqrt(nBasis), sqrt(nBasis))
      }else{
        dim.plot = c(floor(sqrt(nBasis)+1), floor(sqrt(nBasis)+1))
        reduction = floor((prod(dim.plot) - nBasis) / dim.plot[1])
        dim.plot = dim.plot - c(reduction, 0)
      }
      par(mfrow=dim.plot, mar=c(3,3,2,1))
      
      if(is.null(color)){color = brewer.pal(max(nBasis,3), "Set2")}
      
      colSeq = colorRampPalette(c(color[1], color[2]))(24)
      
      colSeq = brewer.pal(9,"YlOrRd")
      
      for(i in 1:nBasis){
        
        idxBeta = (0:(nSites-1))*nBasis+i
        
        zBasis = SPTMresobj$GibbsOut$theta$betaH[,idxBeta,1]
        
        zBasisIQ = diff(apply(zBasis,2,function(x)quantile(x,c(0.25,0.75))))
        
        zBasisIdx <- cut(colMeans(zBasis), breaks = seq(quantile(zBasis,0.01), quantile(zBasis,0.99), len = 10),
                         include.lowest = TRUE)
        
        zBasisColor = colSeq[zBasisIdx]
        
        xyBasis= SPTMresobj$GibbsOut$coords[,2:3]
        
        plot(xyBasis, col = zBasisColor, bg = zBasisColor, pch=19, main = substitute(beta[i], list(i = i)))
        
        yArr =  zBasisIQ / max(zBasisIQ) / 100  / diff(range(xyBasis[,2]))
        
        arrows(xyBasis[,1] + 0.01*diff(range(xyBasis[,1])), xyBasis[,2] - yArr/2,
               xyBasis[,1] + 0.01*diff(range(xyBasis[,1])), xyBasis[,2] + yArr/2, code = 0,
               col = zBasisColor)
        
        legend("bottomright", levels(zBasisIdx), col = colSeq, pch = 19, cex=0.75, pt.cex=0.75, bty='n')
      }
      
    }
    
    if(type == "alpha"){
      
      if(floor(sqrt(nBasis*nComp)) == sqrt(nBasis*nComp)){
        dim.plot = c(sqrt(nBasis*nComp), sqrt(nBasis*nComp))
      }else{
        dim.plot = c(floor(sqrt(nBasis*nComp)+1), floor(sqrt(nBasis*nComp)+1))
        reduction = floor((prod(dim.plot) - nBasis*nComp) / dim.plot[1])
        dim.plot = dim.plot - c(reduction, 0)
      }
      par(mfrow=dim.plot, mar=c(3,3,2,1))
      
      if(is.null(color)){color = brewer.pal(min(max(max(unlist(lapply(SPTMresobj$GibbsOut$basis$idx, length))),3),12), "Set3")}
      
      for(i in 1:nBasis){
        for(j in 1:nComp){
          plot(-100,-1,
               main = paste0(SPTMresobj$GibbsOut$basis$tempPeriod[i], " - Component ", j),
               xlim = c(0,nrow(alphaH)), ylim = range(alphaH[,SPTMresobj$GibbsOut$basis$idx[[i]],j,1]))
          for(k in SPTMresobj$GibbsOut$basis$idx[[i]]){
            points(alphaH[,k,j,1], type="l", col  = color[k %% 12 + 1], lwd=2)
          }
        }
      }
    }
    
    if(type == "beta"){
      
      if(floor(sqrt(nSites)) == sqrt(nSites)){
        dim.plot = c(sqrt(nSites), sqrt(nSites))
      }else{
        dim.plot = c(floor(sqrt(nSites)+1), floor(sqrt(nSites)+1))
        reduction = floor((prod(dim.plot) - nSites) / dim.plot[1])
        dim.plot = dim.plot - c(reduction, 0)
      }
      par(mfrow=dim.plot, mar=c(3,3,2,1))
      
      if(is.null(color)){color = brewer.pal(min(max(max(unlist(lapply(SPTMresobj$GibbsOut$basis$idx, length))),3),12), "Set3")}
      
      for(i in 1:nSites){
        idxSite = 1:nBasis + (i-1)*nBasis
        plot(-100,-1,
             main = ID[i],
             xlim = c(0,nrow(betaH)), ylim = range(betaH[keepRun,idxSite,1]))
        for(j in idxSite){
          points(betaH[,j,1], type="l", col  = color[j %% 12 + 1], lwd=2)
        }
      }
      
    }
    
    if(type == "hyperparameter"){
      if(floor(sqrt(nHyper)) == sqrt(nHyper)){
        dim.plot = c(sqrt(nHyper), sqrt(nHyper))
      }else{
        dim.plot = c(floor(sqrt(nHyper)+1), floor(sqrt(nHyper)+1))
        reduction = floor((prod(dim.plot) - nHyper) / dim.plot[1])
        dim.plot = dim.plot - c(reduction, 0)
      }
      par(mfrow=dim.plot, mar=c(3,3,2,1))
      
      if(is.null(color)){color = brewer.pal(12, "Set3")}
      
      plot(keepRun,sigma2H[keepRun,,1], type="l", 
           col  = adjustcolor(color[1],alpha.f = 0.3), lwd=2, main = expression(sigma^2), 
           ylim = range(sigma2H[keepRun,,]))
      points(keepRun,cumsum(sigma2H[keepRun,,1])/1:length(keepRun), type="l", col  = adjustcolor(color[1],alpha.f = 0.8), lwd=3)
      if(nBatch>1){
        for(nb in 2:nBatch){
          points(keepRun,sigma2H[keepRun,,nb], type="l", col  = adjustcolor(color[1], alpha.f = 0.3), lwd=2)
          points(keepRun,cumsum(sigma2H[keepRun,,nb])/1:length(keepRun), type="l", col  = adjustcolor(color[1],alpha.f = 0.8), lwd=3)
        }
      }
      plot(keepRun,tauH[keepRun,,1], type="l", col  = adjustcolor(color[2],alpha.f = 0.5), lwd=2, 
           main = expression(tau), ylim = range(tauH[keepRun,,]))
      points(keepRun,cumsum(tauH[keepRun,,1])/1:length(keepRun), type="l", col  = adjustcolor(color[2],alpha.f = 0.8), lwd=3)
      if(nBatch>1){
        for(nB in 2:nBatch){
          points(keepRun,tauH[keepRun,,nB], type="l", col  = adjustcolor(color[2], alpha.f = 0.5), lwd=2)
          points(keepRun,cumsum(tauH[keepRun,,nB])/1:length(keepRun), type="l", col  = adjustcolor(color[2],alpha.f = 0.8), lwd=3)
        }
      }
      
      plot(keepRun,phiH[keepRun,,1], type="l", col  = adjustcolor(color[3],alpha.f = 0.3), lwd=2, 
           main = expression(phi), ylim = range(phiH[keepRun,,]))
      points(keepRun,cumsum(phiH[keepRun,,1])/1:length(keepRun), type="l", col  = adjustcolor(color[3],alpha.f = 0.8), lwd=3)
      if(nBatch>1){
        for(nB in 2:nBatch){
          points(keepRun,phiH[keepRun,,nB], type="l", col  = adjustcolor(color[3], alpha.f = 0.3), lwd=2)
          points(keepRun,cumsum(phiH[keepRun,,nB])/1:length(keepRun), type="l", col  = adjustcolor(color[3],alpha.f = 0.8), lwd=3)
        }
      }
      
      if(SPTMresobj$model %in% "spatialMixture"){
        plot(keepRun,rhoH[keepRun,,1], type="l", col  = adjustcolor(color[4],alpha.f = 0.3), lwd=2, 
             main = expression(rho), ylim = range(rhoH))
        points(keepRun,cumsum(rhoH[keepRun,,1])/1:length(keepRun), type="l", col  = adjustcolor(color[4],alpha.f = 0.8), lwd=3)
        if(nBatch>1){
          for(nB in 2:nBatch){
            points(keepRun,rhoH[keepRun,,nB], type="l", col  = adjustcolor(color[4], alpha.f = 0.3), lwd=2)
            points(keepRun,cumsum(rhoH[keepRun,,nB])/1:length(keepRun), type="l", col  = adjustcolor(color[4],alpha.f = 0.8), lwd=3)
          }
        }
      }
      
      if(SPTMresobj$tempRE == "corr"){
        
        plot(tauTH, type="l", col  = color[2], lwd=2, main = expression(tau[time]))
        
        plot(phiTH, type="l", col  = color[3], lwd=2, main = expression(phi[time]))
      }
      if(SPTMresobj$tempRE == "gp"){
        parName = c(expression(k[l][-var]),expression(k[l][-range]),
                    expression(k[sea][-var]),expression(k[sea][-shape]),expression(k[sea][-amp]),expression(k[sea][-freq]))
        for(par in 1:ncol(thetaGPH[keepRun,,1])){
        plot(keepRun,thetaGPH[keepRun,par,1], type="l", col  = adjustcolor(color[par+3],alpha.f = 0.3), lwd=2, 
             main = parName[par], ylim = range(thetaGPH[keepRun,par,]))
        points(keepRun,cumsum(thetaGPH[keepRun,par,1])/1:length(keepRun), type="l", col  = adjustcolor(color[par+3],alpha.f = 0.8), lwd=3)
        if(nBatch>1){
          for(nB in 2:nBatch){
            points(keepRun,thetaGPH[keepRun,par,nB], type="l", col  = adjustcolor(color[par+3], alpha.f = 0.3), lwd=2)
            points(keepRun,cumsum(thetaGPH[keepRun,par,nB])/1:length(keepRun), type="l", col  = adjustcolor(color[par+3],alpha.f = 0.8), lwd=3)
          }
        }
        }
      }
    }
    
    if(type == "cluster"){
      
      dimS = ceiling((nSites +4)/4 )
      
      par(mfrow=c(dimS, dimS), mar=c(2,2,1,1))
      
      LO = matrix(NA,dimS,dimS)
      LO[,1:2] = 1:(2*dimS)
      LO[1:2,3:dimS] = (2*dimS):(4*dimS-5) + 1
      LO[which(LO>nSites, arr.ind = T)] = nSites+2
      LO[3:dimS, 3:dimS] = nSites+1
      layout(LO)
      
      if(is.null(color)){color = brewer.pal(min(max(max(unlist(lapply(SPTMresobj$GibbsOut$basis$idx, length))),3),12), "Set3")}
      
      xPol = cos(2*pi * (c((1:nComp),1)-1) / nComp + pi/2)
      yPol = sin(2*pi * (c((1:nComp),1)-1) / nComp + pi/2)
      
      for(i in 1:nSites){
        
        prob = apply(zH[,,i,1],2,sum) / Nrun
        
        prob = prob / max(prob) * 0.66
        
        prob = c(prob, prob[1])
        
        # plot(prob*xPol, prob*yPol, xlim = c(-1.1,1.1), ylim = c(-1.1,1.1), type="l", axes = F, bty = 'n',
        #      col = color[which.max(apply(zH[,,i,1],2,sum) / Nrun)])
        # 
        # points(xPol, yPol, col = color[1:nComp])
        
        pHat = apply(piH[keepRun,,i,1],2,mean)
        pHatq = apply(piH[keepRun,,i,1],2,quantile,prob = c(0.1,0.9))
        plot(1:nComp, pHat, type="p", xlim = c(0.5, nComp+.5), col = adjustcolor(color[1:nComp], alpha.f = 1),
             ylim = c(-0.01,1.01), lwd = 2, main = ID[i], cex.main = 0.75, cex.axis=0.6)
        #points(1:nComp, apply(piH[,,i,1],2,mean), col = color[which.max(apply(zH[,,i,1],2,sum) / Nrun)])
        #arrows(1:nComp, apply(piH[,,i,1],2,quantile, 0.05), 1:nComp, apply(piH[,,i,1],2,quantile, 0.95), code = 0,
        #       col = color[which.max(apply(zH[,,i,1],2,sum) / Nrun)])
        arrows(1:nComp, pHatq[1,], 1:nComp, pHatq[2,], code = 0,
               col = color[1:nComp],#color[which.max(apply(zH[,,i,1],2,sum) / Nrun)], 
               lwd=2, lty=2)
        abline(h = c(0,1))
        for(jj in 1:nBatch){
        pHat = apply(piH[keepRun,,i,jj],2,mean)
        pHatq = apply(piH[keepRun,,i,jj],2,quantile,prob = c(0.1,0.9))
        points(1:nComp, pHat, col = adjustcolor(color[1:nComp], alpha.f = 1),
             lwd = 2)
        #points(1:nComp, apply(piH[,,i,1],2,mean), col = color[which.max(apply(zH[,,i,1],2,sum) / Nrun)])
        #arrows(1:nComp, apply(piH[,,i,1],2,quantile, 0.05), 1:nComp, apply(piH[,,i,1],2,quantile, 0.95), code = 0,
        #       col = color[which.max(apply(zH[,,i,1],2,sum) / Nrun)])
        arrows(1:nComp, pHatq[1,], 1:nComp, pHatq[2,], code = 0,
               col = color[1:nComp],#color[which.max(apply(zH[,,i,1],2,sum) / Nrun)], 
               lwd=2, lty=2)
        }
      }
      
      xyBasis= SPTMresobj$GibbsOut$coords[match(ID, SPTMresobj$GibbsOut$coords$site),2:3]
      
      zBasisColor = color[apply(apply(piH[keepRun,,,1],2:3,mean),2,which.max)]
      
      plot(xyBasis, col = zBasisColor, bg = zBasisColor, pch=19, main = "clustering")
      
      for(jj in 1:nBatch){
        zBasisColor = color[apply(apply(piH[keepRun,,,jj],2:3,mean),2,which.max)]
        xyB = xyBasis + matrix(c(0.01*diff(range(xyBasis[,1])),0), ncol = 2, nrow = nrow(xyBasis), byrow=T)
        points(xyB, col = zBasisColor, bg = zBasisColor, pch=19, main = "clustering")
      }
      
      legend("bottomright", paste("Cluster", 1:(nComp)), col = color[1:(nComp)], pch = 19, cex=0.75, pt.cex=0.75, bty='n')
      
    }
    
    if(type == "map"){
      
     
      par(mfrow=c(1,1), mar=c(2,1,1,1))
      
      coords$comp = apply(apply(zH[,,,1], 2:3,mean),2,which.max)
      
      if(is.null(color)){
        color = brewer.pal(min(max(max(unlist(lapply(SPTMresobj$GibbsOut$basis$idx, length))),3),12), "Set3")
        }
     
      plot(Cells[[1]], col = adjustcolor("grey", alpha.f = 0.3))
      if(!is.null(SpatialPolygon)){
        plot(SpatialPolygon, add=T, col=adjustcolor("blue", alpha.f = 0.25), lwd=2)
      }
      points(coords[,2:3], col = color[coords$comp], pch=19, cex=1)
    }
    
    if(type == "covariates"){
      if(floor(sqrt(nSpTcov)) == sqrt(nSpTcov)){
        dim.plot = c(sqrt(nSpTcov), sqrt(nSpTcov))
      }else{
        dim.plot = c(floor(sqrt(nSpTcov)+1), floor(sqrt(nSpTcov)+1))
        reduction = floor((prod(dim.plot) - nSpTcov) / dim.plot[1])
        dim.plot = dim.plot - c(reduction, 0)
      }
      par(mfrow=dim.plot, mar=c(3,3,2,1))
      
      if(is.null(color)){color = brewer.pal(12, "Set3")}
      
      for(i in 1:nSpTcov){
        plot(keepRun,gammaH[keepRun,i,1], type="l", col  = adjustcolor(color[i],0.35), lwd=2, main = bquote(gamma[.(i)]), ylim = quantile(gammaH[,i,], c(0.001,0.999)))
        points(keepRun,cumsum(gammaH[keepRun,i,1])/(1:length(keepRun))*0.9 + 0.1*gammaH[keepRun,i,1], type="l", lwd=2, col=color[i])
        abline(h = 0, lty = 1)
        
        for(nb in 2:nBatch){
          points(keepRun,gammaH[keepRun,i,nb], type="l", col  = adjustcolor(color[i],0.35), lwd=2)
          points(keepRun,cumsum(gammaH[keepRun,i,nb])/(1:length(keepRun))*0.9 + 0.1*gammaH[keepRun,i,nb], type="l", lwd=2, col=color[i])
        }
        abline(h = quantile(gammaH[keepRun,i,], c(0.05,0.95)), lty = 2, col = "grey", lwd=2)
      }
      
    }
    
}

"sptmRes2mcmc" <- 
function(SPTMresobj, idxRem=NULL){
  
    idxKeep = which(names(SPTMresobj$GibbsOut$theta) %in% c("alphaH", "betaH","sigma2H", "tauH", "phiH","rhoH"))
    
    sumTheta = lapply(SPTMresobj$GibbsOut$theta[idxKeep], function(x) prod(dim(x)[-c(1,length(dim(x)) )]))
    
    chainLength = nrow(SPTMresobj$GibbsOut$theta$betaH)
    nPar = sum(unlist(sumTheta))
    nBatch = dim(SPTMresobj$GibbsOut$theta$betaH)[3]
    
    mcmcChain = array(NA, dim = c(chainLength,nPar,nBatch))
    rowJ = 0
    for(j in idxKeep){
      dimJ = dim(SPTMresobj$GibbsOut$theta[[j]])
      if(prod(dimJ[-c(1, length(dimJ))])>0){
        rowJ = max(rowJ) + 1:prod(dimJ[-c(1, length(dimJ))])
        temp = SPTMresobj$GibbsOut$theta[[j]]
        colnames(temp) = paste(names(sumTheta)[j],1:ncol(SPTMresobj$GibbsOut$theta[[j]]))
        mcmcChain[,rowJ,] = temp
      }
    }
    
    if(!is.null(idxRem)){mcmcChain = mcmcChain[,-idxRem,]}
    eval(parse(text = paste0("combinedchain = mcmc.list(", paste0("mcmc(mcmcChain[,,",1:nBatch,"])", collapse = ","),")")))
  return(combinedchain)
}

"MurraySchemePotts" <- 
  function(beta.init, range.beta = c(0, 1), n.colors=2, Coords = NULL, Nb = NULL, col.obs, N.run=100, print = FALSE){
  
  beta.hist = matrix(0, nrow = N.run, ncol = length(beta.init))
  beta.hist[1,] = beta.init
  
  zz = as.matrix(col.obs)
  
  n.obs = length(zz)
  
  rate = 0
  
  beta = as.matrix(beta.init)
  
  if(!is.null(Coords)){Coords = as.matrix(Coords)}
  
  if(is.null(Nb)){bds = getBonds(Coords, NN = 4, th = 2)
  }else{bds = Nb}
  
  col.sum = getCanStat(bds, zz, n.colors)
  
  for(t in 2:N.run){
    
    if(print == TRUE){if(t/(N.run/10) == round(t/(N.run/10))){print(t)}}
    
    # Simulate beta.p ~ Q(.|beta)
    
    if(length(beta)>1){
      beta.p = runif(length(beta), sapply(beta, function(x) max(9*x/10,range.beta[1])),
                     sapply(beta, function(x) min(11*x/10,range.beta[2])))
    } else {
      beta.p = runif(1, max(9*beta/10,range.beta[1]),min(11*beta/10,range.beta[2]))
    }
    
    # Simulating z.p ~ f(.|beta.p) using Swendsen-Wang Alg. Need to reach stationary distribution
    
    beta.temp = rep(beta.p, n.colors)
    
    res = SwendsenWangAlg(beta.temp, Neigh_Bonds = bds, col.obs = zz, Nrun = 20, sampling = "uniform")
    
    z.p = res$col.new
    
    
    # Calculate probability of acceptance
    
    
    col.sum.p = getCanStat(bds, z.p, n.colors)
    
    p.num = sapply(11/10*beta.p, function(x) min(x,range.beta[2])) - 9 * beta.p/10
    
    p.denom = sapply(11/10*beta, function(x) min(x,range.beta[2])) - 9 * beta/10
    
    p.move =  exp(sum((beta.p - beta) * (col.sum - col.sum.p))) * prod(p.num / p.denom)
    
    p.move = min(p.move, 1, na.rm=T)
    
    # Exchange with probability p.move
    
    if(runif(1) < p.move){
      beta = beta.p
      rate = rate + 1
    }
    
    beta.hist[t,] = beta
    
  }
  
  RET = list(beta.old = beta.init, beta.new = beta, beta.sample = beta.hist,
             rate = rate)
  return(RET)
  }

"SwendsenWangAlg" <- 
  function(beta, Coord = NULL, Neigh_Bonds = NULL, col.obs, Nrun = 10 ,sampling = "else"){
  
  n.vert = length(col.obs)
  
  n.colors = length(beta)
  
  # Get neighboorhood bonds
  
  if(is.null(Neigh_Bonds)){Neigh_Bonds = getBonds(Coord, NN = 4, th = 2)}
  
  # Select bonds
  col.temp = col.obs
  
  z.p = loopSW(Bds = Neigh_Bonds, Cols = col.obs, ncolors = n.colors, Nrun = Nrun, Betas = beta)
  
  RET = list(col.new = z.p, col.old = col.obs)
  return(RET)
  }

"comPts" <- 
  function(tiles_list){
  bds = NULL
  for(i in 1:(length(tiles_list)-1)){
    temp.tile = tiles_list[[i]]
    temp.x = temp.tile$x
    temp.y = temp.tile$y
    for(j in (i+1):length(tiles_list)){
      temptemp.tile = tiles_list[[j]]
      dist.tiles = dist(x = cbind(c(temp.x,temptemp.tile$x), c(temp.y,temptemp.tile$y)))
      if(min(dist.tiles)==0){
        bds = rbind(bds, c(i,j))
      }
    }
  }
  return(bds)
}


"marginalPlot" <- 
  function(SPTMresobj, par = "hyperparameter", color = NULL, keepRun = NULL){
    # Based on Gelfand, Smith, Lee (1992)
    if(SPTMresobj$model == "noMixture"){
      
      alphaH = SPTMresobj$GibbsOut$theta$alphaH
      betaH = SPTMresobj$GibbsOut$theta$betaH
      gammaH = SPTMresobj$GibbsOut$theta$gammaH
      sigma2H = SPTMresobj$GibbsOut$theta$sigma2H
      tauH = SPTMresobj$GibbsOut$theta$tauH
      phiH = SPTMresobj$GibbsOut$theta$phiH
      tauTH = SPTMresobj$GibbsOut$theta$tauTH
      phiTH = SPTMresobj$GibbsOut$theta$phiTH
      
      
      y = SPTMresobj$GibbsOut$y
      basis = SPTMresobj$GibbsOut$basis
      ID = SPTMresobj$GibbsOut$ID
      cov = SPTMresobj$GibbsOut$cov
      coords = SPTMresobj$GibbsOut$coords[match(SPTMresobj$GibbsOut$coords$site, ID),]
      basisSplines = do.call(cbind,basis$splines)
      priors = SPTMresobj$GibbsOut$priors
      DistMat = as.matrix(dist(coords[match(coords$site, ID) ,2:3]))
      
      tempBasis = SPTMresobj$GibbsOut$basis$tempBasis
      tempPeriod = SPTMresobj$GibbsOut$basis$tempPeriod
      
      nBasis = length(basis$idx)
      nSites = length(ID)
      nSplines = ncol(basisSplines)
      nCov = ncol(cov)
      nObs = nrow(y)
      Nrun = nrow(alphaH)
      nHyper = 3 + 2 * as.numeric(SPTMresobj$tempRE == "corr")
      nSpTcov = ncol(gammaH)
      nBatch = dim(alphaH)[3]
      
      if(is.null(keepRun)){
        keepRun = 1:Nrun
      }
      
      if(type == "alpha"){
        
        if(floor(sqrt(nBasis)) == sqrt(nBasis)){
          dim.plot = c(sqrt(nBasis), sqrt(nBasis))
        }else{
          dim.plot = c(floor(sqrt(nBasis)+1), floor(sqrt(nBasis)+1))
          reduction = floor((prod(dim.plot) - nBasis) / dim.plot[1])
          dim.plot = dim.plot - c(reduction, 0)
        }
        par(mfrow=dim.plot, mar=c(3,3,2,1))
        
        if(is.null(color)){color = brewer.pal(min(max(max(unlist(lapply(SPTMresobj$GibbsOut$basis$idx, length))),3),12), "Set3")}
        
        for(i in 1:nBasis){
          plot(-100,-1,
               main = SPTMresobj$GibbsOut$basis$tempPeriod[i],
               xlim = c(0,nrow(alphaH)), ylim = range(alphaH[,SPTMresobj$GibbsOut$basis$idx[[i]],]))
          for(j in SPTMresobj$GibbsOut$basis$idx[[i]]){
            for(nb in 1:nBatch){
              points(alphaH[,j,nb], type="l", col  = color[j %% 12 + 1], lwd=0.5)
            }
          }
        }
        
      }
      
      if(type == "beta"){
        
        if(floor(sqrt(nSites)) == sqrt(nSites)){
          dim.plot = c(sqrt(nSites), sqrt(nSites))
        }else{
          dim.plot = c(floor(sqrt(nSites)+1), floor(sqrt(nSites)+1))
          reduction = floor((prod(dim.plot) - nSites) / dim.plot[1])
          dim.plot = dim.plot - c(reduction, 0)
        }
        par(mfrow=dim.plot, mar=c(3,3,2,1))
        
        if(is.null(color)){
          color = brewer.pal(min(max(max(unlist(lapply(SPTMresobj$GibbsOut$basis$idx, length))),3),12), "Set3")}
        
        for(i in 1:nSites){
          idxSite = 1:nBasis + (i-1)*nBasis
          plot(-100,-1,
               main = ID[i],
               xlim = c(0,nrow(betaH)), ylim = range(betaH[,idxSite,1]))
          for(j in idxSite){
            points(betaH[,j,1], type="l", col  = color[j %% 12 + 1], lwd=2)
          }
        }
        
      }
      
      if(type == "hyperparameter"){
        if(floor(sqrt(nHyper)) == sqrt(nHyper)){
          dim.plot = c(sqrt(nHyper), sqrt(nHyper))
        }else{
          dim.plot = c(floor(sqrt(nHyper)+1), floor(sqrt(nHyper)+1))
          reduction = floor((prod(dim.plot) - nHyper) / dim.plot[1])
          dim.plot = dim.plot - c(reduction, 0)
        }
        par(mfrow=dim.plot, mar=c(3,3,2,1))
        
        if(is.null(color)){color = brewer.pal(12, "Set3")}
        
        # Sigma2
        pSigma = pPhi = pTau = array(NA, dim = c(Nrun, 200,nBatch))
        
        betaTilde = simple_triplet_zero_matrix(nrow = nrow(y), ncol = dim(alphaH)[2])
        Mu = y$obs
        seqTau = seq(min(tauH[keepRun,,]), max(tauH[keepRun,,]), length.out = 200)
        seqSig = seq(min(sigma2H[keepRun,,]), max(sigma2H[keepRun,,]), length.out = 200)
        for(nb in 1:nBatch){
          for(t in keepRun){
            betaS = matrix(betaH[t,,nb], ncol=nSites, byrow=F)
            for(nT in 1:nBasis){
              betaTilde[,basis$idx[[nT]]] = t(betaS[rep(nT, length(basis$idx[[nT]])),match(y$ID ,ID)])
            }
            Btilde = basisSplines * betaTilde
            Yp = matprod_simple_triplet_matrix(Btilde,alphaH[t,,nb])
            if(nSpTcov>0){Ypp = Yp + SpTcov %*% gammaH[t,,nb]}else{Ypp = Yp}
            sigma2PostA = priors$sigma2$a0 + nObs / 2
            sigma2PostB = priors$sigma2$b0 + sum((Mu - Ypp)^2, na.rm=T) / 2
            pSigma[t,,nb] = dinvgamma(seqSig, shape = sigma2PostA, rate=sigma2PostB)
            
            iSigmaBT = kronecker(diag(nBasis),solve(exp(- phiH[t,,nb] * DistMat)))
            tauPostA = priors$tau$a0 + nSites*nBasis / 2
            tauPostB = priors$tau$b0 + t(as.matrix(betaH[t,,nb])) %*% iSigmaBT %*% as.matrix(betaH[t,,nb])/2
            pTau[t,,nb] = dinvgamma(seqTau,tauPostA, tauPostB)
        }
        }
        
        plot(seqSig,colMeans(pSigma[,,1] ,na.rm=T),
             type="l", col  = adjustcolor(color[1],alpha.f = 0.5), lwd=2, 
             main = expression(pi[sigma^2]^y), cex=0.5)
        if(nBatch>1){
          for(nB in 2:nBatch){
            points(seqSig,colMeans(pSigma[,,nB] ,na.rm=T),
                   type="l", col  = adjustcolor(color[1], alpha.f = 0.5), lwd=2)
          }
        }
        
        
        plot(seqTau,colMeans(pTau[,,1], na.rm=T),
             type="l", col  = adjustcolor(color[2],alpha.f = 0.5), lwd=2, 
             main = expression(pi[tau]^y), cex=0.5)
        if(nBatch>1){
          for(nB in 2:nBatch){
            points(seqTau,colMeans(pTau[,,nB], na.rm=T),
                   type="l", col  = adjustcolor(color[2], alpha.f = 0.5), lwd=2)
          }
        }
        
        
        if(SPTMresobj$tempRE == "corr"){
          
          plot(tauTH, type="l", col  = color[2], lwd=2, main = expression(tau[time]))
          
          plot(phiTH, type="l", col  = color[3], lwd=2, main = expression(phi[time]))
        }
        
      }
      
      if(type == "gamma"){
        if(floor(sqrt(nSpTcov)) == sqrt(nSpTcov)){
          dim.plot = c(sqrt(nSpTcov), sqrt(nSpTcov))
        }else{
          dim.plot = c(floor(sqrt(nSpTcov)+1), floor(sqrt(nSpTcov)+1))
          reduction = floor((prod(dim.plot) - nSpTcov) / dim.plot[1])
          dim.plot = dim.plot - c(reduction, 0)
        }
        par(mfrow=dim.plot, mar=c(3,3,2,1))
        
        if(is.null(color)){color = brewer.pal(12, "Set3")}
        
        # Sigma2
        pGamma = array(NA, dim = c(length(keepRun), length(keepRun),nBatch,nSpTcov))
        
        betaTilde = simple_triplet_zero_matrix(nrow = nrow(y), ncol = dim(alphaH)[2])
        Mu = y$obs
        tWW = crossprod(SpTcov, SpTcov)
        itWW = solve(tWW)
        iSigmaGamma = priors$gamma$s0 * diag(dim(gammaH)[2])
        for(nb in 1:nBatch){
          for(t in keepRun){
            betaS = matrix(betaH[t,,nb], ncol=nSites, byrow=F)
            for(nT in 1:nBasis){
              betaTilde[,basis$idx[[nT]]] = t(betaS[rep(nT, length(basis$idx[[nT]])),match(y$ID ,ID)])
            }
            Btilde = basisSplines * betaTilde
            Yp = matprod_simple_triplet_matrix(Btilde,alphaH[t,,nb])
            iSigmaObsG = 1/sigma2H[t,,nb]* tWW
            
            
            gammaPostSd =solve(iSigmaGamma + iSigmaObsG)
            gammaPostMu = gammaPostSd %*% (iSigmaObsG %*% (itWW %*% crossprod(SpTcov, (Mu - Yp))))
           
            for(j in 1:nSpTcov){
              muJ = gammaPostMu[j] - gammaPostSd[j,j]*(iSigmaGamma + iSigmaObsG)[j,-j] %*% matrix(gammaH[t,-j,nb]-gammaPostMu[-j], ncol=1)
              pGamma[t-min(keepRun)+1,,nb,j] = dnorm(gammaH[keepRun,j,nb], muJ , gammaPostSd[j,j])
            }
          }
        }
        
        for(nG in 1:nSpTcov){
        plot(sort(gammaH[keepRun,nG,1]),colMeans(pGamma[,,1,nG])[sort(gammaH[keepRun,nG,1], index.return=T)$ix],
             type="l", col  = adjustcolor(color[1],alpha.f = 0.5), lwd=2, 
             main = expression(pi[gamma]^y), cex=0.5)
        if(nBatch>1){
          for(nB in 2:nBatch){
            points(sort(gammaH[keepRun,nG,nB]),colMeans(pGamma[,,nB,nG])[sort(gammaH[keepRun,nG,nB], index.return=T)$ix],
                   type="l", col  = adjustcolor(color[1], alpha.f = 0.5), lwd=2)
          }
        }
        }
        

        
        
      }
      
      
    }
  }

"createSpTdata" <- 
  function(min.date = as.Date("01-01-2009", format = "%d-%m-%Y"),
           max.date = as.Date("31-12-2010", format = "%d-%m-%Y"),
           byDate = "month",
           nSite = 8,
           nLandUse = 2,
           nBasis = 2,
           betaTrue = NULL,
           parameters = NULL,
           missingMeasures = list(status = F, ratio = 0, type="MNAR")
           ){

    date = seq(min.date,max.date,by = byDate)
    seasonality = array(NA, dim = c(length(date),2,nLandUse))
    shift = sample(-6:6, size = nLandUse, replace = F)
    if(byDate == "month"){
      for(i in 1:nLandUse){
        seasonality[,,i] = cbind(cos(((as.numeric(format(date, "%m")) + shift[i])*2*3.146 )/ 12),sin(((as.numeric(format(date, "%m")) + shift[i]))*2*3.146 / 12)^2)
      }
      }
    if(byDate=="week"){
      for(i in 1:nLandUse){
      seasonality = cbind(cos(((as.numeric(format(date, "%W")) + shift[i])*2*3.146) / 52),sin(((as.numeric(format(date, "%W")) + shift[i])*2*3.146) / 52)^2)
      }
      }
    
    if(is.null(betaTrue)){
      betaTrue = matrix(runif((nBasis+1)*nLandUse,-2,2),ncol=nBasis+1)
    }
    if(is.null(parameters)){
      parameters = list(sigma = 0.1, tau = 0.5, phi = 1)
    }

    locations = matrix(20*runif(nSite*2),ncol=2)
    LU.comp = dirmult::rdirichlet(nSite,alpha = rep(0.5,nLandUse))
    LU.true = t(apply(LU.comp,1,function(x) rmultinom(1,1,x)))
    
    sim.raw.cov = data.frame(ID = as.factor(1:nSite), Longitude = locations[,1], Latitude = locations[,2])
    LU.df = as.data.frame(LU.comp)
    names(LU.df) <- paste0("LU",1:nLandUse)
    LU.df$ID = 1:nSite
    
    sim.raw.cov = merge(sim.raw.cov, LU.df, by = "ID")
    sim.raw.obs = data.frame(date = rep(date, each= nSite), ID = as.factor(rep(1:nSite, length(date))))
    sim.raw.obs$obs = NA
    
    X.true = LU.true
    D = as.matrix(dist(sim.raw.cov[,c("Longitude","Latitude")]))
    for(i in 1:nSite){
      sim.raw.obs$obs[sim.raw.obs$ID == i] = betaTrue[which.max(X.true[i,]),1] + 
        seasonality[,1:nBasis, which.max(X.true[i,])] %*% matrix(betaTrue[which.max(X.true[i,]),2:(2+nBasis-1)], ncol=1)
    }
    
    SpatialRF = rmvnorm(1, mean=rep(0,nSite), sigma = parameters$tau*exp(-parameters$phi*D))
    
    for(j in date){
      sim.raw.obs$obs[sim.raw.obs$date == j] = sim.raw.obs$obs[sim.raw.obs$date == j] + SpatialRF + rnorm(length(diag(D)),0,parameters$sigma)
    }
    if(missingMeasures$status){
      if(missingMeasures$type == "MNAR"){
        keepDates = sort(sample(x = 1:length(date), size = round(length(date)*(1-missingMeasures$ratio))))
        sim.raw.obs = sim.raw.obs[sim.raw.obs$date %in% date[keepDates],]}
      if(missingMeasures$type == "MAR"){
        keepObs = sort(sample(x = 1:nrow(sim.raw.obs), size = round(nrow(sim.raw.obs)*(1-missingMeasures$ratio))))
        sim.raw.obs = sim.raw.obs[keepObs,]}
    }
    #sim.raw.obs$ID = as.character(sim.raw.obs$ID)
    #sim.raw.cov$ID = as.character(sim.raw.cov$ID)
    return(list(df = sim.raw.obs, X = sim.raw.cov, coordinates = locations, luComp = LU.df, 
                betaTrue = betaTrue, grTrue = LU.true, parTrue = parameters))
  }
  

"wrap_ll" <-
  function(theta,y,kernelList,xList){
  idx.p = 0
  param = NULL
  for(i in 1:length(kernelList)){
    paramTemp=NULL
    par = attr(kernelList, "parameters")[[i]]
    idx.p = max(idx.p) + seq_len(length(par))
    for(j in 1:length(idx.p)){
      assign(par[j], (theta[idx.p][j]))
    }
    eval(parse(text=paste0("paramTemp = c(paramTemp,c(list(",par,"=",par,")))")))
    param = c(param, list(paramTemp))
  }
  
  paramList = c(param, list(sigma_noise = (theta[length(theta)])))
  if(sum(theta<0)>=1){ll = -Inf}else{ll = logl_gp(paramList, y, kernelList = kernelList, xList=xList)}
  return(ll)
}

"wrap_dll" <-
  function(theta,y,kernelList,xList){
    idx.p = 0
    param = NULL
    for(i in 1:length(kernelList)){
      paramTemp=NULL
      par = attr(kernelList, "parameters")[[i]]
      idx.p = max(idx.p) + seq_len(length(par))
      for(j in 1:length(idx.p)){
        assign(par[j], (theta[idx.p][j]))
      }
      eval(parse(text=paste0("paramTemp = c(paramTemp,c(list(",par,"=",par,")))")))
      param = c(param, list(paramTemp))
    }
    
    paramList = c(param, list(sigma_noise = (theta[length(theta)])))
    if(sum(theta<0)>=1){dll = rep(0, length(unlist(paramList))-1);}else{dll = matrix(dlogl_gp(paramList, y, kernelList = kernelList, xList=xList), nrow=1)}
    return(dll)
  }

# kernelList = list(k_seasonal)
# attr(kernelList, "name") <- c("seasonal")
# attr(kernelList, "parameters") <- list(c("q3","q4", "qs", "f"))
# attr(kernelList, "type") <- c("temporal")
# theta = c(2,1/5,1/5,1/60 ,0.2)
# paramList = parToList(theta, kernelList = kernelList)
# xList  = list(matrix(runif(20,0,10), ncol = 1))
# y = cos(xList[[1]]*4)+sqrt(abs(xList[[1]])) + rnorm(length(xList[[1]]),0,0.2)
# 
# logl_gp(paramList = paramList, y = y, kernelList = kernelList, xList = xList)
# wrap_ll(theta = theta, y = y, kernelList = kernelList, xList = xList)
# wrap_dll(theta = theta, y = y, kernelList = kernelList, xList = xList)
# dlogl_gp(paramList = paramList, y = y, kernelList = kernelList, xList = xList)
# 
# grad(func = wrap_ll, x = theta,  y = y, kernelList = kernelList, xList = xList)

"logl_gp" <- 
  function(paramList, y, kernelList,xList){
  y = c(y)
  
  Ky = K_wrap(paramList, y, kernelList,xList)
  
  ll = -1/2 * t(y) %*% solve(Ky) %*% y - 1/2 * as.numeric(determinant(Ky, logarithm = TRUE)$modulus) - length(y) /2 * log(2*pi)
  return(ll)
  }

"dlogl_gp" <- 
  function(paramList, y, kernelList,xList){
    nParam = length(unlist(paramList))-1
    y = c(y)
    
    Ky = K_wrap(paramList, y, kernelList,xList)
    
    dKy = dK_wrap(paramList, y, kernelList, xList)
    
    dll = NULL
    for(i in 1:nParam){
      dll_tmp = -1/2 * sum(diag(solve(Ky) %*% dKy[[i]])) + 1/2 * t(y) %*% solve(Ky) %*% dKy[[i]] %*% solve(Ky) %*% y
      dll = c(dll, dll_tmp)
    }
    return(dll)
  }

"K_wrap" <-
  function(paramList, y, kernelList,xList){
  y = c(y)
  Kf = matrix(0, ncol = length(y), nrow = length(y))
  for(i in 1:length(kernelList)){
    x = xList[[i]]
    if(ncol(x) == 1){
      Ktemp = outer(c(x),c(x), kernelList[[i]], param = paramList[[i]])
    }else{
      xd1 = x[rep(1:nrow(x), nrow(x)),]
      xd2 = x[rep(1:nrow(x), each = nrow(x)),]
      xdd = cbind(xd1,xd2)
      Ktemp = matrix(apply(xdd,1, function(x) kernelList[[i]](x[1:2],x[3:4],param = paramList[[i]])), ncol = nrow(x))
    }
    Kf = Kf + Ktemp
  }
  Ky =  Kf + paramList$sigma_noise * diag(length(y))
  
  return(Ky)
}

"dK_wrap" <-
  function(paramList, y, kernelList,xList){
    y = c(y)
    dKy = NULL
    for(i in 1:length(kernelList)){
      x = xList[[i]]
      diffx = outer(c(x),c(x), "-")
      Ktemp = outer(c(x),c(x), kernelList[[i]], param = paramList[[i]])
      if(attr(kernelList, "name")[i] == "long term"){
        dKf = list(1 / paramList[[i]]$q1 * Ktemp, -diffx^2*Ktemp)
      }
      if(attr(kernelList, "name")[i] == "seasonal"){
        dKf = list(1 / paramList[[i]]$q3 * Ktemp, 
                   -1/2*diffx^2*Ktemp,
                   -2*sin(3.14*paramList[[i]]$f*diffx)^2*Ktemp,
                   -2*paramList[[i]]$qs*2*3.14*diffx*sin(paramList[[i]]$f*3.14*diffx)*cos(paramList[[i]]$f*3.14*diffx)*Ktemp)
      }
      dKy =  c(dKy, dKf)
    }
    return(dKy)
  }

"optifix" <-
  function(par, fixed, fn, gr = NULL, ..., method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"), lower = -Inf, upper = Inf, control = list(), hessian = FALSE){
  force(fn)
  force(fixed)
  .npar=length(par)
  .fixValues = par[fixed]
  .parStart = par[!fixed]
  .fn <- function(par,...){
    .par = rep(NA,sum(!fixed))
    .par[!fixed] = par
    .par[fixed] = .fixValues
    fn(.par,...) }
  if(!is.null(gr)){
    .gr <- function(par,...){
      .gpar = rep(NA,sum(!fixed))
      .gpar[!fixed] = par
      .gpar[fixed] = .fixValues
      gr(.gpar,...)[!fixed] }
  }else{ .gr <- NULL }
  
  .opt = optim(.parStart,.fn,.gr,...,method=method,lower=lower,control=control,hessian=hessian)
  .opt$fullpars = rep(NA,sum(!fixed))
  .opt$fullpars[fixed]=.fixValues
  .opt$fullpars[!fixed]=.opt$par
  .opt$fixed = fixed
  return(.opt)
}

"parToList" <-
  function(theta, kernelList){
  idx.p = 0
  param = NULL
  for(i in 1:length(kernelList)){
    paramTemp=NULL
    par = attr(kernelList, "parameters")[[i]]
    idx.p = max(idx.p) + seq_len(length(par))
    for(j in 1:length(idx.p)){
      assign(par[j], (theta[idx.p][j]))
    }
    eval(parse(text=paste0("paramTemp = c(paramTemp,c(list(",par,"=",par,")))")))
    param = c(param, list(paramTemp))
  }
  
  param = c(param, list(sigma_noise = (theta[length(theta)])))
  names(param)[1:length(kernelList)] <- attr(kernelList, "name")
  return(param)
}

"k_longterm" <- 
  function(t,tp, param = list(q1=1,q2=1)){
  param$q1 * exp(-(t-tp)^2 / param$q2^2)
}

"k_seasonal" <- 
  function(t,tp,param = list(q3=1,q4=1,qs=1,f=1)){
  param$q3 * exp( -2*sin(3.14*param$f*(t-tp))^2 / param$qs^2) * exp(-1/2*(t-tp)^2 / param$q4^2)
  }

"GPpred" <-
  function(xd,x,y,param, kernel){
  
  y = matrix(y, ncol = 1)
  
  Kpx = matrix(0, nrow = nrow(xd[[1]]), ncol = nrow(x[[1]]))
  Kxx = matrix(0, nrow = nrow(x[[1]]), ncol = nrow(x[[1]]))
  Kpp = matrix(0, nrow = nrow(xd[[1]]), ncol = nrow(xd[[1]]))
  
  for(i in 1:length(kernel)){
    
    if(ncol(x[[i]])==1){
      Kpx = Kpx + outer(c(xd[[i]]),c(x[[i]]),kernel[[i]],param = param[[i]])
      Kxx = Kxx + outer(c(x[[i]]),c(x[[i]]),kernel[[i]],param = param[[i]])
      Kpp = Kpp + outer(c(xd[[i]]),c(xd[[i]]),kernel[[i]],param = param[[i]])
    }
    if(ncol(x[[i]])>1){
      xd1 = xd[[i]][rep(1:nrow(xd[[i]]), nrow(xd[[i]])),]
      xd2 = xd[[i]][rep(1:nrow(xd[[i]]), each = nrow(xd[[i]])),]
      x1 = x[[i]][rep(1:nrow(x[[i]]), nrow(x[[i]])),]
      x2 = x[[i]][rep(1:nrow(x[[i]]), each = nrow(x[[i]])),]
      
      xdd = cbind(xd1,xd2)
      xdx = cbind(xd[[i]][rep(1:nrow(xd[[i]]), nrow(x[[i]])),],x[[i]][rep(1:nrow(x[[i]]), each = nrow(xd[[i]])),])
      xx = cbind(x1,x2)
      Kpx = Kpx + matrix(apply(xdx,1, function(x) kernel[[i]](x[1:2],x[3:4],param = param[[i]])), ncol = nrow(x[[i]]))
      Kxx = Kxx +matrix(apply(xx,1, function(x) kernel[[i]](x[1:2],x[3:4],param = param[[i]])), ncol = nrow(x[[i]]))
      Kpp = Kpp + matrix(apply(xdd,1, function(x) kernel[[i]](x[1:2],x[3:4],param = param[[i]])), ncol = nrow(xd[[i]]))
    }
    
  }
  
  
  mu_p = Kpx %*% ginv(Kxx) %*% y
  sigma_p = Kpp - Kpx %*% ginv(Kxx) %*% t(Kpx)
  
  wd = list(mp = mu_p, sp = sigma_p)
  
  class(wd) <- "gppred"
  attr(wd, "type") <- attr(kernel, "type")
  
  return(wd)
}

"yHatFunction" <-
  function(alpha, beta, betaI, gamma, z, SPTMresobj){
  
  y = SPTMresobj$GibbsOut$y
  basis = SPTMresobj$GibbsOut$basis
  basisSplines = do.call(cbind,basis$splines)
  ID = SPTMresobj$GibbsOut$ID
  cov = SPTMresobj$GibbsOut$cov
  SpTcov = SPTMresobj$GibbsOut$sptcov
  
  nBasis = length(basis$idx)
  nSites = length(ID)
  nCov = ncol(cov)
  nObs = nrow(y)
   
  Xtilde = simple_triplet_zero_matrix(nrow = nObs, ncol = nSites*nCov*nBasis)
  Iter = diag(nSites)[match(y$ID, ID),]
  yHat = matrix(nrow = nrow(beta), ncol = nrow(y))
  
  for(t in 1:nrow(beta)){

  for(si in 1:nSites){
    idxC = ((si-1)*nBasis+1):(si*nBasis)
    xt = NULL
    for(jj in 1:(nBasis)){
      if(!is.null(z)){fTime = as.matrix(basisSplines[,basis$idx[[jj]]]) %*% alpha[t,basis$idx[[jj]], which.max(z[,si,t])]}
      if(is.null(z)){fTime = as.matrix(basisSplines[,basis$idx[[jj]]]) %*% alpha[t,basis$idx[[jj]]]}
      xt = cbind(xt,fTime * cov)
    }
    Xtilde[which(y$ID == ID[si]),idxC] = xt[which(y$ID == ID[si]),]
  }
  if(!is.null(z)){SpTcov = t(z[,match(y$ID,ID),t])}
  Yp = crossprod_simple_triplet_matrix(t(Xtilde), beta[t,]) + SpTcov %*% gamma[t,]
  yHat[t,] = Yp + Iter %*% betaI[t,]
  }
  return(yHat)
}

"samplePar" <- 
  function(n, m, s, constrained = F, A = NULL, e= NULL){
    sam = NULL
    if(!constrained){
      sam = rmvnorm(n, mean = m, sigma = s, method = 'chol')
    }
    if(constrained){
      for(i in 1:n){
      L = t(chol(solve(s)))
      zi = rmvnorm(1,rep(0, nrow(s)), diag(1, nrow(s)))
      v = solve(t(L),t(zi))
      xi = m + v
      Vnk = solve(s) %*% t(A)
      Wkk = A %*% Vnk
      Ukn = solve(Wkk) %*% t(Vnk)
      ci = A %*% xi - e
      sam = rbind(sam,matrix(xi - t(Ukn) %*% ci, nrow=1))
      }
    }
    return(sam)
  }

'dConstrainedPar' <-
  function(theta, m, s, log=T, A = NULL, e = NULL){
    muAlphaC = m - solve(s) %*% t(A) %*% solve(A %*% solve(s) %*% t(A)) %*% (A %*% m - e)
    sigAlphaC = solve(s) - solve(s) %*% t(A) %*% solve(A %*% solve(s) %*% t(A)) %*% A %*% solve(s)
    EVD = eigen(sigAlphaC)
    iEV = 1 / EVD$values
    iEV[EVD$values < 1e-300] = 0
    SigM = EVD$vectors %*% diag(iEV) %*% t(EVD$vectors)
    prAlphaT =  NULL
    for(i in 1:nrow(theta)){
      prAlphaT   = c(prAlphaT,-(ncol(theta)-1)/2 * log(2*pi) -
        1/2*t(theta[i,] - muAlphaC)%*% SigM %*%  (theta[i,]- muAlphaC) -
        1/2*sum(log(EVD$values[EVD$values > 1e-300])))
    }
    return(prAlphaT)
    
  }

'likelihoodSimple' <-
  function(y, yHat, sigma2){
  logL = NULL
  for(i in 1:nrow(yHat)){
    logL = c(logL, sum(dnorm(yHat[i,], y, sigma2[i], log = T)))
  }
  return(logL)
}

'margLikelihood' <-
  function(SPTMresobj, nSample = 100, keepRun = NULL, c = 1){
    
    nObs = nrow(SPTMresobj$GibbsOut$y)
    Nrun = nrow(SPTMresobj$GibbsOut$theta$betaH)
    nBasis = length(SPTMresobj$GibbsOut$basis$idx)
    priors = SPTMresobj$GibbsOut$priors
    
    if(is.null(keepRun)){
      keepRun = round(0.8*Nrun):Nrun
    }
    
    if(dim(SPTMresobj$GibbsOut$theta$alphaH)[3] == 1){

    idxMAP = which.max(SPTMresobj$GibbsOut$logPostDist$Post[keepRun,,1])  
    
    parAlpha = list(m = SPTMresobj$GibbsOut$theta$alphaH[keepRun[idxMAP],,1,1], 
                      s = apply(SPTMresobj$GibbsOut$theta$alphaH[keepRun,,1,1],2,sd))
    
    parBeta = list(m = (SPTMresobj$GibbsOut$theta$betaH[keepRun[idxMAP],,1]),
                   s = apply(SPTMresobj$GibbsOut$theta$betaH[keepRun,,1],2,sd))
    parBetaI = list(m = (SPTMresobj$GibbsOut$theta$betaIH[keepRun[idxMAP],,1]), 
                    s = apply(SPTMresobj$GibbsOut$theta$betaIH[keepRun,,1],2,sd))
    parGamma = list(m = (SPTMresobj$GibbsOut$theta$gammaH[keepRun[idxMAP],,1]), 
                    s = apply(SPTMresobj$GibbsOut$theta$gammaH[keepRun,,1],2,sd))
    parSigma2 = list(m = (SPTMresobj$GibbsOut$theta$sigma2H[keepRun[idxMAP],1,1]), s = 1)
    
    basisSplines = do.call(cbind,SPTMresobj$GibbsOut$basis$splines)
    
    A = matrix(0,nrow = nBasis, ncol = ncol(basisSplines))
    for(cA in 1:nBasis){
      A[cA,SPTMresobj$GibbsOut$basis$idx[[cA]]] = colMeans(basisSplines)[SPTMresobj$GibbsOut$basis$idx[[cA]]]
    }

    alphaTilde = samplePar(n = nSample, m = parAlpha$m, s = diag(parAlpha$s)*c, constrained = T, 
                           A = A, e = matrix(0, nrow = nBasis, ncol = 1))
    
    betaTilde = samplePar(n = nSample, m = parBeta$m, s = diag(parBeta$s)*c, constrained = F)
    betaITilde = samplePar(n = nSample, m = parBetaI$m, s = diag(parBetaI$s)*c, constrained = T, 
                           A = matrix(1, nrow=1, ncol = length(parBetaI$m)), e = matrix(0, nrow = 1, ncol = 1))
    gammaTilde = samplePar(n = nSample, m = parGamma$m, s = diag(parGamma$s)*c, constrained = F)
    sigma2Tilde = rgamma(nSample, parSigma2$m/c,parSigma2$s/c)
    
    z = NULL

    yHat = yHatFunction(alpha = alphaTilde, 
                        beta = betaTilde, 
                        betaI = betaITilde, 
                        gamma = gammaTilde, 
                        z = z, 
                        SPTMresobj = SPTMresobj)
    
    logL = likelihoodSimple(SPTMresobj$GibbsOut$y$obs, yHat, sigma2Tilde)

    dAlpha = dConstrainedPar(theta = alphaTilde, 
                             m = t(matrix(priors$alpha$m0, nrow=1, ncol = ncol(alphaTilde))), 
                             s = diag(priors$alpha$s0, length(parAlpha$m)),  
                             A = A, e = matrix(0, nrow = nBasis, ncol = 1))
    
    dPrior =mvtnorm::dmvnorm(betaTilde, mean = rep(priors$beta$m0, length(parBeta$m)),sigma = diag(priors$beta$s0, length(parBeta$m)),  log = T) + 
      mvtnorm::dmvnorm(gammaTilde, mean = rep(priors$gamma$m0, length(parGamma$m)),sigma = diag(priors$gamma$s0, length(parGamma$m)),  log = T) + 
      dConstrainedPar(theta = betaITilde, m = t(matrix(0, nrow=1, ncol = ncol(betaITilde))), s = diag(ncol(betaITilde)),  
                      A = matrix(1, nrow=1, ncol = length(parBetaI$m)), e = matrix(0, nrow = 1, ncol = 1)) + 
      dAlpha
      
      dIAlpha = dConstrainedPar(theta = alphaTilde, m = parAlpha$m, s = diag(parAlpha$s)*c,  
                                A = A, e = matrix(0, nrow = nBasis, ncol = 1))
    
    dImp = mvtnorm::dmvnorm(betaTilde, mean = parBeta$m, sigma = diag(parBeta$s)*c,  log = T) + 
      mvtnorm::dmvnorm(gammaTilde, m  =  parGamma$m, sigma = diag(parGamma$s)*c, log = T) + 
      dConstrainedPar(theta = betaITilde, m = parBetaI$m, s = diag(parBetaI$s)*c,  
                      A = matrix(1, nrow=1, ncol = length(parBetaI$m)), e = matrix(0, nrow = 1, ncol = 1)) + 
      dIAlpha
    
    df = data.frame(dPrior = dPrior, dImp = dImp, logL = logL, dPost = dPrior + logL - dImp)
    return(df)
    }
    
    if(dim(SPTMresobj$GibbsOut$theta$alphaH)[3] > 1){
      
      idxMAP = which.max(SPTMresobj$GibbsOut$logPostDist$Post[keepRun,,1])  
      
      
      
      
      parAlpha = list(m = SPTMresobj$GibbsOut$theta$alphaH[keepRun[idxMAP],,,1], 
                      s = apply(SPTMresobj$GibbsOut$theta$alphaH[keepRun,,,1], 2:3, sd))
      
      parBeta = list(m = (SPTMresobj$GibbsOut$theta$betaH[keepRun[idxMAP],,1]),
                     s = apply(SPTMresobj$GibbsOut$theta$betaH[keepRun,,1],2,sd))
      parBetaI = list(m = (SPTMresobj$GibbsOut$theta$betaIH[keepRun[idxMAP],,1]), 
                      s = apply(SPTMresobj$GibbsOut$theta$betaIH[keepRun,,1],2,sd))
      parGamma = list(m = (SPTMresobj$GibbsOut$theta$gammaH[keepRun[idxMAP],,1]), 
                      s = apply(SPTMresobj$GibbsOut$theta$gammaH[keepRun,,1],2,sd))
      parSigma2 = list(m = (SPTMresobj$GibbsOut$theta$sigma2H[keepRun[idxMAP],1,1]), s = 1)
      
      basisSplines = do.call(cbind,SPTMresobj$GibbsOut$basis$splines)
      
      A = matrix(0,nrow = nBasis, ncol = ncol(basisSplines))
      for(cA in 1:nBasis){
        A[cA,SPTMresobj$GibbsOut$basis$idx[[cA]]] = colMeans(basisSplines)[SPTMresobj$GibbsOut$basis$idx[[cA]]]
      }

        alphaTilde = array(NA, dim = c(nSample,length(parAlpha$m[,1]), dim(SPTMresobj$GibbsOut$theta$alphaH)[3]))
        for(i in 1:dim(SPTMresobj$GibbsOut$theta$alphaH)[3]){
          alphaTilde[,,i] = samplePar(n = nSample, m = parAlpha$m[,i], 
                                      s = cov(SPTMresobj$GibbsOut$theta$alphaH[keepRun,,i,1])*c, constrained = F, 
                                      A = A, e = matrix(0, nrow = nBasis, ncol = 1))
        }
      
      #  browser()
      #  A %*% SPTMresobj$GibbsOut$theta$alphaH[keepRun[idxMAP],,,1]
        
      betaTilde = samplePar(n = nSample, m = parBeta$m, s = diag(parBeta$s)*c, constrained = F)
      betaITilde = samplePar(n = nSample, m = parBetaI$m, s = diag(parBetaI$s)*c, constrained = T, 
                             A = matrix(1, nrow=1, ncol = length(parBetaI$m)), e = matrix(0, nrow = 1, ncol = 1))
      gammaTilde = samplePar(n = nSample, m = parGamma$m, s = diag(parGamma$s)*c, constrained = F)

      sigma2Tilde = rgamma(nSample, parSigma2$m/c,parSigma2$s/c)
      
      
        zProb = apply(SPTMresobj$GibbsOut$theta$zH[keepRun[idxMAP],,,1],2:3, mean)
        z = replicate(nSample,apply(t(zProb), 1, function(x) rmultinom(1,1,x)))
      
      yHat = yHatFunction(alpha = alphaTilde, 
                          beta = betaTilde, 
                          betaI = betaITilde, 
                          gamma = gammaTilde, 
                          z = z, 
                          SPTMresobj = SPTMresobj)

      logL = likelihoodSimple(SPTMresobj$GibbsOut$y$obs, yHat, sigma2Tilde)
      
        dAlpha = 0
        for(j in 1:dim(SPTMresobj$GibbsOut$theta$alphaH)[3]){
          #dAlpha = dAlpha + 
            #dConstrainedPar(theta = alphaTilde[,,j], m = t(matrix(priors$alpha$m0, nrow=1, ncol = ncol(basisSplines))), 
            #                                s = diag(priors$alpha$s0, ncol(basisSplines)),  
            #                                A = A, e = matrix(0, nrow = nBasis, ncol = 1))
          dAlpha = dAlpha + mvtnorm::dmvnorm(alphaTilde[,,j], mean = t(matrix(priors$alpha$m0, nrow=1, ncol = ncol(basisSplines))), 
                                    sigma = diag(priors$alpha$s0, ncol(basisSplines)),  log = T)
        }
      
        dBetaP = mvtnorm::dmvnorm(betaTilde, mean = rep(priors$beta$m0, length(parBeta$m)),sigma = diag(priors$beta$s0, length(parBeta$m)),  log = T)
        
      dPrior = dBetaP + 
        mvtnorm::dmvnorm(gammaTilde, mean = rep(priors$gamma$m0, length(parGamma$m)),sigma = diag(priors$gamma$s0, length(parGamma$m)),  log = T) + 
        dConstrainedPar(theta = betaITilde, m = t(matrix(0, nrow=1, ncol = ncol(betaITilde))), s = diag(ncol(betaITilde)),  
                        A = matrix(1, nrow=1, ncol = length(parBetaI$m)), e = matrix(0, nrow = 1, ncol = 1)) + 
        dAlpha
      
      
        dIAlpha = 0
        for(j in 1:dim(SPTMresobj$GibbsOut$theta$alphaH)[3]){
          # dIAlpha = dIAlpha + dConstrainedPar(theta = alphaTilde[,,j], m = parAlpha$m[,j], s = diag(parAlpha$s[,j])*c,  
          #                                     A = A, e = matrix(0, nrow = nBasis, ncol = 1))
          dIAlpha = dIAlpha + mvtnorm::dmvnorm(alphaTilde[,,j], mean = parAlpha$m[,j], 
                                      sigma = cov(SPTMresobj$GibbsOut$theta$alphaH[keepRun,,j,1])*c, log=T)
        }
      
      dImp = mvtnorm::dmvnorm(betaTilde, mean = parBeta$m, sigma = diag(parBeta$s)*c,  log = T) + 
       mvtnorm::dmvnorm(gammaTilde, m  =  parGamma$m, sigma = diag(parGamma$s)*c, log = T) + 
        dConstrainedPar(theta = betaITilde, m = parBetaI$m, s =diag(parBetaI$s)*c,  
                        A = matrix(1, nrow=1, ncol = length(parBetaI$m)), e = matrix(0, nrow = 1, ncol = 1)) + 
        dIAlpha

      df = data.frame(dPrior = dPrior, dImp = dImp, logL = logL, dPost = dPrior + logL - dImp)
      return(df)
    }
    
  }

'get_legend' <-
  function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
