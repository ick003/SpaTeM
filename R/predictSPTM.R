"predictSPTM" <- 
  function(SPTMresobj, transform = "none", posterior = F, newdata = NULL, keepRun = NULL){
    
    if(transform == "none"){
      fT = function(x){x}
    }
    if(transform == "log"){
      fT = function(x){exp(x)}
    }
    
    
    alphaH = SPTMresobj$GibbsOut$theta$alphaH
    thetaGPH = SPTMresobj$GibbsOut$theta$thetaGPH
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
    splinesType = SPTMresobj$GibbsOut$basis$type
    
    nBasis = length(basis$idx)
    nSites = length(ID)
    nSplines = ncol(basisSplines)
    nCov = ncol(cov)
    nObs = nrow(y)
    Nrun = nrow(betaH)
    nHyper = 3 + 2 * as.numeric(SPTMresobj$tempRE == "corr")
    nComp = ncol(SPTMresobj$GibbsOut$priors$pi$alpha0)
    nSpTCov = ncol(SpTcov)
    if(is.null(keepRun)){
      keepRun = round(Nrun/2):Nrun
    }
    Xtilde = simple_triplet_zero_matrix(nrow = nObs, ncol = nSites*nCov*nBasis)
    Iter = diag(nSites)[match(y$ID, ID),]
    
    if(SPTMresobj$tempRE == "gp"){
      
      if(is.null(newdata)){
        
        newBetaH = NULL
        
        
        if(posterior){
          
          YpredH = colMeans(yHat[keepRun,,1])
          YpredHN = yHat[keepRun,,1]
          
          YpredHN = YpredHN + rnorm(length(keepRun)*nObs, 0, mean(sqrt(sigma2H[keepRun,,1])))
          
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
        #browser()
        
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
        

        # Predict temporal signal #
        

        Mu = y$obs
        
        d11 = as.matrix(dist(coords[,2:3]))
        d12 = as.matrix(dist(rbind(coords[,2:3],newLoc)))[1:nSites,(nSites+1):(nSites + nNewSites)]
        d22 = as.matrix(dist(newLoc))
        
        newBetaH = matrix(NA,length(keepRun), ncol =  nBasis*nNewSites)
        
        colnames(newBetaH) <- rep(rownames(newLoc), each = nBasis)
        
        MuA = Mu - SpTcov %*% gammaH[keepRun[1],,1]
        
        for(tt in 2:length(keepRun)){
          
          t = keepRun[tt]
          zCurr = matrix(zH[keepRun[tt-1],,,1], nrow=nComp)
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
          
          MuCorr = diag(betaH[t-1,match(y$ID, levels(y$ID))*nBasis,1]^(-1)) %*% (MuA - betaH[t-1,match(y$ID, levels(y$ID))*nBasis-1,1])
          
          for(si in 1:nNewSites){
            idxZ = which(y$ID %in% ID[which(apply(zCurr,2,which.max) == newZ[si])])
            if(length(idxZ)==0){next}
            param = parToList(c(thetaGPH[t,,1],sigma2H[t-1,1,1]), SPTMresobj$GibbsOut$kernelList)
            xList = list(matrix(y$date[idxZ], ncol=1),matrix(y$date[idxZ], ncol=1))
            xPred = list(matrix(newTime, ncol=1),matrix(newTime, ncol=1))

            tGP = GPpred(xd = xPred, x = xList,y = c(MuCorr[idxZ]), param = param, kernel = kernelList)
                
            ft = tGP$mp

            fTime = ft
            xt = cbind(1,fTime * newCov)
            idxC = ((si-1)*nBasis+1):(si*nBasis)
            Xtilde[which(nID == newID),] = xt[which(nID == newID),]
            
          }
          
          #Yp = crossprod_simple_triplet_matrix(t(Xtilde), betaH[t,,1])
          
          Yp = crossprod_simple_triplet_matrix(t(Xtilde), newBetaH[tt,]) + gammaH[t-1,newZ,1]#SpTcov %*% gammaH[t-1,,1] #
          yHat = Yp #+ betaIH[t,match(newID, SPTMresobj$GibbsOut$ID),1]
          YppredHN[tt,] = yHat + rnorm(length(newTime),0,(sigma2H[t,,1]))
          YppredH[tt,] = yHat
          #browser()
          MuA = Mu - SpTcov %*% gammaH[t,,1]
        }
        
        YpredH = colMeans(YppredH, na.rm=T)
        YpredH025 = apply(YppredHN,2,quantile, 0.025, na.rm=T)
        YpredH975 = apply(YppredHN,2,quantile, 0.975, na.rm=T ) 
        YpredHN = colMeans(YppredHN, na.rm=T)
        
      }
    }else{
      if(is.null(newdata)){
      
      newBetaH = NULL
      
      
      if(posterior){
        
        YpredH = colMeans(yHat[keepRun,,1])
        YpredHN = yHat[keepRun,,1]
        
        YpredHN = YpredHN + rnorm(length(keepRun)*nObs, 0, mean(sqrt(sigma2H[keepRun,,1])))
        
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
      
      #browser()
      
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
      
      if(splinesType == "penalized"){
        for(i in 1:length(tempPeriod)){
          j = j +1 
          #browser()
          if(j == 1){
            dat = data.frame(tTemp = as.numeric(as.Date(paste("2014-",format(newTime,"%m-%d"), sep = ""))))
            #BBT <- smoothCon(s(tTemp,k=nSplines[j]),data=dat,knots=NULL)[[1]]
            
            BBT <- PredictMat(basis$splines[[i]],dat)
            #BBT = splineFun(tTemp, df = nSplines[j],
            #                Boundary.knots = c(min(tTemp)-diff(range(tTemp))/(1*nSplines[j]),max(tTemp)+diff(range(tTemp))/(1*nSplines[j])))
          }
          if(j == 2){
            dat = data.frame(tTemp = as.numeric(newTime))
            BBT <- PredictMat(basis$splines[[i]],dat)
            #BBT <- smoothCon(s(tTemp,k=nSplines[j]),data=dat,knots=NULL)[[1]]
            #BBT = splineFun(t, df = nSplines[j], 
            #                Boundary.knots = c(min(t)-diff(range(t))/(1*nSplines[j]),max(t)+diff(range(t))/(1*nSplines[j])))
          }
          if(j == 3){
            BBT = matrix(rep(1, length(t)), ncol=1)
          }
          
          newBB = c(newBB,list(BBT))
        }
        
        #newsplines = do.call(cbind,lapply(newBB, function(x) x$X))
        newsplines = do.call(cbind,newBB)
      }else{
      
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
      }
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
    }
    return(list(Yp = YpredH, YpwN = YpredHN, newB = newBetaH, CI025 = YpredH025, CI975 = YpredH975))
    
  }