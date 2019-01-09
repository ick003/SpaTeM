"SPTMData" <- 
  function(df.obs, tempBasis = "re", tempPeriod = "%m", nSplines = 3, splinesType = "poly"){
  
  SiteID = levels(df.obs$ID)
  
  t = df.obs$date
  
  BB = NULL
  idxT = NULL
  j=0
  maxIdx = 0
  
  
  if(splinesType == "penalized"){
    for(i in tempPeriod){
      j = j +1 
        if(j == 1){
          dat = data.frame(tTemp = as.numeric(as.Date(paste("2014-",format(t,"%m-%d"), sep = ""))))
          BBT <- smoothCon(s(tTemp,k=nSplines[j]),data=dat,knots=NULL)[[1]]
          #BBT = splineFun(tTemp, df = nSplines[j],
          #                Boundary.knots = c(min(tTemp)-diff(range(tTemp))/(1*nSplines[j]),max(tTemp)+diff(range(tTemp))/(1*nSplines[j])))
        }
        if(j == 2){
          dat = data.frame(tTemp = as.numeric(t))
          BBT <- smoothCon(s(tTemp,k=nSplines[j]),data=dat,knots=NULL)[[1]]
          #BBT = splineFun(t, df = nSplines[j], 
          #                Boundary.knots = c(min(t)-diff(range(t))/(1*nSplines[j]),max(t)+diff(range(t))/(1*nSplines[j])))
        }
        if(j == 3){
          BBT = matrix(rep(1, length(t)), ncol=1)
        }

      BB = c(BB,list(BBT))
      idxt = 1:ncol(BBT$X) + maxIdx #(ncol(BB) - ncol(BBT))
      maxIdx = max(idxt)
      idxT = c(idxT,list(idxt))
    }
  }else{
    if(splinesType == "poly"){
      splineFun = function(...) bs(intercept = T, ...)
    }
    if(splinesType == "cubic"){
      splineFun = function(...) ns(intercept  =T, ...)
    }
  for(i in tempPeriod){
    j = j +1 
    if(tempBasis == "bs"){
      if(j == 1){
        tTemp = as.Date(paste("2014-",format(t,"%m-%d"), sep = ""))
        BBT = splineFun(tTemp, df = nSplines[j],
                 Boundary.knots = c(min(tTemp)-diff(range(tTemp))/(1*nSplines[j]),max(tTemp)+diff(range(tTemp))/(1*nSplines[j])))
      }
      if(j == 2){
        BBT = splineFun(t, df = nSplines[j], 
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
  }
  
  sptm.data = list(obs.data = df.obs, basis = BB, list.idx = idxT, tempBasis = tempBasis, tempPeriod = tempPeriod, splinesType = splinesType)
  class(sptm.data) <- 'sptm'
  return(sptm.data)
}
"SPTModel" <- 
  function(df.sptm, df.lu = NULL, coordinates, cov = NULL, SpTcov = NULL){
  
  
  sptm.model = list(data = df.sptm$obs.data,
                    basis = list(splines = df.sptm$basis, idx = df.sptm$list.idx, tempBasis = df.sptm$tempBasis, tempPeriod = df.sptm$tempPeriod, type = df.sptm$splinesType),
                    lu = df.lu,
                    coord = coordinates,
                    covariates = cov,
                    SpTcov = SpTcov)
  
  class(sptm.model) <- 'stpmod'
  return(sptm.model)
  
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
