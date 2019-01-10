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
          if(SPTMresobj$GibbsOut$basis$type == "penalized"){
            
            for(k in 1:nComp){
              for(i in 1:nBasis){
                t0 = y$date
                if(basis$tempPeriod[i] == "%m"){
                  t0 = as.Date(paste("2014-",format(t0,"%m-%d"), sep = ""))
                }
                
                if(i == 1){deltaT = "day"}
                if(i == 2){deltaT = "week"}
                
                xPred = seq(min(t0), max(t0), deltaT)
                
                xBasis = PredictMat(basis$splines[[i]],data.frame(tTemp=as.numeric(xPred)))
                
                #xBasis = predict(basis$splines[[i]],xPred)
                
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
          }else{
            for(k in 1:nComp){
            for(i in 1:nBasis){
              t0 = y$date
              if(basis$tempPeriod[i] == "%m"){
                t0 = as.Date(paste("2014-",format(t0,"%m-%d"), sep = ""))
              }
              
              if(i == 1){deltaT = "day"}
              if(i == 2){deltaT = "month"}
              xPred = seq(min(t0), max(t0), deltaT)
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