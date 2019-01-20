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
           mixture = TRUE, N.run = 10000,Blocks =NULL,
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
        res = gibbSamplerG(priors, y, basis = basis, ID = levels(y$ID), coords = coords, cov = cov,SpTcov = SpTcov,Bds = NULL,Blocks =NULL,
                           N.run = N.run, debug, print.res = print.res, nBatch = nBatch,
                           parallel = parallel, nCluster = nCluster)
      }
      
      if(tempRE == 'corr'){
        res = gibbSamplerS2(priors, y, basis = basis, ID = levels(y$ID), coords = coords, cov = cov,Blocks =NULL,
                            N.run = N.run, debug, print.res = print.res,
                            parallel = parallel, nCluster = nCluster)
      }
      
      
      RET = list(GibbsOut = res, model = model, tempRE = tempRE)
      
    }
    
    if(model == 'simpleMixture'){
      if(tempRE == 'notcorr'){
        res = gibbSamplerG(priors, y, basis = basis, ID = levels(y$ID), coords = coords, cov = cov,SpTcov = SpTcov,Bds = NULL,Blocks =NULL,
                           N.run = N.run, debug, print.res = print.res, nBatch = nBatch,
                           parallel = parallel, nCluster = nCluster)
      }
      if(tempRE == 'gp'){
        res = gibbSamplerGP(priors, y, basis = basis, ID = levels(y$ID), coords = coords, cov = cov,SpTcov = SpTcov,Bds = NULL,Blocks =NULL,
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
      
      res = gibbSamplerG(priors, y, basis = basis, ID = levels(y$ID), coords = coords, cov = cov,SpTcov = SpTcov,Bds= Bds,Blocks = Blocks,
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

