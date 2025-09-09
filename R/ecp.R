#' Emulate the conditional posterior
#'
#' Implements the ECP algorithm, which emulates the conditional posterior distribution to approximate the cut posterior distribution
#'
#' @param gamma_train samples of \gamma for which the conditional posterior of \alpha is sampled. These samples are used as inputs to the emultor.
#' @param gamma_pred samples of \gamma used in phase 2 of the ECP algorithm. These samples define the set of conditional posterior distributions of alpha that are used for the final approximation to the cut posterior
#' @param lp_alpha log posterior of \alpha|\gamma,y. The first argument should be \alpha, and additional arguments are given through '...'
#' @param p dimension of \alpha
#' @param init values of \alpha used to initialize MCMC chains targeting the conditional posterior
#' @param m number of samples drawn from each approximate conditional posterior in the final approximation to the cut posterior
#' @param Fdist distribution used to estimate the conditional posterior. Default is 'normal', which generalizes to the multivariate normal for p>1
#' @param emulator emulation algorithm used for estimating the conditional posterior parameters (gp_mean, gp_sample, custom)
#' @param em_fit fit function for emulator only used when emulator='custom'. First input should be X matrix (inputs) and second should be y vector (outputs)
#' @param em_pred Prediction function for emulator only used when emulator='custom'. Should return numeric vector of predictions at input matrix X.
#' @param nsamp number of MCMC samples drawn from each conditional posterior in phase 1
#' @param nburn number of burn in samples for each conditional posterior in phase 1
#' @param nboot number of bootstrap samples for parameter estimation
#' @param ... additional arguments passed to lp_alpha
#' @return Description of the return value.
#' @examples
#' An example problem can be found in examples/diamond_in_box
#' @export
ecp = function(gamma_train,gamma_pred,lp_alpha,
               p = 1, init = rep(0,p), m=1,
               Fdist = 'laplace',emulator='gp_mean',em_fit=NULL,em_pred=NULL,
               nsamp=15000, nburn=10000, nboot=5,...){
  get_psi = function(gamma){
    n_cov_pars = p*(p-1)/2
    if(Fdist == 'normal'){
      psi = matrix(NA,nrow=nboot,ncol=p+p+n_cov_pars) # p means, p variances, n_cov_pars covariances
    } else if(Fdist == 'laplace'){
      psi = matrix(NA,nrow=1,ncol=p+p+n_cov_pars)
    } else if(Fdist %in% c('gamma','beta','weibull')){
      psi = matrix(NA,nboot,2) # 2 parameter distributions, p=1 so no covariance
    } else if(Fdist == 'betacopula'){
      psi = matrix(NA,nboot,2*p + p*(p-1)/2)
    } else{
      psi = matrix(NA,nboot,3) # 3 parameter gamma (shifted gamma)
    }
    r = ncol(psi)

    if(Fdist=='laplace'){
      OPT = optim(init,lp_alpha,hessian=TRUE,control=list(fnscale=-1),
                  method = ifelse(p>1,'Nelder-Mead','Brent'),
                  lower=ifelse(p>1,-Inf,-1e8),upper=ifelse(p>1,Inf,1e8),
                  gamma=gamma,...)
      psi[1,1:p] = OPT$par
      if(p>1){
        cholSigma = chol(chol2inv(chol(-OPT$hessian)))
        psi[1,(p+1):(2*p)] = log(diag(cholSigma))
        psi[1,(2*p+1):r] = cholSigma[upper.tri(cholSigma)]
      } else{
        Sigma = 1/-OPT$hessian # this is just a variance
        psi[1,p+1] = log(Sigma)
      }
    } else{
      chain = adaptive_mcmc(lp_alpha, init,
                            prop_sigma = diag(max(.01,.1^2*init),p),
                            nsamp = nsamp,nburn = nburn, verbose = FALSE,
                            gamma=gamma, ...)
      if(!is.matrix(chain$par_store))
        chain$par_store = as.matrix(chain$par_store)

      if(Fdist=='normal'){
        bs.id = 1:nrow(chain$par_store)
        for(i in 1:nboot){
          if(nboot>1)
            bs.id = sample(1:nrow(chain$par_store),nrow(chain$par_store),replace = TRUE)
          psi[i,1:p] = apply(chain$par_store[bs.id,,drop=F],2,mean) # mean vector for alpha
          if(p>1){
            cholSigma = chol(cov(chain$par_store[bs.id,,drop=F])) # covariance of alpha
            psi[i,(p+1):(2*p)] = log(diag(cholSigma))
            psi[i,(2*p+1):r] = cholSigma[upper.tri(cholSigma)]
          } else {
            psi[i,p+1] = log(var(chain$par_store[bs.id]))
          }
        }
      } else if(Fdist %in% c('weibull','gamma','gamma3','beta')){
        bs.id = 1:nrow(chain$par_store)
        for(i in 1:nboot){
          if(nboot>1)
            bs.id = sample(1:nrow(chain$par_store),nrow(chain$par_store),replace = TRUE)
          init = FdistFuncs$init_func(chain$par_store[bs.id,,drop=F])
          if(Fdist != 'gamma3'){
            psi[i,] = log(get_mle(FdistFuncs$nloglik_F,init,data = chain$par_store[bs.id,,drop=F])$estimate) # all parameters for must be > 0
          } else{
            tmp = get_mle(FdistFuncs$nloglik_F,init,data = chain$par_store[bs.id,,drop=F])$estimate
            psi[i,] = c(log(tmp[1:2]),tmp[3])
          }
        }
      } else if(Fdist == 'betacopula'){
        bs.id = 1:nrow(chain$par_store)
        for(i in 1:nboot){
          if(nboot>1)
            bs.id = sample(1:nrow(chain$par_store),nrow(chain$par_store),replace = TRUE)

          # estimate marginal betas
          for(j in 1:p){
            init = FdistFuncs$init_func(chain$par_store[bs.id,j,drop=F])
            psi[i,((2*j)-1):(2*j)] = log(get_mle(FdistFuncs$nloglik_F,init,data = chain$par_store[bs.id,i,drop=F])$estimate) # all parameters for must be > 0
          }
          # estimate correlations
          u = qnorm(chain$par_store[bs.id,,drop=F])
          C = cor(u)
          psi[i,(2*j+1):ncol(psi)] = C[upper.tri(C)]
        }
      } else{
        stop(paste0('No implementation for F~',Fdist))
      }
    }
    returns = list(psi=psi)
    if(Fdist != 'laplace')
      returns$alpha_chain=chain$par_store
    return(returns)
  }

  L = nrow(gamma_train)
  M = nrow(gamma_pred)
  if(Fdist == 'laplace')
    nboot=1
  if(Fdist %in% c('gamma','gamma3','beta','weibull')){
    if(p>1){
      stop(paste0('F ~ ',Fdist,' is only available for p=1. Use normal or laplace for p>1.'))
    }
    FdistFuncs = Fdist_to_func(Fdist)
  }

  n_cov_pars = p*(p-1)/2
  if(Fdist %in% c('normal','laplace')){
    psi = matrix(NA,nrow=0,ncol=p+p+n_cov_pars) # p means, p variances, n_cov_pars covariances
    psi_mean = matrix(NA,nrow=L,2*p+n_cov_pars)
  } else if(Fdist %in% c('gamma','beta','weibull')){
    psi = matrix(NA,0,2) # 2 parameter distributions, p=1 so no covariance
    psi_mean = matrix(NA,nrow=L,2)
  } else{
    psi = matrix(NA,0,3) # 3 parameter families
    psi_mean = matrix(NA,nrow=L,3)
  }
  r = ncol(psi)
  gamma_train_bs = gamma_train[rep(1:nrow(gamma_train),each=nboot),,drop=F]

  # Step 1. Learn alpha | gamma_train
  # Step 2. Estimate laplace distribution
  phase_1_chains = list()
  time = proc.time()
  for(j in 1:L){
    tmp = get_psi(gamma_train[j,])
    if(Fdist != 'laplace')
      phase_1_chains[[j]] = tmp$alpha_chain
    psi = rbind(psi,tmp$psi)#,nboot,p,Fdist,lp_alpha,...))
    psi_mean[j,] = apply(tmp$psi,2,mean)
  }
  phase_1_time = proc.time() - time
  # Step 3. Fit emulator: gamma_train->distributional parameters
  # Step 4. Predict emulator at gamma_pred
  # currently we are using the mean but we should really be sampling
  time = proc.time()
  if(emulator %in% c('gp_mean','gp_sample')){
    mult = rep(nboot,L)
    gp_list = lapply(1:r, function(i) hetGP::mleHomGP(X = list(X0=gamma_train,Z0=psi_mean[,i],mult=mult), Z = psi[,i]))
    gp_pred = lapply(1:r, function(i) predict(gp_list[[i]], gamma_pred))
    # gp_list = lapply(1:r, function(i) fit_scaled(y=psi[,i],inputs = gamma_train_bs))
    # gp_pred = lapply(1:r, function(i) predictions_scaled(gp_list[[i]],gamma_pred,m=25,joint = F,predvar = T))
    if(emulator == 'gp_sample'){
      gp_pred = sapply(1:r, function(i) rnorm(length(gp_pred[[i]]$mean),mean=gp_pred[[i]]$mean,sd=sqrt(gp_pred[[i]]$sd2 + gp_pred[[i]]$nugs)))
    } else{
      gp_pred = sapply(1:r, function(i) gp_pred[[i]]$mean)
    }
  } else if(emulator=='custom'){
    gp_list = lapply(1:r, function(i) em_fit(gamma_train_bs, psi[,i]))
    gp_pred = sapply(1:r, function(i) em_pred(gp_list[[i]], gamma_pred))
  } else{
    stop('Emulator type is not implemented')
  }
  phase_2_time = proc.time() - time

  # Step 5: sample from M predicted F's
  time = proc.time()
  if(Fdist %in% c('normal','laplace')){
    if(p>1){
      alpha = array(dim=c(M,p,m))
      for(i in 1:M){
        sigma = diag(exp(gp_pred[i,(p+1):(2*p)]))
        if(p>1){
          sigma[upper.tri(sigma)] = gp_pred[i,(2*p+1):r]
        }
        alpha[i,,] = mvnfast::rmvn(m,mu=gp_pred[i,1:p],sigma=sigma,isChol = T)
      }
    } else{
      alpha = array(dim=c(M,m))
      for(i in 1:M){
        alpha[i,] = rnorm(m,mean=gp_pred[i,1:p],sd = sqrt(exp(gp_pred[i,(p+1):(2*p)])))
      }
    }
  } else if(Fdist %in% c('gamma','beta','weibull')){
    alpha = array(dim=c(M,m))
    for(i in 1:M){
      alpha[i,] = FdistFuncs$ralpha(m,exp(gp_pred[i,1]),exp(gp_pred[i,2])) # 2 parameter family emulated on log-scale
    }
  } else{
    # 3 parameter gamma
    alpha = array(dim=c(M,m))
    for(i in 1:M){
      alpha[i,] = sample_gamma3(m,exp(gp_pred[i,1]),exp(gp_pred[i,2]),gp_pred[i,3]) # 2 parameter family emulated on log-scale
    }
  }
  dim(alpha) = c(M*m,p)
  sample_time = proc.time() - time
  return(list(rcut=alpha,time=list(phase_1 = phase_1_time,
                                   phase_2 = phase_2_time,
                                   sample_a = sample_time,
                                   total = phase_1_time+phase_2_time+sample_time),
              phase_1_chains=phase_1_chains,
              psi_hat = psi,
              psi_hat_bs_mean = psi_mean,
              emulators = gp_list,
              emulator_pred = gp_pred,
              gamma_train=gamma_train,
              gamma_pred=gamma_pred))
}

#' Emulate the conditional posterior
#'
#' Implements the ECP algorithm with provided samples from the conditional posterior
#'
#' @param gamma_train samples of \gamma for which the conditional posterior of \alpha is sampled. These samples are used as inputs to the emultor.
#' @param gamma_pred samples of \gamma used in phase 2 of the ECP algorithm. These samples define the set of conditional posterior distributions of alpha that are used for the final approximation to the cut posterior
#' @param alpha_train_list conditional samples \alpha|\gamma for each \gamma in gamma_train
#' @param p dimension of \alpha
#' @param init starting values for Fdist parameters
#' @param m number of samples drawn from each approximate conditional posterior in the final approximation to the cut posterior
#' @param Fdist distributional assumption for \alpha|\gamma
#' @param emulator emulation algorithm used for predicting the conditional posterior parameters at gamma_pred
#' @param em_fit fit function for emulator only used when emulator='custom'. First input should be X matrix (inputs) and second should be y vector (outputs)
#' @param em_pred Prediction function for emulator only used when emulator='custom'. Should return numeric vector of predictions at input matrix X.
#' @param alpha_transform transformation applied to \alpha prior to fitting Fdist. If for example \alpha \in [0,1], 'logit' or 'probit' would be appropriate
#' @return Description of the return value.
#' @export
ecp2 = function(gamma_train,gamma_pred,alpha_train_list,
               p = ncol(alpha_train_list[[1]]), init = rep(0,p), m=1,
               alpha_transform = 'identity',
               Fdist = 'normal',emulator='gp_mean',
               em_fit=NULL,em_pred=NULL,
               nsamp=15000, nburn=10000, nboot=5,...){
  get_psi = function(gamma,alpha){
    if(Fdist == 'normal'){
      n_cov_pars = p*(p-1)/2
      psi = matrix(NA,nrow=nboot,ncol=p+p+n_cov_pars) # p means, p variances, n_cov_pars covariances
    } else if(Fdist %in% c('gamma','beta','weibull')){
      psi = matrix(NA,nboot,2) # 2 parameter distributions, p=1 so no covariance
    } else if(Fdist == 'betacopula'){
      psi = matrix(NA,nboot,2*p + p*(p-1)/2)
    } else{
      psi = matrix(NA,nboot,3) # 3 parameter gamma
    }
    # total number of parameters need to simulate from Fdist
    r = ncol(psi)

    if(Fdist=='normal'){
      bs.id = 1:nrow(alpha)
      for(i in 1:nboot){
        if(nboot>1)
          bs.id = sample(1:nrow(alpha),nrow(alpha),replace = TRUE)
        psi[i,1:p] = apply(alpha[bs.id,,drop=F],2,mean) # mean vector for alpha
        if(p>1){
          # emulate the cholesky factors of the covariance matrix
          cholSigma = chol(cov(alpha[bs.id,,drop=F]))
          psi[i,(p+1):(2*p)] = log(diag(cholSigma))
          psi[i,(2*p+1):r] = cholSigma[upper.tri(cholSigma)]
        } else {
          psi[i,p+1] = log(var(alpha[bs.id]))
        }
      }
    } else if(Fdist %in% c('weibull','gamma','gamma3','beta')){
      bs.id = 1:nrow(alpha)
      for(i in 1:nboot){
        if(nboot>1)
          bs.id = sample(1:nrow(alpha),nrow(alpha),replace = TRUE)
        init = FdistFuncs$init_func(alpha[bs.id,,drop=F])
        if(Fdist != 'gamma3'){
          psi[i,] = log(get_mle(FdistFuncs$nloglik_F,init,data = alpha[bs.id,,drop=F])$estimate) # all parameters for must be > 0
        } else{
          tmp = get_mle(FdistFuncs$nloglik_F,init,data = alpha[bs.id,,drop=F])$estimate
          psi[i,] = c(log(tmp[1:2]),tmp[3])
        }
      }
    } else if(Fdist == 'betacopula'){
      bs.id = 1:nrow(alpha)
      for(i in 1:nboot){
        if(nboot>1)
          bs.id = sample(1:nrow(alpha),nrow(alpha),replace = TRUE)

        # estimate marginal betas
        u <- matrix(NA_real_, nrow = nrow(alpha), ncol = p)
        for(j in 1:p){
          init = FdistFuncs$init_func(alpha[bs.id,j,drop=F])
          beta_pars = log(get_mle(FdistFuncs$nloglik_F,init,data = alpha[bs.id,j,drop=F])$estimate) # all parameters for must be > 0
          psi[i,((2*j)-1):(2*j)] = beta_pars
          u[,j] = pbeta(alpha[bs.id,j,drop=F],exp(beta_pars[1]),exp(beta_pars[2]))
        }
        # estimate correlations
        eps = 1e-8
        u[u <= 0] = eps
        u[u >= 1] = 1 - eps
        z = qnorm(u)
        C = cor(z)
        psi[i,(2*p+1):ncol(psi)] = C[upper.tri(C)]
      }
    } else{
      stop(paste0('No implementation for F~',Fdist))
    }
    return(list(psi=psi))
  }

  L = nrow(gamma_train)
  M = nrow(gamma_pred)
  if(Fdist %in% c('gamma','gamma3','beta','weibull')){
    if(p>1){
      stop(paste0('F ~ ',Fdist,' is only available for p=1. Use normal for p>1.'))
    }
    FdistFuncs = Fdist_to_func(Fdist)
  }
  if(Fdist == 'betacopula'){
    FdistFuncs = Fdist_to_func(Fdist)
  }

  if(Fdist %in% c('normal','betacopula')){
    r = 2*p + p*(p-1)/2
  } else if(Fdist %in% c('gamma','beta','weibull')){
    r = 2
  } else{
    r = 3
  }

  psi = matrix(NA,0,r)
  psi_mean = matrix(NA,L,r)
  gamma_train_bs = gamma_train[rep(1:nrow(gamma_train),each=nboot),,drop=F]

  # Step 1. Learn alpha | gamma_train
  # Already done, passed in as alpha_train_list
  # Step 2. Estimate distributional parameters
  time = proc.time()
  for(j in 1:L){
    tmp = get_psi(gamma_train[j,],alpha_train_list[[j]])
    if(any(is.na(tmp$psi))){
      hello = 1
    }
    psi = rbind(psi,tmp$psi)#,nboot,p,Fdist,lp_alpha,...))
    psi_mean[j,] = apply(tmp$psi,2,mean)
  }
  phase_1_time = proc.time() - time
  # Step 3. Fit emulator: gamma_train->distributional parameters
  # Step 4. Predict emulator at gamma_pred
  # currently we are using the mean but we should really be sampling
  time = proc.time()
  if(emulator %in% c('gp_mean','gp_sample')){
    mult = rep(nboot,L)
    gp_list = lapply(1:r, function(i) hetGP::mleHomGP(X = list(X0=gamma_train,Z0=psi_mean[,i],mult=mult), Z = psi[,i]))
    gp_pred = lapply(1:r, function(i) predict(gp_list[[i]], gamma_pred))
    if(emulator == 'gp_sample'){
      gp_pred = sapply(1:r, function(i) rnorm(length(gp_pred[[i]]$mean),mean=gp_pred[[i]]$mean,sd=sqrt(gp_pred[[i]]$sd2 + gp_pred[[i]]$nugs)))
    } else{
      gp_pred = sapply(1:r, function(i) gp_pred[[i]]$mean)
    }
  } else if(emulator=='custom'){
    # user should specify in em_pred if a sample from the emulator or the emulator mean is given
    gp_list = lapply(1:r, function(i) em_fit(gamma_train_bs, psi[,i]))
    gp_pred = sapply(1:r, function(i) em_pred(gp_list[[i]], gamma_pred))
  } else{
    stop('Emulator type is not implemented')
  }
  phase_2_time = proc.time() - time

  # Step 5: sample from M predicted F's
  time = proc.time()
  if(Fdist == 'normal'){
    if(p>1){
      alpha = array(dim=c(M,p,m))
      for(i in 1:M){
        sigma = diag(exp(gp_pred[i,(p+1):(2*p)]))
        if(p>1){
          sigma[upper.tri(sigma)] = gp_pred[i,(2*p+1):r]
        }
        alpha[i,,] = mvnfast::rmvn(m,mu=gp_pred[i,1:p],sigma=sigma,isChol = T)
      }
    } else{
      alpha = array(dim=c(M,m))
      for(i in 1:M){
        alpha[i,] = rnorm(m,mean=gp_pred[i,1:p],sd = sqrt(exp(gp_pred[i,(p+1):(2*p)])))
      }
    }
  } else if(Fdist %in% c('gamma','beta','weibull')){
    alpha = array(dim=c(M,m))
    for(i in 1:M){
      alpha[i,] = FdistFuncs$ralpha(m,exp(gp_pred[i,1]),exp(gp_pred[i,2])) # 2 parameter family emulated on log-scale
    }
  } else if(Fdist == 'betacopula'){
    alpha = array(dim=c(M,p,m))
    for(i in 1:M){
      beta_fits = lapply(1:p, function(j) exp(gp_pred[i,((2*j)-1):(2*j)]))
      rho = gp_pred[i,(2*p+1):r]
      alpha[i,,] = FdistFuncs$ralpha(m,rho,beta_fits)
    }
  } else{
    # 3 parameter gamma
    alpha = array(dim=c(M,m))
    for(i in 1:M){
      alpha[i,] = sample_gamma3(m,exp(gp_pred[i,1]),exp(gp_pred[i,2]),gp_pred[i,3]) # 2 parameter family emulated on log-scale
    }
  }
  dim(alpha) = c(M*m,p)
  return(list(rcut=alpha,time=list(phase_1 = phase_1_time,
                                   phase_2 = phase_2_time,
                                   total = phase_1_time+phase_2_time),
              psi_hat = psi,
              psi_hat_bs_mean = psi_mean,
              emulators = gp_list,
              emulator_pred = gp_pred,
              gamma_train=gamma_train,
              gamma_pred=gamma_pred,
              alpha_train=alpha_train_list))
}
#' @param Fdist distributional assumption for \alpha|\gamma
