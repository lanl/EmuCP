# univariate variance inflation of samples x
variance_inflation = function(x, delta = .5, n_dens_est = 1000, bounds=c(Inf,Inf)){
  d = density(x, from = bounds[1], to = bounds[2], n=n_dens_est)
  d$y = d$y^delta # power up density
  d$bw
  cdf = cumsum(d$y)/sum(d$y)
  u = runif(length(x))
  samples = approx(cdf, d$x, xout = u, yleft=0, yright = 1)$y
  return(samples)
}

generate_cartesian_product = function (n, p)
{
  args_list <- replicate(p, 1:n, simplify = FALSE)
  result <- do.call(expand.grid, args_list)
  return(result)
}

# bivariate KS between two sample sets
bvks = function(x, y, n_pts=100){
  nx = nrow(x)
  ny = nrow(y)
  p = ncol(x)
  if(p != ncol(y)){
    stop("x and y must be matrices with same number of columns.")
  }
  if(p != 2){
    stop('x and ')
  }
  bounds = matrix(NA, nrow=2, ncol=p)
  for(i in 1:p){
    bounds[,i] = range(c(x[,i], y[,i]))
  }

  ks = -Inf
  for(ii in 1:n_pts){
    for(jj in 1:n_pts){
      z = c(bounds[1,1] + ((ii - 1) / (n_pts - 1)) * diff(bounds[,1]),
             bounds[1,2] + ((jj - 1) / (n_pts - 1)) * diff(bounds[,2]))

      Fx = sum((x[,1] <= z[1]) & (x[,2] <= z[2]))/nx
      Fy = sum((y[,1] <= z[1]) & (y[,2] <= z[2]))/ny
      curr = abs(Fx-Fy)
      if(curr > ks) ks = curr
    }
  }
  return(ks)
}

# univariate sample based estimator of wasserstein distance
wasserstein_distance = function(x, y, p = 1, n_points = 1000) {
  u = seq(0, 1, length.out = n_points)
  qx = quantile(x, probs = u, type = 8)
  qy = quantile(y, probs = u, type = 8)
  Wp = (mean(abs(qx - qy)^p))^(1/p)
  return(Wp)
}

# univariate wasserstein distance from samples x to known normal parameterized by (mu, sd)
wasserstein_distance_normal = function(x, mu, sd, p = 1, n_points = 1000) {
  u = seq(0,1, length.out = n_points+2)[-c(1,n_points+2)]
  qx = quantile(x, probs = u)
  qy = qnorm(u, mu, sd)
  Wp = (mean(abs(qx - qy)^p))^(1/p)
  return(Wp)
}

# sample based estimator of LP distance between 2 univariate CDF's
cdf_dist_lp = function (x, y, p = 1, n_pts = 30){
  if (!is.matrix(x) | !is.matrix(y)) {
    x = matrix(x, ncol = 1)
    y = matrix(y, ncol = 1)
  }
  nx = nrow(x)
  ny = nrow(y)
  q = ncol(x)
  if (q != ncol(y)) {
    stop("x and y must be matrices with same number of columns.")
  }
  bounds = matrix(NA, nrow = 2, ncol = q)
  for (i in 1:q) {
    bounds[, i] = range(c(x[, i], y[, i]))
  }
  indx = generate_cartesian_product(n_pts, q)
  ks = -Inf
  res = 0
  for (i in 1:nrow(indx)) {
    flagx = rep(TRUE, nx)
    flagy = rep(TRUE, ny)
    for (j in 1:q) {
      ii = indx[i, j]
      z = bounds[1, j] + ((ii - 1)/(n_pts - 1)) * diff(bounds[,j])
      flagx = flagx & (x[, j] <= z)
      flagy = flagy & (y[, j] <= z)
    }
    Fx = sum(flagx)/nx
    Fy = sum(flagy)/ny
    curr = abs(Fx - Fy)
    res = res + curr^p/n_pts^q
    ks = max(ks, curr)
  }
  if (p == Inf) {
    return(ks)
  }
  return(res^(1/p))
}

# helper function which defines some useful functions based on the selected Fdist
# returns negative log likelihood function
#         random sample function
#         function to initialize parameter estimates based on a sample from Fdist
Fdist_to_func = function(Fdist){
  switch(Fdist,
         weibull={
           nloglik_F = nloglik_weibull
           ralpha = rweibull
           init_func = get_init_weibull
         },
         gamma={
           nloglik_F = nloglik_gamma
           ralpha = rgamma
           init_func = get_init_gamma
         },
         gamma3={
           nloglik_F = nloglik_gamma3
           ralpha = sample_gamma3
           init_func = get_init_gamma3
         },
         beta={
           nloglik_F = nloglik_beta
           ralpha = rbeta
           init_func = get_init_beta
         },
         betacopula={
           nloglik_F = nloglik_beta
           ralpha = sample_cop
           init_func = get_init_beta
         }
  )
  return(list(nloglik_F=nloglik_F,
              ralpha=ralpha,
              init_func=init_func))
}

# negative log likelihood functions
nloglik_weibull = function(par,data){
  if(any(par<=0)){
    return(Inf)
  }
  -sum(dweibull(data,par[1],par[2],log = T))
}
nloglik_gamma = function(par,data){
  if(any(par<=0)){
    return(Inf)
  }
  -sum(dgamma(data,par[1],par[2],log = T))
}
nloglik_gamma3 = function(params, data) {
  alpha = params[1]   # Shape parameter (α)
  beta  = params[2]   # Rate parameter (β)
  theta = params[3]   # Threshold parameter (θ)

  # Parameter constraints
  if (alpha <= 0 || beta <= 0 || theta >= min(data)) {
    return(-Inf)  # Invalid parameters
  }

  y = data - theta
  if (any(y <= 0)) {
    return(-Inf)  # Data points must be greater than theta
  }

  # Calculate log-likelihood
  # log_lik = length(data) * (alpha * log(beta) - lgamma(alpha)) +
  #   (alpha - 1) * sum(log(y)) -
  #   beta * sum(y)
  loglik = sum(dgamma(y,alpha,beta,log = T))

  return(-loglik)
}
nloglik_beta = function(par,data){
  if(any(par<=0)){
    return(Inf)
  }
  -sum(dbeta(data,par[1],par[2],log = T))
}

# general purpose MLE function based on nlm() optimizer
get_mle = function(lp,par,data,hessian=F){
  OPT = nlm(lp,par,hessian = hessian,data=data)
  return(OPT)
}

# parameter initialization functions to help get_mle with a good starting location
get_init_weibull = function(x){
  weibull_cv = function(k) {
    sqrt(base::gamma(1 + 2/k) / (base::gamma(1 + 1/k)^2) - 1)
  }
  # Function to find k that matches the sample CV
  find_k = function(cv, lower = 0.1, upper = 100000) {
    uniroot(function(k) weibull_cv(k) - cv, lower = lower, upper = upper)$root
  }
  mean_x = mean(x)
  sd_x = sd(x)
  cv = sd_x / mean_x
  k_start = find_k(cv)
  lambda_start = mean_x / base::gamma(1 + 1/k_start)
  return(c(k_start,lambda_start))
}
get_init_beta = function(x){
  mean_x = mean(x)
  var_x = var(x)
  alpha_mom = (mean_x / sqrt(var_x))^2
  beta_mom = mean_x / var_x
  return(c(alpha_mom,beta_mom))
}
get_init_gamma = function(x){
  S2 = var(x)
  mx = mean(x)
  shape = mx^2/S2
  scale = S2/mx
  rate = 1/scale
  return(c(shape,rate))
}
get_init_gamma3 = function(x){
  theta = min(x)
  x = x - min(x)
  S2 = var(x)
  mx = mean(x)
  shape = mx^2/S2
  scale = S2/mx
  rate = 1/scale
  return(c(shape,rate,theta))
}

# custom sampling functions
sample_gamma3 = function(n, alpha, beta, theta) {
  # Ensure parameters are valid
  if(alpha <= 0 || beta <= 0) {
    stop("Shape (alpha) and Scale (beta) parameters must be positive.")
  }

  # Generate standard Gamma samples and shift by theta
  samples = theta + rgamma(n, shape = alpha, rate = beta)

  return(samples)
}
sample_cop = function(n, copula, beta_fits){
  p = length(copula)
  Corr = matrix(1,p,p)
  Corr[upper.tri(Corr)] = Corr[lower.tri(Corr)] = copula
  # standard normal samples with Correlation matrix defined by the copula
  sample_u = pnorm(matrix(mvnfast::rmvn(n,rep(0,p),Corr),nrow=n))
  # transform U back to [0,1]
  sample_x = c()
  for(i in 1:p){
    sample_x = cbind(sample_x,qbeta(sample_u[,i],beta_fits[[i]][1],beta_fits[[i]][2]))
  }
  return(sample_x)
}

adaptive_mcmc = function (logllh_func, init_pars, prop_sigma = NULL,
                          nsamp = 10000, nburn = 1000,
                          adapt_after = 100, # start adaptation after this iteration
                          adapt_every = 50,  # adapt every
                          adapt_pct = .5,    # percentage of previous samples to use in adaptation
                          adapt_stop = nburn,# iteration to stop adaptation
                          verbose = FALSE, ...)
{
  if(is.null(prop_sigma)){
    if(length(init_pars) != 1){
      # estimate a good proposal covariance using optim
      fit = optim(init_pars, logllh_func, control = list(fnscale = -1),
                   hessian = TRUE, ...)
      prop_sigma = solve(-fit$hessian)
    }
    else{
      # start with 10% sd
      prop_sigma = matrix((.1*init_pars)^2)
    }
  }
  p = length(init_pars)
  eps = 1e-8
  prop_sigma = prop_sigma + diag(eps,p)
  par_curr = init_pars
  llh_curr = logllh_func(par_curr, ...)
  par_store = matrix(nrow = nsamp, ncol = p)
  accept = numeric(nsamp)
  print_i = floor(seq(nsamp/10, nsamp, length.out = 10))
  for(i in 1:nsamp){
    par_prop = c(mvnfast::rmvn(1, par_curr, prop_sigma))
    llh_prop = logllh_func(par_prop, ...)
    log_accept_prob = llh_prop - llh_curr
    if(is.nan(log_accept_prob))
      log_accept_prob = -Inf
    if(log(runif(1)) < log_accept_prob){
      par_curr = par_prop
      llh_curr = llh_prop
      accept[i] = 1
    }
    par_store[i,] = par_curr
    if(i > adapt_after & i%%adapt_every == 0 & i < adapt_stop){
      len = floor(i * adapt_pct):i
      N = length(len)
      prop_sigma = cov(par_store[len,,drop=F]) + diag(eps,p)
    }
    if(verbose & i %in% print_i){
      cat('Sampling:', i/nsamp * 100, '%,  acceptance rate: ',round(mean(accept[1:i]),3),'\n')
    }
  }
  par_store = par_store[(nburn+1):nsamp,]
  accept = accept[(nburn+1):nsamp]
  accept_rate = mean(accept)

  returns = list(par_store = par_store, prop_sigma = prop_sigma, acceptance_rate = accept_rate)
  return(returns)
}
