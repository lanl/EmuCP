# some useful emulators we might want easy access to
require(BASS)
ecp_bass = function(x,y,...){
  bass(x,y,...)
}
ecp_bass_PCA = function(x,y,...){
  bassPCA(x,y,...)
}
ecp_bass_mean = function(model,xx){
  yy = predict(model)
  mu = apply(yy,2,mean)
  return(mu)
}
ecp_bass_sample = function(model,xx){
  yy = predict(model)
  # return single sample of the predictive mean
  return(yy[dim(yy)[1],])
}
