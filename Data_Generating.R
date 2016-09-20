Data_Generating = function(K,seednum, plot_mode ){
  # Call libraries 
  
  library("mvtnorm")
  
  # Generate toy data 
  set.seed(seednum)
  dim = 2 ## Dimension, for visulization 
  
  # Generate parameters 
  data = c()
  Z = c()
  N = c()
  for (cls_idx in 1:K){
    n_cls = sample(50:100,1) 
    mean_cls = rnorm(dim,mean=runif(1,-100,100),sd=2)
    cov_cls = runif(2,1,5) * diag(dim)
    data_cls = rmvnorm(n_cls, mean_cls, cov_cls)
    data = rbind(data,data_cls)
    Z_cls = rep(cls_idx,n_cls)
    Z = c(Z,Z_cls)
    N = c(N,n_cls)
  }
  
  if (plot_mode == T){
    color = c("red","blue","green","black","darkgoldenrod","burlywood","aquamarine","chocolate")
    plot(data,col=color[Z])
  }
  
  return (list("data"=data, "Z"=Z, "K"=K,"N"=N))
}