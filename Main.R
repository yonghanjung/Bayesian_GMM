options(error = recover) # setting the error option

library("mvtnorm")
source("Data_Generating.R")

seed_num = 123456
K = 4
D = 2 

Data = Data_Generating(K,seed_num,F)
data = Data$data 
Z_true = Data$Z
N = sum(Data$N)

# Setting prior of NIW 
m0 = t(rep(0,2))  # prior of mu in NIW 
m0 = t(m0)
k0 = 1
S0 = diag(2)
v0 = 5

# Initialize Z 
Z = sample(1:K, N, replace = T)
N_z = table(Z)
alpha = 10 
iter_MCMC = 100

for (MCMC_idx in 1:iter_MCMC){
  for (i in 1:N){
    xi = data[i,]
    zi_val = Z[i]
    Z_remove = Z[-i]
    N_removal_z = c()
    for (k in 1:K){
      N_removal_z = c(N_removal_z, sum(Z_remove==k))
    }
    Z_post = c()
    for (k in 1:K){
      Nk_i = N_removal_z[k]
      # Calculate cond.prior z_i 
      logprob_cond_prior_zi_k = log((Nk_i + alpha/K )/(N+alpha-1))
      
      # Calculate likelihood 
      ## values for computing posterior NIW 
      data_k = data[Z_remove %in% k,]
      # print(list("length_data_k"=dim(data_k)[1],"Nk_i"=N_removal_z[k]))
      kNk = k0 + Nk_i 
      vNk = v0 + Nk_i
      X_k = scale(data_k,scale=T)
      if (length(data_k)/D > 1){
        mNk = (k0*m0 + colSums(data_k))/kNk
      }
      else if (length(data_k)/D == 1){
        mNk = (k0*m0 + data_k)/kNk
        data_k = t(data_k)
      }
      
      else{
        log_post = logprob_cond_prior_zi_k
        next
      }
      
      S_k = (t(data_k) %*% data_k) 
      SNk = S0 + S_k + k0*m0%*%t(m0) - kNk * mNk %*% t(mNk)
      log_likelihood = dmvt(xi,delta=mNk,sigma = (kNk + 1)/(kNk * (vNk-D+1))*SNk, df = vNk-D+1, log = T )
      # calculate post (log)
      log_post = logprob_cond_prior_zi_k + log_likelihood
      Z_post = c(Z_post,exp(log_post))
    }
    # normalize Z_post 
    Z_post = Z_post / sum(Z_post)
    # sample new_k
    new_k = c(1:K)[rmultinom(1,1,Z_post) %in% 1]
    Z[i] = new_k
    
    # print(list("z_length"=length(Z),"mcmc_idx"=MCMC_idx,"data_idx"=i))
  }
  print(MCMC_idx)
}


