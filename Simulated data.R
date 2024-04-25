library(lme4)
library(ggplot2)
library(plm)

simulate_data_LMM<- function(sim, b0, b1, n, m, slopes=FALSE, LMM=FALSE, logistic=FALSE){
  #simulation id, cluster id and observation within cluster id
  sim_id <- rep(1:sim, each = n*m)
  cluster_id <- as.factor(rep(rep(1:m, each = n), sim))
  obs_id <- rep(rep(1:n,m), sim)
  
  #Cluster-specific effect
  u0 <- rep(rnorm(m*sim, 0, 1), each = n)
  
  #Covariate correalted with cluster-specific effect
  x <- rnorm(n*m*sim, mean = 10+u0, sd = 0.4)
  
  if(slopes){
    u1 <- rep(rnorm(m*sim, 1, 1), each = n)
    mean_y <- b0 + b1*x + u0 + u1*x 
  }else{
    mean_y <- b0 + b1*x + u0
  }
  
  y <- rnorm(n*m*sim, mean=mean_y, sd=1)
  
  df <- data.frame(y, sim_id, cluster_id, obs_id, x)
  return(df)
}

transform_data <- function(df, sim, n, m){
  x <- c()
  y <- c()
  A_t = t(contr.poly(rep(1, each=n)))
  for(k in 1:sim){
    for (i in 1:m){
      x=append(x, A_t%*%df[df$cluster_id==i & df$sim_id==k, ]$x)
      y=append(y, A_t%*%df[df$cluster_id==i & df$sim_id==k, ]$y)
    }
  }
  n_t = n-1
  sim_id <- rep(1:sim, each = n_t*m)
  cluster_id <- as.factor(rep(rep(1:m, each = n_t), sim))
  obs_id <- rep(rep(1:n_t,m), sim)
  df_t <- data.frame(y, sim_id, cluster_id, obs_id, x)
  return (df_t)
}

#Function to estimate regression coefficients by using simulated data from simulate_data_LMM()
lmm <- function(sim, b0, b1, n, m, slopes=FALSE){
  beta0=c()
  beta1=c()
  beta0_c=c()
  beta1_c=c()
  beta0_W=c()
  beta1_W=c()
  #simulate data
  df = simulate_data_LMM(sim, b0, b1, n, m, slopes)
  pdata = pdata.frame(df, index="cluster_id")
  #for each simulation, estimate regression coefficients and store them
  for (i in 1:sim){
    if(slopes){
      lmm = lmer(y~ x + (1+x|cluster_id), data = df[df$sim_id == i,])
      lmm_c = plm(y~ x|x, data = pdata[pdata$sim_id == i,], model = "within")
      
      beta1_c=append(beta0_c,coefficients(lmm_c)[1])
      beta1_c=append(beta0_c,coefficients(lmm_c)[2])
      
    }else{
      lmm = lmer(y~ x + (1|cluster_id), data = df[df$sim_id == i,])
      df_transformed = transform_data(df, sim, n, m)
      lmm_c = lm(y~ x, data = df_transformed[df_transformed$sim_id == i,])
    }
    beta0=append(beta0,fixef(lmm)[1])
    beta1=append(beta1,fixef(lmm)[2])
    
    beta1FE=append(beta1FE,coefficients(FEmodel_omitted)[1])
  }
  betas=data.frame(beta0, beta1, beta1FE)
  return(betas)
}

m=5
n=3
b0=-1
b1=1
sim=2

df = simulate_data_LMM(sim, b0, b1, n, m, slopes=FALSE)

df_t = transform_data(df, sim, n, m)

df_t
ggplot(df, aes(x=x, y=y, color=cluster_id, shape=cluster_id))+geom_point()
lmm(sim, b0, b1, n, m, slopes=FALSE)
lmm(sim, b0, b1, n, m, slopes=TRUE)

