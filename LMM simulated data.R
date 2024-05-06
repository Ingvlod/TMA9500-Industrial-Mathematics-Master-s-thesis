library(lme4)
library(plm)
#Packages for plotting
library(ggplot2)
library(tidyverse)
library(ggridges)
library(latex2exp)

simulate_data_LMM<- function(sim, b0, b1, n, m, slopes=FALSE, LMM=FALSE, logistic=FALSE){
  #simulation id, cluster id and observation within cluster id
  sim_id <- rep(1:sim, each = n*m)
  cluster_id <- as.factor(rep(rep(1:m, each = n), sim))
  obs_id <- rep(rep(1:n,m), sim)
  
  #Cluster-specific effect
  u0 <- rep(rnorm(m*sim, 0, 1), each = n)
  
  #Covariate correalted with cluster-specific effect
  x <- rnorm(n*m*sim, mean = 10+u0, sd = 0.4)
  
  #Cluster-mean of covariates
  x_mean_cluster <- c()
  for (i in seq(1, sim*n*m, n)){
    x_mean_cluster <- append(x_mean_cluster, rep(mean(x[i:(i+n-1)]), n))
  }
  
  x_dev_from_mean <- x - x_mean_cluster
  
  if(slopes){
    u1 <- rep(rnorm(m*sim, 1, 1), each = n)
    mean_y <- b0 + b1*x + u0 + u1*x 
  }else{
    mean_y <- b0 + b1*x + u0
  }
  
  y <- rnorm(n*m*sim, mean=mean_y, sd=1)
  
  df <- data.frame(y, sim_id, cluster_id, obs_id, x, x_mean_cluster, x_dev_from_mean)
  return(df)
}

simulate_data_LMM_BW<- function(sim, b0, bW, bB, n, m, slopes=FALSE){
  #simulation id, cluster id and observation within cluster id
  sim_id <- rep(1:sim, each = n*m)
  cluster_id <- as.factor(rep(rep(1:m, each = n), sim))
  obs_id <- rep(rep(1:n,m), sim)
  
  #Cluster-specific effect
  u0 <- rep(rnorm(m*sim, 0, 1), each = n)
  
  #Covariate correalted with cluster-specific effect
  x <- rnorm(n*m*sim, mean = 10+u0, sd = 0.4)
  
  #Cluster-mean of covariates
  x_mean_cluster <- c()
  for (i in seq(1, sim*n*m, n)){
    x_mean_cluster <- append(x_mean_cluster, rep(mean(x[i:(i+n-1)]), n))
  }
  
  x_dev_from_mean <- x - x_mean_cluster
  
  if(slopes){
    u1 <- rep(rnorm(m*sim, 1, 1), each = n)
    mean_y <- b0 + bW*x_dev_from_mean + bB*x_mean_cluster + u0 + u1*x 
  }else{
    mean_y <- b0 + bW*x_dev_from_mean + bB*x_mean_cluster + u0
  }
  
  y <- rnorm(n*m*sim, mean=mean_y, sd=1)
  
  df <- data.frame(y, sim_id, cluster_id, obs_id, x, x_mean_cluster, x_dev_from_mean)
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
lmm <- function(sim, b0, bW, bB, n, m, slopes=FALSE){
  beta0=c()
  beta1=c()
  beta1_c=c()
  beta1_FE=c()
  betaW=c()
  betaB=c()
  #simulate data
  #df = simulate_data_LMM(sim, b0, b1, n, m, slopes)
  df = simulate_data_LMM_BW(sim, b0, bW, bB, n, m, slopes)
  df_transformed = transform_data(df, sim, n, m)
  pdata = pdata.frame(df, index="cluster_id")
  #for each simulation, estimate regression coefficients and store them
  for (i in 1:sim){
    if(slopes){
      lmm = lmer(y~ x + (1+x|cluster_id), data = df[df$sim_id == i,])

      lmm_c = lmer(y~ 0 + x+ (0 + x|cluster_id), data = df_transformed[df_transformed$sim_id == i,])
      
      lmm_FE = plm(y~ x|x, data = pdata[pdata$sim_id == i,], model = "within")
      
      lmm_bw= lmer(y~ 0 + x_dev_from_mean + x_mean_cluster + (1+x|cluster_id), data = df[df$sim_id == i,])

      beta1_c=append(beta1_c,fixef(lmm_c)[1])
    }else{
      lmm = lmer(y~ x + (1|cluster_id), data = df[df$sim_id == i,])
      
      lmm_c = lm(y~ 0 + x, data = df_transformed[df_transformed$sim_id == i,])
      
      lmm_FE = plm(y~ x, data = pdata[pdata$sim_id == i,], model = "within")
      
      lmm_bw= lmer(y~ 0 + x_dev_from_mean + x_mean_cluster + (1|cluster_id), data = df[df$sim_id == i,])

      beta1_c=append(beta1_c,coefficients(lmm_c)[1])
    }
    beta0=append(beta0, fixef(lmm)[1])
    beta1=append(beta1, fixef(lmm)[2])
    beta1_FE=append(beta1_FE, coefficients(lmm_FE)[1])
    betaW=append(betaW, fixef(lmm_bw)[1])
    betaB=append(betaB, fixef(lmm_bw)[2])

  }
  betas=data.frame(beta0, beta1, beta1_c, beta1_FE, betaW, betaB)
  return(betas)
}


m=10
n=3
b0=-1
bW=3
bB=8
sim=100

df=simulate_data_LMM_BW(sim, b0, 5, -10, n, m, slopes=TRUE)
ggplot(df, aes(x=x, y=y, color=cluster_id, shape=cluster_id))+geom_point()


lmm(sim, b0, bW, bB, n, m, slopes=FALSE)
lmm(sim, b0, bW, bB, n, m, slopes=TRUE)
#Random intercept only LMM

betas=lmm(sim, b0, bW, bB, n, m, slopes=FALSE)

histogram_beta<-ggplot(gather(subset(betas, select = -c(1,6))), aes(x = value, y = key, group = key, fill= key))+
  geom_density_ridges2(stat="binline", bins=50, alpha=0.6, scale=5, color = "azure4")+
  #scale_fill_manual(name="", values=c("darkolivegreen2","brown1","cyan3"))+
  scale_fill_brewer(name="Case", palette="Set1")+
  geom_vline(aes(xintercept=bW, color="True beta"), color= "#bf5252", linewidth=0.9, linetype="dashed", alpha=0.8)+
  labs(x=TeX(r'(Estimated $\beta_1$ and $\beta_W$)'),y=TeX(r'(Count)'))+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=15),legend.text = element_text(size=14),legend.title = element_text(size=16))+
  annotate(geom = "text", x=8, y=7, label=TeX(r'($\beta_W=?$)'), size=5, color="#bf5252")
histogram_beta

#Random intercept and slopes LMM

betas=lmm(sim, b0, bW, bB, n, m, slopes=TRUE)

histogram_beta<-ggplot(gather(subset(betas, select = -c(1,6))), aes(x = value, y = key, group = key, fill= key))+
  geom_density_ridges2(stat="binline", bins=50, alpha=0.6, scale=5, color = "azure4")+
  #scale_fill_manual(name="", values=c("darkolivegreen2","brown1","cyan3"))+
  scale_fill_brewer(name="Case", palette="Set1")+
  geom_vline(aes(xintercept=bW, color="True beta"), color= "#bf5252", linewidth=0.9, linetype="dashed", alpha=0.8)+
  labs(x=TeX(r'(Estimated $\beta_1$ and $\beta_W$)'),y=TeX(r'(Count)'))+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=15),legend.text = element_text(size=14),legend.title = element_text(size=16))+
  annotate(geom = "text", x=8, y=7, label=TeX(r'($\beta_W=?$)'), size=5, color="#bf5252")
histogram_beta




