library(simcausal)
library(mclogit)
library(lme4)
library(survival)
#Packages for plotting
library(ggplot2)
library(tidyverse)
library(ggridges)
library(latex2exp)

simulate_data_log<-function(sim, b0, bW, bB, n, m){
  #simulation id, cluster id and observation within cluster id
  sim_id <- rep(1:sim, each = n*m)
  cluster_id <- as.factor(rep(rep(1:m, each = n), sim))
  obs_id <- rep(rep(1:n,m), sim)
  
  #Cluster-specific effect
  u0 <- rep(rnorm(m*sim, 0, 1), each = n)
  
  eta_x <- 1 + 1*u0
  p_x <- exp(eta_x)/(1+exp(eta_x))
  
  #Covariate correalted with cluster-specific effect
  x <- rbern(n*m*sim, p_x)
  
  #Cluster-mean of covariates
  x_mean_cluster <- c()
  for (i in seq(1, sim*n*m, n)){
    x_mean_cluster <- append(x_mean_cluster, rep(mean(x[i:(i+n-1)]), n))
  }
  
  x_dev_from_mean <- x - x_mean_cluster
  
  eta <- b0 + bW*x_dev_from_mean + bB*x_mean_cluster + u0
  p <- exp(eta)/(1+exp(eta))

  y <- rbern(n*m*sim, p)
  
  df <- data.frame(y, sim_id, cluster_id, obs_id, x, x_mean_cluster, x_dev_from_mean)
  return(df)
}

#Function to estimate regression coefficients by using simulated data from simulate_data_LMM()
log <- function(sim, b0, bW, bB, n, m){
  beta0=c()
  beta1=c()
  beta1_c=c()
  betaW=c()
  betaB=c()
  
  #simulate data
  df = simulate_data_log(sim, b0, bW, bB, n, m)
  
  #for each simulation, estimate regression coefficients and store them
  for (i in 1:sim){
    log = glmer(y~ x + (1|cluster_id), data = df[df$sim_id == i,], family = binomial)
      
    log_c = clogit(y~ x +strata(cluster_id), data = df[df$sim_id == i,])
      
    log_bw= glmer(y~  x_dev_from_mean + x_mean_cluster + (1|cluster_id), data = df[df$sim_id == i,], family = binomial)
  
    beta0=append(beta0, fixef(log)[1])
    beta1=append(beta1, fixef(log)[2])
    beta1_c=append(beta1_c,coefficients(log_c))
    betaW=append(betaW, fixef(log_bw)[2])
    betaB=append(betaB, fixef(log_bw)[3])
  }
  betas=data.frame(beta0, beta1, beta1_c, betaW, betaB)
  return(betas)
}


m=50
n=5
b0=-1
bW=1
bB=1
sim=200

set.seed(3456)

#Only random intercept
betas_log=log(sim, b0, bW, bB, n, m)

mean=apply(betas_log[],2,mean)
quant=apply(betas_log[],2,quantile,probs=c(0.025,0.975))
df_qm=data.frame(mean ,"2.5"=quant[1,],"97.5"=quant[2,])
df_qm

histogram_beta<-ggplot(gather(subset(betas_log, select = -c(1,5))), aes(x = value, y = key, group = key, fill= key))+
  geom_density_ridges2(stat="binline", bins=30, alpha=0.6, scale=5, color = "azure4")+
  scale_fill_viridis(name="", discrete=TRUE, option="plasma")+
  geom_vline(aes(xintercept=bW, color="True beta"), color= "black", linewidth=0.9, linetype="dashed", alpha=0.8)+
  labs(x=TeX(r'(Estimated $\beta_1$ and $\beta_W$)'),y=TeX(r'(Count)'))+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=15),legend.text = element_text(size=14),legend.title = element_text(size=16))+
  annotate(geom = "text", x=0.7, y=6.2, label=TeX(r'($\beta_W=1$)'), size=5, color="black")
histogram_beta



