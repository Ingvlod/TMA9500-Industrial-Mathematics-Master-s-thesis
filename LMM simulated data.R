library(lme4)
library(plm)
#Packages for plotting
library(ggplot2)
library(tidyverse)
library(ggridges)
library(latex2exp)
library(viridis)

simulate_data_LMM_BW<- function(sim, b0, bW, bB, n, m, slopes=FALSE){
  #simulation id, cluster id and observation within cluster id
  sim_id <- rep(1:sim, each = n*m)
  cluster_id <- as.factor(rep(rep(1:m, each = n), sim))
  obs_id <- rep(rep(1:n,m), sim)
  
  #Cluster-specific effect
  u0 <- rep(rnorm(m*sim, 0, 1), each = n)
  
  #Covariate correalted with cluster-specific effect
  x <- rnorm(n*m*sim, mean = 1+5*u0, sd = 1)
  
  #Cluster-mean of covariates
  x_mean_cluster <- c()
  
  for (i in seq(1, sim*n*m, n)){
    x_mean_cluster <- append(x_mean_cluster, rep(mean(x[i:(i+n-1)]), n))
  }
  
  x_dev_from_mean <- x - x_mean_cluster
  
  if(slopes){
    u1 <- rep(rnorm(m*sim, 0, 1), each = n)
    mean_y <- b0 + bW*x_dev_from_mean + bB*x_mean_cluster + u0 + u1*x 
  }else{
    mean_y <- b0 + bW*x_dev_from_mean + bB*x_mean_cluster + u0
  }
  
  y <- rnorm(n*m*sim, mean=mean_y, sd=1)
  
  df <- data.frame(y, sim_id, cluster_id, obs_id, x, x_mean_cluster, x_dev_from_mean)
  return(df)
}

#Function to transform data in order to implement CMLE for LMM
transform_data <- function(df, sim, n, m){
  x <- c()
  y <- c()
  x_dev_from_mean <-c()
  A_t = t(contr.poly(rep(1, each=n)))
  for(k in 1:sim){
    for (i in 1:m){
      x=append(x, A_t%*%df[df$cluster_id==i & df$sim_id==k, ]$x)
      y=append(y, A_t%*%df[df$cluster_id==i & df$sim_id==k, ]$y)
      x_dev_from_mean=append(x_dev_from_mean, A_t%*%df[df$cluster_id==i & df$sim_id==k, ]$x_dev_from_mean)
    }
  }
  n_t = n-1
  sim_id <- rep(1:sim, each = n_t*m)
  cluster_id <- as.factor(rep(rep(1:m, each = n_t), sim))
  obs_id <- rep(rep(1:n_t,m), sim)
  df_t <- data.frame(y, sim_id, cluster_id, obs_id, x, x_dev_from_mean)
  return (df_t)
}

#Function to estimate regression coefficients by using simulated data from simulate_data_LMM_BW()
lmm <- function(sim, b0, bW, bB, n, m, slopes=FALSE){
  beta0=c()
  beta1=c()
  beta1_c=c()
  beta1_FE=c()
  betaW=c()
  betaB=c()
  betaW_lm=c()
  
  #simulate data
  df = simulate_data_LMM_BW(sim, b0, bW, bB, n, m, slopes)
  df_transformed = transform_data(df, sim, n, m)
  pdata = pdata.frame(df, index="cluster_id")
  
  pb <- txtProgressBar(min = 1, max = sim, style = 3)
  
  #for each simulation, estimate regression coefficients and store them
  for (i in 1:sim){
    if(slopes){
      lmm = lmer(y~ x + (1+x|cluster_id), data = df[df$sim_id == i,], REML = TRUE)
      
      lm = lm (y~x_dev_from_mean + x_mean_cluster, data = df[df$sim_id == i,])

      lmm_c = lmer(y~ 0 + x_dev_from_mean + (0 + x_dev_from_mean|cluster_id), data = df_transformed[df_transformed$sim_id == i,], REML = TRUE)
      
      lmm_FE = plm(y~x|x, data = pdata[pdata$sim_id == i,], model = "within")
      
      lmm_bw= lmer(y~ x_dev_from_mean + x_mean_cluster + (1+x|cluster_id), data = df[df$sim_id == i,], REML = TRUE)

      beta1_c=append(beta1_c,fixef(lmm_c)[1])
    }else{
      lmm = lmer(y~ x + (1|cluster_id), data = df[df$sim_id == i,], REML = TRUE)
      
      lm = lm (y~x_dev_from_mean + x_mean_cluster, data = df[df$sim_id == i,])
      
      lmm_c = lm(y~ 0 + x, data = df_transformed[df_transformed$sim_id == i,])
      
      lmm_FE = plm(y~ x, data = pdata[pdata$sim_id == i,], model = "within")
      
      lmm_bw= lmer(y~  x_dev_from_mean + x_mean_cluster + (1|cluster_id), data = df[df$sim_id == i,], REML = TRUE)

      beta1_c=append(beta1_c,coefficients(lmm_c)[1])
    }
    beta0=append(beta0, fixef(lmm)[1])
    beta1=append(beta1, fixef(lmm)[2])
    beta1_FE=append(beta1_FE, coefficients(lmm_FE)[1])
    betaW=append(betaW, fixef(lmm_bw)[2])
    betaB=append(betaB, fixef(lmm_bw)[3])
    betaW_lm=append(betaW_lm, coef(lm)[2])
    
    setTxtProgressBar(pb, i)

  }
  betas=data.frame(beta0, beta1, beta1_c, beta1_FE, betaW, betaB, betaW_lm)
  return(betas)
}


m=50
n=3
b0=-1
bW=1
bB=1
sim=100


#Random intercept only LMM
set.seed(3456)

betas_lmm=lmm(sim, b0, bW, bB, n, m, slopes=FALSE)

mean=apply(betas_lmm[],2,mean)
quant=apply(betas_lmm[],2,quantile,probs=c(0.025,0.975))
df_qm=data.frame(mean ,"2.5"=quant[1,],"97.5"=quant[2,])
df_qm

histogram_beta<-ggplot(gather(subset(betas_lmm, select = -c(1,6))), aes(x = value, y = key, group = key, fill= key))+
  geom_density_ridges2(stat="binline", bins=30, alpha=0.6, scale=5, color = "azure4")+
  scale_fill_viridis(name="Fitting methods", discrete=TRUE, option="plasma", labels=c("GLM","Between-within","Conditional likelihood", "Fixed effects estimator" , "Mixed effects approach"), breaks=c("betaW_lm","betaW", "beta1_FE" , "beta1_c", "beta1"))+
  geom_vline(aes(xintercept=bW, color="True beta"), color= "black", linewidth=0.9, linetype="dashed", alpha=0.8)+
  labs(x=TeX(r'(Estimated $\beta_1$ and $\beta_W$)'),y=TeX(r'(Count)'))+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=13),legend.text = element_text(size=12),legend.title = element_text(size=14), axis.text.y = element_blank())+
  annotate(geom = "text", x=0.95, y=5.75, label=TeX(r'($\beta_W=1$)'), size=5, color="black")
histogram_beta

#Random intercept and slopes LMM
betas_lmm_sl=lmm(sim, b0, bW, bB, n, m, slopes=TRUE)

mean=apply(betas_lmm_sl[],2,mean)
quant=apply(betas_lmm_sl[],2,quantile,probs=c(0.025,0.975))
df_qm=data.frame(mean ,"2.5"=quant[1,],"97.5"=quant[2,])
df_qm

histogram_beta<-ggplot(gather(subset(betas_lmm_sl, select = -c(1,6))), aes(x = value, y = key, group = key, fill= key))+
  geom_density_ridges2(stat="binline", bins=25, alpha=0.6, scale=5, color = "azure4")+
  scale_fill_viridis(name="Fitting methods", discrete=TRUE, option="plasma", labels=c("Between-within", "Fixed effects estimator" ,"Conditional likelihood", "Mixed effects approach"), breaks=c("betaW", "beta1_FE" , "beta1_c", "beta1"))+
  geom_vline(aes(xintercept=bW, color="True beta"), color= "black", linewidth=0.9, linetype="dashed", alpha=0.8)+
  labs(x=TeX(r'(Estimated $\beta_1$ and $\beta_W$)'),y=TeX(r'(Count)'))+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=13),legend.text = element_text(size=12),legend.title = element_text(size=14), axis.text.y = element_blank())+
  annotate(geom = "text", x=1.15, y=7.25, label=TeX(r'($\beta_W=1$)'), size=5, color="black")
histogram_beta




