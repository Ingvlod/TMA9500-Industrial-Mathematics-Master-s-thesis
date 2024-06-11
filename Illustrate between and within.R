library(lme4)
#Packages for plotting
library(ggplot2)
library(viridis)
library(latex2exp)

simulate_data_BW<- function( b0, bW, bB, n, m, slopes=FALSE, randomized=FALSE, nonclusterspecifc=FALSE){
  #simulation id, cluster id and observation within cluster id
  cluster_id <- as.factor(rep(1:m, each = n))
  obs_id <- rep(1:n,m)
  
  if(randomized){
    hours <- obs_id
    for (i in 1:m){
      hours_positive = runif(n/2, 0, 3)
      hours[which(cluster_id==i)] = append(hours_positive, -hours_positive)
    }
    x <- rep(7.5, times= n*m) + hours
  }
  else{
    x <- rnorm(n*m, mean = 7.5, sd =3)
  }
  
  #Cluster-mean of covariates
  x_mean_cluster <- c()
  for (i in seq(1, n*m, n)){
    x_mean_cluster <- append(x_mean_cluster, rep(mean(x[i:(i+n-1)]), n))
  }
  
  x_dev_from_mean <- x - x_mean_cluster
  if (nonclusterspecifc){
    u0<- rep(0, n*m)
    u1<- rep(0, n*m)
    mean_y <- b0 + bW*x_dev_from_mean + bB*x_mean_cluster
  }
  else{
    if(slopes){
      #Cluster-specific effect
      u0 <- rep(rnorm(m, 0, 5), each = n)
      u1 <- rep(rnorm(m, 0, 1), each = n)
      mean_y <- b0 + bW*x_dev_from_mean + bB*x_mean_cluster + u0 + u1*x 
    }else{
      #Cluster-specific effect
      u0 <- rep(rnorm(m, 0, 5), each = n)
      u1<- rep(0, n*m)
      mean_y <- b0 + bW*x_dev_from_mean + bB*x_mean_cluster + u0
    }
  }
  
  y <- rnorm(n*m, mean=mean_y, sd=0.5)
  
  #Cluster-mean of covariates
  y_mean_cluster <- c()
  for (i in seq(1, n*m, n)){
    y_mean_cluster <- append(y_mean_cluster, rep(mean(y[i:(i+n-1)]), n))
  }
  
  df <- data.frame(y, y_mean_cluster, cluster_id, obs_id, x, x_mean_cluster, x_dev_from_mean, u0, u1)
  return(df)
}

bw_lines <- function(df){
  cluster_intercept = unique(df$u0)
  cluster_mean_x = unique(df$x_mean_cluster)
  cluster_mean_y = unique(df$y_mean_cluster)#c(mean(df[df$cluster_id==1,]$y),mean(df[df$cluster_id==2,]$y),mean(df[df$cluster_id==3,]$y))
  cluster_id = unique(df$cluster_id)
  cluster_slopes = unique(df$u1)
  intercept = cluster_intercept+b0+bB*cluster_mean_x-bW*cluster_mean_x
  return(data.frame(intercept, cluster_mean_x, cluster_mean_y, cluster_id, cluster_slopes))
}

b0=0
bW=-2
bB=5
n=20
m=10

df = simulate_data_BW(b0, bW, bB, n, m)
df_mean = bw_lines(df)

plot <- ggplot(df, aes(x=x, y=y, color=as.factor(cluster_id)))+
  geom_point(aes(shape="Observed"),show.legend = FALSE, size=3)+
  geom_abline(data=df_mean, aes(intercept=intercept, slope=bW, color=as.factor(cluster_id), linetype="Within"), linewidth=1, show.legend = FALSE)+
  geom_point(data=df_mean, aes(x=cluster_mean_x,y=cluster_mean_y, shape="Cluster-mean", color=as.factor(cluster_id)), size=7)+
  geom_abline( aes(intercept=b0, slope=bB,linetype="Between"), linewidth=1)+
  scale_linetype_manual(name="Effects", values= c( "solid", "dashed"))+
  guides(color="none")+
  scale_shape_manual(name = "Values", values = c(17, 16))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=17),legend.text = element_text(size=16),legend.title = element_text(size=17))+
  labs(x=TeX(r'(Daily workload (hours))'),y=TeX(r'(Well-being)'))
plot


df = simulate_data_BW(b0, bW, bB, n, m, randomized = TRUE)
df_mean = bw_lines(df)


plot <- ggplot(df, aes(x=x, y=y, color=as.factor(cluster_id)))+
  geom_point(aes(shape="Observed"),show.legend = FALSE, size=3)+
  geom_abline(data=df_mean, aes(intercept=intercept, slope=bW, color=as.factor(cluster_id), linetype="Within"), linewidth=1, show.legend = FALSE)+
  geom_point(data=df_mean, aes(x=cluster_mean_x,y=cluster_mean_y, shape="Cluster-mean", color=as.factor(cluster_id)), size=7)+
  scale_linetype_manual(name="Effects", values= c( "dashed"))+
  guides(color="none")+
  scale_shape_manual(name = "Values", values = c(17, 16))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=17),legend.text = element_text(size=16),legend.title = element_text(size=17))+
  labs(x=TeX(r'(Daily workload (hours))'),y=TeX(r'(Well-being)'))
plot


