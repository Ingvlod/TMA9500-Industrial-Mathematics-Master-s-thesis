library(lme4)
#Packages for plotting
library(ggplot2)
library(ggridges)
library(latex2exp)
library(viridis)

load("NurseGrades.RData")
df=NurseGrades
range(df$TotalScore)
range(df$Poeng_kar_vgs)
count=c()
for (i in range(df$Campus)[1]:range(df$Campus)[2]){
  count = append(count, nrow(df[df$Campus==i,]))
}
count
mean(count)
median(count)
df$Dev_from_campus_mean_karakter_vgs <- df$Poeng_kar_vgs - df$Campus_mean_poeng_kar_vgs

lmm_bw=lmer(TotalScore~  Dev_from_campus_mean_karakter_vgs + Campus_mean_poeng_kar_vgs + (1|Campus), data = df)
lmm=lmer(TotalScore~ Poeng_kar_vgs + (1|Campus), data = df)
summary(lmm)
summary(lmm_bw)

Campus_mean_poeng_kar_vgs <- unique(df$Campus_mean_poeng_kar_vgs)
df_slope <- data.frame(Intercepts=coef(lmm)$Campus[,1], Slope=coef(lmm)$Campus$Poeng_kar_vgs, Population_inter=fixef(lmm)[1])
df_slope_bw <- data.frame(Intercepts=coef(lmm_bw)$Campus[,1]+fixef(lmm_bw)[3]*Campus_mean_poeng_kar_vgs-fixef(lmm_bw)[2]*Campus_mean_poeng_kar_vgs, Slope=fixef(lmm_bw)[2], Population_inter=fixef(lmm_bw)[1]+(fixef(lmm_bw)[3]-fixef(lmm_bw)[2])*mean(Campus_mean_poeng_kar_vgs))

plot_lmm <- ggplot(df, aes(x=Poeng_kar_vgs, y=TotalScore))+geom_blank()+
  geom_abline(data=df_slope, aes(intercept=Intercepts, slope=Slope), linetype="dashed")+
  geom_abline(data=df_slope, aes(intercept=Population_inter, slope=Slope), linewidth=1.5, color="red")+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18))+
  ylim(25,90)+
  labs(x=TeX(r'(Average grade, VGS)'),y=TeX(r'(Total score)'))+
  annotate(geom = "text", x=45, y=80, label=TeX(r'($\hat{\beta}_1=1.152$)'), size=7, color="black")
plot_lmm

plot_lmm_bw <- ggplot(df, aes(x=Poeng_kar_vgs, y=TotalScore))+geom_blank()+
  geom_abline(data=df_slope_bw, aes(intercept=Intercepts, slope=Slope), linetype="dashed")+
  geom_abline(data=df_slope_bw, aes(intercept=Population_inter, slope=Slope), linewidth=1.5, color="red")+
  ylim(25,90)+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18))+
  labs(x=TeX(r'(Average grade, VGS)'),y=TeX(r'(Total score)'))+
  annotate(geom = "text", x=45, y=80, label=TeX(r'($\hat{\beta}_W=1.058$)'), size=7, color="black")
plot_lmm_bw


campus_mean <- function(df){
  Campus_mean_final_grade = df$FinalGrade
  Campus_mean_karakter_mat_vgs = df$FinalGrade
  Campus_mean_poeng_kar_vgs = df$FinalGrade
  Campus_mean_TotalScore = df$FinalGrade
  
  for (i in range(df$Campus)[1]:range(df$Campus)[2]){
    Campus_mean_final_grade[which(df$Campus==i)]= mean(df[df$Campus==i, ]$FinalGrade)
    Campus_mean_karakter_mat_vgs[which(df$Campus==i)] = mean(df[df$Campus==i, ]$Karakter_mat_vgs)
    Campus_mean_poeng_kar_vgs[which(df$Campus==i)] = mean(df[df$Campus==i, ]$Poeng_kar_vgs)
    Campus_mean_TotalScore[which(df$Campus==i)] = mean(df[df$Campus==i, ]$TotalScore)
  }
  
  df$Campus_mean_final_grade <- Campus_mean_final_grade
  df$Campus_mean_karakter_mat_vgs <- Campus_mean_karakter_mat_vgs
  df$Campus_mean_poeng_kar_vgs <- Campus_mean_poeng_kar_vgs
  df$Campus_mean_TotalScore <- Campus_mean_TotalScore
  
  df$Dev_from_campus_mean_karakter_vgs <- df$Poeng_kar_vgs - df$Campus_mean_poeng_kar_vgs
  
  return(df)
}


sample_three_stundents <- function(n=1000, df){
  beta1=c()
  betaw=c()
  for (i in 1:n){
    df_sampled <- df[df$Campus==1, ][sample(1:nrow(df[df$Campus==1, ]),3,replace=FALSE),]
    for (i in unique(df$Campus)){
      df_sampled <- rbind(df_sampled, df[df$Campus==i, ][sample(1:nrow(df[df$Campus==i, ]),3,replace=FALSE),])
    }
    df_sampled <- df_sampled[-c(1,2,3),]
    df_sampled <- campus_mean(df_sampled)
    
    lmm_bw=lmer(TotalScore~  Dev_from_campus_mean_karakter_vgs + Campus_mean_poeng_kar_vgs + (1|Campus), data = df_sampled)
    lmm=lmer(TotalScore~ Poeng_kar_vgs + (1|Campus), data = df_sampled)
    beta1 = append(beta1, fixef(lmm)[2])
    betaw = append(betaw, fixef(lmm_bw)[2])
  }
  
  return (data.frame(beta1,betaw))
}

df_beta=sample_three_stundents(1000, df)
mean(df_beta$beta1)
mean(df_beta$betaw)


#One sampling for illustration 

set.seed(6378)
#Sample 3 students from each campus
df_sampled <- df[df$Campus==1, ][sample(1:nrow(df[df$Campus==1, ]),3,replace=FALSE),]
for (i in unique(df$Campus)){
  df_sampled <- rbind(df_sampled, df[df$Campus==i, ][sample(1:nrow(df[df$Campus==i, ]),3,replace=FALSE),])
}
df_sampled<- df_sampled[-c(1,2,3),]
df_sampled <- campus_mean(df_sampled)



lmm_bw=lmer(TotalScore~  Dev_from_campus_mean_karakter_vgs + Campus_mean_poeng_kar_vgs + (1|Campus), data = df_sampled)
lmm=lmer(TotalScore~ Poeng_kar_vgs + (1|Campus), data = df_sampled)
summary(lmm_bw)
summary(lmm)
fixef(lmm)
round(fixef(lmm)[2], 3)

Campus_mean_poeng_kar_vgs <- unique(df_sampled$Campus_mean_poeng_kar_vgs)


df_slope <- data.frame(Intercepts=coef(lmm)$Campus[,1], Slope=coef(lmm)$Campus$Poeng_kar_vgs, Population_inter=fixef(lmm)[1])
df_slope_bw <- data.frame(Intercepts=coef(lmm_bw)$Campus[,1]+fixef(lmm_bw)[3]*Campus_mean_poeng_kar_vgs-fixef(lmm_bw)[2]*Campus_mean_poeng_kar_vgs, Slope=fixef(lmm_bw)[2], Population_inter=fixef(lmm_bw)[1]+(fixef(lmm_bw)[3]-fixef(lmm_bw)[2])*mean(Campus_mean_poeng_kar_vgs))

plot_lmm <- ggplot(df_sampled, aes(x=Poeng_kar_vgs, y=TotalScore))+geom_blank()+
  geom_abline(data=df_slope, aes(intercept=Intercepts, slope=Slope), linetype="dashed")+
  geom_abline(data=df_slope, aes(intercept=Population_inter, slope=Slope), linewidth=1.5, color="red")+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18))+
  ylim(25,90)+
  labs(x=TeX(r'(Average grade, VGS)'),y=TeX(r'(Total score)'))+
  annotate(geom = "text", x=45, y=80, label=TeX(r'($\hat{\beta}_1=1.97$)'), size=7, color="black")
plot_lmm

plot_lmm_bw <- ggplot(df_sampled, aes(x=Poeng_kar_vgs, y=TotalScore))+geom_blank()+
  geom_abline(data=df_slope_bw, aes(intercept=Intercepts, slope=Slope), linetype="dashed")+
  geom_abline(data=df_slope_bw, aes(intercept=Population_inter, slope=Slope), linewidth=1.5, color="red")+
  ylim(25,90)+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18))+
  labs(x=TeX(r'(Average grade, VGS)'),y=TeX(r'(Total score)'))+
  annotate(geom = "text", x=45, y=80, label=TeX(r'($\hat{\beta}_W=1.357$)'), size=7, color="black")
plot_lmm_bw



