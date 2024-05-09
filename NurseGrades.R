#Packages for plotting
library(ggplot2)

load("NurseGrades.RData")
df=NurseGrades
ggplot(df, aes(x=Campus_mean_karakter_mat_vgs, y=Campus_mean_final_grade, color=as.factor(Campus)))+geom_point()
ggplot(df, aes(x=Campus_mean_poeng_kar_vgs, y=Campus_mean_final_grade, color=as.factor(Campus)))+geom_point()
ggplot(df, aes(x=Karakter_mat_vgs, y=FinalGrade, color=as.factor(Campus)))+geom_point()
