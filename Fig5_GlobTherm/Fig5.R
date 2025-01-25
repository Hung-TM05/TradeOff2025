library(dplyr)
library(vegan)
library(lme4)
library(MuMIn)
library(nlme)
library(scales)
library(AICcmodavg)

setwd("/Users/tmhung/Library/CloudStorage/OneDrive-個人/ESSLAB/1. Manuscript/1. Trade-off/0. Nature EE format/Final codes")
All = read.csv('./Globtherm_data.csv', header = TRUE)
All = All %>% tidyr::drop_na(Chelsa_STmax)

#### MAIN ####
pdf("./Fig5.pdf", width=8, height=8, onefile = TRUE)

par(mfrow=c(2,2))
col.list = c('orange', 'peachpuff', 'royalblue', 'lightblue')

for(i in seq(1)){
  if(i == 2){
    df = All %>% filter(Group == 'Bird' | Group == 'Mammal') %>% filter(Lat_max == Lat_min & Long_max == Long_min)
    main_text = 'Endotherm'
  }else if(i == 1){
    df = All %>% filter(Group == 'Amphibian' | Group == 'Reptile' | Group == 'Invertebrate') %>% filter(Lat_max == Lat_min & Long_max == Long_min)
    main_text = 'Ectotherm'
  }
  
  #---- Fig.5a: STmax-CTlimits ----
  print('Now processing panel 1')
  plot(x=df$Chelsa_STmax, y=df$CTmax, las=1, cex= 3, 
       xlab = '', ylab = '',main = main_text,
       xlim = c(-10,50), ylim = c(-20,60),
       type="n", axes = F)
  box()
  # x-axis
  axis(side=1, at=seq(-10,50,10),cex.axis=1.2)
  mtext(side=1, line=3, "Averaged max temperature of warmest month (°C)", col="black", font=1,cex=1.2)
  # y-axis
  axis(side=2, las=2, at=seq(-20,60,10),cex.axis=1.2)
  mtext(side=2, line=3, "Critical thermal limits (°C)", col="black", font=1,cex=1.2)
  
  df_x = df$Chelsa_STmax
  df_y = df$CTmax
  md_x <- lmer(df_y~ df_x+ (1|Group), data = df)
  mydata.x <- data.frame(expand.grid(df_x = seq(min(df_x),max(df_x),0.1)))
  pred.x <- predictSE(md_x,mydata.x,se.fit=TRUE)
  upper.x <- pred.x$fit+1.96*pred.x$se.fit
  lower.x <- pred.x$fit-1.96*pred.x$se.fit
  polygon(x=c(mydata.x$df_x,rev(mydata.x$df_x)),y=c(lower.x,rev(upper.x)),
          col=scales::alpha('grey70', 0.2),border=NA)
  
  sum = summary(md_x)
  slope = round(sum[["coefficients"]][2], 4)
  r = round(r.squaredGLMM(md_x)[1], 4)
  p = round(car::Anova(md_x)[1,3], 4)
  x_pos = sum(-10,50)/2
  y_pos = max(df_y)
  
  if(p < 0.05 && p >= 0.01){
    text(x_pos, y_pos, paste('p=',p,'*, slope=',slope,', r2=',r,sep=''), cex=1)
  }else if(p < 0.01 && p >= 0.001){
    text(x_pos, y_pos, paste('p=',p,'**, slope=',slope,', r2=',r,sep=''), cex=1)
  }else if(p < 0.001){
    text(x_pos, y_pos, paste('p < 0.001***, slope=',slope,', r2=',r,sep=''), cex=1)
  }else{
    text(x_pos, y_pos, paste('p=',p,', slope=',slope,', r2=',r,sep=''), cex=1)
  }
  
  if(p > 0.05) lty = 2 else lty = 1
  lines(x=mydata.x$df_x,y=pred.x$fit,col='black',lwd=1.5,lty=lty)
  points(x=df$CTmax_Tmax,y=df$CTmax,col=scales::alpha(col.list[1], 0.7),pch=16,cex=0.7)
  
  # CTmin
  df_x = df$Chelsa_STmax
  df_y = df$CTmin
  md_x <- lmer(df_y~ df_x+ (1|Group), data = df)
  mydata.x <- data.frame(expand.grid(df_x = seq(min(df_x),max(df_x),0.1)))
  pred.x <- predictSE(md_x,mydata.x,se.fit=TRUE)
  upper.x <- pred.x$fit+1.96*pred.x$se.fit
  lower.x <- pred.x$fit-1.96*pred.x$se.fit
  polygon(x=c(mydata.x$df_x,rev(mydata.x$df_x)),y=c(lower.x,rev(upper.x)),
          col=scales::alpha('grey70', 0.2),border=NA)
  
  sum = summary(md_x)
  slope = round(sum[["coefficients"]][2], 4)
  r = round(r.squaredGLMM(md_x)[1], 4)
  p = round(car::Anova(md_x)[1,3], 4)
  x_pos = sum(-10,50)/2
  y_pos = min(df_y)
  
  if(p < 0.05 && p >= 0.01){
    text(x_pos, y_pos, paste('p=',p,'*, slope=',slope,', r2=',r,sep=''), cex=1)
  }else if(p < 0.01 && p >= 0.001){
    text(x_pos, y_pos, paste('p=',p,'**, slope=',slope,', r2=',r,sep=''), cex=1)
  }else if(p < 0.001){
    text(x_pos, y_pos, paste('p < 0.001***, slope=',slope,', r2=',r,sep=''), cex=1)
  }else{
    text(x_pos, y_pos, paste('p=',p,', slope=',slope,', r2=',r,sep=''), cex=1)
  }
  
  if(p > 0.05) lty = 2 else lty = 1
  lines(x=mydata.x$df_x,y=pred.x$fit,col='black',lwd=1.5,lty=lty)
  points(x=df$Chelsa_STmax,y=df$CTmin,col=scales::alpha(col.list[3], 0.7),pch=16,cex=0.7)
  
  #---- Fig.5b: STmin-CTlimits ----
  print('Now processing panel 2')
  plot(x=df$Chelsa_STmin, y=df$CTmax, las=1, cex= 3, 
       xlab = '', ylab = '',main = main_text,
       xlim = c(-20,30), ylim = c(-20,60),
       type="n", axes = F)
  box()
  # x-axis
  axis(side=1, at=seq(-20,30,10),cex.axis=1.2)
  mtext(side=1, line=3, "Minimum ambient temperature (°C)", col="black", font=1,cex=1.2)
  # y-axis
  axis(side=2, las=2, at=seq(-20,60,10),cex.axis=1.2)
  mtext(side=2, line=3, "Critical thermal limits (°C)", col="black", font=1,cex=1.2)
  
  df_x = df$Chelsa_STmin
  df_y = df$CTmax
  md_x <- lmer(df_y~ df_x+ (1|Group), data = df)
  mydata.x <- data.frame(expand.grid(df_x = seq(min(df_x),max(df_x),0.1)))
  pred.x <- predictSE(md_x,mydata.x,se.fit=TRUE)
  upper.x <- pred.x$fit+1.96*pred.x$se.fit
  lower.x <- pred.x$fit-1.96*pred.x$se.fit
  polygon(x=c(mydata.x$df_x,rev(mydata.x$df_x)),y=c(lower.x,rev(upper.x)),
          col=scales::alpha('grey70', 0.2),border=NA)
  
  sum = summary(md_x)
  slope = round(sum[["coefficients"]][2], 4)
  r = round(r.squaredGLMM(md_x)[1], 4)
  p = round(car::Anova(md_x)[1,3], 4)
  x_pos = sum(-20,30)/2
  y_pos = max(df_y)
  
  if(p < 0.05 && p >= 0.01){
    text(x_pos, y_pos, paste('p=',p,'*, slope=',slope,', r2=',r,sep=''), cex=1)
  }else if(p < 0.01 && p >= 0.001){
    text(x_pos, y_pos, paste('p=',p,'**, slope=',slope,', r2=',r,sep=''), cex=1)
  }else if(p < 0.001){
    text(x_pos, y_pos, paste('p < 0.001***, slope=',slope,', r2=',r,sep=''), cex=1)
  }else{
    text(x_pos, y_pos, paste('p=',p,', slope=',slope,', r2=',r,sep=''), cex=1)
  }
  
  if(p > 0.05) lty = 2 else lty = 1
  lines(x=mydata.x$df_x,y=pred.x$fit,col='black',lwd=1.5,lty=lty)
  points(x=df$Chelsa_STmin,y=df$CTmax,col=scales::alpha(col.list[1], 0.7),pch=16,cex=0.7)
  
  # CTmin
  df_x = df$Chelsa_STmin
  df_y = df$CTmin
  md_x <- lmer(df_y~ df_x+ (1|Group), data = df)
  mydata.x <- data.frame(expand.grid(df_x = seq(min(df_x),max(df_x),0.1)))
  pred.x <- predictSE(md_x,mydata.x,se.fit=TRUE)
  upper.x <- pred.x$fit+1.96*pred.x$se.fit
  lower.x <- pred.x$fit-1.96*pred.x$se.fit
  polygon(x=c(mydata.x$df_x,rev(mydata.x$df_x)),y=c(lower.x,rev(upper.x)),
          col=scales::alpha('grey70', 0.2),border=NA)
  
  sum = summary(md_x)
  slope = round(sum[["coefficients"]][2], 4)
  r = round(r.squaredGLMM(md_x)[1], 4)
  p = round(car::Anova(md_x)[1,3], 4)
  x_pos = sum(-30,30)/2
  y_pos = min(df_y)
  
  if(p < 0.05 && p >= 0.01){
    text(x_pos, y_pos, paste('p=',p,'*, slope=',slope,', r2=',r,sep=''), cex=1)
  }else if(p < 0.01 && p >= 0.001){
    text(x_pos, y_pos, paste('p=',p,'**, slope=',slope,', r2=',r,sep=''), cex=1)
  }else if(p < 0.001){
    text(x_pos, y_pos, paste('p < 0.001***, slope=',slope,', r2=',r,sep=''), cex=1)
  }else{
    text(x_pos, y_pos, paste('p=',p,', slope=',slope,', r2=',r,sep=''), cex=1)
  }
  
  if(p > 0.05) lty = 2 else lty = 1
  lines(x=mydata.x$df_x,y=pred.x$fit,col='black',lwd=1.5,lty=lty)
  points(x=df$Chelsa_STmin,y=df$CTmin,col=scales::alpha(col.list[3], 0.7),pch=16,cex=0.7)
  
  #---- Fig.5c: TTrange-STR ----
  print('Now processing panel 3')
  plot(x=df$Chelsa_STR, y=df$CTrange, las=1, cex= 3, 
       xlab = '', ylab = '',main = main_text,
       xlim = c(0,50), ylim = c(20,60),
       type="n", axes = F)
  box()
  # x-axis
  axis(side=1, at=seq(0,50,10),cex.axis=1.2)
  mtext(side=1, line=3, "Seasonal temperature range (°C)", col="black", font=1,cex=1.2)
  # y-axis
  axis(side=2, las=2, at=seq(20,60,10),cex.axis=1.2)
  mtext(side=2, line=3, "Thermal tolerance range (°C)", col="black", font=1,cex=1.2)
  
  df_x = df$Chelsa_STR
  df_y = df$CTrange
  md_x <- lmer(df_y~ df_x+ (1|Group), data = df)
  mydata.x <- data.frame(expand.grid(df_x = seq(min(df_x),max(df_x),0.1)))
  pred.x <- predictSE(md_x,mydata.x,se.fit=TRUE)
  upper.x <- pred.x$fit+1.96*pred.x$se.fit
  lower.x <- pred.x$fit-1.96*pred.x$se.fit
  polygon(x=c(mydata.x$df_x,rev(mydata.x$df_x)),y=c(lower.x,rev(upper.x)),
          col=scales::alpha('grey70', 0.2),border=NA)
  
  sum = summary(md_x)
  slope = round(sum[["coefficients"]][2], 4)
  r = round(r.squaredGLMM(md_x)[1], 4)
  p = round(car::Anova(md_x)[1,3], 4)
  x_pos = sum(0,50)/2
  y_pos = max(df_y)
  
  if(p < 0.05 && p >= 0.01){
    text(x_pos, y_pos, paste('p=',p,'*, slope=',slope,', r2=',r,sep=''), cex=1)
  }else if(p < 0.01 && p >= 0.001){
    text(x_pos, y_pos, paste('p=',p,'**, slope=',slope,', r2=',r,sep=''), cex=1)
  }else if(p < 0.001){
    text(x_pos, y_pos, paste('p < 0.001***, slope=',slope,', r2=',r,sep=''), cex=1)
  }else{
    text(x_pos, y_pos, paste('p=',p,', slope=',slope,', r2=',r,sep=''), cex=1)
  }
  
  if(p > 0.05) lty = 2 else lty = 1
  lines(x=mydata.x$df_x,y=pred.x$fit,col='black',lwd=1.5,lty=lty)
  points(x=df$Chelsa_STR,y=df$CTrange,col=scales::alpha('#787878', 0.7),pch=16,cex=0.7)
  
  #---- Fig.5d: CTmax-CTmin ----
  print('Now processing panel 4')
  plot(x=df$CTmin, y=df$CTmax, las=1, cex= 3, 
       xlab = '', ylab = '',main = main_text,
       xlim = c(-20,20), ylim = c(20,60), type="n", axes = F)
  box()
  # x-axis
  axis(side=1, at=seq(-20,20,10),cex.axis=1.2)
  mtext(side=1, line=3, "Critical thermal minimum (°C)", col="black", font=1,cex=1.2)
  # y-axis
  axis(side=2, las=2, at=seq(20,60,10),cex.axis=1.2)
  mtext(side=2, line=3, "Critical thermal maximum (°C)", col="black", font=1,cex=1.2)
  
  df_x = df$CTmin
  df_y = df$CTmax
  md_x <- lmer(df_y~ df_x+ (1|Group), data = df)
  mydata.x <- data.frame(expand.grid(df_x = seq(min(df_x),max(df_x),0.1)))
  pred.x <- predictSE(md_x,mydata.x,se.fit=TRUE)
  upper.x <- pred.x$fit+1.96*pred.x$se.fit
  lower.x <- pred.x$fit-1.96*pred.x$se.fit
  polygon(x=c(mydata.x$df_x,rev(mydata.x$df_x)),y=c(lower.x,rev(upper.x)),
          col=scales::alpha('gray70', 0.2),border=NA)
  
  sum = summary(md_x)
  slope = round(sum[["coefficients"]][2], 4)
  r = round(r.squaredGLMM(md_x)[1], 4)
  p = round(car::Anova(md_x)[1,3], 4)
  x_pos = sum(-20,20)/2
  y_pos = max(df_y)
  
  if(p < 0.05 && p >= 0.01){
    text(x_pos, y_pos, paste('p=',p,'*, slope=',slope,', r2=',r,sep=''), cex=1)
  }else if(p < 0.01 && p >= 0.001){
    text(x_pos, y_pos, paste('p=',p,'**, slope=',slope,', r2=',r,sep=''), cex=1)
  }else if(p < 0.001){
    text(x_pos, y_pos, paste('p < 0.001***, slope=',slope,', r2=',r,sep=''), cex=1)
  }else{
    text(x_pos, y_pos, paste('p=',p,', slope=',slope,', r2=',r,sep=''), cex=1)
  }
  
  if(p > 0.05) lty = 2 else lty = 1
  lines(x=mydata.x$df_x,y=pred.x$fit,col='black',lwd=1.5,lty=lty)
  points(x=df$CTmin,y=df$CTmax,col=scales::alpha('#787878', 0.7),pch=16,cex=0.7)
  
}

dev.off()