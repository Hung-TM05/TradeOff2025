setwd("/Users/tmhung/Library/CloudStorage/OneDrive-個人/ESSLAB/2. OpenResource/2. Dataset_Buckley_2022_Thermal_stress/")

library(dplyr)
library(vegan)
library(lme4)
library(MuMIn)
library(scales)
library(AICcmodavg)
# library(ggplot2)

#### functions ####
lm_p_value = function(df_x, df_y, y_pos){
  md_x = lm(df_y~ df_x, data = df)
  sum = summary(md_x)
  slope = round(sum[["coefficients"]][2], 4)
  r = round(r.squaredGLMM(md_x)[1], 4)
  p = round(car::Anova(md_x)[1,4], 4)
  x_pos = sum(x_lim[1],x_lim[2])/2
  y_pos = y_pos
  
  if(p < 0.05 && p >= 0.01){
    text(x_pos, y_pos, paste('p=',p,'*, slope=',slope,', r2=',r,sep=''), cex=1.5)
  }else if(p < 0.01 && p >= 0.001){
    text(x_pos, y_pos, paste('p=',p,'**, slope=',slope,', r2=',r,sep=''), cex=1.5)
  }else if(p < 0.001){
    text(x_pos, y_pos, paste('p < 0.001***, slope=',slope,', r2=',r,sep=''), cex=1.5)
  }else{
    text(x_pos, y_pos, paste('p=',p,', slope=',slope,', r2=',r,sep=''), cex=1.5)
  }
  return(p)
}
lm_poly = function(df_x, df_y, p){
  md_x <- lm(df_y~ df_x, data = df)
  if(p < 0.05){
    mydata.x <- data.frame(expand.grid(df_x = seq(min(df_x)-0.5,max(df_x)+0.5,0.1)))
    conf_interval <- as.data.frame(predict(md_x, newdata=data.frame(x=mydata.x), interval="confidence",
                                           level = 0.95))
    polygon(x=c(mydata.x$df_x,rev(mydata.x$df_x)),y=c(conf_interval$lwr,rev(conf_interval$upr)),
            col=scales::alpha('#787878', 0.2),border=NA)
  }
  return(md_x)
}
lmer_p_value = function(df_x, df_y, y_pos){
  md_x = lmer(df_y~ df_x+ (1|taxa), data = df)
  sum = summary(md_x)
  slope = round(sum[["coefficients"]][2], 4)
  r = round(r.squaredGLMM(md_x)[1], 4)
  p = round(car::Anova(md_x)[1,3], 4)
  x_pos = sum(x_lim[1],x_lim[2])/2
  y_pos = y_pos
  
  if(p < 0.05 && p >= 0.01){
    text(x_pos, y_pos, paste('p=',p,'*, slope=',slope,', r2=',r,sep=''), cex=1.5)
  }else if(p < 0.01 && p >= 0.001){
    text(x_pos, y_pos, paste('p=',p,'**, slope=',slope,', r2=',r,sep=''), cex=1.5)
  }else if(p < 0.001){
    text(x_pos, y_pos, paste('p < 0.001***, slope=',slope,', r2=',r,sep=''), cex=1.5)
  }else{
    text(x_pos, y_pos, paste('p=',p,', slope=',slope,', r2=',r,sep=''), cex=1.5)
  }
  return(p)
}

#### main ####
df_all <- read.csv('./tpcs.csv', head=T)
df_all$TTrange = df_all$CTmax - df_all$CTmin
df_all$Opt_margin = df_all$CTmax - df_all$Topt
df_all$Opt_ctmin = df_all$Topt - df_all$CTmin
df_all2 = df_all %>% filter(taxa != 'photosynthesis') %>% filter(is.na(CTmax)==F & is.na(CTmin)==F)
# df_all2[(df_all2$genus==''),]$genus = 'NA'

pdf("/Users/tmhung/Library/CloudStorage/OneDrive-個人/ESSLAB/1. Manuscript/1. Trade-off/0. Nature EE format/Figs_TIFF/FigS1_Buckey.pdf", width=15, height=5, onefile = TRUE)

par(mfrow=c(1,3))
# col.list = c('orange', 'peachpuff', 'lightblue', 'lightblue')
axis_size = 1.4
text_size = 1

#---- base plot ----
#---- 0. Final Figure ----
### 2. TTrange ####
taxa_ = 'ectotherm'
df = df_all2 %>% filter(taxa != 'plankton' & taxa != 'fish')
    
### a. x-axis: CTmax ####
x_lim = c(20,60,10)
y_lim = c(10,50,10)

plot(x=df$CTmax, y=df$TTrange, las=1, cex= 3, 
     xlab = '', ylab = '', main=taxa_,
     xlim = c(x_lim[1],x_lim[2]), ylim = c(y_lim[1],y_lim[2]),
     type="n", axes = F)
box()
# x-axis
axis(side=1, at=seq(x_lim[1],x_lim[2],x_lim[3]),cex.axis=axis_size)
mtext(side=1, line=3, "Critical thermal maximum (°C)", col="black", font=1,cex=text_size)
# y-axis
axis(side=2, las=2, at=seq(y_lim[1],y_lim[2],y_lim[3]),cex.axis=axis_size)
mtext(side=2, line=3, "Thermal tolerance range (°C)", col="black", font=1,cex=text_size)  

df_x = df$CTmax
df_y = df$TTrange
p = lmer_p_value(df_x, df_y, y_lim[1])

md_x <- lmer(df_y~ df_x+ (1|taxa), data = df)
mydata.x <- data.frame(expand.grid(df_x = seq(min(df_x),max(df_x),0.1)))
pred.x <- predictSE(md_x,mydata.x,se.fit=TRUE)
upper.x <- pred.x$fit+1.96*pred.x$se.fit
lower.x <- pred.x$fit-1.96*pred.x$se.fit
if(p < 0.05){
  polygon(x=c(mydata.x$df_x,rev(mydata.x$df_x)),y=c(lower.x,rev(upper.x)),
          col=scales::alpha('#787878', 0.1),border=NA)
}
if(p > 0.05) lty = 2 else lty = 1
lines(x=mydata.x$df_x,y=pred.x$fit,col='black',lwd=1.5,lty=lty)
points(x=df$CTmax,y=df$TTrange,col=scales::alpha('#787878', 0.7),pch=16,cex=1)

### b. x-axis: CTmin ####
x_lim = c(-10,20,10)
y_lim = c(10,50,10)

plot(x=df$CTmin, y=df$TTrange, las=1, cex= 3, 
     xlab = '', ylab = '', main=taxa_,
     xlim = c(x_lim[1],x_lim[2]), ylim = c(y_lim[1],y_lim[2]),
     type="n", axes = F)
box()
# x-axis
axis(side=1, at=seq(x_lim[1],x_lim[2],x_lim[3]),cex.axis=axis_size)
mtext(side=1, line=3, "Critical thermal minimum (°C)", col="black", font=1,cex=text_size)
# y-axis
axis(side=2, las=2, at=seq(y_lim[1],y_lim[2],y_lim[3]),cex.axis=axis_size)
mtext(side=2, line=3, "Thermal tolerance range (°C)", col="black", font=1,cex=text_size)  

df_x = df$CTmin
df_y = df$TTrange
p = lmer_p_value(df_x, df_y, y_lim[1])

md_x <- lmer(df_y~ df_x+ (1|taxa), data = df)
mydata.x <- data.frame(expand.grid(df_x = seq(min(df_x),max(df_x),0.1)))
pred.x <- predictSE(md_x,mydata.x,se.fit=TRUE)
upper.x <- pred.x$fit+1.96*pred.x$se.fit
lower.x <- pred.x$fit-1.96*pred.x$se.fit
if(p < 0.05){
  polygon(x=c(mydata.x$df_x,rev(mydata.x$df_x)),y=c(lower.x,rev(upper.x)),
          col=scales::alpha('#787878', 0.1),border=NA)
}

if(p > 0.05) lty = 2 else lty = 1
lines(x=mydata.x$df_x,y=pred.x$fit,col='black',lwd=1.5,lty=lty)
points(x=df$CTmin,y=df$TTrange,col=scales::alpha('#787878', 0.7),pch=16,cex=1)

### c. x-axis: Topt ####
x_lim = c(15,40,5)
y_lim = c(10,50,10)

plot(x=df$Topt, y=df$TTrange, las=1, cex= 3, 
     xlab = '', ylab = '', main=taxa_,
     xlim = c(x_lim[1],x_lim[2]), ylim = c(y_lim[1],y_lim[2]),
     type="n", axes = F)
box()
# x-axis
axis(side=1, at=seq(x_lim[1],x_lim[2],x_lim[3]),cex.axis=axis_size)
mtext(side=1, line=3, "Topt (°C)", col="black", font=1,cex=text_size)
# y-axis
axis(side=2, las=2, at=seq(y_lim[1],y_lim[2],y_lim[3]),cex.axis=axis_size)
mtext(side=2, line=3, "Thermal tolerance range (°C)", col="black", font=1,cex=text_size)  

df_x = df$Topt
df_y = df$TTrange
p = lmer_p_value(df_x, df_y, y_lim[1])

md_x <- lmer(df_y~ df_x+ (1|taxa), data = df)
mydata.x <- data.frame(expand.grid(df_x = seq(min(df_x),max(df_x),0.1)))
pred.x <- predictSE(md_x,mydata.x,se.fit=TRUE)
upper.x <- pred.x$fit+1.96*pred.x$se.fit
lower.x <- pred.x$fit-1.96*pred.x$se.fit
if(p < 0.05){
  polygon(x=c(mydata.x$df_x,rev(mydata.x$df_x)),y=c(lower.x,rev(upper.x)),
          col=scales::alpha('#787878', 0.1),border=NA)
}

if(p > 0.05) lty = 2 else lty = 1
lines(x=mydata.x$df_x,y=pred.x$fit,col='black',lwd=1.5,lty=lty)
points(x=df$Topt,y=df$TTrange,col=scales::alpha('#787878', 0.7),pch=16,cex=1)

dev.off()