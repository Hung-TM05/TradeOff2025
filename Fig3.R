setwd('/Users/tmhung/Library/CloudStorage/OneDrive-個人/ESSLAB/1. Manuscript/1. Trade-off/0. Nature EE format/Final codes')
library(dplyr)
library(vegan)
library(lme4)
library(MuMIn)
library(nlme)
library(AICcmodavg)

#---- p-value function (RF: Family, Location) ----
p_value = function(df_x, df_y, xlim, ylim){
  
  md <- lmer(df_y~ df_x+ (1|Family)+ (1|Location), data = df)
  sum = summary(md)
  slope = round(sum[["coefficients"]][2], 4)
  r = round(r.squaredGLMM(md)[1], 2)
  p = round(car::Anova(md)[1,3], 3)
  
  text_pos_x = sum(xlim)/2
  
  if(p < 0.05 && p >= 0.01){
    text(text_pos_x, ylim, paste('slope = ',slope,', r2 = ',r,', p = ',p,'*',sep=''), cex=1)
  }else if(p < 0.01 && p >= 0.001){
    text(text_pos_x, ylim, paste('slope = ',slope,', r2 = ',r,', p = ',p,'**',sep=''), cex=1)
  }else if(p < 0.001){
    text(text_pos_x, ylim, paste('slope = ',slope,', r2 = ',r,', p < 0.001***',sep=''), cex=1)
  }else{
    text(text_pos_x, ylim, paste('slope = ',slope,', r2 = ',r,', p = ',p,sep=''), cex=1)
  }
  
  return(p)
}

#### MAIN ####
df = read.csv('./Empirical_moth_data.csv',header = TRUE)
loc_list = 'All'

## open plot file ####
pdf(file = './Fig3.pdf', width = 8, height = 8, onefile = TRUE)
par(mfrow=c(2,2))
#---- Fig.3a: CTmax - CTmin ----
plot(x=df$CTmin, y=df$CTmax, las=1, cex= 2, 
     xlab = '', ylab = '',
     xlim = c(-6,18), ylim = c(25,50), type="n", axes = F)
box()

axis(side=1, at=seq(-6,18,6),cex.axis=1.2)
mtext(side=1, line=3, "Critical thermal minimum (°C)", col="black", font=1,cex=1.3)
axis(side=2, at=seq(25,50,5),cex.axis=1.2, las=1)
mtext(side=2, line=3, "Critical thermal maximum (°C)", col="black", font=1,cex=1.3)

df_x = df$CTmin
df_y = df$CTmax
md_x <- lmer(CTmax~ CTmin+ Weight+ (1|Family) +(1|Location), data = df)
mydata.x <- data.frame(expand.grid(CTmin = seq(min(df_x),max(df_x),0.1), Weight = mean(df$Weight)))
pred.x <- predictSE(md_x,mydata.x,se.fit=TRUE)
upper.x <- pred.x$fit+1.96*pred.x$se.fit
lower.x <- pred.x$fit-1.96*pred.x$se.fit

df_ = df %>% filter(Location == 'Malaysia')
points(x=df_$CTmin,y=df_$CTmax,col=scales::alpha('#ED784A', 0.7),pch=16,cex=0.7)
df_ = df %>% filter(Location == 'Taiwan')
points(x=df_$CTmin,y=df_$CTmax,col=scales::alpha('#5DAC81', 0.7),pch=16,cex=0.7)
df_ = df %>% filter(Location == 'China')
points(x=df_$CTmin,y=df_$CTmax,col=scales::alpha('#51A8DD', 0.7),pch=16,cex=0.7)

sum = summary(md_x)
slope = round(sum[["coefficients"]][2], 2)
r = round(r.squaredGLMM(md_x)[1], 2)
p = round(car::Anova(md_x)[1,3], 3)

text_pos_x = sum(c(-6,18))/2
ylim = 25
if(p < 0.05 && p >= 0.01){
  text(text_pos_x, ylim, paste('slope = ',slope,', r2 = ',r,', p = ',p,'*',sep=''), cex=1)
}else if(p < 0.01 && p >= 0.001){
  text(text_pos_x, ylim, paste('slope = ',slope,', r2 = ',r,', p = ',p,'**',sep=''), cex=1)
}else if(p < 0.001){
  text(text_pos_x, ylim, paste('slope = ',slope,', r2 = ',r,', p < 0.001***',sep=''), cex=1)
}else{
  text(text_pos_x, ylim, paste('slope = ',slope,', r2 = ',r,', p = ',p,sep=''), cex=1)
}

if(p <= 0.05){
  polygon(x=c(mydata.x$CTmin,rev(mydata.x$CTmin)),y=c(lower.x,rev(upper.x)),
          col=scales::alpha('#787878', 0.1),border=NA)
}
if(p > 0.05) lty = 2 else lty = 1
lines(x=mydata.x$CTmin,y=pred.x$fit,col=scales::alpha('black', 0.8),lwd=1.5,lty=lty)

#---- Fig.3b: TTlimits - Body length ----
xlim_ = c(5,40,5)
ylim_ = c(-10,50,10)
for(loc in loc_list){
  if(loc == 'All') col.list = c('#787878', '#787878', '#787878', '#787878')else 
    col.list = c('orange', 'peachpuff', 'royalblue', 'lightblue')
  
  if(loc == 'Malaysia'){
    df = df %>% filter(Location == 'Malaysia')
    col_ = '#ED784A'
  }
  else if(loc == 'Taiwan'){
    df = df %>% filter(Location == 'Taiwan')
    col_ = '#5DAC81'
  }
  else if(loc == 'China'){
    df = df %>% filter(Location == 'China')
    col_ = '#51A8DD'
  }else if(loc == 'All'){
    df = df
    col_ = '#787878'
  }
  
  plot(x=df$B_length, y=df$TTrange, las=1, cex= 2, 
       xlab = '', ylab = '', main = loc,
       xlim = c(xlim_[1],xlim_[2]), ylim = c(ylim_[1],ylim_[2]), type="n", axes = F)
  box()
  
  axis(side=1, at=seq(xlim_[1],xlim_[2],xlim_[3]),cex.axis=1.4)
  mtext(side=1, line=3, "Body length (mm)", col="black", font=1,cex=1.2)
  axis(side=2, at=seq(ylim_[1],ylim_[2],ylim_[3]),cex.axis=1.4, las=1)
  mtext(side=2, line=3, "Critical thermal limits (°C)", col="black", font=1,cex=1.2)
  
  # CTmax
  df_x = df$B_length
  df_y = df$CTmax
  if(loc == 'All'){ 
    md_x <- lmer(CTmax~ B_length+ (1|Family) +(1|Location), data = df)
  }else{
    md_x <- lmer(CTmax~ B_length+ (1|Family), data = df)
  }
  mydata.x <- data.frame(expand.grid(B_length = seq(min(df_x),max(df_x),0.1)))
  pred.x <- predictSE(md_x,mydata.x,se.fit=TRUE)
  upper.x <- pred.x$fit+1.96*pred.x$se.fit
  lower.x <- pred.x$fit-1.96*pred.x$se.fit
  if(loc == 'All'){
    df_ = df %>% filter(Location == 'Malaysia')
    points(x=df_$B_length,y=df_$CTmax,col=scales::alpha('#ED784A', 0.7),pch=16,cex=0.7)
    df_ = df %>% filter(Location == 'Taiwan')
    points(x=df_$B_length,y=df_$CTmax,col=scales::alpha('#5DAC81', 0.7),pch=16,cex=0.7)
    df_ = df %>% filter(Location == 'China')
    points(x=df_$B_length,y=df_$CTmax,col=scales::alpha('#51A8DD', 0.7),pch=16,cex=0.7)
  }else{
    points(x=df$B_length,y=df$CTmax,col=scales::alpha(col.list[1], 0.7),pch=16,cex=0.7)
  }
  
  if(loc == 'All'){ 
    p_ = p_value(df_x, df_y, c(xlim_[1],xlim_[2]), ylim_[2])
  }else{
    p_ = p_value_loc(df_x, df_y, c(xlim_[1],xlim_[2]), ylim_[2])
  }
  
  if(p_ > 0.05) lty = 2 else{
    lty = 1
    polygon(x=c(mydata.x$B_length,rev(mydata.x$B_length)),y=c(lower.x,rev(upper.x)),
            col=scales::alpha(col.list[2], 0.2),border=NA)
  }
  lines(x=mydata.x$B_length,y=pred.x$fit,col=scales::alpha('black', 0.8),lwd=1.5,lty=lty)
  
  # CTmin
  df_x = df$B_length
  df_y = df$CTmin
  if(loc == 'All'){ 
    md_x <- lmer(CTmin~ B_length+ (1|Family) +(1|Location), data = df)
  }else{
    md_x <- lmer(CTmin~ B_length+ (1|Family), data = df)
  }
  mydata.x <- data.frame(expand.grid(B_length = seq(min(df_x),max(df_x),0.1)))
  pred.x <- predictSE(md_x,mydata.x,se.fit=TRUE)
  upper.x <- pred.x$fit+1.96*pred.x$se.fit
  lower.x <- pred.x$fit-1.96*pred.x$se.fit
  if(loc == 'All'){
    df_ = df %>% filter(Location == 'Malaysia')
    points(x=df_$B_length,y=df_$CTmin,col=scales::alpha('#ED784A', 0.7),pch=16,cex=0.7)
    df_ = df %>% filter(Location == 'Taiwan')
    points(x=df_$B_length,y=df_$CTmin,col=scales::alpha('#5DAC81', 0.7),pch=16,cex=0.7)
    df_ = df %>% filter(Location == 'China')
    points(x=df_$B_length,y=df_$CTmin,col=scales::alpha('#51A8DD', 0.7),pch=16,cex=0.7)
  }else{
    points(x=df$B_length,y=df$CTmin,col=scales::alpha(col.list[3], 0.7),pch=16,cex=0.7)
  }
  
  if(loc == 'All'){ 
    p_ = p_value(df_x, df_y, c(xlim_[1],xlim_[2]), ylim_[1])
  }else{
    p_ = p_value_loc(df_x, df_y, c(xlim_[1],xlim_[2]), ylim_[1])
  }
  
  if(p_ > 0.05) lty = 2 else{
    lty = 1
    polygon(x=c(mydata.x$B_length,rev(mydata.x$B_length)),y=c(lower.x,rev(upper.x)),
            col=scales::alpha(col.list[4], 0.2),border=NA)
  }
  lines(x=mydata.x$B_length,y=pred.x$fit,col=scales::alpha('black', 0.8),lwd=1.5,lty=lty)
}

#---- Fig.3c: TTrange - Body length ----
plot(x=df$B_length, y=df$TTrange, las=1, cex= 2, 
     xlab = '', ylab = '', main = loc,
     xlim = c(5,40), ylim = c(20,50), type="n", axes = F)
box()

axis(side=1, at=seq(5,40,5),cex.axis=1.4)
mtext(side=1, line=3, "Body length (mm)", col="black", font=1,cex=1.3)
axis(side=2, at=seq(20,50,5),cex.axis=1.4, las=1)
mtext(side=2, line=3, "Thermal tolerance range (°C)", col="black", font=1,cex=1.3)

df_x = df$B_length
df_y = df$TTrange
md_x <- lmer(TTrange~ B_length+ (1|Family)+ (1|Location), data = df)
mydata.x <- data.frame(expand.grid(B_length = seq(min(df_x),max(df_x),0.1)))
pred.x <- predictSE(md_x,mydata.x,se.fit=TRUE)
upper.x <- pred.x$fit+1.96*pred.x$se.fit
lower.x <- pred.x$fit-1.96*pred.x$se.fit
if(loc == 'All'){
  df_ = df %>% filter(Location == 'Malaysia')
  points(x=df_$B_length,y=df_$TTrange,col=scales::alpha('#ED784A', 0.7),pch=16,cex=0.7)
  df_ = df %>% filter(Location == 'Taiwan')
  points(x=df_$B_length,y=df_$TTrange,col=scales::alpha('#5DAC81', 0.7),pch=16,cex=0.7)
  df_ = df %>% filter(Location == 'China')
  points(x=df_$B_length,y=df_$TTrange,col=scales::alpha('#51A8DD', 0.7),pch=16,cex=0.7)
}else{
  points(x=df$B_length,y=df$TTrange,col=scales::alpha(col_, 0.7),pch=16,cex=0.7)
}
polygon(x=c(mydata.x$B_length,rev(mydata.x$B_length)),y=c(lower.x,rev(upper.x)),
        col=scales::alpha('#787878', 0.1),border=NA)
p_ = p_value(df_x, df_y, c(5,40), min(df$TTrange))
if(p_ > 0.05) lty = 2 else lty = 1
lines(x=mydata.x$B_length,y=pred.x$fit,col=scales::alpha('black', 0.8),lwd=1.5,lty=lty)

#---- Fig.3d: Body length - STmax ----
plot(x=df$STmax, y=df$B_length, las=1, cex= 2, 
     xlab = '', ylab = '', main = loc,
     xlim = c(10,30), ylim = c(5,40), type="n", axes = F)
box()

axis(side=1, at=seq(10,30,5),cex.axis=1.3)
mtext(side=1, line=3, "Averanged max temperature of warmest month (°C)", col="black", font=1,cex=1)
axis(side=2, at=seq(5,40,5),cex.axis=1.3, las=1)
mtext(side=2, line=3, "Body length (mm)", col="black", font=1,cex=1)

df_x = df$STmax
df_y = df$B_length
md_x <- lmer(B_length~ STmax+ (1|Family)+ (1|Location), data = df)
mydata.x <- data.frame(expand.grid(STmax = seq(min(df_x),max(df_x),0.1)))
pred.x <- predictSE(md_x,mydata.x,se.fit=TRUE)
upper.x <- pred.x$fit+1.96*pred.x$se.fit
lower.x <- pred.x$fit-1.96*pred.x$se.fit

df_ = df %>% filter(Location == 'Malaysia')
points(x=df_$STmax,y=df_$B_length,col=scales::alpha('#ED784A', 0.7),pch=16,cex=0.7)
df_ = df %>% filter(Location == 'Taiwan')
points(x=df_$STmax,y=df_$B_length,col=scales::alpha('#5DAC81', 0.7),pch=16,cex=0.7)
df_ = df %>% filter(Location == 'China')
points(x=df_$STmax,y=df_$B_length,col=scales::alpha('#51A8DD', 0.7),pch=16,cex=0.7)

polygon(x=c(mydata.x$STmax,rev(mydata.x$STmax)),y=c(lower.x,rev(upper.x)),
        col=scales::alpha('#787878', 0.1),border=NA)
p_ = p_value(df_x, df_y, c(10,30), 5)
if(p_ > 0.05) lty = 2 else lty = 1
lines(x=mydata.x$STmax,y=pred.x$fit,col=scales::alpha('black', 0.8),lwd=1.5,lty=lty)

## close plot file ####
dev.off()