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
  slope = round(sum[["coefficients"]][2], 2)
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
df.all = read.csv('./Empirical_moth_data.csv',header = TRUE)
df = df.all %>% tidyr::drop_na(TTrange, W_length, Weight, B_length)
df$Weight_log = log10(df$Weight)

#---- Fig.2 ----
pdf(file = './Fig2.pdf', width = 12, height = 4, onefile = TRUE)

par(mfrow=c(1,3))

#---- Fig.2a: TTlimits - STmax ----
plot(x=df$STmax, y=df$CTmax, las=1, cex= 2, 
     xlab = '', ylab = '',
     xlim = c(14,30), ylim = c(-10,50), type="n", axes = F)
box()

axis(side=1, at=seq(14,30,4),cex.axis=1.4)
mtext(side=1, line=3, "Average max temperature of warmest month (°C)", col="black", font=1,cex=1)
axis(side=2, at=seq(-10,50,10),cex.axis=1.4, las=1)
mtext(side=2, line=3, "Critical thermal limits (°C)", col="black", font=1,cex=1)

# CTmax
df_x = df$STmax
df_y = df$CTmax
md_x <- lmer(CTmax~ STmax+ (1|Family)+ (1|Location), data = df)
mydata.x <- data.frame(expand.grid(STmax = seq(min(df_x),max(df_x),0.1)))
pred.x <- predictSE(md_x,mydata.x,se.fit=TRUE)
upper.x <- pred.x$fit+1.96*pred.x$se.fit
lower.x <- pred.x$fit-1.96*pred.x$se.fit

df_ = df %>% filter(Location == 'Malaysia')
points(x=df_$STmax,y=df_$CTmax,col=scales::alpha('#ED784A', 0.7),pch=16,cex=0.7)
df_ = df %>% filter(Location == 'Taiwan')
points(x=df_$STmax,y=df_$CTmax,col=scales::alpha('#5DAC81', 0.7),pch=16,cex=0.7)
df_ = df %>% filter(Location == 'China')
points(x=df_$STmax,y=df_$CTmax,col=scales::alpha('#51A8DD', 0.7),pch=16,cex=0.7)

polygon(x=c(mydata.x$STmax,rev(mydata.x$STmax)),y=c(lower.x,rev(upper.x)),
        col=scales::alpha('#787878', 0.1),border=NA)
p_ = p_value(df_x, df_y, c(15,35), min(df$CTmax))
if(p_ > 0.05) lty = 2 else lty = 1
lines(x=mydata.x$STmax,y=pred.x$fit,col=scales::alpha('black', 0.8),lwd=1.5,lty=lty)

# CTmin
df_x = df$STmax
df_y = df$CTmin
md_x <- lmer(CTmin~ STmax+ (1|Family)+ (1|Location), data = df)
mydata.x <- data.frame(expand.grid(STmax = seq(min(df_x),max(df_x),0.1)))
pred.x <- predictSE(md_x,mydata.x,se.fit=TRUE)
upper.x <- pred.x$fit+1.96*pred.x$se.fit
lower.x <- pred.x$fit-1.96*pred.x$se.fit

df_ = df %>% filter(Location == 'Malaysia')
points(x=df_$STmax,y=df_$CTmin,col=scales::alpha('#ED784A', 0.7),pch=16,cex=0.7)
df_ = df %>% filter(Location == 'Taiwan')
points(x=df_$STmax,y=df_$CTmin,col=scales::alpha('#5DAC81', 0.7),pch=16,cex=0.7)
df_ = df %>% filter(Location == 'China')
points(x=df_$STmax,y=df_$CTmin,col=scales::alpha('#51A8DD', 0.7),pch=16,cex=0.7)

polygon(x=c(mydata.x$STmax,rev(mydata.x$STmax)),y=c(lower.x,rev(upper.x)),
        col=scales::alpha('#787878', 0.1),border=NA)
p_ = p_value(df_x, df_y, c(15,35), min(df$CTmin))
if(p_ > 0.05) lty = 2 else lty = 1
lines(x=mydata.x$STmax,y=pred.x$fit,col=scales::alpha('black', 0.8),lwd=1.5,lty=lty)

#---- Fig.2b: TTlimits - STmin ----
plot(x=df$STmin, y=df$TTrange, las=1, cex= 2, 
     xlab = '', ylab = '',
     xlim = c(-20,20), ylim = c(-10,50), type="n", axes = F)
box()

axis(side=1, at=seq(-20,20,10),cex.axis=1.4)
mtext(side=1, line=3, "Averaged min temperature of coldest month (°C)", col="black", font=1,cex=1)
axis(side=2, at=seq(-10,50,10),cex.axis=1.4, las=1)
mtext(side=2, line=3, "Critical thermal limits (°C)", col="black", font=1,cex=1)

# CTmax
df_x = df$STmin
df_y = df$CTmax
md_x <- lmer(CTmax~ STmin+ (1|Family)+ (1|Location), data = df)
mydata.x <- data.frame(expand.grid(STmin = seq(min(df_x),max(df_x),0.1)))
pred.x <- predictSE(md_x,mydata.x,se.fit=TRUE)
upper.x <- pred.x$fit+1.96*pred.x$se.fit
lower.x <- pred.x$fit-1.96*pred.x$se.fit

df_ = df %>% filter(Location == 'Malaysia')
points(x=df_$STmin,y=df_$CTmax,col=scales::alpha('#ED784A', 0.7),pch=16,cex=0.7)
df_ = df %>% filter(Location == 'Taiwan')
points(x=df_$STmin,y=df_$CTmax,col=scales::alpha('#5DAC81', 0.7),pch=16,cex=0.7)
df_ = df %>% filter(Location == 'China')
points(x=df_$STmin,y=df_$CTmax,col=scales::alpha('#51A8DD', 0.7),pch=16,cex=0.7)

polygon(x=c(mydata.x$STmin,rev(mydata.x$STmin)),y=c(lower.x,rev(upper.x)),
        col=scales::alpha('#787878', 0.1),border=NA)
p_ = p_value(df_x, df_y, c(-20,20), min(df$CTmax))
if(p_ > 0.05) lty = 2 else lty = 1
lines(x=mydata.x$STmin,y=pred.x$fit,col=scales::alpha('black', 0.8),lwd=1.5,lty=lty)

# CTmin
df_x = df$STmin
df_y = df$CTmin
md_x <- lmer(CTmin~ STmin+ (1|Family)+ (1|Location), data = df)
mydata.x <- data.frame(expand.grid(STmin = seq(min(df_x),max(df_x),0.1)))
pred.x <- predictSE(md_x,mydata.x,se.fit=TRUE)
upper.x <- pred.x$fit+1.96*pred.x$se.fit
lower.x <- pred.x$fit-1.96*pred.x$se.fit

df_ = df %>% filter(Location == 'Malaysia')
points(x=df_$STmin,y=df_$CTmin,col=scales::alpha('#ED784A', 0.7),pch=16,cex=0.7)
df_ = df %>% filter(Location == 'Taiwan')
points(x=df_$STmin,y=df_$CTmin,col=scales::alpha('#5DAC81', 0.7),pch=16,cex=0.7)
df_ = df %>% filter(Location == 'China')
points(x=df_$STmin,y=df_$CTmin,col=scales::alpha('#51A8DD', 0.7),pch=16,cex=0.7)

polygon(x=c(mydata.x$STmin,rev(mydata.x$STmin)),y=c(lower.x,rev(upper.x)),
        col=scales::alpha('#787878', 0.1),border=NA)
p_ = p_value(df_x, df_y, c(-20,20), min(df$CTmin))
if(p_ > 0.05) lty = 2 else lty = 1
lines(x=mydata.x$STmin,y=pred.x$fit,col=scales::alpha('black', 0.8),lwd=1.5,lty=lty)

#---- Fig.2c: TTrange - STR ----
loc_list = 'All'
for(loc in loc_list){
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
  plot(x=df$STR, y=df$TTrange, las=1, cex= 2, 
       xlab = '', ylab = '', main = loc,
       xlim = c(5,35), ylim = c(20,50), type="n", axes = F)
  box()
  
  axis(side=1, at=seq(5,35,5),cex.axis=1.2)
  mtext(side=1, line=3, "Seasonal temperature range (°C)", col="black", font=1,cex=1)
  axis(side=2, at=seq(20,50,5),cex.axis=1.2, las=1)
  mtext(side=2, line=3, "Thermal tolerance range (°C)", col="black", font=1,cex=1)
  
  df_x = df$STR
  df_y = df$TTrange
  if(loc == 'All'){ 
    md_x <- lmer(TTrange~ STR+ (1|Family) +(1|Location), data = df)
  }else{
    md_x <- lmer(TTrange~ STR+ (1|Family), data = df)
  }
  mydata.x <- data.frame(expand.grid(STR = seq(min(df_x),max(df_x),0.1)))
  pred.x <- predictSE(md_x,mydata.x,se.fit=TRUE)
  upper.x <- pred.x$fit+1.96*pred.x$se.fit
  lower.x <- pred.x$fit-1.96*pred.x$se.fit
  
  if(loc == 'All'){
    df_ = df %>% filter(Location == 'Malaysia')
    points(x=df_$STR,y=df_$TTrange,col=scales::alpha('#ED784A', 0.7),pch=16,cex=0.7)
    df_ = df %>% filter(Location == 'Taiwan')
    points(x=df_$STR,y=df_$TTrange,col=scales::alpha('#5DAC81', 0.7),pch=16,cex=0.7)
    df_ = df %>% filter(Location == 'China')
    points(x=df_$STR,y=df_$TTrange,col=scales::alpha('#51A8DD', 0.7),pch=16,cex=0.7)
  }else{
    points(x=df$STR,y=df$TTrange,col=scales::alpha('#787878', 0.7),pch=16,cex=0.7)
  }
  
  if(loc == 'All'){ 
    p_ = p_value(df_x, df_y, c(5,40), min(df$TTrange))
  }else{
    p_ = p_value_loc(df_x, df_y, c(5,40), min(df$TTrange))
  }
  
  if(p_ > 0.05) lty = 2 else{
    lty = 1
    polygon(x=c(mydata.x$STR,rev(mydata.x$STR)),y=c(lower.x,rev(upper.x)),
            col=scales::alpha('#787878', 0.1),border=NA)
  }
  lines(x=mydata.x$STR,y=pred.x$fit,col=scales::alpha('black', 0.8),lwd=1.5,lty=lty)
}

dev.off()
