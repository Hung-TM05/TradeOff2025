library(dplyr)
library(vegan)
library(piecewiseSEM)
library(nlme)
library(lavaan)
library(lme4)

df <- read.csv('Empirical_moth_data.csv',header = TRUE)

## SEM model
m_CTmax = lmer(CTmax ~ STmax+ STmin+ (1|Location)+ (1|Family), data= df)
m_CTmin = lmer(CTmin ~ STmax+ STmin+ B_length+ (1|Location)+ (1|Family), data= df)
m_TTrange = lmer(TTrange ~ CTmin+ CTmax+ (1|Location)+ (1|Family), data= df)
m1 <- psem(m_CTmax, m_CTmin, m_TTrange,
           CTmax %~~% CTmin,
           data = df)
summary(m1)
# plot(m1)
