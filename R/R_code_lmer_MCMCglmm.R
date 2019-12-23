####################################################
# A. Preparation
####################################################
# Note that data generation appears below the analysis section.
# You can use the simulated data table from the supplementary files to reproduce exactly the same results as presented in the paper.
# Set the work directy that is used for rading/saving data tables
# setwd("/Users/R2")
# load R required packages
# If this is done for the first time, it might need to first download and install the package
# install.packages("arm")
library(arm)
# install.packages("lme4")
# the verson 1.0-5
library(lme4)
# install.packages("MCMCglmm")
library(MCMCglmm)
# adding this for vectorisation
library(tidyverse)
library(purrr)
# using this so everybody can read and write files easily
library(here)

####################################################
# B. Analysis
####################################################
# 1. Analysis of body size (Gaussian mixed models)
#---------------------------------------------------
# Clear memory
rm(list = ls())
# Read body length data (Gaussian, available for both sexes)
Data <- read.csv(here("data","BeetlesBody.csv"))
# Fit null model without fixed effects (but including all random effects)
m0 <- lmer(BodyL ~ 1 + (1 | Population) + (1 | Container), data = Data)
# MCMCglmm model
# prior 
prior <- list(R = list(V = 1, nu = 0.002),
               G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                        G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

# model
mm0<-MCMCglmm(BodyL ~ 1 , random = ~ Population + Container, prior = prior, data=Data)

# Fit alternative model including fixed and all random effects
mF <- lmer(BodyL ~ Sex + Treatment + Habitat + (1 | Population) + (1 | Container), data = Data)
# MCMCglmm model
mmF<-MCMCglmm(BodyL ~ Sex + Treatment + Habitat , random = ~ Population + Container, prior = prior, data=Data)

# View model fits for both models
summary(m0)
summary(mm0) # MCMCglmm
summary(mF)
summary(mmF) # MCMCglmm

# Extraction of fitted value for the alternative model
# fixef() extracts coefficents for fixed effects
# mF@pp$X returns fixed effect design matrix
Fixed <- fixef(mF)[2] * mF@pp$X[, 2] + fixef(mF)[3] * mF@pp$X[, 3] + fixef(mF)[4] * mF@pp$X[, 4]
# MCMCglmm (it is probably better to get a posterior distribuiton of R2 rather than getting each varaince component - we do this below as an alternative)
mFixed <- mean(mmF$Sol[,2]) * mmF$X[, 2] + mean(mmF$Sol[, 3]) * mmF$X[, 3] + mean(mmF$Sol[ ,4]) * mmF$X[, 4]
# Calculation of the variance in fitted values
VarF <- var(Fixed)
mVarF<- var(mFixed)
# An alternative way for getting the same result
VarF <- var(as.vector(fixef(mF) %*% t(mF@pp$X)))
mVarF <- var(as.vector(apply(mmF$Sol,2,mean) %*% t(mmF$X)))
# R2GLMM(m) - marginal R2GLMM
# Equ. 26, 29 and 30
# VarCorr() extracts variance components
# attr(VarCorr(lmer.model),'sc')^2 extracts the residual variance
VarF/(VarF + VarCorr(mF)$Container[1] + VarCorr(mF)$Population[1] + attr(VarCorr(mF), "sc")^2)

# MCMCglmm - marginal
mVarF/(mVarF+sum(apply(mmF$VCV,2,mean)))
# alternative with crebile intervals
vmVarF<-numeric(1000)
for(i in 1:1000){
  Var<-var(as.vector(mmF$Sol[i,] %*% t(mmF$X)))
  vmVarF[i]<-Var}
# quicker way to do this
mmF_list <- split(mmF$Sol, seq(nrow(mmF$Sol)))
vmVarF1 <- map(mmF_list, ~.x %*% t(mmF$X)) %>% map_dbl(~var(as.vector(.x)))
near(vmVarF, vmVarF1)
# getting R2m
R2m<-vmVarF/(vmVarF+mmF$VCV[,1]+mmF$VCV[,2]+mmF$VCV[,3])
mean(R2m)
posterior.mode(R2m)
HPDinterval(R2m)
# R2GLMM(c) - conditional R2GLMM for full model
# Equ. 30
(VarF + VarCorr(mF)$Container[1] + VarCorr(mF)$Population[1])/(VarF + VarCorr(mF)$Container[1] + VarCorr(mF)$Population[1] + (attr(VarCorr(mF), "sc")^2))

# MCMCglmm - conditional
(mVarF+sum(apply(mmF$VCV,2,mean)[-3]))/(mVarF+sum(apply(mmF$VCV,2,mean)))
# alternative with crebile intervals
R2c<-(vmVarF+mmF$VCV[,1]+mmF$VCV[,2])/(vmVarF+mmF$VCV[,1]+mmF$VCV[,2]+mmF$VCV[,3])
mean(R2c)
posterior.mode(R2c)
HPDinterval(R2c)

# 2. Analysis of colour morphs (Binomial mixed models)
#---------------------------------------------------

# Clear memory
rm(list = ls())
# Read colour morph data (Binary, available for males only)
Data <- read.csv(here("data" , "BeetlesMale.csv"))

# Fit null model without fixed effects (but including all random effects)
m0 <- glmer(Colour ~ 1 + (1 | Population) + (1 | Container), family = "binomial", data = Data,
            #control=glmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))
            )

# MCMCglmm 
# binomial models are quite special so you need to do a few things: 1) prior and 2) transformation after analyses
# prior 
prior1 <-list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                       G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)),
                      R = list(V=1,fix=1))
# model
mm0 <- MCMCglmm(Colour ~ 1, random = ~ Population + Container, 
                prior=prior1, data = Data, family="categorical")

c2 <- (16 * sqrt(3)/(15 * pi))^2
mm0$Sol1 <- mm0$Sol/sqrt(1+c2 * mm0$VCV[, "units"]) 
mm0$VCV1 <- mm0$VCV/(1+c2 * mm0$VCV[, "units"]) 

# model is singular - probably variance components are unreliable
summary(m0)
# MCMCglmm 
summary(mm0)
summary(mm0$Sol1)
summary(mm0$VCV1)

# Fit alternative model including fixed and all random effects
mF <- glmer(Colour ~ Treatment + Habitat + (1 | Population) + (1 | Container), 
            family = "binomial", data = Data)

# MCMCglmm 
# prior
prior2 <-list(B=list(mu= c(0,0,0), V=gelman.prior(~Treatment + Habitat, data = Data, scale = sqrt(pi^2/3+1))),
              G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                       G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)),
              R = list(V=1,fix=1))

# View model fits for models with fixed effects
mmF <- MCMCglmm(Colour ~ Treatment + Habitat , random = ~ Population + Container, 
                prior=prior2, data = Data, family="categorical")

c2 <- (16 * sqrt(3)/(15 * pi))^2
mmF$Sol1 <- mmF$Sol/sqrt(1+c2 * mmF$VCV[, "units"]) 
mmF$VCV1 <- mmF$VCV/(1+c2 * mmF$VCV[, "units"]) 

# lmer
summary(mF)


# MCMCglmm
summary(mmF)
summary(mmF$Sol1)
summary(mmF$VCV1)

# Extraction of fitted value for the alternative model 
# fixef() extracts coefficents for fixed effects 
# mF@pp$X returns fixed effect design matrix
Fixed <- fixef(mF)[2] * mF@pp$X[, 2] + fixef(mF)[3] * mF@pp$X[, 3]

# MCMCglmm (it is probably better to get a posterior distribuiton of R2 rather than getting each varaince component - we do this below as an alternative)
mFixed <- mean(mmF$Sol1[,2]) * mmF$X[, 2] + mean(mmF$Sol1[, 3]) * mmF$X[, 3] 

# Calculation of the variance in fitted values
VarF <- var(Fixed)
mVarF<- var(mFixed)
# An alternative way for getting the same result
VarF <- var(as.vector(fixef(mF) %*% t(mF@pp$X)))
mVarF <- var(as.vector(apply(mmF$Sol,2,mean) %*% t(mmF$X)))

# MCMCglmm - marginal
mVarF/(mVarF+sum(apply(mmF$VCV1,2,mean)))
# alternative with crebile intervals
vmVarF<-numeric(1000)
for(i in 1:1000){
  Var<-var(as.vector(mmF$Sol1[i,] %*% t(mmF$X)))
  vmVarF[i]<-Var}
# quicker way to do this
mmF_list <- split(mmF$Sol1, seq(nrow(mmF$Sol1)))
vmVarF1 <- map(mmF_list, ~.x %*% t(mmF$X)) %>% map_dbl(~var(as.vector(.x)))
near(vmVarF, vmVarF1)

# R2GLMM(m) - marginal R2GLMM
# see Equ. 29 and 30 and Table 2
VarF/(VarF + VarCorr(mF)$Container[1] + VarCorr(mF)$Population[1] + pi^2/3)

# getting R2m
(mVarF)/(mVarF+sum(apply(mmF$VCV1,2,mean)[-3])+ pi^2/3)
R2m<-vmVarF/(vmVarF+mmF$VCV1[,1]+mmF$VCV1[,2] + pi^2/3)
mean(R2m)
posterior.mode(R2m)
HPDinterval(R2m)


# R2GLMM(c) - conditional R2GLMM for full model
# Equ. 30
(VarF + VarCorr(mF)$Container[1] + 
    VarCorr(mF)$Population[1])/(VarF + VarCorr(mF)$Container[1] + VarCorr(mF)$Population[1] + pi^2/3)

# MCMCglmm - conditional
(mVarF+sum(apply(mmF$VCV1,2,mean)[-3]))/(mVarF+sum(apply(mmF$VCV1,2,mean)[-3])+ pi^2/3)
# alternative with crebile intervals
R2c<-(vmVarF+mmF$VCV1[,1]+mmF$VCV1[,2])/(vmVarF+mmF$VCV1[,1]+mmF$VCV1[,2]+ pi^2/3)
mean(R2c)
posterior.mode(R2c)
HPDinterval(R2c)




# 3. Analysis of fecundity (Poisson mixed models)
#---------------------------------------------------

# Clear memory
rm(list = ls())

# Read fecundity data (Poisson, available for females only)
Data <- read.csv(here("data","BeetlesFemale.csv"))

# Creating a dummy variable that allows estimating additive dispersion in lmer 
# This triggers a warning message when fitting the model
Unit <- factor(1:length(Data$Egg))

# Fit null model without fixed effects (but including all random effects)
m0 <- glmer(Egg ~ 1 + 
              (1 | Population) + (1 | Container) + (1 | Unit), family = "poisson", data = Data)

# MCMCglmm 
# binomial models are quite special so you need to do a few things: 1) prior and 2) transformation after analyses
# prior 
prior1 <-list(G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                       G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)),
              R = list(V = 1, nu = 0.002))
# model
mm0 <- MCMCglmm(Egg ~ 1, random = ~ Population + Container, 
                prior=prior1, data = Data, family="poisson")

# lmer
summary(m0)
# MCMCglmm 
summary(mm0)

# Fit alternative model including fixed and all random effects
mF <- glmer(Egg ~ Treatment + Habitat + 
              (1 | Population) + (1 | Container) + (1 | Unit), family = "poisson", data = Data)


# View model fits for models with fixed effects
mmF <- MCMCglmm(Egg ~ Treatment + Habitat, random = ~ Population + Container, 
                prior=prior1, data = Data, family="poisson")

# lmer
summary(mF)


# MCMCglmm
summary(mmF)

# Extraction of fitted value for the alternative model 
# fixef() extracts coefficents for fixed effects 
# mF@pp$X returns fixed effect design matrix
Fixed <- fixef(mF)[2] * mF@pp$X[, 2] + fixef(mF)[3] * mF@pp$X[, 3]
# MCMCglmm (it is probably better to get a posterior distribuiton of R2 rather than getting each varaince component - we do this below as an alternative)
mFixed <- mean(mmF$Sol[,2]) * mmF$X[, 2] + mean(mmF$Sol[, 3]) * mmF$X[, 3] + mean(mmF$Sol[ ,4])
# Calculation of the variance in fitted values
VarF <- var(Fixed)
mVarF<- var(mFixed)
# An alternative way for getting the same result
VarF <- var(as.vector(fixef(mF) %*% t(mF@pp$X)))
mVarF <- var(as.vector(apply(mmF$Sol,2,mean) %*% t(mmF$X)))
# R2GLMM(m) - marginal R2GLMM 
# see Equ. 29 and 30 and Table 2 
# fixef(m0) returns the estimate for the intercept of null model
VarF/(VarF + VarCorr(mF)$Container[1] + VarCorr(mF)$Population[1] + VarCorr(mF)$Unit[1] + log(1 + 1/exp(as.numeric(fixef(m0)))))
# fixed - see
# Nakagawa, Shinichi, Paul CD Johnson, and Holger Schielzeth. "The coefficient of determination R 2 and intra-class correlation coefficient from generalized linear mixed-effects models revisited and expanded." Journal of the Royal Society Interface 14.134 (2017): 20170213.
Vdis <- log(1 + 1/exp(as.numeric(fixef(m0) + 0.5*(VarCorr(mF)$Container[1] + VarCorr(mF)$Population[1] + VarCorr(mF)$Unit[1]))))

VarF/(VarF + VarCorr(mF)$Container[1] + VarCorr(mF)$Population[1] + VarCorr(mF)$Unit[1] + Vdis)

# MCMCglmm - marginal
mVdis <- log(1 + 1/exp(mm0$Sol[,1] + 0.5*(rowSums(mmF$VCV))))
mVarF/(mVarF+sum(apply(mmF$VCV,2,mean)) + mean(mVdis))
# alternative with crebile intervals
vmVarF<-numeric(1000)
for(i in 1:1000){
  Var<-var(as.vector(mmF$Sol[i,] %*% t(mmF$X)))
  vmVarF[i]<-Var}
# quicker way to do this
mmF_list <- split(mmF$Sol, seq(nrow(mmF$Sol)))
vmVarF1 <- map(mmF_list, ~.x %*% t(mmF$X)) %>% map_dbl(~var(as.vector(.x)))
near(vmVarF, vmVarF1)
# getting R2m
R2m<-vmVarF/(vmVarF+mmF$VCV[,1]+mmF$VCV[,2]+mmF$VCV[,3] + mean(mVdis))
mean(R2m)
posterior.mode(R2m)
HPDinterval(R2m)

# R2GLMM(c) - conditional R2GLMM for full model
# Equ. 30
(VarF + VarCorr(mF)$Container[1] + 
    VarCorr(mF)$Population[1])/(VarF + VarCorr(mF)$Container[1] + 
                                  VarCorr(mF)$Population[1] + VarCorr(mF)$Unit[1] + 
                                  log(1 + 1/exp(as.numeric(fixef(m0)))))

# the same as above
(VarF + VarCorr(mF)$Container[1] + 
    VarCorr(mF)$Population[1])/(VarF + VarCorr(mF)$Container[1] + 
                                  VarCorr(mF)$Population[1] + VarCorr(mF)$Unit[1] + Vdis)

# MCMCglmm - conditional
(mVarF+sum(apply(mmF$VCV,2,mean)[-3]))/(mVarF+sum(apply(mmF$VCV,2,mean))+ mean(mVdis))
# alternative with crebile intervals
R2c<-(vmVarF+mmF$VCV[,1]+mmF$VCV[,2])/(vmVarF+mmF$VCV[,1]+mmF$VCV[,2]+mmF$VCV[,3]+ mVdis)
mean(R2c)
posterior.mode(R2c)
HPDinterval(R2c)