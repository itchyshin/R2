# revised on 6 Nov 2013
#A general and simple method for obtaining R2 from generalized linear mixed-effects models
#Shinichi Nakagawa1,2 and Holger Schielzeth3
#1 National Centre of Growth and Development, Department of Zoology, University of Otago, Dunedin, New Zealand
#2 Department of Behavioral Ecology and Evolutionary Genetics, Max Planck Institute for Ornithology, Seewiesen, Germany
#3 Department of Evolutionary Biology, Bielefeld University, Bielefeld, Germany
#Running head: Variance explained by GLMMs
#Correspondence:
#S. Nakagawa; Department of Zoology, University of Otago, 340 Great King Street, Dunedin, 9054, New Zealand
#Tel:	+64 (0)3 479 5046
#Fax:	+64 (0)3 479 7584
#e-mail: shinichi.nakagawa@otago.ac.nz 


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


####################################################
# B. Analysis
####################################################

# 1. Analysis of body size (Gaussian mixed models)
#---------------------------------------------------

# Clear memory
rm(list = ls())

# Read body length data (Gaussian, available for both sexes)
Data <- read.csv("BeetlesBody.csv")

# Fit null model without fixed effects (but including all random effects)
m0 <- lmer(BodyL ~ 1 + (1 | Population) + (1 | Container), data = Data)

# Fit alternative model including fixed and all random effects
mF <- lmer(BodyL ~ Sex + Treatment + Habitat + (1 | Population) + (1 | Container), data = Data)

# View model fits for both models
summary(m0)
summary(mF)

# Extraction of fitted value for the alternative model
# fixef() extracts coefficents for fixed effects
# mF@pp$X returns fixed effect design matrix
Fixed <- fixef(mF)[2] * mF@pp$X[, 2] + fixef(mF)[3] * mF@pp$X[, 3] + fixef(mF)[4] * mF@pp$X[, 4]

# Calculation of the variance in fitted values
VarF <- var(Fixed)

# An alternative way for getting the same result
VarF <- var(as.vector(fixef(mF) %*% t(mF@pp$X)))

# R2GLMM(m) - marginal R2GLMM
# Equ. 26, 29 and 30
# VarCorr() extracts variance components
# attr(VarCorr(lmer.model),'sc')^2 extracts the residual variance
VarF/(VarF + VarCorr(mF)$Container[1] + VarCorr(mF)$Population[1] + attr(VarCorr(mF), "sc")^2)

# R2GLMM(c) - conditional R2GLMM for full model
# Equ. 30
(VarF + VarCorr(mF)$Container[1] + VarCorr(mF)$Population[1])/(VarF + VarCorr(mF)$Container[1] + VarCorr(mF)$Population[1] + (attr(VarCorr(mF), "sc")^2))

# AIC and BIC needs to be calcualted with ML not REML in body size models
m0ML <- lmer(BodyL ~ 1 + (1 | Population) + (1 | Container), data = Data, REML = FALSE)
mFML <- lmer(BodyL ~ Sex + Treatment + Habitat + (1 | Population) + (1 | Container), data = Data, REML = FALSE)

# View model fits for both models fitted by ML
summary(m0ML)
summary(mFML)


# 2. Analysis of colour morphs (Binomial mixed models)
#---------------------------------------------------

# Clear memory
rm(list = ls())
# Read colour morph data (Binary, available for males only)
Data <- read.csv("BeetlesMale.csv")

# Fit null model without fixed effects (but including all random effects)
m0 <- glmer(Colour ~ 1 + (1 | Population) + (1 | Container), family = "binomial", data = Data)

# Fit alternative model including fixed and all random effects
mF <- glmer(Colour ~ Treatment + Habitat + (1 | Population) + (1 | Container), family = "binomial", data = Data)

# View model fits for both models
summary(m0)
summary(mF)

# Extraction of fitted value for the alternative model 
# fixef() extracts coefficents for fixed effects 
# mF@pp$X returns fixed effect design matrix
Fixed <- fixef(mF)[2] * mF@pp$X[, 2] + fixef(mF)[3] * mF@pp$X[, 3]

# Calculation of the variance in fitted values
VarF <- var(Fixed)

# An alternative way for getting the same result
VarF <- var(as.vector(fixef(mF) %*% t(mF@pp$X)))

# R2GLMM(m) - marginal R2GLMM
# see Equ. 29 and 30 and Table 2
VarF/(VarF + VarCorr(mF)$Container[1] + VarCorr(mF)$Population[1] + pi^2/3)

# R2GLMM(c) - conditional R2GLMM for full model
# Equ. 30
(VarF + VarCorr(mF)$Container[1] + VarCorr(mF)$Population[1])/(VarF + VarCorr(mF)$Container[1] + VarCorr(mF)$Population[1] + pi^2/3)


# 3. Analysis of fecundity (Poisson mixed models)
#---------------------------------------------------

# Clear memory
rm(list = ls())

# Read fecundity data (Poisson, available for females only)
Data <- read.csv("BeetlesFemale.csv")

# Creating a dummy variable that allows estimating additive dispersion in lmer 
# This triggers a warning message when fitting the model
Unit <- factor(1:length(Data$Egg))

# Fit null model without fixed effects (but including all random effects)
m0 <- glmer(Egg ~ 1 + (1 | Population) + (1 | Container) + (1 | Unit), family = "poisson", data = Data)

# Fit alternative model including fixed and all random effects
mF <- glmer(Egg ~ Treatment + Habitat + (1 | Population) + (1 | Container) + (1 | Unit), family = "poisson", data = Data)

# View model fits for both models
summary(m0)
summary(mF)

# Extraction of fitted value for the alternative model 
# fixef() extracts coefficents for fixed effects 
# mF@pp$X returns fixed effect design matrix
Fixed <- fixef(mF)[2] * mF@pp$X[, 2] + fixef(mF)[3] * mF@pp$X[, 3]

# Calculation of the variance in fitted values
VarF <- var(Fixed)

# An alternative way for getting the same result
VarF <- var(as.vector(fixef(mF) %*% t(mF@pp$X)))

# R2GLMM(m) - marginal R2GLMM 
# see Equ. 29 and 30 and Table 2 
# fixef(m0) returns the estimate for the intercept of null model
VarF/(VarF + VarCorr(mF)$Container[1] + VarCorr(mF)$Population[1] + VarCorr(mF)$Unit[1] + log(1 + 1/exp(as.numeric(fixef(m0)))))

# R2GLMM(c) - conditional R2GLMM for full model
# Equ. 30
(VarF + VarCorr(mF)$Container[1] + VarCorr(mF)$Population[1])/(VarF + VarCorr(mF)$Container[1] + VarCorr(mF)$Population[1] + VarCorr(mF)$Unit[1] + log(1 + 
    1/exp(as.numeric(fixef(m0)))))


####################################################
# C. Data generation
####################################################

# 1. Design matrices 
#---------------------------------------------------

# Clear memory
rm(list = ls())

# 12 different populations n = 960
Population <- gl(12, 80, 960)

# 120 containers (8 individuals in each container)
Container <- gl(120, 8, 960)

# Sex of the individuals. Uni-sex within each container (individuals are sorted at the pupa stage)
Sex <- factor(rep(rep(c("Female", "Male"), each = 8), 60))

# Habitat at the collection site: dry or wet soil (four indiviudal from each Habitat in each container)
Habitat <- factor(rep(rep(c("dry", "wet"), each = 4), 120))

# Food treatment at the larval stage: special food ('Exp') or standard food ('Cont')
Treatment <- factor(rep(c("Cont", "Exp"), 480))

# Data combined in a dataframe
Data <- data.frame(Population = Population, Container = Container, Sex = Sex, Habitat = Habitat, Treatment = Treatment)


# 2. Gaussian response: body length (both sexes)
#---------------------------------------------------

# simulation of the underlying random effects (Population and Container with variance of 1.3 and 0.3, respectively)
PopulationE <- rnorm(12, 0, sqrt(1.3))
ContainerE <- rnorm(120, 0, sqrt(0.3))

# data generation based on fixed effects, random effects and random residuals errors
Data$BodyL <- 15 - 3 * (as.numeric(Sex) - 1) + 0.4 * (as.numeric(Treatment) - 1) + 0.15 * (as.numeric(Habitat) - 1) + PopulationE[Population] + ContainerE[Container] + 
    rnorm(960, 0, sqrt(1.2))

# save data (to current work directory)
write.csv(Data, file = "BeetlesBody.csv", row.names = F)


# 3. Binomial response: colour morph (males only)
#---------------------------------------------------

# Subset the design matrix (only males express colour morphs)
DataM <- subset(Data, Sex == "Male")

# simulation of the underlying random effects (Population and Container with variance of 1.2 and 0.2, respectively)
PopulationE <- rnorm(12, 0, sqrt(1.2))
ContainerE <- rnorm(120, 0, sqrt(0.2))

# generation of response values on link scale (!) based on fixed effects and random effects
ColourLink <- with(DataM, 0.8 * (-1) + 0.8 * (as.numeric(Treatment) - 1) + 0.5 * (as.numeric(Habitat) - 1) + PopulationE[Population] + ContainerE[Container])

# data generation (on data scale!) based on binomial distribution
DataM$Colour <- rbinom(length(ColourLink), 1, invlogit(ColourLink))

# save data (to current work directory)
write.csv(DataM, file = "BeetlesMale.csv", row.names = F)


# 4. Poisson response: fecundity (females only)
#---------------------------------------------------

# Subset the design matrix (only females express colour morphs)
DataF <- Data[Data$Sex == "Female", ]

# random effects
PopulationE <- rnorm(12, 0, sqrt(0.4))
ContainerE <- rnorm(120, 0, sqrt(0.05))

# generation of response values on link scale (!) based on fixed effects, random effects and residual errors
EggLink <- with(DataF, 1.1 + 0.5 * (as.numeric(Treatment) - 1) + 0.1 * (as.numeric(Habitat) - 1) + PopulationE[Population] + ContainerE[Container] + rnorm(480, 
    0, sqrt(0.1)))

# data generation (on data scale!) based on Poisson distribution
DataF$Egg <- rpois(length(EggLink), exp(EggLink))

# save data (to current work directory)
write.csv(DataF, file = "BeetlesFemale.csv", row.names = F)




