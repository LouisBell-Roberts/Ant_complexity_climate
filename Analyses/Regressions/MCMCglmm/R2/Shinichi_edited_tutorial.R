########Edited code for calculating R2 - Shinichi Nakagawa
##Useless sections from the original code removed
#28/10/2023


# 1. Analysis of body size (Gaussian mixed models)

# Read body length data (Gaussian, available for both sexes)
Data <- read.csv(here("data","BeetlesBody.csv"))
# Fit null model without fixed effects (but including all random effects)
# prior 
prior <- list(R = list(V = 1, nu = 0.002),
              G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                       G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

# model
mm0<-MCMCglmm(BodyL ~ 1 , random = ~ Population + Container, prior = prior, data=Data)

# Fit alternative model including fixed and all random effects
mmF<-MCMCglmm(BodyL ~ Sex + Treatment + Habitat , random = ~ Population + Container, prior = prior, data=Data)

# View model fits for both models
summary(mm0) # MCMCglmm
summary(mmF) # MCMCglmm

# Extraction of fitted value for the alternative model
# fixef() extracts coefficents for fixed effects
# mF@pp$X returns fixed effect design matrix
# MCMCglmm (it is probably better to get a posterior distribuiton of R2 rather than getting each varaince component - we do this below as an alternative)
mFixed <- mean(mmF$Sol[,2]) * mmF$X[, 2] + mean(mmF$Sol[, 3]) * mmF$X[, 3] + mean(mmF$Sol[ ,4]) * mmF$X[, 4]
# Calculation of the variance in fitted values
mVarF<- var(mFixed)

# An alternative way for getting the same result
mVarF <- var(as.vector(apply(mmF$Sol,2,mean) %*% t(mmF$X)))
# R2GLMM(m) - marginal R2GLMM
# Equ. 26, 29 and 30
# VarCorr() extracts variance components
# attr(VarCorr(lmer.model),'sc')^2 extracts the residual variance

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

# MCMCglmm - conditional
(mVarF+sum(apply(mmF$VCV,2,mean)[-3]))/(mVarF+sum(apply(mmF$VCV,2,mean)))
# alternative with crebile intervals
R2c<-(vmVarF+mmF$VCV[,1]+mmF$VCV[,2])/(vmVarF+mmF$VCV[,1]+mmF$VCV[,2]+mmF$VCV[,3])
mean(R2c)
posterior.mode(R2c)
HPDinterval(R2c)



# 2. Analysis of colour morphs (Binomial mixed models)
#---------------------------------------------------

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

# MCMCglmm 
summary(mm0)
summary(mm0$Sol1)
summary(mm0$VCV1)

# Fit alternative model including fixed and all random effects
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

# MCMCglmm
summary(mmF)
summary(mmF$Sol1)
summary(mmF$VCV1)

# Extraction of fitted value for the alternative model 

# MCMCglmm (it is probably better to get a posterior distribution of R2 rather than getting each variance component - we do this below as an alternative)
mFixed <- mean(mmF$Sol1[,2]) * mmF$X[, 2] + mean(mmF$Sol1[, 3]) * mmF$X[, 3] 

# Calculation of the variance in fitted values
mVarF<- var(mFixed)
# An alternative way for getting the same result
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

# getting R2m
(mVarF)/(mVarF+sum(apply(mmF$VCV1,2,mean)[-3])+ pi^2/3)
R2m<-vmVarF/(vmVarF+mmF$VCV1[,1]+mmF$VCV1[,2] + pi^2/3)
mean(R2m)
posterior.mode(R2m)
HPDinterval(R2m)


# R2GLMM(c) - conditional R2GLMM for full model
# Equ. 30

# MCMCglmm - conditional
(mVarF+sum(apply(mmF$VCV1,2,mean)[-3]))/(mVarF+sum(apply(mmF$VCV1,2,mean)[-3])+ pi^2/3)
# alternative with crebile intervals
R2c<-(vmVarF+mmF$VCV1[,1]+mmF$VCV1[,2])/(vmVarF+mmF$VCV1[,1]+mmF$VCV1[,2]+ pi^2/3)
mean(R2c)
posterior.mode(R2c)
HPDinterval(R2c)






