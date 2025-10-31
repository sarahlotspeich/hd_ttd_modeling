# data generation script for function and metric demonstrations

library(tidyverse)
library(truncnorm)

###############################################################################
# Each variable will be drawn from a censoring-rate guided mixture of two 
# distributions: one for subjects diagnosed during study and one for subjects 
# not diagnosed (censored) during study. Distribution parameters and weights
# will be guided by true values found in ENROLL-HD study data.
###############################################################################
set.seed(2025)
n <- 10
cens.rate <- 0.78 # ENROLL-HD censoring rate after filtration ~ 78%

## covariate-driving censoring likelihoods - NOT actual censoring indicators
r_diag <- function(n, cens.rate) {rbinom(n, 1, cens.rate)}
diag <- r_diag(n, cens.rate)
n0 <- sum(diag) # number of subjects (out of n) censored
n1 <- n - n0 # number of subjects (out of n) observed

###############################################################################
# covariate generation given censoring stratum r_diag
###############################################################################

## age at study start
r_age <- function(n, diag, mu = c(36, 44), sd = c(10, 11), min = 18, max = 80) {
  ages <- c(rep(0,n))
  n0 <- sum(diag)
  n1 <- n - n0
  ages[diag == 0] <- rtruncnorm(n0, a = min, b = max, mean = mu[1], sd = sd[1]) # censored
  ages[diag == 1] <-  rtruncnorm(n1, a = min, b = max, mean = mu[2], sd = sd[2]) # observed diagnosis
  return(ages)
}
age <- r_age(n, diag)

## CAG repeats
r_cag <- function(n, diag, mean = c(42, 44), size = c(29,29)) {
  cag <- c(rep(0,n))
  n0 <- sum(diag)
  n1 <- n - n0
  cag[diag == 0] <- rbinom(n0, size = size[1], prob = (mean[1]-36)/size[1]) # censored
  cag[diag == 1] <- rbinom(n1, size = size[2], prob = (mean[2]-36)/size[2]) # observed diagnosis
  return(cag+36)
}
cag <- r_cag(n, diag)

## sex
r_sex <- function(n, diag, prop = c(0.40, 0.46)) {
  sex <- logical(n)
  n0 <- sum(diag)
  n1 <- n - n0
  sex[diag == 0] <- ifelse(runif(n0) <= prop[1], "male", "female") # censored
  sex[diag == 1] <- ifelse(runif(n1) <= prop[2], "male", "female") # observed diagnosis
  return(factor(sex, labels = c("male", "female")))
}
sex <- r_sex(n, diag)

## DCL
r_dcl <- function(n, diag, prob0 = c(65, 25, 5, 5), prob1 = c(15, 25, 25, 30)) {
  dcl <- c(rep(0,n))
  n0 <- sum(diag)
  n1 <- n - n0
  levels = 0:3
  dcl[diag == 0] <- sample(levels, n0, replace = TRUE, prob = prob0) # censored
  dcl[diag == 1] <- sample(levels, n1, replace = TRUE, prob = prob1) # observed diagnosis
  return(dcl)
}
dcl <- r_dcl(n, diag)

## TMS
r_tms <- function(n, diag, mu = c(2.5, 9.5), sd = c(4, 9), min = 0, max = 124) {
  tms <- c(rep(0, n))
  n0 <- sum(diag)
  n1 <- n - n0
  tms[diag == 0] <- rtruncnorm(n0, a = min, b = max, mean = mu[1], sd = sd[1]) # censored
  tms[diag == 1] <- rtruncnorm(n1, a = min, b = max, mean = mu[2], sd = sd[2]) # observed diagnosis
  return(tms)
}
tms <- r_tms(n, diag)

## SDMT
r_sdmt <- function(n, diag, mu = c(51, 39), sd = c(12, 13), min = 0, max = 110) {
  sdmt <- c(rep(0, n))
  n0 <- sum(diag)
  n1 <- n - n0
  sdmt[diag == 0] <- rtruncnorm(n0, a = min, b = max, mean = mu[1], sd = sd[1]) # censored
  sdmt[diag == 1] <- rtruncnorm(n1, a = min, b = max, mean = mu[2], sd = sd[2]) # observed diagnosis
  return(sdmt)
}
sdmt <- r_sdmt(n, diag)

## Stroop word
r_sword <- function(n, diag, mu = c(94, 79), sd = c(18, 20), min = 0, max = 112) {
  sword <- c(rep(0, n))
  n0 <- sum(diag)
  n1 <- n - n0
  sword[diag == 0] <- rtruncnorm(n0, a = min, b = max, mean = mu[1], sd = sd[1]) # censored
  sword[diag == 1] <- rtruncnorm(n1, a = min, b = max, mean = mu[2], sd = sd[2]) # observed diagnosis
  return(sword)
}
sword <- r_sword(n, diag)

## Stroop color
r_scolor <- function(n, diag, mu = c(73, 62), sd = c(15, 16), min = 0, max = 112) {
  scolor <- c(rep(0, n))
  n0 <- sum(diag)
  n1 <- n - n0
  scolor[diag == 0] <- rtruncnorm(n0, a = min, b = max, mean = mu[1], sd = sd[1]) # censored
  scolor[diag == 1] <- rtruncnorm(n1, a = min, b = max, mean = mu[2], sd = sd[2]) # observed diagnosis
  return(scolor)
}
scolor <- r_scolor(n, diag)


## Stroop interference
r_sint <- function(n, diag, mu = c(45, 36), sd = c(11, 12), min = 0, max = 112) {
  sint <- c(rep(0, n))
  n0 <- sum(diag)
  n1 <- n - n0
  sint[diag == 0] <- rtruncnorm(n0, a = min, b = max, mean = mu[1], sd = sd[1]) # censored
  sint[diag == 1] <- rtruncnorm(n1, a = min, b = max, mean = mu[2], sd = sd[2]) # observed diagnosis
  return(sint)
}
sint <- r_sint(n, diag)

###############################################################################
# design matrix generation
###############################################################################
## combined covariates
n = 1000
cens.rate = 0.78
dat <- data.frame(diag = r_diag(n, cens.rate),
                  age_0 = r_age(n, diag),
                  CAG = r_cag(n, diag),
                  sex = r_sex(n, diag),
                  diagconf = r_dcl(n, diag),
                  motscore = r_tms(n, diag),
                  sdmt1 = r_sdmt(n, diag),
                  swrt1 = r_sword(n, diag),
                  scnt1 = r_scolor(n, diag),
                  sit1 = r_sint(n, diag)
)

dat <- dat %>%
  mutate(motscore2 = motscore^2,
         DCL1 = as.integer(diagconf == 1),
         DCL2 = as.integer(diagconf == 2),
         DCL3 = as.integer(diagconf == 3),
         cagMotscore = CAG*motscore,
         cagAge = CAG*age_0)
design.mat <- model.matrix( ~ 1 + age_0 + sdmt1 + CAG + motscore + motscore2 +
                              cagMotscore + cagAge + swrt1 + scnt1 + sit1 +
                              DCL1 + DCL2 + DCL3, data = dat)

###############################################################################
# outcome generation
###############################################################################
## outcome time-to-diagnosis data generation
beta <- c(
  `(Intercept)` = -3.2,
  age_0 = 0.02,
  sdmt1 = -0.015,
  CAG = 0.08,
  motscore = 0.03,
  motscore2 = 0.0005,
  cagMotscore = 0.002,
  cagAge = 0.001,
  swrt1 = -0.004,
  scnt1 = -0.003,
  sit1 = -0.005,
  DCL1 = 0.15,
  DCL2 = 0.30,
  DCL3 = 0.50
)

eta <- drop(design.mat %*% beta)

shape <- 1.7 # Weibull shape parameter
scale <- 3 # weibull scale parameter

# generate uniform random variable for inverse transform
U <- runif(n)

# generate event times X using the proportional hazards model
X <- scale * (-log(U) / exp(eta))^(1 / shape)

# generate censoring times
C <- rexp(n = n, rate = cens.rate)

# Observed time and censoring indicator
W <- pmin(X, C)
delta <- as.numeric(X <= C)
sum(delta)

toydat <- dat %>%
  mutate(X = X,
         C = C,
         W = W,
         delta = delta)
  