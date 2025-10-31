hd_ttd_modeling
================

R code for the manuscript “Adjusting for Heavy Censoring and
Double-Dipping to Compare Risk Stratification Abilities of Existing
Models for Time to Diagnosis of Huntington Disease” by Grosser, Foes,
Li, Parikh, Garcia, and Lotspeich (2025+)

# Data Generation

We can simulate data from a Huntington disease clinical trial in order
to evaluate various models for time to diagnosis. We generate data to
match observations from the ENROLL-HD longitudinal study of Huntington
disease.

``` r
# source("datagen.R")    # build functions for data generation

n = 1000               # sample size
cens.rate = 0.78       # desired censoring rate (approximate)

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

# define parameters for covariates in outcome model
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
scale <- 3   # weibull scale parameter

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
```

    ## [1] 705

``` r
# combine survival outcome data and covariates for complete dataset
sim_data <- dat %>%
  mutate(X = X,
         C = C,
         W = W,
         delta = delta,
         CAP = age_0 * (CAG - 34)) %>%
  filter(CAG > 40 & CAG < 57)

head(sim_data)
```

    ##   diag    age_0 CAG    sex diagconf  motscore    sdmt1     swrt1    scnt1
    ## 1    1 47.97340  46   male        1 19.815860 36.77801  74.26561 42.91648
    ## 2    0 45.67558  47 female        0  7.909145 32.87784  58.44761 67.41093
    ## 3    0 33.58008  44 female        2  4.650563 45.34272  44.84870 37.00444
    ## 4    1 37.44354  43 female        2  6.889976 54.70478  81.06527 75.93027
    ## 5    1 62.65389  42   male        3  3.713107 34.98918 102.59544 68.03526
    ## 6    1 29.17181  41 female        0  2.594584 44.07967  85.90685 76.16116
    ##       sit1  motscore2 DCL1 DCL2 DCL3 cagMotscore   cagAge         X         C
    ## 1 30.18573 392.668289    1    0    0    911.5295 2206.776 0.2792411 0.4759966
    ## 2 40.12274  62.554569    0    0    0    371.7298 2146.752 0.3971169 1.3083983
    ## 3 38.52677  21.627736    0    1    0    204.6248 1477.523 1.0700101 1.1052819
    ## 4 19.87406  47.471771    0    1    0    296.2690 1610.072 0.5126373 1.6386207
    ## 5 27.26405  13.787160    0    0    1    155.9505 2631.464 0.8494234 0.2988568
    ## 6 33.66016   6.731865    0    0    0    106.3779 1196.044 3.7328927 3.3924873
    ##           W delta      CAP
    ## 1 0.2792411     1 575.6808
    ## 2 0.3971169     1 593.7826
    ## 3 1.0700101     1 335.8008
    ## 4 0.5126373     1 336.9919
    ## 5 0.2988568     0 501.2312
    ## 6 3.3924873     0 204.2027

For simplicity, this simulated data set has no missing data other than
censored times of diagnosis.

# Modeling Time to Diagnosis

``` r
# CAP model: Accelerated failure time model with age at study start and CAG repeat number (and their interaction) as covariates
CAP_model <- survreg(Surv(W, delta) ~ age_0 + CAG:age_0,
                     data = sim_data, dist = "loglogistic")

# MRS model: Cox proportional hazards model with DCL, motor score, Stroop number, Stroop word, Stroop interference, SDMT, CAG repeat number, age at study start as covariates
MRS_model <- coxph(Surv(W, delta) ~ factor(diagconf) + motscore + motscore2 +
                     scnt1 + swrt1 + sit1 + sdmt1 + CAG + age_0 +
                     CAG:motscore + CAG:age_0,
                   data = sim_data, x = TRUE)

# PI model: Cox proportinal hazards model with motor score, SDMT, and CAG repeat number age product as covariates 
PI_model <- coxph(Surv(W, delta) ~ motscore + sdmt1 + CAP,
                  data = sim_data, x = TRUE)

# Langbehn model: parametric survival with baseline age and CAG repeat number as covariates
Langbehn_params <- tryCatch(
fit_langbehn_model(sim_data, decay_vals = c(0.13, 0.14)),
error = function(e) {
  message(paste("Langbehn fit failed on fold", i, ":", e$message))
  return(NULL)
  }
)
```

    ## Best negloglik = 0

# Uno’s C

# ROC

# Sample Enrichment
