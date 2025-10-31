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

We can fit all four discussed models for time to diagnosis. The Langbehn
model does not directly model time to diagnosis, but age at diagnosis,
so some wrapper functions are necessary.

``` r
# source("Langbehn_functions.R")    # build functions for Langbehn model

# CAP model: Accelerated failure time model with age at study start and CAG repeat number (and their interaction) as covariates
CAP_model <- survreg(Surv(W, delta) ~ age_0 + CAG:age_0,
                     data = sim_data, dist = "loglogistic")

# MRS model: Cox proportional hazards model with DCL, motor score, Stroop number, Stroop word, Stroop interference, SDMT, CAG repeat number, age at study start as covariates
MRS_model <- coxph(Surv(W, delta) ~ factor(diagconf) + motscore + motscore2 +
                     scnt1 + swrt1 + sit1 + sdmt1 + CAG + age_0 +
                     CAG:motscore + CAG:age_0,
                   data = sim_data, x = TRUE)

# PI model: Cox proportional hazards model with motor score, SDMT, and CAG repeat number age product as covariates 
PI_model <- coxph(Surv(W, delta) ~ motscore + sdmt1 + CAP,
                  data = sim_data, x = TRUE)

# Langbehn model: parametric survival with baseline age and CAG repeat number as covariates
Langbehn_params <- tryCatch(
fit_langbehn_model(sim_data, decay_vals = c(0.13, 0.14)))
```

    ## Best negloglik = 0

# Uno’s C

``` r
# source("Uno-C-functions.R")

# set random seed
set.seed(2025)
n_splits = 5

# create training-testing folds, with stratification by delta to preserve censoring rate
folds <- vfold_cv(sim_data, v = n_splits, strata = delta)

# summarize censoring information for each training and testing split
fold_censoring_summary <- map_dfr(1:n_splits, function(i) {
  test_data <- assessment(folds$splits[[i]])
  event_rate <- mean(test_data$delta == 1)
  tibble(
    fold = i,
    event_rate = event_rate,
    censoring_rate = 1 - event_rate,
    n = nrow(test_data),
    n_events = sum(test_data$delta == 1)
  )
})
print(fold_censoring_summary)
```

    ## # A tibble: 5 × 5
    ##    fold event_rate censoring_rate     n n_events
    ##   <int>      <dbl>          <dbl> <int>    <int>
    ## 1     1      0.717          0.283   184      132
    ## 2     2      0.717          0.283   184      132
    ## 3     3      0.717          0.283   184      132
    ## 4     4      0.717          0.283   184      132
    ## 5     5      0.716          0.284   183      131

``` r
# preallocate output for models
fitted_models_by_fold <- vector("list", length = n_splits)

for (i in 1:n_splits) {
  # training-testing split i
  split <- folds$splits[[i]]
  
  # training set i
  train_data <- analysis(split)
  
  # testing set i
  test_data  <- assessment(split)

  # CAP model
  CAP_model <- survreg(Surv(W, delta) ~ age_0 + CAG:age_0,
                       data = train_data, dist = "loglogistic")

  # MRS model
  MRS_model <- coxph(Surv(W, delta) ~ factor(diagconf) +
                       motscore + motscore2 +
                       scnt1 + swrt1 + sit1 + sdmt1 + CAG + age_0 +
                       CAG:motscore + CAG:age_0,
                     data = train_data, x = TRUE)

  # PI model
  PI_model <- coxph(Surv(W, delta) ~ motscore + sdmt1 + CAP,
                    data = train_data, x = TRUE)

  Langbehn_params <- tryCatch(
  fit_langbehn_model(train_data, decay_vals = c(0.13, 0.14)),
  error = function(e) {
    message(paste("Langbehn fit failed on fold", i, ":", e$message))
    return(NULL)
    }
  )

  # Linear predictors (negate for AFT to maintain directionality)
  lp_CAP <- -predict(CAP_model, newdata = test_data, type = "lp")
  lp_MRS <-  predict(MRS_model, newdata = test_data, type = "lp")
  lp_PI  <-  predict(PI_model,  newdata = test_data, type = "lp")
  
  # Langbehn-based risk score (CDF of diagnosis in next 5 years)
  risk_Langbehn <- Langbehn_cond_CDF(Langbehn_params,
                                     current_age = test_data$age_0,
                                     time_to_diag = 5,
                                     cag = test_data$CAG)
  
  # Store models, linear predictors, and training and testing data
  fitted_models_by_fold[[i]] <- list(
    CAP = CAP_model,
    MRS = MRS_model,
    PI = PI_model,
    Langbehn = Langbehn_params,
    lp_CAP = lp_CAP,
    lp_MRS = lp_MRS,
    lp_PI  = lp_PI,
    lp_Langbehn = risk_Langbehn,
    train_data = train_data,
    test_data = test_data
  )
}
```

    ## Best negloglik = 0
    ## Best negloglik = 0
    ## Best negloglik = 0
    ## Best negloglik = 0
    ## Best negloglik = 0

``` r
# Extract Langbehn parameter estimates across folds
langbehn_params_all <- map(fitted_models_by_fold, ~ .x$Langbehn) %>%
  do.call(rbind, .) %>%
  as.data.frame()

# Name the parameters for clarity
colnames(langbehn_params_all) <- c("mu_intercept", "mu_scale", "mu_decay",
                                   "sigma_intercept", "sigma_scale", "sigma_decay")

# Show the parameter estimates across folds
# print(langbehn_params_all)

custom_uno_results <- map_dfr(1:5, function(fold_idx) {
  models <- fitted_models_by_fold[[fold_idx]]
  train_data <- models$train_data
  test_data  <- models$test_data

  # Fit KM for censoring distribution from training data
  km_cens <- survfit(Surv(W, 1 - delta) ~ 1, data = train_data)
  G_hat_left <- function(t) {
    s <- summary(km_cens, times = unique(c(0, km_cens$time)), extend = TRUE)
    approx(
      x = s$time,
      y = s$surv,
      xout = t - 1e-8,
      method = "constant",
      f = 0,
      rule = 2
    )$y
  }

  # Prepare output across models and t_star
  bind_rows(
    tibble(model = "Langbehn", lp = models$lp_Langbehn),
    tibble(model = "CAP", lp = models$lp_CAP),
    tibble(model = "MRS", lp = models$lp_MRS),
    tibble(model = "PIN", lp = models$lp_PI)
  ) %>%
    mutate(fold = fold_idx) %>%
    group_by(model, fold) %>%
    group_modify(~ {
      test_data$lp <- .x$lp
      map_dfr(1:8, function(t_star) {
        tibble(
          time = t_star,
          uno_c_custom = compute_uno_c(test_data, G_hat_left, t_star)
        )
      })
    })
})

print(custom_uno_results)
```

    ## # A tibble: 160 × 4
    ## # Groups:   model, fold [20]
    ##    model     fold  time uno_c_custom
    ##    <chr>    <int> <int>        <dbl>
    ##  1 CAP          1     1        0.703
    ##  2 CAP          1     2        0.699
    ##  3 CAP          1     3        0.699
    ##  4 CAP          1     4        0.699
    ##  5 CAP          1     5        0.699
    ##  6 CAP          1     6        0.699
    ##  7 CAP          1     7        0.699
    ##  8 CAP          1     8        0.699
    ##  9 Langbehn     1     1        0.694
    ## 10 Langbehn     1     2        0.690
    ## # ℹ 150 more rows

# ROC

# Sample Enrichment
