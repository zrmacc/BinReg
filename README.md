---
title: "README"
author: "Zachary McCaw"
date: "2018-08-18"
output: 
  html_document: 
    keep_md: TRUE
--- 

# Package Vignette




## Contents

* [Probit](#probit)
* [Robit](#robit)

# Probit

## Model

Suppose $y_{i}\in\{0,1\}$ is a binary outcome, and $z_{i}\in\mathbb{R}^{p}$ is a set of covariates. Let $\eta_{i} = z_{i}'\beta$ denote the linear predictor. Consider the model:

$$
\zeta_{i}|z_{i} \sim N\big(\eta_{i},1) \\
y_{i} = I[\zeta_{i}>0]
$$

## Implementation


```r
# Subjects
n = 1e3;
# Design matrix
Z = matrix(rnorm(n=n*3),nrow=n);
# Regression coefficient
beta = c(1,-1,1);
# Linear predictor
eta = Z %*% beta;
# Probit outcome
y = 1*(rnorm(n=n,mean=eta)>0);
# Recover regression parameters
M = fit.Probit(y=y,Z=Z);
show(M);
```

```
## 5 update(s) performed before tolerance limit. 
## Fitted Probit Model
## Estimated Coefficients:
##   Coeff    Est     SE      L      U        p
## 1    z1  0.976 0.0700  0.839  1.110 3.74e-44
## 2    z2 -1.010 0.0706 -1.140 -0.867 5.66e-46
## 3    z3  1.080 0.0742  0.931  1.220 1.11e-47
```

# Robit

## Model

Suppose $y_{i}\in\{0,1\}$ is a binary outcome, and $z_{i}\in\mathbb{R}^{p}$ is a set of covariates. Let $\eta_{i} = z_{i}'\beta$ denote the linear predictor. Consider the model:

$$
\zeta_{i}|z_{i} \sim t_{\nu}(\eta_{i}) \\
y_{i} = I[\zeta_{i}>0]
$$

## Implementation


```r
# Degrees of freedom
df = 7;
# Robit outcome
y = 1*(rt(n=n,ncp=eta,df=df)>0);
# Recover regression parameters
M = fit.Robit(y=y,Z=Z);
show(M);
```

```
## 5 update(s) performed before tolerance limit. 
## Fitted Robit Model
## Estimated Coefficients:
##   Coeff   Est     SE      L      U        p
## 1    z1  1.09 0.0827  0.927  1.250 1.35e-39
## 2    z2 -1.06 0.0816 -1.220 -0.898 1.75e-38
## 3    z3  1.11 0.0848  0.948  1.280 2.01e-39
```
