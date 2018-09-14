---
title: "README"
author: "Zachary McCaw"
date: "2018-09-14"
output: 
  html_document: 
    keep_md: TRUE
--- 

# Package Vignette




## Contents

* [Introduction](#introduction)
* [Logistic](#logistic)
* [Probit](#probit)
* [Robit](#robit)

# Introduction

This package provides implementations of regression models for a binary outcome. Throughout denote the outcome $y_{i}\in\{0,1\}$, the covariates $x_{i}\in\mathbb{R}^{p}$, and the linear predictor $\eta_{i} = x_{i}'\beta$.

# Logistic

## Model

In logistic regression, the latent variable $z_{i}$ follows a logistic distribution with center $\eta_{i}$ and unit scale:

$$
z_{i}|x_{i} \sim \text{Logistic}\big(\eta_{i},1) \\
y_{i} = I[z_{i}>0]
$$

## Implementation

Below, data are simulated for $n=10^{3}$ subjects. The linear predictor depends on four independent, standard normal covariates. The function `rBinReg` is provided for simulating from binary regression models. The function `fit.BinReg` with `model="logistic"` specifies a logistic model. 


```r
set.seed(100);
library(BinReg);
# Subjects
n = 1e3;
# Design matrix
X = matrix(rnorm(n=4*n),nrow=n);
# Regression coefficient
b = c(1,-0.5,-0.5,0);
# Logistic outcome
y = rBinReg(X,b,model="logistic");
# Estimate regression parameters
M = fit.BinReg(y=y,X=X,model="logistic");
show(M);
```

```
## Objective increment:  75.1 
## Objective increment:  4.5 
## Objective increment:  0.0513 
## Objective increment:  8.84e-06 
## Objective increment:  2.27e-13 
## 4 update(s) performed before tolerance limit.
## 
## Fitted Logistic Model
## Estimated Coefficients:
##   Coeff   Point     SE       L      U        p
## 1    x1  1.0000 0.0849  0.8370  1.170 3.12e-32
## 2    x2 -0.5480 0.0794 -0.7030 -0.392 5.32e-12
## 3    x3 -0.4870 0.0728 -0.6300 -0.344 2.31e-11
## 4    x4  0.0998 0.0757 -0.0487  0.248 1.88e-01
```

# Probit

## Model

In probit regression, the latent variable $z_{i}$ follows a normal distribution with center $\eta_{i}$ and unit scale:

$$
z_{i}|x_{i} \sim N\big(\eta_{i},1) \\
y_{i} = I[z_{i}>0]
$$

## Implementation

The function `fit.BinReg` with `model="probit"` specifies a probit model. 


```r
# Probit outcome
y = rBinReg(X,b,model="probit");
# Estimate regression parameters
M = fit.BinReg(y=y,X=X,model="probit");
show(M);
```

```
## Objective increment:  106 
## Objective increment:  16.9 
## Objective increment:  1.21 
## Objective increment:  0.0102 
## Objective increment:  2.75e-06 
## Objective increment:  3.79e-10 
## 5 update(s) performed before tolerance limit.
## 
## Fitted Logistic Model
## Estimated Coefficients:
##   Coeff   Point     SE       L      U        p
## 1    x1  1.0600 0.0661  0.9310  1.190 5.73e-58
## 2    x2 -0.4670 0.0544 -0.5730 -0.360 9.00e-18
## 3    x3 -0.4770 0.0505 -0.5760 -0.378 3.65e-21
## 4    x4  0.0698 0.0506 -0.0293  0.169 1.67e-01
```

# Robit

## Model

In robit regression, the latent variable $z_{i}$ follows a $t_{\nu}$ distribution with center $\eta_{i}$, unit scale, and specified degrees of freedom $\nu$:

$$
z_{i}|x_{i} \sim t_{\nu}(\eta_{i}) \\
y_{i} = I[z_{i}>0]
$$

## Implementation

The function `fit.BinReg` with `model="robit"` specifies a robit model. 


```r
# Degrees of freedom
df = 7;
# Robit outcome
y = rBinReg(X,b,model="robit",df=df);
# Recover regression parameters
M = fit.BinReg(y=y,X=X,model="robit",df=df);
show(M);
```

```
## Objective increment:  100 
## Objective increment:  14.9 
## Objective increment:  0.802 
## Objective increment:  0.0036 
## Objective increment:  2.85e-07 
## Objective increment:  1.39e-11 
## 5 update(s) performed before tolerance limit.
## 
## Fitted Logistic Model
## Estimated Coefficients:
##   Coeff   Point     SE      L       U        p
## 1    x1  1.0900 0.0753  0.944  1.2400 1.22e-47
## 2    x2 -0.5370 0.0612 -0.657 -0.4170 1.68e-18
## 3    x3 -0.4340 0.0544 -0.540 -0.3270 1.52e-15
## 4    x4 -0.0938 0.0551 -0.202  0.0141 8.85e-02
```
