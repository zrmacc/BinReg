---
title: "README"
author: "Zachary McCaw"
date: "2018-09-16"
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

## Estimation

Below, data are simulated for $n=10^{3}$ subjects. The model matrix `X` includes an intercept and four standard normal covariates. The function `rBinReg` is provided for simulating from binary regression models. The function `Fit.BinReg` with `model="logistic"` specifies a logistic model. 


```r
set.seed(101);
library(BinReg);
# Subjects
n = 1e3;
# Design matrix
X = cbind(1,matrix(rnorm(n=4*n),nrow=n));
# Regression coefficient
b = c(1,-1,2,-1,0);
# Logistic outcome
y = rBinReg(X,b,model="logistic");
# Estimate regression parameters
M = Fit.BinReg(y=y,X=X,model="logistic");
show(M);
```

```
## Objective increment:  163 
## Objective increment:  28.3 
## Objective increment:  3.23 
## Objective increment:  0.0653 
## Objective increment:  3.39e-05 
## Objective increment:  9.55e-12 
## 5 update(s) performed before tolerance limit.
## 
## Fitted logistic Model
## Estimated Coefficients:
##   Coeff   Point     SE      L      U        p
## 1    x1  1.0200 0.1020  0.825  1.220 7.65e-24
## 2    x2 -1.0000 0.1060 -1.210 -0.794 3.28e-21
## 3    x3  2.0400 0.1390  1.760  2.310 2.26e-48
## 4    x4 -1.0800 0.1080 -1.300 -0.874 6.49e-24
## 5    x5  0.0755 0.0931 -0.107  0.258 4.17e-01
```

## Inference

Wald tests for individual coefficients are stored in `M@Coefficients` and displayed using the `show` method. The function `Score.BinReg` performs score tests on subsets of the regression coefficients. The test is specified using a logical vector `L` with as many elements as columns in the model matrix `X`. An element of `L` is set to `TRUE` if the regression coefficient for that column of `X` is fixed under $H_{0}$. An element of `L` is set to `FALSE` if the regression coefficient for that column requires estimation under $H_{0}$. At least one element of `L` must be `TRUE` (i.e. a test must be specified) and at least one element of `L` must be `FALSE` (i.e. a null model must be estimable).

Below, various hypothses are tested on the example data. The first is an overall test that the coefficients for all covariates, other than the intercept, are zero. The null hypothesis is $H_{0}:\beta_{1}=\beta_{2}=\beta_{3}=\beta_{4}=0$, which is false. The second is an individual test that the coefficient for the fourth covariate is zero. The null hypothesis is $H_{0}:\beta_{4}=0$, which is true. The third is a joint test that the coefficient for the second covariate is zero, and the coefficient for the fourth covariate is two. The null hypothesis is $H_{0}:\beta_{2}=0,\beta_{4}=2$, which if false. The fourth is a joint test that the coefficients for the first and third covariates are equal to negative one. The null hypothesis is $H_{0}:\beta_{1}=\beta_{3}=-1$, which is true. All models include an intercept $\beta_{0}$, since in each case `L[1]=FALSE`. 


```r
cat("Test of b1=b2=b3=b4=0:\n");
Score.BinReg(y=y,X=X,L=c(F,T,T,T,T),model="logistic",report=F);
cat("\n");
cat("Test of b4=0:\n");
Score.BinReg(y=y,X=X,L=c(F,F,F,F,T),model="logistic",report=F);
cat("\n");
cat("Test of b2=0, b4=2:\n");
Score.BinReg(y=y,X=X,L=c(F,F,T,F,T),b10=c(0,2),model="logistic",report=F);
cat("\n");
cat("Test of b1=b3=-1:\n");
Score.BinReg(y=y,X=X,L=c(F,T,F,T,F),b10=c(-1,-1),model="logistic",report=F);
cat("\n");
```

```
## Test of b1=b2=b3=b4=0:
##        Score           df            p 
## 4.271668e+02 4.000000e+00 3.745524e-91 
## 
## Test of b4=0:
##     Score        df         p 
## 0.6583181 1.0000000 0.4171544 
## 
## Test of b2=0, b4=2:
##    Score       df        p 
## 1665.101    2.000    0.000 
## 
## Test of b1=b3=-1:
##     Score        df         p 
## 0.6456440 2.0000000 0.7241027
```

# Probit

## Model

In probit regression, the latent variable $z_{i}$ follows a normal distribution with center $\eta_{i}$ and unit scale:

$$
z_{i}|x_{i} \sim N\big(\eta_{i},1) \\
y_{i} = I[z_{i}>0]
$$

## Estimation

The function `Fit.BinReg` with `model="probit"` specifies a probit model. 


```r
# Probit outcome
y = rBinReg(X,b,model="probit");
# Estimate regression parameters
M = Fit.BinReg(y=y,X=X,model="probit");
show(M);
```

```
## Objective increment:  196 
## Objective increment:  50 
## Objective increment:  17.5 
## Objective increment:  3.38 
## Objective increment:  0.185 
## Objective increment:  0.00134 
## Objective increment:  3.19e-06 
## Objective increment:  8.95e-09 
## 7 update(s) performed before tolerance limit.
## 
## Fitted probit Model
## Estimated Coefficients:
##   Coeff   Point     SE      L       U        p
## 1    x1  1.1000 0.0888  0.929  1.2800 2.06e-35
## 2    x2 -1.0800 0.0913 -1.250 -0.8960 5.16e-32
## 3    x3  2.1800 0.1450  1.900  2.4700 1.58e-51
## 4    x4 -1.0500 0.0892 -1.220 -0.8720 9.40e-32
## 5    x5 -0.0456 0.0684 -0.180  0.0884 5.05e-01
```

## Inference

The same hypotheses considered previously are tested for the probit model. 


```r
cat("Test of b1=b2=b3=b4=0:\n");
Score.BinReg(y=y,X=X,L=c(F,T,T,T,T),model="probit",report=F);
cat("\n");
cat("Test of b4=0:\n");
Score.BinReg(y=y,X=X,L=c(F,F,F,F,T),model="probit",report=F);
cat("\n");
cat("Test of b2=0, b4=2:\n");
Score.BinReg(y=y,X=X,L=c(F,F,T,F,T),b10=c(0,2),model="probit",report=F);
cat("\n");
cat("Test of b1=b3=-1:\n");
Score.BinReg(y=y,X=X,L=c(F,T,F,T,F),b10=c(-1,-1),model="probit",report=F);
cat("\n");
```

```
## Test of b1=b2=b3=b4=0:
##         Score            df             p 
##  5.399379e+02  4.000000e+00 1.537766e-115 
## 
## Test of b4=0:
##     Score        df         p 
## 0.4245797 1.0000000 0.5146607 
## 
## Test of b2=0, b4=2:
##    Score       df        p 
## 14293.39     2.00     0.00 
## 
## Test of b1=b3=-1:
##     Score        df         p 
## 0.6870218 2.0000000 0.7092758
```

# Robit

## Model

In robit regression, the latent variable $z_{i}$ follows a $t_{\nu}$ distribution with center $\eta_{i}$, unit scale, and specified degrees of freedom $\nu$:

$$
z_{i}|x_{i} \sim t_{\nu}(\eta_{i}) \\
y_{i} = I[z_{i}>0]
$$

## Estimation

The function `Fit.BinReg` with `model="robit"` specifies a robit model. Note that the degrees of freedom requires specification. The robit model with `df=7` is similar to the logistic model. 


```r
# Degrees of freedom
df = 7;
# Robit outcome
y = rBinReg(X,b,model="robit",df=df);
# Recover regression parameters
M = Fit.BinReg(y=y,X=X,model="robit",df=df);
show(M);
```

```
## Objective increment:  199 
## Objective increment:  52.1 
## Objective increment:  14.8 
## Objective increment:  1.8 
## Objective increment:  0.0339 
## Objective increment:  6.67e-06 
## Objective increment:  5.46e-10 
## 6 update(s) performed before tolerance limit.
## 
## Fitted robit Model with 7 degrees of freedom
## Estimated Coefficients:
##   Coeff   Point     SE       L      U        p
## 1    x1  1.1900 0.1050  0.9860  1.400 6.32e-30
## 2    x2 -1.2400 0.1100 -1.4500 -1.020 3.62e-29
## 3    x3  2.3500 0.1720  2.0100  2.690 1.40e-42
## 4    x4 -1.1900 0.1060 -1.4000 -0.986 2.72e-29
## 5    x5  0.0816 0.0776 -0.0705  0.234 2.93e-01
```

## Inference

The same hypotheses considered previously are tested for the robit model. Note that the degrees of freedom requires specification.


```r
cat("Test of b1=b2=b3=b4=0:\n");
Score.BinReg(y=y,X=X,L=c(F,T,T,T,T),model="robit",df=df,report=F);
cat("\n");
cat("Test of b4=0:\n");
Score.BinReg(y=y,X=X,L=c(F,F,F,F,T),model="robit",df=df,report=F);
cat("\n");
cat("Test of b2=0, b4=2:\n");
Score.BinReg(y=y,X=X,L=c(F,F,T,F,T),b10=c(0,2),model="robit",df=df,report=F);
cat("\n");
cat("Test of b1=b3=-1:\n");
Score.BinReg(y=y,X=X,L=c(F,T,F,T,F),b10=c(-1,-1),model="robit",df=df,report=F);
cat("\n");
```

```
## Test of b1=b2=b3=b4=0:
##         Score            df             p 
##  5.390587e+02  4.000000e+00 2.382912e-115 
## 
## Test of b4=0:
##     Score        df         p 
## 1.1064617 1.0000000 0.2928524 
## 
## Test of b2=0, b4=2:
##    Score       df        p 
## 3598.797    2.000    0.000 
## 
## Test of b1=b3=-1:
##      Score         df          p 
## 5.52201878 2.00000000 0.06322791
```
