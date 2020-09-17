SLR Analysis for Modeling AoA (sweep-dataset)
================
Mohammad Hosseini (<mohammadhosseini@vt.edu>)

  - [1. Regressing Longitudinal Force on
    AoA](#regressing-longitudinal-force-on-aoa)
  - [2. AoA with a Polynomial Kernel](#aoa-with-a-polynomial-kernel)
  - [3. AoA, Only Quadratic Term](#aoa-only-quadratic-term)
  - [4. Comparing Models via ANOVA](#comparing-models-via-anova)

-----

### 1\. Regressing Longitudinal Force on AoA

``` r
load("longForce_AoA_100.rda")
data <- longForce_AoA_100

## Training and testing data sets
n <- nrow(data)
size <- floor(n*0.7)
train.ind <- sample(n, size, replace=FALSE)
train <- data[train.ind,]
test <- data[-train.ind,]

# Running linear regressions
lm.AoA <- lm(train$longitudinal.force ~ train$AoA)

# Plotting the regression lines
plot(train$AoA, train$longitudinal.force, xlab="AoA",
     ylab="Longitudinal Force", col=1, pch=20, main=""); abline(lm.AoA, col="indianred1", lwd=2)
```

![](SLR_AoA_gitdoc_files/figure-gfm/SLRs.long-1.png)<!-- -->

``` r
summary(lm.AoA)
```

    ## 
    ## Call:
    ## lm(formula = train$longitudinal.force ~ train$AoA)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1090.02  -295.10     3.58   291.29  1114.32 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 4165.638      6.801 612.487  < 2e-16 ***
    ## train$AoA     71.050     11.739   6.053 1.58e-09 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 394.2 on 3358 degrees of freedom
    ## Multiple R-squared:  0.01079,    Adjusted R-squared:  0.0105 
    ## F-statistic: 36.64 on 1 and 3358 DF,  p-value: 1.581e-09

<br>

#### 1.1 Checking SLR Assumptions

``` r
par(mfrow = c(2, 2), mai=c(.7,.7,.2,.2))
plot(train$AoA, train$longitudinal.force, xlab="AoA",
     ylab="Longitudinal Force", col=1, pch=20, main=""); abline(lm.AoA, col="indianred1", lwd=2)
## AoA
plot(lm.AoA$fitted, rstudent(lm.AoA),xlab="Fitted values", ylab="Studentized Residuals",
     col=1, main=""); abline(h=0, col=8, lty=2)
qqnorm(rstudent(lm.AoA), pch=20, col="indianred1", main=""); abline(a=0, b=1, lty=2)
hist(rstudent(lm.AoA),freq=FALSE, col="indianred1", xlab="Studentized Residuals", main="")
```

![](SLR_AoA_gitdoc_files/figure-gfm/assumptions.check-1.png)<!-- -->

<br>

### 2\. AoA with a Polynomial Kernel

``` r
# Running linear regression
train$AoA2 <- train$AoA^2
lm.AoA2 <- lm(longitudinal.force ~ AoA + AoA2, data=train)

# Plotting the regression lines
plot(train$AoA, train$longitudinal.force, xlab="AoA", ylab="Longitudinal Force",
     col=1, pch=20, main="")
xgrid <- seq(-1,1,length=100)
ygrid <- lm.AoA2$coef[1] + lm.AoA2$coef[2]*xgrid + lm.AoA2$coef[3]*xgrid^2
lines(xgrid, ygrid, col="indianred1", lwd=2)

library(sjPlot, quietly = TRUE)
```

![](SLR_AoA_gitdoc_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
summary(lm.AoA2)
```

    ## 
    ## Call:
    ## lm(formula = longitudinal.force ~ AoA + AoA2, data = train)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -980.29 -148.86   -5.61  135.13  844.93 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  4535.513      5.661  801.15   <2e-16 ***
    ## AoA            65.959      6.488   10.17   <2e-16 ***
    ## AoA2        -1101.844     12.610  -87.38   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 217.9 on 3357 degrees of freedom
    ## Multiple R-squared:  0.6979, Adjusted R-squared:  0.6977 
    ## F-statistic:  3877 on 2 and 3357 DF,  p-value: < 2.2e-16

<br>

#### 2.1 Checking SLR Assumptions

``` r
par(mfrow = c(2, 2), mai=c(.7,.7,.2,.2))
plot(train$AoA, train$longitudinal.force, xlab="AoA", ylab="Longitudinal Force",
     col=1, pch=20, main="Polynomial Kernel")
xgrid <- seq(-1,1,length=100)
ygrid <- lm.AoA2$coef[1] + lm.AoA2$coef[2]*xgrid + lm.AoA2$coef[3]*xgrid^2
lines(xgrid, ygrid, col="indianred1", lwd=2)
## AoA
plot(lm.AoA2$fitted, rstudent(lm.AoA2),xlab="Fitted values", ylab="Studentized Residuals",
     col=1, main=""); abline(h=0, col=8, lty=2)
qqnorm(rstudent(lm.AoA2), pch=20, col="indianred1", main=""); abline(a=0, b=1, lty=2)
hist(rstudent(lm.AoA2),freq=FALSE, col="indianred1", xlab="Studentized Residuals", main="")
```

![](SLR_AoA_gitdoc_files/figure-gfm/assumptions.check2-1.png)<!-- -->
<br>

### 3\. AoA, Only Quadratic Term

``` r
# Running linear regression
lm.AoA3 <- lm(longitudinal.force ~ AoA2, data=train)

# Plotting the regression lines
plot(train$AoA, train$longitudinal.force, xlab="AoA", ylab="Longitudinal Force",
     col=1, pch=20, main="")
xgrid <- seq(-1,1,length=100)
ygrid <- lm.AoA3$coef[1] + lm.AoA3$coef[2]*xgrid^2
lines(xgrid, ygrid, col="indianred1", lwd=2)
```

![](SLR_AoA_gitdoc_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
summary(lm.AoA3)
```

    ## 
    ## Call:
    ## lm(formula = longitudinal.force ~ AoA2, data = train)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -948.06 -152.44   -5.11  142.69  782.25 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  4535.826      5.747  789.28   <2e-16 ***
    ## AoA2        -1102.995     12.800  -86.17   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 221.2 on 3358 degrees of freedom
    ## Multiple R-squared:  0.6886, Adjusted R-squared:  0.6885 
    ## F-statistic:  7425 on 1 and 3358 DF,  p-value: < 2.2e-16

<br>

#### 3.1 Checking SLR Assumptions

``` r
par(mfrow = c(2, 2), mai=c(.7,.7,.2,.2))
plot(train$AoA, train$longitudinal.force, xlab="AoA", ylab="Longitudinal Force",
     col=1, pch=20, main="AoA-Squared")
xgrid <- seq(-1,1,length=100)
ygrid <- lm.AoA3$coef[1] + lm.AoA3$coef[2]*xgrid^2
lines(xgrid, ygrid, col="indianred1", lwd=2)
## AoA
plot(lm.AoA3$fitted, rstudent(lm.AoA3),xlab="Fitted values", ylab="Studentized Residuals",
     col=1, main=""); abline(h=0, col=8, lty=2)
qqnorm(rstudent(lm.AoA3), pch=20, col="indianred1", main=""); abline(a=0, b=1, lty=2)
hist(rstudent(lm.AoA3),freq=FALSE, col="indianred1", xlab="Studentized Residuals", main="")
```

![](SLR_AoA_gitdoc_files/figure-gfm/assumptions.check3-1.png)<!-- -->

<br>

### 4\. Comparing Models via ANOVA

``` r
anova(lm.AoA2,lm.AoA3)
```

    ## Analysis of Variance Table
    ## 
    ## Model 1: longitudinal.force ~ AoA + AoA2
    ## Model 2: longitudinal.force ~ AoA2
    ##   Res.Df       RSS Df Sum of Sq      F    Pr(>F)    
    ## 1   3357 159392862                                  
    ## 2   3358 164299632 -1  -4906770 103.34 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

-----
