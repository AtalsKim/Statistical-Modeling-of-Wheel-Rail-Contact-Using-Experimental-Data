Regression Analysis, Modeling Creepage
================
Mohammad Hosseini (<mohammadhosseini@vt.edu>)

  - [1. Creepage with a Quadratic
    Polynomial](#creepage-with-a-quadratic-polynomial)
  - [2. Creepage with Square-root](#creepage-with-square-root)
  - [3. Creepage with Squared and Square-root
    Terms](#creepage-with-squared-and-square-root-terms)
  - [4. Locally Linear Model](#locally-linear-model)
  - [5. Comparison across Models](#comparison-across-models)
  - [6. Transformations](#transformations)
  - [7. Conclusion](#conclusion)

-----

## 1\. Creepage with a Quadratic Polynomial

``` r
load("creepage_sweep_100.rda")
data <- creepage_sweep_100

## Training and testing data sets
n <- nrow(data)
size <- floor(n*0.7)
train.ind <- sample(n, size, replace=FALSE)
test.ind <- setdiff(n, train.ind)
train <- data[train.ind,]
test <- data[test.ind,]

# Run linear regression
train$creepage2 <- train$creepage^2
lm.creep2 <- lm(longitudinal.force ~ creepage + creepage2, data=train)
summary(lm.creep2)
```

    ## 
    ## Call:
    ## lm(formula = longitudinal.force ~ creepage + creepage2, data = train)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1272.23  -202.43    11.02   188.02  1149.26 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   317.31      11.89   26.68   <2e-16 ***
    ## creepage     6246.59      28.81  216.81   <2e-16 ***
    ## creepage2   -2090.30      14.30 -146.20   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 274.6 on 4197 degrees of freedom
    ## Multiple R-squared:  0.9622, Adjusted R-squared:  0.9622 
    ## F-statistic: 5.346e+04 on 2 and 4197 DF,  p-value: < 2.2e-16

<br>

### 1.1 Creepage with Polynomial Kernel, Assumptions Check

``` r
par(mfrow = c(2,2), mai=c(.7,.7,.2,.2))
plot(train$creepage, train$longitudinal.force, xlab="Creepage", ylab="Longitudinal Force",
     col=1, pch=20, main="")
xgrid <- seq(0,2,length=100)
ygrid <- lm.creep2$coef[1] + lm.creep2$coef[2]*xgrid + lm.creep2$coef[3]*xgrid^2
lines(xgrid, ygrid, col="seagreen3", lwd=2)

## Original data
long.trans <- (train$longitudinal.force)
creep.trans <- (train$creepage)
creep.trans2 <- creep.trans^2
lm.creep.tsfd <- lm(long.trans ~ creep.trans + creep.trans2)
plot(lm.creep.tsfd$fitted, rstudent(lm.creep.tsfd), xlab="Fitted Values", ylab="Studentized Residuals",
     col=1, main=""); abline(h=0, col=8, lty=2)
qqnorm(rstudent(lm.creep.tsfd), pch=20, col="seagreen3", main=""); abline(a=0, b=1, lty=2)
hist(rstudent(lm.creep.tsfd),freq=FALSE, col="seagreen3", xlab="Studentized Residuals", main="")
```

![](LongitudinalForce_Creepage_Single_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

## 2\. Creepage with Square-root

``` r
# Run linear regression
train$creepage.sqr <- sqrt(train$creepage)
lm.creep.sqr <- lm(longitudinal.force ~ creepage + creepage.sqr, na.rm=TRUE, data=train)
summary(lm.creep.sqr)
```

    ## 
    ## Call:
    ## lm(formula = longitudinal.force ~ creepage + creepage.sqr, data = train, 
    ##     na.rm = TRUE)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1145.72  -187.37    -8.41   165.04  1653.41 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  -1703.59      22.15  -76.93   <2e-16 ***
    ## creepage     -3004.93      33.36  -90.07   <2e-16 ***
    ## creepage.sqr  8990.79      56.62  158.80   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 256 on 4197 degrees of freedom
    ## Multiple R-squared:  0.9672, Adjusted R-squared:  0.9672 
    ## F-statistic: 6.182e+04 on 2 and 4197 DF,  p-value: < 2.2e-16

<br>

### 2.1 Creepage with Squared Term, Assumptions Check

``` r
par(mfrow = c(2,2), mai=c(.7,.7,.2,.2))
plot(train$creepage, train$longitudinal.force, xlab="Creepage", ylab="Longitudinal Force",
     col=1, pch=20, main="")
lines(xgrid, ygrid, col="seagreen3", lwd=2)
plot(lm.creep.sqr$fitted, rstudent(lm.creep.sqr), xlab="Fitted Values", ylab="Studentized Residuals",
     col=1, main=""); abline(h=0, col=8, lty=2)
qqnorm(rstudent(lm.creep.sqr), pch=20, col="seagreen3", main=""); abline(a=0, b=1, lty=2)
hist(rstudent(lm.creep.sqr),freq=FALSE, col="seagreen3", xlab="Studentized Residuals", main="")
```

![](LongitudinalForce_Creepage_Single_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

<br>

## 3\. Creepage with Squared and Square-root Terms

``` r
# Run linear regression
train$creepage.sqr <- sqrt(train$creepage)
lm.creep.sqr2 <- lm(longitudinal.force ~ creepage + creepage.sqr + creepage2, na.rm=TRUE, data=train)
summary(lm.creep.sqr2)
```

    ## 
    ## Call:
    ## lm(formula = longitudinal.force ~ creepage + creepage.sqr + creepage2, 
    ##     data = train, na.rm = TRUE)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1194.89  -152.27   -14.53   145.45  1067.11 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   -934.06      29.20 -31.987   <2e-16 ***
    ## creepage      1014.52     117.50   8.634   <2e-16 ***
    ## creepage.sqr  5281.21     116.19  45.455   <2e-16 ***
    ## creepage2     -966.37      27.36 -35.325   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 224.8 on 4196 degrees of freedom
    ## Multiple R-squared:  0.9747, Adjusted R-squared:  0.9747 
    ## F-statistic: 5.387e+04 on 3 and 4196 DF,  p-value: < 2.2e-16

<br>

### 3.1 Creepage with Squared and Square-root Terms, Assumptions Check

``` r
par(mfrow = c(2,2), mai=c(.7,.7,.2,.2))
plot(train$creepage, train$longitudinal.force, xlab="Creepage", ylab="Longitudinal Force",
     col=1, pch=20, main="")
lines(xgrid, ygrid, col="seagreen3", lwd=2)
plot(lm.creep.sqr2$fitted, rstudent(lm.creep.sqr2), xlab="Fitted Values", ylab="Studentized Residuals",
     col=1, main=""); abline(h=0, col=8, lty=2)
qqnorm(rstudent(lm.creep.sqr2), pch=20, col="seagreen3", main=""); abline(a=0, b=1, lty=2)
hist(rstudent(lm.creep.sqr2),freq=FALSE, col="seagreen3", xlab="Studentized Residuals", main="")
```

![](LongitudinalForce_Creepage_Single_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

<br>

## 4\. Locally Linear Model

### 4.1 Clustering

``` r
library(mclust, quietly = TRUE)
library(sjPlot, quietly = TRUE)
creepage.clusters <- Mclust(data=train[, c(5, 1)], G=2:5, modelNames = mclust.options("emModelNames"))
plot(creepage.clusters$BIC)
```

![](LongitudinalForce_Creepage_Single_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

<br>

There is not much difference in BIC and the slopes of the lines is very
small, so we take the simplest case of two clusters.

``` r
## 2 clusters
creepage.2clusters <- Mclust(data=train[, c(5, 1)], G=2, modelNames = mclust.options("emModelNames"))
mclust2Dplot(data=train[, c(5, 1)], parameters=creepage.2clusters$parameters, z=creepage.2clusters$z,
             xlab="Creepage", ylab="Longitudinal Force", col=c("seagreen3", "royalblue2"))
abline(v=0.6, col=8, lty=2, lwd=2); text(0.35, 5000, "0.6% Creepage", col=2)
```

![](LongitudinalForce_Creepage_Single_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
## Density plot
x <- seq(min(train$creepage), max(train$creepage), length=50) 
y <- seq(min(train$longitudinal.force), max(train$longitudinal.force), length=50) 
xy <- expand.grid(x,y)
xyD <- matrix(dens(modelName=creepage.2clusters$modelName, data=xy, 
                   parameters=creepage.2clusters$parameters), 
                   nrow=length(x), ncol=length(y))
persp(x, y, xyD, theta=-20, xlab="Creepage", ylab="Longitudinal Force", zlab="density")
```

![](LongitudinalForce_Creepage_Single_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

### 4.2 Locally Linear Regression Model

``` r
## Add categories to train
train$creepage.clusters <- as.factor(ifelse(train$creepage < 0.6, "C1", "C2"))

# Run linear regression
lm.creep.clust <- lm(longitudinal.force ~ creepage*creepage.clusters, data=train)
summary(lm.creep.clust)
```

    ## 
    ## Call:
    ## lm(formula = longitudinal.force ~ creepage * creepage.clusters, 
    ##     data = train)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1102.93  -128.82     9.69   142.16   910.77 
    ## 
    ## Coefficients:
    ##                              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)                    118.15      11.31   10.45   <2e-16 ***
    ## creepage                      6168.73      33.06  186.61   <2e-16 ***
    ## creepage.clustersC2           3288.73      17.42  188.78   <2e-16 ***
    ## creepage:creepage.clustersC2 -5305.35      34.58 -153.42   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 219.1 on 4196 degrees of freedom
    ## Multiple R-squared:  0.976,  Adjusted R-squared:  0.9759 
    ## F-statistic: 5.679e+04 on 3 and 4196 DF,  p-value: < 2.2e-16

<br>

### 4.3 Assumptions Check

``` r
par(mfrow = c(2,2), mai=c(.7,.7,.2,.2))
x1.grid <- with(subset(train, creepage.clusters=="C1"), seq(min(creepage), max(creepage), length=100))
x2.grid <- with(subset(train, creepage.clusters=="C2"), seq(min(creepage), max(creepage), length=100))
y1.grid <- with(subset(train, creepage.clusters=="C2"), 
                coef(lm.creep.clust)[1] + coef(lm.creep.clust)[2]*x1.grid)
y2.grid <- with(subset(train, creepage.clusters=="C2"), 
                (coef(lm.creep.clust)[1] + coef(lm.creep.clust)[3]) +
                (coef(lm.creep.clust)[2] + coef(lm.creep.clust)[4])*x2.grid)
plot(train$creepage, train$longitudinal.force, xlab="Creepage", ylab="Longitudinal Force",
     col=ifelse(train$creepage.clusters=="C1", "royalblue2", "seagreen3"), pch=20, main="")
lines(x1.grid, y1.grid, col="indianred2", lwd=2)
lines(x2.grid, y2.grid, col="indianred2", lwd=2)
plot(lm.creep.clust$fitted, rstudent(lm.creep.clust), xlab="Fitted Values", 
     ylab="tStudentized Residuals", col=1, main=""); abline(h=0, col=8, lty=2)
qqnorm(rstudent(lm.creep.clust), pch=20, col="seagreen3", main=""); abline(a=0, b=1, lty=2)
hist(rstudent(lm.creep.clust),freq=FALSE, col="seagreen3", xlab="Studentized Residuals", main="")
```

![](LongitudinalForce_Creepage_Single_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

<br>

## 5\. Comparison across Models

![](LongitudinalForce_Creepage_Single_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->![](LongitudinalForce_Creepage_Single_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

<br>

## 6\. Transformations

To fix the heteroscedasticity, we apply transformations to the
explanatory variables and the response. Log transformation is typically
used on such occasions and the results of applying it to our data set is
shown below.

![](LongitudinalForce_Creepage_Single_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->![](LongitudinalForce_Creepage_Single_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->![](LongitudinalForce_Creepage_Single_files/figure-gfm/unnamed-chunk-12-3.png)<!-- -->

We observe that taking the logarithm of the response, in our case, the
longitudinal force, resolves the non-constant variance; however, it
heavily skews the distribution of the residuals and makes the non-linear
pattern stronger.

<br>

## 7\. Conclusion

Comparing the log-transformation and the quadratic polynomial, the
latter appears to be more plausible since the normality assumption is
met and the non-linearity in the residuals is not as strong as the
log-transformed model. Here we don’t have a strong model satisfying all
the assumptions of the linear regression as it was the case with the
longitudinal and the AoA but we can consider the polynomial (quadratic)
regression model as a decent proxy.

-----
