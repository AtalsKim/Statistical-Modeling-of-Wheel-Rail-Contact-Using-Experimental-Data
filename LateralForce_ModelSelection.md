Model Selection, Mutliple Regression Model for Lateral Force
================
Mohammad Hosseini (<mohammadhosseini@vt.edu>)

-----

## Dataset Preparation & Single Models

``` r
## Load data
load("FinalDataset.rda")

## Divide into training and testing sets
set.seed(98)
n <- nrow(data)
size <- floor(n*0.7)
train.ind <- sample(n, size, replace=FALSE)
train <- data[train.ind,]
test <- data[-train.ind,]
```

<br>

## Multiple Regression

Before running the regression, we check for multicollinearity among the
independent variables.

``` r
## Correlation
c <- cor(data[, -c(1,2)])
print(c)
```

    ##                 Load           AoA     Creepage
    ## Load      1.00000000 -5.751417e-02 1.775960e-03
    ## AoA      -0.05751417  1.000000e+00 4.530105e-20
    ## Creepage  0.00177596  4.530105e-20 1.000000e+00

There is no strong correlation between the pairs of variables.

<br>

``` r
## Construct the model
library(sjPlot, quietly = TRUE)
multiple.reg <- lm(Lateral_Force ~ Load + AoA + I(AoA^3) + Creepage, data = train)
summary(multiple.reg)
```

    ## 
    ## Call:
    ## lm(formula = Lateral_Force ~ Load + AoA + I(AoA^3) + Creepage, 
    ##     data = train)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -2355.35  -577.66   -43.54   449.62  2728.84 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  3.419e+01  1.959e+00   17.45   <2e-16 ***
    ## Load         1.166e-01  2.439e-04  477.88   <2e-16 ***
    ## AoA          3.604e+03  3.428e+00 1051.37   <2e-16 ***
    ## I(AoA^3)    -2.221e+03  5.234e+00 -424.43   <2e-16 ***
    ## Creepage    -3.570e+01  1.370e+00  -26.05   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 917.4 on 1343995 degrees of freedom
    ## Multiple R-squared:  0.6949, Adjusted R-squared:  0.6949 
    ## F-statistic: 7.652e+05 on 4 and 1343995 DF,  p-value: < 2.2e-16

<br>

## Model Selection via BIC

``` r
## Step-wise selection
train2 <- train[,-1]
train2$AoA_Cubed <- train2$AoA^3
multiple.null <- lm(Lateral_Force ~ 1, data = train2)
multiple.full <- lm(Lateral_Force ~ . + .^2, data = train2)
multiple.fwbk <- step(multiple.null, scope=formula(multiple.full),
                      direction="both", k=log(nrow(train2)), trace=0)
summary(multiple.fwbk)
```

    ## 
    ## Call:
    ## lm(formula = Lateral_Force ~ AoA + Load + AoA_Cubed + Creepage + 
    ##     AoA:Load + Load:AoA_Cubed + AoA:AoA_Cubed + AoA:Creepage + 
    ##     AoA_Cubed:Creepage + Load:Creepage, data = train2)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1554.95  -143.53     3.86   142.27  1809.80 
    ## 
    ## Coefficients:
    ##                      Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        -1.162e+02  7.700e-01 -150.88   <2e-16 ***
    ## AoA                 6.783e+02  2.344e+00  289.33   <2e-16 ***
    ## Load                1.650e-01  1.348e-04 1223.63   <2e-16 ***
    ## AoA_Cubed          -5.113e+02  5.790e+00  -88.31   <2e-16 ***
    ## Creepage           -6.755e+01  6.946e-01  -97.25   <2e-16 ***
    ## AoA:Load            7.165e-01  2.916e-04 2457.06   <2e-16 ***
    ## Load:AoA_Cubed     -4.290e-01  4.458e-04 -962.36   <2e-16 ***
    ## AoA:AoA_Cubed       4.788e+02  9.151e-01  523.25   <2e-16 ***
    ## AoA:Creepage       -9.118e+02  2.160e+00 -422.19   <2e-16 ***
    ## AoA_Cubed:Creepage  9.395e+02  3.133e+00  299.88   <2e-16 ***
    ## Load:Creepage      -1.918e-02  1.161e-04 -165.12   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 250.7 on 1343989 degrees of freedom
    ## Multiple R-squared:  0.9772, Adjusted R-squared:  0.9772 
    ## F-statistic: 5.762e+06 on 10 and 1343989 DF,  p-value: < 2.2e-16

<br>

## Manual Selection

``` r
## Construct the model
library(sjPlot, quietly = TRUE)
multiple.reg <- lm(Lateral_Force ~ Load + AoA + AoA_Cubed + Creepage + Load*AoA, data = train2)
summary(multiple.reg)
```

    ## 
    ## Call:
    ## lm(formula = Lateral_Force ~ Load + AoA + AoA_Cubed + Creepage + 
    ##     Load * AoA, data = train2)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1885.07  -223.20   -13.68   197.59  2373.26 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error  t value Pr(>|t|)    
    ## (Intercept) -4.303e+01  7.676e-01   -56.05   <2e-16 ***
    ## Load         1.436e-01  9.601e-05  1495.16   <2e-16 ***
    ## AoA          1.446e+03  1.559e+00   927.73   <2e-16 ***
    ## AoA_Cubed   -2.209e+03  2.049e+00 -1077.90   <2e-16 ***
    ## Creepage    -3.751e+01  5.365e-01   -69.92   <2e-16 ***
    ## Load:AoA     4.567e-01  1.676e-04  2724.48   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 359.2 on 1343994 degrees of freedom
    ## Multiple R-squared:  0.9532, Adjusted R-squared:  0.9532 
    ## F-statistic: 5.477e+06 on 5 and 1343994 DF,  p-value: < 2.2e-16

``` r
confint(multiple.reg, level = 0.95)
```

    ##                     2.5 %        97.5 %
    ## (Intercept)   -44.5338896   -41.5248379
    ## Load            0.1433625     0.1437389
    ## AoA          1442.8787887  1448.9882574
    ## AoA_Cubed   -2213.0089503 -2204.9756393
    ## Creepage      -38.5629900   -36.4599390
    ## Load:AoA        0.4564097     0.4570668

<br>

## Prediction Error

### Model Selected via BIC

``` r
## Prediction Error
library(MLmetrics, quietly = TRUE)
test2 <- test[,-1]
test2$AoA_Cubed <- test2$AoA^3
multiple.fwbk.pred <- predict(multiple.fwbk, newdata=test2)
MSE <- MSE(multiple.fwbk.pred, test2$Lateral_Force)
RMSE <- sqrt(MSE)

pred.errors <- round(c(MSE, RMSE), 2)
names(pred.errors) <- c("MSE", "RMSE")
print(pred.errors)
```

    ##      MSE     RMSE 
    ## 62569.97   250.14

<br>

### Model Selected Manually

``` r
## Prediction Error
multiple.reg.pred <- predict(multiple.reg, newdata=test2)
MSE <- MSE(multiple.reg.pred, test2$Lateral_Force)
RMSE <- sqrt(MSE)

pred.errors <- round(c(MSE, RMSE), 2)
names(pred.errors) <- c("MSE", "RMSE")
print(pred.errors)
```

    ##       MSE      RMSE 
    ## 128530.50    358.51

-----
