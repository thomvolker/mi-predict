# Small simulation script to evaluate coverage of prediction intervals after
# multiple imputation.


# Do for loop in parallel
library(foreach)

# Set seed
set.seed(123)

# Set number of simulations
nsim <- 200

# Create a function that properly calculates prediction intervals
predfunc <- \(x, fit, alpha = 0.05) {
  tval <- qt(1 - alpha/2, fit$df.residual)
  summfit <- summary(fit)
  sigma <- summfit$sigma
  xmat <- model.matrix(~x)
  se <- sqrt(1 + rowSums((xmat %*% summfit$cov.unscaled) * xmat))
  xfit <- xmat %*% fit$coefficients
  return(list(fit = cbind(xfit, xfit, xfit) + (tval * se * sigma) %x% t(c(0, -1, 1)),
              var = se^2 * sigma^2))
}

# Set the sample size
N <- 30
P <- 5
V <- 0.5 + diag(0.5, P)

# Register parallel back-end
cl <- parallel::makeCluster(parallel::detectCores() - 4)
doParallel::registerDoParallel(cl)

# Run for loop in parallel and store output in x
x <- foreach (i = 1:nsim, .combine = 'c') %dopar% {

  # Create training data set
  Xtrain <- matrix(rnorm(N * P), N, P) %*% chol(V)
  b <- sqrt(0.5 / c(1,2,3,4,5) %*% V %*% c(1,2,3,4,5))
  ytrain <- Xtrain %*% (c(1,2,3,4,5)*c(b)) + rnorm(N, 0, sqrt(1-0.5))

  # Fit model
  # train <- data.frame(y = ytrain, Xtrain[,-5])
  # fit <- lm(y ~ ., data = train)

  #preds <- predict(fit, interval = "pred", se = TRUE)

  #preds$fit |> head()

  #predfunc(Xtrain[,-5], fit)$fit |> head()
  #preds$fit |> head()

  # Create large, fully observed test set
  Ntest <- 10000
  Xtest <- matrix(rnorm(Ntest * P), Ntest, P) %*% chol(V)
  ytest <- Xtest %*% (c(1,2,3,4,5)*c(b)) + rnorm(Ntest, 0, sqrt(1-0.5))


  #predint <- predfunc(Xtest[,-5], fit)
  #mean(ytest > predint$fit[,2] & ytest < predint$fit[,3])

  # Create training data with missingness
  Xamp <- mice::ampute(train, prop = 0.5)


  # Set number of imputations and impute missing data with two approaches
  m <- 10
  Ximp1 <- mice::mice(Xamp$amp, m = 1, method = "norm.predict", print = FALSE)
  Ximpm <- mice::mice(Xamp$amp, m = m, method = "norm", print = FALSE)

  # Fit models on imputed data
  fitimp1 <- with(Ximp1, lm(y ~ X1 + X2 + X3 + X4))$analyses
  fitimpm <- with(Ximpm, lm(y ~ X1 + X2 + X3 + X4))$analyses

  # Calculate prediction intervals on imputed sets
  impvals1 <- purrr::map(fitimp1, \(x) {
    predint <- predfunc(Xtest[,-5], x)
    data.frame(
      q = predint$fit[,1],
      u = predint$var
    )
  }) |>
    unlist() |>
    array(dim = c(Ntest, 2, m)) |>
    apply(1, \(x) {
      pooled <- mice::pool.scalar(x[1,], x[2,])
      if (is.na(pooled$b)) {
        pooled$t <- pooled$ubar
        pooled$df <- fit$df.residual
      }
      pooled$qbar + c(-1, 1) * sqrt(pooled$t) * qt(0.975, pooled$df)
    }) |>
    t()

  impvals1m <- purrr::map(fitimpm, \(x) {
    predint <- predfunc(Xtest[,-5], x)
    data.frame(
      q = predint$fit[,1],
      u = predint$var
    )
  }) |>
    unlist() |>
    array(dim = c(Ntest, 2, m)) |>
    apply(1, \(x) {
      pooled <- mice::pool.scalar(x[1,], x[2,])
      if (is.na(pooled$b)) {
        pooled$t <- pooled$ubar
        pooled$df <- fit$df.residual
      }
      pooled$qbar + c(-1, 1) * sqrt(pooled$t) * qt(0.975, pooled$df)
    }) |>
    t()

  # collect output in vector
  c(imp1 = mean(ytest > impvals1[,1] & ytest < impvals1[,2]),
    impm = mean(ytest > impvals1m[,1] & ytest < impvals1m[,2]))
}

# Single imputation
mean(x[2*0:(nsim-1)+1])
# Multiple imputation
mean(x[2*1:nsim])
