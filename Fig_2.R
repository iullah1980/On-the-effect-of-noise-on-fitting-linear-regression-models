library(MASS)
library(reshape2)
library(tidyverse)

# Function to compute estimates and predictions for the following figures
ridgeless_reg_pred <- function(xTrain, xTest, yTrain){
  
  # Calculate regression coefficients using the generalized inverse
  b <- ginv(crossprod(xTrain))%*%crossprod(xTrain, yTrain)
  # The commented out line below adds a tiny value (lambda) to the diagonal of X'X for numerical stability
  # b <- solve(t(xTrain) %*% xTrain + .00000001 * diag(ncol(xTrain))) %*% t(xTrain) %*% yTrain
  
  # Compute the predictions
  predTest <- xTest%*%b
  
}

#################################################
# Figure 2 column (a) 
#################################################

# Simulation parameters
n0 <- 50 # Size of the training data
d0 <- 25 # Number of original variables
k <- 10 # Number of additional variables per each of d0 variables
N <- 10000 # Total number of observations

# Generate dataset
Sigma <- diag(1, d0) # Covariance matrix
DF <- as.data.frame(mvrnorm(N, mu = rep(0, d0), Sigma = Sigma))

# Define true parameters
TPar <- matrix(rnorm(d0, sd = 1))
TPar <- TPar / sqrt((t(TPar) %*% TPar)[1,1]) # Normalized coefficients

# Generate response variable
y <- as.matrix(DF) %*% TPar + rnorm(N, sd = .5)

# Prepare dataset
myData <- data.frame(y = y, DF)


n_cols_myData <- ncol(myData)-1


# Adding noise variables
for (cvar in 1:n_cols_myData){
  noise_var = matrix(rnorm(nrow(myData)*k), ncol = k)
  myData = data.frame(myData,noise_var)
}


# Prepare the evaluation sequence for the number of variables
d_evaluate <- seq(1, d0 * k + d0, 1)

# Prepare data for regression
X <- model.matrix(~ ., data = myData[, -1])
#y <- myData[, 1]

# True coefficients for the model
beta <- c(0, TPar[, 1], rep(0, ncol(myData) - d0 - 1))

# The last 5000 observations will be used for test data and n0 will be chosen from first 5000
xTest = X[5001:N, ]
Truth <- matrix(y[5001:N])


# Perform simulations to calculate bias and variance
res <- replicate(100, {
  train_index <- sample(1:5000, n0, replace = FALSE)
  xTrain <- X[train_index, ]
  yTrain <- y[train_index,]
  
  # Calculate bias for different numbers of variables
  pred <- sapply(d_evaluate, function(d) {
    ridgeless_reg_pred(xTrain[, 1:d , drop=F], xTest[, 1:d, drop=F], yTrain)
  })
  
  pred
  
})


# Process the results
res_df <- melt(res)
colnames(res_df)=c("obs", "d_eval", "rep", "pred")

# Calculate bias and variance from the results
bias_var = res_df %>% 
  group_by(obs, d_eval) %>% 
  summarize(
    bias = mean(pred, na.rm = TRUE),
    variance = sqrt(var(pred, na.rm = TRUE))
  ) %>%
  ungroup() %>%
  mutate(d_eval=rep(d_evaluate, N-5000)) %>%
  arrange(d_eval) %>%
  mutate(Truth=rep(Truth[,1],length(d_evaluate)),
         bias=sqrt((bias-Truth)^2)
         ) 

# Plotting bias
bias_g <- bias_var %>%
  ggplot(aes(x=d_eval, y=bias)) + 
  geom_line(aes(group=obs), color="red", alpha=.05, linewidth=.5) +
  stat_summary(fun=mean, geom="line", lwd=1.5, col="black") +
  coord_cartesian(ylim = c(0, 5)) +
  geom_vline(xintercept=n0, color="gray",lwd=.5) +
  geom_vline(xintercept=d0, color="black",lwd=.5) +
  labs(x = 'Number of explanatory variables (d)', y = 'Bias', title = "(a)") +
  theme(legend.position="none", 
        plot.caption = element_text(hjust = 0),
        axis.text = element_text(size = 12),
        axis.title = element_text(size=14))
bias_g


# Plotting variance
var_g <- bias_var %>%
  ggplot(aes(x=d_eval, y=variance)) + 
  geom_line(aes(group=obs), color="red", alpha=.05, linewidth=.5) +
  stat_summary(fun=mean, geom="line", lwd=1.5, col="black") +
  coord_cartesian(ylim = c(0, 10)) +
  geom_vline(xintercept=n0, color="gray",lwd=.5) +
  geom_vline(xintercept=d0, color="black",lwd=.5) +
  labs(x = 'Number of explanatory variables (d)', y = 'Variance', title = "") +
  theme(legend.position="none", 
        plot.caption = element_text(hjust = 0),
        axis.text = element_text(size = 12),
        axis.title = element_text(size=14))
var_g

# Save bias plots

png(file="Bias_d025n050.png")
par(mar=c(6, 6, 1, 1))
bias_g
dev.off()

png(file="Var_d025n050.png")
par(mar=c(6, 6, 1, 1))
var_g
dev.off()

#################################################
# Figure 2 column (b) 
#################################################

# Simulation parameters
n0 <- 50 # Size of the training data
d0 <- 50 # Number of original variables
k <- 10 # Number of additional variables per each of d0 variables
N <- 10000 # Total number of observations

# Generate dataset
Sigma <- diag(1, d0) # Covariance matrix
DF <- as.data.frame(mvrnorm(N, mu = rep(0, d0), Sigma = Sigma))

# Define true parameters
TPar <- matrix(rnorm(d0, sd = 1))
TPar <- TPar / sqrt((t(TPar) %*% TPar)[1,1]) # Normalized coefficients

# Generate response variable
y <- as.matrix(DF) %*% TPar + rnorm(N, sd = .5)

# Prepare dataset
myData <- data.frame(y = y, DF)


n_cols_myData <- ncol(myData)-1


# Adding noise variables
for (cvar in 1:n_cols_myData){
  noise_var = matrix(rnorm(nrow(myData)*k), ncol = k)
  myData = data.frame(myData,noise_var)
}


# Prepare the evaluation sequence for the number of variables
d_evaluate <- seq(1, d0 * k + d0, 1)

# Prepare data for regression
X <- model.matrix(~ ., data = myData[, -1])
#y <- myData[, 1]

# True coefficients for the model
beta <- c(0, TPar[, 1], rep(0, ncol(myData) - d0 - 1))

# The last 5000 observations will be used for test data and n0 will be chosen from first 5000
xTest = X[5001:N, ]
Truth <- matrix(y[5001:N])

# Perform simulations to calculate bias and variance
res <- replicate(100, {
  train_index <- sample(1:5000, n0, replace = FALSE)
  xTrain <- X[train_index, ]
  yTrain <- y[train_index,]
  
  # Calculate bias for different numbers of variables
  pred <- sapply(d_evaluate, function(d) {
    ridgeless_reg_pred(xTrain[, 1:d , drop=F], xTest[, 1:d, drop=F], yTrain)
  })
  
  pred
  
})

# Process the results
res_df <- melt(res)
colnames(res_df)=c("obs", "d_eval", "rep", "pred")

# Calculate bias and variance from the results
bias_var = res_df %>% 
  group_by(obs, d_eval) %>% 
  summarize(
    bias = mean(pred, na.rm = TRUE),
    variance = sqrt(var(pred, na.rm = TRUE))
  ) %>%
  ungroup() %>%
  mutate(d_eval=rep(d_evaluate,N-5000)) %>%
  arrange(d_eval) %>%
  mutate(Truth=rep(Truth[,1],length(d_evaluate)),
         bias=sqrt((bias-Truth)^2)
  ) 

# Plotting bias
bias_g <- bias_var %>%
  ggplot(aes(x=d_eval, y=bias)) + 
  geom_line(aes(group=obs), color="red", alpha=.05, linewidth=.5) +
  stat_summary(fun=mean, geom="line", lwd=1.5, col="black") +
  coord_cartesian(ylim = c(0, 5)) +
  geom_vline(xintercept=n0, color="gray",lwd=.5) +
  geom_vline(xintercept=d0, color="black",lwd=.5) +
  labs(x = 'Number of explanatory variables (d)', y = 'Bias', title = "(b)") +
  theme(legend.position="none", 
        plot.caption = element_text(hjust = 0),
        axis.text = element_text(size = 12),
        axis.title = element_text(size=14))
bias_g


# Plotting variance
var_g <- bias_var %>%
  ggplot(aes(x=d_eval, y=variance)) + 
  geom_line(aes(group=obs), color="red", alpha=.05, linewidth=.5) +
  stat_summary(fun=mean, geom="line", lwd=1.5, col="black") +
  coord_cartesian(ylim = c(0, 10)) +
  geom_vline(xintercept=n0, color="gray",lwd=.5) +
  geom_vline(xintercept=d0, color="black",lwd=.5) +
  labs(x = 'Number of explanatory variables (d)', y = 'Variance', title = "") +
  theme(legend.position="none", 
        plot.caption = element_text(hjust = 0),
        axis.text = element_text(size = 12),
        axis.title = element_text(size=14))
var_g

# Save bias plots

png(file="Bias_d050n050.png")
par(mar=c(6, 6, 1, 1))
bias_g
dev.off()

png(file="Var_d050n050.png")
par(mar=c(6, 6, 1, 1))
var_g
dev.off()

#################################################
# Figure 2 column (c) 
#################################################

# Simulation parameters
n0 <- 50 # Size of the training data
d0 <- 75 # Number of original variables
k <- 10 # Number of additional variables per each of d0 variables
N <- 10000 # Total number of observations

# Generate dataset
Sigma <- diag(1, d0) # Covariance matrix
DF <- as.data.frame(mvrnorm(N, mu = rep(0, d0), Sigma = Sigma))

# Define true parameters
TPar <- matrix(rnorm(d0, sd = 1))
TPar <- TPar / sqrt((t(TPar) %*% TPar)[1,1]) # Normalized coefficients

# Generate response variable
y <- as.matrix(DF) %*% TPar + rnorm(N, sd = .5)

# Prepare dataset
myData <- data.frame(y = y, DF)


n_cols_myData <- ncol(myData)-1


# Adding noise variables
for (cvar in 1:n_cols_myData){
  noise_var = matrix(rnorm(nrow(myData)*k), ncol = k)
  myData = data.frame(myData,noise_var)
}


# Prepare the evaluation sequence for the number of variables
d_evaluate <- seq(1, d0 * k + d0, 1)

# Prepare data for regression
X <- model.matrix(~ ., data = myData[, -1])
#y <- myData[, 1]

# True coefficients for the model
beta <- c(0, TPar[, 1], rep(0, ncol(myData) - d0 - 1))

# The last 5000 observations will be used for test data and n0 will be chosen from first 5000
xTest = X[5001:N, ]
Truth <- matrix(y[5001:N])


# Perform simulations to calculate bias and variance
res <- replicate(100, {
  train_index <- sample(1:5000, n0, replace = FALSE)
  xTrain <- X[train_index, ]
  yTrain <- y[train_index,]
  
  # Calculate bias for different numbers of variables
  pred <- sapply(d_evaluate, function(d) {
    ridgeless_reg_pred(xTrain[, 1:d , drop=F], xTest[, 1:d, drop=F], yTrain)
  })
  
  pred
  
})


# Process the results
res_df <- melt(res)
colnames(res_df)=c("obs", "d_eval", "rep", "pred")

# Calculate bias and variance from the results
bias_var = res_df %>% 
  group_by(obs, d_eval) %>% 
  summarize(
    bias = mean(pred, na.rm = TRUE),
    variance = sqrt(var(pred, na.rm = TRUE))
  ) %>%
  ungroup() %>%
  mutate(d_eval=rep(d_evaluate,N-5000)) %>%
  arrange(d_eval) %>%
  mutate(Truth=rep(Truth[,1],length(d_evaluate)),
         bias=sqrt((bias-Truth)^2)
  ) 

# Plotting bias
bias_g <- bias_var %>%
  ggplot(aes(x=d_eval, y=bias)) + 
  geom_line(aes(group=obs), color="red", alpha=.05, linewidth=.5) +
  stat_summary(fun=mean, geom="line", lwd=1.5, col="black") +
  coord_cartesian(ylim = c(0, 5)) +
  geom_vline(xintercept=n0, color="gray",lwd=.5) +
  geom_vline(xintercept=d0, color="black",lwd=.5) +
  labs(x = 'Number of explanatory variables (d)', y = 'Bias', title = "(c)") +
  theme(legend.position="none", 
        plot.caption = element_text(hjust = 0),
        axis.text = element_text(size = 12),
        axis.title = element_text(size=14))
bias_g


# Plotting variance
var_g <- bias_var %>%
  ggplot(aes(x=d_eval, y=variance)) + 
  geom_line(aes(group=obs), color="red", alpha=.05, linewidth=.5) +
  stat_summary(fun=mean, geom="line", lwd=1.5, col="black") +
  coord_cartesian(ylim = c(0, 10)) +
  geom_vline(xintercept=n0, color="gray",lwd=.5) +
  geom_vline(xintercept=d0, color="black",lwd=.5) +
  labs(x = 'Number of explanatory variables (d)', y = 'Variance', title = "") +
  theme(legend.position="none", 
        plot.caption = element_text(hjust = 0),
        axis.text = element_text(size = 12),
        axis.title = element_text(size=14))
var_g

# Save bias plots

png(file="Bias_d075n050.png")
par(mar=c(6, 6, 1, 1))
bias_g
dev.off()

png(file="Var_d075n050.png")
par(mar=c(6, 6, 1, 1))
var_g
dev.off()
