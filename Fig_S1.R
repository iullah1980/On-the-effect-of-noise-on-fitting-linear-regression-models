library(tidyverse)
library(MASS)


# Calculate OLS and minimum norm OLS estimates, and evaluate error
ridgeless_reg_error <- function(xTrain, xTest, yTrain, yTest){
  # Calculate regression coefficients using the Moore-Penrose pseudoinverse
  b <- ginv(crossprod(xTrain))%*%crossprod(xTrain, yTrain)
  # b <- solve(crossprod(xTrain)+.00000001*diag(ncol(xTrain)))%*%crossprod(xTrain, yTrain)
  predTrain <- xTrain%*%b
  predTest <- xTest%*%b
  
  lossTrain <- sqrt(mean((yTrain - predTrain[,1])^2))
  lossTest <- sqrt(mean((yTest - predTest[,1])^2))
  list(lossTrain=lossTrain, lossTest=lossTest)
}

# Aditional features are random noise
train_test_rmses <- function(n, k, myData) {
  

  # Randomly shuffle columns except for the first column (assumed to be the response variable)
  myData <- myData[, c(c(1), sample(2:ncol(myData)))]
  
  # Total number of original features excluding the response
  n_cols_myData <- ncol(myData)-1
  
  # Append k noise variables per each predictor to the data frame
  for (cvar in 1:n_cols_myData){
    if (dim(myData)[2] < d0+d0){
      noise_vars = matrix(rnorm(nrow(myData)*k, sd=.5), ncol = k)
      myData = data.frame(myData, noise_vars)
    } else {
      noise_vars = matrix(rnorm(nrow(myData)*k, sd=5), ncol = k)
      myData = data.frame(myData, noise_vars)
    }
  }
  
  # Determine the range of feature dimensions (d) for error evaluation.
  
  # Reduce computation time by batching the noise variables, 
  # rather than incorporating noise variables individually.
  # d_evaluate <- seq(1, n_cols_myData, 1)
  # for (sp in 1:n_cols_myData){
  #   d_evaluate <- c(d_evaluate, n_cols_myData+sp*k)
  # }
  
  d_evaluate <- c(1:ncol(myData))
  
  
  # Prepare the model matrix and response vector
  X <- model.matrix(~., data=myData[, -c(1)])
  y <- myData[, 1]
  
  # Split data into training and testing sets
  test_index <- sample(1:nrow(myData), n, replace = FALSE)
  xTrain = X[test_index, ]
  xTest = X[-test_index, ]
  yTrain = y[test_index]
  yTest = y[-test_index]
  
  # Calculate RMSEs for each dimension of features
  rmses <- map_df(d_evaluate, ~ridgeless_reg_error(xTrain[, 1:.x], xTest[, 1:.x], yTrain, yTest))
  
  # Combine the dimension values with their corresponding RMSEs
  tibble(d_evaluate, rmses)
  
}

##########################################
# Figure S1a
##########################################

# Define simulation parameters
d0 <- 30 # Number of original variables
n0 <- 50 # Size of training data
k <- 10 # Number of additional variables per each original variable
n_rpt <- 100 # Number of repetitions
N <- 10000 # Total sample size for generating data
Sigma <- diag(1, d0) # Covariance matrix for the variables



set.seed(1000)

# Generate dataset
DF <- as.data.frame(mvrnorm(N, mu=rep(0, d0), Sigma=Sigma))
coef <- matrix(rnorm(d0, sd=1)) # Coefficients for the linear model
# Normalize coefficients to unit length for interpretability
TPar <- coef# / sqrt((t(coef) %*% coef)[1,1])
y <- as.matrix(DF) %*% TPar + rnorm(N, sd=0.5) # Generate response variable
myData <- data.frame(y=y, DF) # Combine response and predictors

# Run simulations to calculate RMSEs for different numbers of features
results <- map_df(1:n_rpt, 
                  ~ train_test_rmses(n=n0, k=k, myData=myData), 
                  .id = "runs")



test_error_summary = results %>% 
  group_by(d_evaluate) %>% 
  summarize(
    rmse_train = mean(lossTrain, na.rm = TRUE),
    rmse_test = mean(lossTest, na.rm = TRUE)
  ) %>% 
  pivot_longer(-d_evaluate, names_to = 'lossType', values_to = 'RMSE') %>% 
  mutate(lossType = str_remove(lossType, 'loss')) 


g <- results %>%
  pivot_longer(-c(runs, d_evaluate), names_to = "Type", values_to = "RMSE") %>%
  mutate(Type = str_remove(Type, "loss")) %>% 
  ggplot(aes(x=d_evaluate, y=RMSE)) + 
  geom_line(aes(color=Type, group=interaction(Type, runs)), alpha=.2, linewidth=.5) +
  scale_colour_manual(values=c(2,4)) +
  stat_summary(fun=mean, geom="line", lwd=1.5, aes(group=Type, color=Type)) +
  geom_vline(xintercept=n0, color="gray",lwd=1) + # Training size, using n0
  geom_vline(xintercept=d0, color="black",lwd=1) + # Original var count
  geom_hline(yintercept=sqrt(mean(y^2)), color="black",lwd=1) + # Baseline error
  coord_cartesian(ylim = c(0, 15)) +
  labs(x = 'Number of explanatory variables (d)', y = 'Error', title = "(a)") +
  guides(color = guide_legend(override.aes = list(linewidth = 3, alpha=1)))+
  geom_point(color = 'black', alpha = .5, size = 3,
             data = test_error_summary %>% filter(lossType == 'rmse_test') %>% 
               filter(RMSE == min(RMSE)),
             show.legend = FALSE
  ) +
  theme(legend.position = c(0.9, 0.85), 
        plot.caption = element_text(hjust = 0),
        axis.text = element_text(size = 12),
        axis.title = element_text(size=14))



g # Display the plot


png(file="TripleDd030n050.png")
par(mar=c(6, 6, 1, 1))
g
dev.off()



##########################################
# Figure S1b
##########################################


# Aditional features are random noise
train_test_rmses <- function(n, k, myData=myData){
  
  # Randomly shuffle columns except for the first column (assumed to be the response variable)
  myData <- myData[, c(c(1), sample(2:ncol(myData)))]
  
  # Total number of original features excluding the response
  n_cols_myData <- ncol(myData)-1
  
  # Append k noise variables per each predictor to the data frame
  for (cvar in 1:n_cols_myData){
    if (dim(myData)[2] < d0+d0){
      noise_vars = matrix(rnorm(nrow(myData)*k, sd=.5), ncol = k)
      myData = data.frame(myData, noise_vars)
    } else {
      noise_vars = matrix(rnorm(nrow(myData)*k, sd=5), ncol = k)
      myData = data.frame(myData, noise_vars)
    }
  }
  
  # Determine the range of feature dimensions (d) for error evaluation.
  
  # Reduce computation time by batching the noise variables, 
  # rather than incorporating noise variables individually.
  # d_evaluate <- seq(1, n_cols_myData, 1)
  # for (sp in 1:n_cols_myData){
  #   d_evaluate <- c(d_evaluate, n_cols_myData+sp*k)
  # }
  
  d_evaluate <- c(1:ncol(myData))
  
  # Prepare the response vector
  y <- myData[, 1]
  test_index <- sample(1:nrow(myData), n, replace = FALSE)
  yTrain = y[test_index]
  yTest = y[-test_index]
  
  # Scale the data and prepare the model matrix
  myData <- data.frame(scale(myData[, -c(1)]))
  X <- model.matrix(~., data=myData)
  
  xTrain = X[test_index, ]
  xTest = X[-test_index, ]

  # Calculate RMSEs for each dimension of features
  rmses <- map_df(d_evaluate, ~ridgeless_reg_error(xTrain[, 1:.x], xTest[, 1:.x], yTrain, yTest))
  
  # Combine the dimension values with their corresponding RMSEs
  tibble(d_evaluate, rmses)
}

# Define simulation parameters
d0 <- 30 # Number of original variables
n0 <- 50 # Size of training data
k <- 10 # Number of additional variables per each original variable
n_rpt <- 100 # Number of repetitions
N <- 10000 # Total sample size for generating data
Sigma <- diag(1, d0) # Covariance matrix for the variables



set.seed(1000)

# Generate dataset
DF <- as.data.frame(mvrnorm(N, mu=rep(0, d0), Sigma=Sigma))
coef <- matrix(rnorm(d0, sd=1)) # Coefficients for the linear model
# Normalize coefficients to unit length for interpretability
TPar <- coef# / sqrt((t(coef) %*% coef)[1,1])
y <- as.matrix(DF) %*% TPar + rnorm(N, sd=0.5) # Generate response variable
myData <- data.frame(y=y, DF) # Combine response and predictors

# Run simulations to calculate RMSEs for different numbers of features
results <- map_df(1:n_rpt, 
                  ~ train_test_rmses(n=n0, k=k, myData=myData), 
                  .id = "runs")



test_error_summary = results %>% 
  group_by(d_evaluate) %>% 
  summarize(
    rmse_train = mean(lossTrain, na.rm = TRUE),
    rmse_test = mean(lossTest, na.rm = TRUE)
  ) %>% 
  pivot_longer(-d_evaluate, names_to = 'lossType', values_to = 'RMSE') %>% 
  mutate(lossType = str_remove(lossType, 'loss')) 


g <- results %>%
  pivot_longer(-c(runs, d_evaluate), names_to = "Type", values_to = "RMSE") %>%
  mutate(Type = str_remove(Type, "loss")) %>% 
  ggplot(aes(x=d_evaluate, y=RMSE)) + 
  geom_line(aes(color=Type, group=interaction(Type, runs)), alpha=.2, linewidth=.5) +
  scale_colour_manual(values=c(2,4)) +
  stat_summary(fun=mean, geom="line", lwd=1.5, aes(group=Type, color=Type)) +
  geom_vline(xintercept=n0, color="gray",lwd=1) + # Training size, using n0
  geom_vline(xintercept=d0, color="black",lwd=1) + # Original var count
  geom_hline(yintercept=sqrt(mean(y^2)), color="black",lwd=1) + # Baseline error
  coord_cartesian(ylim = c(0, 15)) +
  labs(x = 'Number of explanatory variables (d)', y = 'Error', title = "(b)") +
  guides(color = guide_legend(override.aes = list(linewidth = 3, alpha=1)))+
  geom_point(color = 'black', alpha = .5, size = 3,
             data = test_error_summary %>% filter(lossType == 'rmse_test') %>% 
               filter(RMSE == min(RMSE)),
             show.legend = FALSE
  ) +
  theme(legend.position = c(0.9, 0.85), 
        plot.caption = element_text(hjust = 0),
        axis.text = element_text(size = 12),
        axis.title = element_text(size=14))



g # Display the plot


png(file="DoubleDd030n050.png")
par(mar=c(6, 6, 1, 1))
g
dev.off()
