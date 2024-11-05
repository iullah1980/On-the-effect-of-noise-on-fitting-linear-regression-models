
# Load necessary libraries
library(MASS)      # For generating multivariate normal samples
library(Matrix)    # For matrix operations
library(reshape2)  # For reshaping data
library(tidyverse) # For data manipulation and plotting


#################################################
# Figure 6 
#################################################

# Initialize parameters
d0 <- 30          # Dimensionality of the feature space
N <- 10000        # Number of observations

# Covariance matrix
# A <- matrix(.5, p, p)
# Sigma <- A^ abs(row(A)-col(A))

Sigma <- matrix(0, d0, d0); diag(Sigma) = 1 

# Generate data
DF <- as.data.frame(mvrnorm(N, mu = rep(0, d0), Sigma = Sigma))

# Generate true parameters and response variable
TPar <- rnorm(d0, sd = 1) # True parameters
y <- as.matrix(DF) %*% TPar + rnorm(N, sd = .5) # Generate response variable with noise


# Prepare the dataset
myData <- data.frame(y = y, DF)

# Create model matrix excluding the response variable
X <- model.matrix(~., data = myData[,-c(1)])
y <- myData[, 1]

# Define grid for regularization parameter, number of samples for training, and number of noise dimensions
reg_par <- seq(-1, 1, length.out = 100)
n0_evaluate <- seq(10, 200, by = 10)
k_evaluate <- seq(0, 200, by = 10)


# Loop through each sample size and noise dimension to evaluate the ridge regression
# It will take more than 20 hours 
no_loop <- sapply(n0_evaluate, function(n0) {
  k_loop <- sapply(k_evaluate, function(k) {
    lambda_loop <- sapply(reg_par, function(lambda) {
      
      # Replicate experiment to average over random selections and noise
      ridge_grid <- replicate(500, {
        noise = matrix(rnorm((ncol(myData) - 1) * k, sd = 1), nrow = k)
        
        # Randomly select training data and optionally add noise
        test_index <- sample(1:nrow(myData), n0, replace = FALSE)
        if (k > 0) {
          xTrain = rbind(X[test_index, -1], noise)
        } else {
          xTrain = X[test_index, -1]
        }
        xTest = X[-test_index, -1]
        yTrain = c(y[test_index], rep(0, k))
        yTest = y[-test_index]
        
        # Perform ridge regression
        b <- solve(t(xTrain) %*% xTrain + (n0 + k) * lambda * diag(d0)) %*% t(xTrain) %*% yTrain
        predTest <- xTest %*% b
        # Calculate RMSE for test data
        lossTest <- yardstick::rmse_vec(yTest, predTest[, 1])
        lossTest
      })
      ridge_grid
    })
    # Cap outliers in the test error
    lambda_loop[lambda_loop > quantile(lambda_loop, probs = .95)] = 
      quantile(lambda_loop, probs = .95)
    # Return the best lambda for this setting of n0 and k
    reg_par[which.min(apply(lambda_loop, 2, mean))]
  })
  k_loop
})


# save(no_loop, file="all_testerr.RData")
# load(file="all_testerr.RData")

# Set column and row names for the results matrix
colnames(no_loop) <- n0_evaluate
rownames(no_loop) <- k_evaluate

# Reshape the results for plotting
results <- melt(no_loop)
colnames(results) <- c("k", "n0", "value")

# Plot the results using ggplot2
g <- ggplot(results, aes(x = n0, y = k, fill = value)) +
  geom_tile(show.legend = TRUE) + #color = "white"
  scale_fill_gradient2(low = "red", high = "blue", mid = "white", 
                       midpoint = 0, limits = c(-.5,.5)) + 
  labs(fill=expression(lambda)) +
  theme_light() +
  xlab(expression(n[0])) + ylab(expression(n-n[0])) +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(5, 205), expand=c(0,0))+
  scale_y_continuous(limits = c(-5, 205), expand=c(0,0))

g

png(file="NegLambdaIdentity.png")
par(mar=c(6,6,1,1))
g
dev.off()


