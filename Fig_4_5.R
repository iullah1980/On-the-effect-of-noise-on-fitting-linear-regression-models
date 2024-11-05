library(MASS) # For generating multivariate normal samples
library(tidyverse) # For data manipulation and visualization

###########################################
#  Functions used to generate data and perform computations for the following figures.
###########################################

# Define a function for ridgeless regression, calculating RMSE for training and test sets
ridgeless_reg <- function(X_train, X_test, y_train, y_test){
  
  # Calculate regression coefficients using the generalized inverse (Moore-Penrose inverse)
  b <- ginv(crossprod(X_train)) %*% crossprod(X_train, y_train)
  # An alternative using a tiny regularization factor to ensure numerical stability
  # b <- solve(crossprod(X_train) + .00000001 * diag(p)) %*% crossprod(X_train, y_train)
  
  # Predict responses for training and test datasets
  predTrain <- X_train %*% b
  predTest <- X_test %*% b
  
  # Calculate Root Mean Squared Error (RMSE) for both training and test predictions
  RMSE_train <- sqrt(mean((y_train - predTrain)^2))
  RMSE_test <- sqrt(mean((y_test - predTest)^2))
  
  # Return RMSEs as a list
  list(RMSE_Train=RMSE_train, RMSE_Test=RMSE_test)
}

# Function to evaluate model performance under augmentation
res <- function(N, n, k){
  # Randomly select training rows
  train_rows <- sample(1:N, n)
  
  # Split the data into training and testing sets
  X_train <- DF[train_rows, ]
  y_train <- y[train_rows, , drop = FALSE]
  
  # Augment X_train with additional observations 
  X_aug <- rbind(X_train, mvrnorm(k, mu = rep(0, p), Sigma = diag(p)))
  y_aug <- rbind(y_train, matrix(sample(0, k, replace = TRUE)))
  
  # Define test set as the data not included in the training set
  X_test <- DF[-train_rows, ]
  y_test <- y[-train_rows, , drop = FALSE]
  
  # Calculate RMSEs for a range of augmentation sizes
  rmses <- map_df(n_evaluate,
                  ~ridgeless_reg(X_aug[1:.x, ],
                                 X_test,
                                 y_aug[1:.x, ],
                                 y_test)
  )
  
  # Return a tibble with n_evaluate and corresponding RMSEs
  tibble(n_evaluate, rmses)
}

####################################################
# Figure 4a
####################################################
# Define simulation parameters
n0 <- 5 # Size of the training data
d0 <- 15 # Number of predictors
N <- 10000 # Total number of observations
Sigma <- diag(1, d0) # Create a covariance matrix

# Generate data matrix DF and response vector y based on true parameters TPar
DF <- mvrnorm(N, mu = rep(0, d0), Sigma = Sigma)

# Generate true regression coefficients
coef <- matrix(rnorm(d0, sd = 1))
TPar <- coef / sqrt((t(coef) %*% coef)[1,1])

# Generate response variable with added noise
y <- as.matrix(DF) %*% TPar + rnorm(N, sd = .5)

# Define a sequence of training sizes to evaluate
n_evaluate <- seq(2, 100, 1)

n_rpt <- 100 # Number of repetitions for the simulation

# Simulate and calculate RMSEs for different training sizes
results <- map_df(1:n_rpt, 
                  ~ res(N = N,
                        n = n0,
                        k = 100), # k represents number of augmented data points
                  .id = "runs")



# Summarize the test errors across different training sizes
test_error_summary = results %>% 
  group_by(n_evaluate) %>% 
  summarize(
    rmse_train = mean(RMSE_Train, na.rm = TRUE),
    rmse_test = mean(RMSE_Test, na.rm = TRUE)
  ) %>% 
  pivot_longer(-n_evaluate, names_to = 'Type', values_to = 'RMSE') %>% 
  mutate(lossType = str_remove(Type, 'RMSE_')) 

# Plot RMSEs for training and test sets across varying training sizes
g <- results %>%
  pivot_longer(-c(runs, n_evaluate), names_to = "Type", values_to = "RMSE") %>%
  mutate(Type = str_remove(Type, "RMSE_")) %>% 
  ggplot(aes(x=n_evaluate, y=RMSE)) + 
  geom_line(aes(color=Type, group=interaction(Type, runs)), alpha=.2, linewidth=.5) +
  scale_colour_manual(values=c(2,4)) +
  stat_summary(fun=mean, geom="line", lwd=1.5, aes(group=Type, color=Type)) +
  geom_vline(xintercept=n0, color="gray",lwd=1) +
  geom_vline(xintercept=d0, color="black",lwd=1) +
  geom_hline(yintercept=sqrt(mean(y^2)), color="black",lwd=1) +
  coord_cartesian(ylim = c(0, 8)) +#max(test_error_summary$RMSE) + .5)) +
  labs(x = 'Number of training samples (n)', y = 'Error', title = "(a)") +
  guides(color = guide_legend(override.aes = list(linewidth = 3, alpha=1)))+
  geom_point(
    color = 'black',
    alpha = .5,
    size = 3,
    data = test_error_summary %>% filter(lossType == 'rmse_test') %>% 
      filter(RMSE == min(RMSE)),
    show.legend = FALSE
  ) +
  theme(legend.position = c(0.9, 0.85), 
        plot.caption = element_text(hjust = 0),
        axis.text=element_text(size=12), 
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=20))
# Display the plot
g

# Save the plot
png(file="Fig_d015n05.png")
par(mar=c(6, 6, 1, 1))
g
dev.off()

####################################################
# Figure 4b
####################################################

n0 <- 15 # Size of the training data

# Simulate and calculate RMSEs for different training sizes
results <- map_df(1:n_rpt, 
                  ~ res(N = N,
                        n = n0,
                        k = 100), # k represents number of augmented data points
                  .id = "runs")



# Summarize the test errors across different training sizes
test_error_summary = results %>% 
  group_by(n_evaluate) %>% 
  summarize(
    rmse_train = mean(RMSE_Train, na.rm = TRUE),
    rmse_test = mean(RMSE_Test, na.rm = TRUE)
  ) %>% 
  pivot_longer(-n_evaluate, names_to = 'Type', values_to = 'RMSE') %>% 
  mutate(lossType = str_remove(Type, 'RMSE_')) 

# Plot RMSEs for training and test sets across varying training sizes
g <- results %>%
  pivot_longer(-c(runs, n_evaluate), names_to = "Type", values_to = "RMSE") %>%
  mutate(Type = str_remove(Type, "RMSE_")) %>% 
  ggplot(aes(x=n_evaluate, y=RMSE)) + 
  geom_line(aes(color=Type, group=interaction(Type, runs)), alpha=.2, linewidth=.5) +
  scale_colour_manual(values=c(2,4)) +
  stat_summary(fun=mean, geom="line", lwd=1.5, aes(group=Type, color=Type)) +
  geom_vline(xintercept=n0, color="gray",lwd=1) +
  geom_vline(xintercept=d0, color="black",lwd=1) +
  geom_hline(yintercept=sqrt(mean(y^2)), color="black",lwd=1) +
  coord_cartesian(ylim = c(0, 8)) +#max(test_error_summary$RMSE) + .5)) +
  labs(x = 'Number of training samples (n)', y = 'Error', title = "(b)") +
  guides(color = guide_legend(override.aes = list(linewidth = 3, alpha=1)))+
  geom_point(
    color = 'black',
    alpha = .5,
    size = 3,
    data = test_error_summary %>% filter(lossType == 'rmse_test') %>% 
      filter(RMSE == min(RMSE)),
    show.legend = FALSE
  ) +
  theme(legend.position = c(0.9, 0.85), 
        plot.caption = element_text(hjust = 0),
        axis.text=element_text(size=12), 
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=20))

# Display the plot
g

png(file="Fig_d015n015.png")
par(mar=c(6, 6, 1, 1))
g
dev.off()

####################################################
# Figure 4c
####################################################

n0 <- 25 # Size of the training data

# Simulate and calculate RMSEs for different training sizes
results <- map_df(1:n_rpt, 
                  ~ res(N = N,
                        n = n0,
                        k = 100), # k represents number of augmented data points
                  .id = "runs")



# Summarize the test errors across different training sizes
test_error_summary = results %>% 
  group_by(n_evaluate) %>% 
  summarize(
    rmse_train = mean(RMSE_Train, na.rm = TRUE),
    rmse_test = mean(RMSE_Test, na.rm = TRUE)
  ) %>% 
  pivot_longer(-n_evaluate, names_to = 'Type', values_to = 'RMSE') %>% 
  mutate(lossType = str_remove(Type, 'RMSE_')) 

# Plot RMSEs for training and test sets across varying training sizes
g <- results %>%
  pivot_longer(-c(runs, n_evaluate), names_to = "Type", values_to = "RMSE") %>%
  mutate(Type = str_remove(Type, "RMSE_")) %>% 
  ggplot(aes(x=n_evaluate, y=RMSE)) + 
  geom_line(aes(color=Type, group=interaction(Type, runs)), alpha=.2, linewidth=.5) +
  scale_colour_manual(values=c(2,4)) +
  stat_summary(fun=mean, geom="line", lwd=1.5, aes(group=Type, color=Type)) +
  geom_vline(xintercept=n0, color="gray",lwd=1) +
  geom_vline(xintercept=d0, color="black",lwd=1) +
  geom_hline(yintercept=sqrt(mean(y^2)), color="black",lwd=1) +
  coord_cartesian(ylim = c(0, 8)) +#max(test_error_summary$RMSE) + .5)) +
  labs(x = 'Number of training samples (n)', y = 'Error', title = "(c)") +
  guides(color = guide_legend(override.aes = list(linewidth = 3, alpha=1)))+
  geom_point(
    color = 'black',
    alpha = .5,
    size = 3,
    data = test_error_summary %>% filter(lossType == 'rmse_test') %>% 
      filter(RMSE == min(RMSE)),
    show.legend = FALSE
  ) +
  theme(legend.position = c(0.9, 0.85), 
        plot.caption = element_text(hjust = 0),
        axis.text=element_text(size=12), 
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=20))

# Display the plot
g


png(file="Fig_d015n025.png")
par(mar=c(6, 6, 1, 1))
g
dev.off()


####################################################
# Figure 4d
####################################################

# Define simulation parameters
n0 <- 5 # Size of the training data
d0 <- 15 # Number of predictors
N <- 10000 # Total number of observations
Sigma <- diag(1, d0) # Create a covariance matrix

# Generate data matrix DF and response vector y based on true parameters TPar
DF <- mvrnorm(N, mu = rep(0, d0), Sigma = Sigma)

# Generate true regression coefficients
coef <- matrix(rnorm(d0, sd = 1))
TPar <- coef# / sqrt((t(coef) %*% coef)[1,1])

# Generate response variable with added noise
y <- as.matrix(DF) %*% TPar + rnorm(N, sd = .5)

# Define a sequence of training sizes to evaluate
n_evaluate <- seq(2, 100, 1)

n_rpt <- 100 # Number of repetitions for the simulation

# Simulate and calculate RMSEs for different training sizes
results <- map_df(1:n_rpt, 
                  ~ res(N = N,
                        n = n0,
                        k = 100), # k represents number of augmented data points
                  .id = "runs")



# Summarize the test errors across different training sizes
test_error_summary = results %>% 
  group_by(n_evaluate) %>% 
  summarize(
    rmse_train = mean(RMSE_Train, na.rm = TRUE),
    rmse_test = mean(RMSE_Test, na.rm = TRUE)
  ) %>% 
  pivot_longer(-n_evaluate, names_to = 'Type', values_to = 'RMSE') %>% 
  mutate(lossType = str_remove(Type, 'RMSE_')) 

# Plot RMSEs for training and test sets across varying training sizes
g <- results %>%
  pivot_longer(-c(runs, n_evaluate), names_to = "Type", values_to = "RMSE") %>%
  mutate(Type = str_remove(Type, "RMSE_")) %>% 
  ggplot(aes(x=n_evaluate, y=RMSE)) + 
  geom_line(aes(color=Type, group=interaction(Type, runs)), alpha=.2, linewidth=.5) +
  scale_colour_manual(values=c(2,4)) +
  stat_summary(fun=mean, geom="line", lwd=1.5, aes(group=Type, color=Type)) +
  geom_vline(xintercept=n0, color="gray",lwd=1) +
  geom_vline(xintercept=d0, color="black",lwd=1) +
  geom_hline(yintercept=sqrt(mean(y^2)), color="black",lwd=1) +
  coord_cartesian(ylim = c(0, 8)) +#max(test_error_summary$RMSE) + .5)) +
  labs(x = 'Number of training samples (n)', y = 'Error', title = "(d)") +
  guides(color = guide_legend(override.aes = list(linewidth = 3, alpha=1)))+
  geom_point(
    color = 'black',
    alpha = .5,
    size = 3,
    data = test_error_summary %>% filter(lossType == 'rmse_test') %>% 
      filter(RMSE == min(RMSE)),
    show.legend = FALSE
  ) +
  theme(legend.position = c(0.9, 0.85), 
        plot.caption = element_text(hjust = 0),
        axis.text=element_text(size=12), 
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=20))
# Display the plot
g

# Save the plot
png(file="Fig_StrongSig_d015n05.png")
par(mar=c(6, 6, 1, 1))
g
dev.off()
####################################################
# Figure 4e
####################################################

n0 <- 15 # Size of the training data

# Simulate and calculate RMSEs for different training sizes
results <- map_df(1:n_rpt, 
                  ~ res(N = N,
                        n = n0,
                        k = 100), # k represents number of augmented data points
                  .id = "runs")



# Summarize the test errors across different training sizes
test_error_summary = results %>% 
  group_by(n_evaluate) %>% 
  summarize(
    rmse_train = mean(RMSE_Train, na.rm = TRUE),
    rmse_test = mean(RMSE_Test, na.rm = TRUE)
  ) %>% 
  pivot_longer(-n_evaluate, names_to = 'Type', values_to = 'RMSE') %>% 
  mutate(lossType = str_remove(Type, 'RMSE_')) 

# Plot RMSEs for training and test sets across varying training sizes
g <- results %>%
  pivot_longer(-c(runs, n_evaluate), names_to = "Type", values_to = "RMSE") %>%
  mutate(Type = str_remove(Type, "RMSE_")) %>% 
  ggplot(aes(x=n_evaluate, y=RMSE)) + 
  geom_line(aes(color=Type, group=interaction(Type, runs)), alpha=.2, linewidth=.5) +
  scale_colour_manual(values=c(2,4)) +
  stat_summary(fun=mean, geom="line", lwd=1.5, aes(group=Type, color=Type)) +
  geom_vline(xintercept=n0, color="gray",lwd=1) +
  geom_vline(xintercept=d0, color="black",lwd=1) +
  geom_hline(yintercept=sqrt(mean(y^2)), color="black",lwd=1) +
  coord_cartesian(ylim = c(0, 8)) +#max(test_error_summary$RMSE) + .5)) +
  labs(x = 'Number of training samples (n)', y = 'Error', title = "(e)") +
  guides(color = guide_legend(override.aes = list(linewidth = 3, alpha=1)))+
  geom_point(
    color = 'black',
    alpha = .5,
    size = 3,
    data = test_error_summary %>% filter(lossType == 'rmse_test') %>% 
      filter(RMSE == min(RMSE)),
    show.legend = FALSE
  ) +
  theme(legend.position = c(0.9, 0.85), 
        plot.caption = element_text(hjust = 0),
        axis.text=element_text(size=12), 
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=20))
# Display the plot
g

png(file="Fig_StrongSig_d015n015.png")
par(mar=c(6, 6, 1, 1))
g
dev.off()

####################################################
# Figure 4f
####################################################

n0 <- 25 # Size of the training data

# Simulate and calculate RMSEs for different training sizes
results <- map_df(1:n_rpt, 
                  ~ res(N = N,
                        n = n0,
                        k = 100), # k represents number of augmented data points
                  .id = "runs")



# Summarize the test errors across different training sizes
test_error_summary = results %>% 
  group_by(n_evaluate) %>% 
  summarize(
    rmse_train = mean(RMSE_Train, na.rm = TRUE),
    rmse_test = mean(RMSE_Test, na.rm = TRUE)
  ) %>% 
  pivot_longer(-n_evaluate, names_to = 'Type', values_to = 'RMSE') %>% 
  mutate(lossType = str_remove(Type, 'RMSE_')) 

# Plot RMSEs for training and test sets across varying training sizes
g <- results %>%
  pivot_longer(-c(runs, n_evaluate), names_to = "Type", values_to = "RMSE") %>%
  mutate(Type = str_remove(Type, "RMSE_")) %>% 
  ggplot(aes(x=n_evaluate, y=RMSE)) + 
  geom_line(aes(color=Type, group=interaction(Type, runs)), alpha=.2, linewidth=.5) +
  scale_colour_manual(values=c(2,4)) +
  stat_summary(fun=mean, geom="line", lwd=1.5, aes(group=Type, color=Type)) +
  geom_vline(xintercept=n0, color="gray",lwd=1) +
  geom_vline(xintercept=d0, color="black",lwd=1) +
  geom_hline(yintercept=sqrt(mean(y^2)), color="black",lwd=1) +
  coord_cartesian(ylim = c(0, 8)) +#max(test_error_summary$RMSE) + .5)) +
  labs(x = 'Number of training samples (n)', y = 'Error', title = "(f)") +
  guides(color = guide_legend(override.aes = list(linewidth = 3, alpha=1)))+
  geom_point(
    color = 'black',
    alpha = .5,
    size = 3,
    data = test_error_summary %>% filter(lossType == 'rmse_test') %>% 
      filter(RMSE == min(RMSE)),
    show.legend = FALSE
  ) +
  theme(legend.position = c(0.9, 0.85), 
        plot.caption = element_text(hjust = 0),
        axis.text=element_text(size=12), 
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=20))
# Display the plot
g

png(file="Fig_StrongSig_d015n025.png")
par(mar=c(6, 6, 1, 1))
g
dev.off()


####################################################
####################################################
# Figure 5
####################################################
####################################################

# Simulation parameters
d0 <- 15 # Number of original predictors
n0 <- 25 # Size of the training dataset
k <- 10  # Number of additional observations per variable
n_rpt <- 1000 # Number of repetitions for simulation

N <- 10000 # Total number of observations
Sigma <- diag(1, d0) # Diagonal covariance matrix for predictors

# Generate dataset
set.seed(1000) # Ensure reproducibility
DF <- as.data.frame(mvrnorm(N, mu = rep(0, d0), Sigma = Sigma))

# Generate true coefficients
coef <- matrix(rnorm(d0, sd = 2))
TPar <- coef / sqrt((t(coef) %*% coef)[1,1]) # Normalize coefficients

# Generate response variable with added noise
y <- as.matrix(DF) %*% TPar + rnorm(N, sd = .5)

# Combine predictors and response in one dataframe
myData <- data.frame(y = y, DF)

# Define a sequence of augmentation sizes to evaluate
n_evaluate <- c(0, 25, 100, 300)


# Simulate and compute estimated coefficients for each augmentation size
bhats <- sapply(n_evaluate, function(n){
  res <- replicate(n_rpt, {
    # Randomly sample training data
    train_index <- sample(1:nrow(myData), n0, replace = FALSE)
    X_train <- myData[train_index, -1]
    y_train <- myData[train_index, 1]
    
    # Augment training data if k is not zero
    if(n == 0){
      X_aug <- X_train
    } else {
      X_aug <- rbind(X_train, mvrnorm(n, mu = rep(0, d0), Sigma = diag(d0)))
    }
    y_aug <- c(y_train, sample(0, n, replace = TRUE))
    
    # Create model matrix and perform regression
    xTrain <- model.matrix(~., data = X_aug)
    yTrain <- y_aug
    
    b <- ginv(crossprod(xTrain)) %*% crossprod(xTrain, yTrain)
    b[,1]
  })
  stack(as.data.frame(t(res)))
})


# Prepare data for plotting
dat <- data.frame(X8 = 1:(d0 + 1), b = c(0, TPar[,1]))
label_text <- data.frame(
  lab = c("n-n[0]==0","n-n[0]==25","n-n[0]==100","n-n[0]==300"),
  k_eval = c("(a)", "(b)", "(c)", "(d)")
)

# Generate boxplot with facets for each augmentation size
g <- data.frame(t(plyr::ldply(bhats, rbind))) %>%
  dplyr::select(-c(X2,X4,X6)) %>%
  pivot_longer(!X8, names_to = "k_eval", values_to = "b")%>% 
  mutate(k_eval = factor(k_eval)) %>%
  mutate(k_eval = recode(k_eval, "X1"="(a)", "X3"="(b)", "X5"="(c)", "X7"="(d)")) %>%
  ggplot(aes(x = X8, y = b, group=X8)) +
  geom_boxplot(outlier.shape = NA, width=.3) + 
  facet_wrap(.~k_eval, nrow=4) +
  geom_point(data = dat,
             color = 'red',
             alpha = 1,
             size = 1,
             show.legend = FALSE
  ) +
  labs(y = expression(hat(beta)), x = "Explanatory Variable") +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    #panel.grid = element_blank()
    panel.grid.major = element_line(colour="gray", size=0.5),
    #panel.grid.major.x = element_blank(),
    #panel.grid.minor.x = element_blank(),
    #panel.grid.minor.y = element_blank()
  ) +
  xlim(0,p+1+.5) + 
  ylim(-.8,.8) +
  geom_text(x=16, y=.58, aes(label = lab), parse = TRUE, data = label_text, inherit.aes = FALSE)

# Display and save the plot
g 

png(file="FigBoxPlotCoefObs.png")
par(mar=c(6, 6, 1, 1))
g
dev.off()

###################################################################
