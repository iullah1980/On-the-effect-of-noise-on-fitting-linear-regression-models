library(tidyverse)
library(MASS)

###########################################
#  Functions used to generate data and perform computations for the following figures.
###########################################

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
    noise_vars = matrix(rnorm(nrow(myData)*k, sd=1), ncol = k)
    myData = data.frame(myData, noise_vars)
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

#################################################
# Figure 1a
#################################################

# Define simulation parameters
d0 <- 25 # Number of original variables
n0 <- 50 # Size of training data
k <- 10 # Number of additional variables per each original variable
n_rpt <- 100 # Number of repetitions
N <- 10000 # Total sample size for generating data
Sigma <- diag(1, d0) # Covariance matrix for the variables


# Generate dataset
DF <- as.data.frame(mvrnorm(N, mu=rep(0, d0), Sigma=Sigma))
coef <- matrix(rnorm(d0, sd=1)) # Coefficients for the linear model
# Normalize coefficients to unit length for interpretability
TPar <- coef / sqrt((t(coef) %*% coef)[1,1])
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
  coord_cartesian(ylim = c(0, 7)) +
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

# Save the plot
png(file="Fig_d025n050.png")
par(mar=c(6, 6, 1, 1))
g
dev.off()


#################################################
# Figure 1b
#################################################

# Define simulation parameters
d0 <- 50 # Number of original variables
n0 <- 50 # Size of training data
k <- 10 # Number of additional variables per each original variable
n_rpt <- 100 # Number of repetitions
N <- 10000 # Total sample size for generating data
Sigma <- diag(1, d0) # Covariance matrix for the variables


# Generate dataset
DF <- as.data.frame(mvrnorm(N, mu=rep(0, d0), Sigma=Sigma))
coef <- matrix(rnorm(d0, sd=1)) # Coefficients for the linear model
# Normalize coefficients to unit length for interpretability
TPar <- coef / sqrt((t(coef) %*% coef)[1,1])
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
  coord_cartesian(ylim = c(0, 7)) +
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

# Save the plot
png(file="Fig_d050n050.png")
par(mar=c(6, 6, 1, 1))
g
dev.off()


#################################################
# Figure 1c
#################################################

# Define simulation parameters
d0 <- 75 # Number of original variables
n0 <- 50 # Size of training data
k <- 10 # Number of additional variables per each original variable
n_rpt <- 100 # Number of repetitions
N <- 10000 # Total sample size for generating data
Sigma <- diag(1, d0) # Covariance matrix for the variables


# Generate dataset
DF <- as.data.frame(mvrnorm(N, mu=rep(0, d0), Sigma=Sigma))
coef <- matrix(rnorm(d0, sd=1)) # Coefficients for the linear model
# Normalize coefficients to unit length for interpretability
TPar <- coef / sqrt((t(coef) %*% coef)[1,1])
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
  coord_cartesian(ylim = c(0, 7)) +
  labs(x = 'Number of explanatory variables (d)', y = 'Error', title = "(c)") +
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

# Save the plot
png(file="Fig_d075n050.png")
par(mar=c(6, 6, 1, 1))
g
dev.off()


#################################################
# Figure 1d
#################################################


# Simulation parameters
d0 <- 25 # Number of original variables
n0 <- 50 # Size of training data
k <- 10 # Number of additional variables per each variable
n_rpt <- 100 # Number of repetitions
N <- 10000 # Total sample size for generating data

Sigma <- diag(1, d0) # Covariance matrix for the variables

# Generate dataset
# set.seed(1000) # Ensuring reproducibility
DF <- as.data.frame(mvrnorm(N, mu = rep(0, d0), Sigma = Sigma))

# Generate true parameters
coef <- matrix(rnorm(d0, sd = 1))
TPar <- coef # No normalization applied here

# Generate response variable
y <- as.matrix(DF) %*% TPar + rnorm(N, sd = .5)

# Prepare dataset
myData <- data.frame(y = y, DF)


# Run simulations and calculate RMSE for training and testing sets
results <- map_df(1:n_rpt, 
                  ~ train_test_rmses(n = n0, 
                                     k = k, 
                                     myData = myData),
                  .id = "runs")


# Summarize test errors
test_error_summary <- results %>% 
  group_by(d_evaluate) %>% 
  summarize(rmse_train = mean(lossTrain, na.rm = TRUE),
            rmse_test = mean(lossTest, na.rm = TRUE)) %>% 
  pivot_longer(-d_evaluate, names_to = 'lossType', values_to = 'RMSE') %>% 
  mutate(lossType = str_remove(lossType, 'rmse_'))


# Visualize RMSEs
g <- results %>%
  pivot_longer(-c(runs, d_evaluate), names_to = "Type", values_to = "RMSE") %>%
  mutate(Type = str_remove(Type, "loss")) %>% 
  ggplot(aes(x=d_evaluate, y=RMSE)) + 
  geom_line(aes(color=Type, group=interaction(Type, runs)), alpha=.2, size=.5) +
  scale_colour_manual(values=c(2,4)) +
  stat_summary(fun=mean, geom="line", lwd=1.5, aes(group=Type, color=Type)) +
  geom_vline(xintercept=n0, color="gray",lwd=1) +
  geom_vline(xintercept=d0, color="black",lwd=1) +
  geom_hline(yintercept=sqrt(mean(y^2)), color="black",lwd=1) +
  coord_cartesian(ylim = c(0, 15)) +#max(test_error_summary$RMSE) + .5)) +
  labs(x = 'Number of explanatory variables (d)', y = 'Error', title = "(d)") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha=1)))+
  geom_point(
    color = 'black',
    alpha = .5,
    size = 3,
    data = test_error_summary %>% filter(lossType == 'test') %>% 
      filter(RMSE == min(RMSE)),
    show.legend = FALSE
  ) +
  theme(legend.position = c(0.9, 0.85), plot.caption = element_text(hjust = 0))

# Display and save the plot
g

# Save the plot
png(file="Fig_StrongSig_d025n050.png")
par(mar=c(6, 6, 1, 1))
g
dev.off()

#################################################
# Figure 1e
#################################################


# Simulation parameters
d0 <- 50 # Number of original variables
n0 <- 50 # Size of training data
k <- 10 # Number of additional variables per each variable
n_rpt <- 100 # Number of repetitions
N <- 10000 # Total sample size for generating data

Sigma <- diag(1, d0) # Covariance matrix for the variables

# Generate dataset
# set.seed(1000) # Ensuring reproducibility
DF <- as.data.frame(mvrnorm(N, mu = rep(0, d0), Sigma = Sigma))

# Generate true parameters
coef <- matrix(rnorm(d0, sd = 1))
TPar <- coef # No normalization applied here

# Generate response variable
y <- as.matrix(DF) %*% TPar + rnorm(N, sd = .5)

# Prepare dataset
myData <- data.frame(y = y, DF)


# Run simulations and calculate RMSE for training and testing sets
results <- map_df(1:n_rpt, 
                  ~ train_test_rmses(n = n0, 
                                     k = k, 
                                     myData = myData),
                  .id = "runs")


# Summarize test errors
test_error_summary <- results %>% 
  group_by(d_evaluate) %>% 
  summarize(rmse_train = mean(lossTrain, na.rm = TRUE),
            rmse_test = mean(lossTest, na.rm = TRUE)) %>% 
  pivot_longer(-d_evaluate, names_to = 'lossType', values_to = 'RMSE') %>% 
  mutate(lossType = str_remove(lossType, 'rmse_'))


# Visualize RMSEs
g <- results %>%
  pivot_longer(-c(runs, d_evaluate), names_to = "Type", values_to = "RMSE") %>%
  mutate(Type = str_remove(Type, "loss")) %>% 
  ggplot(aes(x=d_evaluate, y=RMSE)) + 
  geom_line(aes(color=Type, group=interaction(Type, runs)), alpha=.2, size=.5) +
  scale_colour_manual(values=c(2,4)) +
  stat_summary(fun=mean, geom="line", lwd=1.5, aes(group=Type, color=Type)) +
  geom_vline(xintercept=n0, color="gray",lwd=1) +
  geom_vline(xintercept=d0, color="black",lwd=1) +
  geom_hline(yintercept=sqrt(mean(y^2)), color="black",lwd=1) +
  coord_cartesian(ylim = c(0, 15)) +#max(test_error_summary$RMSE) + .5)) +
  labs(x = 'Number of explanatory variables (d)', y = 'Error', title = "(e)") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha=1)))+
  geom_point(
    color = 'black',
    alpha = .5,
    size = 3,
    data = test_error_summary %>% filter(lossType == 'test') %>% 
      filter(RMSE == min(RMSE)),
    show.legend = FALSE
  ) +
  theme(legend.position = c(0.9, 0.85), plot.caption = element_text(hjust = 0))

# Display and save the plot
g

# Save the plot
png(file="Fig_StrongSig_d050n050.png")
par(mar=c(6, 6, 1, 1))
g
dev.off()

#################################################
# Figure 1f
#################################################


# Simulation parameters
d0 <- 75 # Number of original variables
n0 <- 50 # Size of training data
k <- 10 # Number of additional variables per each variable
n_rpt <- 100 # Number of repetitions
N <- 10000 # Total sample size for generating data

Sigma <- diag(1, d0) # Covariance matrix for the variables

# Generate dataset
# set.seed(1000) # Ensuring reproducibility
DF <- as.data.frame(mvrnorm(N, mu = rep(0, d0), Sigma = Sigma))

# Generate true parameters
coef <- matrix(rnorm(d0, sd = 1))
TPar <- coef # No normalization applied here

# Generate response variable
y <- as.matrix(DF) %*% TPar + rnorm(N, sd = .5)

# Prepare dataset
myData <- data.frame(y = y, DF)


# Run simulations and calculate RMSE for training and testing sets
results <- map_df(1:n_rpt, 
                  ~ train_test_rmses(n = n0, 
                                     k = k, 
                                     myData = myData),
                  .id = "runs")


# Summarize test errors
test_error_summary <- results %>% 
  group_by(d_evaluate) %>% 
  summarize(rmse_train = mean(lossTrain, na.rm = TRUE),
            rmse_test = mean(lossTest, na.rm = TRUE)) %>% 
  pivot_longer(-d_evaluate, names_to = 'lossType', values_to = 'RMSE') %>% 
  mutate(lossType = str_remove(lossType, 'rmse_'))


# Visualize RMSEs
g <- results %>%
  pivot_longer(-c(runs, d_evaluate), names_to = "Type", values_to = "RMSE") %>%
  mutate(Type = str_remove(Type, "loss")) %>% 
  ggplot(aes(x=d_evaluate, y=RMSE)) + 
  geom_line(aes(color=Type, group=interaction(Type, runs)), alpha=.2, size=.5) +
  scale_colour_manual(values=c(2,4)) +
  stat_summary(fun=mean, geom="line", lwd=1.5, aes(group=Type, color=Type)) +
  geom_vline(xintercept=n0, color="gray",lwd=1) +
  geom_vline(xintercept=d0, color="black",lwd=1) +
  geom_hline(yintercept=sqrt(mean(y^2)), color="black",lwd=1) +
  coord_cartesian(ylim = c(0, 15)) +#max(test_error_summary$RMSE) + .5)) +
  labs(x = 'Number of explanatory variables (d)', y = 'Error', title = "(f)") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha=1)))+
  geom_point(
    color = 'black',
    alpha = .5,
    size = 3,
    data = test_error_summary %>% filter(lossType == 'test') %>% 
      filter(RMSE == min(RMSE)),
    show.legend = FALSE
  ) +
  theme(legend.position = c(0.9, 0.85), plot.caption = element_text(hjust = 0))

# Display and save the plot
g

# Save the plot
png(file="Fig_StrongSig_d075n050.png")
par(mar=c(6, 6, 1, 1))
g
dev.off()

####################################################
####################################################
# Figure 3
####################################################
####################################################

# Simulation parameters
d0 <- 50 # Number of original variables
n0 <- 50 # Size of training data
k <- 10  # Number of additional variables per each variable
n_rpt <- 1000 # Number of shuffles

N <- 10000
Sigma <- diag(1, d0) # Covariance matrix for the variables

# Generate dataset
# set.seed(1000)
DF <- as.data.frame(mvrnorm(N, mu = rep(0, d0), Sigma = Sigma))

# Generate true coefficients
coef <- matrix(rnorm(d0, sd = 2))
TPar <- coef / sqrt((t(coef) %*% coef)[1,1])

# Generate response variable
y <- as.matrix(DF) %*% TPar + rnorm(N, sd = .5)

# Prepare dataset
myData <- data.frame(y = y, DF)

# Evaluation points
d_evaluate <- c(25, 75, 100, 300)

# Estimate coefficients
bhats <- sapply(d_evaluate, function(eval_indx){
  
  res <- replicate(n_rpt, {
    
    myData_aug <- myData
    
    n_cols_myData <- ncol(myData)-1
    
    # Randomly select and reorder columns for model inclusion
    permute_col <- sample(2:ncol(myData), min(d0, eval_indx))
    
    ix <- sort(permute_col, index.return=TRUE)$ix
    myData_aug <- myData_aug[, c(1, permute_col)]
    #permute_col[ix]
    
    for (cvar in 1:n_cols_myData){
      noise_vars = matrix(rnorm(nrow(myData_aug)*k, sd=1), ncol = k)
      myData_aug = data.frame(myData_aug, noise_vars)
    }
    
    
    X <- model.matrix(~., data=myData_aug[, -c(1)])
    y <- myData_aug[, 1]
    
    # Training set selection
    train_index <- sample(1:nrow(myData_aug), n0, replace = FALSE)
    xTrain = X[train_index, 1:(eval_indx+1)]
    yTrain = y[train_index]
    
    
    bhat <- numeric(max(d0, eval_indx)+1)
    # Coefficient estimation
    b <- ginv(crossprod(xTrain))%*%crossprod(xTrain, yTrain)
    # b <- solve(crossprod(xTrain)+.00000001*diag(ncol(xTrain)))%*%crossprod(xTrain, yTrain)
    # Store estimated coefficients
    if(eval_indx > d0){
      bhat[c(1,permute_col,(length(permute_col)+2):(eval_indx+1))]=b[,1]
    } else {bhat[c(1,permute_col)]=b[,1]}
    bhat
    
  })
  stack(as.data.frame(t(res)))
})

# Preparing data for plotting
dat <- data.frame(X8=1:(d0+1), b=c(0, TPar[,1]))
label_text  <- data.frame(lab=c("d-d[0]==0","d-d[0]==25","d-d[0]==50","d-d[0]==250"),
                          d_eval = c("(a)","(b)","(c)","(d)"))

# Plotting
g <- data.frame(t(plyr::ldply(bhats, rbind))) %>%
  dplyr::select(-c(X2,X4,X6)) %>%
  mutate(X1 = na_if(X1,0), X3 = na_if(X3,0), X5 = na_if(X5,0), X7 = na_if(X7,0)) %>%
  pivot_longer(!X8, names_to = "d_eval", values_to = "b")%>% 
  mutate(d_eval = factor(d_eval)) %>%
  mutate(d_eval = recode(d_eval, "X1"="(a)", "X3"="(b)", "X5"="(c)", "X7"="(d)")) %>%
  ggplot(aes(x = X8, y = b, group=X8)) +
  geom_boxplot(outlier.shape = NA, width=.3) + 
  facet_wrap(.~d_eval, nrow=4) +
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
    panel.grid.major = element_line(colour="gray", size=0.5),
  ) +
  xlim(0,60.5) + 
  #ylim(-.6,.6) +
  geom_text(x=56, y=.5, aes(label = lab), parse = TRUE, data = label_text, inherit.aes = FALSE)

g 

# Save the plot
png(file="FigBoxPlotCoef.png")
par(mar=c(6, 6, 1, 1))
g
dev.off()

#################################################
#################################################
# Figure S2a
#################################################
#################################################

# Define simulation parameters
n0 <- 50 # Size of training data
k <- 10 # Number of additional variables per each original variable
d0 <- 25 # Number of original variables
n_rpt <- 100 # Number of repetitions
N <- 10000 # Total number of observations
Sigma <- diag(1, d0) # Covariance matrix for the original variables

# Generate dataset
DF <- mvrnorm(N, mu=rep(0, d0), Sigma=Sigma)

# Generate coefficients and response variable
coef <- matrix(rnorm(d0, sd=1))
TPar <- coef / sqrt((t(coef) %*% coef)[1,1])
y <- as.matrix(DF) %*% TPar + rnorm(N, sd=1)

myData <- DF
n_cols_myData <- dim(DF)[2]

for (cvar in 1:n_cols_myData){
  noise_vars = matrix(rnorm(nrow(myData)*k, sd=1), ncol = k)
  myData = data.frame(myData, noise_vars)
}

# Prepare d_evaluate vector for exploring condition numbers across d values
d_evaluate <- seq(1, 275, 1)

# Calculate mean generalized condition number for each value in d_evaluate
mean_CNo <- sapply(d_evaluate, function(d) {
  CNos <- replicate(n_rpt, {
    train_rows <- sample(1:N, n0)
    X_train <- myData[train_rows, 1:d]
    
    SVal <- svd(X_train)$d
    # Generalized condition number
    CNo <- max(SVal) / min(SVal)
    
    CNo
  })
  mean(CNos)
})

g <- data.frame(cbind(mean_CNo,d_evaluate)) %>%
  ggplot(aes(x=d_evaluate, mean_CNo)) + 
  geom_line(alpha=.8, size=2, color="red") +
  labs(x = 'Number of explanatory variables (d)', 
       y = 'Generalized condition number', title = "(a)") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size=14)) +
  xlim(0,275)+
  geom_vline(xintercept=n0, color="gray",lwd=1) +
  geom_vline(xintercept=d0, color="black",lwd=1)

g # Display the plot

# Save the plot
png(file="Fig_d025n050_CNo.png")
par(mar=c(6, 6, 1, 1))
g
dev.off()


#################################################
# Figure S2b
#################################################
# Define simulation parameters
n0 <- 50 # Size of training data
k <- 10 # Number of additional variables per each original variable
d0 <- 50 # Number of original variables
n_rpt <- 100 # Number of repetitions
N <- 10000 # Total number of observations
Sigma <- diag(1, d0) # Covariance matrix for the original variables

# Generate dataset
DF <- mvrnorm(N, mu=rep(0, d0), Sigma=Sigma)

# Generate coefficients and response variable
coef <- matrix(rnorm(d0, sd=1))
TPar <- coef / sqrt((t(coef) %*% coef)[1,1])
y <- as.matrix(DF) %*% TPar + rnorm(N, sd=1)

myData <- DF
n_cols_myData <- dim(DF)[2]

for (cvar in 1:n_cols_myData){
  noise_vars = matrix(rnorm(nrow(myData)*k, sd=1), ncol = k)
  myData = data.frame(myData, noise_vars)
}

# Prepare d_evaluate vector for exploring condition numbers across d values
d_evaluate <- seq(1, n_cols_myData+sp*k, 1)

# Calculate mean generalized condition number for each value in d_evaluate
mean_CNo <- sapply(d_evaluate, function(d) {
  CNos <- replicate(n_rpt, {
    train_rows <- sample(1:N, n0)
    X_train <- myData[train_rows, 1:d]
    
    SVal <- svd(X_train)$d
    # Generalized condition number
    CNo <- max(SVal) / min(SVal)
    
    CNo
  })
  mean(CNos)
})

g <- data.frame(cbind(mean_CNo,d_evaluate)) %>%
  ggplot(aes(x=d_evaluate, mean_CNo)) + 
  geom_line(alpha=.8, size=2, color="red") +
  labs(x = 'Number of explanatory variables (d)', 
       y = 'Generalized condition number', title = "(b)") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size=14)) +
  xlim(0,275)+
  geom_vline(xintercept=n0, color="gray",lwd=1) +
  geom_vline(xintercept=d0, color="black",lwd=1)

g # Display the plot

# Save the plot

png(file="Fig_d050n050_CNo.png")
par(mar=c(6, 6, 1, 1))
g
dev.off()

#################################################
# Figure S2c
#################################################

# Define simulation parameters
n0 <- 50 # Size of training data
k <- 10 # Number of additional variables per each original variable
d0 <- 75 # Number of original variables
n_rpt <- 100 # Number of repetitions
N <- 10000 # Total number of observations
Sigma <- diag(1, d0) # Covariance matrix for the original variables

# Generate dataset
DF <- mvrnorm(N, mu=rep(0, d0), Sigma=Sigma)

# Generate coefficients and response variable
coef <- matrix(rnorm(d0, sd=1))
TPar <- coef / sqrt((t(coef) %*% coef)[1,1])
y <- as.matrix(DF) %*% TPar + rnorm(N, sd=1)

myData <- DF
n_cols_myData <- dim(DF)[2]

for (cvar in 1:n_cols_myData){
  noise_vars = matrix(rnorm(nrow(myData)*k, sd=1), ncol = k)
  myData = data.frame(myData, noise_vars)
}

# Prepare d_evaluate vector for exploring condition numbers across d values
d_evaluate <- seq(1, n_cols_myData+sp*k, 1)

# Calculate mean generalized condition number for each value in d_evaluate
mean_CNo <- sapply(d_evaluate, function(d) {
  CNos <- replicate(n_rpt, {
    train_rows <- sample(1:N, n0)
    X_train <- myData[train_rows, 1:d]
    
    SVal <- svd(X_train)$d
    # Generalized condition number
    CNo <- max(SVal) / min(SVal)
    
    CNo
  })
  mean(CNos)
})

g <- data.frame(cbind(mean_CNo,d_evaluate)) %>%
  ggplot(aes(x=d_evaluate, mean_CNo)) + 
  geom_line(alpha=.8, size=2, color="red") +
  labs(x = 'Number of explanatory variables (d)', 
       y = 'Generalized condition number', title = "(c)") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size=14)) +
  xlim(0,275)+
  geom_vline(xintercept=n0, color="gray",lwd=1) +
  geom_vline(xintercept=d0, color="black",lwd=1)

g # Display the plot

# Save the plot
png(file="Fig_d075n050_CNo.png")
par(mar=c(6, 6, 1, 1))
g
dev.off()

