
library(tidyverse)
library(MASS)

who <- read.csv("C:\\Users\\iulla\\Downloads\\Stanford-AI-Alignment-Double-Descent-Tutorial-master\\Stanford-AI-Alignment-Double-Descent-Tutorial-master\\data\\Life Expectancy Data.csv")

# Clean data: drop non-numeric and missing
df <- who %>% 
  dplyr::select(-Country, -Year, -Status) %>%
  drop_na()

# Define X and y
X <- scale(dplyr::select(df, -Life.expectancy))
y <- scale(df$Life.expectancy, center = TRUE, scale = FALSE)

# Ridgeless regression function
ridgeless_reg <- function(X_train, X_test, y_train, y_test){
  b <- MASS::ginv(crossprod(X_train)) %*% crossprod(X_train, y_train)
  predTrain <- X_train %*% b
  predTest <- X_test %*% b
  RMSE_train <- sqrt(mean((y_train - predTrain)^2))
  RMSE_test <- sqrt(mean((y_test - predTest)^2))
  list(RMSE_Train = RMSE_train, RMSE_Test = RMSE_test)
}

# Evaluation function for WHO
evaluate_who <- function(X, y, fixed_train_size, n_evaluate){
  n_total <- nrow(X)
  d <- ncol(X)
  train_rows <- sample(1:n_total, fixed_train_size)
  
  X_train <- X[train_rows, ]
  y_train <- y[train_rows]
  
  max_aug <- max(n_evaluate) - fixed_train_size
  X_aug <- rbind(X_train, MASS::mvrnorm(n = max_aug, mu = rep(0, d), Sigma = diag(d)))
  y_aug <- c(y_train, sample(0, max_aug, replace = TRUE))
  
  X_test <- X[-train_rows, ]
  y_test <- y[-train_rows]
  
  map_df(n_evaluate, ~{
    res <- ridgeless_reg(X_aug[1:.x, ], X_test, y_aug[1:.x], y_test)
    tibble(n_evaluate = .x, RMSE_Train = res$RMSE_Train, RMSE_Test = res$RMSE_Test)
  })
}

# Parameters
fixed_train_size <- 80
n_total <- nrow(X)
n_evaluate <- seq(2, 200, by = 2)
num_repeats <- 50

# Run simulation
#set.seed(123)
results <- map_df(
  1:num_repeats,
  ~ evaluate_who(X, y, fixed_train_size, n_evaluate),
  .id = "runs"
)

# Summarise RMSE
summary_df <- results %>%
  group_by(n_evaluate) %>%
  summarise(
    rmse_train = mean(RMSE_Train, na.rm = TRUE),
    rmse_test = mean(RMSE_Test, na.rm = TRUE)
  ) %>%
  pivot_longer(-n_evaluate, names_to = "Type", values_to = "RMSE") %>%
  mutate(lossType = str_remove(Type, "rmse_"))

# Plot
g <- results %>%
  pivot_longer(-c(runs, n_evaluate), names_to = "Type", values_to = "RMSE") %>%
  mutate(Type = str_remove(Type, "RMSE_")) %>%
  ggplot(aes(x = n_evaluate, y = RMSE)) +
  geom_line(aes(color = Type, group = interaction(Type, runs)), alpha = .2, linewidth = .5) +
  scale_colour_manual(values = c(2, 4)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.5, aes(group = Type, color = Type)) +
  geom_vline(xintercept = fixed_train_size, color = "gray", linewidth = 1) +
  geom_vline(xintercept = ncol(X), color = "black", linewidth = 1) +
  geom_hline(yintercept = sqrt(mean(y^2)), color = "black", linewidth = 1) +
  coord_cartesian(ylim = c(0, 20)) +
  labs(x = 'Number of training samples (n)', y = 'Error', title = "(WHO example, fixed train size = 50)") +
  guides(color = guide_legend(override.aes = list(linewidth = 3, alpha = 1))) +
  geom_point(
    color = 'black',
    alpha = .5,
    size = 3,
    data = summary_df %>% filter(lossType == 'test') %>% filter(RMSE == min(RMSE)),
    show.legend = FALSE
  ) +
  theme(
    legend.position = c(0.9, 0.85),
    plot.caption = element_text(hjust = 0),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 20)
  )

# Show plot
print(g)

