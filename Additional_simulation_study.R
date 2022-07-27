rm(list = ls())

source("core_functions/core.R", echo = TRUE)
Rcpp::sourceCpp("core_functions/core.cpp")

# ----------------------------------------------
# PRELIMINARY GRAPHS
# ----------------------------------------------

# Simulated datasets
n <- 100
TT <- 50
time <- (1:TT) / TT
dataset_std <- matrix(0, n, TT)

clusters <- rep(1:4, each = 25)

# FUNCTIONS DEFINITION
set.seed(123)
f1 <- 1 - 2 * time
f2 <- 0.5 * cos(2 * pi * time) + 0.5 * sin(2 * pi * time)
f3 <- -1 + 2 * time^4
f4 <- 0.5 * cos(4 * pi * time) + 0.5 * sin(4 * pi * time)

par(mfrow = c(2, 2))
plot(f1, type = "l")
plot(f2, type = "l")
plot(f3, type = "l")
plot(f4, type = "l")
par(mfrow = c(1, 1))

# -----------------------------------------------
# Hyperparameter settings
# -----------------------------------------------

L <- 4
B <- cbind(1, time, time^4, cos(2 * pi * time), sin(2 * pi * time), cos(2 * pi * 2 * time), sin(2 * pi * 2 * time))
indexB <- matrix(c(
  TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE,
  TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE,
  TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE,
  TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE
), 4, 7, byrow = TRUE)

cc <- c(1, 1, 1, 1)
Hl <- c(5, 5, 5, 5)
alpha <- c(1, 1, 1, 1)

Sigma <- list(diag(10, 2), diag(10, 3), diag(10, 2), diag(10, 3))
a_sigma <- 1
b_sigma <- 1

# -----------------------------------------------
# HIGH SIGNAL SCENARIO
# -----------------------------------------------

SE <- 0.1

dataset_std <- matrix(0, n, TT)

for (i in 1:n) {
  if (1 <= i && i <= 25) dataset_std[i, ] <- f1
  if (26 <= i && i <= 50) dataset_std[i, ] <- f2
  if (51 <= i && i <= 75) dataset_std[i, ] <- f3
  if (76 <= i && i <= 100) dataset_std[i, ] <- f4
}

set.seed(123)
dataset_std <- dataset_std + matrix(rnorm(n * TT, 0, SE), n, TT)
# dataset_std <-  dataset_std[sample(1:n),]

# EXECUTION OF THE ALGORITHM
library(foreach)

# library(doMC)      # For UNIX machines

registerDoSEQ()

# Select different seeds
set.seed(555)
seeds <- sample(1:2^15, 5) # Convergence occurs always.

lower_bounds <- foreach(seeds = seeds, .combine = rbind) %dopar% {
  set.seed(seeds)
  out <- EFDMP(
    X = dataset_std, Hl = Hl, L = L, B = B, indexB = indexB, time = time, cc = cc,
    alpha = alpha, Sigma = Sigma, a_sigma = a_sigma, b_sigma = b_sigma, verbose = TRUE, maxiter = 1500
  )
  Rsquared <- 1 - sum((dataset_std - out$prediction)^2) / sum((dataset_std - mean(dataset_std))^2)
  Kunique <- length(unique(out$cluster))
  c(seeds, out$lowerbound, Rsquared, Kunique)
}

set.seed(lower_bounds[which.max(lower_bounds[, 2]), 1])
fit <- EFDMP(
  X = dataset_std, Hl = Hl, L = L, B = B, indexB = indexB, time = time, cc = cc,
  alpha = alpha, Sigma = Sigma, a_sigma = a_sigma, b_sigma = b_sigma, verbose = TRUE, maxiter = 1500
)

lower_bounds
table(fit$Fclass)
table(fit$cluster)
table(clusters, fit$Fclass)
table(clusters, fit$cluster)

tab1 <- table(clusters, fit$cluster)
library(xtable)
print(xtable(tab1[, colSums(tab1) > 0], format = "latex"), booktabs = TRUE)

library(ggplot2)
data_plot <- reshape2::melt(dataset_std)
data_plot$Cluster <- factor(data_plot$Var1)
levels(data_plot$Cluster) <- fit$cluster
data_plot$Cluster <- factor(as.numeric(data_plot$Cluster))

selected_cluster <- as.numeric(unique(fit$cluster))
data_plot2 <- reshape2::melt(fit$pred_cluster[selected_cluster, ])
p1 <- ggplot(data = data_plot, aes(x = Var2 / TT, y = value, group = as.factor(Var1), col = Cluster)) +
  geom_line(alpha = 0.3) +
  theme_bw() +
  scale_colour_brewer(palette = "Set1") +
  ylab("y") +
  xlab("Time") +
  geom_line(data = data_plot2, size = 1.1, aes(x = Var2 / TT, y = value, col = as.factor(Var1))) +
  ggtitle("Small variance scenario") +
  theme(legend.position = "none")
p1


# -----------------------------------------------
# LOW SIGNAL SCENARIO
# -----------------------------------------------

SE <- 1.5

dataset_std <- matrix(0, n, TT)

for (i in 1:n) {
  if (1 <= i && i <= 25) dataset_std[i, ] <- f1
  if (26 <= i && i <= 50) dataset_std[i, ] <- f2
  if (51 <= i && i <= 75) dataset_std[i, ] <- f3
  if (76 <= i && i <= 100) dataset_std[i, ] <- f4
}

set.seed(123)
dataset_std <- dataset_std + matrix(rnorm(n * TT, 0, SE), n, TT)

# EXECUTION OF THE ALGORITHM

# Select different seeds
set.seed(555)
seeds <- sample(1:2^15, 5) # Convergence occurs always.

lower_bounds <- foreach(seeds = seeds, .combine = rbind) %dopar% {
  set.seed(seeds)
  out <- EFDMP(
    X = dataset_std, Hl = Hl, L = L, B = B, indexB = indexB, time = time, cc = cc,
    alpha = alpha, Sigma = Sigma, a_sigma = a_sigma, b_sigma = b_sigma, verbose = TRUE, maxiter = 1500
  )
  Rsquared <- 1 - sum((dataset_std - out$prediction)^2) / sum((dataset_std - mean(dataset_std))^2)
  Kunique <- length(unique(out$cluster))
  return(c(seeds, out$lowerbound, Rsquared, Kunique))
}

set.seed(lower_bounds[which.max(lower_bounds[, 2]), 1])
fit <- EFDMP(
  X = dataset_std, Hl = Hl, L = L, B = B, indexB = indexB, time = time, cc = cc,
  alpha = alpha, Sigma = Sigma, a_sigma = a_sigma, b_sigma = b_sigma, verbose = TRUE, maxiter = 1500
)

lower_bounds
table(fit$Fclass)
table(fit$cluster)
table(clusters, fit$Fclass)
table(clusters, fit$cluster)

tab2 <- table(clusters, fit$cluster)
library(xtable)
print(xtable(tab2[, colSums(tab2) > 0], format = "latex"), booktabs = TRUE)

# Due to labeling issues, the errors must be computed by hand.

library(ggplot2)
data_plot <- reshape2::melt(dataset_std)
data_plot$Cluster <- factor(data_plot$Var1)
levels(data_plot$Cluster) <- fit$cluster
data_plot$Cluster <- factor(as.numeric(data_plot$Cluster))

selected_cluster <- as.numeric(unique(fit$cluster))
data_plot2 <- reshape2::melt(fit$pred_cluster[selected_cluster, ])
p2 <- ggplot(data = data_plot, aes(x = Var2 / TT, y = value, group = as.factor(Var1), col = Cluster)) +
  geom_line(alpha = 0.2) +
  theme_bw() +
  scale_colour_brewer(palette = "Set1") +
  ylab("y") +
  xlab("Time") +
  geom_line(data = data_plot2, size = 1.1, aes(x = Var2 / TT, y = value, col = as.factor(Var1))) +
  ggtitle("High variance scenario") +
  theme(legend.position = "none")
p2

library(gridExtra)
grid.arrange(p1, p2, ncol = 1)
# ggsave("/Volumes/Macfiles/Google Drive/University/Lavori/Flight_route_segmentation/BA_tex/img/crossval.pdf",width=10,height=8)
