rm(list = ls())

source("core_functions/core.R", echo = TRUE)
Rcpp::sourceCpp("core_functions/core.cpp")

# Preliminary quantities -----------------------------------------
n <- 100
TT <- 50
time <- (1:TT) / TT

B <- cbind(
  1, time, time^4,
  cos(2 * pi * time),
  sin(2 * pi * time),
  cos(2 * pi * 2 * time),
  sin(2 * pi * 2 * time)
)
indexB <- matrix(c(
  TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE,
  TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE,
  TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE,
  TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE
), 4, 7, byrow = TRUE)


# Specification of the true curves ----------------------------
set.seed(1234)
f_curves <- matrix(0, nrow = 16, ncol = TT)
par(mfrow = c(1, 1))
plot(time, rep(0, TT), type = "l", lty = "dotted", ylim = c(-5, 5))

f_curves[1, ] <- 1 - 2 * time
f_curves[2, ] <- -1 + time
f_curves[3, ] <- 1 - time
f_curves[4, ] <- -1 + 2 * time
f_curves[5, ] <- -1 + 2 * time^4
f_curves[6, ] <- 1 - 2 * time^4
f_curves[7, ] <- -1.5 + 4 * time^4
f_curves[8, ] <- 1.5 - 4 * time^4
f_curves[9, ] <- 0.5 * cos(2 * pi * time) + 0.5 * sin(2 * pi * time)
f_curves[10, ] <- 0.25 * cos(2 * pi * time) + 0.75 * sin(2 * pi * time)
f_curves[11, ] <- 0.75 * cos(2 * pi * time) + 0.25 * sin(2 * pi * time)
f_curves[12, ] <- sin(2 * pi * time)
f_curves[13, ] <- 0.5 * cos(4 * pi * time) + 0.5 * sin(4 * pi * time)
f_curves[14, ] <- 0.25 * cos(4 * pi * time) + 0.75 * sin(4 * pi * time)
f_curves[15, ] <- 0.75 * cos(4 * pi * time) + 0.25 * sin(4 * pi * time)
f_curves[16, ] <- sin(4 * pi * time)

for (i in 1:4) {
  for (j in 1:4) {
    lines(time, f_curves[4 * (i - 1) + j, ], col = 4 * (i - 1) + j)
  }
}

# Generation of the dataset ---------------------
SE <- 0.3
dataset_signal <- matrix(0, n, TT)
clusters <- rep(0, n)
set.seed(100)

for (i in 1:n) {
  id <- sample(1:16, 1) # Uniform probabilities
  clusters[i] <- id
  dataset_signal[i, ] <- f_curves[id, ]
}

set.seed(123)
dataset_high_signal <- dataset_signal + matrix(rnorm(n * TT, 0, SE), n, TT)

library(ggplot2)
data_plot <- reshape2::melt(dataset_high_signal)
data_plot2 <- reshape2::melt(f_curves)
p0 <- ggplot(data = data_plot, aes(x = Var2, y = value, group = as.factor(Var1))) +
  geom_line(alpha = 0.2, col = "gray") +
  theme_bw() +
  scale_colour_brewer(palette = "Set1") +
  ylab("y") +
  xlab("Time") +
  geom_line(data = data_plot2, size = 0.7, aes(x = Var2, y = value, group = as.factor(Var1))) +
  theme(legend.position = "none")
p0

# Estimation E-FDMP----------------------------------------

library(splines)
L <- 4

# List of hyperparameters
cc <- c(1, 1, 1, 1)
alpha <- c(1, 1, 1, 1)
Sigma <- list(diag(10, 2), diag(10, 2), diag(10, 3), diag(10, 3))
a_sigma <- 1
b_sigma <- 1
Hl <- c(20, 20, 20, 20)

# Execution of the algorithm ---------------------------------
library(doMC)
library(foreach)
registerDoMC(6)

# Select different seeds
set.seed(555)
seeds <- sample(1:2^15, 50) # Convergence occurs always.

lower_bounds <- foreach(seeds = seeds, .combine = rbind) %dopar% {
  set.seed(seeds)
  out <- EFDMP(
    X = dataset_high_signal, Hl = Hl, L = L, B = B, indexB = indexB, time = time, cc = cc,
    alpha = alpha, Sigma = Sigma, a_sigma = a_sigma, b_sigma = b_sigma, verbose = TRUE, maxiter = 1500, tol = 1e-12
  )
  Rsquared <- 1 - sum((dataset_high_signal - out$prediction)^2) / sum((dataset_high_signal - mean(dataset_high_signal))^2)
  Kunique <- length(unique(out$cluster))
  c(seeds, out$lowerbound, Rsquared, Kunique)
}

set.seed(lower_bounds[which.max(lower_bounds[, 2]), 1])
fit <- EFDMP(
  X = dataset_high_signal, Hl = Hl, L = L, B = B, indexB = indexB, time = time, cc = cc,
  alpha = alpha, Sigma = Sigma, a_sigma = a_sigma, b_sigma = b_sigma, verbose = TRUE, maxiter = 1500
)

# Display some useful quantities
lower_bounds
length(unique(fit$cluster))
table(fit$Fclass)
table(fit$cluster)
table(clusters, fit$Fclass)
table(clusters, fit$cluster)

# Plot generation
data_plot <- reshape2::melt(dataset_high_signal)
data_plot$Cluster <- factor(data_plot$Var1)

selected_cluster <- as.numeric(unique(fit$cluster))
data_plot2 <- fit$pred_cluster[selected_cluster, ]
rownames(data_plot2) <- selected_cluster
data_plot2 <- reshape2::melt(data_plot2)
data_plot2$Var3 <- fit$Findex[data_plot2$Var1]
p1 <- ggplot(data = data_plot, aes(x = Var2, y = value, group = as.factor(Var1))) +
  geom_line(alpha = 0.2, col = "gray") +
  theme_bw() +
  scale_colour_brewer(palette = "Set1") +
  ylab("y") +
  xlab("Time") +
  geom_line(data = data_plot2, size = 0.7, aes(x = Var2, y = value, group = as.factor(Var1), col = as.factor(Var3))) +
  theme(legend.position = "none")
p1

# Mixture parametric ------------------------------------------------

L <- 1
B <- cbind(
  1, time, time^4,
  cos(2 * pi * time),
  sin(2 * pi * time),
  cos(2 * pi * 2 * time),
  sin(2 * pi * 2 * time)
)
indexB <- matrix(rep(TRUE, 7), 1, 7, byrow = TRUE)

# List of hyperparameters
L <- 1
cc <- 1
alpha <- 1
Sigma <- list(diag(rep(10, 7)))
a_sigma <- 1
b_sigma <- 1
Hl <- 80

# Execution of the algorithm ---------------------------------
library(doMC)
library(foreach)
registerDoMC(6)

# Select different seeds
set.seed(555)
seeds <- sample(1:2^15, 50) # Convergence occurs always.

lower_bounds2 <- foreach(seeds = seeds, .combine = rbind) %dopar% {
  set.seed(seeds)
  out <- EFDMP(
    X = dataset_high_signal, Hl = Hl, L = L, B = B, indexB = indexB, time = time, cc = cc,
    alpha = alpha, Sigma = Sigma, a_sigma = a_sigma, b_sigma = b_sigma, verbose = TRUE, maxiter = 1500
  )
  Rsquared <- 1 - sum((dataset_high_signal - out$prediction)^2) / sum((dataset_high_signal - mean(dataset_high_signal))^2)
  Kunique <- length(unique(out$cluster))
  c(seeds, out$lowerbound, Rsquared, Kunique)
}

set.seed(lower_bounds2[which.max(lower_bounds2[, 2]), 1])
fit2 <- EFDMP(
  X = dataset_high_signal, Hl = Hl, L = L, B = B, indexB = indexB, time = time, cc = cc,
  alpha = alpha, Sigma = Sigma, a_sigma = a_sigma, b_sigma = b_sigma, verbose = TRUE, maxiter = 1500
)

# Display some useful quantities
lower_bounds2
length(unique(fit2$cluster))
table(fit2$Fclass)
table(fit2$cluster)
table(clusters, fit2$Fclass)
table(clusters, fit2$cluster)

data_plot <- reshape2::melt(dataset_high_signal)
data_plot$Cluster <- factor(data_plot$Var1)

selected_cluster <- as.numeric(unique(fit2$cluster))
data_plot2 <- fit2$pred_cluster[selected_cluster, ]
rownames(data_plot2) <- selected_cluster
data_plot2 <- reshape2::melt(data_plot2)
data_plot2$Var3 <- fit2$Findex[data_plot2$Var1]
p2 <- ggplot(data = data_plot, aes(x = Var2, y = value, group = as.factor(Var1))) +
  geom_line(alpha = 0.2, col = "gray") +
  theme_bw() +
  scale_colour_brewer(palette = "Set1") +
  ylab("y") +
  xlab("Time") +
  geom_line(data = data_plot2, size = 0.7, aes(x = Var2, y = value, group = as.factor(Var1))) +
  theme(legend.position = "none")
p2

# Mixture splines --------------------------------------------

library(splines)
L <- 1

B <- bs(time, df = 20, intercept = T)
indexB <- matrix(rep(TRUE, 20), 1, 20, byrow = TRUE)

# List of hyperparameters
L <- 1
cc <- 1
alpha <- 1
Sigma <- list(diag(rep(0.5, 20)))
a_sigma <- 1
b_sigma <- 1
Hl <- 80

# Execution of the algorithm ---------------------------------
library(doMC)
library(foreach)
registerDoMC(6)

# Select different seeds
set.seed(555)
seeds <- sample(1:2^15, 50) # Convergence occurs always.

lower_bounds3 <- foreach(seeds = seeds, .combine = rbind) %dopar% {
  set.seed(seeds)
  out <- EFDMP(
    X = dataset_high_signal, Hl = Hl, L = L, B = B, indexB = indexB, time = time, cc = cc,
    alpha = alpha, Sigma = Sigma, a_sigma = a_sigma, b_sigma = b_sigma, verbose = TRUE, maxiter = 1500
  )
  Rsquared <- 1 - sum((dataset_high_signal - out$prediction)^2) / sum((dataset_high_signal - mean(dataset_high_signal))^2)
  Kunique <- length(unique(out$cluster))
  c(seeds, out$lowerbound, Rsquared, Kunique)
}

set.seed(lower_bounds3[which.max(lower_bounds3[, 2]), 1])
fit3 <- EFDMP(
  X = dataset_high_signal, Hl = Hl, L = L, B = B, indexB = indexB, time = time, cc = cc,
  alpha = alpha, Sigma = Sigma, a_sigma = a_sigma, b_sigma = b_sigma, verbose = TRUE, maxiter = 1500
)

# Display some useful quantities
lower_bounds3
length(unique(fit3$cluster))
table(fit3$Fclass)
table(fit3$cluster)
table(clusters, fit3$Fclass)
table(clusters, fit3$cluster)

data_plot <- reshape2::melt(dataset_high_signal)
data_plot$Cluster <- factor(data_plot$Var1)

selected_cluster <- as.numeric(unique(fit3$cluster))
data_plot2 <- fit3$pred_cluster[selected_cluster, ]
rownames(data_plot2) <- selected_cluster
data_plot2 <- reshape2::melt(data_plot2)
data_plot2$Var3 <- fit3$Findex[data_plot2$Var1]
p3 <- ggplot(data = data_plot, aes(x = Var2, y = value, group = as.factor(Var1))) +
  geom_line(alpha = 0.2, col = "gray") +
  theme_bw() +
  scale_colour_brewer(palette = "Set1") +
  ylab("y") +
  xlab("Time") +
  geom_line(data = data_plot2, size = 0.7, aes(x = Var2, y = value, group = as.factor(Var1))) +
  theme(legend.position = "none")
p3

# Final results --------------------------------
tab <- cbind(cbind(summary(apply(dataset_signal - fit$prediction, 1, function(x) sqrt(mean(x^2) - mean(x)^2))),
summary(apply(dataset_signal - fit2$prediction, 1, function(x) sqrt(mean(x^2) - mean(x)^2))),
summary(apply(dataset_signal - fit3$prediction, 1, function(x) sqrt(mean(x^2) - mean(x)^2)))),

cbind(summary(apply(dataset_signal - fit$prediction, 1, function(x) mean(abs(x)))),
summary(apply(dataset_signal - fit2$prediction, 1, function(x) mean(abs(x)))),
summary(apply(dataset_signal - fit3$prediction, 1, function(x) mean(abs(x))))))

fit_kmeans <- kmeans(dataset_high_signal, 16, nstart = 20)

library(mcclust)
mcclust::vi.dist(clusters, fit_kmeans$cluster)
mcclust::vi.dist(clusters, as.numeric(fit$cluster))
mcclust::vi.dist(clusters, as.numeric(fit2$cluster))
mcclust::vi.dist(clusters, as.numeric(fit3$cluster))

ggsave("../sim0.pdf", p0, width=10,height=6)
ggsave("../sim1.pdf", p1, width=10,height=6)
ggsave("../sim2.pdf", p2, width=10,height=6)
ggsave("../sim3.pdf", p3, width=10,height=6)