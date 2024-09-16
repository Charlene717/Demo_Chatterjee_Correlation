##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)

#### Load Packages ####
if(!require('ggplot2')) {install.packages('ggplot2'); library(ggplot2)}
if(!require('GGally')) {install.packages('GGally'); library(GGally)}
if(!require('XICOR')) {install.packages('XICOR'); library(XICOR)} # Chatterjee Correlation
if(!require('gridExtra')) {install.packages('gridExtra'); library(gridExtra)}

#### Generate simulated data ####
set.seed(123) # Keep results reproducible
x1 <- rnorm(100)
y1 <- x1 + rnorm(100, sd=0.1) # Highly linearly correlated

x2 <- rnorm(100)
y2 <- 0.5 * x2 + rnorm(100, sd=0.2) # Moderately linearly correlated

x3 <- rnorm(100)
y3 <- rnorm(100) # Almost no correlation

x4 <- rnorm(100)
y4 <- -0.5 *x4 + rnorm(100, sd=0.1) # Inversely correlated


x5 <- seq(0, 1, length.out = 100)
y5 <- 10 * (x5^11) + rnorm(100, sd = 0.1) # 11th-degree Highly Non-Linear Positive

x6 <- seq(0, 1, length.out = 100)
y6 <- 10 * (x6^7) + rnorm(100, sd = 0.1) # 7th-degree Highly Non-Linear Positive

x7 <- seq(0, 1, length.out = 100)
y7 <- 10 * (x7^7) + rnorm(100, sd = 1) # 7th-degree Moderate Non-linear Positive

x8 <- seq(0, 1, length.out = 100)
y8 <- -10 * (x8^7) + rnorm(100, sd = 0.5)  # 7th-degree Highly Non-linear Negative

x9 <- seq(-10, 10, length.out = 100)
y9 <- x9^2 + rnorm(100, sd=10) # Moderate Quadratic relationship

x10 <- seq(-10, 10, length.out = 100)
y10 <- x10^2 + rnorm(100, sd=1) # Highly Quadratic relationship

x11 <- seq(-10, 10, length.out = 200)
y11 <- sin(x11) + rnorm(200, sd=0.3) # Moderate Sine wave relationship

x12 <- seq(-10, 10, length.out = 200)
y12 <- sin(x12) + rnorm(200, sd=0.1) # Highly Sine wave relationship

# Store the data in a list
data_list <- list(
  list(x = x1, y = y1, title = "Highly Linear Positive"),
  list(x = x2, y = y2, title = "Moderate Linear Positive"),
  list(x = x3, y = y3, title = "No Correlation"),
  list(x = x4, y = y4, title = "Highly Linear Negative"),
  list(x = x5, y = y5, title = "11th-degree Highly Non-Linear Positive"),
  list(x = x6, y = y6, title = "7th-degree Highly Non-linear Positive"),
  list(x = x7, y = y7, title = "7th-degree Moderate Non-linear Positive"),
  list(x = x8, y = y8, title = "7th-degree Highly Non-linear Negative"),
  list(x = x9, y = y9, title = "Moderate Quadratic Relationship"),
  list(x = x10, y = y10, title = "Highly Quadratic Relationship"),
  list(x = x11, y = y11, title = "Moderate Sine Wave Relationship"),
  list(x = x12, y = y12, title = "Highly Sine Wave Relationship")
)

#### Modified_Chatterjee ####
modified_chatterjee_A <- function(X, Y, h = NULL) {
  n <- length(X)
  
  if(length(Y) != n) {
    stop("X 和 Y 的长度必须相同")
  }
  
  # 步骤 1：计算 Y 的秩
  R <- rank(Y, ties.method = "average")
  
  # 步骤 2：计算权重矩阵 W
  # 如果未指定 h，使用默认的带宽（Silverman's rule of thumb）
  if(is.null(h)) {
    h <- 1.06 * sd(X) * n^(-1/5)
  }
  
  # 创建距离矩阵
  X_diff <- outer(X, X, "-")
  
  # 计算权重矩阵 W，使用高斯核函数
  W <- exp(- (X_diff^2) / (2 * h^2))
  
  # 步骤 3：计算加权秩差的绝对值之和
  R_diff <- outer(R, R, "-")
  S <- sum(W * abs(R_diff))
  
  # 步骤 4：计算权重的总和
  W_sum <- sum(W)
  
  # 步骤 5：计算改良的 Chatterjee 相关系数
  D_max <- n - 1
  xi_global <- 1 - S / (D_max * W_sum)
  
  return(xi_global)
}


modified_chatterjee_D <- function(x, y) {
  n <- length(x)
  if (length(y) != n) {
    stop("x and y must have the same length.")
  }
  
  # 计算秩
  rank_x <- rank(x)
  rank_y <- rank(y)
  
  # 计算秩的距离矩阵
  A <- as.matrix(dist(rank_x, method = "euclidean"))
  B <- as.matrix(dist(rank_y, method = "euclidean"))
  
  # 中心化距离矩阵
  A_row_mean <- rowMeans(A)
  A_col_mean <- colMeans(A)
  A_mean <- mean(A)
  A_centered <- A - outer(A_row_mean, rep(1, n)) - outer(rep(1, n), A_col_mean) + A_mean
  
  B_row_mean <- rowMeans(B)
  B_col_mean <- colMeans(B)
  B_mean <- mean(B)
  B_centered <- B - outer(B_row_mean, rep(1, n)) - outer(rep(1, n), B_col_mean) + B_mean
  
  # 计算距离协方差和方差
  dCov2 <- sum(A_centered * B_centered) / (n^2)
  dVarX2 <- sum(A_centered^2) / (n^2)
  dVarY2 <- sum(B_centered^2) / (n^2)
  
  # 计算改良的 Chatterjee 相关系数
  xi_improved <- sqrt(dCov2) / sqrt(sqrt(dVarX2) * sqrt(dVarY2))
  
  return(xi_improved)
}


combined_correlation <- function(x, y, alpha = 0.5) {
  n <- length(x)
  if (length(y) != n) {
    stop("x and y must have the same length.")
  }
  
  # 计算原始的 Chatterjee 相关系数
  chatterjee_corr <- function(x, y) {
    n <- length(x)
    rank_x <- rank(x, ties.method = "average")
    rank_y <- rank(y, ties.method = "average")
    order_x <- order(rank_x)
    rank_y_ordered <- rank_y[order_x]
    S <- sum(abs(diff(rank_y_ordered)))
    xi_n <- 1 - (3 * S) / (n^2 - 1)
    return(xi_n)
  }
  xi_chatterjee <- chatterjee_corr(x, y)
  
  # 计算距离相关系数
  if(!require('energy')) {install.packages('energy'); library(energy)}
  xi_dcor <- dcor(x, y)
  
  # 组合相关系数
  xi_combined <- alpha * xi_chatterjee + (1 - alpha) * xi_dcor
  
  return(list(
    xi_chatterjee = xi_chatterjee,
    xi_dcor = xi_dcor,
    xi_combined = xi_combined
  ))
}


#### Calculate correlations ####
# Define a function to calculate and display Pearson, Spearman, and Chatterjee correlations
calculate_correlations <- function(x, y) {

  # chatterjee_corr_A <- modified_chatterjee_A(x, y)
  chatterjee_corr_D <- modified_chatterjee_D(x, y)
  chatterjee_corr_E <- combined_correlation(x, y)$xi_combined
  chatterjee_corr <- xicor(x, y)
  
  pearson_corr <- cor(x, y, method = "pearson")
  spearman_corr <- cor(x, y, method = "spearman")
  
  return(c(chatterjee_corr_D, chatterjee_corr_E, chatterjee_corr, pearson_corr, spearman_corr))
}

#### Visualization ####
# Define the plotting function
plot_data <- function(x, y, title, correlations) {

  chatterjee_corr_D <- round(correlations[1], 2)
  chatterjee_corr_E <- round(correlations[2], 2)
  chatterjee_corr <- round(correlations[3], 2)
  pearson_corr <- round(correlations[4], 2)
  spearman_corr <- round(correlations[5], 2)
  
  ggplot(data = data.frame(x = x, y = y), aes(x = x, y = y)) +
    geom_point(color = 'blue', size = 2) +
    ggtitle(paste0(title, "\nChatterjee_D: ", chatterjee_corr_D,
                   " | Chatterjee_E: ", chatterjee_corr_E,
                   " | Chatterjee: ", chatterjee_corr,
                   "\nPearson: ", pearson_corr,
                   " | Spearman: ", spearman_corr)) +
    theme_minimal() +
    
    # Set title, axis labels, and axis tick size
    theme(
      plot.title = element_text(size = 9),            # Title font size
      axis.title.x = element_text(size = 14),          # X-axis title size
      axis.title.y = element_text(size = 14),          # Y-axis title size
      axis.text.x = element_text(size = 12),           # X-axis tick text size
      axis.text.y = element_text(size = 12),           # Y-axis tick text size
      
      # Add a black border to the plot
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    labs(x = "X", y = "Y")  # Add labels for X and Y axes
}

# Plot each dataset and show correlations
plot_list <- list()
for (data in data_list) {
  correlations <- calculate_correlations(data$x, data$y)
  plot_list[[data$title]] <- plot_data(data$x, data$y, data$title, correlations)
}

# Combine all plots
grid.arrange(grobs = plot_list, ncol = 4)
