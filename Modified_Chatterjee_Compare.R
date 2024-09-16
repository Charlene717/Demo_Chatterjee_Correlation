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


# 加载包
library(infotheo)
library(Rfast)

# install.packages("mpmi")
library(mpmi)


# 改良的 Chatterjee 相关系数函数，结合互信息
modified_chatterjee_B <- function(X, Y, alpha = 0.5, method = "knn", k = 10, num_bins = NULL) {
  n <- length(X)
  
  if(length(Y) != n) {
    stop("X 和 Y 的长度必须相同")
  }
  
  # 检查 alpha 是否在 [0, 1] 范围内
  if(alpha < 0 || alpha > 1) {
    stop("参数 alpha 必须在 [0, 1] 范围内")
  }
  
  # 步骤 1：计算原始的 Chatterjee 相关系数
  chatterjee_corr <- function(X, Y) {
    n <- length(X)
    order_X <- order(X)
    Y_ordered <- Y[order_X]
    R <- rank(Y_ordered, ties.method = "average")
    S <- sum(abs(diff(R)))
    xi <- 1 - (3 * S) / (n^2 - 1)
    return(xi)
  }
  
  xi <- chatterjee_corr(X, Y)
  
  # 步骤 2：计算互信息
  if(method == "discrete") {
    # 如果未指定 num_bins，使用 Sturges' 公式
    if(is.null(num_bins)) {
      num_bins <- ceiling(1 + log2(n))
    }
    
    # 离散化数据（指定包名）
    X_discrete <- infotheo::discretize(X, disc = "equalfreq", nbins = num_bins)
    Y_discrete <- infotheo::discretize(Y, disc = "equalfreq", nbins = num_bins)
    
    # 计算互信息
    mi <- infotheo::mutinformation(X_discrete, Y_discrete)
    
    # 计算熵
    hx <- infotheo::entropy(X_discrete)
    hy <- infotheo::entropy(Y_discrete)
  }
  else if(method == "knn") {
    # 使用 mpmi 包的基于 k-NN 的方法
    # 计算互信息
    mi <- knn_mi(X, Y, k = k)
    
    # 计算熵
    hx <- knn_entropy(X, k = k)
    hy <- knn_entropy(Y, k = k)
  }
  else {
    stop("未知的方法。请使用 'discrete' 或 'knn'")
  }
  
  # 计算标准化互信息
  nmi <- 2 * mi / (hx + hy)
  
  # 处理特殊情况
  if(is.na(nmi) || is.infinite(nmi)) {
    nmi <- 0
  }
  
  # 将 NMI 限制在 [0, 1] 范围内
  nmi <- max(0, min(1, nmi))
  
  # 步骤 3：计算结合的相关系数
  xi_combined <- alpha * xi + (1 - alpha) * nmi
  
  # 返回结果
  return(list(
    xi_chatterjee = xi,
    mutual_information = mi,
    normalized_mutual_information = nmi,
    xi_combined = xi_combined
  ))
}


kernel_chatterjee <- function(x, y, kernel = "gaussian", sigma = NULL, degree = 3, coef0 = 1) {
  # 检查输入
  if (length(x) != length(y)) {
    stop("x and y must have the same length.")
  }
  n <- length(x)
  
  # 定义核函数
  kernel_function <- function(u, v) {
    if (kernel == "gaussian") {
      # 如果 sigma 未指定，使用数据的中位数距离
      if (is.null(sigma)) {
        sigma <- median(dist(c(u, v)))
      }
      return(exp(- (outer(u, v, "-")^2) / (2 * sigma^2)))
    } else if (kernel == "linear") {
      return(outer(u, v, "*"))
    } else if (kernel == "polynomial") {
      return((outer(u, v, "*") + coef0)^degree)
    } else {
      stop("Unsupported kernel type.")
    }
  }
  
  # 计算核矩阵
  KX <- kernel_function(x, x)
  KY <- kernel_function(y, y)
  
  # 计算中心化矩阵 H
  H <- diag(n) - matrix(1, n, n) / n
  
  # 中心化核矩阵
  KX_centered <- H %*% KX %*% H
  KY_centered <- H %*% KY %*% H
  
  # 计算迹
  numerator <- sum(KX_centered * KY_centered)  # 等价于 trace(A %*% B)
  denominator <- sqrt(sum(KX_centered * KX_centered) * sum(KY_centered * KY_centered))
  
  # 计算改良的 Chatterjee 相关系数
  xi_kernel <- numerator / denominator
  
  return(xi_kernel)
}

global_chatterjee <- function(x, y, gamma = NULL) {
  n <- length(x)
  if (length(y) != n) {
    stop("x and y must have the same length.")
  }
  
  # 计算秩
  rank_x <- rank(x, ties.method = "average")
  rank_y <- rank(y, ties.method = "average")
  
  # 计算距离矩阵
  dist_x <- as.matrix(dist(x, method = "euclidean"))
  dist_y <- as.matrix(dist(y, method = "euclidean"))
  
  # 计算权重矩阵
  if (is.null(gamma)) {
    gamma <- 1 / median(dist_x[upper.tri(dist_x)])
  }
  w <- exp(-gamma * dist_x)
  
  # 计算秩差的绝对值矩阵
  rank_diff <- abs(outer(rank_y, rank_y, "-"))
  
  # 计算分子和分母
  numerator <- sum(w * rank_diff)
  denominator <- sum(w) * (n - 1)
  
  # 计算改良的 Chatterjee 相关系数
  xi_global <- 1 - numerator / denominator
  
  return(xi_global)
}


improved_chatterjee <- function(x, y, gamma = NULL) {
  n <- length(x)
  if (length(y) != n) {
    stop("x and y must have the same length.")
  }
  
  # 计算秩
  rank_x <- rank(x, ties.method = "average")
  rank_y <- rank(y, ties.method = "average")
  
  # 计算距离
  dist_x <- abs(x - mean(x))
  
  # 计算权重
  if (is.null(gamma)) {
    gamma <- 1 / median(abs(dist_x - median(dist_x)))
  }
  w <- exp(-gamma * dist_x)
  
  # 对权重进行归一化
  w <- w / sum(w)
  
  # 计算加权平均值
  mean_rx <- sum(w * rank_x)
  mean_ry <- sum(w * rank_y)
  
  # 计算加权协方差和方差
  cov_w <- sum(w * (rank_x - mean_rx) * (rank_y - mean_ry))
  var_rx <- sum(w * (rank_x - mean_rx)^2)
  var_ry <- sum(w * (rank_y - mean_ry)^2)
  
  # 计算改良的 Chatterjee 相关系数
  xi_improved <- cov_w / sqrt(var_rx * var_ry)
  
  return(xi_improved)
}

distance_chatterjee <- function(x, y) {
  n <- length(x)
  if (length(y) != n) {
    stop("x and y must have the same length.")
  }
  
  # 计算距离矩阵
  D_X <- as.matrix(dist(x))
  D_Y <- as.matrix(dist(y))
  
  # 双中心化距离矩阵
  J <- diag(n) - matrix(1, n, n) / n
  A <- -0.5 * J %*% (D_X^2) %*% J
  B <- -0.5 * J %*% (D_Y^2) %*% J
  
  # 计算迹
  numerator <- sum(A * B)
  denominator <- sqrt(sum(A * A) * sum(B * B))
  
  # 计算距离相关系数
  xi_distance <- numerator / denominator
  
  return(xi_distance)
}


improved_chatterjee <- function(x, y) {
  n <- length(x)
  if (length(y) != n) {
    stop("x and y must have the same length.")
  }
  
  # 计算秩
  rank_x <- rank(x, ties.method = "average")
  rank_y <- rank(y, ties.method = "average")
  
  # 初始化 L_i
  L <- numeric(n)
  
  for (i in 1:n) {
    # 计算 Y_i 与其他 Y 的绝对差值
    diff_y <- abs(y[i] - y[-i])
    # 找到最近邻的索引
    nn_index <- which.min(diff_y)
    if (nn_index >= i) nn_index <- nn_index + 1  # 调整索引，因为 y[-i]
    
    # 计算 L_i
    L[i] <- abs(rank_x[i] - rank_x[nn_index]) + 1
  }
  
  # 计算改进的 Chatterjee 相关系数
  xi_hat <- 1 - (sum(L - 1)) / (n * (n - 1) / 2)
  
  return(xi_hat)
}


improved_chatterjee <- function(x, y) {
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

#### Calculate correlations ####
# Define a function to calculate and display Pearson, Spearman, and Chatterjee correlations
calculate_correlations <- function(x, y) {
  chatterjee_corr_A <- modified_chatterjee_A(x, y)
  # chatterjee_corr_B <- modified_chatterjee_B(x, y,method = "discrete")$xi_combined
  
  # chatterjee_corr_B <- kernel_chatterjee(x, y, kernel = "gaussian")
  # chatterjee_corr_B <-  global_chatterjee(x, y)
  # chatterjee_corr_B <- improved_chatterjee(x, y)
  # chatterjee_corr_B <- distance_chatterjee(x, y)
  # chatterjee_corr_B <- improved_chatterjee(x, y)
  chatterjee_corr_B <- improved_chatterjee(x, y)
  
  chatterjee_corr <- xicor(x, y)
  
  return(c(chatterjee_corr_A, chatterjee_corr_B, chatterjee_corr))
}

#### Visualization ####
# Define the plotting function
plot_data <- function(x, y, title, correlations) {
  chatterjee_corr_A <- round(correlations[1], 2)
  chatterjee_corr_B <- round(correlations[2], 2)
  chatterjee_corr <- round(correlations[3], 2)
  
  ggplot(data = data.frame(x = x, y = y), aes(x = x, y = y)) +
    geom_point(color = 'blue', size = 2) +
    ggtitle(paste0(title, "\nChatterjee_A: ", chatterjee_corr_A,
                   " | Chatterjee_B: ", chatterjee_corr_B,
                   " | Chatterjee: ", chatterjee_corr)) +
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
