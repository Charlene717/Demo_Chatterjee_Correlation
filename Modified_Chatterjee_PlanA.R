# 改良的 Chatterjee 相关系数函数
modified_chatterjee <- function(X, Y, h = NULL) {
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




# 示例数据
X <- c(3, 1, 4, 2, 5)
Y <- c(2, 5, 1, 4, 3)

# 计算改良的 Chatterjee 相关系数
xi_global <- modified_chatterjee(X, Y)
print(xi_global)

ggplot(data = data.frame(x = X, y = Y), aes(x = X, y = Y)) +
  geom_point(color = 'blue', size = 2)

#################################################################
#
x7 <- seq(0, 1, length.out = 100)
y7 <- 10 * (x7^7) + rnorm(100, sd = 1) # 7th-degree Moderate Non-linear Positive
xi_global_7 <- modified_chatterjee(x7, y7)
print(xi_global_7)

ggplot(data = data.frame(x = x7, y = y7), aes(x = x7, y = y7)) +
  geom_point(color = 'blue', size = 2)

#
x12 <- seq(-10, 10, length.out = 200)
y12 <- sin(x12) + rnorm(200, sd=0.1) # Highly Sine wave relationship

xi_global_12 <- modified_chatterjee(x12, y12)
print(xi_global_12)

ggplot(data = data.frame(x = x12, y = y12), aes(x = x12, y = y12)) +
  geom_point(color = 'blue', size = 2)

#
x11 <- seq(-10, 10, length.out = 200)
y11 <- sin(x11) + rnorm(200, sd=0.3) # Moderate Sine wave relationship

xi_global_11 <- modified_chatterjee(x11, y11)
print(xi_global_11)

ggplot(data = data.frame(x = x11, y = y11), aes(x = x11, y = y11)) +
  geom_point(color = 'blue', size = 2)


#
x1 <- rnorm(100)
y1 <- x1 + rnorm(100, sd=0.1) # Highly linearly correlated

xi_global_1 <- modified_chatterjee(x1, y1)
print(xi_global_1)

ggplot(data = data.frame(x = x1, y = y1), aes(x = x1, y = y1)) +
  geom_point(color = 'blue', size = 2)


#
x2 <- rnorm(100)
y2 <- 0.5 * x2 + rnorm(100, sd=0.2) # Moderately linearly correlated

xi_global_2 <- modified_chatterjee(x2, y2)
print(xi_global_2)

ggplot(data = data.frame(x = x2, y = y2), aes(x = x2, y = y2)) +
  geom_point(color = 'blue', size = 2)
