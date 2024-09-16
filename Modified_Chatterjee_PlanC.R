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



# 生成示例数据
set.seed(123)
n <- 100
x <- runif(n)
y <- sin(2 * pi * x) + rnorm(n, sd = 0.1)

# 计算基于核的 Chatterjee 相关系数
xi_kernel <- kernel_chatterjee(x, y, kernel = "gaussian")
cat("基于核的 Chatterjee 相关系数（高斯核）：", xi_kernel, "\n")

# 比较原始的 Chatterjee 相关系数
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

xi_original <- chatterjee_corr(x, y)
cat("原始的 Chatterjee 相关系数：", xi_original, "\n")

