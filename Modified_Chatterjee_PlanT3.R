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


set.seed(789)
x3 <- rnorm(100)
y3 <- rnorm(100)

# 计算改进的 Chatterjee 相关系数
xi_improved <- improved_chatterjee(x3, y3)
cat("改进的 Chatterjee 相关系数（独立数据）：", xi_improved, "\n")
