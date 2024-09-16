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



# 生成二次关系数据
set.seed(123)
n <- 100
x <- runif(n, -1, 1)
y <- x^2 + rnorm(n, sd = 0.1)

# 计算原始 Chatterjee 相关系数
xi_original <- chatterjee_corr(x, y)
cat("原始的 Chatterjee 相关系数：", xi_original, "\n")

# 计算改良的全局 Chatterjee 相关系数
xi_global <- global_chatterjee(x, y)
cat("改良的全局 Chatterjee 相关系数：", xi_global, "\n")





# 生成正弦波关系数据
set.seed(456)
n <- 100
x <- runif(n, 0, 2 * pi)
y <- sin(x) + rnorm(n, sd = 0.1)

# 计算原始 Chatterjee 相关系数
xi_original <- chatterjee_corr(x, y)
cat("原始的 Chatterjee 相关系数：", xi_original, "\n")

# 计算改良的全局 Chatterjee 相关系数
xi_global <- global_chatterjee(x, y)
cat("改良的全局 Chatterjee 相关系数：", xi_global, "\n")



