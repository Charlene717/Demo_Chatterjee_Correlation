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


# 生成独立的随机数据
set.seed(789)
x3 <- rnorm(100)
y3 <- rnorm(100)

# 计算原始 Chatterjee 相关系数
xi_original <- chatterjee_corr(x3, y3)
cat("原始的 Chatterjee 相关系数：", xi_original, "\n")

# 计算改良的 Chatterjee 相关系数
xi_improved <- improved_chatterjee(x3, y3)
cat("改良的 Chatterjee 相关系数：", xi_improved, "\n")
