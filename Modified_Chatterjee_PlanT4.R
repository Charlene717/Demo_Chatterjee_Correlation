improved_chatterjee <- function(x, y, gamma = NULL) {
  n <- length(x)
  if (length(y) != n) {
    stop("x and y must have the same length.")
  }
  
  # 计算秩
  rank_x <- rank(x, ties.method = "average")
  rank_y <- rank(y, ties.method = "average")
  
  # 计算 Y 的秩差绝对值矩阵
  D_Y <- abs(outer(rank_y, rank_y, "-"))
  
  # 计算 X 的距离矩阵
  D_X <- abs(outer(x, x, "-"))
  
  # 计算相似度矩阵 S_X
  if (is.null(gamma)) {
    gamma <- 1 / median(D_X[upper.tri(D_X)])
  }
  S_X <- exp(-gamma * D_X)
  
  # 计算加权平均秩差
  numerator <- sum(S_X * D_Y)
  denominator <- sum(S_X)
  
  bar_D_Y <- numerator / denominator
  
  # 计算改良的 Chatterjee 相关系数
  xi_global <- 1 - (bar_D_Y / (n - 1))
  
  return(xi_global)
}


set.seed(789)
x3 <- rnorm(100)
y3 <- rnorm(100)

# 计算改良的 Chatterjee 相关系数
xi_global <- improved_chatterjee(x3, y3)
cat("改良的 Chatterjee 相关系数（独立数据）：", xi_global, "\n")
