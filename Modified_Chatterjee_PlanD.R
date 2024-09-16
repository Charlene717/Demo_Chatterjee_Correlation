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

set.seed(789)
x3 <- rnorm(100)
y3 <- rnorm(100)

# 计算改良的 Chatterjee 相关系数
xi_improved <- improved_chatterjee(x3, y3)
cat("改良的 Chatterjee 相关系数（独立数据）：", xi_improved, "\n")

