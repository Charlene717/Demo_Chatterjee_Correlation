# 如果尚未安装 energy 包，请先安装并加载
if(!require('energy')) {install.packages('energy'); library(energy)}



# 定义组合相关系数的函数
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

# 测试独立的随机数据
set.seed(789)
x3 <- rnorm(100)
y3 <- rnorm(100)

# 计算组合相关系数
result <- combined_correlation(x3, y3)
cat("Chatterjee 相关系数（独立数据）：", result$xi_chatterjee, "\n")
cat("距离相关系数（独立数据）：", result$xi_dcor, "\n")
cat("组合相关系数（独立数据）：", result$xi_combined, "\n")
