# 如果尚未安装 energy 包，请先安装并加载
if(!require('energy')) {install.packages('energy'); library(energy)}



# 定义组合相关系数的函数
combined_correlation <- function(x, y, alpha = 0.5, threshold = 0.8, k = 10, method = "spearman") {
  n <- length(x)
  if (length(y) != n) {
    stop("x and y must have the same length.")
  }
  
  # 計算原始的 Chatterjee 相關系數
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
  
  # 計算選定的相關系數
  if (method == "pearson") {
    xi_corr <- abs(cor(x, y, method = "pearson"))
  } else if (method == "spearman") {
    xi_corr <- abs(cor(x, y, method = "spearman"))
  } else if (method == "dcor") {
    if (!require('energy')) { install.packages('energy'); library(energy) }
    xi_corr <- dcor(x, y)
  } else {
    stop("Invalid method. Use 'pearson', 'spearman', or 'dcor'.")
  }
  
  # 增強權重的函數，當超過閾值時使用二次指數增長
  weight_adjustment <- function(corr_value, threshold, k) {
    if (corr_value > threshold) {
      return(exp(k * (corr_value - threshold)^2))
    } else {
      return(1)
    }
  }
  
  # 根據 Chatterjee 和選定的相關系數調整權重
  alpha_adjusted_chatterjee <- alpha * weight_adjustment(xi_chatterjee, threshold, k)
  alpha_adjusted_corr <- (1 - alpha) * weight_adjustment(xi_corr, threshold, k)
  
  # 標準化權重使其和為1
  total_weight <- alpha_adjusted_chatterjee + alpha_adjusted_corr
  alpha_adjusted_chatterjee <- alpha_adjusted_chatterjee / total_weight
  alpha_adjusted_corr <- alpha_adjusted_corr / total_weight
  
  # 組合相關系數
  xi_combined <- alpha_adjusted_chatterjee * xi_chatterjee + alpha_adjusted_corr * xi_corr
  
  return(list(
    xi_chatterjee = xi_chatterjee,
    xi_corr = xi_corr,
    alpha_adjusted_chatterjee = alpha_adjusted_chatterjee,
    alpha_adjusted_corr = alpha_adjusted_corr,
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
