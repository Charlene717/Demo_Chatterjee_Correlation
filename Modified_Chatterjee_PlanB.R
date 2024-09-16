# 安装必要的包
# install.packages(c("infotheo", "Rfast"))

# 加载包
library(infotheo)
library(Rfast)

# install.packages("mpmi")
library(mpmi)


# 改良的 Chatterjee 相关系数函数，结合互信息
modified_chatterjee_mi <- function(X, Y, alpha = 0.5, method = "knn", k = 10, num_bins = NULL) {
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




# 安装并加载必要的包
# install.packages("mpmi")
# install.packages("infotheo")
library(mpmi)
library(infotheo)

# 示例数据
set.seed(123)
n <- 100
X <- rnorm(n)
Y <- sin(X) + rnorm(n, sd = 0.1)

# 使用 k-NN 方法计算改良的 Chatterjee 相关系数
result_knn <- modified_chatterjee_mi(X, Y, alpha = 0.5, method = "knn", k = 10)
print("使用 k-NN 方法的结果：")
print(result_knn)

# 使用离散化方法计算改良的 Chatterjee 相关系数
result_discrete <- modified_chatterjee_mi(X, Y, alpha = 0.5, method = "discrete", num_bins = 10)
print("使用离散化方法的结果：")
print(result_discrete)


