# 檢查並安裝所需的套件
if(!require('ggplot2')) {install.packages('ggplot2'); library(ggplot2)}
if(!require('GGally')) {install.packages('GGally'); library(GGally)}
if(!require('XICOR')) {install.packages('XICOR'); library(XICOR)} # Chatterjee Correlation
if(!require('gridExtra')) {install.packages('gridExtra'); library(gridExtra)}


# 生成模擬數據
set.seed(123) # 保持結果可重複
x1 <- rnorm(100)
y1 <- x1 + rnorm(100, sd=0.1) # 高度線性相關

x2 <- rnorm(100)
y2 <- 0.5 * x2 + rnorm(100, sd=0.5) # 中等線性相關

x3 <- rnorm(100)
y3 <- rnorm(100) # 幾乎無相關

x4 <- rnorm(100)
y4 <- -x4 + rnorm(100, sd=0.1) # 反向相關

x5 <- seq(-10, 10, length.out = 100)
y5 <- x5^2 + rnorm(100, sd=1) # 二次曲線

x6 <- seq(-10, 10, length.out = 100)
y6 <- sin(x6) + rnorm(100, sd=0.5) # 正弦波

# 將數據存入列表
data_list <- list(
  list(x = x1, y = y1, title = "Highly Linear Positive"),
  list(x = x2, y = y2, title = "Moderate Linear Positive"),
  list(x = x3, y = y3, title = "No Correlation"),
  list(x = x4, y = y4, title = "Highly Linear Negative"),
  list(x = x5, y = y5, title = "Quadratic Relationship"),
  list(x = x6, y = y6, title = "Sine Wave Relationship")
)

# 定義函數來計算並顯示 Pearson, Spearman 和 Chatterjee correlation
calculate_correlations <- function(x, y) {
  pearson_corr <- cor(x, y, method = "pearson")
  spearman_corr <- cor(x, y, method = "spearman")

  chatterjee_corr <- xicor(x, y)
  # chatterjee_corr <- calculateXI(x, y)
  
  return(c(pearson_corr, spearman_corr, chatterjee_corr))
}


# 計算並顯示結果
for (data in data_list) {
  correlations <- calculate_correlations(data$x, data$y)
  cat("Dataset:", data$title, "\n")
  cat("Pearson Correlation: ", round(correlations[1], 2), "\n")
  cat("Spearman Correlation: ", round(correlations[2], 2), "\n")
  cat("Chatterjee Correlation: ", round(correlations[3], 2), "\n\n")
}


# 定義繪圖函數
plot_data <- function(x, y, title, correlations) {
  pearson_corr <- round(correlations[1], 2)
  spearman_corr <- round(correlations[2], 2)
  chatterjee_corr <- round(correlations[3], 2)
  
  ggplot(data = data.frame(x = x, y = y), aes(x = x, y = y)) +
    geom_point(color = 'blue', size = 2) +
    ggtitle(paste0(title, "\nPearson: ", pearson_corr,
                   " | Spearman: ", spearman_corr,
                   " | Chatterjee: ", chatterjee_corr)) +
    theme_minimal()
}

# 繪製每個數據集的圖表並顯示相關性
plot_list <- list()
for (data in data_list) {
  correlations <- calculate_correlations(data$x, data$y)
  plot_list[[data$title]] <- plot_data(data$x, data$y, data$title, correlations)
}

# 組合所有圖表
grid.arrange(grobs = plot_list, ncol = 3)
