##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)

#### Load Packages ####
if(!require('ggplot2')) {install.packages('ggplot2'); library(ggplot2)}
if(!require('GGally')) {install.packages('GGally'); library(GGally)}
if(!require('XICOR')) {install.packages('XICOR'); library(XICOR)} # Chatterjee Correlation
if(!require('gridExtra')) {install.packages('gridExtra'); library(gridExtra)}

#### Generate simulated data ####
set.seed(123) # Keep results reproducible
x1 <- rnorm(100)
y1 <- x1 + rnorm(100, sd=0.1) # Highly linearly correlated

x2 <- rnorm(100)
y2 <- 0.5 * x2 + rnorm(100, sd=0.2) # Moderately linearly correlated

x3 <- rnorm(100)
y3 <- rnorm(100) # Almost no correlation

x4 <- rnorm(100)
y4 <- -0.5 *x4 + rnorm(100, sd=0.1) # Inversely correlated


x5 <- seq(0, 1, length.out = 100)
y5 <- 10 * (x5^11) + rnorm(100, sd = 0.1) # 11th-degree Highly Non-Linear Positive

x6 <- seq(0, 1, length.out = 100)
y6 <- 10 * (x6^7) + rnorm(100, sd = 0.1) # 7th-degree Highly Non-Linear Positive

x7 <- seq(0, 1, length.out = 100)
y7 <- 10 * (x7^7) + rnorm(100, sd = 1) # 7th-degree Moderate Non-linear Positive

x8 <- seq(0, 1, length.out = 100)
y8 <- -10 * (x8^7) + rnorm(100, sd = 0.5)  # 7th-degree Highly Non-linear Negative

x9 <- seq(-10, 10, length.out = 100)
y9 <- x9^2 + rnorm(100, sd=10) # Moderate Quadratic relationship

x10 <- seq(-10, 10, length.out = 100)
y10 <- x10^2 + rnorm(100, sd=1) # Highly Quadratic relationship

x11 <- seq(-10, 10, length.out = 200)
y11 <- sin(x11) + rnorm(200, sd=0.3) # Moderate Sine wave relationship

x12 <- seq(-10, 10, length.out = 200)
y12 <- sin(x12) + rnorm(200, sd=0.1) # Highly Sine wave relationship

# Store the data in a list
data_list <- list(
  list(x = x1, y = y1, title = "Highly Linear Positive"),
  list(x = x2, y = y2, title = "Moderate Linear Positive"),
  list(x = x3, y = y3, title = "No Correlation"),
  list(x = x4, y = y4, title = "Highly Linear Negative"),
  list(x = x5, y = y5, title = "11th-degree Highly Non-Linear Positive"),
  list(x = x6, y = y6, title = "7th-degree Highly Non-linear Positive"),
  list(x = x7, y = y7, title = "7th-degree Moderate Non-linear Positive"),
  list(x = x8, y = y8, title = "7th-degree Highly Non-linear Negative"),
  list(x = x9, y = y9, title = "Moderate Quadratic Relationship"),
  list(x = x10, y = y10, title = "Highly Quadratic Relationship"),
  list(x = x11, y = y11, title = "Moderate Sine Wave Relationship"),
  list(x = x12, y = y12, title = "Highly Sine Wave Relationship")
)



#### Calculate correlations ####
# Define a function to calculate and display Pearson, Spearman, and Chatterjee correlations
calculate_correlations <- function(x, y) {
  pearson_corr <- cor(x, y, method = "pearson")
  spearman_corr <- cor(x, y, method = "spearman")
  
  chatterjee_corr <- xicor(x, y)
  
  return(c(pearson_corr, spearman_corr, chatterjee_corr))
}

#### Visualization ####
# Define the plotting function
plot_data <- function(x, y, title, correlations) {
  pearson_corr <- round(correlations[1], 2)
  spearman_corr <- round(correlations[2], 2)
  chatterjee_corr <- round(correlations[3], 2)
  
  ggplot(data = data.frame(x = x, y = y), aes(x = x, y = y)) +
    geom_point(color = 'blue', size = 2) +
    ggtitle(paste0(title, "\nPearson: ", pearson_corr,
                   " | Spearman: ", spearman_corr,
                   " | Chatterjee: ", chatterjee_corr)) +
    theme_minimal() +
    
    # Set title, axis labels, and axis tick size
    theme(
      plot.title = element_text(size = 11),            # Title font size
      axis.title.x = element_text(size = 14),          # X-axis title size
      axis.title.y = element_text(size = 14),          # Y-axis title size
      axis.text.x = element_text(size = 12),           # X-axis tick text size
      axis.text.y = element_text(size = 12),           # Y-axis tick text size
      
      # Add a black border to the plot
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    labs(x = "X", y = "Y")  # Add labels for X and Y axes
}

# Plot each dataset and show correlations
plot_list <- list()
for (data in data_list) {
  correlations <- calculate_correlations(data$x, data$y)
  plot_list[[data$title]] <- plot_data(data$x, data$y, data$title, correlations)
}

# Combine all plots
grid.arrange(grobs = plot_list, ncol = 4)
