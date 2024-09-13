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

x5 <- seq(-10, 10, length.out = 100)
y5 <- x5^2 + rnorm(100, sd=10) # Moderate Quadratic relationship

x6 <- seq(-10, 10, length.out = 100)
y6 <- x6^2 + rnorm(100, sd=1) # Highly Quadratic relationship

x7 <- seq(-10, 10, length.out = 200)
y7 <- sin(x7) + rnorm(200, sd=0.3) # Moderate Sine wave relationship

x8 <- seq(-10, 10, length.out = 200)
y8 <- sin(x8) + rnorm(200, sd=0.1) # Highly Sine wave relationship

# Store the data in a list
data_list <- list(
  list(x = x1, y = y1, title = "Highly Linear Positive"),
  list(x = x2, y = y2, title = "Moderate Linear Positive"),
  list(x = x3, y = y3, title = "No Correlation"),
  list(x = x4, y = y4, title = "Highly Linear Negative"),
  list(x = x5, y = y5, title = "Moderate Quadratic Relationship"),
  list(x = x6, y = y6, title = "Highly Quadratic Relationship"),
  list(x = x7, y = y7, title = "Moderate Sine Wave Relationship"),
  list(x = x8, y = y8, title = "Highly Sine Wave Relationship")
)


#### Calculate correlations ####
# Define a function to calculate and display Pearson, Spearman, and Chatterjee correlations
calculate_correlations <- function(x, y) {
  pearson_corr <- cor(x, y, method = "pearson")
  spearman_corr <- cor(x, y, method = "spearman")

  chatterjee_corr <- xicor(x, y)
  # chatterjee_corr <- calculateXI(x, y)
  
  return(c(pearson_corr, spearman_corr, chatterjee_corr))
}


# Calculate and display results
for (data in data_list) {
  correlations <- calculate_correlations(data$x, data$y)
  cat("Dataset:", data$title, "\n")
  cat("Pearson Correlation: ", round(correlations[1], 2), "\n")
  cat("Spearman Correlation: ", round(correlations[2], 2), "\n")
  cat("Chatterjee Correlation: ", round(correlations[3], 2), "\n\n")
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

