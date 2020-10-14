# # R FOR STATISTICS

# Generate and visualise a distribution!
# Generate a vector of 1000 normally distributed values with mean 10 and standard deviation 5.
num_vector <- rnorm(n=1000,mean=10, sd=5)
num_vector

# Inspect the output of the function for that vector
summary(num_vector)

# Compute the mean and standard deviation for those values
mean(num_vector)
sd(num_vector)

# Compute the deciles (i.e. 10 evenly spaced quantiles) for those values.
quantile(num_vector, probs = seq(from = 0, to = 1, by = 0.1))

# Visualise the distribution of those values as a histogram. You can use base R or ggplot2.
hist(num_vector, breaks = 50)

num_df <- data.frame(NAME = num_vector)
head(num_df)
library(tidyverse)
install.packages("ggthemes")
library(ggthemes)

num_df %>%
    ggplot(aes(x=NAME)) +
    geom_histogram() +
    theme_economist() +
    labs(title = "1000 normally distributed numbers", x = "Number")

# Visualise as vertical lines on the histogram: the mean, median, one standard deviation, and one median absolute deviation.
mean_vector <- mean(num_df$NAME)

num_df %>%
    ggplot(aes(x=NAME)) +
    geom_histogram() +
    theme_economist() +
    labs(title = "1000 normally distributed numbers", x = "Number") +
    geom_vline(xintercept = mean_vector)

# Generate a new vector with a lot more values (e.g., one million). Draw again a histogram. How does the distribution compare with more data points?
num_vector2 <- rnorm(n=1e6,mean=10, sd=5)
num_vector2
hist(num_vector2, breaks = 1000)


# Query distributions and probabilities
# For the standard normal distribution N(u = 0, variance = 1)


# Plot the cumulative distribution function in the range .
# Plot the inverse cumulative distribution function for quantiles in 0.01 increment. Plot the density function in the range .
# What is the probability of observing a value greater than 2?
#     What is the probability of observing a value between -2 and 2?
#     What is the probability of observing a value more extreme than -2 or 2?
