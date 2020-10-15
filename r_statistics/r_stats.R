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
# Plot the cumulative distribution function in the range [-5,5].
x_index = seq(-5,5,0.001) # set up x coordinates
dst <- tibble(
    x_coordinate = x_index,
    cd_normal_distribution = pnorm(q = seq(-5,5,0.001), mean=0, sd=1)
)
plot(dst)

dst %>%
    ggplot(aes(x=x_coordinate, y=cd_normal_distribution)) +
    geom_point() +
    theme_economist()

# Plot the inverse cumulative distribution function for quantiles in 0.01 increment.
x_prob = seq(0,1, 0.01) # set up x coordinates
dst <- tibble(
    x_coordinate = x_prob,
    inverse_cd_normal_distribution = qnorm(p = x_prob, mean=0, sd=1)
)
plot(dst)

dst %>%
    ggplot(aes(x=x_prob, y=inverse_cd_normal_distribution)) +
    geom_point() +
    theme_economist()

# Plot the density function in the range [-5,5].
x_index = seq(-5,5,0.001) # set up x coordinates
dst <- tibble(
    x_coordinate = x_index,
    density_normal_distribution = dnorm(x = x_index, mean=0, sd=1),
    random_data = rnorm(n=x_index)
)
plot(dst)

dst %>%
    ggplot(aes(x=x_index, y=density_normal_distribution)) +
    geom_point() +
    theme_economist()
#  What is the probability of observing a value greater than 2?
1 - pnorm(2)

#     What is the probability of observing a value between -2 and 2?
prob2 <- pnorm(q = c(2))
probneg2 <- pnorm(q = c(-2))
prob_between_twos <- prob2-probneg2
prob_between_twos

#     What is the probability of observing a value more extreme than -2 or 2?
1 - prob_between_twos

# We didn't do this
# Compute an Empirical Cumulative Distribution Function
# Use the ecdf() function to compute the empirical cumulative distribution function for the variable Sepal.Lenght in the iris data set.
# Use the plot() function to visualise the empirical cumulative distribution function.
# Use the knots() function on the ecdf() output and compare this with the list of unique values for the variable Sepal.Length


# Statistical tests
# The data set gives the measurements in centimeters of the variables sepal
# length and width and petal length and width, respectively, for 50 flowers from each of 3 species of iris

# Use the summary() function to view some information about each column.
summary(iris)

# Visualise the distribution of Sepal.Length stratified by species
iris %>%
    ggplot(aes(x=Sepal.Length)) +
    geom_histogram(colour = "black") +
    facet_wrap( ~ Species, ncol = 1)

# Is Sepal.Length normally distributed? Overall? Within each species?
# Use Shapiro-Wilk
shapiro.test(iris$Sepal.Length)

for (species_name in levels(iris$Species)){
    print(species_name)
    print(shapiro.test(iris$Sepal.Length[iris$Species==species_name]))
}
shapiro.test(iris$Sepal.Length[iris$Species=='setosa'])
shapiro.test(iris$Sepal.Length[iris$Species=='versicolor'])
shapiro.test(iris$Sepal.Length[iris$Species=='virginica'])

iris %>%
    group_by(Species) %>%

    summarise(shapiro=shapiro.test(Sepal.Length)$p.value)

# Is there a significant variation of Sepal.Length between various species?
anova <- aov (formula = Sepal.Length ~ Species, data = iris)
summary(anova)

test1 <- t.test(iris$Sepal.Length[iris$Species=='setosa'], iris$Sepal.Length[iris$Species=='versicolor'])
test1

test2 <- t.test(iris$Sepal.Length[iris$Species=='virginica'], iris$Sepal.Length[iris$Species=='versicolor'])
test2

test3 <- t.test(iris$Sepal.Length[iris$Species=='setosa'], iris$Sepal.Length[iris$Species=='virginica'])
test3


#  The ChickWeight data set measures the impact of different diets on the early growth of chicks.
# Fit a linear model to measure the effect of Time and Diet on weight in the ChickWeight data set.
summary(ChickWeight)
formula <- formula(weight ~ Diet + Time)
fit <- lm(formula = formula, data = ChickWeight)
summary(fit)
#Coefficient tells us the average difference between


#Let's plot the data
ChickWeight %>%
    ggplot(aes(x=Time, y=weight, colour=Diet))+
    geom_point()+
    geom_smooth(method='lm')

# Interaction effect
formula2 <- formula(weight ~ Diet * Time)
fit2 <- lm(formula = formula2, data = ChickWeight)
summary(fit2)

ChickWeight$Diet <- relevel(ChickWeight$Diet, "3") # diet 3 becomes reference
fit2 <- lm(formula = formula2, data = ChickWeight)
summary(fit2)

levels(ChickWeight$Diet)
# Which diet leads to the fastest increase in body weight?
#Diet3


# Multiple test correction

log_counts <- read.csv("data/logcounts.csv", row.names = 1)
# View(log_counts)
cell_md <- read.csv("data/cell_metadata.csv",row.names=1)
# View(cell_md)

#
gene1 <- data.frame("log_count" = as.numeric(log_counts[1,]), infection = cell_md$Infection)
test_result<- t.test(log_count~infection, gene1)
str(test_result)
test_result[['p.value']]

diff_exp <- function(gene_index, matrix, groups){
    gene_row <- data.frame("log_count" = as.numeric(matrix[gene_index,]), infection = groups)
    test_result<- t.test(log_count~infection, gene_row)
    return(test_result[['p.value']])
}

diff_exp(3,log_counts,cell_md$Infection)
# apply(log_counts, 1, diff_exp, matrix = log_counts, groups = cell_md$Infection)
p_values <- vapply(seq(1,nrow(log_counts)), diff_exp, numeric(1), matrix = log_counts, groups = cell_md$Infection)
names(p_values) <- rownames(log_counts)
p_values

padj_values <- p.adjust(p_values, method = "holm")
hist(padj_values)

head(sort(padj_values))

#OVER-REPRESENTATION ANALYSIS

human_go_bp <- read.csv("data/human_go_bp.csv")
# View(log_counts)
go_info <- read.csv("data/go_info.csv")
# View(cell_md)

#Make list
#Take first pathway
#Turn first pathway's code into vapply
#Do vapply
## For each, build table with Fisher's test, how many belong in the pathway or not (over-rep)
