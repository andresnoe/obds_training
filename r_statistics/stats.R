
#Maisha
#Generate a vector of 1000 normally distributed values with mean 10 and standard deviation 5.
#Inspect the output of the summary() function for that vector.
#Compute the mean and standard deviation for those values
#Compute the deciles (i.e. 10 evenly spaced quantiles) for those values.
#Visualise the distribution of those values as a histogram. You can use base  of ggplot2.
#Visualise as vertical lines on the histogram: the mean, median, one standard deviation, and one median absolute deviation.
#Generate a new vector with a lot more values (e.g., one million). Draw again a histogram. How doesthe distribution compare with more data points?

num_vector <- rnorm(1000, mean = 10, sd = 5)
summary(num_vector)
mean_vector <- mean(num_vector)
sd(num_vector)
quantile(num_vector, probs = seq(0, 1, 0.1))

hist(num_vector, breaks = 20)

#ggplot
vector_df <- data.frame(NAME = num_vector)
head(vector_df)
ggplot(data = vector_df, aes(x=NAME)) +
    geom_histogram() +
    geom_vline(xintercept = mean_vector)

#generate a larger vector
large_vector <- rnorm(1e6, mean = 10, sd = 5)
summary(large_vector)
hist(large_vector, breaks = 100)

#-----------------------------------
#-----------------------------------
#-----------------------------------



#------------------Phil-------------------#
# For the standard normal distribution (mean = 0, sd = 1)
# 1. Plot the cumulative distribution function in the range [-5, 5]
# 2. Plot the inverse cumulative distribution function for quantiles in 0.01 increment.
# 3. Plot the density function in the range [-5, 5]
# 4. What is the probability of observing a value greater than 2?
# 5. What is the probability of observing a value between -2 and 2?
# 6. What is the probability of observing a value more extreme than -2 or 2?
#-----------------------------------------#

## /Q1/
# set up the x coordinate
x_index = seq(-5, 5, length.out = 1000)

# generate the required function of the distribution, and input it as column in a tidyverse table (ie tibble)
dst <- tibble(
    x_coordinate = x_index,
    cdf_normal_distribution = pnorm(q = x_index, mean = 0, sd = 1),
)

# visualize using ggplot
dst %>%
    ggplot(aes(x=x_coordinate, y=cdf_normal_distribution)) +
    geom_point()

## /Q2/
# set up the x coordinate for probability
x_prob = seq(0, 1, 0.01)# length.out = 1000)

dst <- tibble(
    x_prob = x_prob,
    inverse_cdf_normal = qnorm(p = x_prob, mean=0, sd=1)
)

dst %>%
    ggplot(aes(x=x_prob, y=inverse_cdf_normal)) +
    geom_point()

## /Q3/
# set up the x coordinate
x_index = seq(-5, 5, length.out = 1000)
# generate the required function of the distribution, and input it as column in a tidyverse table (ie tibble)
dst <- tibble(
    x_coordinate = x_index,
    density_normal_distribution = dnorm(x = x_index, mean = 0, sd = 1)
)
# plot using ggplot
dst %>%
    ggplot(aes(x=x_coordinate, y=density_normal_distribution)) +
    geom_point()

## /Q4/
1 - pnorm(2)

## /Q5/
pnorm(2) - pnorm(-2)

## /Q6/
1 - (pnorm(2) - pnorm(-2))

#-----------------------------------
#-----------------------------------
#-----------------------------------

#Jen - Iris data set, normal distribution t-test and ANOVA

summary(iris)

#visualise distribution

ggplot(iris, aes(Sepal.Length)) +
    geom_histogram(color = "black") +
    facet_wrap(~Species, ncol = 1)

shapiro.test(iris$Sepal.Length)

#compare  between species

species_setosa <- subset(iris, Species == "setosa", select = Sepal.Length)
species_versicolor <- subset(iris, Species == "versicolor", select = Sepal.Length)
species_virginica <- subset(iris, Species == "virginica", select = Sepal.Length)


shapiro.test(species_setosa$Sepal.Length)
shapiro.test(species_versicolor$Sepal.Length)
shapiro.test(species_virginica$Sepal.Length)



anova<- aov(formula = Sepal.Length ~ Species, data = iris)
summary(anova)

levels(iris$Species)

ttest1<- t.test(species_setosa$Sepal.Length, species_versicolor$Sepal.Length)
str(ttest1)

ttest2<- t.test(species_setosa$Sepal.Length, species_virginica$Sepal.Length)
ttest2

ttest3<- t.test(species_versicolor$Sepal.Length, species_virginica$Sepal.Length)
ttest3

#-----------------------------------
#-----------------------------------
#-----------------------------------

Linear Regression Exercise (sophie):

#Fit a linear model to measure the effect of time & diet on weight
ChickWeight
head(ChickWeight)


#first formula is running 2 separate univariate linear model,

formula <- formula(weight ~ Diet + Time)
lm_chick <- lm(formula = formula, data = ChickWeight)
summary(lm_chick)
#coeffi estiamte is an average

ggplot(ChickWeight, aes(x= Time, y = weight, colour = Diet)) +
   geom_point() +
   geom_smooth(method = "lm") +
   geom_abline(slope = 8.7505, intercept = 10.9244)
#abline drawring average across diets

#taking diet and time instead of separately, now we're running multivariate
formula2 <- formula(weight ~ Diet * Time)
lm_chick <- lm(formula = formula2, data = ChickWeight)
summary(lm_chick)

#on avg all gain 6.8g over time. diet2 gives extra 1.7g on top etc etc

ChickWeight$Diet <- relevel(ChickWeight$Diet, "3")
lm_chick <- lm(formula = formula2, data = ChickWeight)
summary(lm_chick)
levels(ChickWeight$Diet) #all diets defined by diet 3 (reference)

ggplot(ChickWeight, aes(x= Time, y = weight, colour = Diet)) +
   geom_point() +
   geom_smooth(method = "lm") +
   geom_abline(slope = 11.4229, intercept = 18.2503) +
   geom_abline(slope = (11.4229 - 4.5811), intercept = (18.2503 + 12.6807))




#-----------------------------------
#-----------------------------------
#-----------------------------------

# Andrés - AT THE END (PASTE ABOVE THIS ^^)
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


## Extras
#-----------------
# over representation analysis - unfinished
#-----------------
sig_de <- padj_values[padj_values<0.05]
go_info <- read.csv('data/go_info.csv')
human_go <- read.csv('data/human_go_bp.csv')

# map go_id to significant genes in tidy form
ora <- lapply(
    names(sig_de),
    function(x) filter(human_go, ensembl_gene_id==x)
    ) %>%
    bind_rows %>%
    group_by(go_id) %>%
    summarise(count=n()) %>%
    arrange(-count) %>%
    left_join(go_info, by=c('go_id'='GOID'))
head(ora, 10)

# go_id      count TERM
# GO:0060337    16 type I interferon signaling pathway
# GO:0045071    14 negative regulation of viral genome replication
# GO:0030593     8 neutrophil chemotaxis
# GO:0048661     8 positive regulation of smooth muscle cell proliferation
# GO:0060333     8 interferon-gamma-mediated signaling pathway
# GO:0070098     8 chemokine-mediated signaling pathway
# GO:0071346     8 cellular response to interferon-gamma
# GO:0032755     7 positive regulation of interleukin-6 production
# GO:0071347     7 cellular response to interleukin-1
# GO:0002548     6 monocyte chemotaxis
