character_vector <- c("apple", "table", "penny", "cherry", "door")
names(character_vector) <- c("green", "brown", "gold", "red", "white")
character_vector

character_vector["green"]
character_vector["gold"]
character_vector[length(character_vector)]

character_vector

num <- c(6,3,2,18)
bool <- c(TRUE,FALSE,TRUE,TRUE)
sum(num[bool])
num*bool

num_2 <- c(3,4,12)

num+num_2

quality <- factor(c('low','high', 'medium','high', 'low'),
                    levels = c('low', 'medium', 'high'),
                    ordered = TRUE)

mylist <- list(character_vector,bool,quality,num)
mylist
subset_list <- mylist[c(2,4)]
mylist[[3]][[2]]

names(mylist) <- c("objects","answers", "heat", "numbers")
mylist
mylist$numbers


numeric_vector <- 1:10
lap_output <- lapply(numeric_vector,function(x) {x^2})
sap_output <- sapply(numeric_vector,function(x) {x^2})
lap_output
sap_output

num1 <- c(1,2,3,4,5)
num2 <- c(6,7,8,9,10)
log1 <- c(TRUE,FALSE,TRUE,FALSE)
log2 <- c(FALSE,TRUE,FALSE,TRUE)

mylist2 <- list(num1,num2,log1,log2)
class(mylist2)
mylist2

lapply(mylist2,sum)
sapply(mylist2,sum)


rep(num1,each=3)
sapply_list <- sapply(mylist2,function(x){rep(x, each=3)})
sapply_unlist <- unlist(sapply_list)
sapply_unlist

str(sapply_list)
summary(sapply_list)


mymatrix <- matrix(seq(from=2,to=100,by=2),nrow=5,byrow=TRUE)
mymatrix

sum_squares <- apply(mymatrix,1,function(x){sum(x^2)})
sum_squares

min_max <- apply(mymatrix,2,function(x){c(min(x),max(x))})
min_max
min(mymatrix[,1])

mymatrix2 <- matrix(seq(from=1,to=60,by=1),nrow=10, ncol=6,byrow=TRUE)
mymatrix2
transposed <- t(mymatrix2)
dim(mymatrix)
dim(transposed)

joined_matrix <- rbind(transposed,mymatrix)
dim(joined_matrix)
df <- as.data.frame(joined_matrix)
df
class(df)
matrix_as_df_as_list <- list(df)
class(df_as_list)

getwd()
setwd("/Users/andresnoe")
coding_gene_region <- read.table("coding_gene_region.bed", header=FALSE)
str(coding_gene_region)
head(coding_gene_region)
tail(coding_gene_region)
colnames(coding_gene_region) <- c("chr","start","end","ensemblID","score","strand")
head(coding_gene_region)
class(coding_gene_region)


coding_gene_region$length <- coding_gene_region$end-coding_gene_region$start
ordered_coding_gene_region <- coding_gene_region[order(coding_gene_region$length, decreasing = TRUE),]
ordered_coding_gene_region

ordered_coding_gene_region[30,3]

ordered_coding_gene_region[,2]
ordered_coding_gene_region$start

max_chrom <- ordered_coding_gene_region$chr[ordered_coding_gene_region$length==max(ordered_coding_gene_region$length)]
max_chrom

tail(ordered_coding_gene_region, 50)
head(ordered_coding_gene_region[ordered_coding_gene_region$length>=100001 & ordered_coding_gene_region$length<=200000,])
length_regions <- ordered_coding_gene_region[ordered_coding_gene_region$length>=100001 & ordered_coding_gene_region$length<=200000,]
length_regions
write.table(length_regions,'length_regions.tsv',sep = "\t", col.names = TRUE,row.names = FALSE,quote=FALSE)


order <- ordered_coding_gene_region
dim(order[order$score == 100,])
order$score[order$strand=="+" &
              order$length>200000 &
              (order$chr=="chr4" | order$chr=="chr17")] <- 100
dim(order[order$score == 100,])

str(ordered_coding_gene_region)
new_row <- data.frame(c("chr42"),
                         c(42L),
                         c(84L),
                         c("ENSG00000174489"),
                         c(0L),
                         c("+"),
                         c(42L))
str(new_row)
colnames(new_row) <- c("chr","start","end","ensemblID","score","strand","length")
new_row
new_data_frame <- rbind(ordered_coding_gene_region, new_row)
tail(new_data_frame)

new_data_frame <- new_data_frame[,-5]

ranges <- apply(ordered_coding_gene_region, 2, function(x){range(x, na.rm = TRUE, finite = TRUE)})
ranges

colours_vector <- c("red","orange","purple","yellow","pink","blue")
nchar(colours_vector[[1]])

if length(colours_vector[[x]]) in colours_vector{
  
}

for (i in colours_vector){
  if (nchar(i) == 4){
    print(i)
  }
}

if 
for i in colours_vector{
  
  colours_vector[c(print_positions)]
  }
