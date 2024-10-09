# SLE777
## Part 1: Gene Expression and Growth Data Analysis

### 1.1 Reading the file and setting gene identifiers as row names
#download required files 
download.file("https://raw.githubusercontent.com/ghazkha/Assessment4/refs/heads/main/gene_expression.tsv", destfile = "gene_expression.tsv")

download.file("https://raw.githubusercontent.com/ghazkha/Assessment4/refs/heads/main/growth_data.csv", destfile = "growth_data.csv")

#reading the file into R

gene_expression <- read.delim("gene_expression.tsv", header = TRUE, row.names = 1)

#displaying the first six rows

head(gene_expression, n = 6)


### 1.2 Adding a new column for the mean of other columns
#calculating the mean of expression values across the columns
rowMeans(gene_expression) 
#give the mean of expression values a name
meanofrow <- rowMeans(gene_expression) 
#adding new column into table
gene_expression$Mean_Expression <- meanofrow 
#show a table of values for the first six genes
table_of_values <- head(gene_expression, 6)  
print(table_of_values)  


### 1.3 List the 10 genes with the highest mean expression

#sort by Mean_Expression and list top 10 in a table 
top_genes <- gene_expression[order(-gene_expression$Mean_Expression), ]
table_of_top10 <- top_genes[1:10, ]
View(table_of_top10)

### 1.4 Determine the number of genes with a mean <10

#count and save number of genes with Mean_Expression < 10 as a value
low_expression_genes <- sum(gene_expression$Mean_Expression < 10)
low_expression_genes


### 1.5  Plot Histogram of Mean Expression

#install plot package
install.packages("ggplot2")
library(ggplot2)
#Apply a log10 transformation to the Mean Expression for better visualization
#Added +1 to avoid log of zero
ggplot(gene_expression, aes(x = log10(Mean_Expression + 1))) +  
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black") +
  labs(title = "Histogram of Log10 Mean Gene Expression", x = "Log10 Mean Expression", y = "Frequency") +  theme_minimal()

