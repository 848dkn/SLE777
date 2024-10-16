# SLE777
## Part 1: Gene Expression and Growth Data Analysis
## Introduction
In modern data science, the ability to process, analyze, and interpret complex datasets is essential for drawing meaningful insights. This report involves building proficiency in data wrangling, performing mathematical operations, and creating visualizations using R, one of the most powerful tools for statistical computing. Additionally, mastering collaboration platforms like GitHub is crucial for effective teamwork, version control, and sharing code.
## Packages Use:
```ggplot2```: Use for data visualization (histograms, box plots).
```dplyr```: Use for data manipulation, grouping, summarizing, and filtering.
```tidyr```: Use for reshaping data.
## Dataset:
``` gene_expression.tsv ```: A TSV file containing gene expression levels across various genes.
``` growth_data.csv ```: A CSV file containing tree circumference measurements at different time points (2005, 2010, 2015, 2020) and at different sites (northeast and southwest).

# Objectives of the Analysis
Compute an average expression for each gene and exhibit the distribution of protein expression levels overall to provide a summary of gene expression.
Detect Highly Expressed Genes: Develop a list comprising the top ten markers with the highest average expression values; these genes might be significant factors in the physiological procedure under consideration.

 
 

## Code Summary and Output
The code specifically looks to understand the distribution of gene expression from a given dataset across genes and tries to find key patterns. The following is what the above code accomplishes.
Read and Prepare the Gene Expression Data The following code reads in gene expression data from a tab-separated file where gene identifiers are set to row names.
The code computes the mean gene expression for each gene (row) across all of these conditions(columns). It then adds a new column Mean_Expression to its dataset, which has this calculated value.
Top Genes: The code ranks the dataset using values from Mean_Expression and generates a list of 10 genes with maximum mean expression.
Low-Expression Genes: It tells you a count of genes with an average expression value <10, which might indicate either non-expressed or lowly expressed under experimental conditions.
Plots the Data: To better handle a broad spectrum of expression values, a log10 transformation is performed to develop a histogram that illustrates the range of the mean protein expression values.

## Purpose of the Code
### Installing a Package: R, like a lot of languages used for programming, gets its novel functionality and features from other libraries or packages. One of the most commonly employed R utilities for complex and adaptable representations of data is the ggplot2 package.

#Explanation
The ggplot2 package has been added and pulled from the CRAN (Comprehensive R Archive Network) online resource for R packages via the `install.packages("ggplot2")` command. R retrieves the ggplot2 package on CRAN, installs it on the machine where you are, and ultimately connects to it when it is executed. Usually, only in an instance that the package isn't already installed or once per session is this step important. You can load and use the bundle in current and future R sessions after it is successfully installed.

### ggplot2 Package
ggplot2 is a common R tool that facilitates the creation of a wide range of visualizations by adhering to the grammar of graphics, a lucid and consistent design framework. Charts like scatter plots, line graphs, and histograms are simple to make with ggplot2. Also, it offers you a plethora of personalization choices to precisely adjust things such as colors, labels, and scales to match your preferred visual style.


### 1.1 Reading the file and setting gene identifiers as row names

```ruby
#download required files and packages
install.packages("ggplot2")
install.packages("dplyr")
install.packages("tidyr")
download.file("https://raw.githubusercontent.com/ghazkha/Assessment4/refs/heads/main/gene_expression.tsv", destfile = "gene_expression.tsv")

download.file("https://raw.githubusercontent.com/ghazkha/Assessment4/refs/heads/main/growth_data.csv", destfile = "growth_data.csv")

#reading the file into R

gene_expression <- read.delim("gene_expression.tsv", header = TRUE, row.names = 1)

#displaying the first six rows
first6row <- head(gene_expression, n = 6)
knitr::kable(first6row)
```
#### Input: 

The ```download.file``` function is used to download files from a URL with the argument ```destfile``` specifies the destination path where the file will be saved locally. 
The command ```read.delim``` reads a tab-separated file (.tsv) into a data frame. In this command, the argument ```header = TRUE``` and ```name =1```  is use to indicate that the first row contains the column names row and specify that the first column ```gene identifiers``` will be used as the row names of the data frame.
```knitr::kable()```function from the knitr package is used for generating dynamic reports in R Markdown similar to ```print()``` in R script.

#### Output:
The gene expression data is read into R as a data frame named ```gene_expression``` use for data manipulation. A table showing the first 6 genes.

### 1.2 Adding a new column for the mean of other columns

```ruby
#calculating the mean of expression values across the columns
rowMeans(gene_expression) 
#give the mean of expression values a name
meanofrow <- rowMeans(gene_expression) 
#adding new column into table
gene_expression$Mean_Expression <- meanofrow 
#show a table of values for the first six genes
table_of_values <- head(gene_expression, 6)  
knitr::kable(table_of_values)  
````
#### Input: 
Command ```rowMeans()``` computes the mean of each row, averaging expression values across the different sample columns. The mean values are saved into a variable to be added as a new column. The ```head()``` command is used to display the first six rows of the updated data frame for quick inspection. 
```knitr::kable()```function from the knitr package is used for generating dynamic reports in R Markdown similar to ```print()``` in R script.
#### Output: 
A table showing the first six genes, including their expression values and the new mean expression column.


### 1.3 List the 10 genes with the highest mean expression

```ruby
#sort by Mean_Expression and list top 10 in a table 
top_genes <- gene_expression[order(-gene_expression$Mean_Expression), ]
table_of_top10 <- top_genes[1:10, ]
knitr::kable(table_of_top10)
```
#### Input: 
The ```order()```command sorts the data frame in descending order. The first 10 rows are sorted from the data frame as the top 10 genes.
```knitr::kable()```function from the knitr package is used for generating dynamic reports in R Markdown similar to ```print()``` in R script.
#### Output: 
A table of the top 10 genes with the highest mean expression values.


### 1.4 Determine the number of genes with a mean <10

```ruby
#count and save number of genes with Mean_Expression < 10 as a value
low_expression_genes <- sum(gene_expression$Mean_Expression < 10)
knitr::kable(low_expression_genes)

```
#### Input: 
Command ```sum()``` counts the number of genes in column ```Mean_Expression``` that are lower than 10
```knitr::kable()```function from the knitr package is used for generating dynamic reports in R Markdown similar to ```print()``` in R script.
### Output: 
The count of genes where the mean expression is less than 10.


### 1.5  Plot Histogram of Mean Expression

```ruby
#plot package
library(ggplot2)
#Apply a log10 transformation to the Mean Expression for better visualization
#Added +1 to avoid log of zero
ggplot(gene_expression, aes(x = log10(Mean_Expression + 1))) +  
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black") +
  labs(title = "Histogram of Log10 Mean Gene Expression", x = "Log10 Mean Expression", y = "Frequency") +  theme_minimal()
```
#### Input: 
Command ```ggplot()``` creates a histogram, with ```x = log10(Mean_Expression)```. 
Command ```geom_histogram()``` creates the histogram, and various parameters (```binwidth```, ```fill```, ```color```) control the appearance of the plot.
Command ```labs()``` adds titles and axis labels, and ```theme_minimal()``` applies a minimalistic style to the plot (this is optional).
#### Output: 
A histogram visualizing the distribution of log-transformed mean gene expression values.


### 1.6 Import Tree Circumference Growth Data and display column names
```ruby
#reading the CSV file
growth_data <- read.csv("growth_data.csv")
# view column names in console
colnames(growth_data)
```
#### Input: 
Command ```read.csv()``` reads the CSV file into an R data frame, making the data ready for analysis. 
#### Output:
A data frame named ```growth_data```, containing the information from the CSV file.

### 1.7 Calculate the mean and standard deviation of tree circumference

```ruby
#load necessary libraries

library(dplyr)

#calculate the mean and standard deviation of tree circumferences at the start (2005) and end (2020) of the study for both control and treatment sites
summary_stats <- growth_data %>%
  group_by(Site) %>%
  summarise(
    Mean_Circumf_2005 = mean(Circumf_2005_cm, na.rm = TRUE),
    SD_Circumf_2005 = sd(Circumf_2005_cm, na.rm = TRUE),
    Mean_Circumf_2020 = mean(Circumf_2020_cm, na.rm = TRUE),
    SD_Circumf_2020 = sd(Circumf_2020_cm, na.rm = TRUE)  )
#print the summary statistics
knitr::kable(summary_stats)
```
#### Input: 
Command ```library()``` load the required dplyr package. 
Command ```group_by()``` groups the data by the Site column
Command ```summarise()``` calculates summary statistics for each group of data from the command mean() and sd(), with the ```na.rm = TRUE``` argument ensuring that missing values (NA) are ignored during these calculations.
```knitr::kable()```function from the knitr package is used for generating dynamic reports in R Markdown similar to ```print()``` in R script.
#### Output: 
A new data frame name ```summary_stats``` containing the mean and standard deviation of tree circumferences at both sites in 2005 and 2020.


### 1.8 Box plot of tree circumference at the start and end of the study at both sites
 
```ruby
# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
#gather the data for Circumf_2005_cm and Circumf_2020_cm into long format
long_data <- growth_data %>%
  select(Site, Circumf_2005_cm, Circumf_2020_cm) %>%
  pivot_longer(cols = starts_with("Circumf"), 
               names_to = "Year", 
               values_to = "Circumference")
#create the box plot
ggplot(long_data, aes(x = Year, y = Circumference, fill = Site)) +
  geom_boxplot() +
  labs(title = "Box Plot of Tree Circumference at Start and End of Study At Both Sites",
       x = "Year",
       y = "Tree Circumference (cm)") +
  scale_x_discrete(labels = c("Circumf_2005_cm" = "2005", "Circumf_2020_cm" = "2020")) +  theme_minimal()
```  
#### Input: 
Command ```library()``` load the required dplyr, tidyr and ggplot2 package
Command ```pivot_longer()``` convert the wide-format data (with separate columns for each year) into long format that just has one column for the year and one for the circumference. The data frame created is turned into a plot using command ```ggplot()```  and ```geom_boxplot()``` which turn it into a boxplot with various commands to edit it including ```fill = Site``` differentiates the two sites using different colors, ```labs()``` adds title and axis labels and ```scale_x_discrete()``` name the axis lable.
#### Ouput:
A box plot comparing tree circumferences at both sites in 2005 and 2020.

### 1.9 Calculation of the mean growth over the last 10 years at each site

```ruby
#calculate the growth for each tree
growth_data <- growth_data %>%
  mutate(Growthover10years = Circumf_2020_cm - Circumf_2010_cm)  #mutate add a new column that contain calculation

#calculate the mean growth at each site
mean_growth <- growth_data %>%
  group_by(Site) %>%
  summarise(Mean_Growth = mean(Growthover10years))  # Compute mean growth
#print the mean growth
knitr::kable(mean_growth)
```
#### Input: 
The pipe operator ```%>%``` takes the ```growth_data``` data frame and passes it as input to the next function in the chain. 
Command ```mutate()``` adds a new column to the data frame that has a difference between the values of the start (2005) and end (2020) of the study. 
Command ```group_by()``` groups the data by site.
Command ```summarise(Mean_Growth = mean())``` summarises the average 10-year growth for each site.
```knitr::kable()```function from the knitr package is used for generating dynamic reports in R Markdown similar to ```print()``` in R script.
#### Output: 
A summary table mean_growth that shows the mean growth over 10 years at each site.


### 1.10 Perform a t-test to estimate the p-value
```ruby
#perform t-test to compare the growth between the two sites
t_test_result <- t.test(Growthover10years ~ Site, data = growth_data)

# Extract the t-test results into a data frame
t_test_summary <- data.frame(
  Statistic = t_test_result$statistic,
  P_Value = t_test_result$p.value,
  Conf_Interval_Lower = t_test_result$conf.int[1],
  Conf_Interval_Upper = t_test_result$conf.int[2],
  Mean_Difference = t_test_result$estimate[1] - t_test_result$estimate[2]
)

# Print the t-test results using knitr::kable
knitr::kable(t_test_summary, caption = "T-Test Results: Growth Over 10 Years Between Sites")

```
#### Input: 
Command ```t.test()``` performs a t-test to determine the significant difference in 10-year growth between the two sites.
The ```~ Site``` grouping variable specifies that the growth differences should be compared between the sites.
```knitr::kable()```function from the knitr package is used for generating dynamic reports in R Markdown similar to ```print()``` in R script.
#### Output: 
The results of a t-test comparing the 10-year growth between the two sites

## Part 2: Examining Biological Sequence Diversity
##Introduction
This report analyses the genomic sequence diversity between Escherichia coli (E. coli) and Mycoplasma hyopneumoniae strain ES-2 (GCA_004768725). The key genomic features are compared by analysing coding DNA sequences (CDS), such as the number and length of coding sequences, nucleotide and amino acid compositions, and codon usage bias. Visualisations, such as boxplots and bar graphs, effectively highlight these differences. Furthermore, over- and under-represented k-mers in protein sequences are analysed, exploring how evolutionary pressures may have influenced the genomes of these two bacteria.
 
##Packages Use
```Biostrings```: Use for handling DNA and protein sequences.
```R.utils```: Use for file manipulation, memory-efficient data handling, and working with compressed files.
```seqinr```: Use for sequence input/output and codon usage calculations.
```ggplot2```: Use for creating visualizations such as histograms, box plots, and bar plots.
##Dataset
```E. coli CDS```: Coding DNA sequences for E. coli downloaded from the Ensembl Bacteria database.
```Mycoplasma hyopneumoniae CDS```: Coding DNA sequences for Mycoplasma hyopneumoniae strain ES-2 downloaded from the Ensembl Bacteria database.


```ruby
### Install and load packages
install.packages("Biostrings")
install.packages("R.utils")
install.packages("ggplot2")
install.packages("seqinr")
library(Biostrings)
library("R.utils")
library("ggplot2")
library("seqinr")
```
### 2.1 Download and count Coding Sequences

```ruby
```{r download_count_cds, echo=FALSE, eval=FALSE}
#Download coding sequences for Mesomycoplasma hyopneumoniae
mh_url <- "https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-59/fasta/bacteria_40_collection/mesomycoplasma_hyopneumoniae_gca_004768725/cds/Mesomycoplasma_hyopneumoniae_gca_004768725.ASM476872v1.cds.all.fa.gz"
download.file(mh_url, destfile = "mesomycoplasma_cds.fa.gz")
gunzip("mesomycoplasma_cds.fa.gz")
unzipmycoplasma <- seqinr::read.fasta("mesomycoplasma_cds.fa")

#Download coding sequences for E. coli
ecoli_url <- "http://ftp.ensemblgenomes.org/pub/bacteria/release-53/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/cds/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.cds.all.fa.gz"
download.file(ecoli_url, destfile = "ecoli_cds.fa.gz")
gunzip("ecoli_cds.fa.gz")
unzipecoli <- seqinr::read.fasta("ecoli_cds.fa")

#Read the sequences from the unzipped files
Ecoli_sequences <- Biostrings::readDNAStringSet("ecoli_cds.fa")
Mycoplasma_sequences <- Biostrings::readDNAStringSet("mesomycoplasma_cds.fa")

#Count the number of sequences in both FASTA files
num_sequences <- data.frame(Organism = c("E. coli", "Mycoplasma hyopneumoniae"),
                            Num_CDS = c(length(Ecoli_sequences), length(Mycoplasma_sequences)))

# Print the result using knitr::kable
knitr::kable(num_sequences)
```
#### Input:
The ```download.file``` function is used to download FASTA files containing the coding sequences (CDS) for Escherichia coli and Mycoplasma hyopneumoniae. The downloaded files are compressed in ```.fa.gz``` format and unzipped using the ```gunzip``` command. The unzipped CDS were then read into R using the ```seqinr::read.fasta``` and ```Biostrings::readDNAStringSet``` functions for handling biological sequence data. Finally, the ```length``` function counts the number of coding sequences in each organism.

#### Output:
A table showing the number of coding sequences for both organisms, E. coli and Mycoplasma hyopneumoniae.

### 2.2 Calculate total coding DNA for both organisms

```ruby
```{r calculate_total, echo=FALSE, eval=FALSE}
#Calculate total coding DNA 
total_coding <- data.frame(Organism = c("E. coli", "Mycoplasma hyopneumoniae"),
                           Total_Coding_DNA = c(sum(width(Ecoli_sequences)), sum(width(Mycoplasma_sequences))))
total_coding #Present table for total coding DNA


#Calculate total length of coding sequences for E. coli and Mycoplasma hyopneumoniae
total_length_ecoli <- sum(width(Ecoli_sequences))
total_length_ecoli
total_length_mycoplasma <- sum(width(Mycoplasma_sequences))

#Create a table with total coding DNA length
total_coding_length <- data.frame(
  Organism = c("E. coli", "Mycoplasma hyopneumoniae"),
  Total_Coding_Length = c(total_length_ecoli, total_length_mycoplasma)
)


# Print the result using knitr::kable
knitr::kable(total_coding_length)


```
#### Input: The ```width``` function is applied to both E. coli and Mycoplasma hyopneumoniae sequences to calculate the length (in base pairs) of each coding sequence (CDS). The ```sum``` function is then used to compute the total length of all coding DNA for each organism. The results are stored in a ```data.frame``` showing the total coding DNA length for both organisms.
#### Output: A table displaying the total length of coding DNA (in base pairs) for E. coli and Mycoplasma hyopneumoniae.


### 2.3 Calculate the length of all coding sequences and box plot


```ruby
```{r calculate_length, echo=FALSE, message=FALSE, warning=FALSE}
#Calculate the length of all coding sequences for both organisms
lengths_ecoli <- width(Ecoli_sequences)  # Length of CDS of E. coli
lengths_mycoplasma <- width(Mycoplasma_sequences)  # Length of CDS of Mycoplasma hyopneumoniae

#Create a data frame for the lengths
sequence_lengths <- data.frame(
  Organism = rep(c("E. coli", "Mycoplasma hyopneumoniae"), 
                 times = c(length(lengths_ecoli), length(lengths_mycoplasma))),
  Length = c(lengths_ecoli, lengths_mycoplasma)
)

#Load ggplot2 for visualization
library(ggplot2)

#Create boxplot
library(ggplot2)
ggplot(sequence_lengths, aes(x = Organism, y = Length)) +
  geom_boxplot(fill = c("#FF9999", "#99CCFF")) +
  labs(title = "Boxplot of Coding Sequence Lengths", x = "Organism", y = "Sequence Length (bp)") +
  theme_minimal()

#Calculate mean and median for E. coli
mean_length_ecoli <- mean(lengths_ecoli)
median_length_ecoli <- median(lengths_ecoli)

#Calculate mean and median for Mycoplasma hyopneumoniae
mean_length_mycoplasma <- mean(lengths_mycoplasma)
median_length_mycoplasma <- median(lengths_mycoplasma)

#Create a table for mean and median lengths
mean_median_lengths <- data.frame(
  Organism = c("E. coli", "Mycoplasma hyopneumoniae"),
  Mean_Length = c(mean_length_ecoli, mean_length_mycoplasma),
  Median_Length = c(median_length_ecoli, median_length_mycoplasma)
)

# Print the result using knitr::kable
knitr::kable(mean_median_lengths)

```
#### Input: The ```width``` function is used to calculate the length (in base pairs) of each coding sequence (CDS) for both E. coli and Mycoplasma hyopneumoniae. The lengths are then stored in a ```data.frame``` that includes the organism's name and the corresponding sequence lengths. Using ```ggplot2```, a boxplot is created to visualize the distribution of coding sequence lengths for both organisms.
Additionally, the ```mean``` and ```median``` functions are used to calculate the average and median lengths of coding sequences for both organisms.

#### Output: A boxplot of the coding sequence lengths for E. coli and Mycoplasma hyopneumoniae and a table displaying the mean and median coding sequence lengths for both organisms.
 

### 2.4 The frequency of DNA bases in the total coding sequences for both organisms

```ruby
#Calculate base frequencies for E. coli
ecoli_base_frequencies <- alphabetFrequency(Ecoli_sequences, baseOnly = TRUE)

#Calculate base frequencies for Mycoplasma hyopneumoniae
mycoplasma_base_frequencies <- alphabetFrequency(Mycoplasma_sequences, baseOnly = TRUE)

ecoli_frequencies <- colSums(ecoli_base_frequencies)
mycoplasma_frequencies <- colSums(mycoplasma_base_frequencies)


#Create a data frame for base frequencies
base_frequencies <- data.frame(
  Organism = rep(c("E. coli", "Mycoplasma hyopneumoniae"), each = 5),
  Base = rep(c("A", "T", "G", "C", "other")),
  Frequency = c(ecoli_frequencies, mycoplasma_frequencies))

print(base_frequencies)
# Print the result using knitr::kable
knitr::kable(base_frequencies)
```
#### Input: The ```alphabetFrequency``` function calculates DNA and amino acid base frequencies in E. coli and Mycoplasma hyopneumoniae coding sequences. These frequencies are aggregated and visualized using ```colSums``` and ```ggplot2```.
#### Output: A data frame of the frequency of DNA bases, amino acids and bar plot of nucleotide and amino acids frequencies between E. coli and Mycoplasma hyopneumoniae.


### Calculate frequency for total protein sequence

```ruby
#Translate DNA sequences to protein sequences
ecoli_protein_sequences <- Biostrings::translate(Ecoli_sequences)
mycoplasma_protein_sequences <- Biostrings::translate(Mycoplasma_sequences)

#Calculate amino acid frequencies for E. coli
ecoli_amino_frequencies <- alphabetFrequency(ecoli_protein_sequences)

#Calculate amino acid frequencies for Mycoplasma hyopneumoniae
mycoplasma_amino_frequencies <- alphabetFrequency(mycoplasma_protein_sequences)

#Create a data frame for amino acid frequencies
amino_frequencies <- data.frame(
  Organism = rep(c("E. coli", "Mycoplasma hyopneumoniae"), each = 62),
  Amino_Acid = rep(letters[1:31], times = 2),  # Assuming 20 amino acids
  Frequency = c(colSums(ecoli_amino_frequencies), colSums(mycoplasma_amino_frequencies)))


print(colSums(ecoli_amino_frequencies))
print(colSums(mycoplasma_amino_frequencies))

# Print the result using knitr::kable
knitr::kable(colSums(ecoli_amino_frequencies))

# Print the result using knitr::kable
knitr::kable((mycoplasma_amino_frequencies))

#Plot nucleotide frequencies
ggplot(base_frequencies, aes(x = Base, y = Frequency, fill = Organism)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(title = "Nucleotide Frequency Comparison", x = "Base", y = "Frequency")
```

```ruby
#Plot amino acid frequencies
ggplot(amino_frequencies, aes(x = Amino_Acid, y = Frequency, fill = Organism)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(title = "Amino Acid Frequency Comparison", x = "Amino Acid", y = "Frequency")
```

### 2.5 codon usage table and quantify the codon usage

```ruby
#Calculate codon usage for E. coli
ecoli_codon_usage <- uco(unzipecoli)

#Calculate codon usage for Mycoplasma hyopneumoniae
mycoplasma_codon_usage <- uco(unzipmycoplasma)
class(Mycoplasma_sequences)  # Kiểm tra loại đối tượng
mycoplasma_vector <- as.character(Mycoplasma_sequences)

#Create a table for codon usage
codon_usage <- data.frame(
  Codon = names(ecoli_codon_usage),
  Ecoli_Usage = ecoli_codon_usage,
  Mycoplasma_Usage = mycoplasma_codon_usage
)

print(ecoli_codon_usage)
print(mycoplasma_codon_usage)
codon_usage

# Print the result using knitr::kable
knitr::kable(ecoli_codon_usage)

# Print the result using knitr::kable
knitr::kable(mycoplasma_codon_usage)

```

#### Input: The ```uco``` function is used to calculate codon usage for both E. coli and Mycoplasma hyopneumoniae. Coding sequences from both organisms were converted into nucleotide vectors, enabling the computation of codon usage statistics. A comparative data frame printed using ```knitr: :kable``` is subsequently generated to analyse codon usage differences between the two species.
#### Output: A codon usage table displays each codon's frequency for E. coli and Mycoplasma hyopneumoniae.



### 2.6 Comparison of Over and Under Represented Protein Sequence k-mers (3-5) in M. hyopneumoniae and E. coli

### Mycoplasma hyopneumoniae k-mer
```ruby
#Function to calculate k-mer frequencies for protein sequences
calculate_kmers <- function(sequences, k) {
  kmers <- sapply(sequences, function(seq) {
    # Convert sequence to character string
    seq <- as.character(seq)
    # Extract all possible k-mers from the sequence
    kmer_list <- sapply(1:(nchar(seq) - k + 1), function(i) {
      substr(seq, i, i + k - 1)
    })
    return(kmer_list)
  })
  
  # Flatten the list of k-mers and calculate frequencies
  kmers <- unlist(kmers)
  kmer_freqs <- table(kmers)
  return(kmer_freqs)
}
```

```ruby
#Calculate 3-mer, 4-mer, and 5-mer frequencies for Mycoplasma hyopneumoniae protein sequences
mycoplasma_kmers_3 <- calculate_kmers(mycoplasma_protein_sequences, 3)
mycoplasma_kmers_4 <- calculate_kmers(mycoplasma_protein_sequences, 4)
mycoplasma_kmers_5 <- calculate_kmers(mycoplasma_protein_sequences, 5)

#Combine the k-mer frequencies into one
mycoplasma_kmers <- c(mycoplasma_kmers_3, mycoplasma_kmers_4, mycoplasma_kmers_5)

#Sort and get top 10 over- and under-represented k-mers
top_10_kmers_mycoplasma <- head(sort(mycoplasma_kmers, decreasing = TRUE), 10)
bottom_10_kmers_mycoplasma <- head(sort(mycoplasma_kmers, decreasing = FALSE), 10)

#Print top 10 over- and under-represented k-mers
print(top_10_kmers_mycoplasma)
print(bottom_10_kmers_mycoplasma)

# Print the result using knitr::kable
knitr::kable(top_10_kmers_mycoplasma)
# Print the result using knitr::kable
knitr::kable(bottom_10_kmers_mycoplasma)
```

### Ecoli k-mer

```ruby
#Calculate k-mers for E. coli protein sequences
ecoli_kmers_3 <- calculate_kmers(ecoli_protein_sequences, 3)
ecoli_kmers_4 <- calculate_kmers(ecoli_protein_sequences, 4)
ecoli_kmers_5 <- calculate_kmers(ecoli_protein_sequences, 5)

#Combine k-mers into one table for E. coli
ecoli_kmers <- c(ecoli_kmers_3, ecoli_kmers_4, ecoli_kmers_5)

#Get top 10 over- and under-represented k-mers in E. coli
top_10_kmers_ecoli <- head(sort(ecoli_kmers, decreasing = TRUE), 10)
bottom_10_kmers_ecoli <- head(sort(ecoli_kmers, decreasing = FALSE), 10)

#Print results
print(top_10_kmers_ecoli)
print(bottom_10_kmers_ecoli)

# Print the result using knitr::kable
knitr::kable(top_10_kmers_ecoli)
# Print the result using knitr::kable
knitr::kable(bottom_10_kmers_ecoli)

``` 


## Plot mycoplasma vs ecoli k-mer

```ruby
# Combine top and bottom k-mers from both organisms into a single data frame
combined_kmers <- data.frame(
  Kmer = c(names(top_10_kmers_mycoplasma), names(bottom_10_kmers_mycoplasma),
           names(top_10_kmers_ecoli), names(bottom_10_kmers_ecoli)),
  Frequency = c(top_10_kmers_mycoplasma, bottom_10_kmers_mycoplasma,
                top_10_kmers_ecoli, bottom_10_kmers_ecoli),
  Organism = rep(c("Mycoplasma hyopneumoniae", "Mycoplasma hyopneumoniae", 
                   "E. coli", "E. coli"), each = 10),
  Representation = rep(c("Top 10", "Bottom 10"), times = 4)
)
```

```ruby
# Plot the k-mers
ggplot(combined_kmers, aes(x = Kmer, y = Frequency, fill = Organism)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Representation) +
  labs(title = "Comparison of Over and Under K-mers", 
       x = "K-mers", y = "Frequency") +
  theme_minimal()
```

#### Input: The ```calculate_kmers``` function determines k-mer frequencies (k=3, 4, 5) within the protein sequences. A k-mer represents a substring of length k extracted from the sequence. The function identified all possible k-mers and calculated their frequencies using the ```table``` function. K-mer frequencies are computed for Mycoplasma hyopneumoniae and E. coli protein sequences, resulting in 3-mer, 4-mer, and 5-mer frequency tables. To identify over- and under-represented k-mers, the tables are sorted in decreasing and increasing order, yielding the top 10 k-mers for each category.The "print" outputs the raw data to the console, while ```knitr::kable``` formats the same data into clean tables for better presentation in reports.
#### Output: The top 10 over- and under-represented k-mers for both organisms are listed, and a bar plot of E. coli and Mycoplasma hyopneumoniae compares these top and bottom k-mers for both species.
