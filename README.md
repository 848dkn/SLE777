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

### 1.5  Plot Histogram of Mean Expression
#install plot package
install.packages("ggplot2")
library(ggplot2)
#Apply a log10 transformation to the Mean Expression for better visualization
#Added +1 to avoid log of zero
ggplot(gene_expression, aes(x = log10(Mean_Expression + 1))) +  
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black") +
  labs(title = "Histogram of Log10 Mean Gene Expression", x = "Log10 Mean Expression", y = "Frequency") +  theme_minimal()
### 1.6 Import Tree Circumference Growth Data and display column names
#reading the CSV file
growth_data <- read.csv("growth_data.csv")
# view column names in console
colnames(growth_data)
### 1.7 Calculate the mean and standard deviation of tree circumference
#load necessary libraries
install.packages("dplyr")
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
print(summary_stats)

### 1.8 Box plot of tree circumference at the start and end of the study at both sites
# Load necessary libraries
install.packages("tidyr")
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
### 1.9 Calculation of the mean growth over the last 10 years at each site
#calculate the growth for each tree
growth_data <- growth_data %>%
  mutate(Growthover10years = Circumf_2020_cm - Circumf_2010_cm)  #mutate add a new column that contain calculation

#calculate the mean growth at each site
mean_growth <- growth_data %>%
  group_by(Site) %>%
  summarise(Mean_Growth = mean(Growthover10years))  # Compute mean growth
#print the mean growth
print(mean_growth)
### 1.10 Perform a t-test to estimate the p-value that the 10-year growth is different between the two sites
#perform t-test to compare the growth between the two sites
t_test_result <- t.test(Growthover10years ~ Site, data = growth_data)
#print the t-test results
print(t_test_result)

## Part 2: Examining Biological Sequence Diversity
### Install and load packages
install.packages("Biostrings")
install.packages("R.utils")
install.packages("ggplot2")
install.packages("seqinr")
library(Biostrings)
library("R.utils")
library("ggplot2")
library("seqinr")
###2.1 Download and count Coding Sequences

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

num_sequences


### 2.2 Calculate total coding DNA for both organisms
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

total_coding_length

### 2.3 Calculate the length of all coding sequences and box plot

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

mean_median_lengths

###2.4 The frequency of DNA bases in the total coding sequences for both organisms
# Calculate base frequencies for E. coli
ecoli_base_frequencies <- alphabetFrequency(Ecoli_sequences, baseOnly = TRUE)

# Calculate base frequencies for Mycoplasma hyopneumoniae
mycoplasma_base_frequencies <- alphabetFrequency(Mycoplasma_sequences, baseOnly = TRUE)

ecoli_frequencies <- colSums(ecoli_base_frequencies)
mycoplasma_frequencies <- colSums(mycoplasma_base_frequencies)


# Create a data frame for base frequencies
base_frequencies <- data.frame(
  Organism = rep(c("E. coli", "Mycoplasma hyopneumoniae"), each = 5),
  Base = rep(c("A", "T", "G", "C", "other")),
  Frequency = c(ecoli_frequencies, mycoplasma_frequencies))

print(base_frequencies)

#Calculate frequency for total protein sequence
# Translate DNA sequences to protein sequences
ecoli_protein_sequences <- Biostrings::translate(Ecoli_sequences)
mycoplasma_protein_sequences <- Biostrings::translate(Mycoplasma_sequences)

# Calculate amino acid frequencies for E. coli
ecoli_amino_frequencies <- alphabetFrequency(ecoli_protein_sequences)

# Calculate amino acid frequencies for Mycoplasma hyopneumoniae
mycoplasma_amino_frequencies <- alphabetFrequency(mycoplasma_protein_sequences)

# Create a data frame for amino acid frequencies
amino_frequencies <- data.frame(
  Organism = rep(c("E. coli", "Mycoplasma hyopneumoniae"), each = 62),
  Amino_Acid = rep(letters[1:31], times = 2),  # Assuming 20 amino acids
  Frequency = c(colSums(ecoli_amino_frequencies), colSums(mycoplasma_amino_frequencies)))


print(colSums(ecoli_amino_frequencies))
print(colSums(mycoplasma_amino_frequencies))



# Plot nucleotide frequencies
ggplot(base_frequencies, aes(x = Base, y = Frequency, fill = Organism)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(title = "Nucleotide Frequency Comparison", x = "Base", y = "Frequency")

# Plot amino acid frequencies
ggplot(amino_frequencies, aes(x = Amino_Acid, y = Frequency, fill = Organism)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(title = "Amino Acid Frequency Comparison", x = "Amino Acid", y = "Frequency")


###2.5 codon usage table and quantify the codon usage
# Calculate codon usage for E. coli
ecoli_codon_usage <- uco(unzipecoli)
# Calculate codon usage for Mycoplasma hyopneumoniae
mycoplasma_codon_usage <- uco(unzipmycoplasma)
class(Mycoplasma_sequences)  # Kiểm tra loại đối tượng
mycoplasma_vector <- as.character(Mycoplasma_sequences)

# Create a table for codon usage
codon_usage <- data.frame(
  Codon = names(ecoli_codon_usage),
  Ecoli_Usage = ecoli_codon_usage,
  Mycoplasma_Usage = mycoplasma_codon_usage
)

print(ecoli_codon_usage)
print(mycoplasma_codon_usage)
codon_usage


###2.6 Comparison of Over and Under Represented Protein Sequence k-mers (3-5) in Mycoplasma hyopneumoniae and E. coli

## Mycoplasma hyopneumoniae k-mer
# Function to calculate k-mer frequencies for protein sequences
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

# Calculate 3-mer, 4-mer, and 5-mer frequencies for Mycoplasma hyopneumoniae protein sequences
mycoplasma_kmers_3 <- calculate_kmers(mycoplasma_protein_sequences, 3)
mycoplasma_kmers_4 <- calculate_kmers(mycoplasma_protein_sequences, 4)
mycoplasma_kmers_5 <- calculate_kmers(mycoplasma_protein_sequences, 5)

# Combine the k-mer frequencies into one
mycoplasma_kmers <- c(mycoplasma_kmers_3, mycoplasma_kmers_4, mycoplasma_kmers_5)

# Sort and get top 10 over- and under-represented k-mers
top_10_kmers_mycoplasma <- head(sort(mycoplasma_kmers, decreasing = TRUE), 10)
bottom_10_kmers_mycoplasma <- head(sort(mycoplasma_kmers, decreasing = FALSE), 10)

# Print top 10 over- and under-represented k-mers
print(top_10_kmers_mycoplasma)
print(bottom_10_kmers_mycoplasma)


## Ecoli k-mer

# Calculate k-mers for E. coli protein sequences
ecoli_kmers_3 <- calculate_kmers(ecoli_protein_sequences, 3)
ecoli_kmers_4 <- calculate_kmers(ecoli_protein_sequences, 4)
ecoli_kmers_5 <- calculate_kmers(ecoli_protein_sequences, 5)

# Combine k-mers into one table for E. coli
ecoli_kmers <- c(ecoli_kmers_3, ecoli_kmers_4, ecoli_kmers_5)

# Get top 10 over- and under-represented k-mers in E. coli
top_10_kmers_ecoli <- head(sort(ecoli_kmers, decreasing = TRUE), 10)
bottom_10_kmers_ecoli <- head(sort(ecoli_kmers, decreasing = FALSE), 10)

# Print results
print(top_10_kmers_ecoli)
print(bottom_10_kmers_ecoli)

##Plot mycoplasma vs ecoli k-mer

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

# Plot the k-mers
ggplot(combined_kmers, aes(x = Kmer, y = Frequency, fill = Organism)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Representation) +
  labs(title = "Comparison of Over and Under K-mers", 
       x = "K-mers", y = "Frequency") +
  theme_minimal()


```ruby
require 'redcarpet'
markdown = Redcarpet.new("Hello World!")
puts markdown.to_html
```
