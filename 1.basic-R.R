# Lesson 1. Basic R Script for scRNA-seq Course

# This script provides an introduction to R programming basics
# for students learning single-cell RNA sequencing analysis.

# 1. Basic Arithmetic
# R can be used as a simple calculator
2 + 3  # Addition
10 - 4  # Subtraction
5 * 2  # Multiplication
20 / 4  # Division
2^3  # Exponentiation (2 to the power of 3)

# 2. Variable Assignment
# Variables allow us to store and manipulate data
x <- 10  # Assign the value 10 to x
x  # Print the value of x

y <- 5  # Assign the value 5 to y
y  # Print the value of y

z <- x + y  # Add x and y, and assign the result to z
print(z)  # Print the value of z using the print() function

# 3. Basic Data Types
# R has several basic data types. Here are the three most common:
num <- 10.5  # Numeric (for numbers)
class(num)  # Check the class (type) of num

char <- "Hello"  # Character (for text)
class(char)  # Check the class of char

log <- TRUE  # Logical (for TRUE/FALSE values)
class(log)  # Check the class of log

# 4. Vectors
# Vectors are one-dimensional arrays that can hold multiple values of the same type
numbers <- c(1, 2, 3, 4, 5)  # Create a numeric vector
numbers  # Print the vector

fruits <- c("apple", "banana", "cherry")  # Create a character vector
fruits  # Print the vector

bool_vec <- c(TRUE, FALSE, TRUE, TRUE)  # Create a logical vector
bool_vec  # Print the vector

# Creating a sequence of numbers
seq_numbers <- 1:10  # Create a vector with numbers from 1 to 10
seq_numbers  # Print the vector

# 5. Logical Operators
# These are used for comparisons and combining conditions
a <- 5
b <- 7

a < b   # Is a less than b?
a <= b  # Is a less than or equal to b?
a > b   # Is a greater than b?
a >= b  # Is a greater than or equal to b?
a == b  # Is a equal to b?
a != b  # Is a not equal to b?

# Combining conditions
(a > 3) | (b > 10)  # Is a > 3 OR b > 10? (| means OR)
(a > 3) & (b < 10)  # Is a > 3 AND b < 10? (& means AND)

# 6. Retrieving Elements from Vectors
numbers <- c(10, 20, 30, 40, 100)  # Create a new numeric vector
numbers[3]  # Get the third element
numbers[c(1, 3, 5)]  # Get the first, third, and fifth elements
numbers[2:4]  # Get elements from index 2 to 4

# Using logical vectors for indexing
numbers[c(TRUE, FALSE, TRUE, FALSE, TRUE)]  # Select elements where TRUE appears

# Conditional selection
numbers[numbers > 30]  # Select elements greater than 30

# Let's break down this expression:
numbers > 30  # This creates a logical vector
numbers[numbers > 30]  # Use this logical vector to select elements

# Combining conditions
numbers[numbers > 25 & numbers < 50]  # Select elements between 25 and 50

# Step-by-step breakdown:
high_expr <- numbers > 25  # Which elements are greater than 25?
high_expr

low_expr <- numbers < 50  # Which elements are less than 50?
low_expr

high_and_low_expr <- high_expr & low_expr  # Combine conditions with AND (&)
high_and_low_expr

numbers[high_and_low_expr]  # Use this logical vector to select elements

# Creating and using named vectors
gene_expr <- c(gene1 = 100, gene2 = 200, gene3 = 150, gene4 = 300)
print(gene_expr)

gene_expr["gene2"]  # Access element by name

gene_expr[c("gene1", "gene3")]  # Access multiple elements by name

# 7. Changing and Removing Elements
numbers[2] <- 1000  # Change the second element to 1000
numbers

numbers <- numbers[-3]  # Remove the third element
numbers

numbers <- c(numbers, 60, 70)  # Add new elements to the end of the vector
numbers

gene_expr["gene5"] <- 250  # Add a new named element
gene_expr

gene_expr["gene2"] <- 180  # Change an existing element by name
gene_expr

# 8. Installing and Loading Packages
# Packages extend R's functionality. Here's how to install and load them:
# install.packages("ggplot2")  # Install ggplot2 package (uncomment to run)
# library(ggplot2)  # Load ggplot2 package

# For Bioconductor packages:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")  # Install DESeq2 package

# library(DESeq2)  # Load DESeq2 package

# 9. Getting Help
# R provides built-in help documentation
?mean  # Get help for the mean() function
help(mean)  # Another way to get help for mean()

# 10. Data Frames
# Data frames are table-like structures for storing data
df <- data.frame(
  sample = c("sample1", "sample2", "sample3", "sample4"),
  expression = c(100, 200, 150, 300),
  condition = c("control", "treatment", "control", "treatment")
)
print(df)  # Print the data frame

# Getting information about the data frame
dim(df)  # Dimensions of the data frame
nrow(df)  # Number of rows
ncol(df)  # Number of columns

colnames(df)  # Column names

# Accessing data in a data frame
df$sample  # Access the 'sample' column
df[["sample"]]  # Another way to access the 'sample' column
df["sample"]  # Yet another way (returns a data frame)

# Subsetting a data frame
df[1:2, 2:3]  # Select rows 1-2 and columns 2-3
df[1, ]  # Select the first row and all columns
df[df$condition == "control", c("sample", "condition")]  # Select rows where condition is "control", and only the "sample" and "condition" columns

# Adding a new column
df$log_expression <- log2(df$expression)  # Add a new column with log2 transformed expression
df

# Adding a new row
new_row <- data.frame(sample = "sample5", expression = 250, condition = "control", log_expression = log2(250))
df <- rbind(df, new_row)  # Add the new row to the data frame
df

# Changing values in a data frame
df[df$condition == "control", "condition"] <- "C"  # Change "control" to "C"
df[df$condition == "treatment", "condition"] <- "T"  # Change "treatment" to "T"
print(df)

# 11. Matrices
# Matrices are 2-dimensional structures where all elements are of the same type
m <- matrix(1:12, nrow = 3, ncol = 4)  # Create a 3x4 matrix
print(m)

# Adding row and column names
rownames(m) <- c("gene1", "gene2", "gene3")
colnames(m) <- c("cell1", "cell2", "cell3", "cell4")
print(m)

# Matrix operations
t(m)  # Transpose the matrix
m * 2  # Multiply each element by 2

# Apply functions to rows or columns
colSums(m)  # Sum of each column
rowMeans(m)  # Mean of each row

# Subsetting matrices
m[1, ]  # First row
m[, 2]  # Second column
m[1:2, 3:4]  # Submatrix (rows 1-2, columns 3-4)
m[c("gene1","gene3"), c(2,4)]  # Select specific rows and columns by name

# 12. Factors
# Factors are used to represent categorical data
symptoms <- factor(c("mild", "severe", "mild", "moderate", "moderate"))
print(symptoms)

levels(symptoms)  # Get the levels of the factor

# Changing factor levels
levels(symptoms) <- c("M", "Mo", "S")  # Change levels to abbreviations
print(symptoms)

# Using factors in data analysis
counts <- table(symptoms)  # Count occurrences of each level
counts

# Visualizing factor data
barplot(counts, main="Symptom Severity", xlab="Severity Level")

# Reordering factor levels
symptoms_reordered <- factor(symptoms, levels = c("S", "Mo", "M"))
print(symptoms_reordered)

counts_reordered <- table(symptoms_reordered)
barplot(counts_reordered, main="Symptom Severity (Reordered)", xlab="Severity Level")

# 13. Lists
# Lists can contain elements of different types
patient <- list(
  id = "PT001",
  age = 45,
  symptoms = c("fever", "cough"),
  test_results = data.frame(
    test = c("PCR", "Antibody"),
    result = c("Positive", "Negative")
  )
)
print(patient)

# Accessing list elements
patient$id  # Access by name using $
patient[["age"]]  # Access by name using [[]]
patient[[3]]  # Access by index

# Adding to a list
patient$medication <- c("Aspirin", "Cough Syrup")
print(patient)

# Accessing nested list elements
patient$test_results$result

# 14. If-Else Statements
# Used for conditional execution of code
gene_expression <- 100

if (gene_expression > 50) {
  print("High expression")
} else {
  print("Low expression")
}

# Multiple conditions
cd4_expression <- 80
cd8_expression <- 20

if (cd4_expression > 50 & cd8_expression < 30) {
  print("This cell is likely a CD4+ T cell")
} else if (cd8_expression > 50 & cd4_expression < 30) {
  print("This cell is likely a CD8+ T cell")
} else {
  print("Cell type is uncertain")
}

# 15. Basic For Loops
# Used for repeating operations
for (i in 1:5) {
  print(i)
}

# Example: Calculating mean expression for multiple genes
genes <- c("GENE1", "GENE2", "GENE3", "GENE4")
expression_values <- list(
  GENE1 = c(10, 20, 15, 25),
  GENE2 = c(50, 60, 55, 65),
  GENE3 = c(5, 8, 6, 7),
  GENE4 = c(100, 120, 110, 130)
)
print(expression_values)

for (gene in genes) {
  mean_expression <- mean(expression_values[[gene]])
  cat("Mean expression of", gene, ":", mean_expression, "\n")
}

# 16. Working with File Paths
# Important for reading/writing files
current_dir <- getwd()  # Get current working directory
print(current_dir)

# setwd("/path/to/your/directory")  # Set working directory (uncomment and modify as needed)

# Construct file paths
results_dir <- "results"
csv_file <- file.path(results_dir, "gene_expression.csv")
print(csv_file)

# Create directory if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
  print(paste("Created directory:", results_dir))
}

# 17. Writing and Reading CSV Files
# Create a sample data frame
gene_data <- data.frame(
  gene_id = c("GENE1", "GENE2", "GENE3", "GENE4"),
  expression_level = c(10.5, 20.3, 15.7, 30.2),
  p_value = c(0.001, 0.05, 0.01, 0.001)
)
print(gene_data)

# Write the data frame to a CSV file
write.csv(gene_data, file = csv_file, row.names = FALSE)
print(paste("Data written to CSV file:", csv_file))

# Read the CSV file we just created
read_data <- read.csv(csv_file)
print(read_data)

# Read a specific number of rows
gene_data_head <- read.csv(csv_file, nrows = 2)
print("First 2 rows of data:")
print(gene_data_head)

# 18. Writing and Reading RDS Files
# RDS files can store any R object
rds_file <- file.path(results_dir, "gene_expression.rds")
saveRDS(gene_data, file = rds_file)
print(paste("Data written to RDS file:", rds_file))

read_rds_data <- readRDS(rds_file)
print(read_rds_data)

# Install Seurat package (uncomment to run)
# install.packages("Seurat")

# This script provides a foundation for R programming.
# In the next lessons, we'll explore the Seurat package for scRNA-seq analysis.
