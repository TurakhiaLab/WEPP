suppressPackageStartupMessages({
  library("SQUAREM")
  library(turboEM)
})

## Start measuring time
start_time <- Sys.time()

## INITIALIZE the variables
setwd("/home/AD.UCSD.EDU/pgangwar/usher")
#file_name <- "my_vcf_mismatch_matrix.csv"
file_name <- "test.csv"
alpha <- 0.005
read_len <- 150

## READ the Input file
# Read the CSV file as a matrix
matrix_data <- t(as.matrix(read.csv(file_name, header = FALSE, sep = ",")))
# Convert the matrix to an array with each sub-array representing one row
mismatch_matrix <- array(matrix_data, dim = c(1, nrow(matrix_data), ncol(matrix_data)))
#Dropping the extra dimension
mismatch_matrix <- drop(mismatch_matrix)
#Interchanging rows and columns to match the input file
mismatch_matrix <- aperm(mismatch_matrix, c(2, 1))
#Get numer of haplotypes
num_hap <- dim(mismatch_matrix)[2]


## Objective Function -> Calculate (-)LOG-LIKELIHOOD
loglik <- function(prob, prob_read_from_haplotype) {
  #A*B happens in R by multiplying each element of a column of matrix A with each element of the vector B
  weighted_row_sum <- rowSums(t(t(prob_read_from_haplotype) * prob))
  -sum(log(weighted_row_sum))
}


## EM algorithm function
em <- function(prob, prob_read_from_haplotype) {
  weighted_row_sum <- rowSums(t(t(prob_read_from_haplotype) * prob))
  prob_row_element <- sweep(prob_read_from_haplotype, 2, prob, "*")
  prob_new <- colSums(sweep(prob_row_element, 1, weighted_row_sum, "/"))
  prob_new / dim(prob_read_from_haplotype)[1]
}


##EVALUATE
#Generate random probabilities b/w 0 and 1. Normalize initial probabilities
prob0 <- runif(num_hap)
prob0 <- prob0 / sum(prob0)

# Calculate prob_read_from_haplotype outside the EM algorithm loop
prob_read_from_haplotype <- (alpha^mismatch_matrix) * ((1-alpha)^(read_len-mismatch_matrix))

f1 <- turboem(
  par = prob0, fixptfn = em, objfn = loglik, method = "squarem", 
  prob_read_from_haplotype = prob_read_from_haplotype, parallel = F,
  control.run = list(
    convtype = "parameter", tol = 1.0e-7,
    stoptype = "maxtime", maxtime = 14400
  )
)

#Normalize the results with absolute values
f1$par <- abs(f1$par)
f1$par <- f1$par / sum(f1$par)
# Convert the result to a single column
result <- matrix(f1$par, nrow = length(f1$par), ncol = 1)
print(result)

# Stop measuring time
print(Sys.time() - start_time)