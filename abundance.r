library("SQUAREM")

## READ the Input file
# Set the working directory to the location where the CSV file is located
setwd("/home/AD.UCSD.EDU/pgangwar/usher")
# Read the CSV file as a matrix
matrix_data <- t(as.matrix(read.csv("test.csv", header = FALSE, sep = ",")))
# Convert the matrix to an array with each sub-array representing one row
mismatch_matrix <- array(matrix_data, dim = c(1, nrow(matrix_data), ncol(matrix_data)))

## INITIALIZE the variables
alpha <- 0.000001
read_len <- 150
prob0 <- c(0.3, 0.5, 0.1, 0.1)
prob_read_from_haplotype <- (alpha^mismatch_matrix) * ((1-alpha)^(read_len-mismatch_matrix))


###Objective Function
##Calculate LOG-LIKELIHOOD
loglik <- function(prob, prob_read_from_haplotype) {
  mat_dimensions <- dim(prob_read_from_haplotype)
  log_liklihood <- 0
  for (i in 1:mat_dimensions[3]) {
    weighted_row_sum <- 0
    for (j in 1:mat_dimensions[2])
      weighted_row_sum <- weighted_row_sum + ((prob_read_from_haplotype[,j,i]) * prob[j])
    log_liklihood <- log_liklihood + log(weighted_row_sum)
  }
  return(-log_liklihood)
}


## EM algorithm function
em <- function(prob, prob_read_from_haplotype) {
  prob_new <- rep(NA_real_, length(prob))
  mat_dimensions<-dim(prob_read_from_haplotype)
  for (k in 1:mat_dimensions[2]) {
    inv_weighted_row_sum <- 0
    #Getting the weighted sum of denominator
    for (i in 1:mat_dimensions[3]) {
      weighted_row_sum <- 0
      for (j in 1:mat_dimensions[2])
        weighted_row_sum <- weighted_row_sum + ((prob_read_from_haplotype[,j,i]) * prob[j])
      inv_weighted_row_sum <- 1 / (weighted_row_sum)
    }
    #Getting the column sum of q
    col_sum <- 0
    for (i in 1:mat_dimensions[3])
      col_sum <- col_sum + prob_read_from_haplotype[,k,i]
    prob_new[k] <- (col_sum * inv_weighted_row_sum * prob[k]) / (mat_dimensions[3])
  }
  return(prob_new)
}


##EVALUATE
print(f1 <- squarem(par = prob0, fixptfn = em, objfn = loglik, control = list(tol = 1.e-08), prob_read_from_haplotype = prob_read_from_haplotype))