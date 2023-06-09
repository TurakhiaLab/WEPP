library("SQUAREM")

## Start measuring time
start_time <- Sys.time()

## INITIALIZE the variables
setwd("/home/AD.UCSD.EDU/pgangwar/usher")
#file_name <- "my_vcf_mismatch_matrix.csv"
file_name <- "test.csv"
alpha <- 0.0000001
read_len <- 150

## READ the Input file
# Read the CSV file as a matrix
matrix_data <- t(as.matrix(read.csv(file_name, header = FALSE, sep = ",")))
# Convert the matrix to an array with each sub-array representing one row
mismatch_matrix <- array(matrix_data, dim = c(1, nrow(matrix_data), ncol(matrix_data)))
num_hap <- dim(mismatch_matrix)[2]

## Objective Function -> Calculate (-)LOG-LIKELIHOOD
loglik <- function(prob, prob_read_from_haplotype) {
  mat_dimensions <- dim(prob_read_from_haplotype)
  log_liklihood <- 0
  #Iteration on Row
  for (i in 1:mat_dimensions[3]) {
    weighted_row_sum <- 0
    #Iteration on Column
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
  #Iteration on Column
  for (k in 1:mat_dimensions[2]) {
    prob_new[k] <- 0
    #Iteration on Row -> Getting the weighted sum of denominator
    for (i in 1:mat_dimensions[3]) {
      weighted_row_sum <- 0
      #Iteration on Column
      for (j in 1:mat_dimensions[2])
        weighted_row_sum <- weighted_row_sum + ((prob_read_from_haplotype[,j,i]) * prob[j])
      prob_new[k] <- prob_new[k] + ((prob_read_from_haplotype[,k,i] * prob[k])/  weighted_row_sum)
    }
    prob_new[k] <- (prob_new[k] / mat_dimensions[3])
  }
  return(prob_new)
}


##EVALUATE
#Generate random probabilities b/w 0 and 1. Normalize initial probabilities
prob0 <- runif(num_hap)
prob0 <- prob0 / sum(prob0)
prob_read_from_haplotype <- (alpha^mismatch_matrix) * ((1-alpha)^(read_len-mismatch_matrix))

f1 <- fpiter(par = prob0, fixptfn = em, objfn = loglik, control = list(tol = 1.e-08), prob_read_from_haplotype = prob_read_from_haplotype)
#Normalize the results
f1$par <- f1$par / sum(f1$par) 
print(f1$par)

#Check Warnings
if (length(warnings()) > 0)
  print(warnings())

# Stop measuring time
print(Sys.time() - start_time)