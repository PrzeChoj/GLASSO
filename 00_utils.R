is.positive.semi.definite.matrix <- function(matrix_of_interest, tolerance = 1e-06) {
  eigenvalues <- eigen(
    matrix_of_interest,
    symmetric = TRUE,
    only.values = TRUE
  )[["values"]]
  
  return(all(eigenvalues >= -tolerance * abs(eigenvalues[1]))) # 1st is the biggest eigenvalue
}

T_function <- function(x, lambda){
  return(sign(x)*max(abs(x)-lambda,0))
}


off_diagonal_part <- function(my_matrix){
  (sum(abs(my_matrix) > 0.0001) - p) / (p * (p-1))
}

frob_norm <- function(matrix1, matrix2){
  diff_matrix <- matrix1 - matrix2
  sum(diff_matrix * diff_matrix)
}


