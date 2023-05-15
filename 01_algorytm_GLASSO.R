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

my_GLASSO <- function(my_cov, lambda){
  stopifnot(lambda > 0)
  stopifnot(is.positive.semi.definite.matrix(my_cov))
  
  p <- ncol(my_cov) # wymiar macierzy
  
  sigma_matrix <- my_cov + diag(lambda, nrow = p) # macierz precyzji
  # Diagonal are already ok
  B <- matrix(numeric(p*p), ncol = p) # macierz zer
  
  continue_outer <- TRUE
  v <- p
  while (continue_outer){ # this will stop if no changes were made to sigma_matrix
    continue_outer <- FALSE # this will be changed to TRUE if anything will change in B or sigma
    v <- ifelse(v == p, 1, v + 1)
    
    continue_inner <- TRUE
    u <- p
    counter_changes_B <- -1 # will be changed in the 1st run of the while
    
    while (continue_inner){
      if (u == p){
        if(counter_changes_B == 0){
          continue_inner <- FALSE # zbieglismy juz
          next
        }
        
        u <- 1
        counter_changes_B <- 0
      }
      else{
        u <- u + 1
      }
      if (u == v){
        next
      }
      
      old_B_u_v <- B[u,v]
      my_sum <- sum(sigma_matrix[u,-v] * B[-v,v])
      B[u,v] <- T_function(my_cov[u,v] - my_sum, lambda)/sigma_matrix[v,v]
      
      if(abs(old_B_u_v - B[u,v]) > 0.001){
        counter_changes_B <- counter_changes_B + 1
        continue_outer <- TRUE
      }
    }
    
    continue_inner <- TRUE
    counter_changes_sigma <- 1
    while (continue_inner){
      if (u == p){
        if(counter_changes_sigma == 0){
          continue_inner <- FALSE # zbieglismy juz
          next
        }
        
        u <- 1
        counter_changes_sigma <- 0
      }
      else{
        u <- u + 1
      }
      if (u == v){
        next
      }
      
      old_sigma_u_v <- sigma_matrix[u,v]
      sigma_matrix[u,v] <- sum(sigma_matrix[u,-v] * B[-v,v])
      
      if(abs(old_sigma_u_v - sigma_matrix[u,v]) > 0.001){
        counter_changes_sigma <- counter_changes_sigma + 1
        continue_outer <- TRUE
      }
    }
  }
  
  # Teraz budujemy macierz K.
  my_K <- matrix(numeric(p*p), nrow=p) #inicjalizujemy
  
  for (v in 1:p) {
    my_K[v,v] <- 1/(sigma_matrix[v,v] - sum(sigma_matrix[v,-v] * B[-v,v]))
  }
  
  for (v in 1:p) {
    for (u in 1:p) {
      if (u == v){
        next
      }
      else{
        my_K[u,v] <- -B[u,v] * my_K[v,v]
      }
    }
  }
  
  my_K
}


