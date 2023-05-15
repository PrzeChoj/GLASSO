library(glmnet)

# Ponizej to implementacja z zajec. Niestety algorytm rozbiega
z_zajec_GLASSO <- function(my_cov, lambda, maxiter = 100){
  stopifnot(lambda > 0)
  stopifnot(is.positive.semi.definite.matrix(my_cov))
  
  p <- ncol(my_cov) # wymiar macierzy
  
  sigma_matrix <- my_cov + diag(lambda, nrow = p) # macierz precyzji
  # Diagonal are already ok
  B <- matrix(numeric(p*p), ncol = p) # macierz zer
  
  continue_outer <- TRUE
  counter_outer <- 0
  v <- p
  while (continue_outer){ # this will stop if no changes were made to sigma_matrix
    counter_outer <- counter_outer + 1
    continue_outer <- FALSE # this will be changed to TRUE if anything will change in B or sigma
    v <- ifelse(v == p, 1, v + 1)
    
    continue_inner <- TRUE
    counter_inner <- 0
    u <- p
    counter_changes_B <- -1 # will be changed in the 1st run of the while
    
    while (continue_inner){
      counter_inner <- counter_inner + 1
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
      B[v,u] <- B[u,v]
      
      ic(abs(old_B_u_v - B[u,v]))
      if(abs(old_B_u_v - B[u,v]) > 0.001){
        counter_changes_B <- counter_changes_B + 1
        continue_outer <- TRUE
      }
      
      if(counter_inner > maxiter){
        continue_inner <- FALSE
      }
    }
    
    continue_inner <- TRUE
    counter_inner <- 0
    counter_changes_sigma <- 1
    while (continue_inner){
      counter_inner <- counter_inner + 1
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
      sigma_matrix[v,u] <- sigma_matrix[u,v]
      
      if(abs(old_sigma_u_v - sigma_matrix[u,v]) > 0.001){
        counter_changes_sigma <- counter_changes_sigma + 1
        continue_outer <- TRUE
      }
      
      if(counter_inner > maxiter){
        continue_inner <- FALSE
      }
    }
    
    if(counter_outer > maxiter){
      continue_outer <- FALSE
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

# Ponizej to implementacja z papiera z uzyciem wbudowanej LASSO solver. Algorytm dziala zgodnie z oczekiwaniami
z_papiera_GLASSO <- function(my_cov, lambda, maxiter = 100, t = 0.001, verbose = FALSE){
  stopifnot(lambda > 0)
  stopifnot(is.positive.semi.definite.matrix(my_cov))
  
  stop_treshold <- t * mean(abs(my_cov[row(my_cov) != col(my_cov)]))
  
  
  p <- ncol(my_cov) # wymiar macierzy
  
  sigma_matrix <- my_cov + diag(lambda, nrow = p) # macierz precyzji
  # Diagonal are already ok
  B <- matrix(numeric(p*p), ncol = p) # macierz zer
  
  
  
  v <- p
  old_sigma_matrix <- sigma_matrix
  stop_crit <- FALSE
  num_iter <- 0
  while(!stop_crit && (num_iter < maxiter)){
    v <- ifelse(v == p, 1, v + 1)
    
    tryCatch(x <- expm::sqrtm(sigma_matrix[-v,-v]),
             error = function(e) {browser()}) # sometimes this surprisingly happened
    y <- Matrix::solve(x, my_cov[-v,v])
    
    # Rozwiaz LASSO i odczytaj wynik
    coefs <- unclass(coef(glmnet(x, y, lambda = lambda, intercept = FALSE)))
    my_beta <- numeric(p-1)
    for (i in 1:(p-1)) {
      if(i %in% attr(coefs, "i")){
        my_beta[i] <- attr(coefs, "x")[attr(coefs, "i") == i]
      }
    }
    
    B[-v, v] <- my_beta # rozwiazania trzymam w kolumnie
    
    sigma_matrix[-v,v] <- sigma_matrix[-v,-v] %*% my_beta
    sigma_matrix[v,-v] <- t(sigma_matrix[-v,v])
    
    if(v == p){
      num_iter <- num_iter + 1
      
      diff_sigma <- sum(abs(old_sigma_matrix - sigma_matrix))
      if(verbose){
        print(paste0("Roznica w sigmie w ostatniej iteracji: ", diff_sigma))
      }
      if(diff_sigma < stop_treshold){ # kryt stopu zgodne z papierem
        stop_crit <- TRUE
      }
      old_sigma_matrix <- sigma_matrix
    }
  }
  
  
  # Teraz budujemy macierz K.
  my_K <- matrix(numeric(p*p), nrow=p) # inicjalizujemy
  
  for (v in 1:p) {
    my_K[v,v] <- 1/(sigma_matrix[v,v] - sum(sigma_matrix[v,] * B[,v])) # B[v,v] == 0
    my_K[-v,v] <- -B[-v,v] * my_K[v,v]
    my_K[v,-v] <- my_K[-v,v] # na podstawie https://github.com/MGallow/GLASSOO/blob/master/src/glasso.cpp, komentarz w linii 141
  }
  
  my_K
}

