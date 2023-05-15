# Pierwszy zbior
if (!require("huge", quietly = TRUE)) {
  install.packages("huge")
}
library(huge)
data(stockdata)
data1 <- log(stockdata$data[2:1258, ] / stockdata$data[1:1257, ])
colnames(data1) <- stockdata$info[, 1]




# Drugi zbior
if (!require("visNetwork", quietly = TRUE)) {
  install.packages("visNetwork")
}
library(visNetwork)

get_edges <- function(p, plot_points = TRUE){
  cordinates <- data.frame(matrix(runif(p*2), ncol = 2))
  if(plot_points){
    plot(cordinates, type="p", pch=20)
  }
  
  from <- numeric(0)
  to <- numeric(0)
  
  num_of_edges <- numeric(p)
  
  progressBar_counter <- 1
  progressBar <- utils::txtProgressBar(min = 0, max = (p-1)*(p-1)/2, initial = 1)
  for (i in 1:(p-1)) {
    progressBar_counter <- progressBar_counter + (p-i)
    utils::setTxtProgressBar(progressBar, progressBar_counter)
    for (j in (i+1):p){
      if (num_of_edges[i] == 4 || num_of_edges[j] == 4){
        next
      }
      
      if (runif(1) < dnorm(dist(cordinates[c(i,j),])*sqrt(p))){
        num_of_edges[i] <- num_of_edges[i] + 1
        num_of_edges[j] <- num_of_edges[j] + 1
        from <- c(from, i)
        to <- c(to, j)
        
        if(plot_points){
          lines(c(cordinates[i,1], cordinates[j,1]), c(cordinates[i,2], cordinates[j,2]))
        }
      }
    }
  }
  close(progressBar)
  
  data.frame(from = from, to = to)
}

macierz_kowariancji <- function(edges_list, p) {
  M <- matrix(0, nrow = p, ncol = p)
  diag(M) <- rep(1, p)
  edges <- t(edges_list)
  n_edges <- length(edges) / 2
  
  for (i in 1:n_edges) {
    M[edges[1, i], edges[2, i]] <- 0.245
    M[edges[2, i], edges[1, i]] <- 0.245
  }
  
  cov2cor(solve(M))
}


get_data2 <- function(n, p, plot_points = TRUE){
  nodes <- data.frame(id = 1:p)
  edges <- get_edges(p, plot_points = plot_points)
  
  sigma_matrix <- macierz_kowariancji(edges, p)
  
  MASS::mvrnorm(n, mu = numeric(p), Sigma = sigma_matrix)
}







