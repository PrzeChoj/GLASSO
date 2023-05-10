# Pierwszy zbior
if (!require("huge", quietly = TRUE)) {
  install.packages("huge")
}
library(huge)
data(stockdata)
x <- log(stockdata$data[2:1258, ] / stockdata$data[1:1257, ])
colnames(x) <- stockdata$info[, 1]
x




# Drugi zbior
if (!require("visNetwork", quietly = TRUE)) {
  install.packages("visNetwork")
}
library(visNetwork)




# Generate graph:
p <- 400 # 500 --->>> 25 seconds; 1000 --->>> 170 seconds
nodes <- data.frame(id = 1:p)

get_edges <- function(p){
  cordinates <- data.frame(matrix(runif(p*2), ncol = 2))
  plot(cordinates, type="p", pch=20)
  
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
        
        lines(c(cordinates[i,1], cordinates[j,1]), c(cordinates[i,2], cordinates[j,2]))
      }
    }
  }
  close(progressBar)
  
  edges <- data.frame(from = from, to = to)
  
  edges
}

edges <- get_edges(p)

# TODO(Pomyslec nad tym:)
# Nie jestem tego pewny z 2 powodow:
# 1. Napisałem *sqrt(p), bo tak mi pasowało, ale w PDF-ie było /sqrt(p)
# 2. Pisali, że robili jakies usuwanie, zeby byly 4 krawedzie, a ja po prostu
    # nie generowalem nastepnych, bo to by bylo bardzo czasochlonne.





