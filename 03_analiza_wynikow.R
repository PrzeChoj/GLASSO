source("./01_algorytm_GLASSO.R")
source("./02_wczytanie_danych.R")

data1

# Generate data2:
set.seed(1234)
p <- 50 # 500 --->>> 25 seconds; 1000 --->>> 170 seconds
n <- floor(p*4/5)
data2 <- get_data2(n, p, plot_points = TRUE)





lambda <- 0.01
K <- z_papiera_GLASSO(cov(data2), lambda, verbose = FALSE)
(sum(abs(K) > 0.0001) - p) / (p * (p-1)) # Jak duzo rzeczy poza diagonala jest niezerowych

