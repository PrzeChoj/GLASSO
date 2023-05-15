source("./01_algorytm_GLASSO.R")
source("./02_wczytanie_danych.R")

data1
# Generate graph:
p <- 40 # 500 --->>> 25 seconds; 1000 --->>> 170 seconds
data2 <- get_data2(p, plot_points = TRUE)



set.seed(234567)
p <- 5
lambda <- 0.3
my_cov <- matrix(rnorm(p*p), ncol=p)
my_cov <- my_cov %*% t(my_cov)
my_cov
solve(my_cov)

K <- z_papiera_GLASSO(my_cov, lambda)

K