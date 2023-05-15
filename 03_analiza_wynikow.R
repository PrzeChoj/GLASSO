source("./01_algorytm_GLASSO.R")
source("./02_wczytanie_danych.R")

data1

# Generate data2:
set.seed(1234)
p <- 50 # 500 --->>> 25 seconds; 1000 --->>> 170 seconds
n <- floor(p*4/5)
data2 <- get_data2(n, p, plot_points = TRUE)

data2_true_cov <- data2[[2]]
data2 <- data2[[1]]

# Jak duzo rzeczy poza diagonala jest niezerowych
off_diagonal_part(data2_true_cov) # 0.17
frob_norm(cov(data2), data2_true_cov) # 64


set.seed(1234)
lambda <- 0.03
# install.packages("icecream")
K <- z_papiera_GLASSO(cov(data2), lambda, verbose = TRUE) # verbose = TRUE wymaga pakietu icecream
off_diagonal_part(K)

frob_norm(solve(K), data2_true_cov) # 10, czyli duzo mniej niz 64


gout <- glasso::glasso(cov(data2), lambda)
frob_norm(gout$w, data2_true_cov) # 51, czyli troche mniej niz 64




set.seed(1234)
lambda <- 0.03

data1_cov <- cov(data1)
#K <- z_papiera_GLASSO(data1_cov, lambda, verbose = TRUE) # niestety baaaaardzo dlugo sie liczy... :<

gout <- glasso::glasso(data1_cov, lambda)




