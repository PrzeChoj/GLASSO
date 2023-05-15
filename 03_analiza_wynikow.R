source("./01_algorytm_GLASSO.R")
#source("./02_wczytanie_danych.R")

set.seed(234567)
p <- 50
lambda <- 0.1
my_cov <- matrix(rnorm(p*p), ncol=p)
my_cov <- my_cov %*% t(my_cov)
my_cov

K <- z_papiera_GLASSO(my_cov,lambda)

K