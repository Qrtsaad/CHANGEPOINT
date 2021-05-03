devtools::install_github("Qrtsaad/CHANGEPOINT")
library(CHANGEPOINT)

best_beta(c(0,1,1,2,2,3))

cost_gauss(c(0,0,1,1,2,2,3,4,5))
cost_poiss(c(1,2,3))

data_generator(20, chpts = c(10), means = c(0,10), type = "gauss")
data_generator(25, chpts = c(10,20), means = c(2,0,1), type = "gauss")


myOP(c(0,0,1,1,0,0,0), beta = 0.00001)
myOP(c(rnorm(50, mean = 0, sd = 0.1), rnorm(50, mean = 10, sd = 0.1)))
myOP(data_generator(25, chpts = c(10,20), means = c(20,0,20), type = "gauss"), beta = 5)

myPELT(c(0,0,1,1,0,0,0), beta = 0.00001)
myPELT(c(rnorm(50, mean = 0, sd = 0.1), rnorm(50, mean = 10, sd = 0.1)))
myPELT(data_generator(25, chpts = c(10,20), means = c(20,0,20), type = "gauss"), beta = 5)
