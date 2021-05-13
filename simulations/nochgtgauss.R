## Library
devtools::install_github("Qrtsaad/CHANGEPOINT", force = TRUE)
library(CHANGEPOINT)
library("parallel")
library(ggplot2)

##
cores <- detectCores()
cores

#### Initialisation de la séquence
seq <- seq(100,1000,100)


#### Calcul des coûts pour OP et PELT
resOP1 <- NULL
resPELT1 <- NULL

for(i in seq)
{
  print(i)
  data <- data_generator(n = i, type = "gauss")
  plot(data)
  tOP <- system.time(myOP(data, cost = "gauss"))[3]
  print(tOP)
  resOP1 <- c(resOP1, tOP)
  tPELT <- system.time(myPELT(data, cost = "gauss"))[3]
  print(tPELT)
  resPELT1 <- c(resPELT1, tPELT)
}


#### Affichage des vecteurs
seq
resOP1
resPELT1


#### Plot des graphes
n <- seq
OP <- resOPg
PELT <- resPELTg

## First plot
xyplot(OP + PELT ~ n, ylab = "Time", main = "OP vs PELT : Gaussian model without changepoint", type = "l", auto.key = list(points = F,lines = T), par.settings = list(superpose.line = list(col = c("red","blue"))))

## Second plot
xyplot(log(OP) + log(PELT) ~ log(n), ylab = "log(Time)", main = "OP vs PELT : Gaussian model without changepoint", type = "l", auto.key = list(points = F,lines = T), par.settings = list(superpose.line = list(col = c("red","blue"))))
