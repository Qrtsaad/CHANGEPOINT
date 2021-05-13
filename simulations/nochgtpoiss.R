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
resOP2 <- NULL
resPELT2 <- NULL

for(i in seq)
{
  print(i)
  data <- data_generator(n = i, type = "poisson")
  plot(data)
  tOP <- system.time(myOP(data, cost = "poisson"))[3]
  print(tOP)
  resOP2 <- c(resOP2, tOP)
  tPELT <- system.time(myPELT(data, cost = "poisson"))[3]
  print(tPELT)
  resPELT2 <- c(resPELT2, tPELT)
}


#### Affichage des vecteurs
seq
resOP2
resPELT2


#### Plot des graphes
n <- seq
OP <- resOPp
PELT <- resPELTp


## First plot
xyplot(OP + PELT ~ n, ylab = "Time", main = "OP vs PELT : Poisson model without changepoint", type = "l", auto.key = list(points = F,lines = T), par.settings = list(superpose.line = list(col = c("red","blue"))))

## Second plot
xyplot(log(OP) + log(PELT) ~ log(n), ylab = "log(Time)", main = "OP vs PELT : Poisson model without changepoint", type = "l", auto.key = list(points = F,lines = T), par.settings = list(superpose.line = list(col = c("red","blue"))))

