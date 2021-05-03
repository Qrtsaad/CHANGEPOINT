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
## First plot
df <- data.frame(seq,resOP2,resPELT2)

ggplot(df, aes(seq)) + # basic graphical object
  geom_line(aes(y=resOP2), colour="blue") +  # first layer
  geom_line(aes(y=resPELT2), colour="red") + # second layer
  ylab("Time") +
  xlab("n")

## Second plot
dfb <- data.frame(log(seq),log(resOP2),log(resPELT2))

ggplot(dfb, aes(log(seq))) + # basic graphical object
  geom_line(aes(y=log(resOP2)), colour="blue") +  # first layer
  geom_line(aes(y=log(resPELT2)), colour="red") + # second layer
  ylab("log(Time)") +
  xlab("log(n)")

