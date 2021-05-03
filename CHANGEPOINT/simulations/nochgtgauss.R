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
## First plot
df <- data.frame(seq,resOP1,resPELT1)

ggplot(df, aes(seq)) + # basic graphical object
  geom_line(aes(y=resOP1), colour="blue") +  # first layer
  geom_line(aes(y=resPELT1), colour="red") + # second layer
  ylab("Time") +
  xlab("n")

## Second plot
dfb <- data.frame(log(seq),log(resOP1),log(resPELT1))

ggplot(dfb, aes(log(seq))) + # basic graphical object
  geom_line(aes(y=log(resOP1)), colour="blue") +  # first layer
  geom_line(aes(y=log(resPELT1)), colour="red") + # second layer
  ylab("log(Time)") +
  xlab("log(n)")

