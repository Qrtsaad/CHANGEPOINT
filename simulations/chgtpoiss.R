## Library
devtools::install_github("Qrtsaad/CHANGEPOINT", force = TRUE)
library(CHANGEPOINT)
library("parallel")
library(ggplot2)


##
cores <- detectCores()
cores


### FUNCTIONS TO CREATE VECTORS FOR DATA_GENERATOR PARAMS


make_chpts <- function(n)
{
  if (n<0){stop('n must be non-negative')}
  else if (n <= 100){res <- seq(from = 0, to = n-1, by = n/2)}
  else if (n <= 200){res <- seq(from = 0, to = n-1, by = n/4)}
  else if (n <= 500){res <- seq(from = 0, to = n-1, by = n/5)}
  else if (n <= 1000){res <- seq(from = 0, to = n-1, by = n/10)}
  else {res <- seq(from = 0, to = n-1, by = n/20)}
  return (res)
}

make_means <- function(n)
{
  res <- NULL
  tmp <- 0
  for (i in 1:(n+1)){
    rand <- sample(0:10, 1)
    while (rand == tmp){rand <- sample(1:10, 1)}
    tmp <- rand
    res <- c(res,rand)
  }
  return (res)
}



#### Initialisation de la séquence
seq <- seq(100,1000,100)


#### Calcul des coûts pour OP et PELT
reschgtOP2 <- NULL
reschgtPELT2 <- NULL

for(i in seq)
{
  print(i)
  chgt <- make_chpts(i)
  print(chgt)
  moys <- make_means(length(chgt))
  data <- data_generator(n = i, chpts = chgt, means = moys, type = "poisson")
  plot(data)
  tOP <- system.time(myOP(data, cost = "poisson"))[3]
  print(tOP)
  reschgtOP2 <- c(reschgtOP2, tOP)
  tPELT <- system.time(myPELT(data, cost = "poisson"))[3]
  print(tPELT)
  reschgtPELT2 <- c(reschgtPELT2, tPELT)
}



#### Affichage des vecteurs
seq
reschgtOP2
reschgtPELT2


#### Plot des graphes
## First plot
df <- data.frame(seq,reschgtOP2,reschgtPELT2)

ggplot(df, aes(seq)) + # basic graphical object
  geom_line(aes(y=reschgtOP2), colour="blue") +  # first layer
  geom_line(aes(y=reschgtPELT2), colour="red") + # second layer
  ylab("Time") +
  xlab("n")

## Second plot
dfb <- data.frame(log(seq),log(reschgtOP2),log(reschgtPELT2))

ggplot(dfb, aes(log(seq))) + # basic graphical object
  geom_line(aes(y=log(reschgtOP2)), colour="blue") +  # first layer
  geom_line(aes(y=log(reschgtPELT2)), colour="red") + # second layer
  ylab("log(Time)") +
  xlab("log(n)")




library(lattice)
n <- seq
OP <- reschgtOP2
PELT <- reschgtPELT2
xyplot(OP + PELT ~ n, ylab = "Time", type = "l", auto.key = list(points = F,lines = T), par.settings = list(superpose.line = list(col = c("red","blue"))))

