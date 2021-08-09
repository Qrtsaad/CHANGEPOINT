## Library
devtools::install_github("Qrtsaad/CHANGEPOINT", force = TRUE)
library(CHANGEPOINT)
library(pracma)
library(caret)
library(lattice)
library(ggplot2)

### Nombre d'itérations en fonction de la tolérance

myseq <- 0:15
log10seq <- 10^-myseq

matABC <- rbind(
  c(100,1,10),
  c(100,1,100),

  c(1,100,10),
  c(1,100,100),

  c(1000,1,100),
  c(1000,1,1000),
  c(1000,10,100),

  c(5000,1,1000),
  c(5000,10,5000),

  c(100000,10,10000),
  c(100000,10,100000),

  c(1000000,10,10000),
  c(1000000,10,100000),
  c(1000000,10,1000000)
)

for(t in 1:nrow(matABC))
{
  Newton <- NULL
  Secant <- NULL
  Muller <- NULL

  a <- matABC[t,1]
  b <- matABC[t,2]
  c <- matABC[t,3]

  for(i in log10seq)
  {
    Newton <- c(Newton, searchzero(A = a, B = b, C = c, tol = i, method = "newton")$iterations)
    Secant <- c(Secant, searchzero(A = a, B = b, C = c, tol = i, method = "secante")$iterations)
    Muller <- c(Muller, searchzero(A = a, B = b, C = c, tol = i, method = "muller")$iterations)
  }

  mymain <- paste("Comparaison des méthodes avec ", paste("A = ", a, "B = ", b, "C = ", c), sep = ": ")
  p <- xyplot(Newton + Secant + Muller ~ myseq, xlab = "tol (10^-n)", ylab = "it", main = mymain, type = "l", auto.key = list(points = F,lines = T), par.settings = list(superpose.line = list(col = c("black","red","blue"))))
  print(p)
}

### Cp en fonction de la tolérance
myseq <- 0:15
log10seq <- 10^-myseq

matABC <- rbind(
  c(100,1,10),
  c(100,1,100),

  c(1,100,10),
  c(1,100,100),

  c(1000,1,100),
  c(1000,1,1000),
  c(1000,10,100),

  c(5000,1,1000),
  c(5000,10,5000),

  c(100000,10,10000),
  c(100000,10,100000),

  c(1000000,10,10000),
  c(1000000,10,100000),
  c(1000000,10,1000000)
)

for(t in 1:nrow(matABC))
{
  Newton <- NULL
  Secant <- NULL
  Muller <- NULL

  a <- matABC[t,1]
  b <- matABC[t,2]
  c <- matABC[t,3]

  for(i in log10seq)
  {
    Newton <- c(Newton, searchzero(A = a, B = b, C = c, tol = i, method = "newton")$complexity)
    Secant <- c(Secant, searchzero(A = a, B = b, C = c, tol = i, method = "secante")$complexity)
    Muller <- c(Muller, searchzero(A = a, B = b, C = c, tol = i, method = "muller")$complexity)
  }

  mymain <- paste("Comparaison des méthodes avec ", paste("A = ", a, "B = ", b, "C = ", c), sep = ": ")
  p <- xyplot(Newton + Secant + Muller ~ myseq, xlab = "tol (10^-n)", ylab = "cp", main = mymain, type = "l", auto.key = list(points = F,lines = T), par.settings = list(superpose.line = list(col = c("black","red","blue"))))
  print(p)
}
