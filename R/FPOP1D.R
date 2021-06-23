#' Functional Pruned Optimal Partitioning 1D v2
#'
#' @description Functional Pruned Optimal Partitioning 1D Algorithm
#'
#' @param data vector of data points
#' @param cost a number
#' @param beta a number
#'
#' @return a vector of changepoints, global cost, a matrix with labels, values and positions, a vector of tau
#' @export
#'
#' @examples
#' myFPOP1D(c(0,0,1,1,0,0,0), beta = 0.00001)
#' myFPOP1D(c(rnorm(50, mean = 0, sd = 0), rnorm(50, mean = 10, sd = 0)), beta = 1)
myFPOP1Dv2 <- function(data, beta = best_beta(data))
{
  n <- length(data)
  tau <- rep(0,n)
  tau[1] <- 1 #si on met s > 1 sur le backtracking on obtient changepoint = numeric(0)
  #tau[1] <- 0 # résoud un probleme mais ?


  #### min
  #old version : pb d'indices avec m[i-1] for i in 1:length(L)
  #mi <- rep(0,n) #vecteur qui stocke le min de chaque itération
  #mi[1] <- 0
  #new version : on créé un vecteur de taille n+1 avec mi[1] = -beta qui est en réalité m_0 => décalagé d'indice facile à gerer (voir autres modifs dans étape 9 et ajout des valeurs de mi[t] qui devient mi[t+1])
  #rq : ici on donne la valeur de m_0 = - beta "à la main" et on donne m_1 = 0 via rep(0,n+1)
  mi <- rep(0,n+1)
  mi[1] <- -beta


  #### vecteurs preprocessing
  csd <- cumsum(data) #vecteur qui stocke la somme des valeurs
  csd2 <- cumsum(data^2) #vecteur qui stocke la somme du carré des valeurs

  #######

  v <- matrix(nrow = 0, ncol = 4)
  colnames(v) <- c("label1","label2","value","position")
  v <- rbind(v, c(1,1,0,data[1]))

  for (t in 1:(n-1))
  {
    # STEP 1 "FIRST" NEW DATA
    v <- rbind(v, c(t+1, t+1, mi[t+1] + beta, data[t+1]))

    # STEP 2 UPDATE POINTS
    for (i in 1:t)
    {
      if (v[i,1] == v[i,2])
      {
        if(i == 1)
        {
          Vart <- csd2[t+1] - (csd[t+1])^2/(t+1)
          mut <- csd[t+1]/(t+1)
        }
        else
        {
          Vart <- (csd2[t+1] - csd2[i-1]) - (csd[t+1] - csd[i-1])^2/(t+1-i+1)
          mut <- (csd[t+1] - csd[i-1])/(t+1-i+1)
        }

        v[i,3] <- mi[i] + beta + Vart
        v[i,4] <- mut
      }
      else if (v[i,1] != v[i,2])
      {
        v[i,3] <- v[i,3] + (data[t+1] - v[i,4])^2
      }
    }

    # STEP 3 : NEW QUAD

    L <- unique(v[,1])
    lenL <- length(L) - 1
    #print("L")
    #print(L)

    for (i in 1:lenL)
    {
      j <- L[i]

      if(j == 1)
      {
        Vjtp1 <- csd2[t+1] - (csd[t+1])^2/(t+1)
        VR <- csd2[t]/t - (csd[t]/t)^2
        mujt <- (csd[t])/t
        mujtp1 <- (csd[t+1])/(t+1)
      }

      else if (j != 1)
      {
        mujt <- (csd[t]-csd[j-1])/(t-j+1)
        mujtp1 <- (csd[t+1]-csd[j-1])/(t+1-j+1)
        Vjtp1 <- (csd2[t+1] - csd2[j-1]) - (csd[t+1] - csd[j-1])^2/(t+1-j+1)

        if (j<t)
        {
          VR <- (csd2[t] - csd2[j-1])/(t-j+1) - (csd[t] - csd[j-1])^2
          #VR <- (csd2[t] - csd2[j-1])/(t-j+1) - ((csd[t] - csd[j-1])/(t-j+1))^2
        }

        else
        {
          VR <- (csd2[j-1] - csd2[t-1])/(j-(t+1)) - (csd[j-1] - csd[t-1])^2
        }
      }

      mjtp1 <- mi[j] + beta + Vjtp1

      delta <- abs(mujt - mujtp1)

      R2 <- (mi[t+1] - mi[j])/(t+1-j) - VR
      R <- sqrt(abs(R2))
      #R <- (mi[t+1] - mi[j])/(t+1-j) - VR
      #print(R)

      q1 <- (t+1-j+1)*(delta - R)^2 + mjtp1
      theta1 <- mujt - R
      q2 <- (t+1-j+1)*(delta + R)^2 + mjtp1
      theta2 <- mujt + R

      v <- rbind(v,c(j,t+1,q1,theta1))
      v <- rbind(v,c(t+1,j,q2,theta2))
    }



    # STEP 4 : SORT BY VALUES
    v <- v[order(v[,3]),]

    # STEP 5 : PRUNINGS

    # FIRST

    #I <- v[1,1]
    #i <- 2
    #while (i < length(v[,1]) & v[i,3] < v[1,1] + beta)
    #{
    #  l <- v[i,1]
    #  p <- v[i,4]
    #  ql <- qin(p,data,mi,beta,l,t+1)
    #  print("ql")
    #  print(ql)

    #  for (j in I)
    #  {
    #    qj <- qin(p,data,mi,beta,j,t+1)
    #    if (qj < ql)
    #    {
    #      if (v[i,1] != v[i,2])
    #      {
    #        v <- v[-i,]
    #      }

    #    }

    #    else
    #    {
    #      I <- unique(I,v[i,1],v[i,2])
    #      i <- i+1
    #    }
    #  }
    #}
    # => pb boucle infinie, comment bien gerer si il existe j \in I tq ...

    I <- v[1,1]
    i <- 2
    while (i < length(v[,1]) & v[i,3] < v[1,1] + beta)
    {
      l <- v[i,1]
      p <- v[i,4]
      ql <- qin(p,data,mi,beta,l,t+1)
      #print("ql")
      #print(ql)

      if ((qin(p,data,mi,beta,I,t+1) < ql) & (v[i,1] != v[i,2]))
      {
        v <- v[-i,]
      }

      else
      {
        I <- unique(I,v[i,1],v[i,2])
        i <- i+1
      }
    }

    # SECOND

    #    for (i in 1:length(v[1,]))
    #    {
    #      if (v[i,3] > v[1,3] + beta)
    #      {
    #        v <- v[-i,]
    #      }
    #    }

    #=> PB: Pour un beta petit pruning enlève tout

    v <- v[v[,3] < v[1,3] + beta,]
    

    # THIRD

    for (i in v[,1])
    {
      if (v[i,1] == v[i,2])
      {
        if ((i %in% I) == FALSE)
        {
          v <- v[-i,]
        }
      }
    }



    #PRINT MATRIX
    #afficher matrice pour voir
    #print(v)

    tau[t+1] <- v[1,1]

    #new version
    mi[t+2] <- v[1,3]



  }


  # affichage des quadratiques
  #for (t in 1:length(data))
  #{
  #  X <- v[,4]
  #  for (i in 1:t)
  #  {
  #    X <- cbind(X,Qn(v[,4],data,mi,0.5,i))
  #  }



    # on convertit le tableau en dataframe
  #  X <- data.frame(X)


    # on utilise la fonction melt pour pouvoir afficher plusieurs courbes
  #  X.melted = melt(X, id='X')

    # on affiche les courbes
  #  p <- ggplot(data = X.melted, aes(x = X, y = value, color = variable)) + geom_line()
  #  print(p)
  #}


  #affichage du min des quadratiques
  #for (i in 1:length(data))
  #{
  #  print(plotminQuad(v[,4],data,mi,beta,i))
  #}

  #tau => backtracking -> vérifier qu'on a ce qu'on veut.
  s <- tau[n]
  P <- tau[n]

  while (s>1)
  {
    P <- c(P, tau[s-1])
    s <- tau[s-1]
  }

  P <- rev(P)[-1]
  P <- append(P,1,0)

  #P
  

  return(list(tau = tau, changepoints = P, vec = v, globalCost = mi[n+1] - length(P)*beta))
}
