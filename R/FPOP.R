#' Functional Pruned Optimal Partitioning 1D
#'
#'
#' @description Functional Pruned Optimal Partitioning 1D Algorithm
#'
#' @param data vector of data points
#' @param cost a number
#' @param beta a number
#'
#' @return a vector of changepoints, global cost
#' @export
#'
#' @examples
#' myFPOP1D(c(0,0,1,1,0,0,0), beta = 0.00001)
#' myFPOP1D(c(rnorm(50, mean = 0, sd = 0), rnorm(50, mean = 10, sd = 0)), beta = 1)
myFPOP <- function(data, beta = best_beta(data))
{
  n <- length(data)


  tau <- rep(0,n)
  tau[1] <- 1


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
  #exemples d'utilisation :
  #(t-j+1)*Vit(data,j,t) remplacé par (csd2[t] - csd2[j-1]) - (csd[t] - csd[j-1])^2/(t-j+1)
  #Vit(data,j,t) remplacé par (csd2[t] - csd2[j-1])/(t-j+1) - ((csd[t]-csd[j-1])/(t-j+1))^2


  v <- matrix(nrow = 0, ncol = 4)
  colnames(v) <- c("label1","label2","value","position")
  v <- rbind(v, c(1,1,0,data[1]))


  for (t in 2:n)
  {
    # line 5
    v <- rbind(v, c(t, t, v[1,3] + beta, data[t]))

    #ligne 7 update
    for (i in 2:t)
    {
      if (v[i,1] == v[i,2])
      {
        #v[i,3] <- v[1,3] + beta + (t-i+1)*Vit(data,i,t)
        v[i,3] <- v[1,3] + beta + (csd2[t] - csd2[i-1]) - (csd[t] - csd[i-1])^2/(t-i+1)
        v[i,4] <- mean(data[i:t])
      }
      else if (v[i,1] != v[i,2])
      {
        v[,3] <- v[,3] + (data[t] - v[,4])^2
      }
    }



    # ligne 9 => new elements (intersection quadratique(new)/quadratique(old))

    L <- unique(v[,1])
    lenL <- length(L) - 1


    for (i in 1:lenL)
    {
      j <- L[i]


      if (j == 1)
      {
        Vjt <- (csd2[t] - csd2[1])/t - ((csd[t]-csd[1])/t)^2
        Vjtm1 <- (csd2[t-1] - csd2[1])/t - ((csd[t-1]-csd[1])/t)^2
      }

      else
      {
        Vjt <- (csd2[t] - csd2[j-1])/(t-j+1) - ((csd[t]-csd[j-1])/(t-j+1))^2
        Vjtm1 <- (csd2[t-1] - csd2[j-1])/(t-j) - ((csd[t-1]-csd[j-1])/(t-j))^2
      }


      #old version
      #mjt <- (t-j+1)*Vit(data,j,t) + mi[j-1] + beta
      #new version mi[j] = m_{j-1}
      #mjt <- (t-j+1)*Vit(data,j,t) + mi[j] + beta
      # on remplace Vit
      #Vjt <- (csd2[t] - csd2[j-1])/(t-j+1) - ((csd[t]-csd[j-1])/(t-j+1))^2
      mjt <- (t-j+1)*Vjt + mi[j] + beta


      # Delta
      #delta <- abs(mean(data[j:t]) - mean(data[j:(t+1)]))
      delta <- abs(mean(data[j:(t-1)]) - mean(data[j:t]))

      # R
      #old version
      #R <- (mi[t-1] - mi[j-1])/(t-j) - Vit(data,j,t-1)
      #new version mi[t] = m_{t-1} et mi[j] = m_{j-1}
      #R <- (mi[t] - mi[j])/(t-j) - Vit(data,j,t-1)
      # on remplace Vit
      #Vjtm1 <- (csd2[t-1] - csd2[j-1])/(t-1-j+1) - ((csd[t-1]-csd[j-1])/(t-1-j+1))^2
      R <- (mi[t] - mi[j])/(t-j) - Vjtm1


      q1 <- (t-j+1)*(delta - R)^2 + mjt
      q2 <- (t-j+1)*(delta + R)^2 + mjt
      theta1 <- mean(data[j:(t-1)]) - R
      theta2 <- mean(data[j:(t-1)]) + R


      v <- rbind(v,c(j,t,q1,theta1))
      v <- rbind(v,c(t,j,q2,theta2))

    }


    #ligne 10 (tri selon les valeurs)
    v <- v[order(v[,3]),]

    I <- v[1,1]
    i <- 2
    while (i < length(v[,1]) & v[i,3] < v[1,1] + beta)
    {
      l <- v[i,1]
      p <- v[i,4]
      ql <- qin(p,data,mi,beta,l,t)
      print("ql")
      print(ql)

      if ((qin(p,data,mi,beta,I,t) < ql) & (v[i,1] != v[i,2]))
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











    tau[t] <- v[1,1]
    #old version
    #mi[t] <- v[1,3]
    #new version
    mi[t+1] <- v[1,3]

    #PRINT MATRIX
    #afficher matrice pour voir
    #print(v)

    #points
    #plot(v[,4],v[,3], xlab = "position", ylab = "value")
    #abline(v = data[t], col = "lightgray")
    #par(new = FALSE)

  }


  #tau => backtracking -> vérifier qu'on a ce qu'on veut.
  s <- tau[n]
  P <- tau[n]

  while (s>1)
  {
    P <- c(P, tau[s-1])
    s <- tau[s-1]
  }
  P <- rev(P)[-1] - 1 #on fait -1 pour gerer les indices sinon ça donne les resultats des commentaires ci-dessous

  P


  # affichage des quadratiques
  #for (t in P)
  #{
  #    X <- v[,4]
  #    for (i in 1:t)
  #    {
  #      X <- cbind(X,Qn(v[,4],data,mi,0.5,i))
  #    }



  # on convertit le tableau en dataframe
  #    X <- data.frame(X)


  # on utilise la fonction melt pour pouvoir afficher plusieurs courbes
  #    X.melted = melt(X, id='X')

  # on affiche les courbes
  #    p <- ggplot(data = X.melted, aes(x = X, y = value, color = variable)) +
  #      geom_line()
  #    print(p)
  #}

  #affichage du min des quadratiques
  #for (i in P)
  #{
  #  print(plotminQuad(v[,4],data,mi,beta,i))
  #}

  return(list(min = mi, tau = tau, changepoints = P, vec = v, globalCost = as.numeric(v[1,3]) - length(P)*beta))
}
