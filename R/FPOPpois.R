#' Functional Pruned Optimal Partitioning 1D poisson
#'
#'
#' @description Functional Pruned Optimal Partitioning 1D Algorithm for poisson cost function
#'
#' @param data vector of data points
#' @param cost a number
#' @param beta a number
#'
#' @return a vector of changepoints, global cost
#' @export
#'
#' @examples
#' downupFPOPpois(c(rpois(10, lambda = 1), rpois(10, lambda = 50)))
downupFPOPpois <- function(data, beta = best_beta(data), affiche = FALSE)
{

  sizev <- 0

  n <- length(data)
  tau <- rep(0, n)
  tau[1] <- 1
  mi <- rep(0, n + 1)
  mi[1] <- - beta

  #### vecteurs preprocessing
  csd <- cumsum(data)
  csd2 <- cumsum(data^2)
  csd <- c(0, csd)
  csd2 <- c(0, csd2)

  ## pour avoir mean(data[1:t]) on fait csd[t+1]/t
  ## pour avoir ((t+1)-i+1)V_{i:t+1} on fait csd2[t+2]-csd2[i] - (csd[t+2]-csd[i])^2/((t+1)-i+1)

  #######

  v <- matrix(nrow = 0, ncol = 4)
  colnames(v) <- c("label1", "label2", "value", "position")

  #1er élément
  v <- rbind(v, c(1, 1, 0, data[1]))

  for (t in 1:(n - 1))
  {
    #####
    # STEP 1 ADD NEW DATA
    #####
    v <- rbind(v, c(t + 1, t + 1, v[1, 3] + beta, data[t + 1]))

    #####
    # STEP 2 UPDATE POINTS (values)
    #####
    for (i in 1:(nrow(v)-1)) ### CORRIGE (enlever le dernier ajouté à l'instant)
    {
      if (v[i,1] == v[i,2])
      {
        s <- v[i,1] #### CORRIGE
        mut <- (csd[t+2] - csd[s]) / ((t+1)-s+1)
        v[i,3] <- mi[s] + beta + cost_poiss(data[s:t+1])
        v[i,4] <- mut
      }
      else
      {
        v[i,3] <- v[i,3] + v[i,4] - data[t+1]*log(v[i,4])
      }
    }

    #####
    # STEP 3 : NEW INTERSECTIONS
    #####
    indices <- unique(v[,1]) #### indices = indices des quadratiques présentes
    nb_old_indices <- length(indices) - 1

    for (i in 1:nb_old_indices) #new index = t+1
    {
      j <- indices[i]

      A <- -(t+1 - j - 1)
      B <- sum(data[j:t+1])
      C <- -sum(log(factorial(data[j:t+1])))
      delta <- (B/A)*exp(-C/A)
      m_ti <- mi[j] + beta

      if(is.nan(delta) == FALSE)
      {
        Wdeltap <- lambertWp(delta)
        Wdeltan <- lambertWn(delta)

        if(is.nan(Wdeltap) == FALSE)
        {
          thetap <- exp(-Wdeltap - (C/A))
          sgammap <- sum(thetap - data[j:t+1]*log(thetap))
          qp <- m_ti + sgammap
          v <- rbind(v, c(j, t+1, qp, thetap))
        }

        if(is.nan(Wdeltan) == FALSE)
        {
          thetan <- exp(-Wdeltan - (C/A))
          sgamman <- sum(thetan - data[j:t+1]*log(thetan))
          qn <- m_ti + sgamman
          v <- rbind(v, c(t+1, j, qn, thetan))
        }
      }
    }

    #####
    # STEP 4 : SORT BY VALUES
    #####
    v <- v[order(v[,3]),]

    #####
    # STEP 5 : PRUNING
    #####

    # SECOND :
    #v <- v[v[,3] < v[1,3] + beta,]


    # FIRST :

    I <- v[1,1]
    i <- 2
    while(i<nrow(v))
    {
      l <- v[i,1]
      p <- v[i,4]

      qI <- NULL
      for (j in I)
      {
        qj <- mi[j] + beta + sum(p - data[j:t+1]*log(p))
        qI <- c(qI,qj)
      }

      ql <- mi[l] + beta + sum(p - data[l:t+1]*log(p))

      if((is.element('TRUE',qI < ql) == TRUE) & (v[i,1]!=v[i,2]))
      {
        v <- v[-i,]
      }
      else
      {
        I <- c(I,v[i,1],v[i,2])
        I <- unique(I)
        i <- i + 1
      }
    }


    # THIRD :


    lab <- v[v[,1] == v[,2],]
    lab <- lab[is.element(lab[,1],I),] #on garde que les lignes du type l,l,v,p avec l \in I
    inter <- v[v[,1] != v[,2],] #partie de v avec les lignes du type l,l',v,p
    v <- rbind(lab,inter) #on associe les deux
    rownames(v) <- NULL




    mi[t+2] <- v[1,3]
    tau[t+1] <- v[1,1]

    if (affiche == TRUE) # pour vérifier les valeurs à chaque itération d'un test simple
    {
      print(v)
    }

    sizev <- c(sizev, nrow(v))

  }


  ########
  # BACKTRACKING
  ########

  s <- tau[n]
  P <- tau[n]

  while (s>1)
  {
    P <- c(P, tau[s-1])
    s <- tau[s-1]
  }

  P <- rev(P)[-1] - 1

  P

  return(list(vecofsizev = sizev, changepoints = P, globalCost = v[1,3] - length(P)*beta))
}
