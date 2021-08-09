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
#' downupFPOP(c(rnorm(50, mean = 0, sd = 1), rnorm(50, mean = 100, sd = 1)))
downupFPOP <- function(data, cost = "gauss", beta = best_beta(data), affiche = FALSE)
{
  allowed.cost <- c("gauss", "poisson", "negbin")
  if(!cost %in% allowed.cost){stop('type must be one of: ', paste(allowed.cost, collapse=", "))}

  if (cost == "gauss") {cost_f <- cost_gauss}
  else if (cost == "poisson") {cost_f <- cost_poiss}
  else if (cost == "negbin") {cost_f <- cost_negbin}

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
        Vart <- (csd2[t+2] - csd2[s]) - (csd[t+2] - csd[s])^2 / ((t+1)-s+1)
        mut <- (csd[t+2] - csd[s]) / ((t+1)-s+1)
        v[i,3] <- mi[s] + beta + Vart
        v[i,4] <- mut
      }
      else
      {
        v[i,3] <- v[i,3] + (data[t+1] - v[i,4])^2
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

      mujt <- (csd[t+1] - csd[j])/(t-j+1)
      mujtp1 <- (csd[t+2] - csd[j])/((t+1)-j+1)
      delta <- abs(mujt - mujtp1)

      V_itm1 <- (csd2[t+1] - csd2[j])/(t+1-j+1) - (csd[t+1] - csd[j])^2/(t+1-j+1)^2

      R2 <- (mi[t+1] - mi[j])/(t+1-j) - V_itm1 ###ALWAYS > 0

      if (R2 > 0)
      {
        R <- sqrt(R2)

        theta1 <- (csd[t+1] - csd[j])/(t-j+1) - R
        theta2 <- (csd[t+1] - csd[j])/(t-j+1) + R

        m_ti <- mi[j] + beta + (csd2[t+2] - csd2[j] - (csd[t+2] - csd[j])^2 / (t+1-j+1))

        q1 <- (t+1-j+1)*(delta - R)^2 + m_ti
        q2 <- (t+1-j+1)*(delta + R)^2 + m_ti

        v <- rbind(v, c(j, t+1, q1, theta1))
        v <- rbind(v, c(t+1, j, q2, theta2))
      }
    }



    #####
    # STEP 5 : SORT BY VALUES
    #####
    v <- v[order(v[,3]),]

    #####
    # STEP 6 : PRUNING
    #####

    # SECOND :

    v <- v[v[,3] < v[1,3] + beta,]

    # FIRST :

    I <- v[1,1]
    i <- 2
    while(i<nrow(v))
    {
      l <- v[i,1]
      p <- v[i,4]

      qI <- NULL
      #qIbis <- NULL
      for (j in I)
      {
        Vjtp1 <- (csd2[t+2] - csd2[j])/(t+1-j+1) - (csd[t+2] - csd[j])^2/(t+1-j+1)^2
        qj <- mi[j] + beta + (t+1-j+1)*((p - (csd[t+2]-csd[j])/(t+1-j+1) )^2 + Vjtp1)

        qI <- c(qI,qj)
      }

      Vltp1 <- (csd2[t+2] - csd2[l])/(t+1-l+1) - (csd[t+2] - csd[l])^2/(t+1-l+1)^2
      ql <- mi[l] + beta + (t+1-l+1)*((p - (csd[t+2]-csd[l])/(t+1-l+1) )^2 + Vltp1)

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
    #min <- min(v[,3])
    #print(min)
    tau[t+1] <- v[1,1]

    if (affiche == TRUE) # pour vérifier les valeurs à chaque itération d'un test simple
    {
      print(v)
    }

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

  return(list(vsize = nrow(v), changepoints = P, globalCost = v[1,3] - length(P)*beta))
}
