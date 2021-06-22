qin <- function(theta, data, min, beta, i, n)
{
  #mi[i] et pas mi[i-1] car déja décalage d'index pour le vecteur mi
  return(min[i] + beta + (n-i+1)*theta^2 - 2*sum(data[i:n])*theta + sum(data[i:n]^2))
}



Qn <- function(theta, data, min, beta, n)
{
  if(n == 1)
  {
    return(qin(theta, data, min, beta, 1, 1))
  }
  else
  {
    mini <- qin(theta, data, min, beta, 1, n)
    for (i in 2:n)
    {
      val <- qin(theta,data,min,beta,i,n)
      if(val < mini)
      {
        mini <- val
      }
    }
    return (mini)
  }
}



plotminQuad <- function(theta,data,mi,beta,t)
{
  myQs <- NULL
  myMins <- NULL

  for (i in 1:t)
  {
    myQs <- cbind(myQs, Qn(theta,data,mi,beta,i))
  }

  for(i in 1:length(theta))
  {
    myMins[i] <- min(myQs[i,])
  }

  myMins

  ggplot(data.frame(x=theta),aes(x=theta)) + geom_line(aes(y = myMins))
}
