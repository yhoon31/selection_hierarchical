#BH procedure
BH <- function(p_values, alpha) 
{
  ##inputs
  #p_values : vector of p values
  #alpha : target level
  
  ##output
  #rejection threshold
  
  p_sorted <- sort(p_values)
  m <- length(p_values)
  if(min(p_sorted*m/(1:m)) <= alpha)
  {
    threshold <- p_sorted[max(which(p_sorted <= (1:m) / m * alpha))]
  }
  else
  {
    threshold <- 0
  }
  return(threshold)
  
}

#e-BH procedure
eBH <- function(e_values,alpha)
{
  ##inputs
  #e_values : vector of e values
  #alpha : target level
  
  ##output
  #rejection threshold
  
  m <- length(e_values)
  e_sorted <- sort(e_values,decreasing=T)
  e_adjusted <- (1:m)*e_sorted/m
  if(max(e_adjusted) < 1/alpha)
  {
    threshold <- Inf
  }
  else
  {
    threshold <- e_sorted[max(which(e_adjusted >= 1/alpha))]
  }
  return(threshold)
}


#function for generating group sizes
f_N_const <- function(K,N)
{
  return(rep(K,N))
}

f_N_pois <- function(K,lambda)
{
  return(2+rpois(K,lambda))
}