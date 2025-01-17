source("functions.R")

##------------------------------------------------------------------------------

library(MASS)
library(randomForest)

p_G <- 10
p_X <- 20

K_train <- 100
K <- 200
K_test <- 20 

#data geneeration
set.seed(1234)
A <- matrix(rnorm(p_X*p_G),nrow=p_X)
beta_1 <- runif(p_X,0,1)/p_X

G_train <- matrix(runif(K_train*p_G,-5,5),nrow=K_train)
X_train <- list()
Y_train_1 <- list()
Y_train_0 <- list()
lambda <- 10 

X_train_mat <- matrix(nrow=0,ncol=p_X)
G_train_mat <- matrix(nrow=0,ncol=p_G)
N_k_train <-  f_N_pois(K_train,lambda)
for(i in 1:K_train)
{
  G_train_mat <- rbind(G_train_mat,rep(1,N_k_train[i])%*%t(G_train[i,]))
  X_train[[i]] <- matrix(runif(N_k_train[i]*p_X),nrow=N_k_train[i])
  X_train_mat <- rbind(X_train_mat,X_train[[i]])
  Y_train_1[[i]] <- sapply(1:N_k_train[i], function(j) rbeta(1,1+sum(beta_1*X_train[[i]][j,]),1-sum(beta_1*X_train[[i]][j,])))
  Y_train_0[[i]] <- sapply(1:N_k_train[i], function(j) rbeta(1,1-sum(beta_1*X_train[[i]][j,]),1+sum(beta_1*X_train[[i]][j,])))
  
}

#fit random forest regression
rf <- randomForest(X_train_mat,unlist(Y_train_1))

alpha_set <- seq(0.05,0.25,by=0.025)

#run procedure
tn <- 500
FDP_vec <- matrix(0,nrow=tn,ncol=length(alpha_set))
power_vec <- matrix(0,nrow=tn,ncol=length(alpha_set))
FDP_vec_U <- matrix(0,nrow=tn,ncol=length(alpha_set))
power_vec_U <- matrix(0,nrow=tn,ncol=length(alpha_set))
FDP_p_vec <- matrix(0,nrow=tn,ncol=length(alpha_set))
power_p_vec <- matrix(0,nrow=tn,ncol=length(alpha_set))


pb <- txtProgressBar(min = 0, max = tn, style = 3)

for(t in 1:tn)
{
  G <- matrix(runif(K*p_G,-5,5),nrow=K)
  G_test <- matrix(runif(K_test*p_G,-5,5),nrow=K_test)
  N_k <- f_N_pois(K,lambda)
  N_k_test <- f_N_pois(K_test,lambda)
  X <- list()
  Y_1 <- list()
  Y_0 <- list()
  X_test <- list()
  Y_test_1 <- list()
  Y_test_0 <- list()
  
  for(i in 1:K)
  {
    X[[i]] <- matrix(runif(N_k[i]*p_X),nrow=N_k[i])
    Y_1[[i]] <- sapply(1:N_k[i], function(j) rbeta(1,1+sum(beta_1*X[[i]][j,]),1-sum(beta_1*X[[i]][j,])))
  }
  for(i in 1:K_test)
  {
    X_test[[i]] <- matrix(runif(N_k_test[i]*p_X),nrow=N_k_test[i])
    Y_test_1[[i]] <- sapply(1:N_k_test[i], function(j) rbeta(1,1+sum(beta_1*X_test[[i]][j,]),1-sum(beta_1*X_test[[i]][j,])))
    Y_test_0[[i]] <- sapply(1:N_k_test[i], function(j) rbeta(1,1-sum(beta_1*X_test[[i]][j,]),1+sum(beta_1*X_test[[i]][j,])))
  }
  
  score_cal_1 <- lapply(1:K,function(i) Y_1[[i]]-predict(rf,X[[i]]))
  score_test_0 <- lapply(1:K_test,function(i) Y_test_0[[i]]-predict(rf,X_test[[i]]))
  
  ind <- sapply(1:K, function(k) sample(1:N_k[k],1))
  
  score_cal_1_sample <- sapply(1:K, function(k) score_cal_1[[k]][ind[k]])

  FDP_hat <- function(t,j)
  {
    numer <- sum((score_cal_1_sample < t))+1
    denom <- max(1,sum(sapply(score_test_0,function(v) sum(v < t))[-j]))
    FDP_t <- (numer/denom)*(sum(N_k_test)/(K+1))
    return(FDP_t)
  }
  FDP_hat <- Vectorize(FDP_hat)
  
  T_rej <- rep(0,K_test)
  
  t_v <- c(score_cal_1_sample,unlist(score_test_0))
  t_set <- sort(t_v)
  FDP_t_set <- matrix(0,nrow=K_test,ncol=length(t_set))
  
  for(j in 1:K_test)
  {
    FDP_t_set[j,] <- FDP_hat(t_set,j)
  }
  
  
  l <- 1
  for(alpha in alpha_set)
  {
    alpha <- alpha_set[l]
    alpha_tilde <- 0.9*alpha
    T_rej <- rep(0,K_test)
    for(j in 1:K_test)
    {
      if(min(FDP_t_set[j,]) <= alpha_tilde )
      {
        T_rej[j] <- max(t_set[which(FDP_t_set[j,] <= alpha_tilde )])
      } else{
        T_rej[j] <- -Inf
      } 
    }
    
    
    U <- runif(1)
    e_list <- lapply(1:K_test, function(j) (score_test_0[[j]] < T_rej[j])*(K+1)/(( sum((score_cal_1_sample < T_rej[j]))+1)))
    e_list_U <- lapply(1:K_test, function(j) (score_test_0[[j]] < T_rej[j])*(K+1)/(U*( sum((score_cal_1_sample < T_rej[j]))+1)))
    
    e_threshold <- eBH(unlist(e_list),alpha)
    e_threshold_U <- eBH(unlist(e_list_U),alpha)
    
    p_list <- lapply(score_test_0, function(v) (colSums(sapply(v, function(x) (x>score_cal_1_sample)))+1)/(K+1))
    p_threshold <- BH(unlist(p_list),alpha)
    
    R <- sum(unlist(sapply(1:K_test, function(k) e_list[[k]]>= e_threshold)))
    V <- sum(unlist(sapply(1:K_test, function(k) (e_list[[k]]>= e_threshold)*(Y_test_1[[k]]<=Y_test_0[[k]]))))
    tr_pos <- sum(unlist(sapply(1:K_test, function(k) (e_list[[k]]>= e_threshold)*(Y_test_1[[k]]>Y_test_0[[k]]))))
    
    R_U <- sum(unlist(sapply(1:K_test, function(k) e_list_U[[k]]>= e_threshold_U)))
    V_U <- sum(unlist(sapply(1:K_test, function(k) (e_list_U[[k]]>= e_threshold_U)*(Y_test_1[[k]]<=Y_test_0[[k]]))))
    tr_pos_U <- sum(unlist(sapply(1:K_test, function(k) (e_list_U[[k]]>= e_threshold_U)*(Y_test_1[[k]]>Y_test_0[[k]]))))
    
    pos <- sum(unlist(sapply(1:K_test, function(k) (Y_test_1[[k]]>Y_test_0[[k]]))))
    
    R_p <- sum(unlist(sapply(1:K_test, function(k) p_list[[k]]<= p_threshold)))
    V_p <- sum(unlist(sapply(1:K_test, function(k) (p_list[[k]]<=p_threshold)*(Y_test_1[[k]]<=Y_test_0[[k]]))))
    tr_pos_p <- sum(unlist(sapply(1:K_test, function(k) (p_list[[k]]<= p_threshold)*(Y_test_1[[k]]>Y_test_0[[k]]))))

    FDP_vec[t,l] <- V/max(R,1)
    power_vec[t,l] <- tr_pos/max(pos,1)
    FDP_vec_U[t,l] <- V_U/max(R_U,1)
    power_vec_U[t,l] <- tr_pos_U/max(pos,1)
    FDP_p_vec[t,l] <- V_p/max(R_p,1)
    power_p_vec[t,l] <- tr_pos_p/max(pos,1)
    
    
    l <- l+1
  }
  
  
  setTxtProgressBar(pb, t)
}


#plots
col1 <- rgb(0.098, 0.463, 0.824, alpha=0.2)
col2 <- rgb(0.827, 0.184, 0.184, alpha=0.2)

FDP_mean <- colMeans(FDP_vec)
FDP_se <- apply(FDP_vec,2,sd)/sqrt(tn)
FDP_mean_U <- colMeans(FDP_vec_U)
FDP_se_U <- apply(FDP_vec_U,2,sd)/sqrt(tn)
FDP_p_mean <- colMeans(FDP_p_vec)
FDP_p_se <- apply(FDP_p_vec,2,sd)/sqrt(tn)

power_mean <- colMeans(power_vec)
power_se <- apply(power_vec,2,sd)/sqrt(tn)
power_mean_U <- colMeans(power_vec_U)
power_se_U <- apply(power_vec_U,2,sd)/sqrt(tn)
power_p_mean <- colMeans(power_p_vec)
power_p_se <- apply(power_p_vec,2,sd)/sqrt(tn)

setwd("C:/Wharton/selection for hierarchically structured data")
write.table(cbind(FDP_mean,FDP_se,FDP_mean_U,FDP_se_U,FDP_p_mean,FDP_p_se),"FDP_M_200_lam_10_ite")
write.table(cbind(power_mean,power_se,power_mean_U,power_se_U,power_p_mean,power_p_se),"power_M_200_lam_10_ite")

col1 <- rgb(0.098, 0.463, 0.824, 1)
col2 <- rgb(0.827, 0.184, 0.184, 1)

par(mfrow=c(1,2), cex.lab = 2, cex.main=2, cex.axis=2,mar=c(3,3,3,3), xpd=NA, oma = c(3, 3, 1, 10), mgp = c(3.5, 1.25, 0))
matplot(alpha_set,cbind(FDP_mean,FDP_mean_U,FDP_p_mean),type="l",lwd=2,xlab=expression(alpha),ylab="FDR",ylim=c(0,0.3),col=c(col1,col1,col2),lty=c(1,2,1))
arrows(alpha_set, FDP_mean-FDP_se, alpha_set, FDP_mean+FDP_se, length=0.03, angle=90, code=3,lwd=2,col=col1)
arrows(alpha_set, FDP_mean_U-FDP_se_U, alpha_set, FDP_mean_U+FDP_se_U, length=0.03, angle=90, code=3,lwd=2,col=col1)
arrows(alpha_set, FDP_p_mean-FDP_p_se, alpha_set, FDP_p_mean+FDP_p_se, length=0.03, angle=90, code=3,lwd=2,col=col2)
segments(0.05,0.05,0.25,0.25,lty=3,col="black",lwd=2)


matplot(alpha_set,cbind(power_mean,power_mean_U,power_p_mean),type="l",lwd=2,xlab=expression(alpha),ylab="Power",ylim=c(0,1),col=c(col1,col1,col2),lty=c(1,2,1))
arrows(alpha_set, power_mean-power_se, alpha_set, power_mean+power_se, length=0.03, angle=90, code=3,lwd=2,col=col1)
arrows(alpha_set, power_mean_U-power_se_U, alpha_set, power_mean_U+power_se_U, length=0.03, angle=90, code=3,lwd=2,col=col1)
arrows(alpha_set, power_p_mean-power_p_se, alpha_set, power_p_mean+power_p_se, length=0.03, angle=90, code=3,lwd=2,col=col2)

legend1="e-BH"
legend2="U-eBH"
legend3="p-BH"

legend(0.28,0.78, legend = c(legend1, legend2, legend3), col = c(col1,col1,col2),lty=c(1,2,2),lwd=2,box.lty=0,cex=1.7,bty="n",)


