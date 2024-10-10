source("functions.R")

##Testing individual nulls & group-global nulls

library(MASS)
library(randomForest)

p_G <- 10
p_X <- 20

sd_X <- 3
sd_Y <- 1

K_train <- 100
K <- 200
K_test <- 50

set.seed(1234)
A <- matrix(rnorm(p_X*p_G),nrow=p_X)
beta_1 <- runif(p_X)
beta_2 <- runif(p_G)

set.seed(1)
G_train <- matrix(runif(K_train*p_G,-5,5),nrow=K_train)

X_train <- list()
Y_train <- list()
lambda <- 5

X_train_mat <- matrix(nrow=0,ncol=p_X)
G_train_mat <- matrix(nrow=0,ncol=p_G)
N_k_train <-  f_N_pois(K_train,lambda)
for(i in 1:K_train)
{
  G_train_mat <- rbind(G_train_mat,rep(1,N_k_train[i])%*%t(G_train[i,]))
  X_train[[i]] <- mvrnorm(N_k_train[i], mu=A%*%G_train[i,], Sigma=sd_X*diag(p_X))
  X_train_mat <- rbind(X_train_mat,X_train[[i]])
  Y_train[[i]] <- sapply(1:N_k_train[i], function(j) rnorm(1,mean=sum(beta_1*X_train[[i]][j,])+log(abs(sum(beta_2*G_train[i,]))), sd=sd_Y*sqrt(sum(X_train[[i]][j,]^2))/p_X))
}


rf <- randomForest(cbind(G_train_mat,X_train_mat),unlist(Y_train))

plot(X_train_mat[,1],unlist(Y_train))

c <- 20

alpha_set <- seq(0.05,0.25,by=0.025)

tn <- 500

FDP_vec <- matrix(0,nrow=tn,ncol=length(alpha_set))
power_vec <- matrix(0,nrow=tn,ncol=length(alpha_set))
FDP_indiv_vec <- matrix(0,nrow=tn,ncol=length(alpha_set))
power_indiv_vec <- matrix(0,nrow=tn,ncol=length(alpha_set))
FDP_group_vec <- matrix(0,nrow=tn,ncol=length(alpha_set))
power_group_vec <- matrix(0,nrow=tn,ncol=length(alpha_set))

pb <- txtProgressBar(min = 0, max = tn, style = 3)

for(t in 1:tn)
{
  G <- matrix(runif(K*p_G,-5,5),nrow=K)
  G_test <- matrix(runif(K_test*p_G,-5,5),nrow=K_test)
  N_k <- f_N_pois(K,lambda)
  N_k_test <- f_N_pois(K_test,lambda)

  X <- lapply(1:K,function(i) mvrnorm(N_k[i], mu=A%*%G[i,], Sigma=sd_X*diag(p_X)))
  Y <- lapply(1:K, function(i) sapply(1:N_k[i], function(j) rnorm(1,mean=sum(beta_1*X[[i]][j,])+log(abs(sum(beta_2*G[i,]))), sd=sd_Y*sqrt(sum(X[[i]][j,]^2))/p_X)))
  
  X_test <- lapply(1:K_test,function(i) mvrnorm(N_k_test[i], mu=A%*%G_test[i,], Sigma=sd_X*diag(p_X)))
  Y_test <- lapply(1:K_test, function(i) sapply(1:N_k_test[i], function(j) rnorm(1,mean=sum(beta_1*X_test[[i]][j,])+log(abs(sum(beta_2*G_test[i,]))), sd=sd_Y*sqrt(sum(X_test[[i]][j,]^2))/p_X)))

  score_cal <- lapply(1:K,function(i) Y[[i]]-predict(rf,cbind(rep(1,N_k[i])%*%t(G[i,]),X[[i]])))
  score_cal_hat <- lapply(1:K,function(i) c-predict(rf,cbind(rep(1,N_k[i])%*%t(G[i,]),X[[i]])))
  score_test_hat <- lapply(1:K_test,function(i) c-predict(rf,cbind(rep(1,N_k_test[i])%*%t(G_test[i,]),X_test[[i]])))
  
  ind <- sapply(1:K, function(k) sample(1:N_k[k],1))
  
  score_cal_hat_sample <- sapply(1:K, function(k) score_cal_hat[[k]][ind[k]])
  score_cal_sample <- sapply(1:K, function(k) score_cal[[k]][ind[k]])
  Y_sample <- sapply(1:K, function(k) Y[[k]][ind[k]])
  
  FDP_hat <- function(t,j)
  {
    numer <- sum((score_cal_hat_sample < t)*(Y_sample <= c))+1
    denom <- max(1,sum(sapply(score_test_hat,function(v) sum(v < t))[-j]))
    FDP_t <- (numer/denom)*(sum(N_k_test)/(K+1))
    return(FDP_t)
  }
  FDP_hat <- Vectorize(FDP_hat)
  
  T_rej <- rep(0,K_test)
  
  t_v <- c(score_cal_hat_sample,unlist(score_test_hat))
  t_set <- sort(t_v[which((t_v<20)*(t_v>-20)==1)])
  
  FDP_t_set <- t(sapply(1:K_test, function(j) FDP_hat(t_set,j)))
  
  l <- 1
  for(alpha in alpha_set)
  {
    alpha_tilde <- 0.9*alpha
    T_rej <- rep(0,K_test)
    for(j in 1:K_test)
    {
      if(min(FDP_t_set[j,]) <= alpha_tilde)
      {
        T_rej[j] <- max(t_set[which(FDP_t_set[j,] <= alpha_tilde)])
      } else{
        T_rej[j] <- -Inf
      } 
    }
    
    e_list <- lapply(1:K_test, function(j) (score_test_hat[[j]] < T_rej[j])*(K+1)/(( sum((score_cal_hat_sample < T_rej[j])*(Y_sample <= c))+1)))
  
    e_group <- lapply(e_list,mean)
    
    e_threshold <- eBH(c(unlist(e_list),unlist(e_group)),alpha)

    null_global <- sapply(Y_test,function(v) prod(v <= c))
    tr_global <- 1-null_global
  
    R <- sum(unlist(sapply(1:K_test, function(k) e_list[[k]]>= e_threshold)))
    V <- sum(unlist(sapply(1:K_test, function(k) (e_list[[k]]>= e_threshold)*(Y_test[[k]]<=c))))
    R_group <- sum(e_group>= e_threshold)
    V_group <- sum((e_group>= e_threshold)*null_global)
    tr_pos <- sum(unlist(sapply(1:K_test, function(k) (e_list[[k]]>= e_threshold)*(Y_test[[k]]>c))))
    pos <- sum(unlist(sapply(1:K_test, function(k) (Y_test[[k]]>c))))
    tr_pos_group <- sum((e_group>= e_threshold)*tr_global)
    pos_group <- sum(tr_global)
    
    FDP_vec[t,l] <- (V+V_group)/max(R+R_group,1)
    power_vec[t,l] <- (tr_pos+tr_pos_group)/max(pos+pos_group,1)
    FDP_indiv_vec[t,l] <- V/max(R,1)
    power_indiv_vec[t,l] <- tr_pos/max(pos,1)
    FDP_group_vec[t,l] <- V_group/max(R_group,1)
    power_group_vec[t,l] <- tr_pos_group/max(pos_group,1)
    
    
    l <- l+1
  }
  setTxtProgressBar(pb, t)
  
}

FDP_summary <- cbind(colMeans(FDP_vec),apply(FDP_vec,2,sd)/sqrt(tn),colMeans(FDP_group_vec),apply(FDP_group_vec,2,sd)/sqrt(tn),colMeans(FDP_indiv_vec),apply(FDP_indiv_vec,2,sd)/sqrt(tn))
print(round(FDP_summary,4))
power_summary <- cbind(colMeans(power_vec),apply(power_vec,2,sd)/sqrt(tn),colMeans(power_group_vec),apply(power_group_vec,2,sd)/sqrt(tn),colMeans(power_indiv_vec),apply(power_indiv_vec,2,sd)/sqrt(tn))
print(round(power_summary,4))


col1 <- rgb(0.098, 0.463, 0.824, 1)
col2 <- rgb(0.980, 0.686, 0.549, 1)
col3 <- rgb(0.749, 0.686, 0.824, 1)


par(mfrow=c(1,2), cex.lab = 2, cex.main=2, cex.axis=1.5,mar=c(2,3,2,3), xpd=NA, oma = c(6, 4, 1, 1), mgp = c(3.5, 1.25, 0))
matplot(alpha_set,FDP_summary[,c(1,3,5)],type="l",lwd=2,xlab="",ylab="FDR",ylim=c(0,0.25),col=c(col1,col2,col3),lty=c(1,1,1))
arrows(alpha_set, FDP_summary[,1]-FDP_summary[,2], alpha_set, FDP_summary[,1]+FDP_summary[,2], length=0.02, angle=90, code=3,lwd=2,col=col1)
arrows(alpha_set, FDP_summary[,3]-FDP_summary[,4], alpha_set, FDP_summary[,3]+FDP_summary[,4], length=0.02, angle=90, code=3,lwd=2,col=col2)
arrows(alpha_set, FDP_summary[,5]-FDP_summary[,6], alpha_set, FDP_summary[,5]+FDP_summary[,6], length=0.02, angle=90, code=3,lwd=2,col=col3)
segments(0.05,0.05,0.25,0.25,lty=3,col="black",lwd=2)

matplot(alpha_set,power_summary[,c(1,3,5)],type="l",lwd=2,xlab=expression(alpha),ylab="Power",ylim=c(0,1),col=c(col1,col2,col3),lty=c(1,1,1))
arrows(alpha_set, power_summary[,1]-power_summary[,2], alpha_set, power_summary[,1]+power_summary[,2], length=0.02, angle=90, code=3,lwd=2,col=col1)
arrows(alpha_set, power_summary[,3]-power_summary[,4], alpha_set, power_summary[,3]+power_summary[,4], length=0.02, angle=90, code=3,lwd=2,col=col2)
arrows(alpha_set, power_summary[,5]-power_summary[,6], alpha_set, power_summary[,5]+power_summary[,6], length=0.02, angle=90, code=3,lwd=2,col=col3)

legend1="Overall"
legend2="Group"
legend3="Individual"
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", legend = c(legend1, legend2, legend3), col = c(col1,col2,col3),lty=c(1,1,1),lwd=2,box.lty=0,cex=1.7,bty="n",horiz = T,xpd = TRUE)






##Testing individual nulls & group-mean nulls with constant group sizes

p_G <- 10
p_X <- 20

sd_X <- 3
sd_Y <- 1

K_train <- 100
K <- 200
K_test <- 200

set.seed(1234)
A <- matrix(rnorm(p_X*p_G),nrow=p_X)
beta_1 <- runif(p_X)
beta_2 <- runif(p_G)

set.seed(1)
G_train <- matrix(runif(K_train*p_G,-5,5),nrow=K_train)

X_train <- list()
Y_train <- list()
N <- 10

X_train_mat <- matrix(nrow=0,ncol=p_X)
G_train_mat <- matrix(nrow=0,ncol=p_G)
N_k_train <-  f_N_const(N,K_train)
for(i in 1:K_train)
{
  G_train_mat <- rbind(G_train_mat,rep(1,N_k_train[i])%*%t(G_train[i,]))
  X_train[[i]] <- mvrnorm(N_k_train[i], mu=A%*%G_train[i,], Sigma=sd_X*diag(p_X))
  X_train_mat <- rbind(X_train_mat,X_train[[i]])
  Y_train[[i]] <- sapply(1:N_k_train[i], function(j) rnorm(1,mean=sum(beta_1*X_train[[i]][j,])+log(abs(sum(beta_2*G_train[i,]))), sd=sd_Y*sqrt(sum(X_train[[i]][j,]^2))/p_X))
}


rf <- randomForest(cbind(G_train_mat,X_train_mat),unlist(Y_train))

c <- 20

alpha_set <- seq(0.05,0.25,by=0.025)

tn <- 500

FDP_vec <- matrix(0,nrow=tn,ncol=length(alpha_set))
power_vec <- matrix(0,nrow=tn,ncol=length(alpha_set))
FDP_indiv_vec <- matrix(0,nrow=tn,ncol=length(alpha_set))
power_indiv_vec <- matrix(0,nrow=tn,ncol=length(alpha_set))
FDP_group_vec <- matrix(0,nrow=tn,ncol=length(alpha_set))
power_group_vec <- matrix(0,nrow=tn,ncol=length(alpha_set))

pb <- txtProgressBar(min = 0, max = tn, style = 3)

for(t in 1:tn)
{
  G <- matrix(runif(K*p_G,-5,5),nrow=K)
  G_test <- matrix(runif(K_test*p_G,-5,5),nrow=K_test)
  N_k <- f_N_const(N,K)
  N_k_test <- f_N_const(N,K_test)
  
  X <- lapply(1:K,function(i) mvrnorm(N_k[i], mu=A%*%G[i,], Sigma=sd_X*diag(p_X)))
  Y <- lapply(1:K, function(i) sapply(1:N_k[i], function(j) rnorm(1,mean=sum(beta_1*X[[i]][j,])+log(abs(sum(beta_2*G[i,]))), sd=sd_Y*sqrt(sum(X[[i]][j,]^2))/p_X)))
  
  X_test <- lapply(1:K_test,function(i) mvrnorm(N_k_test[i], mu=A%*%G_test[i,], Sigma=sd_X*diag(p_X)))
  Y_test <- lapply(1:K_test, function(i) sapply(1:N_k_test[i], function(j) rnorm(1,mean=sum(beta_1*X_test[[i]][j,])+log(abs(sum(beta_2*G_test[i,]))), sd=sd_Y*sqrt(sum(X_test[[i]][j,]^2))/p_X)))
  
  score_cal <- lapply(1:K,function(i) Y[[i]]-predict(rf,cbind(rep(1,N_k[i])%*%t(G[i,]),X[[i]])))
  score_cal_hat <- lapply(1:K,function(i) c-predict(rf,cbind(rep(1,N_k[i])%*%t(G[i,]),X[[i]])))
  score_test_hat <- lapply(1:K_test,function(i) c-predict(rf,cbind(rep(1,N_k_test[i])%*%t(G_test[i,]),X_test[[i]])))
  
  score_cal_group_hat <- sapply(score_cal_hat,mean)
  score_test_group_hat <- sapply(score_test_hat,mean)
  
  ind <- sapply(1:K, function(k) sample(1:N_k[k],1))
  
  score_cal_hat_sample <- sapply(1:K, function(k) score_cal_hat[[k]][ind[k]])
  score_cal_sample <- sapply(1:K, function(k) score_cal[[k]][ind[k]])
  Y_sample <- sapply(1:K, function(k) Y[[k]][ind[k]])
  Y_cal_mean <- sapply(Y,mean)
  
  FDP_hat <- function(t,j)
  {
    numer <- sum((score_cal_hat_sample < t)*(Y_sample <= c))+1
    denom <- max(1,sum(sapply(score_test_hat,function(v) sum(v < t))[-j]))
    FDP_t <- (numer/denom)*(sum(N_k_test)/(K+1))
    return(FDP_t)
  }
  FDP_hat <- Vectorize(FDP_hat)
  
  t_v <- c(score_cal_hat_sample,unlist(score_test_hat))
  t_set <- sort(t_v[which((t_v<20)*(t_v>-20)==1)])
  
  FDP_t_set <- t(sapply(1:K_test, function(j) FDP_hat(t_set,j)))
  
  FDP_group_hat <- function(t,j)
  {
    numer <- sum((score_cal_group_hat < t)*(Y_cal_mean <= c))+1
    denom <- max(1,sum((score_test_group_hat < t)[-j])+1)
    FDP_t <- (numer/denom)*(K_test/(K+1))
    return(FDP_t)
  }
  FDP_group_hat <- Vectorize(FDP_group_hat)
  
  FDP_t_set_group <- t(sapply(1:K_test, function(j) FDP_group_hat(t_set,j)))
  
  
  l <- 1
  for(alpha in alpha_set)
  {
    alpha_tilde <- 0.9*alpha
    
    T_rej <- rep(0,K_test)
    T_rej_group <- rep(0,K_test)
    for(j in 1:K_test)
    {
      if(min(FDP_t_set[j,]) <= alpha_tilde )
      {
        T_rej[j] <- max(t_set[which(FDP_t_set[j,] <= alpha_tilde)])
      } else{
        T_rej[j] <- -Inf
      } 
      
      if(min(FDP_t_set_group[j,]) <= alpha_tilde )
      {
        T_rej_group[j] <- max(t_set[which(FDP_t_set_group[j,] <= alpha_tilde)])
      } else{
        T_rej_group[j] <- -Inf
      } 
    }
    

    e_list <- lapply(1:K_test, function(j) (score_test_hat[[j]] < T_rej[j])*(K+1)/( sum((score_cal_hat_sample < T_rej[j])*(Y_sample <= c))+1))
    
    e_group <- sapply(1:K_test, function(j) (score_test_group_hat[[j]] < T_rej_group[j])*(K+1)/(sum((score_cal_group_hat < T_rej_group[j])*(Y_cal_mean <= c))+1))

    e_threshold <- eBH(c(unlist(e_list),e_group),alpha)
    
    null_global <- sapply(Y_test,function(v) prod(v <= c))
    tr_global <- 1-null_global
    
    R <- sum(unlist(sapply(1:K_test, function(k) e_list[[k]]>= e_threshold)))
    V <- sum(unlist(sapply(1:K_test, function(k) (e_list[[k]]>= e_threshold)*(Y_test[[k]]<=c))))
    R_group <- sum(e_group>= e_threshold)
    V_group <- sum((e_group>= e_threshold)*null_global)
    tr_pos <- sum(unlist(sapply(1:K_test, function(k) (e_list[[k]]>= e_threshold)*(Y_test[[k]]>c))))
    pos <- sum(unlist(sapply(1:K_test, function(k) (Y_test[[k]]>c))))
    tr_pos_group <- sum((e_group>= e_threshold)*tr_global)
    pos_group <- sum(tr_global)
    
    FDP_vec[t,l] <- (V+V_group)/max(R+R_group,1)
    power_vec[t,l] <- (tr_pos+tr_pos_group)/max(pos+pos_group,1)
    FDP_indiv_vec[t,l] <- V/max(R,1)
    power_indiv_vec[t,l] <- tr_pos/max(pos,1)
    FDP_group_vec[t,l] <- V_group/max(R_group,1)
    power_group_vec[t,l] <- tr_pos_group/max(pos_group,1)
    
    l <- l+1
  }
  setTxtProgressBar(pb, t)
  
}

FDP_summary <- cbind(colMeans(FDP_vec),apply(FDP_vec,2,sd)/sqrt(tn),colMeans(FDP_group_vec),apply(FDP_group_vec,2,sd)/sqrt(tn),colMeans(FDP_indiv_vec),apply(FDP_indiv_vec,2,sd)/sqrt(tn))
print(round(FDP_summary,4))
power_summary <- cbind(colMeans(power_vec),apply(power_vec,2,sd)/sqrt(tn),colMeans(power_group_vec),apply(power_group_vec,2,sd)/sqrt(tn),colMeans(power_indiv_vec),apply(power_indiv_vec,2,sd)/sqrt(tn))
print(round(power_summary,4))

col1 <- rgb(0.098, 0.463, 0.824, 1)
col2 <- rgb(0.980, 0.686, 0.549, 1)
col3 <- rgb(0.749, 0.686, 0.824, 1)


par(mfrow=c(1,2), cex.lab = 2, cex.main=2, cex.axis=1.5,mar=c(2,3,2,3), xpd=NA, oma = c(6, 4, 1, 1), mgp = c(3.5, 1.25, 0))
matplot(alpha_set,FDP_summary[,c(1,3,5)],type="l",lwd=2,xlab="",ylab="FDR",ylim=c(0,0.25),col=c(col1,col2,col3),lty=c(1,1,1))
arrows(alpha_set, FDP_summary[,1]-FDP_summary[,2], alpha_set, FDP_summary[,1]+FDP_summary[,2], length=0.02, angle=90, code=3,lwd=2,col=col1)
arrows(alpha_set, FDP_summary[,3]-FDP_summary[,4], alpha_set, FDP_summary[,3]+FDP_summary[,4], length=0.02, angle=90, code=3,lwd=2,col=col2)
arrows(alpha_set, FDP_summary[,5]-FDP_summary[,6], alpha_set, FDP_summary[,5]+FDP_summary[,6], length=0.02, angle=90, code=3,lwd=2,col=col3)
segments(0.05,0.05,0.25,0.25,lty=3,col="black",lwd=2)

matplot(alpha_set,power_summary[,c(1,3,5)],type="l",lwd=2,xlab=expression(alpha),ylab="Power",ylim=c(0,1),col=c(col1,col2,col3),lty=c(1,1,1))
arrows(alpha_set, power_summary[,1]-power_summary[,2], alpha_set, power_summary[,1]+power_summary[,2], length=0.02, angle=90, code=3,lwd=2,col=col1)
arrows(alpha_set, power_summary[,3]-power_summary[,4], alpha_set, power_summary[,3]+power_summary[,4], length=0.02, angle=90, code=3,lwd=2,col=col2)
arrows(alpha_set, power_summary[,5]-power_summary[,6], alpha_set, power_summary[,5]+power_summary[,6], length=0.02, angle=90, code=3,lwd=2,col=col3)

legend1="Overall"
legend2="Group"
legend3="Individual"
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", legend = c(legend1, legend2, legend3), col = c(col1,col2,col3),lty=c(1,1,1),lwd=2,box.lty=0,cex=1.7,bty="n",horiz = T,xpd = TRUE)




##Testing individual nulls & group-mean nulls with random group sizes

p_G <- 10
p_X <- 20

sd_X <- 3
sd_Y <- 1

K_train <- 100
K <- 200
K_test <- 20

set.seed(1234)
A <- matrix(rnorm(p_X*p_G),nrow=p_X)
beta_1 <- runif(p_X)
beta_2 <- runif(p_G)

G_train <- matrix(runif(K_train*p_G,-5,5),nrow=K_train)

X_train <- list()
Y_train <- list()
lambda <- 5

X_train_mat <- matrix(nrow=0,ncol=p_X)
G_train_mat <- matrix(nrow=0,ncol=p_G)
N_k_train <-  f_N_pois(K_train,lambda)
for(i in 1:K_train)
{
  G_train_mat <- rbind(G_train_mat,rep(1,N_k_train[i])%*%t(G_train[i,]))
  X_train[[i]] <- mvrnorm(N_k_train[i], mu=A%*%G_train[i,], Sigma=sd_X*diag(p_X))
  X_train_mat <- rbind(X_train_mat,X_train[[i]])
  Y_train[[i]] <- sapply(1:N_k_train[i], function(j) rnorm(1,mean=sum(beta_1*X_train[[i]][j,])+log(abs(sum(beta_2*G_train[i,]))), sd=sd_Y*sqrt(sum(X_train[[i]][j,]^2))/p_X))
}


rf <- randomForest(cbind(G_train_mat,X_train_mat),unlist(Y_train))

c <- 20

alpha_set <- seq(0.05,0.25,by=0.025)

tn <- 500

FDP_vec <- matrix(0,nrow=tn,ncol=length(alpha_set))
power_vec <- matrix(0,nrow=tn,ncol=length(alpha_set))
FDP_indiv_vec <- matrix(0,nrow=tn,ncol=length(alpha_set))
power_indiv_vec <- matrix(0,nrow=tn,ncol=length(alpha_set))
FDP_group_vec <- matrix(0,nrow=tn,ncol=length(alpha_set))
power_group_vec <- matrix(0,nrow=tn,ncol=length(alpha_set))

pb <- txtProgressBar(min = 0, max = tn, style = 3)

for(t in 1:tn)
{
  G <- matrix(runif(K*p_G,-5,5),nrow=K)
  G_test <- matrix(runif(K_test*p_G,-5,5),nrow=K_test)
  N_k <- f_N_pois(K,lambda)
  N_k_test <- f_N_pois(K_test,lambda)
  
  X <- lapply(1:K,function(i) mvrnorm(N_k[i], mu=A%*%G[i,], Sigma=sd_X*diag(p_X)))
  Y <- lapply(1:K, function(i) sapply(1:N_k[i], function(j) rnorm(1,mean=sum(beta_1*X[[i]][j,])+log(abs(sum(beta_2*G[i,]))), sd=sd_Y*sqrt(sum(X[[i]][j,]^2))/p_X)))
  
  X_test <- lapply(1:K_test,function(i) mvrnorm(N_k_test[i], mu=A%*%G_test[i,], Sigma=sd_X*diag(p_X)))
  Y_test <- lapply(1:K_test, function(i) sapply(1:N_k_test[i], function(j) rnorm(1,mean=sum(beta_1*X_test[[i]][j,])+log(abs(sum(beta_2*G_test[i,]))), sd=sd_Y*sqrt(sum(X_test[[i]][j,]^2))/p_X)))
  
  score_cal <- lapply(1:K,function(i) Y[[i]]-predict(rf,cbind(rep(1,N_k[i])%*%t(G[i,]),X[[i]])))
  score_cal_hat <- lapply(1:K,function(i) c-predict(rf,cbind(rep(1,N_k[i])%*%t(G[i,]),X[[i]])))
  score_test_hat <- lapply(1:K_test,function(i) c-predict(rf,cbind(rep(1,N_k_test[i])%*%t(G_test[i,]),X_test[[i]])))
  
  score_cal_group_hat <- sapply(score_cal_hat,mean)
  score_test_group_hat <- sapply(score_test_hat,mean)
  
  ind <- sapply(1:K, function(k) sample(1:N_k[k],1))
  
  score_cal_hat_sample <- sapply(1:K, function(k) score_cal_hat[[k]][ind[k]])
  score_cal_sample <- sapply(1:K, function(k) score_cal[[k]][ind[k]])
  Y_sample <- sapply(1:K, function(k) Y[[k]][ind[k]])
  Y_cal_mean <- sapply(Y,mean)
  
  FDP_hat <- function(t,j)
  {
    numer <- sum((score_cal_hat_sample < t)*(Y_sample <= c))+1
    denom <- max(1,sum(sapply(score_test_hat,function(v) sum(v < t))[-j]))
    FDP_t <- (numer/denom)*(sum(N_k_test)/(K+1))
    return(FDP_t)
  }
  FDP_hat <- Vectorize(FDP_hat)
  
  
  t_v <- c(score_cal_hat_sample,unlist(score_test_hat))
  t_set <- sort(t_v[which((t_v<20)*(t_v>-20)==1)])
  
  FDP_t_set <- t(sapply(1:K_test, function(j) FDP_hat(t_set,j)))
  
  score_hat_group_r <- lapply(2:max(N_k), function(r) sapply(score_cal_hat[which(N_k >= r)], function(v) mean(v[1:r])))
  Y_mean_r <- lapply(2:max(N_k), function(r) sapply(Y[which(N_k >= r)], function(v) mean(v[1:r])))
  FDP_group_hat_1 <- function(t)
  {
    return(sapply(2:max(N_k), function(r) (sum((score_hat_group_r[[r-1]] < t)*(Y_mean_r[[r-1]] <= c))+1)/length(which(N_k >= r))))
  }
  FDP_group_hat_1 <- Vectorize(FDP_group_hat_1)
  
  FDP_group_hat <- function(t,j)
  {
    if(N_k_test[j] <= max(N_k))
    {
      FDP_t <- FDP_group_hat_1(t)[N_k_test[j]-1,] * (K_test)/max(1,sum((score_test_group_hat < t)))
    }
    else
    {
      FDP_t <- (K_test)/max(1,sum((score_test_group_hat < t)))
    }
    return(FDP_t)
  }
  FDP_group_hat <- Vectorize(FDP_group_hat)
  
  FDP_t_set_group <- t(sapply(1:K_test, function(j) FDP_group_hat(t_set,j)))
  
  l <- 1
  for(alpha in alpha_set)
  {
    alpha_tilde <- 0.9*alpha
    
    T_rej <- rep(0,K_test)
    T_rej_group <- rep(0,K_test)
    for(j in 1:K_test)
    {
      if(min(FDP_t_set[j,]) <= alpha_tilde )
      {
        T_rej[j] <- max(t_set[which(FDP_t_set[j,] <= alpha_tilde)])
      } else{
        T_rej[j] <- -Inf
      } 
      
      if(min(FDP_t_set_group[j,]) <= alpha_tilde )
      {
        T_rej_group[j] <- max(t_set[which(FDP_t_set_group[j,] <= alpha_tilde)])
      } else{
        T_rej_group[j] <- -Inf
      } 
    }
    
    e_list <- lapply(1:K_test, function(j) (score_test_hat[[j]] < T_rej[j])*(K+1)/(sum((score_cal_hat_sample < T_rej[j])*(Y_sample <= c))+1))
    
    e_group <- sapply(1:K_test, function(j) ifelse(N_k_test[j] <= max(N_k),(score_test_group_hat[[j]] < T_rej_group[j])/(FDP_group_hat_1(T_rej_group[j])[N_k_test[j]-1]),(score_test_group_hat[[j]] < T_rej_group[j])/U))
    
    e_threshold <- eBH(c(unlist(e_list),e_group),alpha)
    
    null_group <- sapply(Y_test,function(v) mean(v) <= c)
    tr_group <- 1-null_group
    
    R <- sum(unlist(sapply(1:K_test, function(k) e_list[[k]]>= e_threshold)))
    V <- sum(unlist(sapply(1:K_test, function(k) (e_list[[k]]>= e_threshold)*(Y_test[[k]]<=c))))
    R_group <- sum(e_group>= e_threshold)
    V_group <- sum((e_group>= e_threshold)*null_group)
    tr_pos <- sum(unlist(sapply(1:K_test, function(k) (e_list[[k]]>= e_threshold)*(Y_test[[k]]>c))))
    pos <- sum(unlist(sapply(1:K_test, function(k) (Y_test[[k]]>c))))
    tr_pos_group <- sum((e_group>= e_threshold)*tr_group)
    pos_group <- sum(tr_group)
    
    FDP_vec[t,l] <- (V+V_group)/max(R+R_group,1)
    power_vec[t,l] <- (tr_pos+tr_pos_group)/max(pos+pos_group,1)
    FDP_indiv_vec[t,l] <- V/max(R,1)
    power_indiv_vec[t,l] <- tr_pos/max(pos,1)
    FDP_group_vec[t,l] <- V_group/max(R_group,1)
    power_group_vec[t,l] <- tr_pos_group/max(pos_group,1)
    
    
    
    l <- l+1
  }
  

  if(anyNA(FDP_vec)==T) {break}
  setTxtProgressBar(pb, t)
  
}


FDP_summary <- cbind(colMeans(FDP_vec),apply(FDP_vec,2,sd)/sqrt(tn),colMeans(FDP_group_vec),apply(FDP_group_vec,2,sd)/sqrt(tn),colMeans(FDP_indiv_vec),apply(FDP_indiv_vec,2,sd)/sqrt(tn))
print(round(FDP_summary,4))
power_summary <- cbind(colMeans(power_vec),apply(power_vec,2,sd)/sqrt(tn),colMeans(power_group_vec),apply(power_group_vec,2,sd)/sqrt(tn),colMeans(power_indiv_vec),apply(power_indiv_vec,2,sd)/sqrt(tn))
print(round(power_summary,4))


col1 <- rgb(0.098, 0.463, 0.824, 1)
col2 <- rgb(0.980, 0.686, 0.549, 1)
col3 <- rgb(0.749, 0.686, 0.824, 1)


par(mfrow=c(1,2), cex.lab = 2, cex.main=2, cex.axis=1.5,mar=c(2,3,2,3), xpd=NA, oma = c(6, 4, 1, 1), mgp = c(3.5, 1.25, 0))
matplot(alpha_set,FDP_summary[,c(1,3,5)],type="l",lwd=2,xlab="",ylab="FDR",ylim=c(0,0.25),col=c(col1,col2,col3),lty=c(1,1,1))
arrows(alpha_set, FDP_summary[,1]-FDP_summary[,2], alpha_set, FDP_summary[,1]+FDP_summary[,2], length=0.02, angle=90, code=3,lwd=2,col=col1)
arrows(alpha_set, FDP_summary[,3]-FDP_summary[,4], alpha_set, FDP_summary[,3]+FDP_summary[,4], length=0.02, angle=90, code=3,lwd=2,col=col2)
arrows(alpha_set, FDP_summary[,5]-FDP_summary[,6], alpha_set, FDP_summary[,5]+FDP_summary[,6], length=0.02, angle=90, code=3,lwd=2,col=col3)
segments(0.05,0.05,0.25,0.25,lty=3,col="black",lwd=2)

matplot(alpha_set,power_summary[,c(1,3,5)],type="l",lwd=2,xlab=expression(alpha),ylab="Power",ylim=c(0,1),col=c(col1,col2,col3),lty=c(1,1,1))
arrows(alpha_set, power_summary[,1]-power_summary[,2], alpha_set, power_summary[,1]+power_summary[,2], length=0.02, angle=90, code=3,lwd=2,col=col1)
arrows(alpha_set, power_summary[,3]-power_summary[,4], alpha_set, power_summary[,3]+power_summary[,4], length=0.02, angle=90, code=3,lwd=2,col=col2)
arrows(alpha_set, power_summary[,5]-power_summary[,6], alpha_set, power_summary[,5]+power_summary[,6], length=0.02, angle=90, code=3,lwd=2,col=col3)

legend1="Overall"
legend2="Group"
legend3="Individual"
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", legend = c(legend1, legend2, legend3), col = c(col1,col2,col3),lty=c(1,1,1),lwd=2,box.lty=0,cex=1.7,bty="n",horiz = T,xpd = TRUE)


