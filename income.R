source("functions.R")

library(dplyr)

X <- read.table('ca_features.txt')
Y <- read.table('ca_labels.txt')
muhat <- read.table('estimates.txt')

ind_1 <- 9:26
ind_2 <- 36:44
ind_3 <- 45:69

group_1 <- as.factor(apply(X,1,function(x) which(x[ind_1]==1)))
group_2 <- as.factor(apply(X,1,function(x) which(x[ind_2]==1)))
group_3 <- as.factor(apply(X,1,function(x) ifelse(max(x[ind_3])==1,which(x[ind_3]==1),0)))

groups <- as.data.frame(cbind(group_1,group_2,group_3))

Y_list <- split(Y,groups[c("group_1","group_2","group_3")])

length(Y_list)
group_sizes <- unname(sapply(Y_list,function(y) dim(y)[1]))
Y_list_filtered <- Y_list[which(group_sizes >= 10)]
Y_list_filtered <- lapply(Y_list_filtered, function(y) y$V1)

length(Y_list_filtered)

K_cal <- 650
K_test <- 200
ind_shuffle <- sample(1:length(Y_list_filtered))

Y_cal <- Y_list_filtered[ind_shuffle[1:K_cal]]
Y_test <- Y_list_filtered[ind_shuffle[(K_cal+1):(K_cal+K_test)]]


muhat_list <- split(muhat,groups[c("group_1","group_2","group_3")])
muhat_list_filtered <- muhat_list[which(group_sizes >= 10)]
muhat_list_filtered <- lapply(muhat_list_filtered, function(y) y$V1)
muhat_cal <- muhat_list_filtered[ind_shuffle[1:K_cal]]
muhat_test <- muhat_list_filtered[ind_shuffle[(K_cal+1):(K_cal+K_test)]]



N_k <- unname(sapply(Y_cal,length))
N_k_test <- unname(sapply(Y_test,length))
c <- 0.5

score_cal <- lapply(1:K_cal,function(i) Y_cal[[i]]-muhat_cal[[i]])
score_cal_hat <- lapply(1:K_cal,function(i) c-muhat_cal[[i]])
score_test_hat <- lapply(1:K_test,function(i) c-muhat_test[[i]])

set.seed(1234)
alpha_set <- seq(0.05,0.25,by=0.025)

tn <- 500
FDP_vec <- matrix(0,nrow=tn,ncol=length(alpha_set))
power_vec <- matrix(0,nrow=tn,ncol=length(alpha_set))
FDP_vec_U <- matrix(0,nrow=tn,ncol=length(alpha_set))
power_vec_U <- matrix(0,nrow=tn,ncol=length(alpha_set))
FDP_p_vec <- matrix(0,nrow=tn,ncol=length(alpha_set))
power_p_vec <- matrix(0,nrow=tn,ncol=length(alpha_set))


pb <- txtProgressBar(min = 0, max = tn, style = 3)

for(i in 1:tn)
{
  ind <- sapply(1:K_cal, function(k) sample(1:N_k[k],1))
  score_cal_hat_sample <- sapply(1:K_cal, function(k) score_cal_hat[[k]][ind[k]])
  Y_sample <- sapply(1:K_cal, function(k) Y_cal[[k]][ind[k]])
  FDP_hat <- function(t,j)
  {
    numer <- sum((score_cal_hat_sample < t)*(Y_sample <= c))+1
    denom <- max(1,sum(sapply(score_test_hat,function(v) sum(v < t))[-j]))
    FDP_t <- (numer/denom)*(sum(N_k_test)/(K_cal+1))
    return(FDP_t)
  }
  FDP_hat <- Vectorize(FDP_hat)
  T_rej <- rep(0,K_test)
  t_set <- seq(-0.5,0.5,0.025)
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
    e_list <- lapply(1:K_test, function(j) (score_test_hat[[j]] < T_rej[j])*(K_cal+1)/(( sum((score_cal_hat_sample < T_rej[j])*(Y_sample <= c))+1)))
    e_list_U <- lapply(1:K_test, function(j) (score_test_hat[[j]] < T_rej[j])*(K_cal+1)/(U*( sum((score_cal_hat_sample < T_rej[j])*(Y_sample <= c))+1)))
    e_threshold <- eBH(unlist(e_list),alpha)
    e_threshold_U <- eBH(unlist(e_list_U),alpha)
    p_list <- lapply(score_test_hat, function(v) (colSums(sapply(v, function(x) (x>score_cal_hat_sample)*(Y_sample <= c)))+1)/(K_cal+1))
    p_threshold <- BH(unlist(p_list),alpha)
    
    R <- sum(unlist(sapply(1:K_test, function(k) e_list[[k]]>= e_threshold)))
    V <- sum(unlist(sapply(1:K_test, function(k) (e_list[[k]]>= e_threshold)*(Y_test[[k]]<=c))))
    tr_pos <- sum(unlist(sapply(1:K_test, function(k) (e_list[[k]]>= e_threshold)*(Y_test[[k]]>c))))
    
    R_U <- sum(unlist(sapply(1:K_test, function(k) e_list_U[[k]]>= e_threshold_U)))
    V_U <- sum(unlist(sapply(1:K_test, function(k) (e_list_U[[k]]>= e_threshold_U)*(Y_test[[k]]<=c))))
    tr_pos_U <- sum(unlist(sapply(1:K_test, function(k) (e_list_U[[k]]>= e_threshold_U)*(Y_test[[k]]>c))))
    
    pos <- sum(unlist(sapply(1:K_test, function(k) (Y_test[[k]]>c))))
    R_p <- sum(unlist(sapply(1:K_test, function(k) p_list[[k]]<= p_threshold)))
    V_p <- sum(unlist(sapply(1:K_test, function(k) (p_list[[k]]<=p_threshold)*(Y_test[[k]]<=c))))
    tr_pos_p <- sum(unlist(sapply(1:K_test, function(k) (p_list[[k]]<= p_threshold)*(Y_test[[k]]>c))))
    
    FDP_vec[i,l] <- V/max(R,1)
    power_vec[i,l] <- tr_pos/max(pos,1)
    FDP_vec_U[i,l] <- V_U/max(R_U,1)
    power_vec_U[i,l] <- tr_pos_U/max(pos,1)
    FDP_p_vec[i,l] <- V_p/max(R_p,1)
    power_p_vec[i,l] <- tr_pos_p/max(pos,1)
    
    l <- l+1
  }
  
  setTxtProgressBar(pb, i)  
}

FDP_e <- colMeans(FDP_vec)
FDP_e_U <- colMeans(FDP_vec_U)
FDP_p <- colMeans(FDP_p_vec)
FDP_e_se <- apply(FDP_vec,2,sd)/sqrt(tn)
FDP_e_U_se <- apply(FDP_vec_U,2,sd)/sqrt(tn)
FDP_p_se <- apply(FDP_p_vec,2,sd)/sqrt(tn)

FDP <- cbind(FDP_e,FDP_e_se,FDP_e_U,FDP_e_U_se,FDP_p,FDP_p_se )

power_e <- colMeans(power_vec)
power_e_U <- colMeans(power_vec_U)
power_p <- colMeans(power_p_vec)
power_e_se <- apply(power_vec,2,sd)/sqrt(tn)
power_e_U_se <- apply(power_vec_U,2,sd)/sqrt(tn)
power_p_se <- apply(power_p_vec,2,sd)/sqrt(tn)

power <- cbind(power_e,power_e_se,power_e_U,power_e_U_se,power_p,power_p_se )

##----------------------------------------------------------------------

FDP <- read.table("real_FDP.txt")
power <- read.table("real_power.txt")

col1 <- rgb(0.098, 0.463, 0.824, 1)
col2 <- rgb(0.827, 0.184, 0.184, 1)

par(mfrow=c(1,2), cex.lab = 2, cex.main=2, cex.axis=1.5,mar=c(2,2,2,2), xpd=NA, oma = c(4, 2, 1, 1), mgp = c(3.5, 1.25, 0))
matplot(alpha_set,FDP[,c(1,3,5)],type="l",lwd=2,xlab="",ylab="",ylim=c(0,0.3),col=c(col1,col1,col2),lty=c(1,2,1))
title(main="FDP",line=1)
arrows(alpha_set, FDP[,1]-FDP[,2], alpha_set, FDP[,1]+FDP[,2], length=0.02, angle=90, code=3,lwd=2,col=col1)
arrows(alpha_set, FDP[,3]-FDP[,4], alpha_set, FDP[,3]+FDP[,4], length=0.02, angle=90, code=3,lwd=2,col=col1)
arrows(alpha_set, FDP[,5]-FDP[,6], alpha_set, FDP[,5]+FDP[,6], length=0.02, angle=90, code=3,lwd=2,col=col2)
segments(0.05,0.05,0.25,0.25,lty=3,col="black",lwd=2)

matplot(alpha_set,power[,c(1,3,5)],type="l",lwd=2,xlab="",ylab="",ylim=c(0,1),col=c(col1,col1,col2),lty=c(1,2,1))
title(main="Empirical power",line=1)
arrows(alpha_set, power[,1]-power[,2], alpha_set, power[,1]+power[,2], length=0.02, angle=90, code=3,lwd=2,col=col1)
arrows(alpha_set, power[,3]-power[,4], alpha_set, power[,3]+power[,4], length=0.02, angle=90, code=3,lwd=2,col=col1)
arrows(alpha_set, power[,5]-power[,6], alpha_set, power[,5]+power[,6], length=0.02, angle=90, code=3,lwd=2,col=col2)

legend1 <- "e-BH"
legend2 <- "U-eBH"
legend3 <- "p-BH"
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", legend = c(legend1, legend2,legend3), col = c(col1,col1,col2),lty=c(1,2,1),lwd=2,box.lty=0,cex=1.7,bty="n",horiz = T,xpd = TRUE)


