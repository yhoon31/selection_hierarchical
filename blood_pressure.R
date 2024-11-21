source("functions.R")

library(readxl)
set.seed(1234)

#load data
data <- read_excel("real_1.xlsx")

#delete samples with missing observations
na.omit(data)
data <- data[-c(31, 48, 49, 123, 145, 165, 174, 211, 215, 278, 279, 283, 285, 337, 338, 358, 392, 394, 403, 404),]


data$`CV Risk` <- as.factor(as.numeric(factor(data$`CV Risk`)))
data$`Education Level` <- as.factor(as.numeric(factor(data$`Education Level`)))
data$`BMI Classification` <- as.factor(data$`BMI Classification`)
data$Gender <- as.factor(data$Gender)


G <- data$`Site Number`
X <- data[,c(3:5,7:14)]
Y <- data$`SBP reduction`
n <- dim(data)[1]

#generate counterfactual outcomes
X_Y <- cbind(X,Y)
lm1 <- lm(Y~.,X_Y)
lm2 <- lm(Y~.,X_Y)

lm1$coefficients <- c(-55,-1,-2,-7,-2,0.5,0,-2,0,-0.5,0,-0.2,-4,2,10,-3,0,-3,0,4,1,2,0,1,0.5)
lm2$coefficients <- c(-55,1,-8,-3,-2,0.3,0,-2,0,0.5,0,0.3,-6,3,7,-4,0,-1,0,7,-1,1.5,0,0.5,1)

Y1 <- predict(lm1,newdata=X)+rnorm(n,0,1)
Y0 <- predict(lm2,newdata=X)+rnorm(n,0,1)
par(mfrow=c(1,2))
hist(Y1)
hist(Y0)
mean(Y1 > Y0)

groups <- unique(G)
K_total <- length(groups)

X_list <- lapply(1:K_total,function(c) X[G==groups[c],])
Y1_list <- lapply(1:K_total,function(c) Y1[G==groups[c]])
Y0_list <- lapply(1:K_total,function(c) Y0[G==groups[c]])

#data splitting
group_shuffle <- sample(1:length(groups))

K_train <- 5
K_cal <- 20
K_test <- 4

ind_train <- 1:K_train
ind_cal <- (K_train+1):(K_train+K_cal)
ind_test <- (K_train+K_cal+1):(K_train+K_cal+K_test)

X_cal <- X_list[group_shuffle[ind_cal]]
Y_cal <- Y1_list[group_shuffle[ind_cal]]
X_test <- X_list[group_shuffle[ind_test]]
Y_test_1 <- Y1_list[group_shuffle[ind_test]]
Y_test_0 <- Y0_list[group_shuffle[ind_test]]
N_k <- sapply(Y_cal,length)
N_k_test <- sapply(Y_test_0,length)

ind_train_df <- 1:length(which(G <= 5))
X_train_df <- X[ind_train_df,]
Y_train_vec <- Y1[ind_train_df]

#fit random forest regression
library(randomForest)
rf <- randomForest(X_train_df,Y_train_vec)
score_cal_1 <- lapply(1:length(Y_cal),function(j) Y_cal[[j]] - predict(rf,newdata=X_cal[[j]]))
score_test_0 <- lapply(1:length(Y_test_1),function(j) Y_test_0[[j]] - predict(rf,newdata=X_test[[j]]))


set.seed(1234)
alpha_set <- seq(0.05,0.4,by=0.025)

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
  score_cal_1_sample <- sapply(1:K_cal, function(k) score_cal_1[[k]][ind[k]])
  Y_sample <- sapply(1:K_cal, function(k) Y_cal[[k]][ind[k]])
  FDP_hat <- function(t,j)
  {
    numer <- sum((score_cal_1_sample < t))+1
    denom <- max(1,sum(sapply(score_test_0,function(v) sum(v < t))[-j]))
    FDP_t <- (numer/denom)*(sum(N_k_test)/(K_cal+1))
    return(FDP_t)
  }
  FDP_hat <- Vectorize(FDP_hat)
  T_rej <- rep(0,K_test)
  t_set <- seq(-10,10,0.05)
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
    e_list <- lapply(1:K_test, function(j) (score_test_0[[j]] < T_rej[j])*(K_cal+1)/(( sum((score_cal_1_sample < T_rej[j]))+1)))
    e_list_U <- lapply(1:K_test, function(j) (score_test_0[[j]] < T_rej[j])*(K_cal+1)/(U*( sum((score_cal_1_sample < T_rej[j]))+1)))
    e_threshold <- eBH(unlist(e_list),alpha)
    e_threshold_U <- eBH(unlist(e_list_U),alpha)
    p_list <- lapply(score_test_0, function(v) (colSums(sapply(v, function(x) (x>score_cal_1_sample)))+1)/(K_cal+1))
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


#summary of results
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

#plots
col1 <- rgb(0.098, 0.463, 0.824, 1)
col2 <- rgb(0.827, 0.184, 0.184, 1)

par(mfrow=c(1 ,2), cex.lab = 2, cex.main=2, cex.axis=1.5,mar=c(2,2,2,2), xpd=NA, oma = c(4, 2, 1, 1), mgp = c(3.5, 1.25, 0))
matplot(alpha_set,FDP[,c(1,3,5)],type="l",lwd=2,xlab="",ylab="",ylim=c(0,0.45),col=c(col1,col1,col2),lty=c(1,2,1))
title(main="FDP",line=1)
arrows(alpha_set, FDP[,1]-FDP[,2], alpha_set, FDP[,1]+FDP[,2], length=0.02, angle=90, code=3,lwd=2,col=col1)
arrows(alpha_set, FDP[,3]-FDP[,4], alpha_set, FDP[,3]+FDP[,4], length=0.02, angle=90, code=3,lwd=2,col=col1)
arrows(alpha_set, FDP[,5]-FDP[,6], alpha_set, FDP[,5]+FDP[,6], length=0.02, angle=90, code=3,lwd=2,col=col2)
segments(0.05,0.05,0.4,0.4,lty=3,col="black",lwd=2)

matplot(alpha_set,power[,c(1,3,5)],type="l",lwd=2,xlab="",ylab="",ylim=c(0,0.7),col=c(col1,col1,col2),lty=c(1,2,1))
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


