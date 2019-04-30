###STAT225 Project###
library(mvtnorm)
library(MCMCpack)

##Data Manipulation##
#read data
abalone <- read.table("abalone.txt",header = FALSE, sep = ",",dec = ".")
#preprocessdata
abalone <- subset(abalone,V1=="M",select=V2:V9)
colnames(abalone) <- c("X1","X2","X3","X4","X5","X6","X7","Y")
#split into training and test set
test.sample <- sample(1:nrow(abalone),100)
abalone.test <- abalone[test.sample,]
abalone.train <- abalone[-test.sample,]
#setting some universal variables
n <- nrow(abalone.train)
k <- 8
#some manipulation to make testing easier
ntest <- nrow(abalone.test)
y_true <- abalone.test$Y
x.test <- cbind(rep(1,ntest),abalone.test[,1:7])
colnames(x.test)[1] <- "INTERCEPT"


##Bayesian linear regression model##
#fit regular linear regression as start
lm_fit <- lm(Y~.,data=abalone.train)
#obtain parameter estimates: beta(vector),sigma^2,v_beta
beta_hat <- coef(lm_fit)#beta(vector)
v_beta <- vcov(lm_fit)#covariance of beta
s_square <- mean(resid(lm_fit)^2)#error variance
#setting up for sampling
nloop <- 5000
lm_beta <- array(NA,c(nloop,k))
colnames(lm_beta) <- c("INTERCEPT","X1","X2","X3","X4","X5","X6","X7")
lm_sigma_square <- array(NA,nloop)
lm_ypred <- array(NA,c(nloop,ntest))
#MC sampling
for(loop in 1:nloop){
  lm_sigma_square[loop] <- (n-k)*s_square/rchisq(1,n-k)#inv-chi posterior
  lm_beta[loop,] <- rmvnorm(1,beta_hat,v_beta*lm_sigma_square[loop])#mvnormal posterior
  pred_mean <- as.matrix(x.test)%*%lm_beta[loop,]#calculate mean vector for predictive
  for (i in 1:ntest){
    lm_ypred[loop,i] <- rnorm(1,pred_mean[i],lm_sigma_square)
  }
}
#use mean of predictive draw as estimate
lm_ypred_mean <- array(NA,ntest)
for (i in 1:ntest){lm_ypred_mean[i] <- mean(lm_ypred[,i])}
#inference on beta
beta_hat
summary(cbind(lm_beta,lm_sigma_square))
#check model fit
plot(y_true,lm_ypred_mean)
plot(y_true,y_true - lm_ypred_mean)

#Finite Mixture Model
H <- 6
W <- as.matrix(abalone.train)
#function for calculating inverse of sum of square matrix as inv-Wishart parameter
sum_square_matrix <- function(x,mu){
  S <- matrix(rep(0,64),8,8)
  for(i in 1:8){
    for (j in 1:8) {
      S[i,j] <- sum((x[,i]-mu[i])*(x[,j]-mu[j]))
    }
  }
  solve(S)
}
#setting initial values for sampling
pi_0 <- rep(1/H,H)
mu_0 <- array(NA,k)
for(j in 1:k) mu_0[j] <- mean(abalone.train[,j])
sigma_0 <- cov(abalone.train)
invS_0 <- sum_square_matrix(abalone.train,mu_0)
#Gibbs sampling
Nloop <- 200
pi_h <- array(NA,c(Nloop,H))
z <- array(NA,c(Nloop,n,H))
theta <- rep(0,H)
n_h <- rep(0,H)
invS_h <- array(NA,c(k,k))
hit <- rep(0,n)
w_bar <- rep(0,k)
sigma_h <- array(NA,c(Nloop,H,k,k))
mu_h <- array(NA,c(Nloop,H,k))
mix_beta <- array(NA,c(Nloop,H,k))
mix_sigma <- array(NA,c(Nloop,H))
pi_x <- array(NA,c(ntest,H))
mix_ypred <- array(NA,c(Nloop,ntest))

for (i in 1:n) {
  for (h in 1:H) {
    theta[h] <- pi_0[h]*dmvnorm(W[i,],mu_0,sigma_0)
  }
  z[1,i,] <- rmultinom(1,size = 1, prob = theta)
  n_h <- n_h + z[1,i,]
}
pi_h[1,] <- rdirichlet(1,rep(1,H)+n_h)
for (h in 1:H) {
  hit <- z[1,,h]
  for(j in 1:k) w_bar[j] <- mean(abalone.train[hit==1,j])
  invS_h <- sum_square_matrix(abalone.train[hit==1,],w_bar)
  sigma_h[1,h,,] <- riwish(n_h[h],invS_h)
  mu_h[1,h,] <- rmvnorm(1,w_bar,sigma_h[1,h,,]/n_h[h])
  temp_fit <- lm(Y~.,data = abalone.train[hit==1,])
  mix_beta[1,h,] <- coef(temp_fit)
  mix_sigma[1,h] <- mean(resid(temp_fit)^2)
}
for (i in 1:ntest) {
  temp_sum <- 0
  mix_ypred[1,i] <- 0
  for (h in 1:H) {
    pred_mean <- as.matrix(x.test)%*%mix_beta[1,h,]
    mix_ypred[1,i] <- mix_ypred[1,i]+rnorm(1,pred_mean[i],mix_sigma[1,h])*pi_h[1,h]#*dmvnorm(abalone.test[i,-8],mu_h[1,h,-8],sigma_h[1,h,-8,-8])
    #temp_sum <- temp_sum + dmvnorm(abalone.test[i,-8],mu_h[1,h,-8],sigma_h[1,h,-8,-8])*pi_h[1,h]
  }
  #mix_ypred[1,i] <- mix_ypred[1,i] / temp_sum
}

for (m in 2:Nloop) {
  n_h <- rep(0,H)
  for (i in 1:n) {
    for (h in 1:H) {
      theta[h] <- pi_h[m-1,h]*dmvnorm(W[i,],mu_h[m-1,h,],sigma_0)
    }
    z[m,i,] <- rmultinom(1,size = 1, prob = theta)
    n_h <- n_h + z[m,i,]
  }
  pi_h[m,] <- rdirichlet(1,rep(1,H)+n_h)
  for (h in 1:H) {
    hit <- z[m,,h]
    for(j in 1:k) w_bar[j] <- mean(abalone.train[hit==1,j])
    invS_h <- sum_square_matrix(abalone.train[hit==1,],w_bar)
    sigma_h[m,h,,] <- riwish(n_h[h],invS_h)
    mu_h[m,h,] <- rmvnorm(1,w_bar,sigma_h[m,h,,]/n_h[h])
    temp_fit <- lm(Y~.,data = abalone.train[hit==1,])
    mix_beta[m,h,] <- coef(temp_fit)
    mix_sigma[m,h] <- mean(resid(temp_fit)^2)
  }
  for (i in 1:ntest) {
    temp_sum <- 0
    mix_ypred[m,i] <- 0
    for (h in 1:H) {
      pred_mean <- as.matrix(x.test)%*%mix_beta[m,h,]
      mix_ypred[m,i] <- mix_ypred[m,i]+rnorm(1,pred_mean[i],mix_sigma[1,h])*pi_h[m,h]#*dmvnorm(abalone.test[i,-8],mu_h[m,h,-8],sigma_h[m,h,-8,-8])
      #temp_sum <- temp_sum + dmvnorm(abalone.test[i,-8],mu_h[m,h,-8],sigma_h[m,h,-8,-8])*pi_h[m,h]
    }
    #mix_ypred[m,i] <- mix_ypred[m,i] / temp_sum
  }
}

par(mfrow=c(2,3))
plot(y_true,mix_ypred[1,])
plot(y_true,mix_ypred[2,])
plot(y_true,mix_ypred[3,])
plot(y_true,mix_ypred[1,]-y_true)
plot(y_true,mix_ypred[2,]-y_true)
plot(y_true,mix_ypred[3,]-y_true)
