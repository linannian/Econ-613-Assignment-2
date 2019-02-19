getwd()
setwd('Documents/R programming/Econ 613/assignment2/')

## Exercise 1
set.seed(2018)
X1 <- runif(n = 10000,min = 1,max = 3)
X2 <- rgamma(n = 10000,shape = 3,scale = 2)
X3 <- rbinom(n = 10000,size = 1,prob = 0.3)
eps <- rnorm(n = 10000,mean = 2,sd = 1)
k <- as.matrix(cbind(X1,X2,X3,eps))
Y <- matrix(0,ncol = 1, nrow = 10000)
k <- cbind(k,Y)
colnames(k)[5] <- 'Y'
for (i in 1:10000) {
  k[i,5] <- 0.5 + 1.2*k[i,1] - 0.9*k[i,2] + 0.1*k[i,3] + k[i,4]
}
Y_mean <- mean(k[,5])
k <- cbind(k,Y)
colnames(k)[6] <- 'ydum'
for (i in 1:10000){
  if(k[i,5] > Y_mean){
    k[i,6] <- 1
  }
}

## Exercise 2
corr_y_x1 <- cor(k[,1],k[,5]) #0.2158, quite different from 1.2
corr_y_x1

Y <- k[,5]
X <- as.matrix(cbind(1,X1,X2,X3))

getOLSCoef <- function(X,Y){
  coef  <- solve(t(X) %*% X) %*% t(X) %*% Y
  return(coef)
}

getOLSSe <- function(X,Y){
  Coef <- getOLSCoef(X,Y)
  resid <- Y -  X %*% Coef
  sigma2 <- (t(resid) %*% resid) / (nrow(X) - ncol(X))
  vcov_beta <- c(sigma2) *(solve(t(X) %*% X))
  StdErr_OLS <- as.matrix(sqrt(diag(vcov_beta)))
  return(StdErr_OLS)
}

Coef_OLS <- getOLSCoef(X,Y)
StdErr_OLS <- getOLSSe(X,Y)

regressor <- cbind(X,Y)
stdErr_boot <- matrix(0,ncol = 4,nrow = 49)
colnames(stdErr_boot) <- c('const_SE','X1_SE','X2_SE','X3_SE')
holder <- matrix(0,ncol = 5, nrow = 10000)

for (i in 1:49) {
  resample <- sample(1:10000,10000,replace = T)
  for (j in 1:10000){
    holder[j,] <- regressor[resample[j],]
  }
  se <- getOLSSe(holder[,1:4],holder[,5])
  stdErr_boot[i,] <- t(se)
}
stdErr_boot_49 <- as.matrix(colMeans(stdErr_boot))
stdErr_boot_49

stdErr_boot_large <- matrix(0,ncol = 4,nrow = 499)
colnames(stdErr_boot_large) <- c('const_SE','X1_SE','X2_SE','X3_SE')
for (i in 1:499) {
  resample <- sample(1:10000,10000,replace = T)
  for (j in 1:10000){
    holder[j,] <- regressor[resample[j],]
  }
  se <- getOLSSe(holder[,1:4],holder[,5])
  stdErr_boot_large[i,] <- t(se)
}
stdErr_boot_499 <- as.matrix(colMeans(stdErr_boot_large))
stdErr_boot_499

SE <- cbind(StdErr_OLS,stdErr_boot_49,stdErr_boot_499)
colnames(SE) <- c('OLS','49_boot','499_boot')
SE

## Exercise 3
# return log-likelihood function
loglik_probit <- function(X,Y,beta){
  norm_cdf <- pnorm(X %*% beta)
  n <- length(Y)
  k <- length(beta)
  f <- sum(Y*log(norm_cdf)) + sum((1-Y)*log(1-norm_cdf))
  f <- -f
  return(f)
}

## gradient function for probit_loglik function
loglik_probit_deriv <- function(X,Y,beta){
  norm_cdf <- pnorm(X %*% beta)
  norm_pdf <- dnorm(X %*% beta)
  n <- length(Y)
  k <- length(beta)
  gradient <- t(matrix(rep(norm_pdf/norm_cdf,k),nrow=n)*X) %*% Y - 
    t(matrix(rep(norm_pdf/(1-norm_cdf),k),nrow=n)*X) %*% (1-Y)
  gradient <- -1*gradient
  return(gradient)
}

y_dum <- k[,6]
# Implement gradient descent
b <- Coef_OLS ## set a start value for anything
diff <- 1
alpha <- 0.000001

while (diff > 1e-9) {
  ll1 <- loglik_probit(X,y_dum,b)
  g <- loglik_probit_deriv(X,y_dum,b)
  bn <- b - alpha*g
  ll2 <- loglik_probit(X,y_dum,bn)
  diff <- abs(ll1 - ll2)/abs(ll2)
  b <- bn
  print(t(b))
}
b # 2.913662 1.312044 -0.9267062 0.2236294
## b is the answer. From the regression in question 4, 
## this is fairly close to the right answer. if we are talking about
## the coefficient 1.2,-0.9,0.1, b is also very close

## Exercise 4
# logit
logit_fit <- glm(formula = k[,6]~ X[,2] + X[,3] + X[,4],family = binomial(link = 'logit'))
summary(logit_fit) #all coefficient are significant at 99% confidence level
beta_logit <- summary(logit_fit)$coefficient[,1]
beta_logit ## 5.2166223 2.3847400 -1.6706475 0.3916688 

loglik_logit <- function(X,Y,beta){
  loglik <- -sum(-Y*log(1 + exp(-(X%*%beta))) - (1-Y)*log(1 + exp(X%*%beta)))
  return(loglik)
}

optim_logit <- optim(par = Coef_OLS,loglik_logit,X = X, Y=y_dum)
optim_logit$par ## 5.216315 2.382636 -1.669837 0.3902768 quite close

# probit
probit_fit <- glm(formula = k[,6]~ X[,2] + X[,3] + X[,4],family = binomial(link = 'probit'))
summary(probit_fit) #all coefficient are significant at 99% confidence level
beta_probit <- summary(probit_fit)$coefficient[,1]
beta_probit ## 2.9284027 1.3096814 -0.9283100 0.2228331 

optim_probit <- optim(par = Coef_OLS,loglik_probit,X = X, Y=y_dum)
optim_probit$par ## 2.928325 1.309655 -0.9282772 0.2228664, very close to the last one

# linear model
linear_fit <- lm(formula = k[,6]~ + X[,2] + X[,3] + X[,4])
summary(linear_fit) #all coefficient are significant at 99% confidence level
coef_linear <- getOLSCoef(X,y_dum)
coef_linear # 0.8932198 0.1417872 -0.102673 0.02826532

## Exercise 5
# Marginal effects of probit & logit model
margin_mean_probit <- function(X, b){
  m_p <- dnorm(X%*%b) %*% b
  mean_probit <- colMeans(m_p)
  return(mean_probit)
}

margin_mean_logit <- function(X,b){
  m_l <- (exp(X%*%b)/(1+exp(X%*%b))^2) %*% b
  mean_logit <- colMeans(m_l)
  return(mean_logit)
}

marginal_probit_in_mean <- margin_mean_probit(X,beta_probit)
marginal_probit_in_mean

marginal_logit_in_mean <- margin_mean_logit(X,beta_logit)
marginal_logit_in_mean

# delta method
cov_matrix_probit <- vcov(probit_fit)
cov_matrix_logit <- vcov(logit_fit)

m_beta_probit <- as.matrix(beta_probit,nrow=4)
counter_matrix_probit <- matrix(0,nrow = 4,ncol = 4)

for (i in 1:10000){
  phi <- as.numeric(dnorm(X[i,] %*% m_beta_probit))
  iden <- diag(4)
  xb <- as.numeric(X[i,] %*% m_beta_probit) 
  bx <- m_beta_probit %*% X[i,]
  counter_matrix_probit <- counter_matrix_probit + (phi * (iden - xb*bx))
}
Jacobian_probit <- counter_matrix_probit/10000
delta_vcov_probit <- Jacobian_probit %*% cov_matrix_probit %*% t(Jacobian_probit)
std_delta_probit <- diag(delta_vcov_probit)^0.5
std_delta_probit ## 0.0090222694 0.0042899569 0.0003706689 0.0054517609 

m_beta_logit <- as.matrix(beta_logit,nrow=4)
counter_matrix_logit <- matrix(0,nrow = 4,ncol = 4)
for (i in 1:10000){
  f <- exp(as.numeric(X[i,] %*% m_beta_logit))/(1+exp(as.numeric(X[i,] %*% m_beta_logit)))
  bx <- m_beta_logit %*% X[i,]
  iden <- diag(4)
  counter_matrix_logit <- counter_matrix_logit + f*(1-f)*(iden + (1-2*f) * bx)
}
Jacobian_logit <- counter_matrix_logit/10000
delta_vcov_logit <- Jacobian_logit %*% cov_matrix_logit %*% t(Jacobian_logit)
std_delta_logit <- diag(delta_vcov_logit)^0.5
std_delta_logit ## 0.008976668 0.004277820 0.000364945 0.005442773 

## 2nd way to estimate Jacobian
Jacobian_probit2 <- matrix(0,ncol = 4,nrow = 4)
for (i in 1:4){
  beta_plus <- beta_probit
  beta_minus <- beta_probit
  beta_plus[i] <- beta_plus[i]+0.00001
  beta_minus[i] <- beta_minus[i]-0.00001
  mean1 <- margin_mean_probit(X,beta_plus)
  mean2 <- margin_mean_probit(X,beta_minus)
  Jacobian_probit2[,i] <- (mean1-mean2)/(2*0.00001)
}
delta_vcov_probit2 <- Jacobian_probit2 %*% cov_matrix_probit %*% t(Jacobian_probit2)
std_delta_probit2 <- diag(delta_vcov_probit2)^0.5
std_delta_probit2 # 0.0090222694 0.0042899569 0.0003706689 0.0054517609 Exactly same to the last one

Jacobian_logit2 <- matrix(0,ncol = 4, nrow = 4)
for (i in 1:4){
  beta_plus_l <- beta_logit
  beta_minus_l <- beta_logit
  beta_plus_l[i] <- beta_plus_l[i]+0.00001
  beta_minus_l[i] <- beta_minus_l[i]-0.00001
  mean1_l <- margin_mean_logit(X,beta_plus_l)
  mean2_l <- margin_mean_logit(X,beta_minus_l)
  Jacobian_logit2[,i] <- (mean1_l-mean2_l)/(2*0.00001)
}
delta_vcov_logit2 <- Jacobian_logit2 %*% cov_matrix_logit %*% t(Jacobian_logit2)
std_delta_logit2 <- diag(delta_vcov_logit2)^0.5
std_delta_logit2 # 0.008976668 0.004277820 0.000364945 0.005442773 

## bootstrap method
regressor_p_l <- cbind(X,y_dum)
holder_boot_probit <- matrix(0,nrow = 499, ncol = 4)
for (i in 1:499){
  resample <- sample(1:10000,10000,replace = T)
  for (j in 1:10000){
    holder[j,] <- regressor_p_l[resample[j],]
  }
  fit <- glm(formula = holder[,5] ~ holder[,2]  + holder[,3] + holder[,4],family = binomial(link = 'probit'))
  Beta <- summary(fit)$coefficient[,1]
  marg <- dnorm(X%*%Beta) %*% Beta
  holder_boot_probit[i,] <- colMeans(marg)
}
std_boot_probit <- diag(var(holder_boot_probit))^0.5
std_boot_probit ## 0.0091417027 0.0043695428 0.0003929125 0.0057003705

holder_boot_logit <- matrix(0,nrow = 499, ncol = 4)
for (i in 1:499){
  resample <- sample(1:10000,10000,replace = T)
  for (j in 1:10000){
    holder[j,] <- regressor_p_l[resample[j],]
  }
  fit_logi <- glm(formula = holder[,5] ~ holder[,2]  + holder[,3] + holder[,4],family = binomial(link = 'logit'))
  Beta_logi <- summary(fit_logi)$coefficient[,1]
  marg_logi <- (exp(X%*%Beta_logi)/(1+exp(X%*%Beta_logi))^2) %*% Beta_logi
  holder_boot_logit[i,] <- colMeans(marg_logi)
}
std_boot_logit <- diag(var(holder_boot_logit))^0.5
std_boot_logit ## 0.0094798995 0.0045467522 0.0003553502 0.0053360342
