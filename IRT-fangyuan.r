# IRT project

# Instruction:
  # change n and k in *generate parameters* part and fit a new model
  # change the item number and student number in *RQR checking* part to conduct graphical checking for the fitted model

## prepare environment
library(tidyverse)
library(rstan)
library(KernSmoothIRT)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# date simulation

## generate parameters
set.seed(1)
n <- 526  # number of test takers  526/1052
k <- 30     # number of items
alpha <- rlnorm(k,0,0.3)
beta <- rlnorm(k,0,0.3)
theta <- rnorm(n,0,1)
gam <- rbeta(k,5,23)
if(n==526){
  testn <- 126
}else{
  testn <- 352
}
testset <- sample(n,testn)    # index of test set


## calculate probability using simulated parameters and generate response
  ## function IRT: 
          ## input: alpha, beta, theta, gamma
          ## output: probability of correct response
          ## - For Rasch Model, set alpha = 1 and gamma = 0
          ## - For 2PL Model, set gamma = 0

# model definition
IRT <- function(alpha,beta,theta,gam){
  gam + (1-gam)*exp(alpha*(theta-beta))/(1+exp(alpha*(theta-beta)))
}
# prepare empty matrix
Rasch_matrix <- matrix(nrow=n,ncol=k) # store probability
Rasch_response <- matrix(nrow=n,ncol=k) # store response
PL2_matrix <- matrix(nrow=n,ncol=k)
PL2_response <- matrix(nrow=n,ncol=k)
PL3_matrix <- matrix(nrow=n,ncol=k)
PL3_response <- matrix(nrow=n,ncol=k)
# fill matrix
for(i in 1:n){
  for(j in 1:k){
    # Rasch
    p <- IRT(1,beta[j],theta[i],0)
    Rasch_matrix[i,j] <- p
    if(runif(1) <= p){
      Rasch_response[i,j] <- 1
    }else{Rasch_response[i,j] <- 0}
    # 2PL
    p <- IRT(alpha[j],beta[j],theta[i],0)
    PL2_matrix[i,j] <- p
    if(runif(1) <= p){
      PL2_response[i,j] <- 1
    }else{PL2_response[i,j] <- 0}
    # 3PL
    p <- IRT(alpha[j],beta[j],theta[i],0.1)  # assume gamma = 0.1 for all items
    PL3_matrix[i,j] <- p
    if(runif(1) <= p){
      PL3_response[i,j] <- 1
    }else{PL3_response[i,j] <- 0}
  }
}

# real data
#library(TAM)
#data(data.sim.rasch)
#dat <- data.sim.rasch
#n <- 2000    # number of test takers
#k <- 40
#testset <- sample(n,1500)


## estimate parameters using Rstan
  ## Each of the model fitting functions will return a list of two stanfit objects, corresponding to in-sample case and out-of-sample case.  
  
  #### The running time of this part is about 45-50 minutes.

Rasch_est <- function(n,k,testset,response){
  model_Rasch <- "
    data{
      int<lower=0> n; 
      int<lower=0> k; 
      int<lower=0,upper=1> y[n,k];
    }
    parameters {
      vector[k] beta;
      vector[n] theta;
    }
    model {
      theta ~ normal(0,1);
      beta ~ normal(0,1);
      
      for(i in 1:n){
        for(j in 1:k){
          y[i,j] ~ bernoulli_logit(theta[i]-beta[j]);
        }
      }
    }
  "
  # in-sample
  stan.data <- list(n=n,k=k,y=response)
  fit <- stan(model_code = model_Rasch,
              data = stan.data,
              iter = 2000,
              chains = 4,
              cores = 2)
  
  # out-of-sample
  stan.data2 <- list(n=n-length(testset),k=k,y=response[-testset,])
  fit2 <- stan(model_code = model_Rasch,
               data = stan.data2,
               iter = 2000,
               chains = 4,
               cores = 2)
  
  return(c(fit,fit2))
}
PL2_est <- function(n,k,testset,response){
  model_2PL <- "
    data{
      int<lower=0> n; 
      int<lower=0> k; 
      int<lower=0,upper=1> y[n,k];
    }
    parameters {
      vector<lower=0>[k] alpha;
      vector[k] beta;
      vector[n] theta;
    }
    model {
      theta ~ normal(0,1);
      alpha ~ lognormal(0,0.3);
      beta ~ normal(0,1);
      
      for(i in 1:n){
        for(j in 1:k){
          y[i,j] ~ bernoulli_logit(alpha[j]*(theta[i]-beta[j]));
        }
      }
    }
  "
  # in-sample
  stan.data <- list(n=n,k=k,y=response)
  fit <- stan(model_code = model_2PL,
              data = stan.data,
              iter = 2000,
              chains = 4,
              cores = 2)
  
  # out-of-sample
  stan.data2 <- list(n=n-length(testset),k=k,y=response[-testset,])
  fit2 <- stan(model_code = model_2PL,
               data = stan.data2,
               iter = 2000,
               chains = 4,
               cores = 2)
  
  return(c(fit,fit2))
}
PL3_est <- function(n,k,gam,testset,response){
  model_3PL <- "
    data{
      int<lower=0> n; 
      int<lower=0> k; 
      int<lower=0,upper=1> y[n,k];
    }
    parameters {
      vector<lower=0>[k] alpha;
      vector[k] beta;
      vector[n] theta;
      vector<lower=0,upper=1>[k] gamma; //item pseudo-guessing
    }
    model {
      theta ~ normal(0,1);
      alpha ~ lognormal(0,0.3);
      beta ~ normal(0,1);
      gamma ~ beta(5,23);
      
      for(i in 1:n){
        for(j in 1:k){
          real p;
          p = inv_logit(alpha[j]*(theta[i]-beta[j]));
          y[i,j] ~ bernoulli(gamma[j]+(1-gamma[j])*p);
        }
      }
    }
  "
  # in-sample
  stan.data <- list(n=n,k=k,gam=gam,y=response)
  fit <- stan(model_code = model_3PL,
              data = stan.data,
              iter = 2000,
              chains = 4,
              cores = 2)
  
  # out-of-sample
  stan.data2 <- list(n=n-length(testset),k=k,y=response[-testset,])
  fit2 <- stan(model_code = model_3PL,
               data = stan.data2,
               iter = 2000,
               chains = 4,
               cores = 2)
  
  return(c(fit,fit2))
}

# prepare data for RQR checking
result_Rasch <- Rasch_est(n,k,testset,PL2_response)
result_2PL <- PL2_est(n,k,testset,PL2_response)
result_3PL <- PL3_est(n,k,gam,testset,PL2_response)
#result_Rasch <- Rasch_est(n,k,testset,dat)
#result_2PL <- PL2_est(n,k,testset,dat)
#result_3PL <- PL3_est(n,k,gam,testset,dat)



## RQR checking
  ## In this part, I set 1 - Rasch, 2 - 2PL, 3 - 3PL
  ## For Item check, input item number, stan outputs of the three models generated in the functions above.
  ## For person check, input student number, stan outputs of the three models generated in the functions above.

item_checking <- function(item_number, testset, Rasch, PL2, PL3, response1,response2,response3){
  item <- item_number
  testn <- length(testset)
  is1 <- rstan::extract(Rasch[[1]],permuted = TRUE)
  os1 <- rstan::extract(Rasch[[2]],permuted = TRUE)
  is2 <- rstan::extract(PL2[[1]],permuted = TRUE)
  os2 <- rstan::extract(PL2[[2]],permuted = TRUE)
  is3 <- rstan::extract(PL3[[1]],permuted = TRUE)
  os3 <- rstan::extract(PL3[[2]],permuted = TRUE)
  
  # RQR calculation function
  # input: probability, response / output: RQR
  RQR <- function(p,o){
    if(o==0){
      rqr <- map(p,function(x){qnorm(x)})
    }else{
      rqr <- map(p,function(x){qnorm(x)})
    }
    unlist(rqr)
  }
  
  # in-sample
  # prepare sampling parameters
  beta1 <- is1$beta[2001:dim(is1$beta)[1],]
  theta1 <- is1$theta[2001:dim(is1$theta)[1],]
  alpha2 <- is2$alpha[2001:dim(is2$alpha)[1],]
  beta2 <- is2$beta[2001:dim(is2$beta)[1],]
  theta2 <- is2$theta[2001:dim(is2$theta)[1],]
  alpha3 <- is3$alpha[2001:dim(is3$alpha)[1],]
  beta3 <- is3$beta[2001:dim(is3$beta)[1],]
  theta3 <- is3$theta[2001:dim(is3$theta)[1],]
  gamma3 <- is3$gamma[2001:dim(is3$gamma)[1],]
  # calculate RQR
  rqr1_in_sample <- c()
  rqr2_in_sample <- c()
  rqr3_in_sample <- c()
  for(i in 1:testn){
    # Rasch
    p <- IRT(rep(1,2000),beta1[,item],theta1[,testset[i]],rep(0,2000))
    rqr1_in_sample <- c(rqr1_in_sample,RQR(p,response1[testset[i],item]))
    # 2PL
    p <- IRT(alpha2[,item],beta2[,item],theta2[,testset[i]],rep(0,2000))
    rqr2_in_sample <- c(rqr2_in_sample,RQR(p,response2[testset[i],item]))
    # 3PL
    p <- IRT(alpha3[,item],beta3[,item],theta3[,testset[i]],gamma3[,item])
    rqr3_in_sample <- c(rqr3_in_sample,RQR(p,response3[testset[i],item]))
  }
  
  # out-of-sample
  # prepare sampling parameters
  beta1 <- os1$beta[2001:dim(os1$beta)[1],]
  alpha2 <- os2$alpha[2001:dim(os2$alpha)[1],]
  beta2 <- os2$beta[2001:dim(os2$beta)[1],]
  alpha3 <- os3$alpha[2001:dim(os3$alpha)[1],]
  beta3 <- os3$beta[2001:dim(os3$beta)[1],]
  gamma3 <- os3$gamma[2001:dim(os3$gamma)[1],]
  # calculate RQR
  # rank in all students
  z_score1 <- scale(rank(apply(response1,1,sum)))
  z_score2 <- scale(rank(apply(response2,1,sum)))
  z_score3 <- scale(rank(apply(response3,1,sum)))
  rqr1_out_of_sample <- c()
  rqr2_out_of_sample <- c()
  rqr3_out_of_sample <- c()
  for(i in 1:testn){
    # Rasch
    p <- IRT(rep(1,2000),beta1[,item],rep(z_score1[i],2000),rep(0,2000))
    rqr1_out_of_sample <- c(rqr1_out_of_sample, RQR(p,response1[testset[i],item]))
    # 2PL
    p <- IRT(alpha2[,item],beta2[,item],rep(z_score2[i],2000),rep(0,2000))
    rqr2_out_of_sample <- c(rqr2_out_of_sample, RQR(p,response2[testset[i],item]))
    # 3PL
    p <- IRT(alpha3[,item],beta3[,item],rep(z_score3[i],2000),gamma3[,item])
    rqr3_out_of_sample <- c(rqr3_out_of_sample, RQR(p,response3[testset[i],item]))
  }
  
  # plot 
  par(mfrow=c(1,2))
  qqnorm(rqr1_in_sample,main = paste("item",as.character(item),"in-sample"),ylab="Rasch",pch=20,cex=0.5)
  qqnorm(rqr1_out_of_sample,main = paste("item",as.character(item),"out-of-sample"),ylab="Rasch",pch=20,cex=0.5)
  par(mfrow=c(1,2))
  qqnorm(rqr2_in_sample,main = paste("item",as.character(item),"in-sample"),ylab="2PL",pch=20,cex=0.5)
  qqnorm(rqr2_out_of_sample,main = paste("item",as.character(item),"out-of-sample"),ylab="2PL",pch=20,cex=0.5)
  par(mfrow=c(1,2))
  qqnorm(rqr3_in_sample,main = paste("item",as.character(item),"in-sample"),ylab="3PL",pch=20,cex=0.5)
  qqnorm(rqr3_out_of_sample,main = paste("item",as.character(item),"out-of-sample"),ylab="3PL",pch=20,cex=0.5)
}
student_checking <- function(student_number, k, Rasch, PL2, PL3, response1,response2,response3){
  student <- student_number
  is1 <- rstan::extract(Rasch[[1]],permuted = TRUE)
  os1 <- rstan::extract(Rasch[[2]],permuted = TRUE)
  is2 <- rstan::extract(PL2[[1]],permuted = TRUE)
  os2 <- rstan::extract(PL2[[2]],permuted = TRUE)
  is3 <- rstan::extract(PL3[[1]],permuted = TRUE)
  os3 <- rstan::extract(PL3[[2]],permuted = TRUE)
  
  # RQR calculation function
  # input: probability, response / output: RQR
  RQR <- function(p,o){
    if(o==0){
      rqr <- map(p,function(x){qnorm(x)})
    }else{
      rqr <- map(p,function(x){qnorm(x)})
    }
    unlist(rqr)
  }
  
  # in-sample
  # prepare sampling parameters
  beta1 <- is1$beta[2001:dim(is1$beta)[1],]
  theta1 <- is1$theta[2001:dim(is1$theta)[1],]
  alpha2 <- is2$alpha[2001:dim(is2$alpha)[1],]
  beta2 <- is2$beta[2001:dim(is2$beta)[1],]
  theta2 <- is2$theta[2001:dim(is2$theta)[1],]
  alpha3 <- is3$alpha[2001:dim(is3$alpha)[1],]
  beta3 <- is3$beta[2001:dim(is3$beta)[1],]
  theta3 <- is3$theta[2001:dim(is3$theta)[1],]
  gamma3 <- is3$gamma[2001:dim(is3$gamma)[1],]
  # calculate RQR
  rqr1_in_sample <- c()
  rqr2_in_sample <- c()
  rqr3_in_sample <- c()
  for(i in 1:k){
    # Rasch
    p <- IRT(rep(1,2000),beta1[,i],theta1[,student],rep(0,2000))
    rqr1_in_sample <- c(rqr1_in_sample,RQR(p,response1[student,i]))
    # 2PL
    p <- IRT(alpha2[,i],beta2[,i],theta2[,student],rep(0,2000))
    rqr2_in_sample <- c(rqr2_in_sample,RQR(p,response2[student,i]))
    # 3PL
    p <- IRT(alpha3[,i],beta3[,i],theta3[,student],gamma3[,i])
    rqr3_in_sample <- c(rqr3_in_sample,RQR(p,response3[student,i]))
  }
  
  # out-of-sample
  # prepare sampling parameters
  beta1 <- os1$beta[2001:dim(os1$beta)[1],]
  alpha2 <- os2$alpha[2001:dim(os2$alpha)[1],]
  beta2 <- os2$beta[2001:dim(os2$beta)[1],]
  alpha3 <- os3$alpha[2001:dim(os3$alpha)[1],]
  beta3 <- os3$beta[2001:dim(os3$beta)[1],]
  gamma3 <- os3$gamma[2001:dim(os3$gamma)[1],]
  # calculate RQR
  # rank in all students
  z_score1 <- scale(rank(apply(response1,1,sum)))
  z_score2 <- scale(rank(apply(response2,1,sum)))
  z_score3 <- scale(rank(apply(response3,1,sum)))
  rqr1_out_of_sample <- c()
  rqr2_out_of_sample <- c()
  rqr3_out_of_sample <- c()
  for(i in 1:k){
    # Rasch
    p <- IRT(rep(1,2000),beta1[,i],rep(z_score1[student],2000),rep(0,2000))
    rqr1_out_of_sample <- c(rqr1_out_of_sample, RQR(p,response1[student,i]))
    # 2PL
    p <- IRT(alpha2[,i],beta2[,i],rep(z_score2[student],2000),rep(0,2000))
    rqr2_out_of_sample <- c(rqr2_out_of_sample, RQR(p,response2[student,i]))
    # 3PL
    p <- IRT(alpha3[,i],beta3[,i],rep(z_score3[student],2000),gamma3[,i])
    rqr3_out_of_sample <- c(rqr3_out_of_sample, RQR(p,response3[student,i]))
  }
  
  # plot 
  par(mfrow=c(1,2))
  qqnorm(rqr1_in_sample,main = paste("student",as.character(student),"in-sample"),ylab="Rasch",pch=20,cex=0.5)
  qqnorm(rqr1_out_of_sample,main = paste("student",as.character(student),"out-of-sample"),ylab="Rasch",pch=20,cex=0.5)
  par(mfrow=c(1,2))
  qqnorm(rqr2_in_sample,main = paste("student",as.character(student),"in-sample"),ylab="2PL",pch=20,cex=0.5)
  qqnorm(rqr2_out_of_sample,main = paste("student",as.character(student),"out-of-sample"),ylab="2PL",pch=20,cex=0.5)
  par(mfrow=c(1,2))
  qqnorm(rqr3_in_sample,main = paste("student",as.character(student),"in-sample"),ylab="3PL",pch=20,cex=0.5)
  qqnorm(rqr3_out_of_sample,main = paste("student",as.character(student),"out-of-sample"),ylab="3PL",pch=20,cex=0.5)
}

  ## input different item number and student number to do graphical checking
item_checking(1,testset,result_Rasch,result_2PL,result_3PL,PL2_response,PL2_response,PL2_response)
student_checking(testset[1],k,result_Rasch,result_2PL,result_3PL,PL2_response,PL2_response,PL2_response)
#item_checking(1,testset,result_Rasch,result_2PL,result_3PL,dat,dat,dat)
#student_checking(testset[1],k,result_Rasch,result_2PL,result_3PL,dat,dat,dat)




## Kernel Smoothing Model checking
   ## not finished yet
ksModel <- ksIRT(PL2_response,key = 1,format = 1,kernel = "gaussian",bandwidth = "Silverman")
#plot(ksModel,items=c(1))
item_ooc_y <- ksModel[["OCC"]][seq(1,k*2,2),4:54]
plot(ksModel[["expectedscores"]],item_ooc_y[1,],type="l",lwd=2)

item_plot <- function(item_number, testset, Rasch, PL2, PL3, response1,response2,response3){
  item <- item_number
  testn <- length(testset)
  k = dim(response1)[2]
  is1 <- rstan::extract(Rasch[[1]],permuted = TRUE)
  os1 <- rstan::extract(Rasch[[2]],permuted = TRUE)
  is2 <- rstan::extract(PL2[[1]],permuted = TRUE)
  os2 <- rstan::extract(PL2[[2]],permuted = TRUE)
  is3 <- rstan::extract(PL3[[1]],permuted = TRUE)
  os3 <- rstan::extract(PL3[[2]],permuted = TRUE)
  
  # RQR calculation function
  # input: probability, response / output: RQR
  RQR <- function(p,o){
    if(o==0){
      rqr <- map(p,function(x){qnorm(x)})
    }else{
      rqr <- map(p,function(x){qnorm(x)})
    }
    unlist(rqr)
  }
  
  # in-sample
  # prepare sampling parameters
  beta1 <- is1$beta[dim(is1$beta)[1],]
  theta1 <- is1$theta[dim(is1$theta)[1],]
  alpha2 <- is2$alpha[dim(is2$alpha)[1],]
  beta2 <- is2$beta[dim(is2$beta)[1],]
  theta2 <- is2$theta[dim(is2$theta)[1],]
  alpha3 <- is3$alpha[dim(is3$alpha)[1],]
  beta3 <- is3$beta[dim(is3$beta)[1],]
  theta3 <- is3$theta[dim(is3$theta)[1],]
  gamma3 <- is3$gamma[dim(is3$gamma)[1],]
  # calculate RQR
  p1_in_sample <- matrix(nrow=testn,ncol=k)
  p2_in_sample <- matrix(nrow=testn,ncol=k)
  p3_in_sample <- matrix(nrow=testn,ncol=k)
  for(i in 1:testn){
    for(j in 1:k){
      # Rasch
      p1_in_sample[i,j] <- IRT(1,beta1[j],theta1[testset[i]],0)
      # 2PL
      p2_in_sample[i,j] <- IRT(alpha2[j],beta2[j],theta2[testset[i]],0)
      # 3PL
      p3_in_sample[i,j] <- IRT(alpha3[j],beta3[j],theta3[testset[i]],gamma3[j])
    }
  }

  ksModel <- ksIRT(response1,key = 1,format = 1,kernel = "gaussian",bandwidth = "Silverman")
  #plot(ksModel,items=c(item))
  item_ooc_y <- ksModel[["OCC"]][seq(1,dim(response1)[2]*2,2),4:54]
  plot(ksModel[["expectedscores"]],item_ooc_y[item,],type="l",lwd=2,xlab="expected score",ylab="probability")
  points(rowSums(p1_in_sample),p1_in_sample[,item],pch = 15, cex=0.5,col="green")
  points(rowSums(p2_in_sample),p2_in_sample[,item],pch = 15, cex=0.5,col="red")
  points(rowSums(p3_in_sample),p3_in_sample[,item],pch = 15, cex=0.5,col="blue")
  title(paste("item",as.character(item)))
  legend(0,1,c("Rasch","2PL","3PL"), lwd = c(3,1), col=c("green","red","blue"), y.intersp=1.5)
}
#item_plot(9,testset,result_Rasch,result_2PL,result_3PL,dat,dat,dat)
item_plot(7,testset,result_Rasch,result_2PL,result_3PL,PL2_response,PL2_response,PL2_response)

