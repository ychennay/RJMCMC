#derived from Ai Jialin's code in his master's thesis work:
#http://www1.maths.leeds.ac.uk/~voss/projects/2011-RJMCMC/AiJialin.pdf

readgeno=function(fname){
  r=readLines (fname)
  x=matrix(as.numeric(unlist(strsplit (r,split =""))),nrow=length(r), byrow=T)
  return (x)
}

library("boot")
getwd()
setwd("/Users/ychen/Desktop")
genotype= readgeno("geno.geno")
genotype = t(genotype)
dimensions <- dim(genotype)

readgeno=function(fname){
  r=readLines (fname)
  x=matrix(as.numeric(unlist(strsplit (r,split =""))),nrow=length(r), byrow=T)
  return (x)
}


# disaster number function m(.)
count.SNPs <- function(s1, s2) {
  sum_vector <- list()
  for (i in 1:dimensions[2]){
    sum_vector[i] <-   sum(genotype[,i] >= s1 & genotype[,i] < s2)
  }
  return(unlist(sum_vector))}

k0 <- 6
s0 <- c(1,1100,1200,1550,1600,1606,1614)
h0 <- c(1,10,20,30,40,50)

lambda = 7 #this value is in the generation of prior distributions of likelihood for k

# value of Gamma(a,b), used for the prior distribution of heights
a=3 
b=1

birth <- function(change_points){return(min(1,lambda/(change_points+1)))}
death <- function(change_points){return(min(1,change_points/lambda))}


Move <- function(n, h, s, k) {
  H <- list()
  S <- list()
  K <- vector(length=n) # K is a vector the length of the number of iterations
  K[1] <- k
  H[[1]] <- h
  S[[1]] <- s
  
  for (i in 2:n){
    position_prob <- ifelse(k<=1, 0, 0.5 * (1 - birth(k-1) - death(k-1)))
    #if (k<=1) {position_prob <- 0} else {position_prob <- 0.5 * (1 - birth(k-1) - death(k-1))}
    
    height_prob <-( 1 - birth(k-1) - death(k-1) - position_prob) #probability of height is
    #defined as 1 - probability of birth from the previous model - probability of death from
    #the previous model - position probability
    
    type <- runif(1) #a randomly generated number between 0 and 1
    
    if (type>1-height_prob) {#from 1 to k is for Height h_1 to h_k
      j <- sample(1:k,size=1) #pick a number between 1 and K, with K being the total number of
      #available numbers at the current iteration
      u <- runif(1,-0.5,0.5) #generate a u variable
      h_tilde <- h[j] * exp(u) #use the u varible to generate a proposal for a new height
      U <- runif (1) # M-H algorithm
      if (U < alpha.hmove(s,h,h_tilde,j)) {
        h[j] <- h_tilde #accept the update
      }
    }
    
    if (type<=1-height_prob && type>1-height_prob-position_prob) {#from 2 to k is for Position s_2 to s_k
      j <- sample(1:(k-1), size=1) + 1 #pick a number between 1:k
      s_tilde <- runif(1,s[j-1],s[j+1]) #produce a random number from the previous model start to the proposed model's end
      U <- runif (1)
      if (U < alpha.smove(s,h,s_tilde,j)) { #if alpha is> U / M-H algorithm
        s[j] <- s_tilde # update s
      }
    }
    ########################################################################################################################################
    if (type>=birth(k-1) && type<=birth(k-1)+death(k-1)) {#from 1 to k-1 is for death of steps d_2 to d_k
      #conditional is true when type is a probability greater than birth but less than both birth
      #and death
      j <- sample(1:(k-1), size=1) # pick a number between 1 and k-1
      r <- (s[j+2] - s[j+1]) / (s[j+2] - s[j])
      hj_prime <- h[j]^(1-r) * h[j+1]^r  #exp((1-r) * log(h[j]) + r * log(h[j+1]))
      U <- runif(1)
      if (U < alpha.death(s,h,hj_prime,j,k)){#Metropolis Hastings Algorithm
        k <- k - 1
        h[j] <- hj_prime
        h <- h[-(j+1)]
        s <- s[-(j+1)]
      }
    }
    
    ########################################################################################################################################
    
    if (type<=birth(k-1)){#from 1 to k is for birth of steps b_1 to b_k
      # if the type variable generated is less than the probability of birth from prev model
      j <- sample(1:k,size=1)
      s_star <- runif(1,s[j],s[j+1]) # generate a proposed position between the current model
      # and the previous model position
      u <- runif(1)
      r <- exp(log(s[j+1] - s_star) - log(s[j+1] - s[j])) # ratio calculation for proposed position
      h1_prime = h[j] * exp(r * (log(u) - log(1-u)))
      h2_prime = h[j] * exp((1-r) * (log(1-u) - log(u)))
      U <- runif(1)
      if (U < alpha.birth(s,h,s_star,h1_prime,h2_prime,j,k)){ # if birth happens
        s <- c(s[1:j], s_star, s[(j+1):(k+1)]) # the s_star proposed positio is inserted into the
        # middle of the list of positions
        if (j > 1) { # if j is greater than 1, the first model
          left <- h[1:(j-1)] #  populate left variable with h from 1: j-1
        } else {# if j is one, then there is nothing to the left!
          left <- c()
        }
        if (j < k) { # assuming that j is less than k
          right <- h[(j+1):k] # populate the right variable with the new inputs
        } else {
          right <- c()# if j is greater than k, then there is nothing to the right
        }
        h <- c(left, h1_prime, h2_prime, right) # concat the left, right, and h primes
        k <- k + 1 # update k, since we have added a model
      }
    }
    
    
    #updates
    K[i] <- k
    H[[i]] <- h
    S[[i]] <- s
  }
  
  #here the algorithm finally exits from the for loop and returns results
  return(list(h=H,k=K,s=S))
}
X <- Move(20000, h0, s0, k0) #################################################################################
K <- X$k
H <- X$h
S <- X$s

h_list = list()
for (i in 1:3999){
  h_list[i] = H[[i]][1]
}

h_list2 = list()
for (i in 1:3999){
  h_list2[i] = H[[i]][2]
}

s_list = list()
for (i in 1:3999){
  s_list[i] = S[[i]][1]
}

B <- unlist(s_list)
A <- unlist(s_list)
C <- unlist(s_list)

unlist(h_list)


hist(K,breaks=seq(0.5, max(K)+0.5, by=0.5)+0.20, main = "Posterior Distribution of SNP Clusters, N = 20000")
plot(K, type='l', main="Trace Plot for K Convergence, N = 20000", xlab="Iterations", ylab="K Value", xlim = c(0,20000))
plot(S)
hist(S)
plot(unlist(h_list), type='l', main="Plot Trace of Height for Clusters", ylab="Height", xlab="Iterations")
lines(unlist(h_list2), type='l', col='red')

plot(C, type='l', col='blue', main="Plot Trace of Position for Clusters", ylab="Position", xlab="Iterations")
lines(A, col = 'pink')
lines(B, col = 'orange')

write.csv(H,"Heights.csv")
non.null.list <- lapply(H, Filter, f = Negate(is.null))
