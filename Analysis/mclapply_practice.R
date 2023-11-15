# Create a folder in your directory called testruns

library(parallel)

# Write a single function which contains your full pipeline
# The first argument will be the thing you want to iterate over, i.e., an index
fcn.test <- function(x,n){
  set.seed(x)
  x.df <- matrix(runif(n), nrow = 10)
  save(file = paste0(save_path, "testruns/test_res_", x, ".dat"), x.df)
}

n <- 10000
s1 <- system.time(mapply(function(i) fcn.test(x=i,n=n), 1:5))
s2 <- system.time(mclapply(1:5, function(i) fcn.test(x=i, n=n), mc.cores = 5))



# Fancier version: now iterate over different probabilities instead of chains

# Create list of probabilities to iterate over
p1.site <- c(0.5, 0.7, 0.9)
p2.hab.c <- seq(0.2, 0.8, by = 0.2)
p3.hab.r <- seq(0.1, 0.4, by = 0.2)
p4.h <- seq(0.2, 0.4, by = 0.2)
p5.c <- c(0.01, 0.05, 0.1)
p6.r <- c(0.01, 0.05, 0.1)

prob.list <- list() 
prob.list[[1]] <- c(1, 0.75, 0.5, 0.45, 0.25, 0.1, 0.05)
counter <- 2
# First exercise, just create nested loops that create a grid of probabilities
for(p1 in p1.site){
  for(p2 in p2.hab.c[p2.hab.c<p1]){
    for(p3 in p3.hab.r[p3.hab.r< p2]){
      for(p4 in p4.h[p4.h<p3]){
        for(p5 in p5.c[p5.c<p4]){
          for(p6 in p6.r[p6.r<p5]){
            prob.list[[counter]] <- c(1, p1, p2, p3, p4, p5, p6)
            counter <- counter + 1
          }
        }
      }
    }
    
  }
}


# Write a single function which contains your full pipeline
# The first argument will be the thing you want to iterate over, i.e., an index

fcn.test.2 <- function(x,n,n.cv,prob.list){
  p <- prob.list[[x]]
  means.df <- matrix(NA, nrow = n.cv, ncol = length(p))
  for(rr in 1:n.cv){
    x.df <- t(matrix(rbinom(n*length(p),1,p), nrow = length(p)))
    means.df[rr,] <- colMeans(x.df) 
  }
  colMeans(means.df)
  #save(file = paste0(save_path, "testruns/test_res2_", x, ".dat"), x.df)
}

n.cv <- 2
n <- 1000
t(mapply(function(i) fcn.test.2(x=i,n=n, n.cv,prob.list), 1:5))

s1 <- system.time(mapply(function(i) fcn.test.2(x=i,n=n, prob.list), 1:5))
s2 <- system.time(mclapply(1:5, function(i) fcn.test.2(x=i, n=n, prob.list, savemeans), mc.cores = 5))

