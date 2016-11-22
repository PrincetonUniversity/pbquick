library(poibin)
library(poisbinom)
library(microbenchmark)


phis <- runif(1000)
eval.vec <- 0:1000

identical(round(dpoibin(eval.vec,phis), 6), round(dpoisbinom(eval.vec,phis), 6))
  
microbenchmark(dpoibin(eval.vec,phis)
              ,dpoisbinom(eval.vec,phis),
               times = 5L)
