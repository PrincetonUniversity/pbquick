library(poibin)
library(pbquick)
library(microbenchmark)


phis <- runif(1000)
eval.vec <- 0:1000

identical(round(dpoibin(eval.vec,phis), 6), round(dpbquick(eval.vec,phis), 6))
  
microbenchmark(dpoibin(eval.vec,phis)
              ,dpbquick(eval.vec,phis),
               times = 5L)
