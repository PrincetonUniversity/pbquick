library(poibin)
library(poisbinom)
library(microbenchmark)


set.seed(831213)
phis <- runif(5000)
eval.vec <- 0:5000

#Test dpoisbinom
cat("Test probability mass function:\n")
identical(round(dpoibin(eval.vec,phis), 6), round(dpoisbinom(eval.vec,phis, FALSE), 6))
microbenchmark(dpoibin(eval.vec,phis)
              ,dpoisbinom(eval.vec, phis, FALSE),
               times = 5L)

#Test ppoisbinom
cat("Test cumulative distribution function:\n")
identical(round(ppoibin(eval.vec,phis), 6), round(ppoisbinom(eval.vec,phis, TRUE, FALSE), 6))
microbenchmark(ppoibin(eval.vec,phis)
               ,ppoisbinom(eval.vec, phis, TRUE, FALSE),
               times = 5L)

#Test qpoisbinom
cat("Test quantile function:\n")
identical(round(qpoibin(phis,phis), 6), round(qpoisbinom(phis,phis), 6))
microbenchmark(qpoibin(phis,phis)
               ,qpoisbinom(phis, phis),
               times = 5L)
