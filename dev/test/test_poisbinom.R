library(poibin)
library(poisbinom)
library(microbenchmark)


set.seed(831213)
phis <- runif(5001)
eval.vec <- 0:5001

## Test dpoisbinom
cat("\n")
cat("Test probability mass function:\n")
cat("   Results identical to 'poibin' package?")
identical(round(dpoibin(eval.vec,phis), 6), round(dpoisbinom(eval.vec,phis, FALSE), 6))
cat("\n")
cat("Benchmark test for computing PMF:\n")
microbenchmark(dpoibin(eval.vec,phis)
              ,dpoisbinom(eval.vec, phis, FALSE),
               times = 5L)

## Test ppoisbinom
cat("\n")
cat("Test cumulative distribution function:\n")
cat("   Results identical to 'poibin' package?")
identical(round(ppoibin(eval.vec,phis), 6), round(ppoisbinom(eval.vec,phis, TRUE, FALSE), 6))
cat("\n")
cat("Benchmark test for computing CDF:\n")
microbenchmark(ppoibin(eval.vec,phis)
               ,ppoisbinom(eval.vec, phis, TRUE, FALSE),
               times = 5L)

## Test qpoisbinom
cat("\n")
cat("Test quantile function:\n")
cat("   Results identical to 'poibin' package?")
identical(round(qpoibin(phis,phis), 6), round(qpoisbinom(phis,phis), 6))
cat("\n")
cat("Benchmark test for computing quantile function:\n")
microbenchmark(qpoibin(phis,phis)
               ,qpoisbinom(phis, phis),
               times = 5L)

## Test rpoisbinom
cat("\n")
cat("Test random number generator:\n")
cat("   Results identical to 'poibin' package?")
set.seed(45692); out.rpoibin <- as.integer(rpoibin(400, phis))
set.seed(45692); out.rpoisbinom <- rpoisbinom(400, phis)
identical(out.rpoibin, out.rpoisbinom)
cat("\n")
cat("Benchmark test for generating random draws:\n")
microbenchmark(rpoibin(1000, phis)
              ,rpoisbinom(1000, phis),
               times = 5L)
     
