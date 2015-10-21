load("MonocarpSimData.Rdata")

LU <- matrix(NA, 10000, 2)
for (k in 1:10000) {
  sampleSizes <- sample(sim.data$size, 1000, replace = FALSE)
  LU[k, ] <- range(sampleSizes)
}
cat(apply(LU, 2, mean), "\n")

quantile(LU[, 1])
quantile(LU[, 2])