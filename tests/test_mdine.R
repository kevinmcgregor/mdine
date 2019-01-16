#Test mdine

set.seed(893497)
n=50
n.spec = 5

x = cbind(1, rbinom(n, 1, 0.5))
w = matrix(rnorm(n*n.spec), n, n.spec)

beta = rbind(c(5, 4, 3, 2, 1), c(-1, 1, 0, 1, 2))

exp_part = cbind(exp(x%*%beta+w), 1)
probs = exp_part/rowSums(exp_part)

rd = floor(rnorm(n, 50000, 5000))

counts = matrix(0, n, n.spec+1)
for (i in 1:n) {
  counts[i,] = rmultinom(1, rd[i], probs[i,])
}

library(mdine)

md_test = mdine(Y = counts, X = x, Z = x[,2], mc.cores = 4, iter = 500)


