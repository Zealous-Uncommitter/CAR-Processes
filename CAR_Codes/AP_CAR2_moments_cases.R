set.seed(123)

NewCAR2 <- function(a1, a2, b0, N, K, noise, Y0 = 0, Y1 = 0)
{
  h <- 1 / K
  
  Y <- numeric(N * K + 100001)
  
  Y[1] <- Y0
  Y[2] <- Y1
  
  for(i in 2:(length(Y)-1))
  {
    Y[i+1] <-
      2 * Y[i] -
      Y[i-1] -
      a1 * h * (Y[i] - Y[i-1]) -
      a2 * h^2 * Y[i] +
      b0 * h * noise[i]
  }
  
  Y[100001:length(Y)]
}

# Parameters
N = 500
K = 10000
M = 100
sigma = 1
mu = 3
a1 = 8
a2 = 4
b0 = 1
s = 1

cases <- data.frame(a1 = c(2, 3, 1.5,  4, 0.5, 2,  2, 4),
                    a2 = c(0.9, 1, 0.5,  8, 2, 5,  1, 4))

for (i in 1:nrow(cases)) {
  
  ###Simulate and Sample Process###
  
  a1 <- cases$a1[i]; a2 <- cases$a2[i]
  Z <- rnorm(N*K + 100001, mean = mu*(1/K), sd = sigma*sqrt(1/K))
  Y <- NewCAR2(a1, a2, 1, N, K, noise = Z, Y0 = 0, Y1 = 0)   
  relation <- K / M
  
  sample_index <- seq(
    1,
    length(Y),
    by = relation
  )
  Ysample <- Y[sample_index]
  
  ###Theoretical moment values###
  
  cases$theoretical_mean[i] <- mu/a2
  cases$theoretical_variance[i] <- eta^2/(2*a2*a1)
  l1 <- (-a1 + sqrt(as.complex(a1^2-(4*a2))))/2
  l2 <- (-a1 - sqrt(as.complex(a1^2-(4*a2))))/2
  
  if (l1 != l2) {cases$theoretical_covariance[i] <- Re(eta^2/(2*a2*a1*(l1-l2)) * (-l2*exp(l1*s)+l1*exp(l2*s)))}
  if (l1 == l2) {cases$theoretical_covariance[i] <- Re(eta^2 * (-l2 * (s * exp(l2 * s)) + exp(l2 * s)) / (-4 * l2^3))}
  
  ###Simulated moment values###
  
  k <- s * M   # lag s in time units = s*M positions in Ysample
  n <- length(Ysample)
  cases$simulated_mean[i] <- mean(Ysample)
  cases$simulated_variance[i] <- var(Ysample)
  cases$simulated_covariance[i] <- cov(Ysample[1:(n-k)], Ysample[(1+k):n])
}

print(round(cases, 4))


