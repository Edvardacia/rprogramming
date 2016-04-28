library("HiddenMarkov")
library("quantmod")

Pi <- matrix(c(1/2, 1/2, 0, 1/3, 1/3, 1/3, 0, 1/2, 1/2), byrow=TRUE, nrow=3)
delta <- c(0, 1, 0)
x <- dthmm(NULL, Pi, delta, "norm", list(mean=c(1, 6, 3), sd=c(0.5, 1, 0.5)))
x <- simulate(x, nsim=1000)
# use above parameter values as initial values
y <- BaumWelch(x)
print(summary(y))
print(logLik(y))
hist(residuals(y))
# check parameter estimates
print(sum(y$delta))
print(y$Pi %*% rep(1, ncol(y$Pi)))

# m*m matrix, veroitnost perehoda odnorodnoy skritoy markovskoy sepi

Pi <- matrix(c(1/2, 1/2, 0, 0, 0, 1/3, 1/3, 1/3, 0, 0, 0, 1/3, 1/3, 1/3, 0, 0, 0, 1/3, 1/3, 1/3, 0, 0, 0, 1/2, 1/2), byrow=TRUE, nrow=5)

delta <- c(0, 1, 0, 0, 0)
lambda <- c(1, 4, 2, 5, 3)
m <- nrow(Pi)
x <- dthmm(NULL, Pi, delta, "pois", list(lambda=lambda), discrete=TRUE)
x <- simulate(x, nsim=2000)
states <- Viterbi(x)
states <- factor(states, levels=1:m)

p <- matrix(NA, nrow = m, ncol = m)

for (j in 1:m) { 
  a <- (x$y==j) 
  p[j,] <- table(states[a])/sum(a) }

print(p)

#------ Local Decoding ------
# locally decode at i=100
print(which.max(Estep(x$x, Pi, delta, "pois", list(lambda=lambda))$u[100,]))
#---------------------------------------------------
# simulate a beta HMM
Pi <- matrix(c(0.8, 0.2,
               0.3, 0.7),
             byrow=TRUE, nrow=2)
delta <- c(0, 1)
y <- seq(0.01, 0.99, 0.01)
plot(y, dbeta(y, 2, 6), type="l", ylab="Density", col="blue")
points(y, dbeta(y, 6, 2), type="l", col="red")
n <- 100
x <- dthmm(NULL, Pi, delta, "beta",
           list(shape1=c(2, 6), shape2=c(6, 2)))
x <- simulate(x, nsim=n)
# colour denotes actual hidden Markov state
plot(1:n, x$x, type="l", xlab="Time", ylab="Observed Process")
points((1:n)[x$y==1], x$x[x$y==1], col="blue", pch=15)
points((1:n)[x$y==2], x$x[x$y==2], col="red", pch=15)
states <- Viterbi(x)
# mark the wrongly predicted states
wrong <- (states != x$y)
points((1:n)[wrong], x$x[wrong], pch=1, cex=2.5, lwd=2)






Pi <- matrix(c(1/2, 1/2, 0, 0, 0, 1/3, 1/3, 1/3, 0, 0, 0, 1/3, 1/3, 1/3, 0, 0, 0, 1/3, 1/3, 1/3, 0, 0, 0, 1/2, 1/2), byrow=TRUE, nrow=5)
delta <- c(0, 1, 0, 0, 0)
lambda <- c(1, 4, 2, 5, 3)
m <- nrow(Pi)
x <- dthmm(NULL, Pi, delta, "pois", list(lambda=lambda), discrete=TRUE)

bull1 = rnorm ( 100, 1, 5 )
bear = rnorm ( 100, -5, 1 )
bull2 = rnorm (100, 4, 4 )
x <- c( bull1, bear, bull2 )
states <- Viterbi(x)


