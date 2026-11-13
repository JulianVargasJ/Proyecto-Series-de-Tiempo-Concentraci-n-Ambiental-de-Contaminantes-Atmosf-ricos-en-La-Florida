#############################
## Algoritmo de Innovación ##
#############################

algoritmo.innov = function(X, ar = NULL, ma = NULL){
N = length(X)
Theta = matrix(0, ncol = N, nrow = N)
rho = ARMAacf(ar = ar, ma = ma, lag.max = N)[2:(N+1)]
X.hat = c()
v = rep(0,N)
v0 = var(X) 
X.hat[1] = 0 
Theta[1,1] = rho[1]
v[1] = var(X)-Theta[1,1]^2*v0
for(n in 2:N){
Theta[n,n] = rho[n]
for(k in 1:(n-1)){
j = 0:(k-1)
Theta[n,n-k] = (rho[n-k]*v0-sum(c(Theta[k, k-j])*Theta[n, n-j]*c(v0,v[j])))/v[k]
}
v[n] = var(X) - sum(Theta[n,n-0:(n-1)]^2*c(v0,v[1:(n-1)]))
}

for(n in 1:N-1){
X.hat[n+1] = sum(Theta[n,1:n] * (X[n+1-1:n]-X.hat[n+1-1:n]))
}
X.hat
}