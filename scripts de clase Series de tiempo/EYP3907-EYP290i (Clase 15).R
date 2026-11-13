##############
## Clase 15 ##
##############

## Capitulo 5: Metodos de Estimacion (parte 2)
## - Metodo de momentos: Yule-Walker
## - Maxima verosimilitud: exacta (AR(1) a mano) e innovaciones
##   (Durbin-Levinson para un ARMA(p,q) general)
## - Whittle: la verosimilitud construida desde el periodograma
## - Inferencia: normalidad asintotica y test z de los coeficientes
## (continuacion de la Clase 14: AIC/BIC, media y bandas de Bartlett)

## Requiere: LSTS, forecast

## La serie de trabajo es la misma de la Clase 14
X <- LSTS::malleco
n <- length(X)

##################################################################
## Paso 1: metodo de momentos --> Yule-Walker                   ##
##################################################################

## Caso AR(p): igualar gamma(h) teorica y muestral para h = 1,...,p
## entrega el sistema lineal G phi = g

ar.yule.walker <- function(X, p = 1){
  g <- (acf(X, type = "covariance", plot = FALSE, lag.max = p)$acf[,,1])[2:(p+1)]
  G <- matrix(NA, ncol = p, nrow = p)
  g.aux <- (acf(X, type = "covariance", plot = FALSE, lag.max = p)$acf[,,1])[1:p]
  aux2 <- c(rev(g.aux)[seq_len(p-1)], g.aux)  ## seq_len evita el caso p = 1
  for(i in 1:p){
    G[i,] <- aux2[(p-(i-1)):((p-(i-1)) + p - 1)]
  }
  list(ar = c(solve(G) %*% g))
}

p <- 10
M <- cbind(ar.yule.walker(X, p = p)$ar,
           ar(X, aic = FALSE, order.max = p)$ar,
           ar.yw(X, aic = FALSE, order.max = p)$ar)
colnames(M) <- c("ar.yule.walker", "ar", "ar.yw")
rownames(M) <- paste("ar", 1:p, sep = "")
round(M, 4)
## nuestra funcion coincide con las implementaciones de R

##################################################################
## Paso 2: maxima verosimilitud                                 ##
##################################################################

## (a) AR(1): verosimilitud EXACTA escrita a mano
##     y1 ~ N(0, sigma^2/(1-phi^2));  yt | y(t-1) ~ N(phi y(t-1), sigma^2)
loglik <- function(x, serie){
  phi <- x[1]; sigma <- x[2]
  Y <- serie; n <- length(Y)
  -(sum(log(dnorm(Y[2:n], mean = phi*Y[1:(n-1)], sd = sigma))) +
    log(dnorm(Y[1], mean = 0, sd = sigma/sqrt(1-phi^2))))
}

fit1 <- nlminb(objective = loglik, start = c(0,1), serie = X - mean(X),
               lower = c(-0.99, 1e-6), upper = c(0.99, Inf))
fit2 <- ar(X - mean(X), aic = FALSE, order.max = 1)
fit3 <- forecast::Arima(X - mean(X), order = c(1,0,0), include.mean = FALSE)
M <- rbind(c(fit1$par),
           c(fit2$ar, sd(fit2$resid, na.rm = TRUE)),
           c(fit3$coef, sqrt(fit3$sigma2)))
colnames(M) <- c("ar1", "sigma"); rownames(M) <- c("nlminb", "ar", "Arima")
round(M, 4)

## (b) Durbin-Levinson para un ARMA(p,q) GENERAL: la verosimilitud SIN
##     determinantes ni inversas. La recursion solo necesita la ACF
##     teorica del modelo (ARMAacf) y entrega, para cada t, el predictor
##     a un paso X.hat[t] y la varianza (relativa a gamma(0)) de su
##     error, nu[t]. Los phi[t,t] son la PACF teorica del modelo.
DurbinLevinson <- function(serie, ar = NULL, ma = NULL){
  X <- serie; N <- length(X)
  rho <- ARMAacf(ar = ar, ma = ma, lag.max = N)[-1]   ## rho(1),...,rho(N)
  Phi <- matrix(0, N, N); nu <- numeric(N); X.hat <- numeric(N)
  Phi[1,1] <- rho[1]
  nu[1] <- 1 - Phi[1,1]^2
  X.hat[1] <- 0
  X.hat[2] <- Phi[1,1]*X[1]
  for(n in 2:(N-1)){
    Phi[n,n] <- (rho[n] - sum(Phi[n-1,1:(n-1)]*rho[(n-1):1]))/nu[n-1]
    Phi[n,1:(n-1)] <- Phi[n-1,1:(n-1)] - Phi[n,n]*Phi[n-1,(n-1):1]
    nu[n] <- nu[n-1]*(1 - Phi[n,n]^2)
    X.hat[n+1] <- sum(X[n:1]*Phi[n,1:n])
  }
  ## nu[t] = Var(X[t] - X.hat[t])/gamma(0); para t = 1 es 1 (X.hat[1] = 0)
  list(fitted = X.hat, nu = c(1, nu[1:(N-1)]), pacf = diag(Phi)[1:(N-1)])
}

## Verosimilitud gaussiana exacta via innovaciones e[t] = X[t] - X.hat[t]:
##   L = -N/2 log(2 pi) - 1/2 sum log(gamma0 nu[t]) - 1/2 sum e[t]^2/(gamma0 nu[t])
## gamma(0) se "perfila" (maximiza en forma cerrada): gamma0 = mean(e^2/nu).
## Asi la funcion depende SOLO de los coeficientes ARMA (menos parametros
## que optimizar), y sigma^2 se recupera al final con gamma0 = sigma^2 sum psi_j^2.
loglik.DL <- function(x, serie, p = 1, q = 0){
  ar <- if(p > 0) x[1:p] else NULL
  ma <- if(q > 0) x[p + 1:q] else NULL
  ## fuera de la region estacionaria/invertible la recursion no tiene
  ## sentido: devolvemos un valor enorme para que nlminb se aleje
  if(p > 0 && any(Mod(polyroot(c(1, -ar))) <= 1)) return(1e10)
  if(q > 0 && any(Mod(polyroot(c(1,  ma))) <= 1)) return(1e10)
  fit <- DurbinLevinson(serie, ar = ar, ma = ma)
  e <- serie - fit$fitted; nu <- fit$nu; N <- length(serie)
  gamma0 <- mean(e^2/nu)
  -(-N/2*log(2*pi) - N/2*log(gamma0) - 0.5*sum(log(nu)) - N/2)  ## -logLik
}

sigma2.DL <- function(x, serie, p = 1, q = 0){
  ar <- if(p > 0) x[1:p] else NULL
  ma <- if(q > 0) x[p + 1:q] else NULL
  fit <- DurbinLevinson(serie, ar = ar, ma = ma)
  gamma0 <- mean((serie - fit$fitted)^2/fit$nu)
  gamma0/(1 + sum(ARMAtoMA(ar, ma, lag.max = 1000)^2))   ## gamma0 / sum psi_j^2
}

## Chequeo 1: en el MLE del AR(1) reproduce el logLik de Arima
cbind(DL = -loglik.DL(fit3$coef[1], X - mean(X), p = 1, q = 0),
      Arima = as.numeric(logLik(fit3)))

## Chequeo 2: estimacion ML de varios ARMA(p,q) con nuestra verosimilitud,
## comparada con forecast::Arima (mismos coeficientes, sigma^2 y logLik).
## Punto de partida: Yule-Walker del Paso 1 para la parte AR y 0 para la
## MA (con un ARMA sobreparametrizado la verosimilitud tiene maximos
## locales y partir de 0 puede quedarse en uno de ellos)
Y <- X - mean(X)
for(pq in list(c(1,0), c(1,1), c(2,1))){
  p <- pq[1]; q <- pq[2]
  st <- c(ar.yule.walker(Y, p = p)$ar, rep(0, q))
  f <- nlminb(objective = loglik.DL, start = st, serie = Y, p = p, q = q)
  a <- forecast::Arima(Y, order = c(p,0,q), include.mean = FALSE)
  M <- rbind(DL = c(f$par, sigma2.DL(f$par, Y, p, q), -f$objective),
             Arima = c(a$coef, a$sigma2, a$loglik))
  colnames(M) <- c(names(a$coef), "sigma2", "logLik")
  cat("\nARMA(", p, ",", q, ")\n", sep = ""); print(round(M, 4))
}
## los phi[t,t] de la recursion son la PACF, y nu[t] la varianza del
## error de prediccion a un paso: predecir y evaluar verosimilitud
## son el mismo calculo (asi trabaja arima, via Kalman)
cbind(DL = DurbinLevinson(Y, ar = 0.6, ma = 0.4)$pacf[1:5],
      teorica = ARMAacf(ar = 0.6, ma = 0.4, lag.max = 5, pacf = TRUE)[1:5])

##################################################################
## Paso 3: Whittle --> la verosimilitud desde el periodograma   ##
##################################################################

## Los datos entran SOLO via I(lambda) (Clase 13); el modelo, solo via
## f(lambda). Whittle busca la f_theta que mejor "calza" con I:
##   L(theta) = sum [ log f(lambda_j) + I(lambda_j)/f(lambda_j) ]

whittle <- function(x, serie){
  phi <- x[1]; sigma <- x[2]
  aux <- LSTS::periodogram(serie, plot = FALSE)
  f <- LSTS::spectral.density(ar = phi, sd = sigma, lambda = aux$lambda)
  sum(log(f) + aux$periodogram/f)
}
round(nlminb(objective = whittle, start = c(0,1), serie = X - mean(X),
             lower = c(-0.99, 1e-4), upper = c(0.99, Inf))$par, 4)
## comparar con la ML exacta del Paso 2: practicamente identicos.
## LSTS::LS.whittle implementa esta idea (en version localmente
## estacionaria; con coeficientes constantes es este caso)

## Experimento Monte Carlo: ECM de Whittle vs ML exacta
B <- 200                     ## suba B si hay tiempo en clase
M <- c()
for(i in 1:B){
  Y <- arima.sim(n = 200, model = list(order = c(1,0,0), ar = 0.6))
  phi1 <- nlminb(objective = whittle, start = c(0,1), serie = Y - mean(Y),
                 lower = c(-0.99, 1e-4), upper = c(0.99, Inf))$par[1]
  phi2 <- forecast::Arima(Y, order = c(1,0,0), include.mean = FALSE)$coef[1]
  M <- rbind(M, c(whittle = phi1, ML = phi2))
}
Sesgo <- apply(M, 2, mean) - 0.6
ECM <- Sesgo^2 + apply(M, 2, var)
rbind(Sesgo, ECM)
## mismo sesgo y ECM casi identico: Whittle entrega la calidad de la
## ML exacta a una fraccion del costo computacional (clave en ARFIMA)

##################################################################
## Paso 4: inferencia --> normalidad asintotica y test z        ##
##################################################################

## Teorema: sqrt(n) (theta.hat - theta0) --> N(0, Gamma^-1(theta0)).
## Para el AR(1): Var asintotica = (1 - phi^2)/n

set.seed(3907)
phi <- 0.7; n <- 200
ph <- replicate(1000, {
  Y <- arima.sim(model = list(ar = phi), n = n)
  as.numeric(arima(Y, order = c(1,0,0), include.mean = FALSE)$coef[1])
})
par(bty = "n", las = 1)
hist(ph, breaks = 40, freq = FALSE, col = "gray90", border = "gray60",
     xlab = expression(hat(phi)), main = "")
curve(dnorm(x, phi, sqrt((1-phi^2)/n)), add = TRUE, col = "red", lwd = 2)
c(empirica = sd(ph), teorica = sqrt((1-phi^2)/n))

## En la practica: los s.e. que reporta Arima salen de Gamma^-1(theta),
## y con ellos el test z de significancia de cada coeficiente
fit <- forecast::Arima(X - mean(X), order = c(1,0,0), include.mean = FALSE)
se <- sqrt(diag(fit$var.coef))
z  <- fit$coef/se
cbind(coef = fit$coef, se = se, z = z, p.value = 2*pnorm(-abs(z)))
## |z| > 1.96: el coeficiente es significativo al 5%

## Moraleja:
## - Yule-Walker: el estimador de momentos del AR(p) es un sistema
##   lineal con la ACVF muestral; con el 1/n de acf() siempre tiene
##   solucion causal.
## - Durbin-Levinson evalua la verosimilitud gaussiana exacta solo con
##   predicciones a un paso; Whittle la aproxima desde el periodograma.
## - El MLE es asintoticamente normal y eficiente: de ahi salen los
##   s.e., los IC y los test z que reporta Arima.
