##############
## Clase 14 ##
##############

## Capitulo 5: Metodos de Estimacion (parte 1)
## - Parsimonia: seleccion de orden via AIC/BIC
## - Estimacion de la media: IC correcto para series correlacionadas
## - ACF muestral: bandas de Bartlett y ajuste por momentos del AR(1)
## - Diagnostico de residuos del ajuste
## (continua en la Clase 15: Yule-Walker, maxima verosimilitud,
##  Whittle e inferencia)

## Requiere: LSTS, forecast

########################################################
## Paso 1: parsimonia --> AIC y BIC                   ##
########################################################

## AIC = -2 log L + 2r ; BIC = -2 log L + r log(n), con r = p + q + 1.
## Se prefiere el modelo con el criterio MAS PEQUENO; BIC penaliza mas.

set.seed(3907)
X <- arima.sim(model = list(ar = c(0.6, -0.4)), n = 300)
plot(X)
P <- 0:6
crit <- t(sapply(P, function(p){
  f <- arima(X, order = c(p,0,0), include.mean = FALSE, method = "ML")
  r <- p + 1                       ## p coeficientes + sigma
  c(AIC = -2*as.numeric(logLik(f)) + 2*r,
    BIC = -2*as.numeric(logLik(f)) + r*log(length(X)))
}))
round(cbind(p = P, crit), 2)

par(bty = "n", las = 1, mfrow = c(1,1))
matplot(P, crit, type = "b", pch = 19, lty = 1, lwd = 2,
        col = c("blue", "red"), xlab = "orden p", ylab = "criterio")
legend("topright", legend = c("AIC", "BIC"), col = c("blue","red"),
       pch = 19, lwd = 2, bty = "n")
## ambos criterios recuperan el orden verdadero p = 2; desde ahi cada
## parametro extra mejora tan poco la verosimilitud que la penalizacion
## domina. Esto es lo que auto.arima hace por nosotros (Clase 11):
forecast::auto.arima(X, d = 0, allowmean = FALSE, stepwise = FALSE)

##################################################################
## Paso 2: estimacion de la media --> el IC iid queda CORTO     ##
##################################################################

## Var(y.barra) ~ v/n con v = 2 pi f(0) = sum gamma(h): bajo dependencia
## positiva v > gamma(0) y el IC "iid" con sd(X)/sqrt(n) es muy angosto

set.seed(3907290)
n <- 100
phi <- 0.7
mu <- 1
sigma <- 0.12
X <- mu + arima.sim(n, model = list(order = c(1,0,0), ar = phi), sd = sigma)

par(bty = "n", las = 1)
plot(X, col = "gray", xlab = "")
abline(h = mean(X), lty = 2)
## IC ingenuo (como si fuera iid): sd(X)/sqrt(n)
LI <- mean(X) - qnorm(0.975)*sd(X)/sqrt(n)
LS <- mean(X) + qnorm(0.975)*sd(X)/sqrt(n)
abline(h = c(LI, LS), col = "red", lty = 3)
## IC correcto para el AR(1): v = sigma^2/(1-phi)^2  (= 2 pi f(0))
v <- sigma^2/(1 - phi)^2
LI <- mean(X) - qnorm(0.975)*sqrt(v/n)
LS <- mean(X) + qnorm(0.975)*sqrt(v/n)
abline(h = c(LI, LS), col = "blue", lty = 3)
legend("topleft", legend = c("IC iid (malo)", "IC AR(1)"),
       col = c("red","blue"), lty = 3, bty = "n")

## En la practica phi y sigma se estiman; arima ademas reporta el
## error estandar de mu directamente (var.coef):
mod <- arima(X, order = c(1,0,0)) ## Adelanto de la Clase 15
phi   <- mod$coef[1]
sigma <- sqrt(mod$sigma2)
mu    <- mod$coef[2]
LI <- mu - qnorm(0.975)*sqrt((sigma^2/(1-phi)^2)/n)
LS <- mu + qnorm(0.975)*sqrt((sigma^2/(1-phi)^2)/n)
cbind(LI, LS)                                  ## formula v = 2 pi f(0)
LI <- mu - qnorm(0.975)*sqrt(diag(mod$var.coef))[2]
LS <- mu + qnorm(0.975)*sqrt(diag(mod$var.coef))[2]
cbind(LI, LS)                                  ## var.coef de arima: similar

##################################################################
## Paso 3: ACF muestral --> bandas de Bartlett                  ##
##################################################################

## sqrt(n) (rho.hat(h) - rho(h)) --> N(0, w(h)), con w(h) de Bartlett.
## Bajo RB, w = 1 para todo h: las bandas +-1.96/sqrt(n) de acf()

Bartlett <- function(ma = NULL, ar = NULL, m = 3000, lag.max = 10){
  rho <- ARMAacf(ma = ma, ar = ar, lag.max = lag.max + m)
  j <- 1:m
  w <- c()
  for(h in 1:lag.max){
    w[h] <- sum((rho[abs(h+j)+1] + rho[abs(h-j)+1] - 2*rho[h+1]*rho[j+1])^2)
  }
  w
}
Bartlett(ar = 0, lag.max = 5)   ## RB: w = 1 en todos los rezagos

X <- LSTS::malleco
n <- length(X)
par(bty = "n", las = 1)
acf(X, lag.max = 20, ylim = c(-1,1), main = "", lwd = 2, ci = 0)

## Ajuste AR(1) por emparejamiento de momentos: phi.hat = rho.hat(1)
phi <- acf(X, plot = FALSE)$acf[2,,1]
ACF <- ARMAacf(ar = phi, lag.max = 20)
lines(ACF ~ c(0:20), type = "p", pch = 20, col = "red")

## Banda de confianza de Bartlett alrededor de la ACF del AR(1)
w  <- Bartlett(ma = 0, ar = phi, lag.max = 20)
LI <- ACF - qnorm(0.975)*sqrt(c(0,w)/n)
LS <- ACF + qnorm(0.975)*sqrt(c(0,w)/n)
lines(LI ~ c(0:20), col = "red", lty = 2)
lines(LS ~ c(0:20), col = "red", lty = 2)
## la ACF muestral cae dentro de la banda: el AR(1) es compatible

##################################################################
## Paso 4: diagnostico --> los residuos del AR(1) por momentos  ##
##################################################################

## Fijamos phi.hat = rho.hat(1) (fixed) y miramos si lo que sobra es
## ruido blanco: ACF de residuos y Box-Ljung en bloque (Clase 12)
Z <- forecast::Arima(X - mean(X), order = c(1,0,0), fixed = phi,
                     include.mean = FALSE)$res
LSTS::ts.diag(Z, lag = 20)
## residuos compatibles con blancura: el AR(1) "a ojo" desde la ACF ya
## captura la dependencia del malleco. En la Clase 15 lo estimamos bien:
## Yule-Walker, maxima verosimilitud y Whittle

## Moraleja:
## - AIC/BIC resuelven el compromiso ajuste vs parsimonia; BIC es mas
##   exigente con cada parametro extra.
## - Con datos correlacionados, los IC "iid" para la media son
##   incorrectos: la varianza correcta es v = 2 pi f(0).
## - La ACF muestral es asintoticamente normal con varianza de Bartlett;
##   las bandas +-1.96/sqrt(n) valen SOLO bajo ruido blanco (testean
##   blancura, no son un IC para rho(h) de un proceso cualquiera).
