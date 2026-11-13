##############
## Clase 17 ##
##############

## Capitulo 6: Series de Tiempo No Lineales
## - ARCH(1) simulado: ruido blanco con cuadrados AR(1)
## - Deteccion en datos reales: retornos del Dow Jones

## Requiere: fGarch, quantmod, moments

########################################################
## Paso 1: ARCH(1) simulado                           ##
########################################################

## y[t] = sigma[t] e[t],  sigma[t]^2 = alpha0 + alpha1 y[t-1]^2

set.seed(3907290)
n <- 200
alpha0 <- 0.3
alpha1 <- 0.8
e <- rnorm(n)
y <- e[1]
sigma <- c()
for(t in 2:n){
  sigma[t] <- sqrt(alpha0 + alpha1*y[t-1]^2)
  y[t] <- e[t] * sigma[t]
}

par(mfrow = c(1,1), bty = "n", las = 1)
ts.plot(y, col = "gray40")
## rachas de volatilidad a la vista, con media constante 0

## La serie es RUIDO BLANCO: ACF sin estructura...
acf(y, lag.max = 10, ylim = c(-1,+1), main = "")

## ...pero los CUADRADOS tienen la firma de un AR(1) con parametro alpha1
par(mfrow = c(1,2), bty = "n", las = 1)
acf(y^2, lag.max = 10, ylim = c(-1,+1), main = "")
pacf(y^2, lag.max = 10, ylim = c(-1,+1), xlim = c(0,10), main = "")
## (resultado clave de las slides: y^2 - mu = alpha1 (y^2[t-1] - mu) + nu[t])

## Ajuste ARCH(1) = garch(1,0): recupera (alpha0, alpha1)
library(fGarch)
garchFit(formula = ~ arma(0,0) + garch(1,0), data = y)

## ¿Y si le ponemos un ARMA a la media? Nada es significativo de verdad
garchFit(formula = ~ arma(1,1) + garch(1,0), data = y)
## La aparente significancia del ARMA(1,1) se debe a que los polinomios
## tratan de compensarse (podrian anularse mutuamente)

## Colas pesadas GRATIS: con alpha1 = 0.3 (6to momento finito) la curtosis
## teorica es 3(1-a1^2)/(1-3a1^2) = 3.74 > 3, aun con innovaciones normales
a1 <- 0.3
yk <- garchSim(garchSpec(model = list(omega = 0.5, alpha = a1, beta = 0)),
               n = 100000)$garch
c(curtosis.sim = mean(yk^4)/mean(yk^2)^2,
  teorica = 3*(1-a1^2)/(1-3*a1^2))

##################################################################
## Paso 2: deteccion en datos reales --> Dow Jones              ##
##################################################################

library(quantmod)
getSymbols("^DJI", src = "yahoo", from = "2025-01-01", to = Sys.Date())

n <- dim(DJI)[1]
X <- ts(DJI$DJI.Adjusted)
par(mfrow = c(1,1), bty = "n")
plot(X)

## Rentabilidad diaria
Y <- c(X[2:n]/X[1:(n-1)] - 1)
par(mfrow = c(1,2), bty = "n", las = 1)
ts.plot(Y)
plot(density(na.omit(Y)), main = "")

mean(Y, na.rm = TRUE)
sd(Y, na.rm = TRUE)
moments::kurtosis(Y, na.rm = TRUE)   ## >> 3: colas pesadas

## La receta de deteccion (la P3 de la I1): Box-Ljung en Y y en Y^2
par(mfrow = c(1,2), bty = "n", las = 1)
acf(Y,   lag.max = 10, ylim = c(-1,+1), main = "", na.action = na.pass)
acf(Y^2, lag.max = 10, ylim = c(-1,+1), main = "", na.action = na.pass)
