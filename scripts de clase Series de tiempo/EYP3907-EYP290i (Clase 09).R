#############
## Clase 9 ##
#############

## Capitulo 2: Procesos Lineales
## - Causalidad e invertibilidad (raices del polinomio AR y MA)
## - Memoria corta (ARMA: decaimiento geometrico de la ACF)
##   vs memoria larga (decaimiento hiperbolico)

## Requiere: LSTS, forecast, longmemo

#######################################################
## Paso 1: AR(1) para malleco --> rho(k) = phi^k     ##
#######################################################

## ARMA(p,q):  phi(B) (Y[t] - mu) = theta(B) e[t],  {e[t]} ~ RB(0, sigma^2)
## Un AR(1) es un ARMA(1,0):  Y[t] - mu = phi (Y[t-1] - mu) + e[t]

X <- LSTS::malleco
par(bty = "n", las = 1, mfrow = c(1,1))
plot(X)

## ACF muestral y estimacion de phi por metodo de momentos: rho(1) = phi
rho <- acf(X, lag.max = 10, ylim = c(-1,+1))$acf[-1,,1]
phi <- rho[1]
phi

## ACF teorica del AR(1): decaimiento GEOMETRICO phi^h (memoria corta)
h <- 0:10
lines(phi^h ~ h, type = "p", pch = 20, col = "red")

## Ajuste del AR(1) con el coeficiente fijado en el valor de momentos
mod <- forecast::Arima(X, order = c(1,0,0), fixed = c(phi),
                       include.mean = FALSE)
mod

## plot() de un objeto Arima muestra las RAICES INVERSAS del polinomio AR:
## proceso CAUSAL <=> raices de phi(z) fuera del circulo unitario
##                <=> raices inversas DENTRO del circulo
plot(mod)

##########################################################
## Paso 2: ¿Y si el polinomio AR tiene raices "malas"?  ##
##########################################################

## Ejemplo: phi(z) = 1 + 0.55 z - 1.6 z^2  (phi1 = -0.55, phi2 = 1.6)
phi1 <- -0.55
phi2 <- +1.6

## Raices de phi(z) y sus modulos: si alguna cae DENTRO del circulo
## unitario, el proceso NO es causal (no admite Wold hacia el pasado)
z <- polyroot(c(1, -phi1, -phi2))
z
Mod(z)          ## ambas raices con modulo < 1 --> NO causal
1/Mod(z)        ## raices inversas: quedan FUERA del circulo unitario

## forecast::Arima se puede "forzar" a dibujarlas:
mod <- forecast::Arima(X, order = c(2,0,0), fixed = c(phi1, phi2),
                       include.mean = FALSE, method = "ML")
plot(mod)       ## raices inversas fuera del circulo rojo

## (si Arima reclama "non-stationary AR part", es la misma moraleja:
##  el software exige causalidad para poder estimar/filtrar)

#################################################################
## Paso 3: el caso limite phi = 1 --> paseo aleatorio          ##
#################################################################

## Y[t+1] = Y[t] + e[t+1] es un "AR(1) con phi = 1": raiz unitaria,
## NO estacionario (la clase lo anota sobre la Figura 1)
set.seed(3907)
e <- rnorm(300)
Y <- cumsum(e)

par(bty = "n", las = 1, mfrow = c(1,2))
ts.plot(Y)
acf(Y, main = "")   ## ACF muestral decae lentisimo (casi lineal)

## Si igual ajustamos un AR(1), el phi estimado se pega a 1:
forecast::Arima(Y, order = c(1,0,0), include.mean = FALSE)

## Moraleja: ACF que "no baja" --> sospechar no estacionaridad
## (o memoria larga, Paso 4); un ARMA de memoria corta no la reproduce

####################################################
## Paso 4: memoria larga --> rio Nilo (NileMin)   ##
####################################################

## Memoria corta: sum |gamma(h)| < infinito (ARMA: gamma(h) ~ C phi^h)
## Memoria larga: sum |gamma(h)| = infinito (gamma(h) ~ C h^(2d-1))

library(longmemo)
data(NileMin)
X <- ts(NileMin, start = 622)

par(bty = "n", las = 1, mfrow = c(1,1))
plot(X, ylab = "", ylim = c(800, 1600), col = "gray")
abline(h = mean(X))

## ACF muestral: decaimiento MUY lento (hiperbolico),
## correlaciones significativas hasta rezagos altisimos
acf(X, lag.max = 80, ylim = c(-1,+1), main = "")