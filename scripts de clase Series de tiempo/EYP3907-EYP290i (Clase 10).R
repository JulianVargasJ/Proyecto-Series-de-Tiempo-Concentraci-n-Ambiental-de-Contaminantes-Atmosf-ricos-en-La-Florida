##############
## Clase 10 ##
##############

## Capitulo 2: Funcion de auto-covarianza y auto-correlacion ARMA
## - ACF exacta con ARMAacf(...): Compute Theoretical ACF for an ARMA Process
## - ARMA(1,1): rho(h) = phi^(h-1) rho(1) --> la TASA la gobierna phi,
##   el nivel/signo lo gobierna theta
## - Identificacion de un MA(q): la ACF se corta en el rezago q
## - Yule-Walker "a mano" para malleco

## Requiere: LSTS, forecast

##########################################################
## Paso 1: un caso curioso --> cancelacion de raices    ##
##########################################################

## ARMA(1,1):  (1 - phi B) y[t] = (1 + theta B) e[t]
## Con theta = -phi ambos polinomios comparten la raiz:
## (1 - 0.5 B) y[t] = (1 - 0.5 B) e[t]  -->  y[t] = e[t]  (RUIDO BLANCO)
phi   <- +0.5
theta <- -0.5
ARMAacf(ar = phi, ma = theta, lag.max = 10)   ## ACF nula en todo rezago

## Por eso la definicion de ARMA(p,q) exige que phi(z) y theta(z)
## NO compartan raices: si no, el modelo esta sobre-parametrizado

################################################################
## Paso 2: ACF exacta del ARMA(1,1) --> rho(h) = phi^(h-1) k  ##
################################################################

phi   <- +0.5
theta <- +0.5
ACF <- ARMAacf(ar = phi, ma = theta, lag.max = 10)

par(bty = "n", las = 1, mfrow = c(1,1))
plot(ACF ~ c(0:10), type = "h", xlab = "Lag")
abline(h = 0)

## ¿Es la ACF de un AR(1) puro? NO: phi^h queda por debajo
lines(phi^(0:10) ~ c(0:10), type = "p", pch = 20, col = "red")

## Formula de la clase:  rho(1) = (phi+theta)(1+theta*phi) / (1+2*phi*theta+theta^2)
k <- (phi+theta)*(1+theta*phi) / (1+2*phi*theta + theta^2)
k

## ¿Decae como k^h? Tampoco:
lines(k^c(0:10) ~ c(0:10), type = "p", col = "blue", pch = 20)

## La ACF exacta es  rho(h) = phi^(h-1) * rho(1)  para h >= 1:
## la TASA de decaimiento es phi (como el AR(1)), pero "arranca" en rho(1)
lines(c(1, k*phi^(0:9)) ~ c(0:10), type = "p", col = "darkgreen", pch = 1, cex = 2)
## los circulos verdes calzan exactamente con ARMAacf(...)

#################################################
## Paso 3: galeria de ACF exactas (Figura 2)   ##
#################################################

par(bty = "n", las = 1, mfrow = c(2,2))

## AR(1) con phi < 0: decaimiento geometrico ALTERNANTE
ACF <- ARMAacf(ar = -0.8, lag.max = 10)
plot(ACF ~ c(0:10), type = "h", ylim = c(-1,+1), xlab = "Lag",
     main = expression(AR(1)*":  "*phi == -0.8))
abline(h = 0)

## AR(1) con phi > 0: decaimiento geometrico positivo
ACF <- ARMAacf(ar = +0.8, lag.max = 10)
plot(ACF ~ c(0:10), type = "h", ylim = c(-1,+1), xlab = "Lag",
     main = expression(AR(1)*":  "*phi == +0.8))
abline(h = 0)

## ARMA(1,1), phi = 0.5: la tasa la pone phi ...
ACF <- ARMAacf(ar = +0.5, ma = +0.7, lag.max = 10)
plot(ACF ~ c(0:10), type = "h", ylim = c(-1,+1), xlab = "Lag",
     main = expression(list(phi == +0.5, theta == +0.7)))
abline(h = 0)

## ... y el SIGNO lo pone theta (Clave de la Figura 2)
ACF <- ARMAacf(ar = +0.5, ma = -0.7, lag.max = 10)
plot(ACF ~ c(0:10), type = "h", ylim = c(-1,+1), xlab = "Lag",
     main = expression(list(phi == +0.5, theta == -0.7)))
abline(h = 0)

## Otro ejemplo mixto: phi negativo domina la alternancia
par(bty = "n", las = 1, mfrow = c(1,1))
ACF <- ARMAacf(ar = -0.8, ma = 0.9, lag.max = 10)
plot(ACF ~ c(0:10), type = "h", ylim = c(-1,+1), xlab = "Lag")
abline(h = 0)

####################################################################
## Paso 4: identificacion de un MA(q) --> la ACF se CORTA en q    ##
####################################################################

## Regla practica: si la ACF muestral no es significativa despues del
## rezago q, un MA(q) es un candidato natural (Figura de la clase)
par(bty = "n", las = 1, mfrow = c(1,2))

ACF <- ARMAacf(ma = c(0.7, -0.6), lag.max = 10)
plot(ACF ~ c(0:10), type = "h", ylim = c(-1,+1), xlab = "Lag", main = "MA(2)")
abline(h = 0)

ACF <- ARMAacf(ma = c(0.3, -0.4, 0.7, 0.5), lag.max = 10)
plot(ACF ~ c(0:10), type = "h", ylim = c(-1,+1), xlab = "Lag", main = "MA(4)")
abline(h = 0)

## gamma(h) = 0 para |h| > q: memoria FINITA (a diferencia del AR)

########################################################
## Paso 5: malleco --> AR(1) vs ARMA(1,1)             ##
########################################################

X <- LSTS::malleco

par(bty = "n", las = 1, mfrow = c(1,1))
rho <- acf(X, ylim = c(-1,+1), lag.max = 10)$acf[-1,,1]

## (i) AR(1) por momentos: phi = rho(1) --> ACF teorica phi^h
phi <- rho[1]
lines(phi^(0:10) ~ c(0:10), type = "p", col = "red", pch = 20)

## (ii) ARMA(1,1) por Yule-Walker "a mano" (pizarra de la clase):
## gamma(2) = phi gamma(1)  -->  phi = rho(2)/rho(1)
phi.yw <- rho[2] / rho[1]
phi.yw
## y con rho(h) = phi^(h-1) rho(1) reconstruimos TODA la ACF teorica
## sin necesidad de conocer theta:
lines(c(1, rho[1]*phi.yw^(0:9)) ~ c(0:10), type = "p", col = "darkgreen", pch = 20)

## (iii) ARMA(1,1) por maxima verosimilitud (forecast::Arima)
aux <- forecast::Arima(X - mean(X), order = c(1,0,1), include.mean = FALSE)$coef
aux
ACF <- ARMAacf(ar = aux[1], ma = aux[2], lag.max = 10)
lines(ACF ~ c(0:10), type = "p", col = "blue", pch = 20)

## Moraleja:
## - Rojo (AR(1), phi = rho(1)) y verde/azul (ARMA(1,1)) son muy parecidos
##   aqui: malleco tiene una dependencia corta que un AR(1) ya captura bien.
## - Yule-Walker entrega estimadores por MOMENTOS (rapidos, sin verosimilitud);
##   ML es mas eficiente. Ambos viven de la misma estructura rho(h).
