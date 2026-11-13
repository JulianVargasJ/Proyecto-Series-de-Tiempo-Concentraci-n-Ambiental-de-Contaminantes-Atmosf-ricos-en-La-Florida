##############
## Clase 25 ##
##############

## Capitulo 9: Estacionalidad (cierre)
## - GARMA: un peak de larga memoria donde uno quiera (Gegenbauer),
##   ciclos persistentes que no se apagan ni son deterministas; las
##   manchas solares (sunspot.year, base R) ajustadas por Whittle
## - SARIMAX con el Imacec de la I3: la "X" mas simple es la tendencia
##   (drift); Box-Cox; drift vs d = 1 (sobrediferenciacion); back test 2026

## Requiere: forecast, LSTS, rio
## Datos:    CCNN2018_IMACEC_01_A.xlsx (Banco Central)

##################################################################
## Capitulo 9: GARMA                                            ##
##################################################################

## Viene de la Clase 24 (SARFIMA): mismos n y grilla de frecuencias
n <- 1200
lam <- seq(0.005, pi, length = 3000)

########################################################
## Paso 1: GARMA --> un peak de larga memoria donde   ##
##         uno quiera (Gegenbauer)                    ##
########################################################

## (1 - 2 eta B + B^2)^(-d) = sum_j C_j(d, eta) B^j, polinomios de Gegenbauer
## con la recursion C_0 = 1, C_1 = 2 eta d,
## C_j = 2 eta ((d-1)/j + 1) C_{j-1} - (2 (d-1)/j + 1) C_{j-2}
gegen <- function(d, eta, m){
  C <- numeric(m + 1); C[1] <- 1; C[2] <- 2*eta*d
  for(j in 2:m) C[j+1] <- 2*eta*((d - 1)/j + 1)*C[j] - (2*(d - 1)/j + 1)*C[j-1]
  C
}
d <- 0.3
eta0 <- cos(2*pi/12)                     ## peak en la frecuencia estacional
C <- gegen(d, eta0, 3000)
round(C[1:8], 4)                         ## oscilan y decaen lentamente

set.seed(1)
e <- rnorm(n + length(C))
y.g <- ts(stats::filter(e, C, sides = 1)[-(1:length(C))][1:n], frequency = 12)

par(mfrow = c(2,1))
plot(y.g, xlim = c(0, 20), ylab = "", main = "GARMA(d = 0.3, eta = cos(2 pi/12)): ciclo de 12 que no se apaga")
acf(c(y.g), lag.max = 120, main = "ACF: oscila con periodo 12 y decae como |h|^(2d - 1)")
## compare con un AR(2) con raices complejas (Cap. 3): tambien oscila,
## pero su ACF muere geometricamente

## Las tres ACF TEORICAS lado a lado (d = 0.3, periodo 12):
acf.ma <- function(psi, H)                ## ACF de un MA(infinito) truncado
  sapply(0:H, function(h) sum(psi[1:(length(psi) - h)]*psi[(1 + h):length(psi)]))/sum(psi^2)
H <- 120; h <- 0:H
rho.g  <- acf.ma(gegen(d, eta0, 20000), H)                 ## GARMA
r <- 0.9; rho.ar <- ARMAacf(ar = c(2*r*cos(2*pi/12), -r^2), lag.max = H)  ## AR(2), periodo 12
rho.sf <- numeric(H + 1); rho.sf[1] <- 1; k <- 1:(H %/% 12)          ## SARFIMA (1 - B^12)^(-d)
rho.sf[1 + 12*k] <- gamma(1 - d)/gamma(d)*gamma(k + d)/gamma(k + 1 - d)

par(mfrow = c(3,1), mar = c(3.5,4,2.5,1))
plot(h, rho.g, type = "h", ylim = c(-1,1), ylab = "", main = "GARMA: oscila con periodo 12, envolvente |h|^(2d - 1)")
abline(h = 0); lines(h[-1], 0.75*h[-1]^(2*d - 1), col = "red", lty = 2); lines(h[-1], -0.75*h[-1]^(2*d - 1), col = "red", lty = 2)
plot(h, rho.ar, type = "h", ylim = c(-1,1), ylab = "", main = "AR(2) raices complejas, periodo 12, modulo 0.9: envolvente r^h")
abline(h = 0); lines(h, r^h, col = "red", lty = 2); lines(h, -r^h, col = "red", lty = 2)
plot(h, rho.sf, type = "h", ylim = c(-1,1), ylab = "", main = "SARFIMA: solo en los multiplos de 12, envolvente |h|^(2d - 1)")
abline(h = 0)
round(rbind(h = c(12, 24, 36, 60, 120), GARMA = rho.g[1 + c(12, 24, 36, 60, 120)],
            AR2 = rho.ar[1 + c(12, 24, 36, 60, 120)], SARFIMA = rho.sf[1 + c(12, 24, 36, 60, 120)]), 3)
## en h = 12 GARMA y AR(2) son casi iguales (0.25 vs 0.28): el periodo no
## discrimina. En h = 36 difieren en un orden de magnitud (0.16 vs 0.02):
## lo que discrimina es la persistencia. El SARFIMA es un peine: nada entre
## multiplos de 12, y mas alto en ellos porque concentra ahi toda la varianza

## Densidad espectral: f(lambda) = sigma^2/(2 pi) |2 (cos lambda - eta)|^(-2 d)
f.g <- abs(2*(cos(lam) - eta0))^(-2*d)/(2*pi)
I <- LSTS::periodogram(c(y.g), plot = FALSE)
par(mfrow = c(1,1))
plot(I$lambda, I$periodogram, type = "l", col = "gray60", ylim = c(0, 4),
     xlab = expression(omega), ylab = "", main = "GARMA: periodograma de la simulacion y f teorica")
lines(lam, f.g, col = "blue", lwd = 2)
abline(v = acos(eta0), lty = 3, col = "red")
## eta = 1  --> acos(1) = 0: el ARFIMA del Cap. 2 (larga memoria en la tendencia)
## eta = cos(2 pi/12) --> UN peak estacional; el SARFIMA pone peaks en TODOS
## los armonicos. GARMA: el peak donde uno quiera

## Las tres densidades TEORICAS lado a lado, como con las ACF. Se dividen
## por gamma(0) (varianza 1) para que sean comparables, igual que rho(h):
phi1 <- 2*r*cos(2*pi/12); phi2 <- -r^2
g0.g  <- sum(gegen(d, eta0, 20000)^2)                           ## gamma(0) del GARMA
g0.ar <- (1 - phi2)/((1 + phi2)*((1 - phi2)^2 - phi1^2))        ## gamma(0) del AR(2)
g0.sf <- gamma(1 - 2*d)/gamma(1 - d)^2                          ## gamma(0) del SARFIMA
f.ar <- 1/Mod(1 - phi1*exp(-1i*lam) - phi2*exp(-2i*lam))^2/(2*pi)
f.sf <- Mod(1 - exp(-12i*lam))^(-2*d)/(2*pi)
par(mfrow = c(3,1), mar = c(3.5,4,2.5,1))
plot(lam, f.g/g0.g, type = "l", col = "blue", lwd = 2, ylim = c(0, 2), ylab = "",
     main = "GARMA: UN peak en acos(eta) = 0.52, angosto y diverge")
abline(v = acos(eta0), lty = 3, col = "red")
plot(lam, f.ar/g0.ar, type = "l", col = "blue", lwd = 2, ylim = c(0, 2), ylab = "",
     main = "AR(2) periodo 12, modulo 0.9: peak en la misma frecuencia, pero FINITO y ancho")
abline(v = acos(eta0), lty = 3, col = "red")
plot(lam, f.sf/g0.sf, type = "l", col = "blue", lwd = 2, ylim = c(0, 2), ylab = "", xlab = expression(omega),
     main = "SARFIMA: peaks en 0 y en TODOS los armonicos 2 pi j/12, todos divergen")
abline(v = 2*pi*(0:6)/12, lty = 3, col = "red")
## GARMA y AR(2) coinciden en DONDE esta el peak (el periodo); difieren en
## su forma: infinito y angosto (larga memoria) vs finito y ancho (memoria
## corta). El SARFIMA repite el peak del GARMA en los siete armonicos.
## Es el espejo de las ACF: onda persistente / onda que muere / peine

########################################################
## Paso 2: GARMA con datos --> las manchas solares    ##
########################################################

## sunspot.year (base R): numero de Wolf anual, 1700-1988. El ciclo solar
## de ~11 anios (el del proyecto de rayos cosmicos de Oulu)
Y <- sunspot.year
x <- sqrt(Y) - mean(sqrt(Y))             ## raiz: estabiliza la varianza
par(bty = "n", las = 1, mfrow = c(2,1), mar = c(4,4,2,1))
plot(Y, ylab = "", main = "Manchas solares, 1700-1988")
I <- LSTS::periodogram(c(x), plot = FALSE)
plot(I$lambda, I$periodogram, type = "l", col = "gray40", xlab = expression(omega), ylab = "",
     main = "Periodograma: un solo peak, lejos de 0")
lam0 <- I$lambda[which.max(I$periodogram)]
abline(v = lam0, lty = 3, col = "red")
c(lambda0 = lam0, periodo = 2*pi/lam0)   ## 11.1 anios

## Whittle (Cap. 5) con la densidad GARMA(1,d,0):
## f(lambda) = s2/(2 pi) |2 (cos lambda - eta)|^(-2d) / |1 - phi e^{-i lambda}|^2
f.garma <- function(lam, d, eta, phi)
  abs(2*(cos(lam) - eta))^(-2*d)/Mod(1 - phi*exp(-1i*lam))^2/(2*pi)
whittle <- function(th){
  d <- th[1]; eta <- th[2]; phi <- th[3]
  if(d <= 0 | d >= 0.5 | abs(eta) >= 1 | abs(phi) >= 1) return(1e10)
  f <- f.garma(I$lambda, d, eta, phi)
  s2 <- mean(I$periodogram/f)            ## sigma^2 concentrado
  sum(log(s2*f) + I$periodogram/(s2*f))
}
fit.g <- optim(c(0.3, cos(2*pi/11), 0.3), whittle)
d.hat <- fit.g$par[1]; eta.hat <- fit.g$par[2]; phi.hat <- fit.g$par[3]
c(d = d.hat, eta = eta.hat, periodo = 2*pi/acos(eta.hat), phi = phi.hat)
## d = 0.34, periodo 11.1 anios: larga memoria EN EL CICLO, no en la tendencia

## El competidor de memoria corta: AR(2) con raices complejas (Cap. 3)
fit.ar2 <- arima(x, order = c(2,0,0), include.mean = FALSE)
r <- polyroot(c(1, -coef(fit.ar2)))
c(periodo.AR2 = 2*pi/Arg(r[1]), modulo = Mod(r[1]))   ## tambien ~11 anios...
## ...pero su ACF muere geometricamente; la muestral no:
a <- acf(c(x), lag.max = 55, plot = FALSE)$acf
round(rbind(muestral = a[1 + 11*(1:5)],
            AR2      = ARMAacf(ar = coef(fit.ar2), lag.max = 55)[1 + 11*(1:5)]), 3)
## en el rezago 33 (tres ciclos) la ACF muestral es 0.24 y el AR(2) predice
## 0.003: el ciclo solar NO se olvida como un AR(2). Whittle lo confirma:
w.ar2 <- function(th){
  f <- 1/Mod(1 - th[1]*exp(-1i*I$lambda) - th[2]*exp(-2i*I$lambda))^2/(2*pi)
  s2 <- mean(I$periodogram/f); sum(log(s2*f) + I$periodogram/(s2*f))
}
c(whittle.GARMA = fit.g$value, whittle.AR2 = optim(coef(fit.ar2), w.ar2)$value)
## menor es mejor: GARMA(1,d,0) con 3 parametros le gana al AR(2) por lejos

par(mfrow = c(1,1))
plot(I$lambda, I$periodogram, type = "l", col = "gray70", xlab = expression(omega), ylab = "", xlim = c(0, 2))
s2 <- mean(I$periodogram/f.garma(I$lambda, d.hat, eta.hat, phi.hat))
lines(lam, s2*f.garma(lam, d.hat, eta.hat, phi.hat), col = "blue", lwd = 2)
abline(v = acos(eta.hat), lty = 3, col = "red")
## la densidad ajustada diverge en acos(eta) = 0.57: ahi vive el ciclo

## Moraleja (Cap. 9, GARMA):
## - Gegenbauer pone la larga memoria en acos(eta): ciclos persistentes;
##   eta = 1 (lambda_0 = 0) devuelve el ARFIMA. Un peak, donde uno quiera;
##   el SARFIMA pone el mismo peak en todos los armonicos.
## - Las manchas solares: ciclo de 11 anios con d = 0.34, que un AR(2)
##   del mismo periodo olvida demasiado rapido (rezago 33: 0.24 vs 0.003).
## - SARIMAX (regresion + SARIMA) cierra el capitulo, a continuacion; su
##   parte X tiene teoria propia: Cap. 10, la proxima clase.

##################################################################
## Capitulo 9: SARIMAX                                          ##
##################################################################

## El Imacec de la I3: serie original (indice 2018 = 100), ene-1996 a
## sep-2026. SARIMAX = regresion + SARIMA en una sola verosimilitud;
## aqui la parte "X" es la mas simple de todas: la tendencia t (drift)
Data <- rio::import("CCNN2018_IMACEC_01_A.xlsx", skip = 2)
head(Data)
Y <- ts(Data[,2], start = c(1996, 1), frequency = 12)
par(bty = "n", las = 1)
plot(Y, ylab = "", main = "Imacec, serie original")

########################################################
## Paso 3: la transformacion --> estabilizar la       ##
##         varianza antes de modelar                  ##
########################################################

## Coeficiente de variacion por anio: si la amplitud estacional crece con
## el nivel, el CV es estable en nivel y lo que hay que estabilizar es
## la escala (log o Box-Cox)
cov <- function(x) sd(x)/mean(x)
plot(aggregate(Y ~ trunc(time(Y)), FUN = cov), ylab = "", main = "CV por anio, escala original")
## ~0.04-0.06 todos los anios: la dispersion es proporcional al nivel

lambda <- forecast::BoxCox.lambda(Y)        ## Guerrero: lambda = -0.58
Y.bc <- forecast::BoxCox(Y, lambda = lambda)
head(forecast::InvBoxCox(Y.bc, lambda = lambda))  ## la vuelta: recupera Y
plot(Y.bc, ylab = "", main = "Box-Cox(Y), lambda = -0.58")
plot(aggregate(Y.bc ~ trunc(time(Y.bc)), FUN = cov), ylim = c(0, 0.01), ylab = "", main = "CV por anio, Box-Cox")
plot(aggregate(log(Y) ~ trunc(time(log(Y))), FUN = cov), ylab = "", main = "CV por anio, log")
## log (lambda = 0) ya deja el CV plano en ~0.01; Box-Cox lo aplana mas
## pero en una escala sin interpretacion. Con lambda = -0.58 seguimos
## para ver que forecast::Arima maneja la transformacion por nosotros:
## estima en la escala Box-Cox y devuelve fitted y pronosticos en la
## escala original (biasadj = FALSE: es la mediana, no la media)

########################################################
## Paso 4: SARIMAX con drift --> regresion sobre t    ##
##         mas dinamica ARMA                          ##
########################################################

## include.drift = TRUE agrega la columna xreg = t: es la regresion
## y_t = a + b t + eta_t del Cap. 1, ahora con eta_t ARMA
mod <- forecast::Arima(Y, lambda = lambda, include.drift = TRUE)   ## order = (0,0,0): eta_t ruido blanco
mod
plot(Y, ylab = "", main = "Tendencia lineal, errores ruido blanco")
lines(mod$fitted, col = "red")
## la recta no sigue la estacionalidad ni los ciclos: todo queda en el
## residuo (sigma2 = 5.3e-05 en escala Box-Cox)

mod <- forecast::Arima(Y, order = c(1,0,1), lambda = lambda, include.drift = TRUE)
mod
plot(Y, ylab = "", main = "Tendencia lineal, errores ARMA(1,1)")
lines(mod$fitted, col = "red")
## phi = 0.99: la dinamica se come la tendencia, el AR esta al borde de
## la raiz unitaria --> la serie no es estacionaria alrededor de una
## recta, hay que diferenciar. Y falta la estacionalidad: D = 1

## Entrenamiento hasta dic-2025; los 9 meses de 2026 quedan fuera para
## validar (el mismo back test del Cap. 9)
Y.tr <- window(Y, end = c(2025, 12))
Y.te <- window(Y, start = c(2026, 1))

mod1 <- forecast::Arima(Y.tr, order = c(1,0,1), seasonal = c(0,1,0), lambda = lambda, include.drift = TRUE)
mod1
## (1 - B^12) Y.bc = drift + ARMA(1,1): phi = 0.89, theta = 0.04 (casi
## AR(1)), drift = 0.0002/mes. Tras la diferencia estacional el drift
## es el crecimiento medio anual (12 x 0.0002 en escala Box-Cox)
plot(Y, ylab = "", main = "(1,0,1)(0,1,0) con drift")
lines(mod1$fitted, col = "red")
pre1 <- forecast::forecast(mod1, h = 90, level = 95)
plot(pre1); lines(Y, lwd = 1); lines(pre1$fitted, col = "red", lwd = 1)
plot(pre1, xlim = c(2024, 2027), ylim = c(100, 130)); lines(Y, lwd = 1); lines(pre1$fitted, col = "red", lwd = 1)
## el pronostico repite el patron estacional del ultimo anio desplazado
## por el drift; los 9 puntos de 2026 (negro) caen TODOS bajo la media
## pronosticada, pero dentro de la banda

########################################################
## Paso 5: d = 1 + D = 1 + ARMA(1,1) --> sin drift,   ##
##         y la senal de sobrediferenciacion          ##
########################################################

## Alternativa: en vez de drift, diferencia regular (d = 1). Con d = 1 y
## D = 1 forecast::Arima no admite drift (la constante se va con la
## segunda diferencia), por eso include.drift = FALSE
mod2 <- forecast::Arima(Y.tr, order = c(1,1,1), seasonal = c(0,1,0), lambda = lambda, include.drift = FALSE)
mod2
## theta = -0.99: raiz unitaria en la parte MA. Es la P1(f) de la I3:
## (1 - B^12) ya elimino la tendencia lineal, y (1 - B) encima la
## sobrediferencia; el MA(1) con theta -> -1 esta "deshaciendo" esa
## diferencia. sigma2 casi igual que mod1 (2.29e-06 vs 2.27e-06) y
## AICc peor (-3530 vs -3541): el drift era la forma correcta
plot(Y, ylab = "", main = "(1,1,1)(0,1,0): theta = -0.99, sobrediferenciado")
lines(mod2$fitted, col = "red")
pre2 <- forecast::forecast(mod2, h = 90, level = 95)
plot(pre2); lines(Y, lwd = 1); lines(pre2$fitted, col = "red", lwd = 1)
plot(pre2, xlim = c(2024, 2027), ylim = c(100, 130)); lines(Y, lwd = 1); lines(pre2$fitted, col = "red", lwd = 1)
## la banda de mod2 se abre mas rapido: con d = 1 la varianza del error
## de prediccion no converge (Cap. 8), con drift y d = 0 si

## El competidor natural: el airline de la I3, (0,1,1)(0,1,1)
mod3 <- forecast::Arima(Y.tr, order = c(0,1,1), seasonal = c(0,1,1), lambda = lambda)
c(AICc.drift = mod1$aicc, AICc.d1 = mod2$aicc, AICc.airline = mod3$aicc)
## el airline gana por ~95 unidades: el SMA(1) con Theta = -0.68 recoge
## la estacionalidad que (0,1,0) deja en el residuo
sapply(list(drift = mod1, d1 = mod2, airline = mod3),
       function(m) round(Box.test(residuals(m), lag = 12, type = "Ljung")$statistic, 1))
## Q(12) = 81.5, 81.7 y 31.3: ninguno pasa (12 g.l., critico 21.0); los
## pulsos de 2020 siguen ahi (LSTS::Box.Ljung.Test(residuals(mod3), lag = 24)
## para verlo rezago a rezago)

## Validacion fuera de muestra, ene-sep 2026
bt <- function(m){ p <- forecast::forecast(m, h = length(Y.te), level = 95)
  c(RMSE = sqrt(mean((Y.te - p$mean)^2)), MAPE = 100*mean(abs(Y.te - p$mean)/Y.te),
    cobertura = mean(Y.te >= p$lower & Y.te <= p$upper)) }
round(rbind(drift = bt(mod1), d1 = bt(mod2), airline = bt(mod3)), 2)
## MAPE 2.4-3.0% y cobertura 100% en los tres; el mejor en AICc (airline)
## no es el mejor fuera de muestra (d1), y la diferencia es chica: 2026
## fue peor que lo que cualquier modelo estimado hasta 2025 proyectaba
## (la contraccion de la noticia), y las tres bandas la contienen

## Moraleja (Cap. 9, SARIMAX):
## - La "X" mas simple es t: include.drift = TRUE es la regresion sobre
##   la tendencia con errores ARMA, estimada junto con la dinamica.
## - Con D = 1 el drift es el crecimiento medio anual; agregar d = 1
##   encima sobrediferencia y se delata en theta -> -1 (P1(f) de la I3).
## - lambda = Box-Cox se estima aparte (Guerrero); Arima transforma,
##   ajusta y destransforma solo. Con log se pierde poco y se gana
##   interpretacion.
## - Comparar modelos: AICc dentro de muestra y back test fuera; cuando
##   discrepan, mirar la magnitud de la diferencia antes que el orden.
