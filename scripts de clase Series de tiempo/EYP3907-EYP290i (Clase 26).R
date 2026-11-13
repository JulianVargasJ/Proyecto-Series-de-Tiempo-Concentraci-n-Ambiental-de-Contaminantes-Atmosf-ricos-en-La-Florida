##############
## Clase 26 ##
##############

## Capitulo 9 / 10: del SARIMA al SARIMAX, paso a paso con el Imacec
## - auto.arima() como punto de partida, no de llegada: que cambia segun
##   como se le entregue la serie (sin ts, con Box-Cox, con ts)
## - Depurar el modelo a mano: fijar en cero lo no significativo (fixed)
##   y agregar rezagos MA aislados donde el diagnostico lo pide
## - Diagnostico en cada paso: coeftest, raices, ts.diag (Box-Ljung
##   rezago a rezago)
## - Validacion: pronostico de 2026 contra lo observado
##
## Requiere: forecast, rio, lubridate, lmtest, LSTS
## Datos (carpeta "Data Clase 26"): CCNN2018_IMACEC_01_A.xlsx

####################
## Modelo SARIMAX ##
####################

setwd("Data Clase 26")
dir()

###################
## Modelo IMACEC ##
###################

## Imacec (indice 2018 = 100), serie original y no minero, ene-1996 a
## sep-2026. skip = 2: el archivo del Banco Central trae titulo + linea
## en blanco antes de los encabezados
IMACEC <- rio::import("CCNN2018_IMACEC_01_A.xlsx", skip = 2)
colnames(IMACEC) <- c("periodo","imacec","imacec_nm")
IMACEC$periodo <- as.Date(IMACEC$periodo)
head(IMACEC)
IMACEC$year <- lubridate::year(IMACEC$periodo)
IMACEC$month <- lubridate::month(IMACEC$periodo)
head(IMACEC)
tail(IMACEC)

## Opción 1: SARIMA
## Entrenamiento hasta dic-2025 (360 meses); los 9 meses de 2026 quedan
## fuera para validar el pronostico al final
Y <- ts(IMACEC$imacec, start = c(IMACEC$year[1],IMACEC$month[1]), frequency = 12)
Y.train <- window(Y, end = c(2025,12))
Y.test <- window(Y, start = c(2026,1))

########################################################
## Paso 1: auto.arima() --> la misma serie, tres      ##
##         respuestas segun como se la entreguemos    ##
########################################################

mod <- forecast::auto.arima(c(Y.train)) ## Sin atributos ts()
mod
## c() borra la frecuencia: auto.arima no sabe que es mensual y no puede
## proponer parte estacional. Entrega ARIMA(1,1,2) con drift y toda la
## estacionalidad queda en el residuo (sigma2 = 15.6)

mod <- forecast::auto.arima(Y.train, lambda = "auto", allowdrift = T, d = 0)
mod
## Con Box-Cox (lambda = -0.32, Guerrero) y forzando d = 0: la tendencia
## la toma el drift y la estacionalidad D = 1. Resultado
## ARIMA(3,0,2)(0,1,1)[12] con drift. Ojo: sigma2 y AICc estan en la
## escala transformada, NO se comparan con los otros dos

mod <- forecast::auto.arima(Y.train) ## con atributos ts()
mod
lmtest::coeftest(mod) ## ar1 no significativo
## En escala original: ARIMA(3,0,1)(0,1,1)[12] con drift, sigma2 = 2.25,
## AICc = 1284.5. d = 0 + D = 1 + drift: la forma correcta de la Clase 25
## (el drift = 0.204 puntos de indice por mes, ~2.4 por anio, es el
## crecimiento medio tras la diferencia estacional).
## auto.arima minimiza AICc; no mira significancia: ar1 = 0.088 con
## z = 1.08 sobra. Es un punto de partida, no el modelo final

########################################################
## Paso 2: recuperar los ordenes y reajustar con      ##
##         forecast::Arima() --> control total        ##
########################################################

## mod$arma guarda los ordenes en ESTE orden: p, q, P, Q, s, d, D
mod$arma
p <- mod$arma[1]
q <- mod$arma[2]
P <- mod$arma[3]
Q <- mod$arma[4]
s <- mod$arma[5]
d <- mod$arma[6]
D <- mod$arma[7]
drift <- FALSE
if(sum(names(coef(mod)) == "drift") == 1){drift = TRUE}

## fixed: un valor por coeficiente, en el orden en que Arima los reporta:
## ar1..arp, ma1..maq, sar, sma, drift (y despues xreg si hubiera).
## NA = estimar, 0 = fijo en cero. Aqui: ar1 = 0 y los otros 5 libres
mod <- forecast::Arima(Y.train, order = c(p,d,q), seasonal = c(P,D,Q), include.drift = drift, fixed = c(0,rep(NA,5)))
mod
lmtest::coeftest(mod)
plot(mod)
## Sin ar1 todos los coeficientes son significativos y el AICc baja a
## 1283.7. plot(mod): raices inversas AR y MA dentro del circulo unitario
## (estacionario e invertible)

tsdiag(mod)
LSTS::ts.diag(c(mod$res), lag = 24)
## tsdiag (base) llega solo al rezago 10; ts.diag muestra residuos, ACF y
## el valor-p de Box-Ljung rezago a rezago hasta 24. Aqui la ACF residual
## sale de la banda en los rezagos 6, 9, 13, 19 y 24 y Box-Ljung rechaza
## desde el rezago 6: los residuos NO son ruido blanco, falta dinamica

########################################################
## Paso 3: agregar rezagos MA aislados --> uno a la   ##
##         vez, mirando el diagnostico                ##
########################################################

## Estrategia: el primer rezago que sale de la banda es el 6 --> subir q
## a 6 y dejar ma2..ma5 fijos en cero (un MA(6) "con hoyos"). Reajustar,
## mirar significancia y diagnostico, y repetir con el siguiente rezago

q <- 6
fixed <- c(0,NA,NA, NA,0,0,0,0,NA, NA, NA)   ## ar1=0, ar2, ar3 | ma1, ma2..ma5=0, ma6 | sma1 | drift
mod <- forecast::Arima(Y.train, order = c(p,d,q), seasonal = c(P,D,Q), include.drift = drift, fixed = fixed)
mod
lmtest::coeftest(mod)
plot(mod)
LSTS::ts.diag(c(mod$res), lag = 24)
## ma6 = 0.070 (z = 2.6), AICc = 1280.0. El rezago 6 queda resuelto;
## siguen fuera el 9 y el 24 y Box-Ljung rechaza desde el rezago 10

q <- 9
fixed <- c(0,NA,NA, NA,0,0,0,0,NA,0,0,NA, NA, NA)   ## se agrega ma9
mod <- forecast::Arima(Y.train, order = c(p,d,q), seasonal = c(P,D,Q), include.drift = drift, fixed = fixed)
mod
lmtest::coeftest(mod)
plot(mod)
LSTS::ts.diag(c(mod$res), lag = 24)
## ma9 = 0.054 con valor-p 0.10 (al limite), AICc = 1279.2. Box-Ljung ya
## solo rechaza en los rezagos 20 y 24; en la ACF queda el 13

q <- 13
fixed <- c(0,NA,NA, NA,0,0,0,0,NA,0,0,NA,0,0,0,NA, NA, NA)   ## se prueba ma13
mod <- forecast::Arima(Y.train, order = c(p,d,q), seasonal = c(P,D,Q), include.drift = drift, fixed = fixed)
mod
lmtest::coeftest(mod)
plot(mod)
LSTS::ts.diag(c(mod$res), lag = 24)
## ma13 = -0.046 NO significativo (valor-p 0.14) y el AICc no mejora
## (1279.1): se descarta. Aparece el rezago 11 en la ACF

q <- 11
fixed <- c(0,NA,NA, NA,0,0,0,0,NA,0,0,NA,0,NA, NA, NA)   ## ma13 fuera, entra ma11
mod <- forecast::Arima(Y.train, order = c(p,d,q), seasonal = c(P,D,Q), include.drift = drift, fixed = fixed)
mod
lmtest::coeftest(mod)
plot(mod)
LSTS::ts.diag(c(mod$res), lag = 24)
## ma11 = -0.098 (z = -2.7) y de paso ma9 se vuelve claramente
## significativo (0.140, z = 3.0). AICc = 1274.1, diez unidades mejor que
## el auto.arima. Box-Ljung no rechaza en NINGUN rezago hasta 24
## (valor-p 0.52 en el 12 y 0.06 en el 24): residuos ruido blanco

## Dos intentos mas, para confirmar que no falta nada:
q <- 11
Q <- 2
fixed <- c(0,NA,NA, NA,0,0,0,0,NA,0,0,NA,0,NA, NA,NA, NA)   ## un segundo SMA
mod <- forecast::Arima(Y.train, order = c(p,d,q), seasonal = c(P,D,Q), include.drift = drift, fixed = fixed)
mod
lmtest::coeftest(mod)
plot(mod)
LSTS::ts.diag(c(mod$res), lag = 24)
## sma2 = -0.027 no significativo (valor-p 0.67), AICc sube a 1276.1: no

q <- 20
Q <- 1
fixed <- c(0,NA,NA, NA,0,0,0,0,NA,0,0,NA,0,NA,0,0,0,0,0,0,0,0,NA, NA, NA)   ## ma20
mod <- forecast::Arima(Y.train, order = c(p,d,q), seasonal = c(P,D,Q), include.drift = drift, fixed = fixed)
mod
lmtest::coeftest(mod)
plot(mod)
LSTS::ts.diag(c(mod$res), lag = 24)
## ma20 = 0.012 no significativo (valor-p 0.67), AICc 1276.1: tampoco

########################################################
## Paso 4: el modelo final y su pronostico contra     ##
##         lo observado en 2026                       ##
########################################################

## ARIMA(3,0,11)(0,1,1)[12] con drift, con ar1 = 0 y solo ma1, ma6, ma9
## y ma11: 8 coeficientes estimados
q <- 11
Q <- 1
fixed <- c(0,NA,NA, NA,0,0,0,0,NA,0,0,NA,0,NA, NA, NA)
mod <- forecast::Arima(Y.train, order = c(p,d,q), seasonal = c(P,D,Q), include.drift = drift, fixed = fixed)
mod
lmtest::coeftest(mod)
plot(mod)
LSTS::ts.diag(c(mod$res), lag = 12)

pre <- forecast::forecast(mod, h = 12)
plot(pre, xlim = c(2024, 2027), ylim = c(100,140))
lines(Y)
## Y completa (negro) sobre el pronostico: los 9 meses de 2026 quedan
## TODOS bajo la media pronosticada (la contraccion de la economia), y
## 8 de 9 dentro de la banda del 95% (enero queda apenas fuera)
n.te <- length(Y.test)
round(cbind(obs = Y.test, pred = pre$mean[1:n.te], lo95 = pre$lower[1:n.te, 2], hi95 = pre$upper[1:n.te, 2]), 1)
c(MAPE = 100*mean(abs(Y.test - pre$mean[1:n.te])/Y.test),
  cobertura = mean(Y.test >= pre$lower[1:n.te, 2] & Y.test <= pre$upper[1:n.te, 2]))
## MAPE = 3.06% y cobertura 89%

## Cuanto se gano fuera de muestra respecto del punto de partida?
mod0 <- forecast::auto.arima(Y.train)
pre0 <- forecast::forecast(mod0, h = n.te)
c(MAPE.auto = 100*mean(abs(Y.test - pre0$mean)/Y.test), AICc.auto = mod0$aicc, AICc.final = mod$aicc)
## MAPE 3.15% vs 3.06%: diez unidades de AICc y residuos limpios dentro
## de muestra compran una decima fuera de muestra. (El airline con
## Box-Cox de la Clase 25 daba 2.85%.) El error de 2026 es de NIVEL: la
## economia crecio menos que el drift historico, y eso no lo arregla
## ningun rezago MA

## Moraleja:
## - auto.arima() depende de lo que recibe (frecuencia, lambda, d) y
##   elige por AICc sin mirar significancia: es el punto de partida.
## - fixed permite sacar lo que sobra y poner MA aislados donde el
##   diagnostico lo pide; el orden del vector es ar, ma, sar, sma, drift.
## - Un cambio a la vez: reajustar, coeftest, ts.diag. Lo que no es
##   significativo ni baja el AICc se devuelve (ma13, sma2, ma20).
## - Cuidado con perseguir rezagos: ma6, ma9 y ma11 limpian la ACF, pero
##   no tienen lectura economica y el pronostico de 2026 casi no mejora.
##   Lo que falta no es mas ARMA, es la "X" (pulsos 2020, calendario, PIB).