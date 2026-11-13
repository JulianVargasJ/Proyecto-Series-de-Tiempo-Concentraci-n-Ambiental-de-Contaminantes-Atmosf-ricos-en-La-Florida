##############
## Clase 22 ##
##############

## Capitulo 9: Estacionalidad
## - SARMA teorico: ACF/PACF con estructura en multiplos de s
## - Densidad espectral con peaks en los armonicos estacionales
## - El modelo airline sobre co2: nsdiffs, ajuste, diagnostico
## - auto.arima estacional y prediccion
## - Back test del airline: pronosticar a medida que llega la data
##
## Requiere: forecast, LSTS

##################################################################
## Capitulo 9: Estacionalidad                                   ##
##################################################################

########################################################
## Paso 1: SARMA teorico --> ACF y PACF               ##
########################################################

## Ejemplo 1: SAR(1)_12 puro --> Phi = 0.7
ACF <- ARMAacf(ar = c(rep(0,11), 0.7), lag.max = 120)
PACF <- ARMAacf(ar = c(rep(0,11), 0.7), lag.max = 120, pacf = TRUE)
par(bty = "n", las = 1, mfrow = c(1,2))
plot(ACF ~ c(0:120), type = "h", ylim = c(-1,1), xlab = "Lag")
abline(h = 0); lines(ACF ~ c(0:120), type = "p", pch = 20, col = "red")
plot(PACF ~ c(1:120), type = "h", ylim = c(-1,1), xlab = "Lag")
abline(h = 0); lines(PACF ~ c(1:120), type = "p", pch = 20, col = "red")
## ACF: picos SOLO en 12, 24, 36... decayendo como 0.7^k;
## PACF: un solo palo en 12 --> la identificacion es "mirar multiplos de s"

## Ejemplo 2: SARMA(1,0)(1,0)[12] --> phi = 0.5, Phi = 0.7
## la forma multiplicativa (1 - 0.5B)(1 - 0.7B^12) expandida:
ar.exp <- c(0.5, rep(0,10), 0.7, -0.5*0.7)
ACF <- ARMAacf(ar = ar.exp, lag.max = 120)
PACF <- ARMAacf(ar = ar.exp, lag.max = 120, pacf = TRUE)
par(bty = "n", las = 1, mfrow = c(1,2))
plot(ACF ~ c(0:120), type = "h", ylim = c(-1,1), xlab = "Lag")
abline(h = 0); lines(ACF ~ c(0:120), type = "p", pch = 20, col = "red")
plot(PACF ~ c(1:120), type = "h", ylim = c(-1,1), xlim = c(0,120), xlab = "Lag")
abline(h = 0); lines(PACF ~ c(1:120), type = "p", pch = 20, col = "red")
## regular x estacional: decaimiento corto ALREDEDOR de cada multiplo de 12

########################################################
## Paso 2: densidad espectral --> peaks en 2 pi j/12  ##
########################################################

fdens <- function(lam, ar = numeric(), sar = numeric(), s = 12, s2 = 1){
  den <- sapply(lam, function(l)
    Mod(1 - sum(ar*exp(-1i*l*seq_along(ar))))^2 *
    Mod(1 - sum(sar*exp(-1i*l*s*seq_along(sar))))^2)
  s2/(2*pi)/den
}
lam <- seq(0.001, pi, length = 2000)
par(bty = "n", las = 1, mfrow = c(1,1))
plot(lam, fdens(lam, sar = 0.8), type = "l", col = "blue", lwd = 2,
     xlab = expression(omega), ylab = expression(f(omega)))
abline(v = 2*pi*(1:6)/12, lty = 3, col = "gray")
## peaks EXACTAMENTE en los armonicos estacionales 2 pi j/12 del Cap. 4:
## la teoria detras de los peaks del periodograma del co2

########################################################
## Paso 3: el modelo airline sobre co2                ##
########################################################

X <- co2                       ## Mauna Loa, mensual (base R)
par(bty = "n", las = 1, mfrow = c(1,1))
plot(X, ylab = "", xlab = "")

## ¿cuantas diferencias? Primero la estacional, luego la regular
forecast::nsdiffs(X)                 ## D = 1 (medida de fuerza estacional, STL)
forecast::ndiffs(diff(X, lag = 12))  ## d = 1, evaluada tras aplicar (1 - B^12)
## nsdiffs NO es un test de raiz unitaria: recomienda (1-B^12) cuando la
## estacionalidad es fuerte, determinista o estocastica. Y (1-B^12) elimina
## de paso cualquier tendencia lineal. Por eso la pregunta por la
## diferencia regular se hace DESPUES; si se hace antes, ndiffs responde
## a una tendencia que igual iba a desaparecer. auto.arima hace esto mismo.
## En el co2 ambos ordenes dan d = D = 1 (la tendencia es mas que lineal);
## el ejemplo simulado que sigue muestra cuando el orden si importa

## Ejemplo: tendencia lineal + estacionalidad deterministas + AR(1)
set.seed(3907)
n <- 240; tt <- 1:n
y <- ts(0.05*tt + 5*cos(2*pi*tt/12) + arima.sim(list(ar = 0.6), n),
        frequency = 12)
c(ndiffs.y = forecast::ndiffs(y), nsdiffs.y = forecast::nsdiffs(y),
  ndiffs.diff12 = forecast::ndiffs(diff(y, lag = 12)))
## ndiffs(y) = 1 y nsdiffs(y) = 1 --> (1-B)(1-B^12): sobrediferenciado
f1 <- forecast::Arima(y, order = c(1,1,1), seasonal = c(0,1,1))
## nsdiffs(y) = 1 y luego ndiffs(diff(y,12)) = 0 --> (1-B^12) + drift
f2 <- forecast::Arima(y, order = c(1,0,1), seasonal = c(0,1,1),
                      include.drift = TRUE)
round(rbind(f1 = c(coef(f1), drift = NA), f2 = coef(f2)), 3)
c(AICc.f1 = f1$aicc, AICc.f2 = f2$aicc)
## f1: theta = -1.000 (raiz unitaria en la parte MA, la huella de la
## sobrediferenciacion). f2: drift = 0.049 (verdadero 0.05) y menor AICc

## el periodograma como detector (las figuras de las I3):
aux <- LSTS::periodogram(diff(X), plot = FALSE)
par(bty = "n", las = 1)
plot(aux$periodogram ~ aux$lambda, type = "l", col = "gray40",
     xlab = expression(omega), ylab = "")
abline(v = 2*pi*(1:6)/12, lty = 3, col = "red")
## los peaks estacionales siguen ahi tras (1-B): falta (1-B^12)
round(2*pi/aux$lambda[rev(order(aux$periodogram))[1:2]])  ## periodos top

aux2 <- LSTS::periodogram(diff(diff(X), lag = 12), plot = FALSE)
plot(aux2$periodogram ~ aux2$lambda, type = "l", col = "gray40",
     xlab = expression(omega), ylab = "")
abline(v = 2*pi*(1:6)/12, lty = 3, col = "red")
## con (1-B)(1-B^12) los peaks desaparecen: listo para un SARMA corto

## El modelo airline: SARIMA(0,1,1)(0,1,1)[12]
fit <- forecast::Arima(X, order = c(0,1,1), seasonal = c(0,1,1))
summary(fit)
## dos parametros generan estructura en rezagos 1, 12 y 13 (parsimonia)

## diagnostico: residuos blancos, sin estructura en multiplos de 12
LSTS::Box.Ljung.Test(residuals(fit), lag = 24)
par(bty = "n", las = 1)
acf(c(residuals(fit)), lag.max = 36, main = "")

########################################################
## Paso 4: auto.arima estacional y prediccion         ##
########################################################

fit2 <- forecast::auto.arima(X)      ## recorre (p,d,q)(P,D,Q) con AICc
summary(fit2)
## comparar AICc con el airline: suelen quedar muy cerca

pre <- forecast::forecast(fit, h = 24, level = 0.95)
plot(pre, main = "", bty = "n", las = 1, xlim = c(1990, 2000))
## el pronostico hereda nivel, tendencia y patron estacional; la banda
## crece sin estabilizarse (dos raices unitarias: regular y estacional)

##################################################################
## Paso 5: back test --> pronosticar a medida que llega la data ##
##################################################################

## Igual que con el IPC (Cap. 8): nos paramos en cada mes k como si
## fuera "hoy", pronosticamos 12 meses con el airline y comparamos con
## lo que realmente ocurrio.

REESTIMAR   <- FALSE   ## TRUE: reestima el modelo en cada origen (honesto)
                       ## FALSE: usa los coeficientes del ajuste completo
GUARDAR_GIF <- FALSE   ## TRUE: guarda la animacion en backtest_co2.gif
H <- 12                ## horizonte de pronostico

ord <- forecast::arimaorder(fit)                 ## (p,d,q)(P,D,Q)[s] del Paso 3
origenes <- which(time(X) >= 1990 & seq_along(X) < length(X))

pres <- vector("list", length(origenes))         ## pronosticos por origen
err  <- matrix(NA, length(origenes), H)          ## error a h pasos
pct  <- matrix(NA, length(origenes), H)          ## error porcentual a h pasos
cub  <- matrix(NA, length(origenes), H)          ## ¿cayo dentro de la banda?

for(i in seq_along(origenes)){
  k <- origenes[i]
  X.k <- window(X, end = time(X)[k])
  fit.k <- if(REESTIMAR) forecast::Arima(X.k, order = ord[1:3],
                                         seasonal = ord[4:6]) else
                         forecast::Arima(X.k, model = fit)
  pres[[i]] <- forecast::forecast(fit.k, h = H, level = 95)
  m <- min(H, length(X) - k)                     ## horizontes ya observados
  obs <- X[(k+1):(k+m)]
  err[i, 1:m] <- obs - pres[[i]]$mean[1:m]
  pct[i, 1:m] <- 100*err[i, 1:m]/obs
  cub[i, 1:m] <- obs >= pres[[i]]$lower[1:m] & obs <= pres[[i]]$upper[1:m]
}

## Un cuadro de la animacion: pronostico desde el origen i vs lo real
cuadro <- function(i){
  k <- origenes[i]
  plot(pres[[i]], xlim = c(1988, 1999), ylim = c(345, 370),
       main = sprintf("Origen: %d-%02d", floor(time(X)[k] + 1e-6), cycle(X)[k]),
       bty = "n", las = 1, font.main = 1, ylab = "co2")
  lines(X, col = "gray50")                       ## lo que de verdad ocurrio
  abline(v = time(X)[k], lty = 3)                ## "hoy"
}

## Animacion en el panel de graficos (RStudio)
par(mfrow = c(1,1))
for(i in seq_along(origenes)){
  cuadro(i)
  Sys.sleep(0.15)
}

## Resumen: el error crece con h (Cap. 7) y la cobertura deberia
## rondar el 95% si las bandas estan bien calibradas
res <- data.frame(h = 1:H,
                  RMSE = sqrt(colMeans(err^2, na.rm = TRUE)),
                  MAE  = colMeans(abs(err), na.rm = TRUE),
                  MAPE = colMeans(abs(pct), na.rm = TRUE),   ## en %
                  cobertura = colMeans(cub, na.rm = TRUE))
round(res, 3)
## A diferencia del IPC, aqui las bandas si cumplen: cobertura ~93-95%
## a corto plazo y ~100% desde h = 10 (incluso algo conservadoras).
## El RMSE a 12 meses es ~0.65 ppm sobre un nivel de ~360 (MAPE ~0.16%):
## el patron estacional es tan estable que el airline lo pronostica casi
## sin error. El MAPE es la metrica comparable entre series de distinta
## escala (la que pide el proyecto junto a RMSE y MAE)

par(mfrow = c(1,2), bty = "n", las = 1)
plot(res$h, res$RMSE, type = "b", pch = 20, xlab = "h", ylab = "RMSE")
plot(res$h, res$cobertura, type = "b", pch = 20, ylim = c(0, 1),
     xlab = "h", ylab = "cobertura banda 95%")
abline(h = 0.95, col = "red", lty = 2)
par(mfrow = c(1,1))

## Opcional: guardar la animacion como GIF (requiere magick y pdftools)
if(GUARDAR_GIF){
  pdf("backtest_co2.pdf", 10, 5.6)
  for(i in seq_along(origenes)) cuadro(i)
  dev.off()
  img <- magick::image_read_pdf("backtest_co2.pdf", density = 90)
  magick::image_write(magick::image_animate(img, fps = 5), "backtest_co2.gif")
}

## Moraleja:
## - La firma estacional: ACF/PACF en multiplos de s; f con peaks en
##   2 pi j/s. La identificacion es la de siempre, mirando esos rezagos.
## - (1-B^s) es a la estacionalidad lo que (1-B) a la tendencia, con el
##   mismo riesgo de sobrediferenciar; nsdiffs ayuda.
## - El airline (0,1,1)(0,1,1)[12]: dos parametros, rezagos 1, 12 y 13.
## - El periodograma es el mejor detector de que diferencia falta.
## - El back test del airline sobre co2 cumple la cobertura del 95%:
##   con estacionalidad estable, el SARIMA pronostica bien a 12 meses.
