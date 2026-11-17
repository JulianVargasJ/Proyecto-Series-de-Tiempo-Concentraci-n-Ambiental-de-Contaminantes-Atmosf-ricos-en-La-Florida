##############################################################################
## Proyecto EYP3907 — Concentración de MP2,5 en La Florida (2002-2026)      ##
## Código Hito 2: Limpieza, Exploración, Estacionariedad, Modelo y Diag.    ##
##############################################################################

library(lubridate)
library(forecast)
library(LSTS)
library(tseries)
library(lmtest)
library(moments)
library(psych)

##############################################################################
## SEGMENTO 1: Lectura, Agregación Mensual e Imputación en Escala Log       ##
## - Clase 04: Atributos ts() y agregación temporal con aggregate()         ##
## - Clase 20 (Paso 1 y 3): Interpolación KalmanSmooth() y soporte positivo ##
## - Clase 26: Estructura de fechas con lubridate::year() y month()         ##
##############################################################################

## 1.1 Lectura de datos y formato de calendario (Clase 26)


# MP25 <- read.csv2(
#  "/Users/pilar/Proyecto-Series-de-Tiempo-Concentraci-n-Ambiental-de-Contaminantes-Atmosf-ricos-en-La-Florida/datos_020101_260823.csv",
#  fileEncoding = "Latin1")[, c(1, 3, 4, 5)]

# MP25 <- read.csv2("C:/Users/conjv/Downloads/datos_020101_260823.csv",
#                  dec = ",", encoding = "latin1")[, c(1, 3, 4, 5)]


MP25 <- read.csv2( file.choose(), dec = ",",fileEncoding = "Latin1" )[, c(1, 3, 4, 5)]

colnames(MP25) <- c("fecha_num", "validados", "preliminares", "no_validados")
psych::headTail(MP25)

MP25$periodo <- as.Date(sprintf("%06d", MP25$fecha_num), format = "%y%m%d")
MP25$year    <- lubridate::year(MP25$periodo)
MP25$month   <- lubridate::month(MP25$periodo)
psych::headTail(MP25)

## Serie auxiliar unificada para contrastar qué midió el sensor en meses sin validar
MP25$sensor_bruto <- ifelse(!is.na(MP25$validados), MP25$validados,
                            ifelse(!is.na(MP25$preliminares), MP25$preliminares,
                                   MP25$no_validados))
psych::headTail(MP25)

## 1.2 Agregación mensual preservando el orden temporal (Clase 04)
M.val   <- aggregate(validados ~ month + year, data = MP25,
                     FUN = mean, na.rm = TRUE, na.action = na.pass)
M.bruto <- aggregate(sensor_bruto ~ month + year, data = MP25,
                     FUN = mean, na.rm = TRUE, na.action = na.pass)
M.val$validados[is.nan(M.val$validados)] <- NA

X.na    <- ts(M.val$validados, start = c(2002, 1), frequency = 12)
X.bruto <- ts(M.bruto$sensor_bruto, start = c(2002, 1), frequency = 12)

sum(is.na(X.na))           ## 14 meses faltantes (ene-sep 2013 y ene-may 2017)
time(X.na)[is.na(X.na)]

## 1.3 Separación Entrenamiento / Validación (20 % más reciente = validación)
## Se hace antes de imputar porque el suavizador de Kalman usa pasado y futuro, lo que filtraria

n_total <- length(X.na)
n_val   <- round(0.20 * n_total) ## 59 meses (~20%: oct-2021 a ago-2026)
n_train <- n_total - n_val ## 237 meses (~80%: ene-2002 a sep-2021)

X.na_train <- ts(X.na[1:n_train], start = start(X.na), frequency = 12)
X_val      <- ts(X.na[(n_train + 1):n_total],
                 start = time(X.na)[n_train + 1], frequency = 12)

c(n_total = n_total, n_train = n_train, n_val = n_val)
c(fin_train = time(X.na_train)[n_train], ini_val = time(X_val)[1])
sum(is.na(X_val))   ## debe ser 0: los faltantes están todos en el entrenamiento

## 1.4 Imputación con Suavizador de Kalman en escala log (Clase 20, Paso 1 y 3)
# solo sobre entrenamiento
logX.na   <- log(X.na_train)
media.log <- mean(logX.na, na.rm = TRUE)
W.na      <- logX.na - media.log

fit.na  <- stats::arima(W.na, order = c(1, 0, 0),
                        seasonal = list(order = c(0, 1, 1), period = 12))
ks      <- stats::KalmanSmooth(W.na, fit.na$model)
imp.log <- drop(ks$smooth %*% fit.na$model$Z) + media.log

## X_train es la serie de entrenamiento limpia que se usara el resto del analisis
X_train <- X.na_train
idx.na  <- which(is.na(X.na_train))
X_train[idx.na] <- exp(imp.log[idx.na])

## 1.5 Gráfico de la serie limpia e imputaciones (Clase 20, Paso 1)
par(bty = "n", las = 1, mfrow = c(1, 1))
plot(ts.union(X_train, X_val), plot.type = "single", col = c("gray40", "darkgreen"),
     lwd = 1.5, xlab = "Año", ylab = expression(mu*g/m^3),
     main = "MP2,5 mensual en La Florida (2002-2026): entrenamiento y validación")
points(time(X_train)[idx.na], X.bruto[idx.na],  pch = 1,  col = "blue", cex = 1.3)
points(time(X_train)[idx.na], X_train[idx.na],  pch = 19, col = "red",  cex = 1.0)
abline(v = time(X_val)[1], lty = 2)
legend("topright",
       legend = c("Entrenamiento (X_train)", "Validación (X_val)",
                  "Imputación (Kalman Smoother)", "Preliminar / no validado"),
       col = c("gray40", "darkgreen", "red", "blue"),
       lty = c(1, 1, NA, NA), pch = c(NA, NA, 19, 1), bty = "n")


##############################################################################
## SEGMENTO 2: Análisis Exploratorio, Box-Cox y Análisis Espectral          ##
## - Clase 04 y 25: Coeficiente de Variación (CV) y Transformación Box-Cox  ##
## - Clase 13 (Paso 3 y 5): Periodograma crudo y suavizado con splines      ##
##############################################################################

## 2.1 Evaluación de heterocedasticidad: Guerrero (~0.49) vs LogLik (0.20) vs Log (0)
lambda_guerrero <- forecast::BoxCox.lambda(X_train, method = "guerrero") ## ~ 0.55 
lambda_loglik   <- forecast::BoxCox.lambda(X_train, method = "loglik")   ## = 0.20
c(Guerrero = lambda_guerrero, LogLik = lambda_loglik)


cov.fn <- function(x) sd(x)/mean(x)
X_cv <- window(X_train, end = c(2020, 12))
X.bc <- forecast::BoxCox(X_cv, lambda = lambda_loglik)
Y_cv <- log(X_cv)
Y <- log(X_train) ## Serie de entrenamiento completa (237 meses) para modelar

par(bty = "n", las = 1, mfrow = c(1, 3))
plot(aggregate(X_cv ~ trunc(time(X_cv)), FUN = cov.fn), type = "o", pch = 20,
     ylim = c(0, 0.6), ylab = "CV", xlab = "Año", main = "CV escala original (lambda = 1)")
plot(aggregate(X.bc ~ trunc(time(X.bc)), FUN = cov.fn), type = "o", pch = 20,
     ylim = c(0, 0.6), ylab = "CV", xlab = "Año", main = "CV Box-Cox (LogLik = 0.20)")
plot(aggregate(Y_cv ~ trunc(time(Y_cv)), FUN = cov.fn), type = "o", pch = 20,
     ylim = c(0, 0.6), ylab = "CV", xlab = "Año", main = "CV escala log (lambda = 0)")

## 2.2 Firma estacional: ACF, PACF y Periodograma suavizado (Clase 13 Paso 5 y Clase 22)
par(bty = "n", las = 1, mfrow = c(1, 3))
acf(c(Y),  lag.max = 48, main = "ACF de log(MP2,5)")
pacf(c(Y), lag.max = 48, main = "PACF de log(MP2,5)")

P0 <- LSTS::periodogram(c(Y), plot = FALSE)
P1 <- LSTS::smooth.periodogram(c(Y), spar = 0.5, plot = FALSE)
plot(P0$lambda, P0$periodogram, type = "h", col = "gray60",
     xlab = expression(omega), ylab = expression(I(omega)), main = "Periodograma de log(MP2,5)")
lines(P1$lambda, P1$smooth.periodogram, col = "blue", lwd = 2)
abline(v = 2*pi*(1:6)/12, lty = 3, col = "red")
legend("topright", legend = c("Periodograma", "Spline suavizado", "Armónicos s=12"),
       col = c("gray60", "blue", "red"), lty = c(1, 1, 3), lwd = c(1, 2, 1), bty = "n")

freq_peak <- P0$lambda[which.max(P0$periodogram)]
c(Frecuencia = freq_peak, Periodo_meses = 2*pi / freq_peak) ## omega = 0.5302 -> 11.85 meses (grilla 2*pi*j/237, j = 20)


##############################################################################
## SEGMENTO 3: Estacionariedad y Orden de Diferenciación                    ##
## - Clase 21 (Paso 2): Test de raíz unitaria ADF y ndiffs (ADF, KPSS, PP)  ##
## - Clase 22 (Paso 3): Primero diferencia estacional (nsdiffs), luego reg. ##
##############################################################################

forecast::nsdiffs(Y)                        ## D = 1 -> aplicamos (1 - B^12)
d12.Y <- diff(Y, lag = 12)

forecast::ndiffs(d12.Y, test = "adf")       ## d = 0
forecast::ndiffs(d12.Y, test = "kpss")      ## d = 0
forecast::ndiffs(d12.Y, test = "pp")        ## d = 0
tseries::adf.test(d12.Y)                    ## ADF = -4.76, p-valor < 0.01 -> estacionaria

par(bty = "n", las = 1, mfrow = c(1, 3))
plot(d12.Y, ylab = "", xlab = "Año", col = "gray40", main = expression((1 - B^12)*log(X[t])))
abline(h = mean(d12.Y), col = "red", lty = 2)
acf(c(d12.Y), lag.max = 48, main = "ACF tras dif. estacional")

aux12 <- LSTS::periodogram(c(d12.Y), plot = FALSE)
plot(aux12$periodogram ~ aux12$lambda, type = "l", col = "gray40",
     xlab = expression(omega), ylab = "", main = "Periodograma tras (1 - B^12)")
abline(v = 2*pi*(1:6)/12, lty = 3, col = "red")


##############################################################################
## SEGMENTO 4: Estimación, Depuración (Clase 26) y Diagnóstico de Residuos  ##
## - Clase 14 (Paso 1) y 15 (Paso 4): AIC/BIC e inferencia asintótica MLE   ##
## - Clase 25 (Paso 4-5) y 26 (Paso 2-3): Drift, fixed y rezagos aislados   ##
## - Clase 03 (Paso 3) y 17-18: Normalidad, outliers y efecto ARCH en res^2 ##
##############################################################################

## 4.1 Ajuste de modelos competidores en escala log (lambda = 0)
mod.airline  <- forecast::Arima(X_train, order = c(0,1,1), seasonal = c(0,1,1), lambda = 0)
mod.sobredif <- forecast::Arima(X_train, order = c(1,1,1), seasonal = c(0,1,1), lambda = 0)
mod.101      <- forecast::Arima(X_train, order = c(1,0,1), seasonal = c(0,1,1),
                                include.drift = TRUE, lambda = 0)

## Depuración 1 (Clase 26, Paso 2): sacamos ma1 por no ser significativo (p = 0.43)
mod.100      <- forecast::Arima(X_train, order = c(1,0,0), seasonal = c(0,1,1),
                                include.drift = TRUE, lambda = 0)

## Depuración 2 (Clase 26, Paso 3): en ts.diag de mod.100, el rezago 15 sale de la
## banda y baja el valor-p de Box-Ljung al límite del 5%. Agregamos ma15 aislado con fixed:
## orden de coeficientes en fixed: ar1 | ma1..ma14=0, ma15 | sma1 | drift
fixed.15 <- c(NA, rep(0, 14), NA, NA, NA)
mod.ma15 <- forecast::Arima(X_train, order = c(1,0,15), seasonal = c(0,1,1),
                            include.drift = TRUE, lambda = 0, fixed = fixed.15)

lmtest::coeftest(mod.100)  ## Modelo parsimonioso (BIC óptimo)
lmtest::coeftest(mod.ma15) ## Modelo depurado con ma15 significativo (AICc óptimo y blancura total)

## 4.2 Tabla comparativa: Criterios de Información y Box-Ljung (en lag 12, 15 y 24)
resumen_mod <- function(m) {
  r <- residuals(m)
  c(AIC    = m$aic,
    AICc   = m$aicc,
    BIC    = m$bic,
    p_BL12 = Box.test(r, lag = 12, type = "Ljung")$p.value,
    p_BL15 = Box.test(r, lag = 15, type = "Ljung")$p.value,
    p_BL24 = Box.test(r, lag = 24, type = "Ljung")$p.value)
}

tabla.modelos <- rbind(
  Airline_011_011   = resumen_mod(mod.airline),
  Sobredif_111_011  = resumen_mod(mod.sobredif),
  SARIMA_101_drift  = resumen_mod(mod.101),
  SARIMA_100_drift  = resumen_mod(mod.100),
  Depurado_MA15     = resumen_mod(mod.ma15)
)
print(round(tabla.modelos, 3))

## 4.3 Diagnóstico de residuos e identificación de valores atípicos
## Puedes evaluar mod.final <- mod.100 (más parsimonioso en BIC) o mod.ma15 (mejor AICc)
mod.final <- mod.100
res       <- c(residuals(mod.final))
res.std   <- res / sqrt(mod.final$sigma2)

## Detección de outliers en entrenamiento (|Z| > 2.5): identifica feb-2014 y mar-2016
idx.out <- which(abs(res.std) > 2.5) ## 3
data.frame(Fecha = sprintf("%d-%02d", M.val$year[idx.out], M.val$month[idx.out]),
           MP25_ugm3 = round(X_train[idx.out], 2),
           Z_residuo = round(res.std[idx.out], 2))

## (a) Causalidad e invertibilidad (Clase 09 y 26)
par(bty = "n", las = 1, mfrow = c(1, 1))
plot(mod.final)

## (b) Comparación de Box-Ljung rezago a rezago: mod.100 vs mod.ma15 (Clase 26, Paso 3)
LSTS::ts.diag(c(residuals(mod.100)),  lag = 24)
LSTS::ts.diag(c(residuals(mod.ma15)), lag = 24)

## (c) Normalidad, curtosis y receta de detección ARCH en res^2 (Clase 03, 17 y 18)
par(bty = "n", las = 1, mfrow = c(2, 2))
hist(res, breaks = 25, freq = FALSE, col = "gray90", border = "white",
     main = "Distribución de residuos vs Normal", xlab = "Residuos")
lines(density(res), lwd = 2)
curve(dnorm(x, mean = mean(res), sd = sd(res)), add = TRUE, col = "red", lwd = 2)

qqnorm(res, pch = 20, col = "gray40", main = "QQ-plot Normal de residuos")
qqline(res, col = "red", lwd = 2)

acf(res,   lag.max = 24, ylim = c(-0.5, 1), lwd = 2, main = expression("ACF de "*hat(e)[t]))
acf(res^2, lag.max = 24, ylim = c(-0.5, 1), lwd = 2, main = expression("ACF de "*hat(e)[t]^2*" (Efecto ARCH)"))

cat(sprintf("Asimetría: %.3f (Normal = 0)\n", moments::skewness(res)))
cat(sprintf("Curtosis:  %.3f (Normal = 3)\n", moments::kurtosis(res)))
c(p_Shapiro_Normalidad = shapiro.test(res)$p.value,
  p_BoxLjung_res_24    = Box.test(res,   lag = 24, type = "Ljung")$p.value,
  p_BoxLjung_ARCH_12   = Box.test(res^2, lag = 12, type = "Ljung")$p.value,
  p_BoxLjung_ARCH_24   = Box.test(res^2, lag = 24, type = "Ljung")$p.value)




##############################################################################
## SEGMENTO 5: Validación Fuera de Muestra, Bandas de Predicción y Back-test##
## - Clase 20 (Paso 2-3): Bandas de predicción asimétricas con lambda = 0   ##
## - Clase 25 (Paso 4-5) y Clase 26 (Paso 4): Pronóstico sobre X_val y bt() ##
## - Clase 21 (Paso 4) y Clase 22 (Paso 5): Back-test móvil a H = 12 pasos  ##
##############################################################################

## Variable con la serie completa (entrenamiento limpio + validación) para graficar
X_completa <- ts(c(X_train, X_val), start = start(X_train), frequency = 12)
h_val      <- length(X_val) ## 59 meses de validación (oct-2021 a ago-2026)

## 5.1 Pronóstico sobre el conjunto de validación con bandas al 95% (Clases 25 y 26)
## Con lambda = 0, forecast() devuelve pronóstico y bandas asimétricas en ug/m3 (Clase 20)
pre.100      <- forecast::forecast(mod.100,      h = h_val, level = 95)
pre.ma15     <- forecast::forecast(mod.ma15,     h = h_val, level = 95)
pre.sobredif <- forecast::forecast(mod.sobredif, h = h_val, level = 95)

## (a) Gráfico general + Zoom a la ventana de validación (Estilo Clase 25 Paso 4 y Clase 26 Paso 4)
par(bty = "n", las = 1, mfrow = c(1, 2))

## Panel izquierdo: serie completa, ajuste in-sample (rojo) y pronóstico con banda 95%
plot(pre.ma15, main = "Pronóstico SARIMA depurado (2002-2026)",
     xlab = "Año", ylab = expression(mu*g/m^3), shadecols = "lightblue")
lines(X_completa, col = "gray30", lwd = 1)
lines(pre.ma15$fitted, col = "red", lwd = 1)
abline(v = time(X_val)[1], lty = 2)

## Panel derecho: Zoom al tramo de validación (2020-2026) para verificar captura de X_val
plot(pre.ma15, xlim = c(2020, 2026.8), ylim = c(5, 65),
     main = "Zoom Validación (oct-2021 a ago-2026)",
     xlab = "Año", ylab = expression(mu*g/m^3), shadecols = "lightblue")
lines(X_train, col = "gray40", lwd = 1.2)
lines(X_val,   col = "darkgreen", lwd = 1.8)
lines( as.numeric(time(X_val)), as.numeric(pre.ma15$mean), col = "blue", lwd = 1.8)
abline(v = time(X_val)[1], lty = 2)
legend("topleft",
       legend = c("Observado validación (X_val)", "Pronóstico (media)", "Banda IC 95%"),
       col = c("darkgreen", "blue", "lightblue"), lty = 1, lwd = c(1.8, 1.8, 8), bty = "n")

## (b) Contraste teórico de bandas (Clase 25, Paso 5): d = 0 con drift vs d = 1 (sobrediferenciado)
## Con d = 0 + drift la varianza del error converge; con d = 1 la banda se abre más rápido
par(bty = "n", las = 1, mfrow = c(1, 1))
plot(X_val, col = "darkgreen", lwd = 2, ylim = c(5, 75),
     xlab = "Año", ylab = expression(mu*g/m^3),
     main = "Bandas de predicción 95% en validación: d = 0 con drift vs d = 1 (sobredif.)")
lines(pre.ma15$mean,  col = "blue", lwd = 2)
lines(pre.ma15$lower, col = "blue", lty = 2, lwd = 1.5)
lines(pre.ma15$upper, col = "blue", lty = 2, lwd = 1.5)
lines(pre.sobredif$lower, col = "red", lty = 3, lwd = 1.5)
lines(pre.sobredif$upper, col = "red", lty = 3, lwd = 1.5)
legend("topleft",
       legend = c("X_val observado", "IC 95% Depurado (d=0 + drift)", "IC 95% Sobredif. (d=1)"),
       col = c("darkgreen", "blue", "red"), lty = c(1, 2, 3), lwd = 2, bty = "n")

## (c) Tabla de métricas en validación (RMSE, MAE, MAPE y Cobertura 95% - Función bt de Clase 25)
bt <- function(m) {
  p <- forecast::forecast(m, h = h_val, level = 95)
  c(RMSE      = sqrt(mean((X_val - p$mean)^2)),
    MAE       = mean(abs(X_val - p$mean)),
    MAPE      = 100 * mean(abs(X_val - p$mean) / X_val),
    cobertura = 100 * mean(X_val >= p$lower & X_val <= p$upper))
}

tabla.validacion <- rbind(
  Airline_011_011  = bt(mod.airline),
  Sobredif_111_011 = bt(mod.sobredif),
  SARIMA_101_drift = bt(mod.101),
  SARIMA_100_drift = bt(mod.100),
  Depurado_MA15    = bt(mod.ma15)
)
print(round(tabla.validacion, 2))
## Cobertura de 94.92% (56 de 59 meses dentro de la banda del 95%): calibración exacta.

## ¿Qué meses de validación se salen de la banda del 95%? (Destaca feb-2023 por incendios)
fuera.banda <- which(X_val < pre.ma15$lower | X_val > pre.ma15$upper)
data.frame(
  Fecha     = sprintf("%d-%02d", M.val$year[n_train + fuera.banda], M.val$month[n_train + fuera.banda]),
  Observado = round(X_val[fuera.banda], 2),
  Predicho  = round(pre.ma15$mean[fuera.banda], 2),
  LI_95     = round(pre.ma15$lower[fuera.banda], 2),
  LS_95     = round(pre.ma15$upper[fuera.banda], 2)
)

##############################################################################
## 5.2 Back-test móvil: pronosticar a 12 meses a medida que llega la data   ##
##     (Ref: Clase 21 Paso 4 y Clase 22 Paso 5)                             ##
##############################################################################

H <- 12 ## Horizonte de 12 meses (1 ciclo anual completo)
origenes <- which(seq_along(X_completa) >= n_train & seq_along(X_completa) < length(X_completa))

pres <- vector("list", length(origenes))
err <- matrix(NA, length(origenes), H)
pct <- matrix(NA, length(origenes), H)
cub <- matrix(NA, length(origenes), H)

for (i in seq_along(origenes)) {
  k <- origenes[i]
  X.k <- window(X_completa, end = time(X_completa)[k])
  ## Evaluamos el modelo preliminar seleccionado sin reestimar (como en Clase 22 Paso 5)
  fit.k <- forecast::Arima(X.k, model = mod.ma15)
  pres[[i]] <- forecast::forecast(fit.k, h = H, level = 95)
  
  m <- min(H, length(X_completa) - k)
  obs <- X_completa[(k + 1):(k + m)]
  err[i, 1:m] <- obs - pres[[i]]$mean[1:m]
  pct[i, 1:m] <- 100 * err[i, 1:m] / obs
  cub[i, 1:m] <- obs >= pres[[i]]$lower[1:m] & obs <= pres[[i]]$upper[1:m]
}

## Función para graficar cada cuadro del avance temporal (Clase 22, Paso 5)
cuadro <- function(i) {
  k <- origenes[i]
  plot(pres[[i]], xlim = c(2019, 2026.8), ylim = c(5, 65),
       main = sprintf("Back-test móvil — Origen: %d-%02d",
                      floor(time(X_completa)[k] + 1e-6), cycle(X_completa)[k]),
       bty = "n", las = 1, font.main = 1, xlab = "Año", ylab = expression(mu*g/m^3))
  lines(X_completa, col = "gray40", lwd = 1.5) ## Lo que realmente ocurrió
  abline(v = time(X_completa)[k], lty = 3)     ## "Hoy"
}

## Animación en RStudio (puedes detenerla con Esc o cambiar Sys.sleep)
par(mfrow = c(1, 1))
for (i in seq_along(origenes)) {
  cuadro(i)
  Sys.sleep(0.10)
}

## Resumen del Back-test por horizonte h = 1...12 (Clase 22, Paso 5)
res.bt <- data.frame(
  h         = 1:H,
  RMSE      = sqrt(colMeans(err^2, na.rm = TRUE)),
  MAE       = colMeans(abs(err), na.rm = TRUE),
  MAPE      = colMeans(abs(pct), na.rm = TRUE),
  cobertura = colMeans(cub, na.rm = TRUE)
)
print(round(res.bt, 3))

## Gráfico del error (RMSE) y cobertura de la banda del 95% según el horizonte h
par(mfrow = c(1, 2), bty = "n", las = 1)
plot(res.bt$h, res.bt$RMSE, type = "b", pch = 20, col = "blue",
     xlab = "Horizonte h (meses)", ylab = "RMSE", main = "Error de pronóstico vs horizonte")
plot(res.bt$h, res.bt$cobertura, type = "b", pch = 20, col = "darkgreen", ylim = c(0, 1),
     xlab = "Horizonte h (meses)", ylab = "Cobertura banda 95%", main = "Calibración del IC 95%")
abline(h = 0.95, col = "red", lty = 2)
par(mfrow = c(1, 1))