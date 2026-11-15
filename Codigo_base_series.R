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

install.packages("psych")
##############################################################################
## SEGMENTO 1: Lectura, Agregación Mensual e Imputación en Escala Log       ##
## - Clase 04: Atributos ts() y agregación temporal con aggregate()         ##
## - Clase 20 (Paso 1 y 3): Interpolación KalmanSmooth() y soporte positivo ##
## - Clase 26: Estructura de fechas con lubridate::year() y month()         ##
##############################################################################

## 1.1 Lectura de datos y formato de calendario (Clase 26)


MP25 <- read.csv2(
  "/Users/pilar/Proyecto-Series-de-Tiempo-Concentraci-n-Ambiental-de-Contaminantes-Atmosf-ricos-en-La-Florida/datos_020101_260823.csv",
  fileEncoding = "Latin1")[, c(1, 3, 4, 5)]


MP25 <- read.csv2("C:/Users/conjv/Downloads/datos_020101_260823.csv",
                  dec = ",", encoding = "latin1")[, c(1, 3, 4, 5)]



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

## 1.3 Separación Entrenamiento / Validación (10 % más reciente = validación)
## Se hace antes de imputar porque el suavizador de Kalman usa pasado y futuro, lo que filtraria

n_total <- length(X.na)
n_val   <- round(0.10 * n_total)          ## ~30 meses (2,5 años)
n_train <- n_total - n_val

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

## 2.1 Evaluación de heterocedasticidad: Guerrero (~0.55) vs LogLik (0.20) vs Log (0)
lambda_guerrero <- forecast::BoxCox.lambda(X_train, method = "guerrero") ## ~ 0.55 
lambda_loglik   <- forecast::BoxCox.lambda(X_train, method = "loglik")   ## = 0.20
c(Guerrero = lambda_guerrero, LogLik = lambda_loglik)

cov.fn <- function(x) sd(x)/mean(x)
X.bc   <- forecast::BoxCox(X_train, lambda = lambda_loglik)
Y      <- log(X_train) ## lambda = 0: menor CV anual (~0.12) e interpretación porcentual (Clase 25)

par(bty = "n", las = 1, mfrow = c(1, 3))
plot(aggregate(X_train ~ trunc(time(X_train)), FUN = cov.fn), type = "o", pch = 20,
     ylim = c(0, 0.6), ylab = "CV", xlab = "Año", main = "CV escala original (lambda = 1)")
plot(aggregate(X.bc ~ trunc(time(X.bc)), FUN = cov.fn), type = "o", pch = 20,
     ylim = c(0, 0.6), ylab = "CV", xlab = "Año", main = "CV Box-Cox (LogLik = 0.20)")
plot(aggregate(Y ~ trunc(time(Y)), FUN = cov.fn), type = "o", pch = 20,
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
c(Frecuencia = freq_peak, Periodo_meses = 2*pi / freq_peak) ## omega = 0.5197 -> 12.09 meses (grilla 2*pi*j/266)


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
tseries::adf.test(d12.Y)                    ## ADF = -5.63, p-valor < 0.01 -> estacionaria

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

## Detección de outliers (|Z| > 3): identifica feb-2023 (incendios forestales)
idx.out <- which(abs(res.std) > 3)
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

