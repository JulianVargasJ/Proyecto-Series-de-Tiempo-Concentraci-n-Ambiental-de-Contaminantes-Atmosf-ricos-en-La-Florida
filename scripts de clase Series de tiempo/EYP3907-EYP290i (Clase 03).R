##############################################################################
## Script Clase 03 — Series financieras en R                                ##
## Rentabilidad diaria de Banco Santander Chile (BSANTANDER.SN)             ##
##############################################################################
## Contenidos:
##   1. Importación de datos (dos vías: quantmod y tidyquant)
##   2. Precio de cierre ajustado y retorno diario
##   3. ¿Son normales los retornos? (histograma, densidad, QQ-plot)
##   4. ¿Son independientes? (dispersión rezagada, elipses, ACF de X y X^2)
##   5. Ajuste de distribuciones (normal, logística, t de Student)
##   6. Dependencia serial y estacionaridad (serie malleco)
##############################################################################

## ---- 0. Paquetes -----------------------------------------------------------
## Se instalan automáticamente los que falten (útil para los alumnos).
paquetes <- c("quantmod", "tidyquant", "car", "corrplot", "moments",
              "fitdistrplus", "lmtest", "MASS", "LSTS")
faltantes <- setdiff(paquetes, rownames(installed.packages()))
if (length(faltantes) > 0) install.packages(faltantes, repos = "https://cloud.r-project.org")

## OJO: tidyquant carga dplyr, que enmascara stats::lag() y stats::filter().
## Por eso más abajo usamos siempre stats::lag() de forma explícita.
library(quantmod)
library(tidyquant)

op <- par(no.readonly = TRUE)   # guardamos los parámetros gráficos originales

## Parámetros generales del ejemplo
ticker <- "BSANTANDER.SN"
desde  <- as.Date("2010-01-01")
hasta  <- Sys.Date()

## ---- 1. Importación de datos ------------------------------------------------
## Vía 1: quantmod::getSymbols() --> objeto xts (serie de tiempo indexada)
getSymbols(ticker, src = "yahoo", from = desde, to = hasta)
head(BSANTANDER.SN)   # nota: el objeto se llama igual que el ticker
tail(BSANTANDER.SN)

## Vía 2: tidyquant::tq_get() --> data.frame "largo" (formato tidy)
## [, -1] elimina la columna symbol; renombramos para no pisar el xts anterior.
santander <- as.data.frame(tq_get(ticker, from = desde, to = hasta))[, -1]
head(santander)
str(santander)        # date, open, high, low, close, volume, adjusted

## ---- 2. Precio de cierre ajustado y retorno diario ---------------------------
## 2.1 Precio de cierre ajustado (corrige dividendos y splits)
par(bty = "n", las = 1)
plot(adjusted ~ date, data = santander, type = "l",
     xaxt = "n", ylim = range(santander$adjusted, na.rm = TRUE),
     xlab = "", ylab = "Precio de cierre ajustado (CLP)",
     main = paste("Banco Santander Chile:", format(desde), "a", format(hasta)))
axis.Date(1, at = pretty(santander$date, n = 8), format = "%Y")
abline(h = axTicks(2), lty = 2, col = "gray")
abline(v = pretty(santander$date, n = 8), lty = 2, col = "gray")

## 2.2 Retorno diario simple: R_t = P_t / P_{t-1} - 1
n <- nrow(santander)
santander$return <- c(NA, santander$adjusted[-1] / santander$adjusted[-n] - 1)
## Alternativa (retorno logarítmico, muy similar para retornos pequeños):
## santander$return <- c(NA, diff(log(santander$adjusted)))

X <- na.omit(santander$return)

par(bty = "n", las = 1, mfrow = c(1, 2))
plot(return ~ date, data = santander, type = "l", col = "gray",
     xaxt = "n", xlab = "", ylab = "Retorno diario",
     main = "Serie de retornos")
axis.Date(1, at = pretty(santander$date, n = 6), format = "%Y")
abline(h = mean(X), lwd = 2)   # media (visualmente ~0)

## ---- 3. ¿Son normales los retornos? ------------------------------------------
## 3.1 Histograma + densidad empírica vs densidad normal
## Los breaks se construyen desde los datos (evita errores si hay retornos
## fuera de un rango fijado a mano).
paso   <- 0.01
lim    <- ceiling(max(abs(X)) / paso) * paso
breaks <- seq(-lim, +lim, by = paso)

hist(X, breaks = breaks, border = "white", freq = FALSE,
     xlim = c(-0.1, +0.1), main = "Distribución de retornos", xlab = "")
lines(density(X, bw = 0.007, n = 2000), lwd = 2)   # bw fijo con fin comparativo
curve(dnorm(x, mean = mean(X), sd = sd(X)),
      from = -1, to = +1, n = 1000, add = TRUE, col = "red", lwd = 2)
legend("topleft", c("Densidad empírica", "Normal"), bty = "n",
       col = c("black", "red"), lwd = 2)

## 3.2 Estadísticos descriptivos
summary(X)
cat(sprintf("Desv. estándar: %.5f\n", sd(X)))
cat(sprintf("Asimetría:      %+.3f   (normal = 0)\n", moments::skewness(X)))
cat(sprintf("Curtosis:       %.3f   (normal = 3) --> colas pesadas\n",
            moments::kurtosis(X)))

## 3.3 QQ-plot normal: si los datos fueran normales, caerían sobre la recta.
## Las desviaciones en los extremos delatan las colas pesadas.
par(bty = "n", las = 1, mfrow = c(1, 1), font.main = 1)
qqnorm(X, col = "gray", pch = 20, xlim = c(-3, 3),
       main = "QQ-plot normal de los retornos")
qqline(X, lwd = 2)

## ---- 4. ¿Son independientes? --------------------------------------------------
## 4.1 Dispersión de X_t vs X_{t+1} con elipses de concentración
X_ts  <- ts(X)
datos <- na.omit(cbind(X_ts, stats::lag(X_ts, -1)))   # stats::lag explícito
colnames(datos) <- c("X_t", "X_t1")

plot(as.data.frame(datos), pch = 20, col = "gray",
     xlab = expression(X[t]), ylab = expression(X[t + 1]),
     main = "Retorno de hoy vs retorno de mañana")
abline(h = axTicks(2), lty = 2, col = "gray")
abline(v = axTicks(1), lty = 2, col = "gray")

## Elipses de concentración a distintos niveles (si hubiera correlación,
## estarían inclinadas; casi circulares ==> correlación ~ 0)
for (nivel in c(0.995, 0.95, 0.50)) {
  car::dataEllipse(as.numeric(datos[, 1]), as.numeric(datos[, 2]),
                   levels = nivel, add = TRUE, plot.points = FALSE,
                   col = "black", lty = 1, lwd = 2, center.pch = 0)
}

## Matriz de correlación visual
corrplot::corrplot(cor(as.data.frame(datos)), method = "ellipse")

## 4.2 ACF: la "firma" de los retornos financieros
## En un ruido iid, ni Z_t ni Z_t^2 muestran autocorrelación.
set.seed(3907290)
Z <- rnorm(2000)
par(mfrow = c(2, 2), bty = "n")
acf(Z,   main = expression(Z[t]),   lwd = 2, ylim = c(-0.5, +1))
acf(Z^2, main = expression(Z[t]^2), lwd = 2, ylim = c(-0.5, +1))

## En Santander: X_t casi incorrelacionado, pero X_t^2 SÍ está correlacionado
## ==> los retornos no son independientes (volatilidad agrupada, efecto ARCH).
m   <- min(length(X), 2000)                # ventana reciente (robusto si n < 2000)
X_m <- tail(X, m)
acf(X_m,   main = expression(X[t]),   lwd = 2, ylim = c(-0.5, +1))
acf(X_m^2, main = expression(X[t]^2), lwd = 2, ylim = c(-0.5, +1))

## ---- 5. Ajuste de distribuciones a los retornos --------------------------------
## MLE con fitdistrplus (normal, logística) y MASS (t de Student ubicación-escala)
ajuste_norm  <- fitdistrplus::fitdist(c(X), distr = "norm",  method = "mle")
ajuste_logis <- fitdistrplus::fitdist(c(X), distr = "logis", method = "mle")
ajuste_t     <- suppressWarnings(MASS::fitdistr(c(X), "t"))  # (m, s, gl)

par1 <- ajuste_norm$estimate
par2 <- ajuste_logis$estimate
par3 <- ajuste_t$estimate
dt_ls <- function(x, m, s, df) dt((x - m) / s, df = df) / s   # densidad t escalada

par(mfrow = c(1, 1), bty = "n", las = 1)
plot(density(X), lwd = 2, main = "Densidad empírica vs distribuciones ajustadas",
     xlab = "Retorno diario")
curve(dnorm(x, mean = par1[1], sd = par1[2]),
      from = -0.2, to = +0.2, n = 1000, add = TRUE, col = "red")
curve(dlogis(x, location = par2[1], scale = par2[2]),
      from = -0.2, to = +0.2, n = 1000, add = TRUE, col = "blue")
curve(dt_ls(x, m = par3[1], s = par3[2], df = par3[3]),
      from = -0.2, to = +0.2, n = 1000, add = TRUE, col = "darkgreen")
legend("topleft", bty = "n", lwd = c(2, 1, 1, 1),
       col = c("black", "red", "blue", "darkgreen"),
       legend = c("Empírica", "Normal", "Logística",
                  sprintf("t de Student (gl = %.1f)", par3[3])))

## Comparación por AIC (menor = mejor): la t de colas pesadas domina
aic <- c(Normal       = ajuste_norm$aic,
         Logística    = ajuste_logis$aic,
         "t de Student" = 2 * 3 - 2 * ajuste_t$loglik)
print(sort(round(aic, 1)))

##############################################################################
## 6. Dependencia serial y estacionaridad                                   ##
##############################################################################
## Estacionaridad (idea): media y varianza constantes, y una dependencia
## serial que no cambia con t. Ejemplo: ancho de anillos de árboles (malleco).
X <- LSTS::malleco
n <- length(X)
par(op); par(bty = "n", las = 1)
plot(X, ylab = "Ancho de anillo", main = "Serie malleco (LSTS)")

## Medias por bloques: si la serie es estacionaria en media, las medias
## de bloques consecutivos deberían fluctuar alrededor de un mismo nivel.
n_bloques <- 8
tam       <- floor(n / n_bloques)             # tamaño de bloque desde los datos
medias    <- sapply(1:n_bloques, function(i) mean(X[1:tam + (i - 1) * tam]))
plot(medias, type = "b", pch = 19, xlab = "Bloque", ylab = "Media del bloque",
     main = "Medias por bloques", ylim = range(X))
abline(h = mean(X), lty = 2)

## ¿Tendencia determinística? Regresión sobre el tiempo
t <- 1:n
modelo <- lm(X ~ t)
summary(modelo)          # pendiente no significativa ==> sin tendencia clara

## ¿Varianza constante? Test de Breusch-Pagan sobre la regresión
lmtest::bptest(modelo)   # H0: homocedasticidad

## ACF global y por tramos: si la estructura de dependencia es estable,
## las ACF de los tres tramos deberían ser similares.
acf(X, lag.max = 10)
tercio <- floor(n / 3)
par(mfrow = c(1, 3))
acf(X[1:tercio],                  lag.max = 10, ylim = c(-1, +1), main = "Tramo 1")
acf(X[(tercio + 1):(2 * tercio)], lag.max = 10, ylim = c(-1, +1), main = "Tramo 2")
acf(X[(2 * tercio + 1):n],        lag.max = 10, ylim = c(-1, +1), main = "Tramo 3")

par(op)   # restauramos los parámetros gráficos originales
