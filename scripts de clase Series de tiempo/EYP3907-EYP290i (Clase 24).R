##############
## Clase 24 ##
##############

## Capitulo 9: Estacionalidad (continuacion)
## - SARFIMA: larga memoria estacional, ACF hiperbolica en multiplos de s
##   y densidad espectral que diverge en 2 pi j/s
##
## Requiere: LSTS

##################################################################
## Capitulo 9: SARFIMA                                          ##
##################################################################

########################################################
## Paso 1: SARFIMA --> larga memoria estacional       ##
########################################################

## (1 - B^s)^(-ds) = sum_j eta_j B^(s j), con los mismos eta_j del ARFIMA
## (Cap. 2) pero en los rezagos s, 2s, 3s, ...
eta <- function(j, d) exp(lgamma(j + d) - lgamma(j + 1) - lgamma(d))
ds <- 0.3; s <- 12
J <- 0:300
psi <- numeric(s*max(J) + 1)
psi[1 + s*J] <- eta(J, ds)              ## solo los multiplos de 12 son no nulos
round(psi[1 + s*(0:5)], 4)              ## 1, 0.300, 0.195, 0.150, ...

## Simulacion: MA(infinito) truncado, como en el ARFIMA
set.seed(3907)
n <- 1200
e <- rnorm(n + length(psi))
y.sf <- ts(stats::filter(e, psi, sides = 1)[-(1:length(psi))][1:n], frequency = s)

plot(y.sf,  bty = "n", col = "gray", las = 1)
## Para comparar: SAR(1)_12 con la misma correlacion en el rezago 12
rho12 <- acf(y.sf, lag.max = 12, plot = FALSE)$acf[13]
set.seed(3907)
y.sar <- arima.sim(list(ar = c(rep(0, 11), rho12)), n = n)

par(bty = "n", las = 1, mfrow = c(2,1), font.main = 1)
acf(c(y.sar), lag.max = 120, main = "SAR(1)[12]: decae geometricamente en 12, 24, 36, ...")
acf(c(y.sf),  lag.max = 120, main = "SARFIMA(ds = 0.3): decae hiperbolicamente")
## en el rezago 12k: rho(12k) ~ C k^(2 ds - 1) (memoria larga) vs Phi^k (corta)
a <- acf(c(y.sf), lag.max = 120, plot = FALSE)$acf
round(rbind(muestral = a[1 + 12*(1:10)],
            hiperbolica = a[13]*(1:10)^(2*ds - 1),
            geometrica  = a[13]^(1:10)), 3)

## Densidad espectral: f(lambda) = sigma^2/(2 pi) |1 - e^{-i 12 lambda}|^(-2 ds)
## diverge en CADA armonico estacional 2 pi j/12 (el ARFIMA solo en 0)
lam <- seq(0.005, pi, length = 3000)
f.sf <- Mod(1 - exp(-1i*s*lam))^(-2*ds)/(2*pi)
I <- LSTS::periodogram(c(y.sf), plot = FALSE)
par(mfrow = c(1,1))
plot(I$lambda, I$periodogram, type = "l", col = "gray60", ylim = c(0, 4),
     xlab = expression(omega), ylab = "", main = "SARFIMA: periodograma de la simulacion y f teorica")
lines(lam, f.sf, col = "blue", lwd = 2)
abline(v = 2*pi*(1:6)/12, lty = 3, col = "red")
## misma escala: el periodograma estima f (sigma^2 = 1) y f diverge en cada
## armonico 2 pi j/12 (el ARFIMA solo divergia en 0)
## mismos peaks que el co2 diferenciado... pero aqui son de larga memoria:
## no se eliminan con una diferencia entera, se fraccionan con ds

## Moraleja:
## - SARFIMA: la estacionalidad tambien puede tener larga memoria; la ACF
##   decae hiperbolicamente en 12, 24, ... y f diverge en cada 2 pi j/s.
## - Es el intermedio entre un SAR(1)_12 (memoria corta, Phi^k) y la raiz
##   unitaria estacional del co2 (ds = 1): persiste, pero revierte.
