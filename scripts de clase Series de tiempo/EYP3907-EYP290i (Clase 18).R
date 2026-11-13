##############
## Clase 18 ##
##############

## Capitulo 6 (continuacion): GARCH
## - Ajuste GARCH(1,1): innovaciones normales vs t de Student
## - La banda +-2 sigma_t que "respira" (figura de las slides)
## - Asimetria: TGARCH sobre el tipo de cambio USD/EUR
## - Bonus: simulacion de escenarios de inversion con garchSim

## Requiere: fGarch, quantmod
## Archivos: Tgarch11.R, d-useu9910.txt (archivo Tsay.zip)

##################################################################
## Paso 1: ajuste GARCH(1,1) --> normal vs t de Student         ##
##################################################################

## Datos de la Clase 17: rentabilidad diaria del Dow Jones
library(fGarch)
library(quantmod)
getSymbols("^DJI", src = "yahoo", from = "2025-01-01", to = Sys.Date())
n <- dim(DJI)[1]
X <- ts(DJI$DJI.Adjusted)
Y <- c(X[2:n]/X[1:(n-1)] - 1)

ts.plot(Y)
acf(Y)
acf(Y^2)

fit1 <- garchFit(formula = ~ arma(0,0) + garch(1,1), data = Y,
                 cond.dist = "norm")
summary(fit1)
## Que mirar: alpha1 + beta1 ~ 0.9x (persistencia, cerca de IGARCH);
## los test de Ljung-Box sobre R^2 (residuos^2) ya no rechazan

fit2 <- garchFit(formula = ~ arma(0,0) + garch(1,1), data = Y,
                 cond.dist = "std")
summary(fit2)
## "shape" = grados de libertad de la t: pequenos --> colas pesadas;
## comparar AIC/BIC de la tabla "Information Criterion Statistics"

## plot(fit2) es interactivo: utiles 1 (serie), 3 (serie con +-2 sd
## condicional), 11 (ACF de residuos^2), 13 (QQ-plot)
## plot(fit2)

##################################################################
## Paso 2: la banda +-2 sigma_t que "respira" (slides)          ##
##################################################################

s <- fit2@sigma.t
par(mfrow = c(1,1), bty = "n", las = 1)
ts.plot(Y, col = "gray40")
lines(+2*s, col = "red")
lines(-2*s, col = "red")
abline(h = c(-2,2)*sd(Y, na.rm = TRUE), col = "blue", lty = 2)
legend("topright", legend = c("+-2 sigma_t (GARCH)", "+-2 sd (incondicional)"),
       col = c("red", "blue"), lty = c(1,2), bty = "n")
## la banda GARCH se abre en las crisis y se cierra en la calma;
## la banda incondicional (azul) es ciega a las rachas

##################################################################
## Paso 3: asimetria --> TGARCH sobre USD/EUR                   ##
##################################################################

## Malas noticias pesan mas: sigma_t^2 salta mas tras una caida que
## tras un alza de igual magnitud (efecto apalancamiento)

source("Tgarch11.R")
Data <- rio::import("d-useu9910.txt")
Data$date <- as.Date(paste(Data[,1], Data[,2], Data[,3], sep = "-"))
n <- dim(Data)[1]
Data$return <- c(NA, Data$rate[2:n]/Data$rate[1:(n-1)] - 1)

par(bty = "n", las = 1, mfrow = c(1,2))
plot(rate ~ date, data = Data, type = "l")
plot(return ~ date, data = Data, type = "l")

X <- ts(Data$return[-1])
mod <- Tgarch11(X)
## mirar gam1 en la tabla: si es significativo y positivo, una caida
## agrega (alpha + gam1) y[t-1]^2 a la varianza y un alza solo alpha y[t-1]^2

## Curva de impacto de noticias con los parametros estimados
## (Tgarch11 devuelve $par = (mu, omega, alpha, gam1, beta))
cf <- mod$par
a0 <- cf[2]; a1 <- cf[3]; g1 <- cf[4]; b1 <- cf[5]
u  <- seq(-3, 3, length = 200)*sd(X)
par(bty = "n", las = 1, mfrow = c(1,1))
plot(u, a0 + a1*u^2 + g1*(u < 0)*u^2 + b1*var(X), type = "l", col = "red",
     lwd = 2, xlab = "shock y[t-1]", ylab = "sigma[t]^2")
lines(u, a0 + (a1 + g1/2)*u^2 + b1*var(X), col = "blue", lwd = 2)
abline(v = 0, lty = 3)
legend("top", legend = c("TGARCH (asimetrica)", "GARCH equivalente (simetrica)"),
       col = c("red","blue"), lwd = 2, bty = "n")

##################################################################
## Paso 4 (bonus): escenarios de inversion con garchSim         ##
##################################################################

## Ahorro mensual con abono, retornos GARCH-t estimados del Dow Jones:
## garchSpec espera 'alpha' y 'beta' como VECTORES sin nombres

spec <- garchSpec(model = list(
  mu    = coef(fit2)[1],
  omega = coef(fit2)[2],
  alpha = unname(coef(fit2)[3]),
  beta  = unname(coef(fit2)[4]),
  shape = coef(fit2)[5]), cond.dist = "std")

final <- c()
for(j in 1:200){
  ini <- 100000
  abono <- 400000
  for(i in 1:(12*5)){
    Z <- garchSim(spec, n = 20)$garch
    ini <- abono + ini * cumprod(Z + 1)[20]
  }
  final[j] <- ini
}

par(mfrow = c(1,1), bty = "n", las = 1)
boxplot(final/1000000, las = 1, pch = 20, ylab = "MM CL$")
Rentabilidad <- mean(final, na.rm = TRUE)/(100000 + abono*12*5) - 1
paste("Rentabilidad Esperada ", round(100*Rentabilidad), "%", sep = "")
## el modelo no solo describe: permite SIMULAR escenarios realistas
## (con rachas de volatilidad y colas pesadas) para decisiones

## Moraleja:
## - ARCH/GARCH generan ruido blanco con dependencia en los cuadrados:
##   la deteccion es la ACF/Box-Ljung de y^2 (la P3 de la I1).
## - alpha1 + beta1 mide la persistencia de la volatilidad; en datos
##   reales suele quedar cerca de 1 (IGARCH).
## - Innovaciones t de Student acomodan las colas pesadas; el propio
##   GARCH gaussiano ya tiene curtosis > 3.
## - TGARCH/EGARCH capturan la asimetria: las malas noticias mueven
##   mas la volatilidad que las buenas.
