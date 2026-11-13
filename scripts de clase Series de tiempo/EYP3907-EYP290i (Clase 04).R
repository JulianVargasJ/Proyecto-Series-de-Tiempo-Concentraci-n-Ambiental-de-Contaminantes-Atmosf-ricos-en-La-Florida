#####################
## Script Clase 04 ##
#####################

## Requiere: forecast, lmtest, LSTS

## Ausencia de Estacionaridad ##

## Heterocedasticidad --> Transformación de Box-Cox

## Ejemplo 1: Monthly Airline Passenger Numbers 1949-1960
AirPassengers
start(AirPassengers)
end(AirPassengers)
frequency(AirPassengers)

## Sin atributos de serie de tiempo
X <- c(AirPassengers)
plot(X)
## El orden temporal ES la información: si lo destruimos, no queda nada
plot(sample(X))

## Con atributos ts(...)
X <- ts(c(AirPassengers), start = c(1949,1), frequency = 12)
X
plot(X)

## Coeficiente de variación por año:
## si sd crece proporcional a la media (cv ~ constante en la serie original),
## la estructura es multiplicativa --> transformar
cv <- function(x){sd(x)/mean(x)}

aggregate(X ~ trunc(time(X)), FUN = mean)
aggregate(X ~ trunc(time(X)), FUN = sd)
M <- aggregate(X ~ trunc(time(X)), FUN = cv)
colnames(M) <- c("year","cv")
plot(cv ~ year, data = M)

## Box-Cox con lambda = 0.5 (raíz cuadrada): ¿estabiliza?
Y <- forecast::BoxCox(X, lambda = 0.5)
M <- aggregate(Y ~ trunc(time(Y)), FUN = cv)
colnames(M) <- c("year","cv")
plot(cv ~ year, data = M)

## lambda "óptimo" (método de Guerrero)
lambda <- forecast::BoxCox.lambda(X)
lambda

Y <- forecast::BoxCox(X, lambda = lambda)
M <- aggregate(Y ~ trunc(time(Y)), FUN = cv)
colnames(M) <- c("year","cv")
plot(cv ~ year, data = M, ylim = c(0.005,0.015))

par(mfrow = c(1,2))
plot(X)
plot(Y)
par(mfrow = c(1,1))

## lambda = 0 (logaritmo): casi igual de bueno e interpretable (contexto económico)
Y <- forecast::BoxCox(X, lambda = 0)
M <- aggregate(Y ~ trunc(time(Y)), FUN = cv)
colnames(M) <- c("year","cv")
plot(cv ~ year, data = M, ylim = c(0.01,0.04))

par(mfrow = c(1,2))
plot(X)
plot(Y)
par(mfrow = c(1,1))

## Modelamiento Tendencial y Estacional
Data <- data.frame(Y = c(AirPassengers), time = c(time(AirPassengers)), t = 1:length(AirPassengers))
Data$Y.bc <- c(forecast::BoxCox(AirPassengers, lambda = lambda))

## Tendencia lineal: regresor "time" (años) vs regresor "t" (1,2,3,...)
## OJO: los coeficientes cambian, pero los valores ajustados son IDÉNTICOS
## (t es transformación lineal de time)
plot(Y ~ time, data = Data, type = "l")
mod_1 <- lm(Y ~ time, data = Data)
Data$fitted_1 <- fitted.values(mod_1)
mod_2 <- lm(Y ~ t, data = Data)
Data$fitted_2 <- fitted.values(mod_2)
lines(fitted_1 ~ time, data = Data, col = "red")
lines(fitted_2 ~ time, data = Data, col = "blue", lty = 2)
cbind(coef(mod_1), coef(mod_2))
max(abs(Data$fitted_1 - Data$fitted_2))

## Lo mismo sobre la serie transformada
plot(Y.bc ~ time, data = Data, type = "l")
mod_1 <- lm(Y.bc ~ time, data = Data)
Data$fitted_1 <- fitted.values(mod_1)
mod_2 <- lm(Y.bc ~ t, data = Data)
Data$fitted_2 <- fitted.values(mod_2)
lines(fitted_1 ~ time, data = Data, col = "red")
lines(fitted_2 ~ time, data = Data, col = "blue", lty = 2)

## Tendencia cuadrática
plot(Y.bc ~ time, data = Data, type = "l")
mod_1 <- lm(Y.bc ~ I(time)+I(time^2), data = Data)
Data$fitted_1 <- fitted.values(mod_1)
mod_2 <- lm(Y.bc ~ I(t) + I(t^2), data = Data)
Data$fitted_2 <- fitted.values(mod_2)
lines(fitted_1 ~ time, data = Data, col = "red")
lines(fitted_2 ~ time, data = Data, col = "blue", lty = 2)

## Estacionalidad con as.factor
## cycle(...) entrega el mes de cada observación (equivale a rep(1:12, 12))
Data$month <- as.factor(cycle(AirPassengers))
plot(Y.bc ~ time, data = Data, type = "l")
mod_3 <- lm(Y.bc ~ I(t) + I(t^2) + month, data = Data)
Data$fitted_3 <- fitted.values(mod_3)
lines(fitted_3 ~ time, data = Data, col = "blue", lty = 1)

## Volver a la escala original: Box-Cox inversa
Data$fitted_3_bcinv <- forecast::InvBoxCox(Data$fitted_3, lambda = lambda)
plot(Y ~ time, data = Data, type = "l")
lines(fitted_3_bcinv ~ time, data = Data, col = "blue", lty = 1)

## Diagnóstico del modelo con mes
Z <- mod_3$res
par(mfrow = c(2,2))
plot(mod_3)                  ## diagnóstico clásico de regresión
par(mfrow = c(1,1))
lmtest::bptest(mod_3)        ## Ho: homocedasticidad
acf(Z)                       ## dependencia serial remanente --> la modelaremos con ARMA (próximas clases)

## Estacionalidad con nls: regresión armónica
## Paso 1: tendencia sin estacionalidad; los residuos concentran el patrón estacional
plot(Y.bc ~ time, data = Data, type = "l")
mod_4 <- lm(Y.bc ~ I(t) + I(t^2), data = Data)
Data$fitted_4 <- fitted.values(mod_4)
lines(fitted_4 ~ time, data = Data, col = "blue", lty = 1)

Z <- mod_4$res
plot(ts(Z))

## Paso 2: ¿qué frecuencias dominan? --> periodograma
aux <- LSTS::periodogram(Z)
i <- rev(order(aux$periodogram))[1:2]
## Periodos detectados (en meses): 12 (ciclo anual) y su armónico 6
2*pi/aux$lambda[i]

## Paso 3: ajustar la regresión armónica con los periodos detectados
per <- round(2*pi/aux$lambda[i])
mod_5 <- nls(Y.bc ~ b0 + b1*t + b2*t^2 +
               alpha1*sin(2*pi*t/12) + beta1*cos(2*pi*t/12) +
               alpha2*sin(2*pi*t/6)  + beta2*cos(2*pi*t/6),
             data = Data,
             start = list(b0 = mean(Data$Y.bc), b1 = 0, b2 = 0,
                          alpha1 = 0, beta1 = 0, alpha2 = 0, beta2 = 0))
summary(mod_5)

plot(Y.bc ~ time, data = Data, type = "l")
lines(fitted(mod_5) ~ time, data = Data, col = "red", lwd = 2)

## Paso 4: comparar armónicos (4 parámetros estacionales) vs mes factor (11)
Data$fitted_5 <- fitted(mod_5)
plot(Y.bc ~ time, data = Data, type = "l")
lines(fitted_3 ~ time, data = Data, col = "blue")
lines(fitted_5 ~ time, data = Data, col = "red", lty = 2)
legend("topleft", legend = c("mes factor (11 par.)", "armónicos 12 y 6 (4 par.)"),
       col = c("blue","red"), lty = c(1,2), bty = "n")

## En escala original
plot(Y ~ time, data = Data, type = "l")
lines(forecast::InvBoxCox(Data$fitted_5, lambda = lambda) ~ time, data = Data, col = "red")

## Paso 5: diagnóstico de los residuos armónicos
Z5 <- Data$Y.bc - Data$fitted_5
par(mfrow = c(1,2))
plot(ts(Z5))
acf(Z5)                      ## aún hay dependencia serial: motivación para los modelos ARMA
par(mfrow = c(1,1))
