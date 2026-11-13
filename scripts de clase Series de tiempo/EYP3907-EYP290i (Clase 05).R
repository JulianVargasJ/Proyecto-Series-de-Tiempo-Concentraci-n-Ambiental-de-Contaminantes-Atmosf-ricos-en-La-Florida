#####################
## Script Clase 05 ##
#####################

## Requiere: forecast, lmtest, LSTS

## Ausencia de Estacionaridad ##

## Heterocedasticidad --> Transformación de Box-Cox

## Ejemplo 1: Monthly Airline Passenger Numbers 1949-1960
Y <- AirPassengers
plot(Y)
lambda <- forecast::BoxCox.lambda(Y)
lambda
Y.bc <- forecast::BoxCox(Y, lambda = lambda)
Y.bc.INV <- forecast::InvBoxCox(Y.bc, lambda = lambda)
cbind(Y,Y.bc,Y.bc.INV)

Data <- data.frame(Y = c(Y), time = c(time(Y)), t = 1:length(Y), Y.bc = c(Y.bc))
head(Data)

## Estacionalidad con nls: regresión armónica
## Paso 1: tendencia sin estacionalidad; los residuos concentran el patrón estacional
plot(Y.bc ~ time, data = Data, type = "l")
mod_4 <- lm(Y.bc ~ I(t) + I(t^2), data = Data)
Data$fitted_4 <- fitted.values(mod_4)
lines(fitted_4 ~ time, data = Data, col = "blue", lty = 1)

Z <- mod_4$res
plot(ts(Z))

## Paso 2: ¿qué frecuencias dominan? --> periodograma (Ver video Transformada de Fourier)
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
## Modelo de referencia con mes como factor (el mod_3 de la Clase 04)
Data$month <- as.factor(cycle(AirPassengers))
mod_3 <- lm(Y.bc ~ I(t) + I(t^2) + month, data = Data)
Data$fitted_3 <- fitted.values(mod_3)

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


##############################################
## Tendencia Estocástica --> Diferenciación ##
##############################################

## Requiere: jsonlite (la API de mindicador.cl entrega JSON)

## Paso 1: descargar el dólar observado (USD/CLP) desde mindicador.cl
## jsonlite::fromJSON() lee la URL directamente: no se necesitan httr/dplyr/purrr
get_dolar_year <- function(year){
  jsonlite::fromJSON(paste0("https://mindicador.cl/api/dolar/", year))$serie
}

## Con respaldo local: si en la sala no hay internet, usa la última descarga
archivo_respaldo <- "dolar_respaldo.rds"
dolar <- tryCatch({
  aux <- rbind(get_dolar_year(2025), get_dolar_year(2026))
  aux$fecha <- as.Date(aux$fecha)
  aux <- aux[order(aux$fecha), c("fecha", "valor")]
  colnames(aux) <- c("fecha", "dolar")
  aux <- subset(aux, fecha >= Sys.Date() - 365)
  saveRDS(aux, archivo_respaldo)
  aux
}, error = function(e){
  message("Sin conexión: usando respaldo local ", archivo_respaldo)
  readRDS(archivo_respaldo)
})

head(dolar)
par(mfrow = c(1,1), bty = "n", las = 1)
plot(dolar ~ fecha, data = dolar, type = "l")

## Paso 2: ¿tendencia determinística o estocástica?
## (días hábiles: tratamos la serie como equiespaciada)
Y <- ts(dolar$dolar)
forecast::ndiffs(Y)   ## grado de integración sugerido por los tests de la Clase 06

## Paso 3: diferenciar una vez --> la serie de cambios diarios
eta <- diff(Y)        ## Y[t] - Y[t-1] = eta[t]
plot(eta)
mean(eta)             ## ¿deriva (drift)?
t.test(eta)$p.value   ## Ho: media cero --> sin deriva significativa
acf(eta)              ## ¿quedó (aprox.) no correlacionada?

## Paso 4: ¿podemos imitar al dólar? --> paseo aleatorio simulado
## Z iid Normal con la misma media y sd que eta, e Y*[t] = Y[1] + cumsum(Z)
set.seed(3907)
Z <- ts(rnorm(length(eta), mean(eta), sd(eta)))
par(mfrow = c(1,2))
plot(Z, main = "ruido simulado")
ts.plot(Y[1] + cumsum(Z), ylab = "", main = "paseo aleatorio simulado")
## A simple vista es indistinguible de un tipo de cambio real:
## la "tendencia" de un paseo aleatorio es 100% estocástica

## Paso 5: lo que la diferenciación NO arregla --> volatilidad (adelanto Cap. 5)
par(mfrow = c(1,2))
acf(Z^2,   main = expression(Z[t]^2~"(simulado)"))
acf(eta^2, main = expression(eta[t]^2~"(dólar)"))
## Los cambios del dólar casi no están correlacionados, pero sus CUADRADOS sí:
## volatilidad agrupada, que el ruido normal iid no puede replicar --> ARCH/GARCH

## Paso 6: sobrediferenciar tiene costo
## Diferenciemos un ruido que YA es estacionario...
par(mfrow = c(1,2))
plot(Z, ylim = c(-30,30), main = "Z")
plot(diff(Z), ylim = c(-30,30), main = expression((1-B)*Z[t]))
var(diff(Z))/var(Z)   ## ~2: duplica la varianza...
acf(diff(Z))          ## ...y crea correlación artificial: -0.5 en el lag 1
par(mfrow = c(1,1))
