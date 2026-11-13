##############
## Clase 21 ##
##############

## Capitulo 8: Procesos No Estacionarios
## - Tendencia determinista vs estocastica (se parecen, se tratan distinto)
## - Raiz unitaria: Dickey-Fuller a mano y su distribucion NO normal
## - ADF y ndiffs en la practica
## - ARIMA con datos reales: IPC
## - Back test: pronosticar a medida que llega la informacion (animacion)
## - Bonus: la dependencia de malleco cambia en el tiempo (proc. LS)

## Requiere: tseries, forecast, LSTS, rio (magick y pdftools, solo para el GIF)
## Archivos: IPC_EMP_2023.xlsx (carpeta Script de Clase)

##################################################################
## Capitulo 8: Procesos No Estacionarios                        ##
##################################################################

########################################################
## Paso 1: dos mecanismos que se parecen              ##
########################################################

set.seed(3907)
n <- 300
e <- rnorm(n)
det <- 0.02*(1:n) + e          ## determinista: f(t) + RB
sto <- cumsum(0.02 + e)        ## estocastica: paseo aleatorio con drift

par(bty = "n", las = 1, mfrow = c(1,2))
plot(det, type = "l", col = "gray40", xlab = "t", ylab = "",
     main = "f(t) + RB")
abline(a = 0, b = 0.02, col = "red", lwd = 2)
plot(sto, type = "l", col = "gray40", xlab = "t", ylab = "",
     main = "paseo aleatorio con drift")
abline(a = 0, b = 0.02, col = "red", lwd = 2)
## MISMAS innovaciones: en (a) la serie vuelve a la recta (shocks
## transitorios); en (b) se aleja sin castigo (shocks permanentes).
## Var: constante en (a), t*sigma2 en (b)

########################################################
## Paso 2: Dickey-Fuller a mano                       ##
########################################################

## H0: phi = 1 (raiz unitaria) vs H1: phi < 1
DF.stat <- function(y){
  n <- length(y)
  y1 <- y[1:(n-1)]; y2 <- y[2:n]
  ph <- sum(y2*y1)/sum(y1^2)
  s <- sqrt(sum((y2 - ph*y1)^2)/(n-2))
  (ph - 1)/(s/sqrt(sum(y1^2)))
}

## La distribucion de DF bajo H0 NO es Normal: simulemosla
set.seed(1)
DF <- replicate(3000, DF.stat(cumsum(rnorm(200))))
par(bty = "n", las = 1, mfrow = c(1,1))
hist(DF, breaks = 50, freq = FALSE, col = "gray90", border = "gray60",
     xlab = "DF", main = "")
curve(dnorm(x), add = TRUE, col = "blue", lwd = 2, lty = 2)
abline(v = quantile(DF, 0.05), col = "red", lwd = 2)
quantile(DF, c(0.01, 0.05, 0.10))
## el cuantil del 5% es -1.95 (tabla DF), no -1.64 (Normal): usar la
## tabla equivocada = rechazar la raiz unitaria de mas

## En la practica: ACF/PACF, ADF (agrega k rezagos de diff(y)) y ndiffs
par(mfrow = c(2,2), bty = "n", las = 1)
acf(sto, main = "sto"); pacf(sto, main = "sto")  ## ACF muere lento y PACF ~1
                                                 ## en el lag 1: AR(1) con phi ~ 1
acf(det, main = "det"); pacf(det, main = "det")  ## tambien muere lento: la
                                                 ## tendencia "contamina" la ACF

## ADF, H0: raiz unitaria (adf.test incluye constante y tendencia lineal)
tseries::adf.test(sto, k = 1)  ## p-valor 0.65: no se rechaza --> diferenciar
tseries::adf.test(det, k = 5)  ## p-valor < 0.01: se rechaza --> estacionaria
                               ## en torno a la tendencia (no diferenciar)

## ndiffs con tres tests: KPSS (H0: estacionaria), ADF y PP (H0: raiz unitaria)
forecast::ndiffs(ts(sto), test = "kpss")  ## 1
forecast::ndiffs(ts(sto), test = "adf")   ## 1
forecast::ndiffs(ts(sto), test = "pp")    ## 1
forecast::ndiffs(ts(det))                 ## 1 (!) el test por defecto no
                                          ## considera tendencia
## la estocastica pide 1 diferencia; la determinista, ojo: ndiffs con
## test que incluye tendencia deberia decir 0 --> se remueve por regresion
forecast::ndiffs(ts(det), type = "trend")

########################################################
## Paso 3: ARIMA con datos reales --> IPC             ##
########################################################

Data <- rio::import("IPC_EMP_2023.xlsx", skip = 2)   ## titulo + linea en blanco
colnames(Data) <- c("Date", "IPC")
Data$Date <- as.Date(Data$Date)
head(Data); tail(Data)
X <- ts(Data$IPC, start = c(2009,12), frequency = 12)   ## dic-2009 a sep-2026

par(bty = "n", las = 1, mfrow = c(1,1))
plot(X, lwd = 2, ylab = "", xlab = "", xlim = c(2010,2026), xaxt = "n")
axis(1, 2010:2026)

## ACF pegada a 1 y decayendo lineal: la firma del AR(1) con phi ~ 1
par(mfrow = c(1,2), bty = "n", las = 1)
acf(c(X), ylim = c(-1,1), xlim = c(0,12), main = "")
pacf(c(X), ylim = c(-1,1), xlim = c(0,12), main = "")

## ¿Cuantas diferencias?
forecast::ndiffs(X, test = "kpss")  ## 2
forecast::ndiffs(X, test = "adf")   ## 1
forecast::ndiffs(X, test = "pp")    ## 1
## los tests no coinciden: KPSS pide d = 2, ADF y PP piden d = 1

## Candidatos (ojo: el AIC NO es comparable entre distintos d, porque
## cambia la serie que se modela)
forecast::Arima(X, order = c(1,0,0), method = "ML")  ## phi = 0.9999 ~ 1 --> d = 1
forecast::Arima(X, order = c(0,1,1))  ## propuesta ADF y PP
forecast::Arima(X, order = c(0,2,1))  ## propuesta KPSS (y auto.arima): theta = -0.83
forecast::Arima(X, order = c(0,3,1))  ## una de mas: theta = -1, raiz MA
                                      ## unitaria --> sobrediferenciacion

## Sobrediferenciar infla la varianza (los dos errores posibles, Paso 5)
var(diff(X, differences = 1))
var(diff(X, differences = 2))  ## +49% respecto a d = 1
var(diff(X, differences = 3))  ## +207% respecto a d = 2 (+355% respecto a d = 1)
## la regla de la varianza minima apunta a d = 1, pero KPSS y auto.arima
## eligen d = 2: el alza de 2021-2022 hace que la inflacion (d = 1)
## tampoco se vea estacionaria. No hay una unica respuesta "correcta"


## Diferenciar y modelar lo estacionario (Cap. 5 sobre diff(X))
## (izq.) variacion mensual del indice; (der.) su cambio
par(mfrow = c(1,2), bty = "n", las = 1)
plot(diff(X, differences = 1), ylab = "", xlab = "", ylim = c(-1.5,2))
plot(diff(X, differences = 2), ylab = "", xlab = "", ylim = c(-1.5,2))

## en log: diff(log(X)) es la inflacion mensual (aprox. en %/100)
par(mfrow = c(1,2), bty = "n", las = 1)
plot(diff(log(X), differences = 1), ylab = "", xlab = "", ylim = c(-1.5,2)/100)
plot(diff(log(X), differences = 2), ylab = "", xlab = "", ylim = c(-1.5,2)/100)


fit <- forecast::auto.arima(X, seasonal = FALSE)  ## el circuito completo
summary(fit)

## Prediccion: la banda se abre SIN estabilizarse (shocks permanentes)
par(mfrow = c(1,1))
pre <- forecast::forecast(fit, h = 24, level = 0.95)
plot(pre, main = "", bty = "n", las = 1)
## contraste con el Cap. 7: en un proceso estacionario la banda se
## estabiliza en +-z sqrt(gamma(0)); aqui crece sin tope

##################################################################
## Paso 4: back test --> pronosticar a medida que llega la data ##
##################################################################

## Nos paramos en cada mes k como si fuera "hoy": solo conocemos
## X[1:k], pronosticamos 12 meses y comparamos con lo que realmente
## ocurrio. Es la validacion temporal del proyecto, repetida muchas
## veces, y de paso una animacion del pronostico en el tiempo.

REESTIMAR   <- FALSE   ## TRUE: reestima el modelo en cada origen (honesto)
                       ## FALSE: usa los coeficientes del ajuste completo
GUARDAR_GIF <- FALSE   ## TRUE: guarda la animacion en backtest_ipc.gif
H <- 12                ## horizonte de pronostico

ord <- forecast::arimaorder(fit)                 ## orden elegido en el Paso 3
origenes <- which(time(X) >= 2015 & seq_along(X) < length(X))

pres <- vector("list", length(origenes))         ## pronosticos por origen
err  <- matrix(NA, length(origenes), H)          ## error a h pasos
cub  <- matrix(NA, length(origenes), H)          ## ¿cayo dentro de la banda?

for(i in seq_along(origenes)){
  k <- origenes[i]
  X.k <- window(X, end = time(X)[k])
  fit.k <- if(REESTIMAR) forecast::Arima(X.k, order = ord) else
                         forecast::Arima(X.k, model = fit)
  pres[[i]] <- forecast::forecast(fit.k, h = H, level = 95)
  m <- min(H, length(X) - k)                     ## horizontes ya observados
  obs <- X[(k+1):(k+m)]
  err[i, 1:m] <- obs - pres[[i]]$mean[1:m]
  cub[i, 1:m] <- obs >= pres[[i]]$lower[1:m] & obs <= pres[[i]]$upper[1:m]
}

## Un cuadro de la animacion: pronostico desde el origen i vs lo real
cuadro <- function(i){
  k <- origenes[i]
  plot(pres[[i]], xlim = c(2014, max(time(X)) + 1), ylim = c(60, 125),
       main = sprintf("Origen: %d-%02d", floor(time(X)[k] + 1e-6), cycle(X)[k]),
       bty = "n", las = 1, font.main = 1, ylab = "IPC")
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
                  cobertura = colMeans(cub, na.rm = TRUE))
round(res, 3)
## la cobertura queda bajo el 95%: ~86-90% con coeficientes fijos
## (REESTIMAR = FALSE, que ademas usa informacion del futuro) y ~80%
## reestimando: las bandas son demasiado optimistas

par(mfrow = c(1,2), bty = "n", las = 1)
plot(res$h, res$RMSE, type = "b", pch = 20, xlab = "h", ylab = "RMSE")
plot(res$h, res$cobertura, type = "b", pch = 20, ylim = c(0, 1),
     xlab = "h", ylab = "cobertura banda 95%")
abline(h = 0.95, col = "red", lty = 2)

## ¿Donde falla? Error a 12 meses segun el origen: el alza de la
## inflacion 2021-2022 no la anticipa ningun ARIMA ajustado con el pasado
par(mfrow = c(1,1))
plot(time(X)[origenes], err[, H], type = "h", xlab = "origen",
     ylab = "error a 12 meses (real - pronostico)")
abline(h = 0)

## Opcional: guardar la animacion como GIF (requiere magick y pdftools)
if(GUARDAR_GIF){
  pdf("backtest_ipc.pdf", 10, 5.6)
  for(i in seq_along(origenes)) cuadro(i)
  dev.off()
  img <- magick::image_read_pdf("backtest_ipc.pdf", density = 90)
  magick::image_write(magick::image_animate(img, fps = 5), "backtest_ipc.gif")
}

########################################################
## Paso 5: los dos errores posibles                   ##
########################################################

## (a) sobrediferenciar una serie estacionaria (la leccion de la I1):
z <- diff(det)                  ## la determinista NO habia que diferenciarla
c(var(det - 0.02*(1:n)), var(z))     ## la varianza ~ se duplica
acf(z, lag.max = 5, plot = FALSE)$acf[2]   ## y aparece rho(1) ~ -0.5 artificial

## (b) NO diferenciar una raiz unitaria: toda la teoria de los capitulos
## anteriores (estacionaria) queda sin piso: ACF que no muere, regresiones
## espurias, bandas invalidas

##################################################################
## Paso 6 (bonus): la dependencia de malleco cambia en el tiempo##
##################################################################

## Procesos localmente estacionarios (el apellido de LSTS):
## rho(1) en ventanas moviles de 100 anos
X <- LSTS::malleco
n <- length(X)
w <- 100
tt <- w:n
r1 <- sapply(tt, function(i) acf(as.numeric(X[(i-w+1):i]),
                                 plot = FALSE, lag.max = 1)$acf[2])
par(bty = "n", las = 1, mfrow = c(1,1))
plot(c(time(X))[tt], r1, type = "l", col = "blue", lwd = 2,
     xlab = "ano (fin de la ventana)", ylab = expression(hat(rho)(1)))
abline(h = acf(c(X), plot = FALSE, lag.max = 1)$acf[2], lty = 2, col = "red")
## la dependencia se fortalece hacia el presente (~0.55 --> ~0.77): el
## AR(1) "constante" del curso es el promedio de una estructura que
## evoluciona. Modelos para esto: LS-ARMA / LS-ARFIMA (LSTS::LS.whittle)

## Moraleja:
## - Tendencia determinista: regresion. Tendencia estocastica: diferencia.
##   Se parecen a simple vista; DF/ADF y ndiffs ayudan a decidir.
## - La distribucion de DF no es Normal: cuantiles propios (-1.95 al 5%).
## - ARIMA = diferenciar lo justo + ARMA sobre lo estacionario; su banda
##   de prediccion crece sin estabilizarse.
## - El back test (pronosticar con lo que se sabia en cada momento) mide
##   el error real a cada horizonte y la cobertura de las bandas: el
##   modelo no anticipa quiebres como el alza de 2021-2022.
## - Y si los parametros cambian suavemente en el tiempo: procesos
##   localmente estacionarios (LSTS).
