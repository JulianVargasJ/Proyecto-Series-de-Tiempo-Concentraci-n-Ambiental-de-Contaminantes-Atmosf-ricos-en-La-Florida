###################################
## Clase 12                      ##
## Introduccion a Larga Memoria  ##
###################################

## Capitulo 2: Procesos de Larga Memoria
## - Ruido fraccionario: ACF con decaimiento HIPERBOLICO rho(k) ~ k^(2d-1)
## - La firma practica: la media muestral converge LENTO a la media
## - Deteccion (metodos intuitivos): Variance Plot y R/S
## - Estimacion por maxima verosimilitud: fracdiff
## - Simulacion: sesgo y ECM de los estimadores

## Requiere: longmemo, fracdiff

##################################################################
## Paso 1: rio Nilo --> ¿no estacionario o larga dependencia?   ##
##################################################################

library(longmemo)
data(NileMin)
X <- ts(NileMin, start = 622)

par(bty = "n", las = 1, mfrow = c(1,2))
ts.plot(X, ylab = "", ylim = c(800,1600), col = "gray")
abline(h = mean(X), lty = 2)
## la serie "deambula" alrededor de una media fija: rachas largas sobre
## y bajo la media, pero SIN tendencia --> sospecha de larga dependencia

acf(X, main = "", lag.max = 30)
## ACF que decae lentisimo y no se corta: ni AR(p) ni MA(q) calzan

## ACF teorica del ruido fraccionario FN(d), d in (0, 1/2):
## rho(k) = Gamma(1-d) Gamma(k+d) / [Gamma(d) Gamma(1+k-d)] ~ C k^(2d-1)
k <- 1:30
col.d <- c("red", "orange", "blue")
dd <- c(0.30, 0.35, 0.40)
for(i in 1:3){
  d <- dd[i]
  ACF <- gamma(1-d)*gamma(k+d) / (gamma(d)*gamma(1+k-d))
  lines(ACF ~ k, type = "p", col = col.d[i], pch = 20)
}
legend("topright", legend = paste("d =", rev(dd)), pch = 20,
       col = rev(col.d), bty = "n")
## a mayor d, mas lento el decaimiento; el Nilo se parece a un d ~ 0.4

###############################################################
## Paso 2: la firma practica --> la media muestral converge  ##
##         LENTO bajo larga memoria (figura de la clase)     ##
###############################################################

set.seed(3907)
n  <- 3000
fn <- fracdiff::fracdiff.sim(n = n, d = 0.4)$series   ## FN(d = 0.4)
ar <- arima.sim(model = list(ar = 0.5), n = n)        ## AR(1) (memoria corta)
fn <- fn/sd(fn); ar <- ar/sd(ar)                      ## misma escala

par(bty = "n", las = 1, mfrow = c(1,2))
## (a) media muestral acumulada: el AR(1) se "pega" rapido a la media,
##     el FN sigue deambulando
m.fn <- cumsum(fn)/(1:n)
m.ar <- cumsum(ar)/(1:n)
plot(m.fn ~ c(1:n), type = "l", col = "blue", xlab = "n",
     ylab = expression(bar(y)[n]), ylim = range(c(m.fn[-(1:10)], m.ar[-(1:10)])))
lines(m.ar ~ c(1:n), col = "red")
abline(h = 0, lty = 2)
legend("topright", legend = c("FN(d = 0.4)", "AR(1)"),
       col = c("blue","red"), lty = 1, bty = "n")

## (b) varianza de medias por bloques de tamaño m, en escala log-log:
##     memoria corta --> pendiente -1; larga memoria --> pendiente 2d-1
vb <- function(x, m){ k <- floor(length(x)/m)
                      var(colMeans(matrix(x[1:(k*m)], nrow = m))) }
M <- round(exp(seq(log(5), log(300), length = 15)))
v.fn <- sapply(M, vb, x = fn)
v.ar <- sapply(M, vb, x = ar)
plot(log(v.fn) ~ log(M), pch = 20, col = "blue", xlab = "log m",
     ylab = expression("log Var"(bar(y)[m])), ylim = range(log(c(v.fn, v.ar))))
points(log(v.ar) ~ log(M), pch = 20, col = "red")
abline(lm(log(v.fn) ~ log(M)), col = "blue")
abline(lm(log(v.ar) ~ log(M)), col = "red")

## Ojo: con UNA trayectoria la pendiente empirica del FN sale mas negativa
## que la teorica 2d-1 = -0.2 (el estimador tiende a subestimar d); la
## figura de las slides usa las curvas EXACTAS calculadas desde la ACF

########################################################
## Paso 3: deteccion 1 --> Variance Plot              ##
########################################################

## Es exactamente el panel (b) anterior convertido en estimador:
## log Var(ybar_m) ~ c + (2d-1) log m  -->  d = (1 + pendiente)/2

varlm <- function(x, m){
  n <- length(x)
  n1 <- min(trunc(n/2), m + 1)
  s2 <- c(1:(n1 - 1))
  for(k in 2:n1){
    m <- trunc(n/k)
    xbark <- c(1:m)
    for(j in 1:m)
      xbark[j] <- mean(x[((j-1)*k + 1):(j*k)])
    s2[k-1] <- var(xbark, na.rm = TRUE)
  }
  s2
}

varlm1 <- function(aux, m, linety = 2, plot = TRUE){
  x <- log(c(1:m))
  y <- log(aux[1:m])
  fit <- lsfit(x, y)$coef
  b <- mean(x + y)
  if(plot){
    par(bty = "n", font.main = 1, las = 1, mfrow = c(1,1))
    plot(x, y, xlab = expression(log(k)), ylab = expression("log Var"(k)))
    lines(x, -1*x + b, lty = linety)      ## referencia: memoria corta (pendiente -1)
    lines(x, fit[2]*x + fit[1])           ## ajuste: pendiente 2d-1
  }
  (1 + fit[2]) / 2                        ## (1 + pendiente)/2 = d
}

varlm2 <- function(x, m, linety = 2, plot = TRUE){
  ## x: la serie; m: numero de puntos de la regresion (p.ej. m = 10)
  varlm1(varlm(x, m), m, linety, plot = plot)
}

d <- varlm2(x = NileMin, m = 10)
M <- cbind(d); colnames(M) <- "variance plot"; rownames(M) <- "d"
M

## Contraste con una serie de memoria CORTA: d deberia salir ~ 0
varlm2(x = diff(LakeHuron), m = 10)

## Contraste con una serie de memoria CORTA: d deberia salir ~ 0
varlm2(x = diff(LSTS::malleco), m = 10)


##########################################################
## Paso 4: deteccion 2 --> Rescaled Range (R/S), Hurst  ##
##########################################################

## Idea del embalse (slide de la clase): R mide la capacidad necesaria
## para sostener una salida constante dado el caudal acumulado; R/s crece
## como k^(d + 1/2) --> la pendiente de log(R/S) vs log(k) estima d + 1/2

rescaled_range_statistic <- function(y){
  n <- length(y)
  x <- cumsum(y)
  adjusted_x <- x - (1:n) * mean(y)
  R <- max(adjusted_x) - min(adjusted_x)
  R / sd(y)
}

RS <- function(X, plot = TRUE){
  RS <- c(); K <- c()
  n <- length(X)
  for(k in seq(10, 200, 10)){
    for(t in 1:(n-k+1)){
      RS <- c(RS, rescaled_range_statistic(y = X[t:(t+k-1)]))
      K  <- c(K, k)
    }
  }
  if(plot){
    par(bty = "n", font.main = 1, las = 1, mfrow = c(1,1))
    plot(log(RS) ~ log(K), pch = 20, col = "gray")
    abline(lm(log(RS) ~ log(K)))
  }
  lm(log(RS) ~ log(K))$coef[2] - 1/2   ## pendiente - 1/2 = d
}

d <- RS(X)
M <- cbind(d, M); colnames(M)[1] <- "rescaled range"; rownames(M) <- "d"
M

########################################################
## Paso 5: estimacion por ML --> fracdiff             ##
########################################################

## Maxima verosimilitud (aproximada) del modelo ARFIMA(0,d,0)
d <- fracdiff::fracdiff(X)$d
M <- cbind(d, M); colnames(M)[1] <- "fracdiff"; rownames(M) <- "d"
M
## los tres metodos coinciden en d ~ 0.35 - 0.40 para el Nilo

###########################################################
## Paso 6: simulacion --> sesgo y ECM de los estimadores ##
###########################################################

## Simulamos 20 trayectorias ARFIMA(0, d = 0.4, 0) de largo 1000
## y estimamos d con los tres metodos

set.seed(12345)
M <- c()
for(i in 1:20){
  X <- fracdiff::fracdiff.sim(n = 1000, d = 0.4)$series
  d1 <- fracdiff::fracdiff(X)$d
  d2 <- varlm2(X, plot = FALSE, m = 10)
  d3 <- RS(X, plot = FALSE)
  M <- rbind(M, cbind(fracdiff = d1, varplot = d2, RS = d3))
}
M

apply(M, 2, mean)

Sesgo <- apply(M, 2, mean) - 0.4
Sesgo

ECM <- Sesgo^2 + apply(M, 2, var)
ECM

## Moraleja:
## - fracdiff (ML) es el mas preciso (menor ECM), pero asume el modelo.
## - Variance Plot y R/S son graficos y heuristicos: mas sesgados y
##   variables, pero no requieren especificar un modelo y sirven para
##   DETECTAR la larga memoria antes de modelarla.
