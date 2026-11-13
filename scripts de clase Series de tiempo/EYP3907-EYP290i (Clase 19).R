##############
## Clase 19 ##
##############

## Capitulo 7: Prediccion
## - Predictor a un paso con pasado finito: Durbin-Levinson e innovaciones
## - La penalizacion por pasado finito muere rapido: nu_t --> sigma2
## - Predictores aproximados: truncar la expansion MA(inf) / AR(inf)
## - Prediccion a h pasos: sigma2(h) = sigma2 * sum(psi^2) --> gamma(0)

## Requiere: LSTS, forecast
## Archivos: Durbin_Levinson.R, innovation.R

########################################################
## Paso 1: tres caminos al mismo predictor a un paso  ##
########################################################

source("Durbin_Levinson.R")
source("innovation.R")

X <- LSTS::malleco
mod1 <- forecast::Arima(X - mean(X), order = c(0,0,3), include.mean = FALSE)

## (i) Kalman (lo que hace Arima), (ii) Durbin-Levinson, (iii) innovaciones
pre1 <- mod1$fitted + mean(X)
pre2 <- DurbinLevinson(X - mean(X), ar = NULL, ma = mod1$coef)$fitted + mean(X)
pre3 <- algoritmo.innov(X - mean(X), ar = NULL, ma = mod1$coef) + mean(X)

head(cbind(ts(pre1), ts(pre2), ts(pre3)), 12)
tail(cbind(ts(pre1), ts(pre2), ts(pre3)), 5)
## los tres coinciden: son el MISMO predictor optimo con pasado finito.
## Para un MA(3), el algoritmo de innovaciones usa a lo mas 3 innovaciones
## (theta[t,j] = 0 para j > 3): la ventaja de trabajar con terminos
## ortogonales en vez de las observaciones correlacionadas

##################################################################
## Paso 2: la penalizacion por pasado finito muere rapido       ##
##################################################################

## nu_t = gamma(0) prod(1 - phi_jj^2): decae al piso sigma2 en pocos pasos
nu <- DurbinLevinson(X - mean(X), ar = NULL, ma = mod1$coef)$nu
par(bty = "n", las = 1, mfrow = c(1,1))
plot(nu[1:15], type = "b", pch = 20, col = "blue",
     xlab = "t (largo del pasado observado)", ylab = expression(nu[t]))
abline(h = mod1$sigma2, col = "red", lty = 2)
## con 4 o 5 observaciones ya casi se alcanza sigma2 (pasado infinito):
## en memoria corta, "poca" historia basta. En ARFIMA (Cap. 2) este
## decaimiento seria mucho mas lento

########################################################
## Paso 3: predictores aproximados (truncar MA/AR inf)##
########################################################

n <- length(X)

## MA(inf): y.hat[t+1] = sum psi_j (y - y.hat) --> innovaciones "asintotico"
Psi <- ARMAtoMA(ar = NULL, ma = mod1$coef, lag.max = n)
Y <- X - mean(X)
Y.hat <- c(); Y.hat[1] <- 0
for(t in 1:(n-1)) Y.hat[t+1] <- sum((Y[t:1] - Y.hat[t:1])*Psi[1:t])
pre4 <- ts(Y.hat + mean(X))

## AR(inf): y.check[t+1] = sum pi_j y[t+1-j] --> el "predictor aproximado"
## de las slides (phi_tj se reemplaza por pi_j)
Pi <- -ARMAtoMA(ar = -mod1$coef, ma = NULL, lag.max = n)
X.hat <- c(); X.hat[1] <- 0
for(t in 1:(n-1)) X.hat[t+1] <- sum((X[t:1] - mean(X))*Pi[1:t])
pre5 <- ts(X.hat + mean(X))

head(cbind(ts(pre1), ts(pre2), ts(pre3), ts(pre4), ts(pre5)), 5)
tail(cbind(ts(pre1), ts(pre2), ts(pre3), ts(pre4), ts(pre5)), 5)

## ¿Cuanto cuesta la aproximacion? Var(e) sube menos de medio punto:
## el termino r_t de las slides es despreciable en memoria corta
(var(ts(X) - ts(pre4))/var(ts(X) - ts(pre1)) - 1)*100
(var(ts(X) - ts(pre5))/var(ts(X) - ts(pre1)) - 1)*100

##################################################################
## Paso 4: prediccion a h pasos --> sigma2(h) crece a gamma(0)  ##
##################################################################

h <- 1:10

## predict() entrega el se con pasado finito (Kalman); la formula de las
## slides sigma2(h) = sigma2 sum_{j<h} psi_j^2 es el pasado infinito
M <- cbind(predict = predict(mod1, n.ahead = 10)$se,
           MA.inf = sqrt(mod1$sigma2*cumsum(c(1, Psi^2))[h]))
round(M, 5)
## practicamente identicos (n es grande)

## Para un MA(3): psi_j = 0 si j > 3 --> sigma2(h) se ESTANCA en gamma(0)
## desde h = 4: mas alla de la memoria del proceso, el mejor pronostico
## es la media y la incertidumbre es la varianza del proceso
cbind(sigma2.h = mod1$sigma2*cumsum(c(1, Psi[1:9]^2)),
      gamma0 = rep(var(ts(X)), 10))

## Moraleja:
## - El predictor optimo a un paso es la proyeccion sobre el pasado:
##   Durbin-Levinson, innovaciones y Kalman son 3 formas del MISMO calculo.
## - Truncar la expansion AR(inf)/MA(inf) (predictor aproximado) casi no
##   cuesta en memoria corta: r_t es despreciable.
## - sigma2(h) = sigma2 sum(psi^2) crece hasta gamma(0): a largo plazo el
##   mejor pronostico es la media del proceso.
