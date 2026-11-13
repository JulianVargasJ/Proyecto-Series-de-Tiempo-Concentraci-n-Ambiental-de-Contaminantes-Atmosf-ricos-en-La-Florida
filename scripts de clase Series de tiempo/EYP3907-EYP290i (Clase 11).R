##############
## Clase 11 ##
##############

## Capitulo 2: Auto-correlacion PARCIAL (PACF)
## - Idea: correlacion entre y[t] e y[t-h] REMOVIENDO el efecto de los
##   rezagos intermedios y[t-1], ..., y[t-h+1]
## - Herramienta dual de la ACF para identificar:
##     AR(p): la PACF se CORTA en p (la ACF decae)
##     MA(q): la ACF se CORTA en q (la PACF decae)

## Requiere: LSTS, forecast

#####################################################################
## Paso 1: PACF teorica del AR(1) --> se corta en el rezago 1      ##
#####################################################################

phi <- 0.8
ACF  <- ARMAacf(ar = phi, lag.max = 10)
PACF <- ARMAacf(ar = phi, lag.max = 10, pacf = TRUE)

par(bty = "n", las = 1, mfrow = c(1,2))
plot(ACF ~ c(0:10), type = "h", ylim = c(-1,+1), xlab = "Lag", main = "ACF")
abline(h = 0)
plot(PACF ~ c(1:10), type = "h", ylim = c(-1,+1), xlab = "Lag", main = "PACF")
abline(h = 0)

## La ACF decae geometricamente (phi^h: nunca se apaga del todo),
## pero la PACF vale phi en el rezago 1 y CERO despues:
PACF

## ¿Por que? En un AR(1), toda la relacion entre y[t] e y[t-2] pasa
## POR y[t-1]; al remover ese efecto intermedio no queda nada

##################################
## Paso 2: malleco --> ¿AR(1) ? ##
##################################

X <- LSTS::malleco
par(bty = "n", las = 1, mfrow = c(1,2))
acf(X,  lag.max = 10, ylim = c(-1,+1), main = "ACF malleco")
pacf(X, lag.max = 10, ylim = c(-1,+1), main = "PACF malleco")

## Segun ACF --> MA(5), PACF --> AR(1): dejemos que compitan
forecast::auto.arima(X - mean(X), d = 0, max.p = 1, max.q = 5,
                     allowmean = FALSE, stepwise = FALSE)
## gana el AR(1): con ambas lecturas posibles, el criterio de informacion
## prefiere el modelo mas parsimonioso (1 parametro vs 5)

################################################################
## Paso 3: PACF teorica del AR(2) --> se corta en el rezago 2 ##
################################################################

phi <- c(1.2, -0.5)
ACF  <- ARMAacf(ar = phi, lag.max = 10)
PACF <- ARMAacf(ar = phi, lag.max = 10, pacf = TRUE)

par(bty = "n", las = 1, mfrow = c(1,2))
plot(ACF ~ c(0:10), type = "h", ylim = c(-1,+1), xlab = "Lag", main = "ACF")
abline(h = 0)
plot(PACF ~ c(1:10), type = "h", ylim = c(-1,+1), xlab = "Lag", main = "PACF")
abline(h = 0)

## La PACF en el rezago p coincide con el ULTIMO coeficiente del AR:
PACF[2]      ## = phi2 = -0.5
PACF[3:10]   ## = 0: el orden del AR queda "impreso" en la PACF

## Contraste: en un MA(1) es al reves (la ACF se corta, la PACF decae)
theta <- 0.8
round(ARMAacf(ma = theta, lag.max = 10), 3)               ## ACF: corta en 1
round(ARMAacf(ma = theta, lag.max = 10, pacf = TRUE), 3)  ## PACF: decae alternando

#####################################################################
## Paso 4: ¿que es "remover el efecto intermedio"? PACF via lm     ##
#####################################################################

## La PACF en el rezago h es el coeficiente de y[t-h] en la REGRESION
## de y[t] sobre y[t-1], ..., y[t-h]. Verifiquemoslo simulando un AR(2):

set.seed(3907)
n <- 2000
Y <- arima.sim(model = list(ar = c(1.2, -0.5)), n = n)

y0 <- Y[3:n]; y1 <- Y[2:(n-1)]; y2 <- Y[1:(n-2)]

## PACF(1): regresion con un rezago --> coeficiente de y[t-1]
coef(lm(y0 ~ y1))["y1"]
## PACF(2): regresion con dos rezagos --> coeficiente de y[t-2]
coef(lm(y0 ~ y1 + y2))["y2"]

## Comparar con la PACF muestral de R:
pacf(Y, lag.max = 2, plot = FALSE)$acf[,,1]
## practicamente identicas: la PACF ES una secuencia de regresiones

######################################################################
## Paso 5: deteccion empirica del orden p con las bandas de la PACF ##
######################################################################

## Bajo AR(p), las PACF muestrales de rezagos > p son ~ N(0, 1/n):
## las bandas azules de pacf() son +-1.96/sqrt(n)

par(bty = "n", las = 1, mfrow = c(1,2))
acf(Y,  lag.max = 20, ylim = c(-1,+1), main = "ACF muestral")
pacf(Y, lag.max = 20, ylim = c(-1,+1), main = "PACF muestral")

## Lectura: ACF que decae (suavemente, sin cortarse) + PACF con DOS
## barras significativas y el resto dentro de las bandas --> AR(2)

## Ojo con el azar: ~1 de cada 20 rezagos puede escapar de la banda
## aunque el modelo sea correcto (son bandas al 95%)

## Resumen (tabla de identificacion):
##            ACF                  PACF
## AR(p):     decae (geometrica)   se corta en p
## MA(q):     se corta en q        decae
## ARMA(p,q): decae                decae
## --> si ambas decaen sin cortarse, pensar en un ARMA
