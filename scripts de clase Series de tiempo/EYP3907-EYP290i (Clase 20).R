##############
## Clase 20 ##
##############

## Capitulo 7 (continuacion): Prediccion
## - Interpolacion de datos faltantes: suavizador de Kalman
## - Bandas de prediccion; bandas asimetricas con transformacion log

## Requiere: forecast, LSTS

##################################################################
## Capitulo 7 (continuacion): Prediccion                        ##
##################################################################

##################################################################
## Paso 1: interpolacion de datos faltantes                     ##
##################################################################

## Datos y modelo de la Clase 19: malleco con un MA(3)
X <- LSTS::malleco
mod1 <- forecast::Arima(X - mean(X), order = c(0,0,3), include.mean = FALSE)

## Predecir no es solo hacia adelante: un faltante INTERMEDIO se proyecta
## sobre pasado Y futuro. Para un AR(1):
##   y.hat[s] = phi/(1+phi^2) (y[s-1] + y[s+1]),  ECM = sigma2/(1+phi^2)
## (menor que el sigma2 de predecir solo con el pasado)

Y <- X - mean(X)          ## malleco centrada
s <- 150                  ## "perdemos" una observacion intermedia
Y.na <- Y; Y.na[s] <- NA

## arima() acepta NA: el FILTRO de Kalman (solo hacia adelante) omite ese
## instante al evaluar la verosimilitud
fit <- arima(Y.na, order = c(1,0,0), include.mean = FALSE)
phi <- fit$coef[1]

## (i) la formula bilateral del AR(1)
imp1 <- (phi/(1 + phi^2))*(Y[s-1] + Y[s+1])
## (ii) el SUAVIZADOR de Kalman: agrega la pasada hacia atras y entrega
##      E[y_s | toda la muestra] -- la proyeccion sobre pasado y futuro
ks <- KalmanSmooth(Y.na, fit$model)
imp2 <- drop(ks$smooth %*% fit$model$Z)[s]
## (iii) solo filtro (solo pasado): phi * y[s-1]

round(cbind(formula = imp1, smoother = imp2, solo.pasado = phi*Y[s-1],
            valor.real = Y[s]) + mean(X), 4)
## formula y smoother coinciden EXACTO; el que solo usa el pasado queda
## mas lejos: la informacion del futuro reduce el ECM en 1/(1+phi^2)

## grafico zoom con la imputacion
par(bty = "n", las = 1)
plot(c(time(X))[(s-15):(s+15)], X[(s-15):(s+15)], type = "o", pch = 20,
     col = "gray40", xlab = "", ylab = "")
points(c(time(X))[s], imp2 + mean(X), pch = 19, col = "red", cex = 1.3)
points(c(time(X))[s], X[s], pch = 1, col = "blue", cex = 1.6)
legend("topleft", legend = c("imputacion (smoother)", "valor real"),
       col = c("red", "blue"), pch = c(19, 1), bty = "n")

##################################################################
## Paso 2: bandas de prediccion                                 ##
##################################################################

pre <- forecast::forecast(mod1, h = 10, level = 0.95)
## devolvemos la media al pronostico (el modelo se ajusto a X - mean(X))
pre$x <- pre$x + mean(X)
pre$mean <- pre$mean + mean(X)
pre$upper <- pre$upper + mean(X)
pre$lower <- pre$lower + mean(X)
par(bty = "n", las = 1)
plot(pre, main = "", xlim = c(1900, 2010))
## el pronostico converge a la media y la banda se abre hasta
## +- 1.96 sqrt(gamma(0)): la figura de las slides

## verificar: el se implicito en la banda es el de predict()
M <- cbind(banda = (pre$upper - pre$mean)/qnorm(0.975),
           predict = predict(mod1, n.ahead = 10)$se)
round(M, 5)

##################################################################
## Paso 3: bandas asimetricas via transformacion log            ##
##################################################################

## malleco es positiva: modelar en escala log y volver con exp entrega
## bandas asimetricas (mas anchas hacia arriba), lo natural para series
## positivas. En forecast::Arima basta lambda = 0 (Box-Cox log)
mod2 <- forecast::Arima(X, lambda = 0, order = c(0,0,3))
pre <- forecast::forecast(mod2, h = 100, level = 0.95)
plot(pre, main = "", bty = "n", las = 1, xlim = c(1200, 2000))

## la misma banda construida a mano: IC simetrico en escala log,
## exp() en los extremos --> asimetria en la escala original
mod3 <- forecast::Arima(log(X), order = c(0,0,3))
pre3 <- forecast::forecast(mod3, h = 100, level = 0.95)
LI <- exp(pre3$lower); LS <- exp(pre3$upper); pm <- exp(pre3$mean)
lines(pm, col = "red")
lines(LI, col = "red", lty = 2)
lines(LS, col = "red", lty = 2)
## coinciden: lambda = 0 es exactamente modelar log(X) y deshacer con exp

##################################################################
## Paso 4: correccion por metodo delta --> banda simetrica      ##
##################################################################

## Si W = log(X), el metodo delta (Taylor de 1er orden de g(W) = e^W)
## entrega Var(e^W) ~ e^{2 mu_W} Var(W). Con esto la banda en la escala
## original queda SIMETRICA en torno a exp(y.hat):
##   exp(y.hat(h)) +- z * exp(y.hat(h)) * se.log(h)

se.log <- (pre3$upper - pre3$mean)/qnorm(0.975)  ## se en escala log
pm <- exp(pre3$mean)
LI.d <- pm - qnorm(0.975)*sqrt(se.log^2 * pm^2)
LS.d <- pm + qnorm(0.975)*sqrt(se.log^2 * pm^2)

par(bty = "n", las = 1)
plot(pre, main = "", xlim = c(1200, 2000))   ## banda exacta (asimetrica)
lines(LI.d, col = "red")
lines(LS.d, col = "red")
legend("topleft", legend = c("exp de los extremos (asimetrica)",
                             "metodo delta (simetrica)"),
       col = c("blue", "red"), lty = 1, bty = "n")
## para horizontes cortos (se.log pequeno) ambas practicamente coinciden;
## a horizontes largos la delta se abre igual hacia ambos lados e incluso
## puede cruzar el cero -- el precio de forzar la simetria en una serie
## positiva. La asimetrica respeta el soporte; la delta es mas simple de
## reportar

## Moraleja:
## - Interpolacion: un faltante intermedio se proyecta sobre pasado y
##   futuro (suavizador de Kalman); el filtro solo mira hacia adelante.
## - Bandas: y.hat(h) +- z sqrt(sigma2(h)); con series positivas, modelar
##   en log (lambda = 0) y volver con exp --> bandas asimetricas, o
##   corregir por metodo delta si se quiere una banda simetrica.
