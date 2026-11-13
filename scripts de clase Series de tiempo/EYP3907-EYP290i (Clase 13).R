##############
## Clase 13 ##
##############

## Capitulo 3: Analisis Espectral
## - Densidad espectral teorica: LSTS::spectral.density
## - Frecuencias de Fourier y la idea del video (enrollar la senal)
## - Periodograma "a mano" y LSTS::periodogram
## - Insesgado pero NO consistente (simulacion)
## - Suavizado: LSTS::smooth.periodogram (splines, como la Clase 6)

## Requiere: LSTS

########################################################
## Paso 1: densidades espectrales teoricas            ##
########################################################

## f(w) reparte la varianza del proceso entre frecuencias w in [0, pi]:
## gamma(0) = integral de f  -->  la forma de f es la "huella" del modelo

lambda <- seq(0.01, pi, length = 400)

par(bty = "n", las = 1, mfrow = c(2,2))
## Ruido blanco: plana (todas las frecuencias aportan lo mismo)
plot(lambda, LSTS::spectral.density(lambda = lambda), type = "l",
     lwd = 2, col = "blue", xlab = expression(omega), ylab = expression(f(omega)),
     main = "Ruido blanco")
## MA(1): suave, con peso en bajas frecuencias si theta > 0
plot(lambda, LSTS::spectral.density(ma = 0.9, lambda = lambda), type = "l",
     lwd = 2, col = "blue", xlab = expression(omega), ylab = expression(f(omega)),
     main = expression(MA(1)*": "*theta == 0.9))
## AR(1) con phi > 0: peso en w = 0 (persistencia, ciclos largos)
plot(lambda, LSTS::spectral.density(ar = 0.7, lambda = lambda), type = "l",
     lwd = 2, col = "blue", xlab = expression(omega), ylab = expression(f(omega)),
     main = expression(AR(1)*": "*phi == +0.7))
## AR(1) con phi < 0: peso en w = pi (serie que alterna signo)
plot(lambda, LSTS::spectral.density(ar = -0.7, lambda = lambda), type = "l",
     lwd = 2, col = "blue", xlab = expression(omega), ylab = expression(f(omega)),
     main = expression(AR(1)*": "*phi == -0.7))

## ARFIMA: la misma funcion con el argumento d --> diverge en w = 0
## (la caracterizacion espectral de la larga memoria del Cap. 2)
par(bty = "n", las = 1, mfrow = c(1,1))
plot(lambda, LSTS::spectral.density(d = 0.4, lambda = lambda), type = "l",
     lwd = 2, col = "blue", ylim = c(0, 4),
     xlab = expression(omega), ylab = expression(f(omega)))
lines(lambda, LSTS::spectral.density(d = 0.2, lambda = lambda), lwd = 2, col = "orange")
lines(lambda, LSTS::spectral.density(ar = 0.7, lambda = lambda), lwd = 2, col = "red")
legend("topright", legend = c("FN(d = 0.4)", "FN(d = 0.2)", "AR(1)"),
       col = c("blue","orange","red"), lwd = 2, bty = "n")

##########################################################
## Paso 2: enrollar la senal (la idea del video)        ##
##########################################################

## Serie mensual con ciclo anual: n = 120 (10 anos)
n <- 120
t <- 1:n
set.seed(3907)
y <- cos(2*pi*t/12) + rnorm(n, 0, 0.3)
yc <- y - mean(y)

## Multiplicar por exp(-i w t) "enrolla" la serie alrededor del origen;
## el centro de masa m(w) se descentra SOLO si la serie cicla en w
par(bty = "n", las = 1, mfrow = c(1,2), pty = "s")
for(j in c(3, 10)){
  w <- 2*pi*j/n
  z <- yc * exp(-1i*w*t)
  m <- mean(z)
  plot(Re(z), Im(z), pch = 20, col = "gray",
       xlab = "parte real", ylab = "parte imaginaria",
       main = paste0("j = ", j, "  (periodo ", round(n/j,1), " meses)"))
  abline(h = 0, v = 0, lty = 3)
  points(Re(m), Im(m), pch = 19, col = "red", cex = 1.5)
  arrows(0, 0, Re(m), Im(m), col = "red", lwd = 2, length = 0.08)
}
## j = 3: nube centrada (la serie no cicla ahi); j = 10 (periodo 12): el
## ciclo anual empuja todo hacia un lado

##########################################################
## Paso 3: el periodograma "a mano"                     ##
##########################################################

## I(w_j) = (1/2 pi n) [ (sum yc cos)^2 + (sum yc sin)^2 ]
##        = (n/2 pi) |m(w_j)|^2  --> distancia^2 del centro de masa
I <- c()
J <- 1:(n/2)
for(j in J){
  w <- 2*pi*j/n
  I[j] <- (sum(yc*cos(w*t))^2 + sum(yc*sin(w*t))^2)/(2*pi*n)
}

par(bty = "n", las = 1, mfrow = c(1,1), pty = "m")
plot(J, I, type = "h", xlab = "j", ylab = expression(I(omega[j])))
points(10, I[10], pch = 19, col = "red")
## el peak esta en j = 10: 10 ciclos en 120 meses = periodo 12 = ciclo anual

## Comparar con la funcion de la libreria del curso: mismo resultado
P <- LSTS::periodogram(y, plot = FALSE)
cbind(a.mano = I[8:12], LSTS = P$periodogram[8:12])
## identicos (en j = 10 ambos dan 4.4349): misma normalizacion

##################################################################
## Paso 4: insesgado pero NO consistente (simulacion)           ##
##################################################################

## Para RB gaussiano en una frecuencia de Fourier: I(w_j) = f(w_j) chi2_2/2
## --> E[I] = f, pero Var[I] = f^2, SIN importar n

f0 <- 1/(2*pi)   ## f teorica del RB(0,1)
for(n.sim in c(100, 1000, 10000)){
  j <- n.sim/4; w <- 2*pi*j/n.sim; tt <- 1:n.sim
  I.sim <- replicate(1000, {
    e <- rnorm(n.sim)
    abs(sum((e-mean(e))*exp(-1i*w*tt)))^2/(2*pi*n.sim)
  })
  cat("n =", format(n.sim, width = 6),
      " E[I] =", round(mean(I.sim), 4),
      " sd[I] =", round(sd(I.sim), 4),
      "   (f =", round(f0, 4), ")\n")
}
## la media esta siempre en f (insesgado), pero la desviacion estandar
## NO baja al crecer n: mas datos = mas frecuencias, no mejores estimaciones

##################################################################
## Paso 5: suavizado --> LSTS::smooth.periodogram y malleco     ##
##################################################################

X <- LSTS::malleco

## Periodograma crudo y suavizado (splines, con parametro spar, el mismo
## suavizador smooth.spline de la Clase 6). Ojo: las funciones de LSTS
## grafican en ggplot, asi que para SUPERPONER curvas trabajamos con los
## valores que devuelven ($lambda, $periodogram, $smooth.periodogram)
P0 <- LSTS::periodogram(X, plot = FALSE)
P1 <- LSTS::smooth.periodogram(X, spar = 0.7, plot = FALSE)

## Suavizado por PROMEDIO MOVIL (la formula de las slides):
## f-hat(w_j) = sum_{h=-m}^{m} W(h) I(w_{j+h}), con pesos uniformes
## W(h) = 1/(2m+1) --> promediar las 2m+1 ordenadas vecinas
m <- 10
W <- rep(1/(2*m + 1), 2*m + 1)
f.ma <- stats::filter(P0$periodogram, W, sides = 2)

## ¿Y como se compara con el modelo AR(1) que le venimos ajustando?
aux <- forecast::Arima(X - mean(X), order = c(1,0,0), include.mean = FALSE)
phi <- aux$coef[1]
s   <- sqrt(aux$sigma2)

par(bty = "n", las = 1, mfrow = c(1,1))
plot(P0$lambda, P0$periodogram, type = "h", col = "gray",
     xlab = expression(omega), ylab = "")
lines(P0$lambda, f.ma, col = "darkgreen", lwd = 2)
lines(P1$lambda, P1$smooth.periodogram, col = "blue", lwd = 2)
lines(P1$lambda, LSTS::spectral.density(ar = phi, sd = s, lambda = P1$lambda),
      col = "red", lwd = 2)
legend("topright",
       legend = c("periodograma", "promedio movil (slides, m = 10)",
                  "spline (LSTS, spar = 0.7)", "f del AR(1) ajustado"),
       col = c("gray", "darkgreen", "blue", "red"), lwd = c(1,2,2,2), bty = "n")
## Los dos suavizadores cuentan la misma historia y siguen a la f del
## AR(1) ajustado: el promedio movil (formula de las slides) es mas
## rugoso y local; el spline entrega una curva mas lisa. La varianza del
## promedio movil baja en un factor ~ 1/(2m+1) respecto del periodograma

## Moraleja:
## - La densidad espectral y la ACF son la MISMA informacion en dominios
##   distintos; la forma de f es la huella del modelo.
## - El periodograma es insesgado pero ruidoso (no consistente): siempre
##   suavizar antes de leerlo.
## - Frecuencia w_j = 2 pi j/n = "j ciclos completos en la muestra";
##   el peak identifica el periodo dominante (aqui, el ciclo anual).
