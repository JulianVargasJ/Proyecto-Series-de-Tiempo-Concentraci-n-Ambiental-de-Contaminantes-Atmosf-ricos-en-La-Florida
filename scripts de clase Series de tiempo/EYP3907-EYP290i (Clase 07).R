#############
## Clase 6 ##
#############

## Requiere: LSTS

## Modelos Parametricos vs No Parametricos ##

##################################
## Ejemplo 1: Modelo PARAMETRICO ##
##################################

## AR(1) --> Y[t] = phi * Y[t-1] + e[t], {e[t]} ~ RB(0, sigma^2), |phi| < 1
## Todo el proceso queda descrito por un parametro de dimension FINITA:
## theta = (phi, sigma)

X <- LSTS::malleco
par(bty = "n", las = 1)
plot(X)

ACF <- acf(X, lag.max = 10, ylim = c(-1,+1))

## Estimacion por metodo de momentos: en un AR(1), rho(1) = phi
phi <- ACF$acf[2,,1]
h <- 0:10
## ACF teorica del AR(1): rho(h) = phi^h --> con UN parametro reproducimos
## (aproximadamente) toda la estructura de dependencia
lines(phi^h ~ h, type = "p", pch = 20, lwd = 2, col = "red")

## El segundo parametro: Var(Y) = sigma^2/(1-phi^2) --> despejamos sigma
sigma <- sqrt(var(X) * (1 - phi^2))
c(phi = phi, sigma = sigma)

######################################
## Ejemplo 2: Modelo NO PARAMETRICO ##
######################################

## Aqui no se postula una forma funcional: la "curva" vive en un espacio de
## funciones y su complejidad efectiva la eligen los datos

X <- co2
par(bty = "n", las = 1)
plot(X)

## Spline de suavizamiento (spar fijado a mano)
mod <- smooth.spline(x = time(X), y = X, spar = 0.8)
lines(mod, col = "red", lwd = 2)
mod$df   ## grados de libertad EFECTIVOS: la "cantidad de parametros" que uso

## spar elegido automaticamente por validacion cruzada generalizada (GCV)
mod.gcv <- smooth.spline(x = time(X), y = X)
lines(mod.gcv, col = "blue", lwd = 2)
c(spar = mod.gcv$spar, df = mod.gcv$df)
## Con GCV el spline persigue tambien la estacionalidad (df ~ 115!);
## con spar = 0.8 solo captura la macrotendencia (df ~ 10)

## Otro suavizador no parametrico: regresion kernel
lines(ksmooth(time(X), X, kernel = "normal", bandwidth = 2), col = "orange", lwd = 2)
legend("topleft", legend = c("spline spar = 0.8 (df ~ 10)", "spline GCV (df ~ 115)", "kernel normal (bw = 2)"),
       col = c("red", "blue", "orange"), lty = 1, lwd = 2, bty = "n")

## Moraleja:
## - Parametrico:    pocos parametros con interpretacion (phi, sigma); extrapola.
## - No parametrico: flexible, la dimension "efectiva" crece con los datos;
##                   describe muy bien dentro de la muestra, extrapola mal.
## - El suavizador captura la TENDENCIA, pero la dependencia serial de los
##   residuos sigue ahi --> eso lo modelaremos con procesos lineales (ARMA)
