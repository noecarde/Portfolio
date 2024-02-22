
# MÉTODO DE FACTORIZACIÓN DE LA MATRIZ DE COVARIANZAS

# ALGORITMO:

# 1. OBTENER LA FACTORIZACIÓN DE CHOLESKY: Σ = L * L^t 
# 2. SIMULAR Z=(Z1,Z2,...,Zd) i.i.d N(0,1)
# 3. HACER X = mu + LZ
# 4. REPETIR LOS PASOS 2 Y 3 LAS VECES NECESARIAS

# CÓDIGO :

# Función para simular una normal multivariante:

# Define una función llamada simulacion_normal_multivariante que toma 
# tres argumentos: n (número de observaciones a simular), mu (vector de medias) y Sigma (matriz de covarianzas).

# Calcula la factorización de Cholesky de la matriz de covarianza Sigma. 
# La factorización de Cholesky es una descomposición de una matriz definida 
# positiva en el producto de una matriz triangular inferior y su traspuesta.
# chol en R nos da la matriz triangular superior de sigma asi que hacemos su traspuesta
# para conseguir la triangular inferior.

# Aquí,rnorm genera una cantidad n*length(mu) números aleatorios de una distribución normal estándar y matrix crea una matriz 
# con estos valores, donde las filas respresentan el numero de variables aleatorias 
# y las columnas el numero de observaciones

# X <- mu + L %*% Z: Para obtener la matriz X de observaciones de 
# la distribución normal multivariante sumamos el vector de medias mu a la matriz LZ 
# (siendo L la matriz triangular inferior de cholesky) y (Z la matriz de numeros aleatorios de una normal estandar)


simulacion_normal_multivariante <- function(n, mu, Sigma) {
  
  L <- t(chol(Sigma))
  
  Z <- matrix(rnorm(n * length(mu)), nrow = length(mu), ncol = n)
  
  X <- mu + L %*% Z
  
  return(X)
}

# EJEMPLO: 

# Tomamos como ejemplos:

# X1 = N(1,2) y X2 = N(2,3) donde cov(X1,X2) = 1 

# Función de densidad de probabilidad para la normal

mu <- c(1,2)  # Vector de medias
Sigma <- matrix(c(2, 0.014, 0.014, 3), nrow = 2, ncol = 2)  # Matriz de covarianzas

# Establecemos semilla para reproducibilidad

set.seed(100)  

# Haciendo uso de la funcion `simulacion_normal_multivariante` obtenemos la matriz `muestras_normal_multivariante`
# donde su primera fila es la simulación de X1 y su segunda la simulación de X2:

t <- proc.time()
muestras_normal_multivariante_mil <- simulacion_normal_multivariante(1000, mu, Sigma)
proc.time()-t
t <- proc.time()
muestras_normal_multivariante_10mil <- simulacion_normal_multivariante(10000, mu, Sigma)
proc.time()-t
t <- proc.time()
muestras_normal_multivariante_100mil <- simulacion_normal_multivariante(100000, mu, Sigma)
proc.time()-t
t <- proc.time()
muestras_normal_multivariante_millon <- simulacion_normal_multivariante(1000000, mu, Sigma)
proc.time()-t
t <- proc.time()
muestras_normal_multivariante_10millon <- simulacion_normal_multivariante(10000000, mu, Sigma)
proc.time()-t

# Visualizamos la matriz de simulación de ambas v.a. X1 y X2


print(muestras_normal_multivariante)

# Comprobamos gráficamente el ajuste de la densidad original de X1 y X2 a las simulaciones:

par(mfrow = c(1, 2))

hist(muestras_normal_multivariante_10mil[1,], col = "lightblue", main = "Histograma de N(1,2)", xlab = "X1",prob = TRUE, ylim=c(0,0.3))
curve(dnorm(x, mean = 1, sd = sqrt(2)), col = "darkblue", add = TRUE)
hist(muestras_normal_multivariante_10mil[2,], col = "lightpink", main = "Histograma de N(2,3)", xlab = "X2",prob = TRUE,ylim=c(0,0.25))
curve(dnorm(x, mean = 2, sd = sqrt(3)), col = "deeppink", add = TRUE)

par(mfrow = c(1, 1))

par(mfrow = c(1, 3))

hist(muestras_normal_multivariante_mil[2,], col = "lightpink", main = "1000 simulaciones", xlab = "X2",prob = TRUE,ylim=c(0,0.25))
curve(dnorm(x, mean = 2, sd = sqrt(3)), col = "deeppink", add = TRUE)
hist(muestras_normal_multivariante_10mil[2,], col = "lightpink", main = "10.000 simulaciones", xlab = "X2",prob = TRUE,ylim=c(0,0.25))
curve(dnorm(x, mean = 2, sd = sqrt(3)), col = "deeppink", add = TRUE)
hist(muestras_normal_multivariante_100mil[2,], col = "lightpink", main = "100.000 simulaciones", xlab = "X2",prob = TRUE,ylim=c(0,0.25))
curve(dnorm(x, mean = 2, sd = sqrt(3)), col = "deeppink", add = TRUE)

par(mfrow = c(1, 2))

hist(muestras_normal_multivariante_millon[2,], col = "lightpink", main = "1.000.000 de simulaciones", xlab = "X2",prob = TRUE,ylim=c(0,0.25))
curve(dnorm(x, mean = 2, sd = sqrt(3)), col = "deeppink", add = TRUE)
hist(muestras_normal_multivariante_10millon[2,], col = "lightpink", main = "10.000.000 de simulaciones", xlab = "X2",prob = TRUE,ylim=c(0,0.25))
curve(dnorm(x, mean = 2, sd = sqrt(3)), col = "deeppink", add = TRUE)

par(mfrow = c(1, 1))

# Instalar e cargar el paquete e1071
install.packages("e1071")
library(e1071)

# Calcular la asimetría
asimetria <- skewness(muestras_normal_multivariante_10millon[1,])

# Calcular la curtosis
curtosis <- kurtosis(muestras_normal_multivariante_10millon[1,])

# Imprimir resultados
print(paste("Asimetría:", asimetria))
print(paste("Curtosis:", curtosis))

# Calcular la asimetría
asimetria <- skewness(muestras_normal_multivariante_10millon[2,])

# Calcular la curtosis
curtosis <- kurtosis(muestras_normal_multivariante_10millon[2,])

# Imprimir resultados
print(paste("Asimetría:", asimetria))
print(paste("Curtosis:", curtosis))



library("MVN")

NM<-mvn(t(muestras_normal_multivariante_mil))
NM$Descriptives
NM<-mvn(t(muestras_normal_multivariante_10mil))
NM$Descriptives
NM<-mvn(t(muestras_normal_multivariante_100mil))
NM$skew
NM<-mvn(t(muestras_normal_multivariante_millon))
NM$Descriptives
NM<-mvn(t(muestras_normal_multivariante_10millon))
NM$Descriptives

install.packages("plotly")
library(plotly)

# Generar datos de la distribución normal multivariante
set.seed(1)
muestras_normal_multivariante <- simulacion_normal_multivariante(n_simulaciones, mu, Sigma)

# Crear histograma 3D
plot_ly(x = muestras_normal_multivariante_10mil[1,], y = muestras_normal_multivariante_10mil[2,], type = "histogram2d") %>%
  layout(scene = list(xaxis = list(title = "X1"),
                      yaxis = list(title = "X2"),
                      zaxis = list(title = "Frecuencia")))

plot(muestras_normal_multivariante_10mil[1,], muestras_normal_multivariante_10mil[2,],
     main = "Gráfico de dispersión",
     xlab = "X1", ylab = "X2", pch = 16, col = "blue")

# Instala y carga la biblioteca ggplot2 si aún no lo has hecho
# install.packages("ggplot2")
library(ggplot2)

# Supongamos que 'x' es tu vector de datos
x <- rnorm(1000)  # Ejemplo: datos generados aleatoriamente

# Crea un gráfico de densidad
ggplot(data = data.frame(muestras_normal_multivariante[1,]), aes(muestras_normal_multivariante[1,])) +
  geom_density(fill = "skyblue", color = "darkblue", alpha = 0.7) +
  labs(title = "Gráfico de Densidad", x = "Variable X", y = "Densidad")


mu_X1 <- 1
sigma_X1 <- sqrt(2)

mu_X2 <- 2
sigma_X2 <- sqrt(3)

set.seed(100)

muestra_X1 <- rnorm(10000, mean = mu_X1, sd = sigma_X1)
muestra_X2 <- rnorm(10000, mean = mu_X2, sd = sigma_X2)

covarianza <- cov(muestra_X1, muestra_X2)

correlacion <- covarianza/(sigma_X1*sigma_X2)

print(covarianza)
print(correlacion)

# MÉTODO DE LA TABLA GUÍA

# ALGORITMO:

# 1. INICIALIZACIÓN:
# Hacer F1 = p1
# Y desde i = 2 hasta n hacer Fi = Fi − 1 + pi

# 2. CÁLCULO DE LA TABLA GUÍA:
# Hacer g1 = 1, e i = 1
# Desde j = 2 hasta m hacer:
#   a) Mientras (j − 1)/m > Fi hacer i =  i + 1
#   b) gj = i

# 3. SIMULACIÓN MEDIANTE TABLA GUÍA:
# Generar U ∼ U( 0,  1)
# Hacer j = [mU] + 1
# Hacer i = gj
# Mientras U > Fi hacer i = i + 1
# Devolver X = xi

# CÓDIGO :

# El método de la tabla guía ofrece una manera eficiente de generar valores 
# aleatorios de una distribución discreta utilizando una tabla precalculada. 
# El código realiza el cálculo de la tabla guía y luego utiliza esta tabla para asignar 
# valores simulados de manera rápida y eficiente.


simulacion_tabla_guia <- function(x, prob, m, n, as.factor = FALSE) {
  
  # Inicializar tabla y FD
  
  Fx <- cumsum(prob)
  g <- rep(1,m)
  i <- 1
  for(j in 2:m) {
    while (Fx[i] < (j-1)/m) i <- i + 1
    g[j] <- i
  }
  
  # Generar valores
  
  X <- numeric(n)
  U <- runif(n)
  for(j in 1:n) {
    i <- i0 <- g[floor(U[j] * m) + 1]
    while (Fx[i] < U[j]) i <- i + 1
    X[j] <- x[i]
  }
  return(X)
}

# EJEMPLO: 

set.seed(123)

# Parámetros de la distribución binomial
n <- 10  # Número de ensayos
p <- 0.4  # Probabilidad de éxito

# Valores posibles de la variable aleatoria
x <- 0:n

# Aplicar la función de simulación
resultado_simulacion <- simulacion_tabla_guia(x, prob = dbinom(x, n, p), m = 1000, n = 1000, as.factor = FALSE)

# Imprimir los primeros 10 valores simulados
print(table(resultado_simulacion))

# Realizar una inspección cuantitativa de los resultados obtenidos en una simulación es crucial para 
# evaluar la validez y precisión de la simulación. Para ello vamos a hacer 
# Es importante evaluar la validez y precisión de la simulación. Para ello vamos a hacer 

# Comparación con la Distribución Teórica:
plot(ecdf(resultado_simulacion), main = "Distribución empírica simulación de B(10,0.4)")
curve(pbinom(x,10,0.4),add= TRUE, col= "violetred")

# Valores teóricos
teoricos <- dbinom(x, n, p)

# Comparación de las medias
cat("Media teórica:", sum(x * teoricos), "\n")
cat("Media simulada:", mean(resultado_simulacion), "\n")

# Comparación de las desviaciones estándar
cat("Desviación Estándar teórica:", sqrt(sum((x - sum(x * teoricos))^2 * teoricos)), "\n")
cat("Desviación Estándar simulada:", sd(resultado_simulacion), "\n")

# Comparación de las varianzas
cat("Varianza teórica:", sum(n * p * (1 - p)), "\n")
cat("Varianza simulada:",var(resultado_simulacion), "\n")

# MÉTODO DE ALIAS

simulacion_alias <- function(x, prob, n, as.factor = FALSE) {
  # Inicializar tablas
  a <- numeric(length(x))
  q <- prob*length(x)
  low <- q < 1
  high <- which(!low)
  low <- which(low)
  while (length(high) && length(low)) {
    l <- low[1]
    h <- high[1]
    a[l] <- h
    q[h] <- q[h] - (1 - q[l])
    if (q[h] < 1) {
      high <- high[-1]
      low[1] <- h
    } else low <- low[-1]
  } # while
  # Generar valores
  V <- runif(n)
  i <- floor(runif(n)*length(x)) + 1
  X <- x[ ifelse( V < q[i], i, a[i]) ]
  if(as.factor) X <- factor(X, levels = x)
  return(X)
}

# Definir los parámetros de la distribución binomial
n <- 10  # Número de ensayos
p <- 0.4  # Probabilidad de éxito

# Definir los posibles valores de la distribución binomial
x <- 0:n

# Calcular las probabilidades
prob <- dbinom(x, size = n, prob = p)

# Utilizar la función rpmf.alias para simular la distribución binomial
set.seed(123)  # Establecer semilla para reproducibilidad
simul_val <- simulacion_alias(x, prob, n = 1000)

# Comparar resultados con la distribución teórica
table(simul_val)  # Frecuencias relativas simuladas
# Probabilidades teóricas
plot(ecdf(simul_val), main = "Distribución empírica simulación de B(10,0.4)")
curve(pbinom(x,10,0.4),add= TRUE, col= "green")

#teoricos
prob_teoricos <- dbinom(x, size = n, prob = p)

# Comparación de valores teóricos y simulados
cat("Media teórica:", sum(x * prob_teoricos), "\n")
cat("Media simulada:", mean(simul_val), "\n")

cat("Desviación Estándar teórica:", sqrt(sum((x - sum(x * prob_teoricos))^2 * prob_teoricos)), "\n")
cat("Desviación Estándar simulada:", sd(simul_val), "\n")

cat("Varianza teórica:", sum(n * p * (1 - p)), "\n")
cat("Varianza simulada:", var(simul_val), "\n")
