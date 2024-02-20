#NACIMIENTOS TOTALES

#Cargamos paquetes
install.packages('descomponer')
install.packages('nortest')
install.packages('tsoutliers')
install.packages('forecast')
library(tsoutliers)
library(forecast)
library(lmtest)
library(descomponer)
library(nortest)
library(ggplot2)

#Aplicaremos BOX-JENKINS
#Importamos los datos.
#Transformamos el dataset en serie temporal.
datos=ts(nacimientos_totales$Nacimientos, frequency=12)
plot(datos, col = "purple", main = "Gráfico de la Serie Temporal")
#Vemos una tendencia creciente y una componente estacional marcada.

#Estudiamos la existencia de outliers 
outliers.datos=tsoutliers::tso(datos,types="AO",maxit.iloop=10)
outliers.datos
#La serie no contiene outliers

#Realizamos la descomposición de la serie:
des=decompose(datos)
plot(des, col = "purple") 
#Viendo el gráfico de descomposición la serie no es estacionaria.

#Vamos a estudiar el número de diferenciaciones necesarias para que lo sea.

#Calculamos las diferenciaciones estacionales 
ns=nsdiffs(datos) #ESTO ES D EN EL MODELO
ns

datos_dif1= diff(datos,lag=12,ns)
plot(datos_dif1, col="deeppink")
#Calculamos las diferenciaciones simples
n=ndiffs(datos_dif1) #ESTO ES d EN EL MODELO
n
#No es necesaria ninguna dif simple
#La componente estacional ha sido eliminada 

periodograma(datos_dif1)
gperiodograma(datos_dif1)
#Vemos que los picos de densidad del periodograma no aparecen a intervalos ctes por lo que la estacionalidad que existia se ha eliminado.

#Estimamos P, p, Q y q
acf(datos, lag=48, main="ACF para datos", col="purple") #miramos multiplos de 12 entonces probaremos Q=0
pacf(datos, lag=48, main="PACF para datos", col="purple") #miramos multiplos de 12 entonces probaremos P=1
acf(datos_dif1, lag=24,col="pink") #Hay un pico significativo fuera de los intervalos de confianza, lo que nos sugiere probar q=1 
pacf(datos_dif1, lag=24, col="pink") #Hay un pico significativo fuera de los intervalos de confianza, lo que nos sugiere probar que p= 1
#Recordamos que estas estimaciones son orientativas  d=0 y D=1

#Buscamos el mejor modelo:

#SARIMA(1,0,1)x(0,1,0)
fitARIMA1=arima(datos, order=c(1,0,1), seasonal = list(order=c(0,1,0), period=12),method='ML') 
fitARIMA1
coeftest(fitARIMA1) 
confint(fitARIMA1)
#Los coeficientes son significativos porque el p-valor es < que 0,05

#SARIMA(1,0,1)X(1,1,0)      
fitARIMA2=arima(datos, order=c(1,0,1), seasonal = list(order=c(1,1,0), period=12),method='ML') 
fitARIMA2
coeftest(fitARIMA2) 
confint(fitARIMA2)
#Los coeficientes son significativos porque el p-valor es < que 0,05

#SARIMA(1,0,0) X (1,1,0)_12 
fitARIMA3 <- arima(datos, order = c(1, 0, 0), seasonal = list(order = c(1, 1, 0), period = 12), method = 'ML')
fitARIMA3
coeftest(fitARIMA3) 
confint(fitARIMA3)
#Los coeficientes son significativos porque el p-valor es < que 0,05


####################################################################################################
################################# BONDAD DE AJUSTE #################################################
####################################################################################################

#Elegimos empezar con el segundo modelo porque su aic es menor

######################### ANALISIS DE RESIDUOS MODELO SARIMA(1,0,1)X(1,1,0)_12################################

#Verificamos la INDEPENDENCIA de la distribución de los datos con Ljung-Box
residuos2 <- residuals(fitARIMA2)

# Realizar el test de Ljung-Box
ljung_box_test2 <- Box.test(residuos2, lag = 24, type = "Ljung-Box")
ljung_box_test2
#El modelo2 verifica la hipótesis de independecia pq pvalor >0.05

#HOMOCEDASTICIDAD
checkresiduals(fitARIMA2) 
#No hay patrones, realizaremos una transformación funcional más adelante por si pudiera mejorar la homocedasticidad ya presente.


#Estudiamos la NORMALIDAD de los residuos, para ello implementamos el test de Lilliefors (Kolmogorov-Smirnov) ya que tenemos una muestra de datos >50.
lillie.test(fitARIMA2$residuals)
#No existe normalidad de los residuos.
qqnorm(fitARIMA2$residuals) 
qqline(fitARIMA2$residuals, col="deeppink")

######################### ANALISIS DE RESIDUOS MODELO SARIMA(1,0,1)X(0,1,0)_12 ################################

#Vamos a analizar el primer modelo a ver si este cumple todas las hipótesis.
residuos1 <- residuals(fitARIMA1)

# Realizar el test de Ljung-Box
ljung_box_test1 <- Box.test(residuos1, lag = 24, type = "Ljung-Box")
ljung_box_test1
#Verifica la hipótesis de independecia porque pvalor >0.05

checkresiduals(fitARIMA1) 
#No hay patrones, realizaremos una transformación funcional más adelante por si pudiera mejorar la homocedasticidad ya presente.

#Estudiamos la NORMALIDAD de los residuos, para ello implementamos el test de Lilliefors (Kolmogorov-Smirnov) ya que tenemos una muestra de datos >50.

lillie.test(fitARIMA1$residuals)

#No existe normalidad entre los residuos porque pvalor<0.05
qqnorm(fitARIMA1$residuals) 
qqline(fitARIMA1$residuals, col="deeppink")
#Este modelo NO cumple todas las hipótesis.

######################### ANALISIS DE RESIDUOS MODELO SARIMA(1,0,0)X(1,1,0)_12 ################################
#Vamos a analizar el último modelo propuesto.

residuos3 <- residuals(fitARIMA3)

# Realizar el test de Ljung-Box
ljung_box_test3 <- Box.test(residuos3, lag = 24, type = "Ljung-Box")
ljung_box_test3
#El modelo verifica la hipótesis de independecia porque pvalor >0.05

checkresiduals(fitARIMA3) 
#No hay patrones, realizaremos una transformación funcional más adelante por si pudiera mejorar la homocedasticidad ya presente.

#Estudiamos la NORMALIDAD de los residuos, para ello implementamos el test de Lilliefors (Kolmogorov-Smirnov) ya que tenemos una muestra de datos >50.

lillie.test(fitARIMA3$residuals)

#Existe normalidad entre los residuos porque pvalor > 0.05
qqnorm(fitARIMA3$residuals) 
qqline(fitARIMA3$residuals, col="deeppink")

#Este modelo SÍ cumple todas las hipótesis.


####################################################################################################
######################### TRANSFORMACIÓN FUNCIONAL #################################################
####################################################################################################

lambda=forecast::BoxCox.lambda(datos)
lambda
datost=forecast::BoxCox(datos, lambda = lambda)
plot(datost, col="deeppink") 

#Numero de diferenciaciones estacionales
D2 <- nsdiffs(datost)
D2
#Numero de diferenciaciones simples
d2 <-ndiffs(D2)
d2

datos_dif2=diff(datost,lag=12,D2)
periodograma(datos_dif2)
gperiodograma(datos_dif2) 

#HOMOCEDASTICIDAD:
plot(datos, col='deeppink') 
plot(datost, col="purple") 
#Vemos que la transformación no aporta cambios significativos en nuestros datos.

#ELEGIMOS DE P, p, Q, q
acf(datost, lag.max=60, col="deeppink")
pacf(datost, lag.max=60,col="deeppink")
acf(datos_dif2, lag.max=60,col="pink")
pacf(datos_dif2, lag.max=60, col="pink")

#Modelo SARIMA(1,0,0) X (1,1,0)_12
fitARIMA4=arima(datost, order=c(1,0,0), seasonal = list(order=c(1,1,0), period=12),method='ML') 
fitARIMA4
coeftest(fitARIMA4) #sale bien

#Modelo SARIMA(1,0,0) X (0,1,1)_12
fitARIMA5=arima(datost, order=c(1,0,0), seasonal = list(order=c(0,1,1), period=12),method='ML') 
fitARIMA5
coeftest(fitARIMA5) 

######################### ANALISIS DE RESIDUOS POST-TRANSFORMACION ################################

#INDEPENDENCIA DESPUES DE LA TRANSFORMACION, vamos a analizar los dos modelos aunque parece que sería mejor elegir el modelo 5 ya que tiene un aic menor

#La independencia se mantiene para los dos modelos.
residuos4 <- residuals(fitARIMA4)
ljung_box_test4 <- Box.test(residuos4, lag = 24, type = "Ljung-Box")
ljung_box_test4
checkresiduals(fitARIMA4)
#El modelo4 verifica la hipótesis de independecia porque pvalor >0.05
residuos5 <- residuals(fitARIMA5)
ljung_box_test5 <- Box.test(residuos5, lag = 24, type = "Ljung-Box")
ljung_box_test5
checkresiduals(fitARIMA5)
#El modelo5 verifica la hipótesis de independecia porque pvalor >0.05

#NORMALIDAD DESPUES DE LA TRANSFORMACION
lillie.test(fitARIMA4$residuals)
#Los datos NO son normales
qqnorm(fitARIMA4$residuals) 
qqline(fitARIMA4$residuals, col="deeppink") 

lillie.test(fitARIMA5$residuals)
#Los datos NO son normales
qqnorm(fitARIMA5$residuals) 
qqline(fitARIMA5$residuals, col="deeppink") 

#La transformación funcional NO ha funcionado ya que ya teniamos homocedasticidad y además perdememos la normalidad.

#Por lo que concluimos que nuestro mejor modelo es fitarima3

#PREDICCIÓN AÑO 2005 ELIGIENDO EL MODELO INICIAL FITARIMA3
predict(fitARIMA3, n.ahead=12)
val_fut=forecast(fitARIMA3, h=12, level=c(95,99))
plot(val_fut, col="deeppink")
#Se aprecia un aumento, vamos a hacer un gráfico con cuadrícula para confirmarlo.
forecast <- forecast(fitARIMA3, h = 12)
data <- data.frame(
  time = c(time(datos), time(forecast$mean)),
  value = c(datos, forecast$mean),
  Leyenda = c(rep("Real", length(datos)), rep("Predicción", length(forecast$mean)))
)

ggplot(data, aes(x = time, y = value, color = Leyenda)) +
  geom_line() +
  labs(title = "Datos Reales vs Predicciones",
       x = "Tiempo",
       y = "Valor") +
  scale_color_manual(values = c("deeppink", "pink")) +
  theme_minimal()
#Contrastamos el aumento con información real y podemos afirmar que la predicción de aumento es correcta.
