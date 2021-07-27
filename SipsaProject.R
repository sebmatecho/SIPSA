###Cálculo del índice SIPSA###

setwd("G:/Mi unidad/U libertadores/2018.1/Especialización/Trabajo de grado/Archivos finales")
setwd("C:/Users/Seb/Google Drive/U libertadores/2018.1/Especialización/Trabajo de grado/SIPSA semanal")
library(AER)
library(readxl)
library(zoo)
library(vars)
library(car)
library(tseries)
library(strucchange)
library(MTS)
library(urca)
library(lmtest)
library(tsDyn)
library(systemfit)
library(TSA)
library(tseries)
library(lmtest)
library(forecast)
#IPC sale de: 
#https://www.dane.gov.co/index.php/estadisticas-por-tema/precios-y-costos/indice-de-precios-al-consumidor-ipc/grupos-ipc-2012
#Generación de archivo con los nombres de todos los archivos que se consideran

write.table(list.files("C:/Users/santander/Documents/Datos/SIPSA semanal"),"archivos.txt",
            quote = FALSE,sep = "\t",row.names = F,col.names = F)
setwd("G:/Mi unidad/U libertadores/2018.1/Especialización/Trabajo de grado/Archivos finales/ultimo")
#### formatos originales del SIPSA
#Tablas en las que se imprimen las medias y medianas
#write.table("medias_semana.txt",append = T,col.names=F)
write.table("medianas_semana.txt",append = T,col.names=F)
##Buscar los NA. Estos puntos indican el cambio de producto dentro de la base
archivos<-read.table("archivos_ordenado.txt",header = F, stringsAsFactors = FALSE)
for(i in 1:nrow(archivos)){
  print(i)
  file<-read_excel(archivos[i,],col_types = "numeric")
  nas<-which(is.na(file[,4]))
  medias<-c()
  medianas<-c()
  fact<-c()
  j<-1 
  #Hacer el promedio para cada producto
  for(i in 2:(length(nas)-1)){
    if(nas[i+1]-nas[i]>1){
      medias[i]<-sum(file[(nas[i]+1):(nas[i+1]-1),4])/nrow(file[(nas[i]+1):(nas[i+1]-1),4])
    }
    else{
      file[nas[i]+1,4]
    }
    if(nas[i+2] == nas[i+1]+1 & nas[i+1]+1==nas[i]+2)# && nas[i+1]== nas[i+2])
    {
      fact[j]<-i
      j<-j+1
    }
  }
  
  
  fact<-c(0,fact)
  #Hacer el precio promedio para cada tipo de alimento (8 en total)
  medias_semana<-c()
  medianas_semana<-c()
  for(k in 1:7){
  #  medias_semana[k]<-mean(medias[(fact[k]+2):(fact[k+1]-1)])
    medianas_semana[k]<-median(medias[(fact[k]+2):(fact[k+1]-1)])
  }
 # medias_semana[8]<-mean(medias[(fact[8]+2):length(medias)])
  medianas_semana[8]<-median(medias[(fact[8]+2):length(medias)])
  
  #cat(t(round(medias_semana,3)),file="medias_semana.txt",append=TRUE,fill=T,sep="\t")
  cat(t(round(medianas_semana,3)),file="medianas_semana.txt",append=TRUE,fill=T,sep="\t")
  
}

#####Archivos que tienen codificación diferente
nombres<-c("Verduras_Hortalizas","Fruta","Tubérculos",
           "Granos_cereales","Huevos_Lácteos","Carnes",
           "Pescados","Procesados")
nombres2<- c("1.1","1.2","1.3","1.4","1.5","1.6","1.7","1.8") 
#write.table("medianas_semana2formato.txt",append = T,col.names=F)
archivos<-read.table("archivos_ordenado2.txt",header = FALSE, stringsAsFactors = FALSE)
todo<- matrix(NA,nrow=nrow(archivos), ncol= 8)

#cálculo de todo 
## ciclo for para los archivos
for(k in 1:nrow(archivos)){
  # este archivo tiene un formato diferente
  if(k==1){
    ## ciclo for para cada pestaña (categoría alimentaria) dentro del archivo problema 
    for(j in 1:length(nombres)){
      base<- read_excel(archivos[k,], sheet = nombres[j])
      facs<-levels(as.factor(base$Producto))
      medias<- c()
      ##ciclo para los precios de cada producto
      for(i in 1:length(facs)){
        medias[i]<-mean(base$`Precio medio`[which(base$Producto==facs[i])])
      }#ciclo de producto 
      todo[k,j]<- median(medias)
    }#ciclo de cada factor de alimento
  }else{
    #Para los otros tipos de archivo, factor de alimento
    for(h in 1:length(nombres2)){
      base<- read_excel(archivos[k,], sheet = nombres2[h])
      sadasd<- tbl_df(base[-c(1:7),1]) %>% slice() %>%  unlist(., use.names=FALSE)
      facs2<-levels(as.factor(sadasd))
      ##Ciclo para los precios de cada producto 
      for(l in 1:length(facs2)){
        medias[l]<- median(data.matrix(base[which(base[,1]==facs2[l]),5]))
      }## ciclo medio de cada tipo de producto. 
      todo[k,h]<- mean(medias)
    }## fin del tipo de archivos que no están mal 
    cat(t(round(medianas_semana,3)),file="medianas_semana.txt",append=TRUE,fill=T,sep="\t")
    # cat(t(round(todo,3)),file="medianas_semana2formato.txt",append=TRUE,fill=T,sep="\t")
  } #for de archivos
}



setwd("C:/Users/Seb/Google Drive/U libertadores/2018.1/Especialización/Trabajo de grado/Archivos finales")
setwd("G:/Mi unidad/U libertadores/2018.1/Especialización/Trabajo de grado/Archivos finales")

sipsaxreg<-read.table("sipsaxreg.txt")
sipsaxreg2<-read.table("sipsaxreg2.txt")
sipsaxregtodo<-read.table("sipsaxregtodo.txt")
sipsaxregmedias<-read.table("sipsaxregmedias.txt",header=F)
sipsaxregmedias2<-read.table("sipsaxregmedias2.txt",header=F)
base<- read.table("basesipsa.txt",header=T)
inf<-ts(base[,1],start=c(2013,2),freq=12)
sipsa<-ts(base[,2],start=c(2013,2),freq=12)
sipsa_media<-ts(base[,3],start=c(2013,2),freq=12)
inf1<-window(inf,start=c(2013,2),end=c(2017,8))
inf2<-window(inf,start=c(2017,9))
sipsa1<-window(sipsa,start=c(2013,2),end=c(2017,8))
sipsa2<-window(sipsa,start=c(2017,9))
sipsa_medias1<-window(sipsa,start=c(2013,2),end=c(2017,8))
sipsa_medias2<-window(sipsa,start=c(2017,9))
par(mfrow=c(1,1))
postscript("ambas.eps",width = 7, height = 5)
plot.zoo(cbind(zoo(inf), zoo(sipsa)),#,zoo(sipsa_media)), 
         plot.type = "single", 
         col = c("red", "blue","green"),ylim=c(-4,6), 
         main="Contraste entre índices \n SIPSA e Inflación de alimentos",
         xlab="Año",ylab="Variación")
legend("topleft", inset=0.05, legend=c("IPC alimentos", "Índice SIPSA"),#,"SIPSA medias"),
       col=c("red", "blue","green"), lty=1:3, cex=0.8)
dev.off()
cor(inf,sipsa)
cor(inf,sipsa_media)



par(mfrow=c(1,3))
plot(inf1)
acf(inf1,lag.max=48)
pacf(inf1,lag.max=48)
plot(diff(inf1,12))
acf(diff(inf1,12),lag.max=48)
pacf(diff(inf,12),lag.max=48)
postscript("diffdiff.eps",width = 9, height = 3)
par(mfrow=c(1,3))
plot(diff(diff(inf1,12)),main="Datos diferenciados \n en sus componentes \n estacional y no estacional")
abline(h=2*sd(diff(diff(inf1,12))),col=2,lwd=1)
abline(h=-2*sd(diff(diff(inf1,12))),col=2,lwd=1)
acf(diff(diff(inf1,12)),lag.max=48,main="ACF")
pacf(diff(diff(inf1,12)),lag.max=48,main="PACF")
dev.off()
pp.test(inf1, lshort = F)
pp.test(diff(diff(inf1,12)), lshort = F)
pp.test(diff(inf1), lshort = F)
options(digits=5)

adf.test(diff(diff(inf1,12)))

###Modelo basado en SIPSA por medianas
ajuste<-Arima(inf1,order=c(10,1,1),
              fixed=c(0,0,0,0,0,,0,0,0,0,NA,
                     NA,#0,0,0,0,0,0,0,0,0,0,NA,
                      NA,NA,NA),
              season=list(order=c(0,1,1),12), 
              xreg=sipsaxreg[,c(1,2)])


####Validación del modelo basado en SIPSA calculado con medianas
BIC(ajuste)
summary(ajuste)
ajuste$coef[ajuste$coef !=0]/sqrt(diag(ajuste$var.coef))
res<- ajuste$residuals
yfit<-inf1-res
postscript("valarimax.eps",width = 7, height = 7)
par(mfrow=c(2,2))
plot(inf1,col="red",main="Datos originales vs \n valores ajustados",xlab="Tiempo")
legend("topleft", inset=0.05, legend=c("originales", "ajustados"),#,"SIPSA medias"),
       col=c("red", "blue"), lty=1:3, cex=0.6)

lines(yfit,col="blue",lwd=2)
acf(res,main="ACF de residuos",xlab="Rezago")
pacf(res,main="PACF de residuos",xlab="Rezago", ylab="ACF Parcial")
qqPlot(res,main="Qqplot de residuos",xlab="Cuantiles de la Dist. Normal")
dev.off()

#Autocorrelación serial
Box.test(res,type = "Box-Pierce")
#Normalidad en los residuos 
jarque.bera.test(res)
####Evaluación de pronóstico dentro de la muestra
sipsa2_1rezago<- ts(c(0,sipsa2[-length(sipsa2)]),start=c(9,2017),freq=12)
predt<-forecast(ajuste,xreg=sipsaxreg2[,c(1,2)])
inf2prono<-predt$mean
par(mfrow=c(1,1))
plot(predt,main="Pronóstico Fuera de la muestra",xlab="Periodo",ylab="IPC alimentos")
accuracy(inf2prono,inf2)
inf2prono.li <- predt$lower[,2]
inf2prono.ls <- predt$upper[,2]
#plot(inf2prono, type="l", col="blue", lwd=2, 
#     ylim=c(min(inf2prono.li,inf2prono.ls),max(inf2prono.li,inf2prono.ls)))
#lines(inf2prono.li, type="l",col="red", lty=2)
#lines(inf2prono.ls, type="l",col="red", lty=2)
postscript("pronoin.eps",width = 7, height = 5)
plot(inf2prono, type="l", col="blue", lwd=2, 
     main="Pronóstico Fuera de la Muestra",
     xlab="",
     ylab="Variación IPC alimentos", 
     xaxt="n",
     ylim=c(min(inf2prono.li,inf2prono.ls),max(inf2prono.li,inf2prono.ls)))
polygon(c(time(inf2prono.li),rev(time(inf2prono.ls))),
        c(inf2prono.li,rev(inf2prono.ls)),
        col="gray", border=NA)
lines(inf2prono, type="b", col="blue", lwd=2) 
lines(inf2, type="b", col="red", lty=2) 
meses<-c("Abr/18","Mar/18","Feb/18","Ene/18","Dec/17","Nov/17","Oct/17","Sept/17")
axis(1, at=rev(time(inf2prono.ls)),labels=meses, las=2)
legend("topleft", inset=0.03, legend=c("Inflación", "Pronóstico \n (SIPSA)"),
       col=c("red", "blue"), lty=1:2, cex=0.8)
dev.off()


##################Modelos finales#######################

ajuste_final<-Arima(inf, order=c(10,1,1),
                    fixed=c(0,0,0,0,0,0,0,0,0,NA,
                            NA,
                            NA,NA,NA),
                    season=list(order=c(0,1,1),12), 
                    xreg=sipsaxregtodo[-c(nrow(sipsaxregtodo)),c(1,2)])

summary(ajuste_final)
ajuste_final$coef[ajuste_final$coef !=0]/sqrt(diag(ajuste_final$var.coef))
res<- ajuste_final$residuals
yfit<-inf1-res
par(mfrow=c(2,2))
plot(inf1,col="red")
lines(yfit,col="blue",lwd=2)
acf(res)
pacf(res)
qqPlot(res)
#Autocorrelación serial
Box.test(res,type = "Box-Pierce")
#Normalidad en los residuos 
jarque.bera.test(res)

forecast(ajuste_final, xreg=sipsaxregtodo[nrow(sipsaxregtodo),c(1,2)], h=1)
(145.64*(-0.2725109)/100)+145.64 #pronóstico
(145.64*(-1.120084)/100)+145.64 #inferior 80
(145.64*(0.5750621)/100)+145.64 #superior 80
(145.64*(-1.568762)/100)+145.64 #inferior 95
(145.64*(1.02374)/100)+145.64 #superio 95
Pronóstico 145.2431
IC 95% (143.3553;147.131)
IC 80% (144.0087;146.4775)
##################### Series de tiempo multivariadas ##############

par(mfrow=c(1,1))
postscript("acfcruzado.eps",width = 7, height = 5)
ccf(sipsa,inf,type="correlation", main="Correlación cruzada entre variación en los \n precios de alimentos y el índice SIPSA")
dev.off()


adf1 <- ur.df(sipsa, type = "trend", lags = 2)
summary(adf1)
plot(adf1)

adf1 <- ur.df(inf, type = "trend", lags = 2)
summary(adf1)
plot(adf1)

plot(diff(inf), nc = 2, xlab = "")

#######################################
##### PARTE 1. Especificación del 
##### orden VAR. Incorporamos las
##### variables en niveles para poder
##### recoger la estructura de los
##### datos introduciendo tendencias

yt<-cbind(inf,sipsa)
ytinsample<-window(yt,start=c(2013,2),end=c(2017,8))
ytoutsample<-window(yt,start=c(2017,9))
plot(ytinsample)
colnames(ytinsample)<-c("inf","sipsa")
VARselect(ytinsample, lag.max = 8, type = "trend",season=12)$selection
#modelo1 <- MTS::VAR(ytinsample, 2)#, type = "trend")
#refVAR(modelo1,thres=1)
modelo <- vars::VAR(ytinsample, 1, type = "both",season=12)
summary(modelo)


#######################################
##### PARTE 2. Evaluando relaciones de
##### Causalidad entre las variables 
##### de interés

causality(modelo, cause="sipsa")
causality(modelo, cause="inf")

#######################################
##### PARTE 3. Diagnóstico
##### del modelo

est <- residuals(modelo)
postscript("Varvalida.eps",width = 7, height = 7)
par(mfrow=c(2,2))
acf(est[,1],main="Residuales para la varición \n en precios de alimentos")
acf(est[,1]^2,main="Residuos al cuadrado para \n variación en precios de alimentos")
acf(est[,2],main="Residuos para variación \n en índice SIPSA")
acf(est[,2]^2,main="Residuos al cuadrado para \n índice SIPSA")
dev.off()
postscript("Varvalida.eps",width = 7, height = 4)
par(mfrow=c(1,2))
acf(est[,1],main="Residuales para la varición \n en precios de alimentos")
acf(est[,1]^2,main="Residuos al cuadrado para \n variación en precios de alimentos")
acf(est[,2],main="Residuos para variación \n en índice SIPSA")
acf(est[,2]^2,main="Residuos al cuadrado para \n índice SIPSA")
dev.off()
par(mfrow=c(1,1))
postscript("Varvalida2.eps",width = 7, height = 7)
ccf(est[,1],est[,2],main="Correlograma cruzado entre \n residuos del modelo ajustado")
dev.off()


ser11 <- serial.test(modelo, lags.pt = 20, type = "PT.asymptotic")
ser11$serial

norm1 <- normality.test(modelo)
norm1$jb.mul

arch1 <- arch.test(modelo, lags.multi = 5)
arch1$arch.mul

postscript("fittedvar.eps",width = 8, height = 5)
par(mfrow=c(1,2))
plot(ts(modelo$varresult$inf$fitted.values,start=c(2013,2),frequency=12),ylim=c(-1.5,max(inf))
     ,main="Ajustados vs Observados \n modelo VAR. precios de alimentos",
     xlab="Periodo",ylab="Variación",col=2,lty=2)
lines(window(inf,start=c(2013,3),end=c(2017,8)))
legend("topleft", inset=0.01, legend=c("Observados", "Ajustados"),
       col=c("black", "red"), lty=1:2, cex=0.7)

plot(ts(modelo$varresult$sipsa$fitted.values,start=c(2013,2),frequency=12),ylim=c(-1.5,max(sipsa))
     ,main="Ajustados vs Observados \n modelo VAR. SIPSA",
     xlab="Periodo",ylab="Variación",col=2,lty=2)
lines(window(sipsa,start=c(2013,3),end=c(2017,8)))
legend("topleft", inset=0.01, legend=c("Observados", "Ajustados"),
       col=c("black", "red"), lty=1:2, cex=0.7)
dev.off()


adf1 <- ur.df(inf, type = "trend", lags = 2)
summary(adf1)
plot(adf1)

adf2 <- ur.df(diff(inf), type = "drift",lags = 1)
summary(adf2)
plot(adf2)

adf1 <- ur.df(sipsa, type = "trend", lags = 2)
summary(adf1)
plot(adf1)

adf2 <- ur.df(diff(inf), type = "drift",lags = 1)
summary(adf2)
plot(adf2)


plot(modelo,names="sipsa")
plot(modelo,names="inf")


plot(arch1, names = "sipsa")
plot(arch1, names = "inf")


postscript("fluctuacion.eps",width = 7, height = 5)

plot(stability(modelo), nc = 2,main="Estabilidad estructural \n del modelo propuesto"
     , ylab="Fluctuación",xlab="Tiempo")
dev.off()



#######################################
##### PARTE 3. Análisis de 
##### de Cointegración
##### Para tres o más variables. 
vecm <- ca.jo(ytinsample, type = "eigen", ecdet = "none", 
              K =2, spec = "transitory")
summary(vecm)

vecm <- ca.jo(ytinsample, type = "trace", ecdet = "none", K = 2, spec = "transitory")
(vecm.r1 <- cajorls(vecm, r = 1))

var_vecm <- vec2var(vecm,r=1)
var_vecm$datamat
eigen(var_vecm$A$A1)
svec.irf <- irf(var_vecm, response = "inf", n.ahead = 12, boot = TRUE)
plot(svec.irf)

par(mfrow=c(1,1))
fevd <- fevd(var_vecm, n.ahead = 15)

plot(fevd, plot.type="single", col=c("gray","green","red","blue"))


#######################################   
#### Modelo SVECM.
#### Es posible encontrar una 
#### representación estructural para
#### el conjunto de variables bajo
#### análisis.
#### Para la definciión del SVECM
#### se requiere la definición de la
#### matriz B. 
#### Los autores establecen que
#### existe retornos constantes a 
#### escala, así que la
#### productividad es únicamente
#### determinada por los choques
#### del producto. 
#### LP.B1,j = 0 j=2,3,4
#### LP.B1.4 = 0
#### SR.B4,2 = 0

vecm <- ca.jo(ytinsample, type = "trace", ecdet = "trend", K = 2, spec = "transitory",season=12)
jo.results<- summary(vecm)
tsDyn::fitted(vecm,level="original")
var_vecm <- vec2var(vecm,r=1)
lines(ts(fitted(var_vecm)))
plot(ytinsample)
length(fitted(var_vecm))
ser11 <- serial.test(var_vecm, lags.pt = 20, type = "PT.asymptotic")
ser11$serial

extract.summary.lm <- function (model, include.rsquared = TRUE, include.adjrs = TRUE, 
                                include.nobs = TRUE, include.fstatistic = FALSE, include.rmse = TRUE, 
                                ...) 
{
  s <- model;
  names <- rownames(s$coef)
  co <- s$coef[, 1]
  se <- s$coef[, 2]
  pval <- s$coef[, 4]
  rs <- s$r.squared
  adj <- s$adj.r.squared
  n <- length(s$residuals)
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.rsquared == TRUE) {
    gof <- c(gof, rs)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.adjrs == TRUE) {
    gof <- c(gof, adj)
    gof.names <- c(gof.names, "Adj. R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.fstatistic == TRUE) {
    fstat <- s$fstatistic[[1]]
    gof <- c(gof, fstat)
    gof.names <- c(gof.names, "F statistic")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.rmse == TRUE && !is.null(s$sigma[[1]])) {
    rmse <- s$sigma[[1]]
    gof <- c(gof, rmse)
    gof.names <- c(gof.names, "RMSE")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  tr <- createTexreg(coef.names = names, coef = co, se = se, 
                     pvalues = pval, gof.names = gof.names, gof = gof, gof.decimal = gof.decimal)
  return(tr)
}

setMethod("extract",  signature = 'summary.lm', definition = extract.summary.lm)

mysum <- summary(cajorls(vecm)$rlm)
screenreg(list(mysum[[1]], mysum[[2]]))
texreg(list(mysum[[1]], mysum[[2]]), dcolumn = TRUE, booktabs = TRUE,
       use.packages = FALSE, label = "tab:3", caption = "Two linear models.",
       float.pos = "hb!")

plot(ts(vecm@R0[,1],start=c(2013,2),frequency = 12))
lines()
ser11 <- serial.test(modelo, lags.pt = 16, type = "PT.asymptotic")
ser11$serial

str(var_vecm)
jarque.bera.test(var_vecm$resid[,1][-41])

norm1 <- normality.test(var_vecm)
norm1$jb.mul

arch1 <- arch.test(var_vecm)
arch1$arch.mul

vecm.r2<- cajorls(vecm,r=1)
class(jo.results)

slotNames(jo.results)


SR <- matrix(1, nrow = 2, ncol = 2)
SR[1,2]<-NA 
SR[2,1]<-NA
LR <- diag(2)
LR[1,2]<-NA


#vecm <- ca.jo(ytinsample, type = "trace", ecdet = "const", K = 2, spec = "transitory")
#svec <- SVEC(vecm, LR = SR, SR = LR, r = 1, lrtest = FALSE, boot = TRUE, runs = 100)
#summary(svec)


SR <- matrix(NA, nrow = 2, ncol = 2)
SR[1,2]<-0
SR[2,1]<-0

LR <- matrix(NA, nrow = 2, ncol = 2)
LR[2,1]<-0
svec <-SVEC(vecm, LR = LR, SR = SR, r = 1, lrtest = FALSE, boot = TRUE, runs = 100)

salida<-SVAR(modelo, estmethod = "scoring", Amat = LR,Bmat = SR,
             max.iter = 100, maxls = 1000, conv.crit = 1.0e-8)
summary(salida) 

###Qqplots para VEC

plot(vecm)
fitted.values(vecm)
#######################################
##### PARTE 5. Funciones de Impulso
##### Respuesta
#######################################

svec.irf <- irf(salida, response = "inf", n.ahead = 12, boot = TRUE)


plot(svec.irf)
#,main="Impulso de respuesta de la variación\n en los precios de los alimentos (SIPSA)"
#     , xlab="",ylab="Variación")
#######################################
##### PARTE 6. Descomposición de la
##### Varianza de Pronóstico

par(mfrow=c(1,1))
fevd <- fevd(salida, n.ahead = 15)

plot(fevd, plot.type="single", col=c("gray","green","red","blue"))

#     ,ylab="Porcentaje",xlab="Horizonte", main="")



#######################################
##### PARTE 7. Pronóstico
##### el modelo de referencia es el 
##### VECM con r=1


var_vecm <- vec2var(vecm,r=1)
arch.test(var_vecm)
normality.test(var_vecm)
serial.test(var_vecm)
svec.irf <- irf(var_vecm, response = "inf", n.ahead = 12, boot = TRUE)
plot(svec.irf)
fevd <- fevd(var_vecm, n.ahead = 15)
plot(fevd, plot.type="single", col=c("gray","green","red","blue"))
predict(salida)
frc <- predict(modelo,n.ahead=8,ci=0.95) #Pronóstico con VAR  
frc <- predict(var_vecm,n.ahead=8,ci=0.95) #Pronóstico con VECM
plot(frc)
fanchart(frc)
par(mfrow=c(1,1))

options(digits=3)
accuracy(frc$fcst$inf[,1],ytoutsample[,1])
accuracy(frc$fcst$sipsa[,1],ytoutsample[,2])


####Pronósticos VAR#####
frc <- predict(modelo,n.ahead=8,ci=0.95) #Pronóstico con VAR  
inf2pronomedias.li<-ts(frc$fcst$inf[,2],start=c(2017,9),frequency = 12)
inf2pronomedias.ls<-ts(frc$fcst$inf[,3],start=c(2017,9),frequency = 12)
prono_infl_ts<-ts(frc$fcst$inf[,1],start=c(2017,9),frequency = 12)
postscript("pronoinvarinf.eps",width = 7, height = 5)
plot(prono_infl_ts, type="l", col="blue", lwd=2, 
     main="Pronóstico de la inflación dentro de la Muestra \n Modelo VAR",
     xlab="",
     ylab="Inflación de alimentos",
     xaxt="n",
     ylim=c(min(inf2pronomedias.li,inf2pronomedias.ls),max(inf2pronomedias.li,inf2pronomedias.ls)))
meses<-c("Sept/17","Oct/17","Nov/17","Dec/17","Ene/18","Feb/18","Mar/18","Abr/18")
axis(1, at=rev(time(inf2pronomedias.ls)),labels=meses, las=2)
polygon(c(time(inf2pronomedias.li),rev(time(inf2pronomedias.ls))),
        c(inf2pronomedias.li,rev(inf2pronomedias.ls)),
        col="gray", border=NA)

lines(prono_infl_ts, type="b", col="blue", lwd=2) 
lines(ytoutsample[,1], type="b", col="red", lty=2) 

legend("topleft", inset=0.03, legend=c("Inflación", "Pronóstico"),
       col=c("red", "blue"), lty=1:2, cex=0.8)
dev.off()




inf2pronomedias.li<-ts(frc$fcst$sipsa[,2],start=c(2017,9),frequency = 12)
inf2pronomedias.ls<-ts(frc$fcst$sipsa[,3],start=c(2017,9),frequency = 12)
prono_infl_ts<-ts(frc$fcst$sipsa[,1],start=c(2017,9),frequency = 12)
postscript("pronoinvarsipsa.eps",width = 7, height = 5)
plot(prono_infl_ts, type="l", col="blue", lwd=2, 
     main="Pronóstico de SIPSA Dentro de la Muestra \n Modelo VAR",
     xlab="",
     xaxt="n",
     ylab="Inflación de alimentos", 
     ylim=c(min(inf2pronomedias.li,inf2pronomedias.ls),max(inf2pronomedias.li,inf2pronomedias.ls)))
meses<-c("Sept/17","Oct/17","Nov/17","Dec/17","Ene/18","Feb/18","Mar/18","Abr/18")
axis(1, at=rev(time(inf2pronomedias.ls)),labels=meses, las=2)

polygon(c(time(inf2pronomedias.li),rev(time(inf2pronomedias.ls))),
        c(inf2pronomedias.li,rev(inf2pronomedias.ls)),
        col="gray", border=NA)
lines(prono_infl_ts, type="b", col="blue", lwd=2) 
lines(ytoutsample[,2], type="b", col="red", lty=2) 

legend("topleft", inset=0.03, legend=c("SIPSA", "Pronóstico"),
       col=c("red", "blue"), lty=1:2, cex=0.8)
dev.off()


####Pronóstico VECM### 

frc <- predict(var_vecm,n.ahead=8,ci=0.95) #Pronóstico con VAR  
inf2pronomedias.li<-ts(frc$fcst$inf[,2],start=c(2017,9),frequency = 12)
inf2pronomedias.ls<-ts(frc$fcst$inf[,3],start=c(2017,9),frequency = 12)
prono_infl_ts<-ts(frc$fcst$inf[,1],start=c(2017,9),frequency = 12)
postscript("pronoinvarinfvec.eps",width = 7, height = 5)
plot(prono_infl_ts, type="l", col="blue", lwd=2, 
     main="Pronóstico de la inflación dentro de la Muestra \n Modelo VEC",
     xlab="",
     ylab="Inflación de alimentos",
     xaxt="n",
     ylim=c(min(inf2pronomedias.li,inf2pronomedias.ls),max(inf2pronomedias.li,inf2pronomedias.ls)))
meses<-c("Sept/17","Oct/17","Nov/17","Dec/17","Ene/18","Feb/18","Mar/18","Abr/18")
axis(1, at=rev(time(inf2pronomedias.ls)),labels=meses, las=2)
polygon(c(time(inf2pronomedias.li),rev(time(inf2pronomedias.ls))),
        c(inf2pronomedias.li,rev(inf2pronomedias.ls)),
        col="gray", border=NA)

lines(prono_infl_ts, type="b", col="blue", lwd=2) 
lines(ytoutsample[,1], type="b", col="red", lty=2) 

legend("topleft", inset=0.03, legend=c("Inflación", "Pronóstico"),
       col=c("red", "blue"), lty=1:2, cex=0.8)
dev.off()




inf2pronomedias.li<-ts(frc$fcst$sipsa[,2],start=c(2017,9),frequency = 12)
inf2pronomedias.ls<-ts(frc$fcst$sipsa[,3],start=c(2017,9),frequency = 12)
prono_infl_ts<-ts(frc$fcst$sipsa[,1],start=c(2017,9),frequency = 12)
postscript("pronoinvarsipsavec.eps",width = 7, height = 5)
plot(prono_infl_ts, type="l", col="blue", lwd=2, 
     main="Pronóstico de SIPSA Dentro de la Muestra \n Modelo VEC",
     xlab="",
     xaxt="n",
     ylab="Inflación de alimentos", 
     ylim=c(min(inf2pronomedias.li,inf2pronomedias.ls),max(inf2pronomedias.li,inf2pronomedias.ls)))
meses<-c("Sept/17","Oct/17","Nov/17","Dec/17","Ene/18","Feb/18","Mar/18","Abr/18")
axis(1, at=rev(time(inf2pronomedias.ls)),labels=meses, las=2)

polygon(c(time(inf2pronomedias.li),rev(time(inf2pronomedias.ls))),
        c(inf2pronomedias.li,rev(inf2pronomedias.ls)),
        col="gray", border=NA)
lines(prono_infl_ts, type="b", col="blue", lwd=2) 
lines(ytoutsample[,2], type="b", col="red", lty=2) 

legend("topleft", inset=0.03, legend=c("SIPSA", "Pronóstico"),
       col=c("red", "blue"), lty=1:2, cex=0.8)
dev.off()



A1<-matrix(c(0.453,0.150,0.199,0.188),2,2,byrow=T)
A2<-matrix(c(0.0224,-0.0803,0.3812,-0.2310),2,2,byrow=T)
A<-matrix(c(1,-0.2598,0,1),2,2,byrow=T)

A%*%A1
A%*%A2

errores<-matrix(c(27.25,27.38,27.38,105.39),2,2,byrow=T)
errores/100


gg<-c(1.89,1,0.83,3.75,4.31,4.78,0.6,2.99)
sum(round(gg/sum(gg),3))

library(astsa)
xt <- arima.sim(n=200, order=c(1,0,1), ar = 0.6, ma = -0.1,
                  seasonal=list(order=c(1,1,0),sar=0.3),period=4)
xt<-arima.sim(n = 200, list(ar = 0.8897, ma =-0.2279,sar=0.12),freq=4)
arima
plot(xt)
aju<-astsa::sarima(xt,p=1,d=0,q=1,P=1,D=1,Q=0,S=4)
dif<-diff(aju)

