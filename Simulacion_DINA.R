rm(list=ls())

######################################
# Estudio de simulación: Modelo DINA
######################################
dir <- "D:/Yuriko Sosa/Otros/PUCP/Tesis/3 Análisis del modelo"
setwd(dir)
dir.result <- paste0(dir,"/Simulaciones/Items fijos", sep="")
dir.images <- paste0(dir,"/Simulaciones/Items fijos/Convergencia", sep="")

# Cargar paquetes
library(dina)
library(foreign)
library(CDM)

# Matriz Q fija
Q <- read.csv("Matriz_Q_Mate_2S_C1_C6.csv", sep=",", T)
Q <- as.matrix(Q[,-(1:2)])
# Matriz de posibles combinaciones según K
A <- read.csv("matriz_CxK.csv", T, sep=",")
A <- as.matrix(A[,-1])


# Variables fijas
K <- ncol(Q) # número de habilidades
J <- nrow(Q) # número de ítems
chainLength <- 50000 # número de iteraciones
burnin <- 1000 # Quemar las primeras 1000
thin <- 5

# Definiendo los escenarios
nrep <- 50 # número de repeticiones
escenarios <- list(c("100", "rep(.1,J)", "rep(.1,J)", "rep(1/(2^K),2^K)"), # 1
                   c("300", "rep(.1,J)", "rep(.1,J)", "rep(1/(2^K),2^K)"), # 2
c("500", "rep(.1,J)", "rep(.1,J)", "rep(1/(2^K),2^K)"), # 3
c("100", "rep(.2,J)", "rep(.2,J)", "rep(1/(2^K),2^K)"), # 4
c("300", "rep(.2,J)", "rep(.2,J)", "rep(1/(2^K),2^K)"), # 5
c("500", "rep(.2,J)", "rep(.2,J)", "rep(1/(2^K),2^K)"), # 6
c("100", "rep(.1,J)", "rep(.2,J)", "rep(1/(2^K),2^K)"), # 7
c("300", "rep(.1,J)", "rep(.2,J)", "rep(1/(2^K),2^K)"), # 8
c("500", "rep(.1,J)", "rep(.2,J)", "rep(1/(2^K),2^K)"), # 9
c("100", "rep(.2,J)", "rep(.1,J)", "rep(1/(2^K),2^K)"), # 10
c("300", "rep(.2,J)", "rep(.1,J)", "rep(1/(2^K),2^K)"), # 11
c("500", "rep(.2,J)", "rep(.1,J)", "rep(1/(2^K),2^K)"), # 12
c("100", "rep(.1,J)", "rep(.1,J)", "c(rep(0.021,11),rep(0.042,5))"), # 13
c("300", "rep(.1,J)", "rep(.1,J)", "c(rep(0.021,11),rep(0.042,5))"), # 14
c("500", "rep(.1,J)", "rep(.1,J)", "c(rep(0.021,11),rep(0.042,5))"), # 15
c("100", "rep(.2,J)", "rep(.2,J)", "c(rep(0.021,11),rep(0.042,5))"), # 16
c("300", "rep(.2,J)", "rep(.2,J)", "c(rep(0.021,11),rep(0.042,5))"), # 17
c("500", "rep(.2,J)", "rep(.2,J)", "c(rep(0.021,11),rep(0.042,5))"), # 18
c("100", "rep(.1,J)", "rep(.2,J)", "c(rep(0.021,11),rep(0.042,5))"), # 19
c("300", "rep(.1,J)", "rep(.2,J)", "c(rep(0.021,11),rep(0.042,5))"), # 20
c("500", "rep(.1,J)", "rep(.2,J)", "c(rep(0.021,11),rep(0.042,5))"), # 21
c("100", "rep(.2,J)", "rep(.1,J)", "c(rep(0.021,11),rep(0.042,5))"), # 22
c("300", "rep(.2,J)", "rep(.1,J)", "c(rep(0.021,11),rep(0.042,5))"), # 23
c("500", "rep(.2,J)", "rep(.1,J)", "c(rep(0.021,11),rep(0.042,5))")) # 24

######################################################
# Definir objetos
mG_bayes <- matrix(0,J,nrep)
mS_bayes <- matrix(0,J,nrep)
PIoutput_bayes <- matrix(0,2^K,nrep)
mG_EM <- matrix(0,J,nrep)
mS_EM <- matrix(0,J,nrep)
PIoutput_EM <- matrix(0,2^K,nrep)

resbg <- matrix(0,J,nrep); resemg<- matrix(0,J,nrep)
resbs<- matrix(0,J,nrep); resems <- matrix(0,J,nrep)
rmse_bayes_g<-matrix(0,J); rmse_em_g<-matrix(0,J)
rmse_bayes_s<-matrix(0,J); rmse_em_s<-matrix(0,J)
errorG1<-matrix(0,J); errorG2<-matrix(0,J)
errorS1 <-matrix(0,J); errorS2<-matrix(0,J)
resbpi <- matrix(0,2^K,nrep);resempi<- matrix(0,2^K,nrep)
rmse_bayes_PI <- matrix(0,2^K); rmse_em_PI <- matrix(0,2^K)
errorPI1 <-matrix(0,2^K);errorPI2<-matrix(0,2^K)
errorPI1 <-NA;errorPI2<-NA

resumen.guessing <- matrix(0,J,4)
resumen.slipping <- matrix(0,J,4)
resumen.pi <- matrix(0,2^K,4)

######################################################
ptm <- proc.time()

for (i in 2:length(escenarios)){
  ptm <- proc.time()
  N <- eval(parse(text=escenarios[[i]][1]))
  s <- eval(parse(text=escenarios[[i]][2]))
  g <- eval(parse(text=escenarios[[i]][3]))
  PI <- eval(parse(text=escenarios[[i]][4]))
  
  #########################################
  # Simulación
  #########################################
  for (f in 1:nrep){
    
    CL <- c((1:(2^K))%*%rmultinom(n=N,size=1,prob=PI))
    ## Definiendo la matriz de los posibles
    ## perfiles de conocimiento (combinaciones A)
    #A <- rep(0,K)
    #for(j in 1:K){
     # temp <- combn(1:K, m=j)
    #  tempmat <- matrix(0, ncol(temp), K)
     # for(j in 1:ncol(temp)) tempmat[j,temp[,j]] = 1
    #  A <- rbind(A,tempmat)
    #}
    Alphas <- A[CL,] # Matriz de alphas A_{NxK}
    gen <- DINAsim(Alphas,Q,s,g); Y_sim <- gen$Y # Data simulada
    
    # Estimación
    ### Muestreador de Gibbs
    outchain <- DINA_Gibbs(Y_sim, Amat=A, Q, chain_length=chainLength)
    ### Algoritmo EM
    outdin <- din(Y_sim, q.matr = Q, progress=FALSE)
    # Guessing parameter
    mG_bayes[,f] <- apply(outchain$GamS[,seq(burnin, chainLength, by = thin)],1,mean)
    # Slipping parameter
    mS_bayes[,f] <- apply(outchain$SigS[,seq(burnin, chainLength, by = thin)],1,mean)
    # Guessing parameter
    mG_EM[,f] <- outdin$item[,-1]$guess
    # Slipping parameter
    mS_EM[,f] <- outdin$item[,-1]$slip
    # Proporción de alumnos en clases latentes
    PIoutput_bayes[,f] <- apply(outchain$PIs[,seq(burnin, chainLength, by = thin)],1,mean)
    # Proporción de alumnos en clases latentes
    PIoutput_EM[,f] <- outdin$attribute.patt$class.prob
  } # cierre de repetición
  
  ##############################################################
  # Plots de convergencia para parámetros de los ítems G y S
  #D:\Yuriko Sosa\Otros\PUCP\Tesis\3 Análisis del modelo\Simulaciones\Items fijos\Convergencia\Guessing 
  for(j in 1:J){
    png(filename=paste0(dir.images,"/Guessing/cadena_item",j,"_esc",i,".png", sep=""))
    ts.plot(outchain$GamS[j,seq(burnin, chainLength, by = thin)], col="dodgerblue",
            xlab="Iteraciones" , ylab = "Adivinación", main="")
    dev.off()
    png(filename=paste0(dir.images,"/Guessing/acf_item",j,"_esc",i,".png", sep=""))
    acf(outchain$GamS[j,seq(burnin, chainLength, by = thin)], col= "red",
        xlab="Lag", ylab="Autocorrelación", main="")  
    dev.off()
    
    png(filename=paste0(dir.images,"/Slipping/cadena_item",j,"_esc",i,".png", sep=""))
    ts.plot(outchain$SigS[j,seq(burnin, chainLength, by = thin)], col="dodgerblue",
            xlab="Iteraciones" , ylab = "Adivinación", main="")
    dev.off()
    png(filename=paste0(dir.images,"/Slipping/acf_item",j,"_esc",i,".png", sep=""))
    acf(outchain$SigS[j,seq(burnin, chainLength, by = thin)], col= "red",
        xlab="Lag", ylab="Autocorrelación", main="")  
    dev.off()
  }
  
  # Plots de convergencia para PIs
  for(k in 1:2^K){
    png(filename=paste0(dir.images,"/Pi/cadena_hab",k,"_esc",i,".png", sep=""))
    ts.plot(outchain$PIs[k,seq(burnin, chainLength, by = thin)], col="dodgerblue",
            xlab="Iteraciones" , ylab = "Probabilidad", main="")
    dev.off()
    png(filename=paste0(dir.images,"/Pi/acf_hab",k,"_esc",i,".png", sep=""))
    acf(outchain$PIs[k,seq(burnin, chainLength, by = thin)], col= "red",
        xlab="Lag", ylab="Autocorrelación", main="")  
    dev.off()
  }
  
  # Resultados de los escenarios
  for(h in 1:J){
    for(cont in 1:nrep){
      resbg[h,cont] <- (mG_bayes[h,cont] - g[h])^2
      resemg[h,cont] <- (mG_EM[h,cont] - g[h])^2
      resbs[h,cont] <- (mS_bayes[h,cont] - s[h])^2
      resems[h,cont] <- (mS_EM[h,cont] - s[h])^2}
    
    rmse_bayes_g[h] <- sqrt(sum(resbg[h,])/nrep)
    rmse_em_g[h] <- sqrt(sum(resemg[h,])/nrep)
    rmse_bayes_s[h] <- sqrt(sum(resbs[h,])/nrep)
    rmse_em_s[h] <- sqrt(sum(resems[h,])/nrep) 
    errorG1[h] <- mean(mG_bayes[h,]- g[h])
    errorG2[h] <- mean(mG_EM[h,]- g[h])
    errorS1[h] <- mean(mS_bayes[h,]- s[h])
    errorS2[h] <- mean(mS_EM[h,]- s[h])
  }
  resumen.guessing <- cbind(round(rmse_bayes_g,3),round(errorG1,3),
                            round(rmse_em_g,3),round(errorG2,3))
  
  resumen.slipping <- cbind(round(rmse_bayes_s,3),round(errorS1,3),
                            round(rmse_em_s,3),round(errorS2,3)) 
  
  for(it in 1:2^K){
    for(cont in 1:nrep){
      resbpi[it,cont] <- (PIoutput_bayes[it,cont] - PI[it])^2
      resempi[it,cont] <- (PIoutput_EM[it,cont] - PI[it])^2}
    rmse_bayes_PI[it] <- sqrt(sum(resbpi[it,])/nrep)
    rmse_em_PI[it] <- sqrt(sum(resempi[it,])/nrep)
    errorPI1[it] <- mean(PIoutput_bayes[it,] - PI[it])
    errorPI2[it] <- mean(PIoutput_EM[it,] - PI[it])}
  
  resumen.pi <- cbind(round(rmse_bayes_PI,3),round(errorPI1,3),
                      round(rmse_em_PI,3),round(errorPI2,3)) 
  
  # Etiquetas de resumen
  colnames(resumen.guessing) <- c("RMSE_Gibbs","Sesgo_Gibbs",
                                  "RMSE_EM","Sesgo_EM")
  colnames(resumen.slipping) <- c("RMSE_Gibbs","Sesgo_Gibbs",
                                  "RMSE_EM","Sesgo_EM")
  colnames(resumen.pi) <- c("RMSE_Gibbs","Sesgo_Gibbs",
                            "RMSE_EM","Sesgo_EM")
  
  write.csv(resumen.guessing, file=paste(dir.result,"/Guessing_Escenario",i,".csv",sep = ""))
  write.csv(resumen.slipping, file=paste(dir.result,"/Slipping_Escenario",i,".csv", sep = ""))
  write.csv(resumen.pi, file=paste(dir.result,"/PI_Escenario",i,".csv", sep = ""))
  print(proc.time() - ptm)
} 

proc.time() - ptm
