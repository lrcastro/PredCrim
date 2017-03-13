#' Spatio-temporal prediction of crimes in Manaus
#'
#' This function develops an one-step ahead prediction algorithm based on Poisson dynamic generalized
#' linear model (DGLM) for the count of crime in Manaus.
#'
#' @param x A data frame containig the following columns:LAT,LONG,DATA,DIA.SEMANA,
#' HORA,PERIODO,Crime.
#' @param type Character string giving the type of crime to make the predictions. In this version only "ROUBO" is available.
#' @param by1 Character string specifying the periodicity of the counts. In this version only "week" is available.
#' @param k0 Number of periods to use as prior information.
#' @param m0 First posterior moment at t-1.
#' @param c0 Second posterior moment at t-1.
#' @param nMap Number of maps to check the quality of predictions.
#' @param delta Discount factor.
#' @param UsePrior Default is \code{FALSE}, if \code{TRUE} the algorithm will take the informations
#' of the last prediction from \code{crimes1415} data set to use as prior. Incompatible with \code{k0}.
#'
#' @return A map that provides predictions of crime intensity and the mean square error.
#'
#' @references Harrison, Jeff, and Mike West. \emph{Bayesian forecasting & dynamic models}. New York: Springer, 1999.
#' @references TRIANTAFYLLOPOULOS, K. Dynamic z
#' \emph{generalized linear models for non-Gaussian time series forecasting}. arXiv preprint arXiv:0802.0219, 2008.
#'
#' @examples
#' data(crimes1415)
#' predCrimMao(crimes1415,type="ROUBO",by1="week",nMap=2)
#'
#' @import rgeos
#' @import rgdal
#' @import maptools
#' @importFrom  spatstat as.owin
#' @importFrom  spatstat quadratcount
#' @importFrom  spatstat intensity
#' @importFrom  spatstat owin
#' @importFrom  spatstat quadrats
#' @importFrom  spatstat ppp
#' @importFrom  spatstat as.tess
#' @importFrom  spatstat tess
#' @import spdep
#' @import raster
#' @importFrom  plyr ldply
#' @import sp
#' @importFrom MASS fitdistr
#' @export

predCrimMao<-function(x,type="ROUBO",nMap=2,delta=0.89,by1="week",k0=NULL,UsePrior=FALSE,m0=NULL,c0=NULL){

  if(is.null(k0)!=TRUE & UsePrior!=FALSE) stop("UsePrior incompatible with k0.")
  if(UsePrior==TRUE){K0=1}

  data("mao")
  x$DATA<-as.Date(x$DATA.DO.FATO,"%d/%m/%Y")
  x<-x[order(x$DATA),]

  #Selecionar os registros de roubo para os dois primeiros meses da base de dados
  roubo<-x[x$NATUREZA==type,]
  dias<-seq(from=min(x$DATA),to=max(x$DATA),by=by1)
  serie=list() #separar por dia, cada no da lista representa um dia t
  for(i in 1:length(dias)){
    if(i < length(dias)){
      serie[[i]]=roubo[roubo$DATA>=dias[i]&roubo$DATA<dias[i+1],]
    }else{
      serie[[i]]=roubo[roubo$DATA>=dias[i],]
    }
  }

  if(length(unique(serie[[length(serie)]]$DATA))<7){serie=serie[1:(length(serie)-1)]}
  cat(length(serie), "Semanas completas encontradas")

  #transformar as coordenadas de roubos em pontos espaciais para cada t
  spt<-list()
  for(i in 1:length(serie)){
    spt[[i]]<-SpatialPoints(cbind(serie[[i]]$LONG,serie[[i]]$LAT),proj4string=CRS("+proj=longlat"))
  }

  #As coordenadas estao em long/lat necessario converter para UTM
  spt.utm<-list()
  for(i in 1:length(serie)){
    spt.utm[[i]]<-spTransform(spt[[i]],CRS("+init=epsg:32721"))
  }

  #transformar os pontos espaciais em ppp
  spp<-list()
  for(i in 1:length(serie)){
    spp[[i]]<-ppp(spt.utm[[i]]$coords.x1,spt.utm[[i]]$coords.x2,window = as.owin(mao))
  }

  #Criar a grade, e contar o número de ocorrências em cada celula, para cada tempo t
  grid1<-list()
  countdf<-list()
  tt=quadrats(as.owin(mao),nx=30,keepempty = T) #para que mostre todos os pontos inclusive os fora da interceção
  for(i in 1:length(serie)){
    #quad[[i]]<-quadrats(spp[[i]],30,30)
    grid1[[i]]<-quadratcount(spp[[i]],tess=as.tess(tt))
    #plot(grid)
    countdf[[i]]<-as.data.frame(grid1[[i]])
  }

  #################################
  # Criar a grade
  #################################
  grid <- raster(extent(mao)) #criar um raster vazio
  ncol(grid)<-30 #escolher o numero de celulas
  nrow(grid)<-30

  #Transformar o grid para a mesma coordenada do shape
  proj4string(grid)<-proj4string(mao)

  #transformar o raster em um polígono datafrme
  gridpolygon <- rasterToPolygons(grid)
  ##############################
  # Obter a matriz de vizinhaça
  ##############################
  nb<-poly2nb(gridpolygon,queen=TRUE) #cria a estrutura nb
  plot(nb,coordinates(gridpolygon),col="blue") #mostra o grafico da estrutura
  w<-nb2listw(nb,style = 'B') #transforma os pesos em uma lista
  w1<-listw2mat(w) #matriz de proximidade

  ################################
  # Estimação da Matriz V_g
  #################################
  #Estimação dos parâmetros, beta, gama e sigma^2
  #Matrizes Necessarias (ver pag 23 spatial_GEE)
  M<-diag(1/rowSums(w1))
  I<-diag(rep(1,900))


  #Numero de pontos esperados em cada quadrado
  ei<-c()
  e<-list()
  for(t in 1:length(countdf)){
    for(i in 1:nrow(countdf[[1]])){
      ei[i]<-(1*sum(countdf[[t]]$Freq))/900
    }
    e[[t]]<-ei
  }
  #Imagem gravada até aqui!!

  tot=length(serie)

  y.st1=c()
  y.st2=c()
  y.star1=list()
  y.star2=list()
  for(l in 1:tot){
    for(i in 1:900){
      y.st1[i]<-log((countdf[[l]]$Freq[i]+0.5-ifelse(countdf[[l]]$Freq[i]==0,1,countdf[[l]]$Freq[i])/4)/e[[l]][i])
      y.st2[i]<-max(0,log(countdf[[l]]$Freq[i]+0.5)^2 - 1/(countdf[[l]]$Freq[i]+1))
    }
    y.star1[[l]]<-y.st1
    y.star2[[l]]<-y.st2
  }


  if(UsePrior==FALSE){
  #######################################
  #     MLGD Dinâmico Poisson
  #######################################
  if(is.null(k0)) k0 = 5
  j=1
  G = matrix(1)
  F. = matrix(1)
  a = r = f1 = q = s = r = f.st = R = mm = pred1 =  matrix(rep(0,900),nrow = 900,ncol=tot)
  m = matrix(rep(-1,900),nrow = 900,ncol=tot)
  c = matrix(rep(0.03,900),nrow = 900,ncol=tot)
  f.st = q.st = Q.st  = s.t1 = r.t1 =  matrix(rep(0,900),nrow = 900,ncol=tot)
  alfa=0.0001
  beta=0.1
  emv=NULL
  cont0=matrix(NA,nrow = 900,ncol=k0)
  for(i in 1:k0){
    cont0[,i]=countdf[[i]]$Freq
  }

  for(i in 1:900){
    if(any(cont0[i,]==0)){
      r.t1[i,1]=alfa
      s.t1[i,1]=beta
    }else{
      emv=fitdistr(cont0[i,],densfun = "gamma")
      r.t1[i,1] = emv$estimate[1]
      s.t1[i,1] = emv$estimate[2]
    }
  }
  # for(i in 1:900){
  #   r.t1[i,1] = (rowMeans(cont0)[i]+(alfa/k0))/(1+beta/k0)
  # }

  # s.t1[,1] = r.t1[,1]/(k0+beta)
  #
  # r.t1[,1]/s.t1[,1]
  # r.t1=matrix(0,nrow=900,ncol=tot)
  # r.t1[,1]=rep(alfa,900)
  #
  # s.t1=matrix(0,nrow=900,ncol=tot)
  # s.t1[,1]=rep(beta,900)

  f=matrix(rep(0,901),nrow = 901,ncol=tot)
  media = matrix(rep(0,901),nrow = 901,ncol=tot)
  #t=k0
  n=900

  ##Predição para os dados iniciais
  for(t in 2:k0){

    gama=t(as.matrix(y.star1[[t]])-as.matrix(y.star1[[t-1]]))%*%M%*%w1%*%(as.matrix(y.star1[[t]])-as.matrix(y.star1[[t-1]])) /
      (t(as.matrix(y.star1[[t]])-as.matrix(y.star1[[t-1]]))%*%w1%*%M^2%*%w1%*%(as.matrix(y.star1[[t]])-as.matrix(y.star1[[t-1]])) +
         sum(rowSums(w1^2)/colSums(w1^2)*(y.star2[[t]]-y.star1[[t]]^2)))
    #
    sig=(1/n) * (t(as.matrix(y.star1[[t]])-as.matrix(y.star1[[t-1]]))%*%(solve(M)-as.vector(gama)*w1)%*%(as.matrix(y.star1[[t]])-
                                                                                                           as.matrix(y.star1[[t-1]]))+sum(colSums(w1)*(y.star1[[t]]-y.star2[[t]]^2)))

    v_g = solve(diag(1,900,900)-as.vector(gama)*M%*%w1)%*%M
    sigVg = as.vector(sig)*v_g


    y = countdf[[t]]$Freq

    for(i in 1:900){
      a[i,t] = G%*%m[i,t-1]

      R[i,t] = (G%*%c[i,t-1]%*%t(G)+sigVg[i,i])
      #R[i,t] = sigVg[i,i]

      #cat("a",a,"\n")
      #Passo 2 - Previsão a um passo
      f[i,t] = t(F.)%*%a[i,t]
      q[i,t] = t(F.)%*%R[i,t]%*%F.

      #com fator de descontos
      r.t1[i,t] = delta[j]*(r.t1[i,t-1]+y[i])+1-delta[j]
      s.t1[i,t] = delta[j]*(s.t1[i,t-1]+1)
      media[i,t] = r.t1[i,t]/s.t1[i,t]

      #Passo 3 -  Atualização
      f.st[i,t] = digamma(r.t1[i,t]+y[i])-log(s.t1[i,t]+1)
      Q.st[i,t] = trigamma(r.t1[i,t]+y[i])

      m[i,t]  = a[i,t] + R[i,t]%*%F.*(f.st[i,t] - f[i,t])/q[i,t]
      c[i,t] = R[i,t] - R[i,t]%*%F.%*%t(F.)%*%R[i,t]*(1-Q.st[i,t]/q[i,t])/q[i,t]
    }
  }

  m1 = m[1:900,k0]
  m = media =  matrix(rep(0,900),nrow = 900,ncol=tot)
  m[,1] = m1

  c1 = c[1:900,k0]
  c =  matrix(rep(0,900),nrow = 900,ncol=tot)
  c[,1] = c1

  r1 = r.t1[1:900,k0]
  r.t1 =  matrix(rep(0,900),nrow = 900,ncol=tot)
  r.t1[,1] = r1

  s1 = s.t1[1:900,k0]
  s.t1 =  matrix(rep(0,900),nrow = 900,ncol=tot)
  s.t1[,1] = s1

  k0=k0+1
  for(t in k0:tot){
    y = countdf[[t]]$Freq
    for(i in 1:900){
      a[i,t] = G%*%m[i,t-k0+1]

      R[i,t] = (G%*%c[i,t-k0+1]%*%t(G)+sigVg[i,i])
      #R[i,t] = sigVg[i,i]

      #cat("a",a,"\n")
      #Passo 2 - Previsão a um passo
      f[i,t] = t(F.)%*%a[i,t]
      q[i,t] = t(F.)%*%R[i,t]%*%F.

      #com fator de descontos
      r.t1[i,t-k0+2] = delta[j]*(r.t1[i,t-k0+1]+y[i])+1-delta[j]
      s.t1[i,t-k0+2] = delta[j]*(s.t1[i,t-k0+1]+1)
      media[i,t] = r.t1[i,t-k0+2]/s.t1[i,t-k0+2]

      #Passo 3 -  Atualização
      f.st[i,t] = digamma(r.t1[i,t-k0+2]+y[i])-log(s.t1[i,t-k0+2]+1)
      Q.st[i,t] = trigamma(r.t1[i,t-k0+2]+y[i])

      m[i,t-k0+2]  = a[i,t] + R[i,t]%*%F.*(f.st[i,t] - f[i,t])/q[i,t]
      c[i,t-k0+2] = R[i,t] - R[i,t]%*%F.%*%t(F.)%*%R[i,t]*(1-Q.st[i,t]/q[i,t])/q[i,t]
    }

    gama=t(as.matrix(y.star1[[t]])-as.matrix(y.star1[[t-1]]))%*%M%*%w1%*%(as.matrix(y.star1[[t]])-as.matrix(y.star1[[t-1]])) /
      (t(as.matrix(y.star1[[t]])-as.matrix(y.star1[[t-1]]))%*%w1%*%M^2%*%w1%*%(as.matrix(y.star1[[t]])-as.matrix(y.star1[[t-1]])) +
         sum(rowSums(w1^2)/colSums(w1^2)*(y.star2[[t]]-y.star1[[t]]^2)))
    #
    sig=(1/n) * (t(as.matrix(y.star1[[t]])-as.matrix(y.star1[[t-1]]))%*%(solve(M)-as.vector(gama)*w1)%*%(as.matrix(y.star1[[t]])-
                                                                                                           as.matrix(y.star1[[t-1]]))+sum(colSums(w1)*(y.star1[[t]]-y.star2[[t]]^2)))

    v_g = solve(diag(1,900,900)-as.vector(gama)*M%*%w1)%*%M
    sigVg = as.vector(sig)*v_g

  }

  } else{
    data("mRou")
    data("cRou")
    data("rt1Rou")
    data("st1Rou")
    data("SigVgRou")
    j = 1
    G = matrix(1)
    F.= matrix(1)
    a = r = f1 = q = s = r = f.st = R = mm = pred1 =  matrix(rep(0,900),nrow = 900,ncol=tot)
    m = matrix(mRou$V1,nrow = 900,ncol=tot)
    c = matrix(cRou$V1,nrow = 900,ncol=tot)
    f.st = q.st = Q.st  = s.t1 = r.t1 =  matrix(rep(0,900),nrow = 900,ncol=tot)

    r.t1 = matrix(rt1Rou$V1,nrow = 900,ncol=tot)
    s.t1 = matrix(st1Rou$V1,nrow = 900,ncol=tot)

    f=matrix(rep(0,901),nrow = 901,ncol=tot)
    media = matrix(rep(0,901),nrow = 901,ncol=tot)
    #t=k0
    n=900
    sigVg = sigVgRou

    for(t in k0:tot){
      y = countdf[[t]]$Freq
      for(i in 1:900){
        a[i,t] = G%*%m[i,t-k0+1]

        R[i,t] = (G%*%c[i,t-k0+1]%*%t(G)+sigVg[i,i])
        #R[i,t] = sigVg[i,i]

        #cat("a",a,"\n")
        #Passo 2 - Previsão a um passo
        f[i,t] = t(F.)%*%a[i,t]
        q[i,t] = t(F.)%*%R[i,t]%*%F.

        #com fator de descontos
        r.t1[i,t-k0+2] = delta[j]*(r.t1[i,t-k0+1]+y[i])+1-delta[j]
        s.t1[i,t-k0+2] = delta[j]*(s.t1[i,t-k0+1]+1)
        media[i,t] = r.t1[i,t-k0+2]/s.t1[i,t-k0+2]

        #Passo 3 -  Atualização
        f.st[i,t] = digamma(r.t1[i,t-k0+2]+y[i])-log(s.t1[i,t-k0+2]+1)
        Q.st[i,t] = trigamma(r.t1[i,t-k0+2]+y[i])

        m[i,t-k0+2]  = a[i,t] + R[i,t]%*%F.*(f.st[i,t] - f[i,t])/q[i,t]
        c[i,t-k0+2] = R[i,t] - R[i,t]%*%F.%*%t(F.)%*%R[i,t]*(1-Q.st[i,t]/q[i,t])/q[i,t]
      }

      if(t>=2){

      gama=t(as.matrix(y.star1[[t]])-as.matrix(y.star1[[t-1]]))%*%M%*%w1%*%(as.matrix(y.star1[[t]])-as.matrix(y.star1[[t-1]])) /
        (t(as.matrix(y.star1[[t]])-as.matrix(y.star1[[t-1]]))%*%w1%*%M^2%*%w1%*%(as.matrix(y.star1[[t]])-as.matrix(y.star1[[t-1]])) +
           sum(rowSums(w1^2)/colSums(w1^2)*(y.star2[[t]]-y.star1[[t]]^2)))
      #
      sig=(1/n) * (t(as.matrix(y.star1[[t]])-as.matrix(y.star1[[t-1]]))%*%(solve(M)-as.vector(gama)*w1)%*%(as.matrix(y.star1[[t]])-
                                                                                                             as.matrix(y.star1[[t-1]]))+sum(colSums(w1)*(y.star1[[t]]-y.star2[[t]]^2)))

      v_g = solve(diag(1,900,900)-as.vector(gama)*M%*%w1)%*%M
      sigVg = as.vector(sig)*v_g
      }

    }

}

  media=media[,k0:tot]
  Ajuste=media[,-1]
  Pobs=Ppred=erro=NULL
  for(k in 1:(tot-k0)){
    z=cbind(countdf[[k+k0]]$Freq,round(media[1:900,k],0));z#periodo previstp
    erro[k]=sum((z[,2]-z[,1])^2)
    Pobs[k]=sum(z[,1]) #numero de pontos obs
    Ppred[k]=sum(round(z[,2])) #numero de pontos preditos
  }
  errotot=sum(erro)/((tot-k0)*900)

  k=(tot-k0)+1
  for(k in (k-nMap):(k)){

    if((k+k0)>tot){
      z=cbind(rep(NA,900),round(media[1:900,k],0))#periodo previstp

    }else{
      z=cbind(countdf[[k+k0]]$Freq,round(media[1:900,k],0))#periodo previstp
    }

    lambda=z[,2]

    #Obter as coordenadas de cada celula
    i=1
    px=py=list()
    d=list()
    for(i in 1:900){
      d[[i]]=unlist(tt$tiles[i]) #obtem o max e o minimo dos eixos x e y de cada celula
      if(lambda[[i]]==0) next
      seqx=seq(as.numeric(d[[i]][2]),as.numeric(d[[i]][3]))
      seqy=seq(as.numeric(d[[i]][4]),as.numeric(d[[i]][5]))
      px[[i]] = sample(seqx,size = lambda[i],prob=rep(lambda[i]/sum(lambda),length(seqx)))
      py[[i]] = sample(seqy,size = lambda[i],prob=rep(lambda[i]/sum(lambda),length(seqy)))

    }

    pontos = data.frame(x=unlist(px),y=unlist(py))


    loc2 <- ppp(pontos[,1],pontos[,2], window=as.owin(mao))
    cq <- quadratcount(loc2, 30, 30)


    if((k+k0)>tot){

      # plot(quadratcount(spp[[k+k0]],30,30))
      # plot(cq)

      ##Intensidade Observada
      #plot(intensity(grid1[[k+k0]], image=T), main=paste("y","_",k+k0,sep=""), cex.main=0.6)

      ## Intensidade estimada para o número de casos esperados
      par(mfrow=c(1,1))
      plot(intensity(cq, image=T), main=paste("y_pred",k+k0,sep=""), cex.main=0.6)
    }else{
      par(mfrow=c(1,2))

      plot(intensity(grid1[[k+k0]], image=T), main=paste("y","_",k+k0,sep=""), cex.main=0.6)

      ## Intensidade estimada para o número de casos esperados
      plot(intensity(cq, image=T), main=paste("y_pred",k+k0,sep=""), cex.main=0.6)
      #plot(intensity(cq, image=T), main=paste("y_pred",k+k0,sep=""), cex.main=0.6)
    }

    x=readline(prompt="Press [enter] to see the next plot:")
    if(x==""){
      next
    }else{
      break
    }

  }
  return(list(erro=errotot))

}
