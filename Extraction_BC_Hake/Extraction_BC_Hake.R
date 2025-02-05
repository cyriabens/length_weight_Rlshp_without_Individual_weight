#########################################################################################################################
##### Trying to fit length-weight relationship when we have size distribution ~ Total Weight     ########################
##### Extraction of body condition for hake #############################################################################
##### Coded by Bensebaini Meriem Cyria & Gregoire Certain ###############################################################
#########################################################################################################################
#########################################################################################################################
InputDir<-'~/Cyria/Article 1 correction 3/Ready to submission/Extraction_BC_Hake/Input/'

###
### Load packages
###

library(dplyr)
library(parallel)

###
### Load data
###
load(file.path(paste(InputDir,"cubeN.hake.Rdata",sep="")))
load(file.path(paste(InputDir,"cubeB.hake.Rdata",sep="")))
load(file.path(paste(InputDir,"vec.Rdata",sep="")))
### Make cluster
cl <- makeCluster(getOption("cl.cores", 4))

########## The function that look for the best a and b 
getWfromL<-function(ab,vecN,ID_cl=as.numeric(names(vecN))){
  a<-ab[1]
  b<-ab[2]
  return(sum((a*ID_cl^b)*vecN))
}

########## fit length-weight relationship with  least squares method using the getWfromL function
LWresid.fct<-function(ab,NL,B){
  return(sum((parApply(cl=cl,NL,1,getWfromL,ab=ab)-B)^2))}


########## the function that calculates the residual index
LWcond.fct<-function(ab,NL,B){
  return(B-(parApply(cl=cl,NL,1,getWfromL,ab=ab)))}


########## function to randomly select 40 hauls 
boot.by.y.fct<-function(vec){
  row.sel<-NULL
  for(i in unique(vec)){
    row.sel<-c(row.sel,sort(sample(((1:length(vec))[vec==i]),40,replace=F)))
  }
  vec.sel<-rep(F,length(vec))
  vec.sel[row.sel]<-T
  return(vec.sel)
}

### choose nboot and create save object:
nboot<-100
CI_allspecies<-array(NA,dim=c(1,25,nboot),dimnames=list('Hake',unique(vec),1:nboot))

### start the loop for the bootstarp 

  for (n in 1:nboot){ ## PLEASE note that it takes a little time to run !!!!!!
    new.Y<-0          ## For more speed reduce the number of nboot !!!!!
    nmintrawl<-0
    cond<-F
    nwhile<-0
    while(cond==F){
      vecsel<-boot.by.y.fct(vec)
      NL.esp<-cubeN.hake[vecsel,]
      B.esp<-apply(cubeB.hake[vecsel,],1,sum)
      FullTrawl<-apply(NL.esp,1,sum)>0
      new.Y<-vec[vecsel][FullTrawl]
      nmintrawl<-min(table(new.Y))
      nwhile<-nwhile+1
      cond<-(length(unique(new.Y))==length(unique(vec))) & (nmintrawl>=3)
      if(nwhile>1000){cond<-T}
    }
    
    if(nwhile<=1000){
      NL.esp.full<-NL.esp[FullTrawl,]
      B.esp.full<-B.esp[FullTrawl]
      res<-nlminb(start=c(0.01,3), objective=LWresid.fct,NL=NL.esp.full,B=B.esp.full)
      essai<-LWcond.fct(ab=res$par,NL=NL.esp.full,B=B.esp.full)
      new.Y<-vec[vecsel][FullTrawl]
      
      ### Isolating outliers
      
      keep<-essai>quantile(essai,0.025) & essai<quantile(essai,0.975) 
      NL.esp.k<-NL.esp.full[keep,]
      B.esp.k<-B.esp.full[keep]
      res2<-nlminb(start=c(0.01,3), objective=LWresid.fct,NL=NL.esp.k,B=B.esp.k)
      essai2<-LWcond.fct(ab=res2$par,NL=NL.esp.k,B=B.esp.k)
      if (sum(keep)!=0){
        CI_allspecies['Hake',,n]<-tapply(essai2,new.Y[keep],mean)}
    }
  }  


## save data
Cube.boot.stock<-CI_allspecies

## Plot function
plot.ci.fct<-function(Cube,numspe){
  ci.inf<-apply(Cube[numspe,,],c(1),quantile,probs=0.025)
  ci.sup<-apply(Cube[numspe,,],c(1),quantile,probs=0.975)
  med<-apply(Cube[numspe,,],c(1),quantile,probs=0.5)
  plot(dimnames(Cube)[[2]], med,type='b',xlab='Year',ylab='RI'
       ,ylim=range(c(ci.inf,ci.sup)),
       pch=20,lwd=1.3, cex.lab=1.4, cex.axis= 1.4)
  points(dimnames(Cube)[[2]],ci.inf,type="l",lty=3,lwd= 1.3, col='darkgrey')
  points(dimnames(Cube)[[2]],ci.sup,type="l",lty=3,lwd= 1.3, col='darkgrey')
  abline(a=0,b=0,col="blue",lty=1, lwd=1.3)
  title(dimnames(Cube)[[1]][numspe], adj = 0.5, line = 0.3,cex= 1.5)
}

plot.ci.fct(Cube.boot.stock,numspe = 1)
