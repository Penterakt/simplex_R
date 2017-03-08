
library(nnet)

pivot<-function(A1,b1,cst1,wcol1,TS){
  
  T = matrix(unlist(TS[1]),nrow(A1)    )
  z = unlist(TS[2])
  T0= unlist(TS[3])
  z0= unlist(TS[4])
  x = unlist(TS[5])
  
  s = which.is.max(z-cst1)
  
  v=T[,s]
  v1=which(v>0)
  
  a = T0[v1]/v[v1]
 r =v1[ which.is.max(-a) ]
  

  return(list(s,r))
}


SympTab<-function(A1,b1,cst1, wcol1){
  
  
  #macierz kolumn rozpinajacych
  
  
  Ak = A1[,wcol1]
  
  
  
  T = round(solve(Ak)%*%A1,digits =10)
  z = t(T)%*% cst1[wcol1]
  T0=solve(Ak)%*%b1
  z0 = cst1[wcol1]%*%T0
  
  
  
  j=1
  x=rep(0,ncol(A1))
  for (i in wcol1)
  {
    
    x[i] = T0[j] 
    j=j+1
  }
  
  
  return(list(T,z,T0,z0,x))


}


Modyfikacja<-function(A1,b1,cst1,wcol1,TS,Piv){
  
      T = matrix(unlist(TS[1]),nrow(A1)    )
      z = unlist(TS[2])
      T0= unlist(TS[3])
      z0= unlist(TS[4])
      x = unlist(TS[5])
  
      s=unlist(Piv[1])
      r=unlist(Piv[2])
  
  z0 = z0 - (z[s]-cst1[s])*T0[r]/T[r,s]
  
       w=c(1:nrow(A1))[-r]
  
  for(i in w){
    T0[i] = T0[i] - (1/T[r,s])*T[i,s]*T0[r] 
                           }
  
    T0[r] = (1/T[r,s])*T0[r]
  
  
  z = z - ((z[s]-cst1[s])/T[r,s])*T[r,]
  
  
  for(i in c(1:nrow(A1))[-r]) {
    T[i,] = T[i,] - (T[i,s]/T[r,s])*T[r,]
  }
  
  T[r,] = (1/T[r,s])*T[r,]
  
  
  wcol1[r]=s
  
  j=1
  x=rep(0,ncol(A1))
  for (i in  wcol1 )
  {
    
    x[i] = T0[j] 
    j=j+1
  }
  
  
  
  
  
  
  
  
  
  
  
  return(list(T,z,T0,z0,x,wcol1))
  
  
}



codalej<-function(A1,b1,cst1,wcol1,TS){
  
  T = matrix(unlist(TS[1]),nrow(A1)    )
  z = unlist(TS[2])
  T0= unlist(TS[3])
  z0= unlist(TS[4])
  x = unlist(TS[5])
  wcol1 = unlist(TS[6])
  
  stan=1
  v=z-cst1

  if(any(v>1.000000e-12)){  
                 s = which(z-cst1>0)[1]
                          
                   if(all(T[,s]<=0)){stan=2 
                                     s=0}
  }else{stan=3  
         s=0}
  

   return(stan)
}


SympleksI<-function(A1,b1,cst1,wcol1){
  
  TS=SympTab(A1,b1,cst1, wcol1)
     
      
    s00=codalej(A1,b1,cst1,wcol1,TS)
  while(s00==1){
    
    Piv=pivot(A1,b1,cst1,wcol1,TS)
   
    
    TS=Modyfikacja(A1,b1,cst1,wcol1,TS,Piv)
    
    
    T = matrix(unlist(TS[1]),nrow(A1)    )
    z = unlist(TS[2])
    T0= unlist(TS[3])
    z0= unlist(TS[4])
    x = unlist(TS[5])
    wcol1 = unlist(TS[6])
    
    s00=codalej(A1,b1,cst1,wcol1,TS)
    
    
  }
  s00
  
  T = matrix(unlist(TS[1]),nrow(A1)    )
  z = unlist(TS[2])
  T0= unlist(TS[3])
  z0= unlist(TS[4])
  x = unlist(TS[5])
  
  
  
  
  return(list(T,z,T0,z0,x,s00,wcol1))
}




# Dwukrokowy sympleks dla problemu minimum w formie kanonicznej


       #Krok I
DwuSympleksI<-function(A2,b2,c2){
  
  
  A1 = matrix(c(A2,diag(nrow(A2))),nrow(A2))
  cst1 = c(rep(0,ncol(A2)), rep(1,nrow(A2)) )
  b1= b2
  wcol1 = (ncol(A2)+1):(ncol(A2)+nrow(A2))
  
  #co jesli T0 niedodatnie
  TS = SympleksI(A1,b1,cst1,wcol1)
       
  
  
  T = matrix(unlist(TS[1]),nrow(A1)    )
  z = unlist(TS[2])
  T0= unlist(TS[3])
  z0= unlist(TS[4])
  x = unlist(TS[5])
  s00 = unlist(TS[6])
  wcol1 = unlist(TS[7])

  return(list(wcol1,s00) )
  }  
  
  
  DwuSympleksII <- function(A2,b2,c2,wcol1){
    A1 =    A2
    cst1 =  c2
    b1=     b2
    
    wcol1 = unlist(DwuSympleksI(A2,b2,c2)[1] )
    
    
    TS = DwuSympleksI(A1,b1,cst1, wcol1)
    T = matrix(unlist(TS[1]),nrow(A1)    )
    z = unlist(TS[2])
    T0= unlist(TS[3])
    z0= unlist(TS[4])
    x = unlist(TS[5])
    
    return(list(T,z,T0,z0,x))
  }
  
  
     
     
Kanonizacja<-function(A1,b1,cst1,I1,J1){
  n=ncol(A1)
  m=nrow(A1)
  M = A1
  
    #dokladamy Zmienne niedoboru
  cc=cst1
  
  if(is.null(I1)==0){
  for(i in I1){
    v = rep(0,m)
    v[i]= -1
    M=matrix( c(M,v),m    )
    cc=c(cc,0)
  }
                    }
    #dokladamy rownania typu x = u-v
    
    bb=b1
    
    
    if(is.null(J1)==0){
  for(j in J1){
    v = rep(0,ncol(M))
    v[j]=1
  M=rbind(M,v)
    v = rep(0,nrow(M))
    v[nrow(M)]=-1
  
    
  M=cbind(M,v) 
  cc=c(cc,0)
    
    v = rep(0,nrow(M))
    v[nrow(M)]=1
  
  M=cbind(M,v)  
    bb =c(bb,0)
    cc=c(cc,0)
  }

}
    
        
    
    X = list(M,bb,cc)
  
  
    return(X)
  }
  

trzask<-function(A2,b2,cst2,I2,J2){
  n=ncol(A2)
  
  A1 = A2
  b1 = b2
  cst1 = cst2
  I1 = I2
  J1 = J2
  K=Kanonizacja(A1,b1,cst1,I1,J1)
  
  #K--->(A1,b1,cst1)
  
  #Szukamy dobrego wyboru kolumn 
  A2 = matrix(  unlist(K[1]),nrow(A1)+length(J2) )
  b2 = unlist(K[2])
  c2 = unlist(K[3])
    #zeby b>0
   for(i in which(b2<1.0e-12)){A2[i,]=-A2[i,]
       b2[i]=-b2[i]}
  
  L=DwuSympleksI(A2,b2,c2)
  #DwuSympleksI--->(wcol1 , s00)
  
  wcol1 = unlist(L[1])

  
  if(max(wcol1) < ncol(A2)+0.1){ 
    
    
    
  X = SympleksI(A2,b2,c2,wcol1)
  #SympleksI-->(T,z,T0,z0,x,s00,wcol1)
  
  
  x = unlist(X[5])[1:n]
  cost= unlist(X[4])
  s00 = unlist(X[6])  
    
  
  if(s00==2){x=NULL
             cost=c("dowolnie maly")}
  
  
  
  }else{
          x=NULL
          cost = c("nie ma rozw dopuszczalnego")
          }
    
              
  
  
  
  
    
    

  
  
  
  
  return(list(x,cost))  
}



