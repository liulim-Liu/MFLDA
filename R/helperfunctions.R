myfastLDAnonsparse=function(X, Y){
  Y=as.vector(Y)
  nc=max(unique(as.vector(Y)))
  nTime=length(unique(X$time))
  print("start fastLDA iteration")
  
  # reshape the dataset
  # Reshape to wide format: one row per id, columns for each time-variable combo
  X$group <- NULL
  XdataTime <- reshape(X, 
                     timevar = "time", 
                     idvar = "id", 
                     direction = "wide")
  XdataTime$id <- NULL
  
  myX=as.matrix(XdataTime)
  n=dim(myX)[1]
  p=dim(X)[2]-2
  mysvd=svd(t(myX))
  Ux1=mysvd$u
  V=mysvd$v
  W=diag(mysvd$d)
  R=W%*%t(V)
  
  rdata2=cbind(Y,t(R))
  rdata=t(rdata2)
  mrd=aggregate(rdata2[,-1],list(rdata2[,1]),mean)
  mr=rowMeans(rdata[-1,])
  
  
  nc=max(unique(Y))
  C=list()
  for(i in 1:nc)
  {
    C[[i]]=rdata2[rdata2[,1]==i,-1] - matrix(rep(t(mrd[mrd[,1]==i,-1]),times=sum(Y==i)) ,ncol=ncol(rdata2[,-1]),byrow=TRUE)
  }
  C=as.matrix(do.call(rbind,C))
  
  #Swrx=t(C)%*%C /(n-1)
  library(corpcor)
  Swrx = cov.shrink(C)
  
  #Crx=R-t(matrix(rep(mr,times=n),ncol=ncol(R),byrow=TRUE))
  Crx=R-rowMeans(R)
  Srx=Crx%*%t(Crx)/(n-1)
  
  Sbrx=Srx-Swrx
  #Sbrx=Sbrx + t(Sbrx)
  Sbx2=Sbrx
  
  lambda=sqrt(log(p)/n)
  if(n<p){
    Strx2=Swrx + lambda*diag(n)
    #Strx=Strx + t(Strx)
    
  } else if(n >= p){
    Strx2=Swrx
    #Strx=Strx + t(Strx)
  }
  
  sqrtminv2= mysqrtminv(Strx2)$sqrtminv;
  Ux2=Ux1;
  sqrtminvStrxSbrx=Ux2%*%sqrtminv2%*%Sbx2%*%sqrtminv2%*%t(Ux2) #uXSw^(-1/2)SbSw^(-1/2)uX^t
  #print(all.equal(sqrtminv2%*%Sbx2%*%sqrtminv2, t(sqrtminv2%*%Sbx2%*%sqrtminv2)))
  myeigen2=eigs_sym(sqrtminv2%*%Sbx2%*%sqrtminv2,1,which="LM")
  myalphaold2=Ux2%*%myeigen2$vectors
  tildealphamat=myalphaold2/norm(myalphaold2,'2')
  #tildealphamat <- matrix(tildealphamat, nrow = p, ncol = nTime, byrow = FALSE)
  tildelambda=myeigen2$values
  gc()
  
  
  print("end fastLDA iteration")
  result=list(tildealphamat=tildealphamat,tildelambda=tildelambda,sqrtminvmat=sqrtminv2,
              Sbx=Sbx2, Ux=Ux2, Swx=Swrx,SqrtmSwSbSqrtmSw=sqrtminvStrxSbrx);
  return(result)
  
}


myfastLDAnonsparse.ind=function(X, Y){


  #   %--------------------------------------------------------------------------
  #   %myfastLDAnonsparse.R: function to obtain nonsparse solution to sparse functional linear
  # discriminant problem
  # %and to obtain matrix needed in constraints
  # %--------------------------------------------------------------------------
  #
  # X is is in long format with a variable name time
  #first column is ID, second column is time, third column is group, fourth-end is
  #variable names
  # if (is.list(X)) {
  #   Xdata2 = X
  # } else {
  #   Xdata2=list(X)
  # }

  #D = length(Xdata2)

  Y=as.vector(Y)

  nc=max(unique(as.vector(Y)))

  Crxd=list()
  Sbx=list()
  myalphaold1=list()
  myalphaold2=list()
  myalphaoldmat=list()
  rmyalphaoldmat=list()
  sqrtminvmat=list()
  tildealphamat=list()
  tildelambda=list()
  Ux=list()
  Swx=list()
  myeigenT=list()
  sqrtminvStrxSbrx=list()

  Xdata=X
  nTime=length(unique(Xdata$time))
  #nTime=5

  print("start fastLDA iteration")

  for (j in 1:nTime){
    print(j)
    myX=as.matrix(Xdata[Xdata$time==j,-c(1:3)])
    #View(myX)
    Y=as.matrix(Xdata[Xdata$time==j,3])
    n=dim(myX)[1]
    p=dim(myX)[2]

    #myX=as.matrix(Xdata)
    # n=dim(myX)[1]
    # p=dim(myX)[2]
    mysvd=svd(t(myX));
    Ux[[j]]=mysvd$u;
    V=mysvd$v;
    W=diag(mysvd$d)
    R=W%*%t(V)


    rdata2=cbind(Y,t(R))
    rdata=t(rdata2)
    mrd=aggregate(rdata2[,-1],list(rdata2[,1]),mean)

    #print("finish here")

    mr=rowMeans(rdata[-1,])

    #print("finish here 2")

    nc=max(unique(Y))
    C=list()
    for(i in 1:nc)
    {
      C[[i]]=rdata2[rdata2[,1]==i,-1] - matrix(rep(t(mrd[mrd[,1]==i,-1]),times=sum(Y==i)) ,ncol=ncol(rdata2[,-1]),byrow=TRUE)
    }
    C=as.matrix(do.call(rbind,C))
    Swrx=t(C)%*%C /(n-1)

    #print("finish here 3")
    #Crx=R-t(matrix(rep(mr,times=n),ncol=ncol(R),byrow=TRUE))
    Crx=R-rowMeans(R)
    Srx=Crx%*%t(Crx)/(n-1)

    Sbrx=Srx-Swrx
    #Sbrx=Sbrx + t(Sbrx)
    Sbx[[j]]=Sbrx

    lambda=sqrt(log(p)/n)
    if(n<p){
      Strx=Swrx + lambda*diag(n)
      #Strx=Strx + t(Strx)

    } else if(n >= p){
      Strx=Swrx
      #Strx=Strx + t(Strx)
    }

    Swx[[j]]=Strx;


    sqrtminv= mysqrtminv(Strx)$sqrtminv;
    sqrtminvmat[[j]]=sqrtminv;

    #form matrix
    sqrtminvStrxSbrx[[j]]=Ux[[j]]%*%sqrtminv%*%Sbrx%*%sqrtminv%*%t(Ux[[j]]) #uXSw^(-1/2)SbSw^(-1/2)uX^t

    myeigenT[[j]]=eigs_sym(sqrtminv%*%Sbrx%*%sqrtminv,1,which="LM") #pick the first eigenvalue-vector pair

    myalphaold2[[j]]=Ux[[j]]%*%myeigenT[[j]]$vectors
    tildealphamat[[j]]=myalphaold2[[j]]/norm(myalphaold2[[j]],'2')
    tildelambda[[j]]=myeigenT[[j]]$values
    gc()

  }
  print("end fastLDA iteration")
  result=list(tildealphamat=tildealphamat,tildelambda=tildelambda,sqrtminvmat=sqrtminvmat,
              Sbx=Sbx, Ux=Ux, Swx=Swx,SqrtmSwSbSqrtmSw=sqrtminvStrxSbrx);
  return(result)
}


mysqrtminv=function(W){
  #W is symmetric, positive definite
  mysvd=svd(W);
  d=diag(mysvd$d^(-0.5))
  out=mysvd$u%*%d%*%t(mysvd$u)
  result=list(sqrtminv=out)
  return(result)
}

minv=function(X){
  mysvd=svd(X)
  if(length(mysvd$d)==1){
    d=diag(as.matrix(mysvd$d^(-1)))
  }else{
    d=diag(mysvd$d^(-1))
  }
  out=mysvd$v%*%d%*%t(mysvd$u)
  return(out)
}

mfldainner.ind = function(Xtrain,Y,Ux,SqrtmSwSbSqrtmSw,myalphaold, tildelambda,myTau){
#will be used later for time points
  # if (is.list(Xtrain)) {
  #   Xdata = Xtrain
  # } else {
  #   Xdata=list(Xtrain)
  # }
  #Xdata=list(Xtrain)


  # #check data
  # if (is.list(Xdata)) {
  #   D = length(Xdata)
  # } else {
  #   stop("Input data should be a list")
  # }

  print("start inner")
  if (is.list(myTau)) {
    Tau = myTau
  } else {
    Tau=list(myTau)
  }

  nK=length(unique(as.vector(Y))) -1


  tildelambdaT=do.call(rbind,tildelambda)

  print("here 1")

  ## myalphaold is a matrix

  # original: myalphaold=do.call(rbind,myalphaold)
  #View(myalphaold)
  #print(dim(myalphaold))

  #myalphaold <-  lapply(seq_len(ncol(myalphaold)), function(i) myalphaold[,i])
  #myalphaold <- t(t(unlist(myalphaold)))
  #View(myalphaold)
  #print(dim(myalphaold))

  Ux1=as.matrix(do.call(rbind,Ux))

  myhatalpha=list()
  myalphamat=list()

  nTime=length(unique(Xtrain$time))
  #nTime=5
  D=1

  Xdata1=as.matrix(Xtrain[Xtrain$time==1,-c(1:3)])
  p=dim(Xdata1)[2]

  #######
  #myalphaold <- matrix(myalphaold, ncol = 100, byrow = TRUE) #p[[1]][1]
  #myalphaold <- t(myalphaold)
  #View(myalphaold)
  #######

  print("here 2")

  for(d in 1:D){
    #------------------------------------
    # ##solve for SIDA directions
    # #M=1
    # Alphai=Variable(p,nTime) # M is time point
    # Alphai2=vec(Alphai) #vectorize
    # Objx=sum(norm2(Alphai,axis=1))
    #
    # print("here 2")
    # Separation=bdiag(SqrtmSwSbSqrtmSw)
    #
    # print("here 3")
    #
    # # #defining the constraints
    # #            lapply(1:p,tildelambda[[j]]*matrix(1,nrow=p,ncol=1)))
    # tildelambdaT=list()
    # for(j in 1:nTime){
    #   print(j)
    #
    #   Xdata1=as.matrix(Xtrain[Xtrain$time==j,-c(1:3)])
    #   p=dim(Xdata1)[2]
    #   tildelambdaT[[j]]=tildelambda[[j]]*matrix(1,nrow=p,ncol=1)
    # }
    # tildelambdaT2=do.call(rbind,tildelambdaT)
    # print("here 4")
    # gc()
    # constraints=list(norm_inf(sum_entries(abs(Separation %*% myalphaold - tildelambdaT2*Alphai2), axis=1 ))<= Tau[[d]])
    #                                           # as.matrix(Separation)
    # prob=Problem(Minimize(Objx),constraints)
    #
    # result=solve(prob,solver="ECOS")
    # print("here 5")
    #
    # alphai=result$getValue(Alphai)
    #
    # View(alphai)
    # #write.csv(alphai,"old-alphai.csv")

    #-----------------

    #------------------------------------ for (have alphai calculated in each time point)
    #-----------------------------------Alphai=Variable(p,1)
      ##solve for SIDA directions
      #M=1

      Alphai=list()
      print("here 2")

      Separation=SqrtmSwSbSqrtmSw
      tildelambdaT=list()

      for(j in 1:nTime){
        #print(j)
        Alphaij=Variable(p,1)
        Alphaij2=vec(Alphaij) #vectorize
        Objx=sum(norm2(Alphaij,axis=1))

        Xdata1=as.matrix(Xtrain[Xtrain$time==j,-c(1:3)])
        p=dim(Xdata1)[2]
        tildelambdaT[[j]]=tildelambda[[j]]*matrix(1,nrow=p,ncol=1)

        tryCatch(
          expr = {
            constraints=list(norm_inf(sum_entries(abs(Separation[[j]] %*% as.matrix(myalphaold[,j]) - tildelambdaT[[j]]*Alphaij2), axis=1 ))<= Tau[[d]])
            prob=Problem(Minimize(Objx),constraints)
            result=solve(prob,solver="ECOS")
            #View(result)

            alphaj=result$getValue(Alphaij)

            Alphai[[j]]=alphaj
          },
          error = function(e){
            # (Optional)
            # Do this if an error is caught...
            Alphai[[j]]=0
          },
          finally = {
            # (Optional)
            # Do this at the end before quitting the tryCatch structure...
            gc()
          }
        )
      }

      gc()
      print("here 5")

      #View(Alphai)
      #

      alphai = do.call(cbind,Alphai)
      #View(alphai)
      #write.csv(alphai,"new-alphai.csv")

    #----------------- end()
    #----------------- alphavec[[t]]=alphai

    #print(dim(alphai))


    # alphai[abs(alphai) <=10^-5]=0
    alphai[abs(alphai) <=10^-10]=0
    #alphai[abs(alphai) <=10^-8]=0
    #View(alphai)

    if((sum(sum(abs(alphai)))==0)){
      myalpha=alphai
    }else{
      #myalpha=alphai/norm(alphai,"2")
      numerator <- matrix(rep(t(sqrt(colSums(alphai*alphai))),times=p),ncol=ncol(alphai),byrow=TRUE)
      #View(numerator)
      numerator[,colSums(numerator) == 0] = 1
      #View(numerator)
      myalpha=alphai/numerator
    }

    #myalphamat[[ii]]=as.vector(myalpha)

  #hatalpha=matrix(unlist(myalphamat), ncol=length(myalphamat),byrow =FALSE)
  myhatalpha[[d]]=myalpha

  }
  print("end inner")
  gc()
  result=list(hatalpha=myhatalpha)
  return(result)
}

mfldainner = function(Xtrain,Y,Ux,SqrtmSwSbSqrtmSw,myalphaold, tildelambda,myTau){
  
  print("start inner")
  if (is.list(myTau)) {
    Tau = myTau
  } else {
    Tau=list(myTau)
  }
  
  nK=length(unique(as.vector(Y))) -1
  print("here 1")
  
  
  Ux1=Ux
  myhatalpha=list()
  myalphamat=list()
  
  nTime=length(unique(Xtrain$time))
  #nTime=5
  D=1
  
  Xdata1=as.matrix(Xtrain[Xtrain$time==1,-c(1:3)])
  p=dim(Xdata1)[2]
  print("here 2")
  
  for(d in 1:D){
    for(ii in 1:nK){
      print(paste("solveing ",ii, "-th LD"))
      
      if(ii==2){
        alphamat = myalphamat[[1]]
        Xnew = data.frame()
        
        for(j in 1:nTime){
          Xmeta=as.matrix(Xtrain[Xtrain$time==j,1:3])
          Xdata1=as.matrix(Xtrain[Xtrain$time==j,-c(1:3)])
          p=dim(Xdata1)[2]
          
          ProjmX=Xdata1%*%(as.matrix(alphamat)%*%minv(t(as.matrix(alphamat))%*%as.matrix(alphamat)+0.001)%*%t(as.matrix(alphamat))) #0.001*diag(ii-1)
          ProjmX[is.nan(ProjmX)]=0
          Xn1=Xdata1-ProjmX
          Xn=cbind.data.frame(Xmeta,Xn1)
          Xnew = rbind(Xnew, Xn)
        }
        myfastlda=myfastLDAnonsparse(Xnew, Y)
        Ux1=myfastlda$Ux
        SqrtmSwSbSqrtmSw=myfastlda$SqrtmSwSbSqrtmSw
        
        tildelambda = myfastlda$tildelambda
        myalphaold=myfastlda$tildealphamat
        print("succuess solving")
      }
      
      gc()
      Alphai=Variable(p,nTime) 
      Alphai2=vec(Alphai) 
      Objx=sum(norm2(Alphai,axis=1))
      
      Separation=SqrtmSwSbSqrtmSw #PT x PT
      
      tildelambdaT=list()
      for(j in 1:nTime){
        Xdata1=as.matrix(Xtrain[Xtrain$time==j,-c(1:3)])
        p=dim(Xdata1)[2]
        tildelambdaT[[j]]=tildelambda*matrix(1,nrow=p,ncol=1)
        gc()
      }
      tildelambdaT2=do.call(rbind,tildelambdaT)
      gc()
      
      print("inner: solve for constraints")
      tryCatch(
        expr = {
          constraints=list(norm_inf(sum_entries(abs(Separation%*%myalphaold - tildelambdaT2*Alphai2), axis=1 ))<= Tau[[d]])
          prob=Problem(Minimize(Objx),constraints)
          result=solve(prob,solver="ECOS")
          alphai=result$getValue(Alphai) #Px T
          gc()
        },
        error = function(e){
          # (Optional)
          # Do this if an error is caught...
          alphai<<-matrix(0, nrow = p, ncol = nTime)
          #View(alphai)
          #print("catch a error")
        },
        finally = {
          # (Optional)
          # Do this at the end before quitting the tryCatch structure...
          gc()
        }
      )
      
      
      # alphai[abs(alphai) <=10^-5]=0
      alphai[abs(alphai) <=10^-10]=0
      #alphai[abs(alphai) <=10^-8]=0
      #View(alphai)
      
      if((sum(sum(abs(alphai)))==0)){
        myalpha=alphai
      }else{
        #myalpha=alphai/norm(alphai,"2")
        numerator <- matrix(rep(t(sqrt(colSums(alphai*alphai))),times=p),ncol=ncol(alphai),byrow=TRUE)
        #View(numerator)
        numerator[,colSums(numerator) == 0] = 1
        #View(numerator)
        myalpha=alphai/numerator
      }
      
      myalphamat[[ii]]=myalpha
      }
    #myalphamat[[ii]]=as.vector(myalpha)
    #hatalpha=matrix(unlist(myalphamat), ncol=length(myalphamat),byrow =FALSE)
    myhatalpha[[d]]=myalphamat
  }
  print("end inner")
  gc()
  result=list(hatalpha=myhatalpha)
  return(result)
}


