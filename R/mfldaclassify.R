#' @title Classify Test Data Using Multivariate Linear Discriminant Analysis
#'
#' @description This function classifies test data based on learned sparse functional linear discriminant analysis parameters (features discriminant scores). It projects both training and test data onto the discriminant space and predicts the class labels for the test samples.
#'
#' @param hatalpha A matrix or list of matrices containing the estimated coefficients for the linear discriminants obtained from the training data.
#' @param Xtestdata A data frame containing the test data, structured similarly to the training data, with identifiers in the first three columns and features in the subsequent columns.
#' @param Xdata A data frame containing the training data, with the same structure as `Xtestdata`.
#' @param Y A vector of class labels for the training data. It should match the length of the number of unique samples in `Xdata`.
#' @param standardize A logical value indicating whether to standardize the data to have mean zero and variance one for each time point and each variable before classification. Default is TRUE.
#'
#' @return A list containing:
#'   \item{PredictedClass}{A vector of predicted class labels for the test data.}
#'   \item{Predclass.df}{A data frame of predicted class labels for each test sample across different dimensions of the discriminant space.}
#'
#' @details
#' The function performs the following steps:
#' 1. Standardizes the training and test data if the `standardize` parameter is set to TRUE.
#' 2. Projects both the training and test data into the discriminant space using the estimated coefficients.
#' 3. Computes the Euclidean distances between projected test data and the class means in the discriminant space.
#' 4. Assigns class labels to the test samples based on the closest mean in the discriminant space.
#'
#' @examples
#' # Assuming hatalpha contains the discriminant scores, Xtestdata is your test set,
#' # Xdata is your training set, and Y contains the corresponding training labels.
#' classification_result <- mfldaclassify(hatalpha = hatalpha, Xtestdata = Xtestdata, Xdata = Xdata, Y = trainY, standardize = TRUE)
#'
#' # Access predicted class labels
#' predicted_classes <- classification_result$PredictedClass
#'
#' # Access detailed prediction data frame
#' predictions_df <- classification_result$Predclass.df
#'
#' @export
mfldaclassify=function(hatalpha=hatalpha,Xtestdata=Xtestdata,Xdata=Xdata,Y=trainY,standardize=TRUE){

  XdataOrig=Xdata
  XtestdataOrig=Xtestdata
  nk=length(unique(XdataOrig$group))
  nTime=length(unique(XdataOrig$time))
  nsample = length(unique(XdataOrig$id))
  ntest = length(unique(XtestdataOrig$id))
  Xstand=list()
  Xteststand=list()

  #View(XtestdataOrig)

  standardize=TRUE
  if(is.null(standardize)){
    standardize=TRUE
  }

  if(standardize==TRUE){
    for(j in 1:nTime){
      myX=scale(as.matrix(XdataOrig[XdataOrig$time==j,-c(1:3)]), center=TRUE,scale=TRUE)
      Xstand[[j]]=cbind.data.frame(as.matrix(XdataOrig[XdataOrig$time==j,c(1:3)]), myX)
    }
    Xdata=do.call(rbind.data.frame,Xstand)
  }

  if((standardize==TRUE) & (ntest > 1)){
    for(j in 1:nTime){
      myXtest=scale(as.matrix(XtestdataOrig[XtestdataOrig$time==j,-c(1:3)]), center=TRUE,scale=TRUE)
      Xteststand[[j]]=cbind.data.frame(as.matrix(XtestdataOrig[XtestdataOrig$time==j,c(1:3)]), myXtest)
    }
    Xtestdata=do.call(rbind.data.frame,Xteststand)
  } else if (ntest == 1){
    Xtestdata = XtestdataOrig
  }
  #View(Xtestdata)

  if (nk == 2){
    Projtest=list()
    Projtrain=list()
    for(j in 1:nTime){
      Projtest[[j]]=as.matrix(Xtestdata[Xtestdata$time==j,-c(1:3)])%*%hatalpha[,j]
      Projtrain[[j]]=as.matrix(Xdata[Xdata$time==j,-c(1:3)])%*%hatalpha[,j]
    }
    nc= length(unique(as.vector(Y)))
    ntest=dim(Projtest[[1]])[1]
    
  } else if (nk == 3){
    
    ## LD1
    Projtest1=list()
    Projtrain1=list()
    for(j in 1:nTime){
      Projtest1[[j]]=as.matrix(Xtestdata[Xtestdata$time==j,-c(1:3)])%*%hatalpha[[1]][,j]
      Projtrain1[[j]]=as.matrix(Xdata[Xdata$time==j,-c(1:3)])%*%hatalpha[[1]][,j]
    }
    ##LD2
    Projtest2=list()
    Projtrain2=list()
    for(j in 1:nTime){
      Projtest2[[j]]=as.matrix(Xtestdata[Xtestdata$time==j,-c(1:3)])%*%hatalpha[[2]][,j]
      Projtrain2[[j]]=as.matrix(Xdata[Xdata$time==j,-c(1:3)])%*%hatalpha[[2]][,j]
    }
    nc= length(unique(as.vector(Y)))
    ntest=dim(Projtest1[[1]])[1]
  }
  PredclassSeparate=list()

  if (nk == 3){
    for(d in 1:nTime){ #nTime
      
      # LD1
      ProjXtestdatad1=Projtest1[[d]]
      
      ProjXdatad1=Projtrain1[[d]]
      ProjXdata1d1=cbind(Y,ProjXdatad1)
      Projmv1=aggregate(ProjXdata1d1[,-1],list(ProjXdata1d1[,1]),mean)
      
      # LD2
      ProjXtestdatad2=Projtest2[[d]]
      
      ProjXdatad2=Projtrain2[[d]]
      ProjXdata1d2=cbind(Y,ProjXdatad2)
      Projmv2=aggregate(ProjXdata1d2[,-1],list(ProjXdata1d2[,1]),mean)
      
      ProjXtestdatad12 = cbind(ProjXtestdatad1,ProjXtestdatad2)
      Projmv <- cbind(Projmv1[,2], Projmv2[,2])
      
      distv=list()
      for(j in 1: nc){
        #euclidean distance
        dist.j <- c()
        for (r in 1:ntest){
          dist.j[r] <- sqrt(sum((ProjXtestdatad12[r,] - Projmv[j,])^2))
        }
        distv[[j]]=dist.j
      }
      distv=do.call(cbind,distv)
      
      #The following code outputs the assigned class
      predclassX1 <- apply(distv, 1, which.min)
      PredclassSeparate[[d]]=predclassX1
    }
    
    Predclass=PredclassSeparate
    Predclass.df=do.call(rbind.data.frame,PredclassSeparate)
    
  } else if (nk == 2){
    for(d in 1:nTime){ #nTime
      ProjXtestdatad=Projtest[[d]]
      ProjXdatad=Projtrain[[d]]
      ProjXdata1d=cbind(Y,ProjXdatad)
      
      Projmv=aggregate(ProjXdata1d[,-1],list(ProjXdata1d[,1]),mean)
      
      distv=list()
      jrep=list()
      for(j in 1: nc){
        rProjm=matrix( rep(Projmv[j,-1],times= ntest), ncol=ncol(ProjXdatad), byrow=TRUE)
        #euclidean distance
        sqdiff=(ProjXtestdatad-as.numeric(rProjm))^2
        dist1=rowSums(sqdiff)^0.5
        jrep[[j]]=j*rep(1,times=ntest)
        distv[[j]]=dist1
      }
      distv=do.call(cbind,distv)
      dim(distv)=c(nc*nrow(distv),1)
      
      jrep=do.call(cbind,jrep)
      dim(jrep)=c(nc*nrow(jrep),1)
      
      distv=cbind(jrep, distv)
      rdistvX1=matrix(distv[,-1], nrow=ntest,ncol=nc)
      minX1=apply(rdistvX1,1,min)
      minind=which(rdistvX1==minX1, arr.ind=TRUE) #minimum indices
      predclassX1=minind[order(minind[,1]),2]
      PredclassSeparate[[d]]=predclassX1
    }
    
    Predclass=PredclassSeparate
    Predclass.df=do.call(rbind.data.frame,PredclassSeparate)
  }
  
  
  class1 <- c()
  class2 <- c()
  class3 <- c()
  
  for (j in 1: ntest){
    class1[j] <- sum(Predclass.df[,j]==1)
    class2[j] <- sum(Predclass.df[,j]==2)
    if (nk == 3){
      class3[j] <- sum(Predclass.df[,j]==3)
    }
  }
  
  #print(class1)
  #print(class2)
  #print(class3)
  
  if (nk == 2){
    classified.class <- ifelse(class1 > class2, 1, ifelse(class1 == class2, sample(1:2, length(class1), replace = TRUE), 2))
    #print(classified.class)
  } else if (nk == 3){
    classified.class <- ifelse( ((class1 > class2) & (class1 > class3)), 1,
                                ifelse((class2 > class3) & (class2 > class1), 2,
                                       ifelse((class3 > class1) & (class3 > class2), 3,
                                              ifelse((class1 == class2), sample(1:2, length(class1), replace = TRUE),
                                                     ifelse((class1 == class3), sample(c(1,3), length(class1), replace = TRUE),
                                                            ifelse((class3 == class3), sample(c(2,3), length(class2), replace = TRUE),
                                                                   sample(1:3, length(class1), replace = TRUE)))))))
    
    #print(classified.class)
  }
  
  result=list(PredictedClass=classified.class,Predclass.df = Predclass.df)
  return(result)
}
