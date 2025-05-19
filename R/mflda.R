#' @title Multivariate Functional Linear Discriminant Analysis (MFLDA)
#'
#' @description This function performs Multivariate Functional Linear Discriminant Analysis (MFLDA) for classification of multi-class data. It trains the model on provided training data, optionally standardizes the data, and evaluates performance on test data.
#'
#' @param Xtrain A data frame containing the training data with features and identifiers. The first three columns should be identifiers, and subsequent columns should be features.
#' @param Y A vector of class labels for the training data, corresponding to the rows in `Xtrain`.
#' @param Tau Regularization parameter for sLDA; determines the sparsity of the resulting discriminant functions.
#' @param Xtestdata A data frame containing the test data, structured similarly to `Xtrain`. If NULL, training data is used for testing.
#' @param Ytest A vector of class labels for the test data. Required if `Xtestdata` is provided.
#' @param plotIt A logical indicating whether to plot the discriminants and densities. Default is FALSE.
#' @param standardize A logical indicating whether to standardize the training and test data standardized to have mean zero and variance one for each time point and each variable. Default is TRUE.
#' @param maxiteration Maximum number of iterations for convergence in optimization. Default is 20.
#' @param thresh Convergence threshold for the optimization algorithm. Default is 1e-03.
#'
#' @return A list containing:
#'   \item{AverageError}{Average error of the classification.}
#'   \item{Metrics}{A vector of classification performance metrics (accuracy, balanced accuracy, F1 score, precision, recall, and Matthews correlation coefficient).}
#'   \item{TimeSpecificError}{Error for the test data classified based on the provided time points.}
#'   \item{hatalpha}{Estimated coefficients from the sLDA model.}
#'   \item{PredictedClass}{Predicted class labels for the test data.}
#'   \item{myDiscPlot}{Plot object of discriminants and densities, if `plotIt` is TRUE.}
#'   \item{Projection.df}{Data frame of projected values for visualization.}
#'   \item{var.selected}{Indices of selected variables based on sparsity criteria.}
#'   \item{varname.selected}{Names of selected variables based on sparsity criteria.}
#'
#' @details
#' This function implements the MFLDA algorithm by performing the following steps:
#' 1. Checks and preprocesses the data.
#' 2. Optionally standardizes the data.
#' 3. Performs the iterative optimization to obtain the discriminant functions.
#' 4. Classifies the test data based on the learned discriminant functions.
#' 5. Computes various performance metrics and generates plots if requested.
#'
#' @examples
#' # Example usage of mflda function:
#' result <- mflda(Xtrain = train_data, Y = train_labels, Tau = 0.5,
#'                  Xtestdata = test_data, Ytest = test_labels,
#'                  plotIt = TRUE, standardize = TRUE)
#'
#' # Access predicted classes and performance metrics
#' predicted_classes <- result$PredictedClass
#' performance_metrics <- result$Metrics
#'
#' @export
mflda=function(Xtrain=Xtrain,Y=Y,Tau=Tau,Xtestdata=Xtestdata,Ytest=Ytest,
               plotIt=FALSE,standardize=TRUE,maxiteration=20,thresh= 1e-03){

  Xdata1=as.matrix(Xtrain[Xtrain$time==1,-c(1:3)])
  nk=length(unique(Xtrain$group))
  n=length(Y)
  Xdata=list(Xdata1)
  dsizes=lapply(Xdata, function(x) dim(x))
  p=lapply(Xdata, function(x) dim(x)[2])
  #D=length(dsizes)
  D=1

  #check data
  if (is.list(Xdata)) {
    D = length(Xdata)
  } else {
    stop("Input data should be a list")
  }

  #check if testing data are provided. If not, will set training data as testing data.
  if(is.null(Xtestdata)){
    Xtestdata=Xdata
    Ytest=Y
  }


  if(is.null(plotIt)){
    plotIt=FALSE
  }

  if(is.null(standardize)){
    standardize=TRUE
  }

  if(is.null(maxiteration)){
    maxiteration=20
  }


  if(is.null(thresh)){
    thresh=1e-03
  }

  #standardize if true
  #standardize if true
  Xstand=list()
  nTime=length(unique(Xtrain$time))
  if(standardize==TRUE){
    for(j in 1:nTime){
      myX=scale(as.matrix(Xtrain[Xtrain$time==j,-c(1:3)]), center=TRUE,scale=TRUE)
      Xstand[[j]]=cbind.data.frame(as.matrix(Xtrain[Xtrain$time==j,c(1:3)]), myX)
    }
    Xdata=do.call(rbind.data.frame,Xstand)
  }


  #View(Xdata)
  #print("Standarize completed")

  # if(standardize==TRUE){
  #    Xdata=lapply(Xdata,function(x)scale(x,center=TRUE,scale=TRUE))
  #    Xtestdata=lapply(Xtestdata,function(x)scale(x,center=TRUE,scale=TRUE))
  # }

  print("finish standardized")

  nK=length(unique(as.vector(Y))) -1

  #norm function for convergence
  normdiff=function(xnew,xold){
    ndiff=norm(xnew-xold,'f')^2 / norm (xold,'f')^2
  }
  #initialize
  iter=0
  diffalpha=1
  reldiff=1

  #print("finish here")
  mynsparse=myfastLDAnonsparse(Xdata,Y)
  myalpha=mynsparse$tildealphamat
  #myalpha=as.matrix(do.call(cbind,myalpha)) ## change to matrix
  
  if (nk == 3){
    ##############################
    myalpha.temp = list()
    myalpha.temp[[1]] = myalpha
    myalpha.temp[[2]] = myalpha
    myalpha = myalpha.temp
    ##############################
  }
  
  Ux=mynsparse$Ux
  sqrtminvmat=mynsparse$sqrtminvmat
  Sbrx=mynsparse$Sbx
  tildelambda=mynsparse$tildelambda

  ######
  SqrtmSwSbSqrtmSw=mynsparse$SqrtmSwSbSqrtmSw
  ####

  print("finish fastLDA")

  #while convergence is not met
  while(iter < maxiteration && min(reldiff,max(diffalpha))> thresh){
    iter=iter+1
    cat("current iteration is", iter, "\n")

    if (nk == 2){
      if (dim(myalpha)[2] != 1){
        myalphaold=matrix(as.vector(myalpha), ncol = 1)
      } else {
        myalphaold=myalpha
      }
    } else if (nk == 3){
      myalphaold=myalpha
    }
    

    if (nk == 2){
      mymflda=mfldainner(Xdata,Y,Ux,SqrtmSwSbSqrtmSw,myalphaold,tildelambda,Tau)
    } else if (nk == 3){
      mymflda=mfldainner(Xdata,Y,Ux,SqrtmSwSbSqrtmSw,myalphaold[[1]],tildelambda,Tau)
    }
    
    if (nk == 2){
      myalpha=mymflda$hatalpha[[1]][[1]] # matrix
      nz = colSums(myalpha!=0)
      if(sum(nz==0) > 0.1 * dim(myalpha)[2]){
        myalpha=myalphaold
        break
      }
      myalpha2 = myalpha
      myalphaold2 = myalphaold
    } else if (nk ==3){
      myalpha=mymflda$hatalpha[[1]] # matrix
      ##############################
      nz = colSums(myalpha[[1]]!=0)
      if(sum(nz==0) > 0.1 * dim(myalpha[[1]])[2]){
        myalpha=myalphaold
        break
      }
      myalpha2 = myalpha[[1]]
      myalphaold2 = myalphaold[[1]]
      ##############################
    }
    
    tryCatch(
      expr = {
        diffalpha=normdiff(myalpha2, myalphaold2)
        sumnormdiff=norm(myalpha2-myalphaold2,'f')^2
        sumnormold=norm(myalphaold2,'f')^2
        reldiff=sumnormdiff/sumnormold
      },
      error = function(e){
        reldiff = 999
      }
    )
  }
  print("finish convergence")


  myclassified <- rep(1, length(Ytest))
  mfldaerror = 1
  
  #classification
  tryCatch(
    expr = {
      myclassify=mfldaclassify(myalpha,Xtestdata,Xdata,Y)
      mfldaerror=sapply(1:nTime, function(x)  sum(unlist(myclassify$Predclass.df[x,])!=Ytest)/length(Ytest) )
      myclassified <- myclassify$PredictedClass
      print(myclassified)
    },
    error = function(e){
      # (Optional)
      # Do this if an error is caught...
      #View(myalpha)
      mfldaerror = 1
      myclassified <- sample(1:2, length(Ytest), replace = TRUE)
      print(myclassified)
    }

  )

  print(as.factor(myclassified))
  #print(as.factor(Ytest))

  Ytest <- as.factor(Ytest)
  
  lvs <- c("1", "2", "3")
  levels(Ytest) <- lvs[1:nk]
  
  myclassified <- as.factor(c(1:nk, myclassified))
  levels(myclassified) <- lvs[1:nk]
  myclassified <- myclassified[-1]
  myclassified <- myclassified[-1]
  if (nk == 3){
    myclassified <- myclassified[-1]
  }
  
  result.metrics <- confusionMatrix(data = myclassified, reference = Ytest)
  
  if (nk == 3){
    group.prop = table(Ytest)/length(Ytest)
    metrics.blacc <- sum(result.metrics$byClass[,"Balanced Accuracy"]* group.prop)
    metrics.f1 <- sum(result.metrics$byClass[,"F1"]* group.prop)
    metrics.prec <- sum(result.metrics$byClass[,"Precision"]* group.prop)
    metrics.recall <- sum(result.metrics$byClass[,"Recall"]* group.prop)
  } else if (nk ==2){
    metrics.blacc <- result.metrics$byClass["Balanced Accuracy"]
    metrics.f1 <- result.metrics$byClass["F1"]
    metrics.prec <- result.metrics$byClass["Precision"]
    metrics.recall <- result.metrics$byClass["Recall"]
  }
  metrics.acc <- result.metrics$overall["Accuracy"]
  metrics.mcc <- mcc(as.factor(myclassified), as.factor(Ytest))
  
  list.metrics <- c(metrics.acc, metrics.blacc, metrics.f1, metrics.prec, metrics.recall, metrics.mcc)
  print("finish classification")

  #Produce discriminant and correlation plot if plotIt=T
  if(plotIt==TRUE){
    if (nk == 2){
      myalpha_plot = myalpha
    } else if (nk == 3){
      myalpha_plot = myalpha[[1]]
    }
    myDiscPlot <- DiscriminantPlots(Xtestdata = Xtestdata,Ytest = Ytest, myalpha = myalpha_plot, predicted.class = myclassified)
  }else{
    myDiscPlot=NULL
  }



  #View(myDiscPlot$Projection.df)
  gc()

  if(nk == 2){
    hatalpha2 <- myalpha
    selectP.70 <- rowSums(hatalpha2 != 0) >= ncol(hatalpha2)*0.7
    var.name <- colnames(Xtrain)[-c(1:3)][selectP.70]
  } else if (nk == 3){
    hatalpha1 <- myalpha[[1]]
    hatalpha2 <- myalpha[[2]]
    selectP.70.1 <- rowSums(hatalpha1 != 0) >= ncol(hatalpha1)*0.7
    selectP.70.2 <- rowSums(hatalpha2 != 0) >= ncol(hatalpha2)*0.7
    selectP.70 = list(selectP.70.1, selectP.70.2)
    
    var.name.1 <- colnames(Xtrain)[-c(1:3)][selectP.70.1]
    var.name.2 <- colnames(Xtrain)[-c(1:3)][selectP.70.2]
    var.name = list(var.name.1, var.name.2)
  }
  
  
  result=list(AverageError=mean(mfldaerror),Metrics=sum(list.metrics), TimeSpecificError=mfldaerror,hatalpha=myalpha,PredictedClass=myclassified,
              myDiscPlot = myDiscPlot, Projection.df = myDiscPlot$Projection.df, var.selected = selectP.70, varname.selected = var.name)
  return(result)
}
