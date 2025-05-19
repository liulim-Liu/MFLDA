#' @title Cross-Validation for Multivariate Functional Linear Discriminant Analysis (MFLDA)
#'
#' @description This function performs cross-validation for Multivariate Functional Linear Discriminant Analysis (MFLDA).
#' It allows for optional parallel processing and can return various metrics, including accuracy or combined metrics.
#'
#' @param Xdata A data frame containing the training data in long format, with the first three columns for `ID`, `Time`, and `Group`, followed by variable columns.
#' @param Y A vector representing the group labels for the training data.
#' @param plotIt A logical value indicating whether to generate discriminant plots. Default is FALSE.
#' @param metrics.choice A string specifying the metric to optimize. Options are "Accuracy" or "CombinedMetrics". Default is "Accuracy".
#' @param Xtestdata A data frame containing the test data. If NULL, training data is used for testing.
#' @param Ytest A vector representing the group labels for the test data. Required if Xtestdata is provided.
#' @param isParallel A logical value indicating whether to perform parallel processing. Default is TRUE.
#' @param ncores An integer specifying the number of cores to use for parallel processing. Default is NULL, which uses half the available cores.
#' @param nfolds An integer specifying the number of folds for cross-validation. Default is 5.
#' @param ngrid An integer specifying the number of tuning grid values. Default is 8.
#' @param standardize A logical value indicating whether to standardize the data to have mean zero and variance one for each time point and each variable. Default is TRUE.
#' @param maxiteration An integer specifying the maximum number of iterations for optimization. Default is 20.
#' @param thresh A numeric value indicating the convergence threshold. Default is 1e-03.
#'
#' @return A list containing:
#'   \item{CVOut}{A matrix of cross-validation results for different tuning parameter values.}
#'   \item{mfldaerror.test}{The estimated classification error for the test data.}
#'   \item{sidaerror.train}{The estimated classification error for the training data.}
#'   \item{hatalpha}{The estimated alpha coefficients from the MFLDA.}
#'   \item{PredictedClass}{The predicted class labels for the test data.}
#'   \item{var.selected}{The variables selected by the MFLDA.}
#'   \item{varname.selected}{The names of the selected variables.}
#'   \item{optTau}{The optimal tuning parameter values.}
#'   \item{gridValues}{The grid of tuning parameter values used.}
#'   \item{myDiscPlot}{A ggplot object of the discriminant plot if \code{plotIt} is TRUE; otherwise NULL.}
#'   \item{InputData}{The original input data used.}
#'
#' @examples
#' # Example usage of cvMFLDA
#' result <- cvMFLDA(Xdata, Y, plotIt=TRUE, metrics.choice="Accuracy",
#'                    Xtestdata=test_data, Ytest=test_labels,
#'                    nfolds=5, ngrid=8)
#'
#' # Accessing results
#' print(result$mfldaerror.test)
#' print(result$myDiscPlot)
#'
#' @import foreach
#' @import doParallel
#' @import ggplot2
#' @import caret
#' @import Matrix
#' @import RSpectra
#' @import CVXR
#' @import reshape2
#' @import mgcv
#' @import nlme
#' @import lme4
#' @import philentropy
#' @import mltools
#' @import corpcor
#'
#' @export
cvMFLDA=function(Xdata=Xdata,Y=Y,withCov=FALSE,plotIt=FALSE, metrics.choice="Accuracy",
                Xtestdata=NULL,Ytest=NULL,isParallel=TRUE,ncores=NULL,
                gridMethod='RandomSearch',AssignClassMethod='Joint',
                nfolds=5,ngrid=8,standardize=TRUE,maxiteration=20,
                weight=0.5,thresh=1e-03){

  ###################### set-up #########################
  starttimeall=Sys.time()

  XdataOrig=Xdata
  XtestdataOrig=Xtestdata
  YOrig=Y
  YtestOrig=Ytest
  D = 1

  #If testing data are not provided, the default is to use training data
  if(is.null(Xtestdata)){
    Xtestdata=Xdata
    Ytest=Y
  }

  #check inputs for testing data
  ntestsizes=lapply(Xtestdata, function(x) dim(x)[1])

  if(is.null(plotIt)){
    plotIt=FALSE
  }

  if(is.null(standardize)){
    standardize=TRUE
  }


  # #standardize if true (for train)
  Xstand=list()
  nTime=length(unique(Xdata$time))
  if(standardize==TRUE){
    for(j in 1:nTime){
      myX=scale(as.matrix(Xdata[Xdata$time==j,-c(1:3)]), center=TRUE,scale=TRUE)
      Xstand[[j]]=cbind.data.frame(as.matrix(Xdata[Xdata$time==j,c(1:3)]), myX)
    }
    Xdata=do.call(rbind.data.frame,Xstand)
  }


  # #standardize if true (for test)
  Xteststand=list()
  ntestTime=length(unique(Xtestdata$time))
  if(standardize==TRUE){
    for(j in 1:ntestTime){
      myX=scale(as.matrix(Xtestdata[Xtestdata$time==j,-c(1:3)]), center=TRUE,scale=TRUE)
      Xteststand[[j]]=cbind.data.frame(as.matrix(Xtestdata[Xtestdata$time==j,c(1:3)]), myX)
    }
    Xtestdata=do.call(rbind.data.frame,Xteststand)
  }




  if(is.null(gridMethod)){
    gridMethod='RandomSearch'
  }

  if(is.null(AssignClassMethod)){
    AssignClassMethod='Joint'
  }

  if(is.null(isParallel)){
    isParallel=TRUE
  }

  if(is.null(nfolds)){
    nfolds=5
  }

  if(is.null(ngrid)){
    ngrid=8
  }


  if(is.null(maxiteration)){
    maxiteration=20
  }

  if(is.null(weight)){
    weight=0.5
  }

  if(is.null(thresh)){
    thresh=1e-03
  }

  ###################### split the folds ######################################
  set.seed(1234)
  nK=length(unique(as.vector(Y))) -1

  nc=length(unique(as.vector(Y)))
  Nn=mat.or.vec(nc,1)
  foldid=list()
  for(i in 1:nc)
  {
    Nn[i]=sum(Y==i)
    mod1=Nn[i]%%nfolds
    if(mod1==0){
      foldid[[i]]=sample(c(rep(1:nfolds,times=floor(Nn[i])/nfolds)),Nn[i])
    }else if(mod1> 0){
      foldid[[i]]=sample(c(rep(1:nfolds,times=floor(Nn[i])/nfolds), 1:(Nn[i]%%nfolds)),Nn[i])
    }
  }

  foldid=unlist(foldid)

  #print(Y)
  #print(foldid)


  #################### obtain tuning range common to all K #######################
  starttimetune=Sys.time()
  print('Getting tuning grid values')
  myTauvec=mfldatunerange(Xdata,Y,ngrid,standardize)
  #myTauvec = list(list(c(0.07798261, 0.09319861, 0.10841462, 0.12363062, 0.13884663, 0.15406263, 0.16927864, 0.18449465))) ## IBD n100_t50
  #myTauvec = list(list(c(0.06558496, 0.07834018, 0.09109539, 0.10385061, 0.11660582, 0.12936104, 0.14211626, 0.15487147)))
  #myTauvec = list(list(c(0.05336110, 0.06322495, 0.07308880, 0.08295265, 0.09281650, 0.10268036, 0.11254421, 0.12240806))) ## IBD n130_t10
  #0.07802485 0.09253828 0.10705171 0.12156514 0.13607857 0.15059200 0.16510543 0.17961886    ## IBD n90_t10
  endtimetune=Sys.time()
  print('Completed at time')
  print(endtimetune-starttimetune)

  #define the grid
  mygrid=expand.grid(do.call(cbind,myTauvec))
  gridcomb=dim(mygrid)[1]
  gridValues=mygrid

  gc()

  ################### CV ####################################
  starttimeCV=Sys.time()
  CVOut=matrix(0, nfolds, nrow(gridValues))

  #cross validation
  if(isParallel==TRUE){
    cat("Begin", nfolds,"-folds cross-validation", "\n")
    registerDoParallel()
    if(is.null(ncores)){
      ncores=parallel::detectCores()
      ncores=ceiling(ncores/2)}
    cl=makeCluster(ncores)
    registerDoParallel(cl)
    CVOut=matrix(0, nrow(gridValues), nfolds)

    ## .export=c('minv','myfastLDAnonsparse','mysqrtminv','mflda','mfldaclassify','mfldainner', 'mfldatunerange'),
    mycv=foreach(i = 1:nrow(gridValues), .combine='rbind',.export=c('minv','myfastLDAnonsparse','DiscriminantPlots','mysqrtminv','mflda','mfldaclassify','mfldainner', 'mfldatunerange'),
                 .packages=c('readr','caret', 'Matrix', 'RSpectra', 'CVXR','reshape2','mgcv','nlme','lme4','philentropy','ggplot2','mltools','corpcor')) %dopar% {
      myTau=gridValues[i,]

      #cat("Begin CV-fold", i, "\n")

      CVOut[i,]= sapply(1:nfolds, function(j){
        testInd0=which(foldid==j)
        #print(testInd0)

        testInd = subset(Xdata, time == 1)$id[testInd0]
        #print(testInd)
        #View(Xdata)

        testX=subset(Xdata, Xdata$id %in% testInd)
        testY=Y[testInd0]
        trainX=subset(Xdata, !(Xdata$id %in% testInd))
        trainY=Y[-testInd0]

        mymflda=mflda(Xtrain=trainX,Y=trainY,Tau=myTau,Xtestdata=testX,Ytest=testY,
                      plotIt=FALSE,standardize=TRUE,maxiteration=20,thresh= 1e-03)

        if (metrics.choice=="Accuracy"){return(mymflda$AverageError)}
        if (metrics.choice=="CombinedMetrics"){return(mymflda$Metrics)}

      } )
    }
    CVOut=t(mycv)
    stopCluster(cl)
  }else if(isParallel==FALSE){
    cat("Begin", nfolds,"-folds cross-validation", "\n")
    CVOut=matrix(0, nfolds, nrow(gridValues))
    for (j in 1:nfolds){
      testInd0=which(foldid==j)
      #print(testInd0)

      testInd = subset(Xdata, time == 1)$id[testInd0]
      #print(testInd)
      #View(Xdata)

      testX=subset(Xdata, Xdata$id %in% testInd)
      testY=Y[testInd0]
      trainX=subset(Xdata, !(Xdata$id %in% testInd))
      trainY=Y[-testInd0]

      #View(testX)

      #print(testY)

      cat("Begin CV-fold", j, "\n")

      CVOut[j,]= sapply(1:nrow(gridValues), function(itau){
        myTau=gridValues[itau,]
        print(itau)
        mymflda=mflda(Xtrain=trainX,Y=trainY,Tau=myTau,Xtestdata=testX,Ytest=testY,
                      plotIt=FALSE,standardize=TRUE,maxiteration=20,thresh= 1e-03)

        if (metrics.choice=="Accuracy"){return(mymflda$AverageError)}
        if (metrics.choice=="CombinedMetrics"){return(mymflda$Metrics)}
      } )
    }
  }

  View(CVOut)
  gc()

  ############################## results #######################################
  endtimeCV=Sys.time()
  print('Cross-validation completed at time')
  print(endtimeCV-starttimeCV)

  print('Getting Results......')
  #compute average classification error
  if (metrics.choice=="Accuracy"){
    minEorrInd=max(which(colMeans(CVOut, na.rm = TRUE)==min(colMeans(CVOut, na.rm = TRUE))))
    optTau=gridValues[minEorrInd,]
  } else if (metrics.choice=="CombinedMetrics"){
    maxScoreInd=max(which(colMeans(CVOut, na.rm = TRUE)==max(colMeans(CVOut, na.rm = TRUE))))
    optTau=gridValues[maxScoreInd,]
  } else {
    print("please input correct CV target")
  }



  #Apply on testing data
  moptTau= optTau
  print(moptTau)
  mymfldaTest=mflda(Xtrain=Xdata,Y=Y,Tau=moptTau,Xtestdata=Xtestdata,
                    Ytest=Ytest,plotIt=FALSE,standardize,maxiteration,thresh)

  #Apply on training data
  mysidaTrain=mflda(Xtrain=Xdata,Y=Y,Tau=moptTau,Xtestdata=Xdata,
                    Ytest=Y,plotIt=FALSE,standardize,maxiteration,thresh)

  #print out some results
  cat("Estimated Test Classification Error is", mymfldaTest$AverageError, "\n")
  cat("Estimated Train Classification Error is", mysidaTrain$AverageError, "\n")

  endtimeall=Sys.time()
  print("Total time used is")
  print(endtimeall-starttimeall)
  gc()

  ############################## plots #########################################
  #Produce discriminant and correlation plot if plotIt=T
  if(plotIt==TRUE){
    if (nk == 2){
      myalpha_plot = mymfldaTest$hatalpha
    } else if (nk == 3){
      myalpha_plot = mymfldaTest$hatalpha[[1]]
    }
    myDiscPlot <- DiscriminantPlots(Xtestdata = Xtestdata,Ytest = Ytest, myalpha = myalpha_plot,
                                    predicted.class = mymfldaTest$PredictedClass)
  }else{
    myDiscPlot=NULL
  }

  ############################## outputs #######################################
  result=list(CVOut=CVOut,sidaerror.test=mymfldaTest$AverageError,sidaerror.train=mysidaTrain$AverageError,
              hatalpha=mymfldaTest$hatalpha,PredictedClass=mymfldaTest$PredictedClass,
              optTau=moptTau,gridValues=gridValues, AssignClassMethod=AssignClassMethod,
              gridMethod=gridMethod, myDiscPlot = myDiscPlot,
              InputData=XdataOrig)
  class(result)="MFLDA"

  return(result)
}
