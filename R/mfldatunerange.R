#' @title Generate Tuning Range for Multivariate Discriminant Analysis
#'
#' @description This function computes the tuning ranges for the sparsity parameter in multivariate functional linear discriminant analysis (MFLDA) based on provided training data. It standardizes the data if specified and uses a grid search approach to determine optimal parameter values.
#'
#' @param Xtrain A data frame containing the training data. The first three columns should contain identifiers (such as `id`, `time`, and `group`), and subsequent columns should contain the features used for classification.
#' @param Y A vector of class labels corresponding to the training data. It should have the same length as the number of rows in `Xtrain`.
#' @param ngrid An integer specifying the number of grid points to generate for the tuning range. Default is 8.
#' @param standardize A logical value indicating whether to standardize the data to have mean zero and variance one for each time point and each variable. Default is TRUE.
#'
#' @return A list containing:
#'   \item{Tauvec}{A list of sequences of tuning parameters for each class dimension based on the grid search results.}
#'
#' @details
#' The function performs the following steps:
#' 1. It prepares the input data, optionally standardizing it based on the `standardize` parameter.
#' 2. It uses a custom function `myfastLDAnonsparse` to perform the preliminary analysis.
#' 3. It generates a grid of tuning parameters by iteratively refining based on selectivity and sparsity conditions.
#' 4. The final output contains a range of optimal tuning parameters for each class.
#'
#' @examples
#' # Assuming Xtrain is your training data frame and Y contains the corresponding labels
#' result <- mfldatunerange(Xtrain = Xtrain, Y = Y, ngrid = 10, standardize = TRUE)
#'
#' # Access the tuning parameters
#' tuning_params <- result$Tauvec
#'
#' @export
mfldatunerange=function(Xtrain=Xtrain,Y=Y,ngrid=8,standardize=TRUE){

  #check size of each data
  # if (is.list(Xtrain)) {
  #   Xdata = Xtrain
  # } else {
  #   Xdata=list(Xtrain)
  # }
  # dsizes=lapply(Xdata, function(x) dim(x))
  #n=dsizes[[1]][1]
  Xdata1=as.matrix(Xtrain[Xtrain$time==1,-c(1:3)])
  n=length(Y)
  Xdata=list(Xdata1)
  dsizes=lapply(Xdata, function(x) dim(x))
  p=lapply(Xdata, function(x) dim(x)[2])
  #print(length(p))
  #D=length(dsizes)
  D=1

  if(is.null(ngrid)){
    ngrid=8
  }

  if(is.null(standardize)){
    standardize=TRUE
  }


  nK=length(unique(as.vector(Y))) -1

  #standardize if true
  Xstand=list()
  if(standardize==TRUE){
    nTime=length(unique(Xtrain$time))
    for(j in 1:nTime){
      myX=scale(as.matrix(Xtrain[Xtrain$time==j,-c(1:3)]), center=TRUE,scale=TRUE)
      Xstand[[j]]=cbind.data.frame(as.matrix(Xtrain[Xtrain$time==j,c(1:3)]), myX)
    }
  }

  Xstand2=do.call(rbind.data.frame,Xstand)

  #obtain nonsparse solutions
  myfastlda=myfastLDAnonsparse(Xstand2, Y)
  SqrtmSwSbSqrtmSw=myfastlda$SqrtmSwSbSqrtmSw #PT x PT
  myalphaold=myfastlda$tildealphamat #PT x 1
  Ux=myfastlda$Ux #PT x n
  tildelambda=myfastlda$tildelambda #1 x 1
  
  #print(tildelambda)
  #View(SqrtmSwSbSqrtmSw)

  #obtain upper and lower bounds
  Separation=list()
  nTime=length(unique(Xtrain$time))
  Separation=SqrtmSwSbSqrtmSw%*%myalphaold #PT x 1
  
  ## make separation from 4000 x 1 back to matrix 100 x 40
  Separation_m <- matrix(Separation, ncol = p[[1]][1], byrow = TRUE)
  Separation_m <- t(Separation_m)
  #View(Separation_m)

  ## generate the grids
  #grid.tau <- apply(Separation_m,2,max)* (1/sqrt(p[[1]][1]))
  dom <- sqrt(p[[1]][1])
  grid.tau <- apply(Separation_m,2,max)* (1/dom)
  
  grid.10 <- max(grid.tau)
  grid.2 <- min(grid.tau)
  grid.1 <- lapply(1:D, function(x) sqrt(log(p[[1]])/(n*nTime))*grid.2[[x]])
  grid.3.9 <- rep(NA, 7)
  for (i in 1:7){
    grid.3.9[i] <- grid.2 + i * ((grid.10-grid.2)/8)
  }
  grid.tau <- unlist(as.vector(c(grid.1, grid.2, grid.3.9, grid.10)))
  print(grid.tau)

  ## set up the cutoff
  sparsity.low <- 0.1
  sparsity.up <- 0.25
  selectivity <- 0.7

  k = 1  

  while(TRUE){
    print(grid.tau)
    ## calculate the selectivity percantage for each tau
    selectivity.all <- matrix(0, p[[1]][1], length(grid.tau))

    myalphaold=as.matrix(myalphaold) ## change to matrix (for inner function)
    for (i in 1:length(grid.tau)){
      tau <- grid.tau[i]
      #alpha <- quick.mflda(Xtrain=Xtrain,Y=Y,Tau=tau)$hatalpha
      alpha <- mfldainner(Xstand2,Y,Ux,SqrtmSwSbSqrtmSw,myalphaold,tildelambda,tau)$hatalpha[[1]][[1]]
      alpha[alpha == 0] <- NA
      #View(alpha)
      alpha.non.zero <- rowSums(!is.na(alpha))
      selectivity.all[,i] <- alpha.non.zero/nTime
    }

    print(selectivity.all)

    ## find the index of taus which has sparsity needed
    selectivity.all[selectivity.all < selectivity] <- NA
    print(selectivity.all)

    sparsity <- colSums(!is.na(selectivity.all))/p[[1]][1]
    print(sparsity)

    grid.tau.new <- unlist(grid.tau[between(sparsity,sparsity.low,sparsity.up)])
    print(length(grid.tau.new))

    if (length(grid.tau.new) != 0){
      break
    }

    if (k > 10){
      k = 2
    } else {
      k = k + 0.5
    }
    
    #grid.tau <- grid.tau / k #
    print(sparsity > sparsity.up)
    grid.tau <- ifelse(sparsity > sparsity.up, grid.tau * k, grid.tau / k)
    gc()
  }

  ubx <- max(grid.tau.new)
  lbx <- min(grid.tau.new)

  #print(c(ubx, lbx))
  ######################## end way 4 ###########################################

  # lbx=lapply(1:D, function(x) sqrt(log(p[[1]])/(n*nTime))*ubx[[x]])
  # #lbx = 0
  # #ubx=lapply(1:D, function(x) lbx[[x]]*1.2)
  #
  #tuning range for each data
  Taugrid=list()
  cc=lapply(1, function(x1,x2)  cbind(lbx,ubx))

  print(cc)

  cc=as.matrix(do.call(rbind,cc))
  for(d in 1:D){
    Taugrid[[d]]=seq(as.numeric(cc[d,1])/(3*nTime),as.numeric(cc[d,2]),length.out=(ngrid+1))
  }

  #print(Taugrid)

  # up to 25% sparsity

  myperx=lapply(Taugrid, function(x) quantile(x[1:ngrid], c(.1, .15, .2, .25, .35, .45), type=5))#similar to matlab
  myperx2=do.call(rbind,myperx)

  #print(myperx)
  #print(myperx2)


  for(loc in 1:6){
    mTaux=sapply(1:D, function(i) list(t(myperx2[i,loc])))
    myres=mfldainner(Xstand2,Y,Ux,SqrtmSwSbSqrtmSw,myalphaold,tildelambda,mTaux)

    nnz=sapply(1:D, function(i) list(colSums(myres$hatalpha[[i]][[1]]!=0)/dsizes[[i]][2]))
    nnz=cbind(c(do.call(rbind,nnz)))
    #print(nnz)
    if(all(nnz<=0.25)){
      break
    }
  }
  #final grid
  #print(loc)
  Tauvec=sapply(1:D, function(i) list(seq(as.numeric(t(myperx2[i,loc])),as.numeric(ubx[[i]]),len=(ngrid+1))))
  #Tauvec=sapply(1:D, function(i) list(seq(as.numeric(t(myperx2[i,loc])),as.numeric(ubx[[i]]),len=(ngrid+1))))

  Tauvec=sapply(1:D, function(x) list(Tauvec[[x]][1:ngrid]))
  print(Tauvec)

  #result=list(Tauvec=Taugrid)
  result=list(Tauvec=Tauvec)
  return(result)
}


