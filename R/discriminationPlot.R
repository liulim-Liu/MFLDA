#' @title Generate Discriminant Plots for Test Data
#'
#' @description This function generates various discriminant plots for test data based on projections from a linear discriminant analysis.
#' It includes a discrimination plot, density plots of the projections, and performs Kolmogorov-Smirnov tests.
#'
#' @param Xtestdata A data frame containing the test data in long format, with the following columns: `ID` (first column), `Time` (second column), `Group` (third column), Variables (fourth to last columns)
#' @param Ytest A vector representing the group labels for the test data.
#' @param myalpha A matrix of coefficients obtained from the linear discriminant analysis, where each column corresponds to a time point.
#' @param predicted.class A vector of predicted class labels for each observation in the test data.
#'
#' @return A list containing:
#'   \item{Projection.df}{The data frame containing the projections for the test data.}
#'   \item{discriminant.plot}{A ggplot object representing the discrimination plot.}
#'   \item{density.loess}{A ggplot object representing the density plot of the loess fit for projections.}
#'   \item{density.plot}{A ggplot object representing the density plot of row means.}
#'   \item{ks.loess.p}{The p-value from the Kolmogorov-Smirnov test on the loess predictions (if applicable).}
#'   \item{ks.density.p}{The p-value from the Kolmogorov-Smirnov test on the density of row means (if applicable).}
#'
#' @examples
#' # Example usage of DiscriminantPlots
#' result <- DiscriminantPlots(test_data, test_labels, coefficients_matrix, predicted_classes)
#'
#' # Plot the discriminant plot
#' print(result$discriminant.plot)
#'
#' @import ggplot2
#' @import reshape2
#'
#' @export
DiscriminantPlots <- function(Xtestdata,Ytest,myalpha, predicted.class){
  # #standardize
  Xteststand=list()
  ntestTime=length(unique(Xtestdata$time))
  nk=length(unique(Xtestdata$group))

  for(j in 1:ntestTime){
    myX=scale(as.matrix(Xtestdata[Xtestdata$time==j,-c(1:3)]), center=TRUE,scale=TRUE)
    Xteststand[[j]]=cbind.data.frame(as.matrix(Xtestdata[Xtestdata$time==j,c(1:3)]), myX)
  }
  Xtestdata=do.call(rbind.data.frame,Xteststand)
  
  Projtest=list()
  for(j in 1:ntestTime){
    Projtest[[j]]=as.matrix(Xtestdata[Xtestdata$time==j,-c(1:3)])%*%myalpha[,j]
  }
  
  Projtest.df.keep <- as.data.frame(do.call(cbind, Projtest)) #n x t
  

  ########## reshape the projection matrix ####################################
  Projtest.df <- Projtest.df.keep
  Projtest.df <- cbind(newColName = rownames(Projtest.df), Projtest.df)
  rownames(Projtest.df) <- 1:nrow(Projtest.df)

  Projtest.df.melted = melt(Projtest.df, id.vars = 'newColName')
  colnames(Projtest.df.melted) <- c("sample", "var", "output")

  time.ind <- rep(1:ntestTime, each = length(unique(Projtest.df.melted$sample)))
  Projtest.df.melted$time <- time.ind

  group.ind <- rep(Ytest, times = ntestTime)
  Projtest.df.melted$classes <- group.ind
  Projtest.df.melted$classes <- as.factor(Projtest.df.melted$classes)

  #View(Projtest.df.melted)

  ################# discrimination plot #######################################
  colormap <- c('#00bbd6','#faa32b','#a832a4')
  df2 <- aggregate(x=Projtest.df.melted$output, list(Projtest.df.melted$classes, Projtest.df.melted$time), FUN=mean)
  names(df2) <- c("classes","time", "ymean")
  
  p <- ggplot(data = Projtest.df.melted, aes(x = time, y = output, group = sample, color = classes)) +
    geom_line(linewidth = 0.1, alpha = 0.1) +
    geom_smooth(data=df2, method = "loess", se = FALSE, aes(x=time, y=ymean, group = classes), linewidth = 1.5)+
    scale_color_manual(values=colormap[1:nk])+
    theme_classic()+
    xlab("Time") + ylab("Projection") +
    theme(axis.text=element_text(size=15),axis.title=element_text(size=20,face="bold"),
          legend.text=element_text(size=15), legend.title=element_text(size=15,face="bold"))
  
  ################# loess curve density plot ###################################
  Projtest.df.melted.1 <- subset(Projtest.df.melted, classes == 1)
  Projtest.df.melted.2 <- subset(Projtest.df.melted, classes == 2)

  loess_fit_1 <- loess(output ~ time, Projtest.df.melted.1)
  loess_fit_2 <- loess(output ~ time, Projtest.df.melted.2)
  if (nk == 2){
    Projtest.df.melted.dens <- data.frame("loess.pred" = c(predict(loess_fit_1),predict(loess_fit_2)),
                                          "time" = c(1:length(predict(loess_fit_1)),1:length(predict(loess_fit_2))),
                                          "predClasses" = c(rep(1,length(predict(loess_fit_1))),rep(2,length(predict(loess_fit_2))))
    )
  }
  
  if (nk == 3){
    Projtest.df.melted.3 <- subset(Projtest.df.melted, classes == 3)
    loess_fit_3 <- loess(output ~ time, Projtest.df.melted.3)
    Projtest.df.melted.dens <- data.frame("loess.pred" = c(predict(loess_fit_1),predict(loess_fit_2),predict(loess_fit_3)),
                                          "time" = c(1:length(predict(loess_fit_1)),1:length(predict(loess_fit_2)),1:length(predict(loess_fit_3))),
                                          "predClasses" = c(rep(1,length(predict(loess_fit_1))),rep(2,length(predict(loess_fit_2))),rep(3,length(predict(loess_fit_3))))
    )
  }
  
  Projtest.df.melted.dens$loess.pred <- as.numeric(Projtest.df.melted.dens$loess.pred)
  Projtest.df.melted.dens$predClasses <- as.factor(Projtest.df.melted.dens$predClasses)
  
  p3 <- ggplot(Projtest.df.melted.dens, aes(loess.pred, fill = predClasses, colour = predClasses, after_stat(density))) +
    geom_density(trim=TRUE,linewidth = 1, alpha = 0.1) +
    # Apply your custom colormap for fill
    scale_fill_manual(values = colormap) +
    # Apply your custom colormap for color (for the lines)
    scale_color_manual(values = colormap)
  
  if (nk == 2){
    ks.loess <- ks.test(loess.pred ~ predClasses, data = Projtest.df.melted.dens)$p.value
  } else {
    ks.loess <- NULL
  }
  ################# density plot ##############################################
  Projtest.df <- Projtest.df.keep
  Projtest.rowMeans <- rowMeans(Projtest.df)
  
  Projtest.df.rowMeans <- data.frame("rowM" <- Projtest.rowMeans, "predClass" <- predicted.class)
  colnames(Projtest.df.rowMeans) <- c("rowM", "predClass")
  Projtest.df.rowMeans$rowM <- as.numeric(Projtest.df.rowMeans$rowM)
  Projtest.df.rowMeans$predClass <- as.factor(Projtest.df.rowMeans$predClass)
  
  p2 <- ggplot(Projtest.df.rowMeans, aes(rowM, fill = predClass, colour = predClass, after_stat(density))) +
    geom_density(trim=TRUE,linewidth = 1, alpha = 0.1)+
    # Apply your custom colormap for fill
    scale_fill_manual(values = colormap) +
    # Apply your custom colormap for color (for the lines)
    scale_color_manual(values = colormap)+
    theme_classic()+
    xlab("Projection") + ylab("Density") +
    theme(axis.text=element_text(size=15),axis.title=element_text(size=20,face="bold"),
          legend.text=element_text(size=15), legend.title=element_text(size=15,face="bold"))
  
  if (nk == 2){
    ## Kolmogorov-Smirnov Tests
    ks.p <- ks.test(rowM ~ predClass, data = Projtest.df.rowMeans)$p.value
  } else {
    ks.p <- NULL
  }
  
  result = list(Projection.df = Projtest.df, discriminant.plot  = p, density.loess = p3, density.plot = p2,
                ks.loess.p = ks.loess, ks.density.p = ks.p)
  return(result)
}