rm(list = ls())
gc()

home.path <- getwd()
setwd(home.path)

################### Import Librarys ###################
library(readr)
library(Matrix)
library(RSpectra)
library(CVXR)
library(caret)
library(doParallel)
library(ggplot2)
library(reshape2)
library(philentropy)
library(mgcv)
library(nlme)
library(lme4)
library("dplyr")
library(mltools)
library("corpcor")

################### Import Functions ###################
source("R/helperfunctions.R")
source("R/mflda.R")
source("R/mfldatunerange.R")
source("R/mfldaclassify.R")
source("R/cvmflda.R")
source("R/discriminationPlot.R")
source("R/simulation.R")
source("R/spline.R")
source("R/bs.generator.R")

################### Data Simulation ###################
## Binary Class Simulation
data <- simulation(n.class = 2, case.num = 1) ## choose case.num from 1-4
View(data$train)
View(data$test)

## Multi Class Simulation
data <- simulation(n.class = 3, case.num = 5) ## choose case.num from 1-4
View(data$train)
View(data$test)

## b-spline estimation
train.pred <- spline.prediction(data$train)
View(train.pred$pred.df)
test.pred <- spline.prediction(data$test)
View(test.pred$pred.df)

################### MFLDA-3-class ###################
## Read data
XTrain <- get(load("data/simTrain_spline_c3_case4.rda"))
XTest <- get(load("data/simTest_spline_c3_case4.rda"))

trainX <- XTrain
trainX$group <- as.integer(trainX$group)
trainY <- trainX$group[trainX$time == 1]
trainY <- as.integer(trainY)

testX <- XTest
testY <- testX$group[testX$time == 1]
testX$group <- as.integer(testX$group)
testY <- as.integer(testY)

## Tuning parameter
cvmod <- cvMFLDA(Xdata=trainX,Y=trainY, metrics.choice="CombinedMetrics") ## Accuracy or CombinedMetrics
cvmod$optTau
gc()

## mflda with known hyper-parameter
myTau=cvmod$optTau 
#myTau=50 ## can choose your own
mflda.result <- mflda(Xtrain=trainX,Y=trainY,Tau=myTau,Xtestdata=testX,Ytest=testY,plotIt=FALSE,standardize=TRUE,maxiteration=20,thresh= 1e-03)

## find the predicted classes
mflda.result$PredictedClass

## find the first discriminant scores
mflda.result$hatalpha[[1]]

## find the second discriminant scores
mflda.result$hatalpha[[2]]

## find the selected variables
mflda.result$varname.selected

## plot the discriminant plots for the first discriminant vector
myplot <- DiscriminantPlots(Xtestdata=testX,Ytest=testY,
myalpha=mflda.result$hatalpha[[1]], predicted.class = mflda.result$PredictedClass)
myplot$discriminant.plot
myplot$density.loess
myplot$density.plot

################### MFLDA-2-class ###################
## Read data
XTrain <- get(load("data/simTrain_spline_c2_case4.rda"))
XTest <- get(load("data/simTest_spline_c2_case4.rda"))

## Rest is the similar to the 3-class problem
