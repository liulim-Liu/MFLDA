#' @title Simulate Data for Classification Problems
#'
#' @description This function simulates datasets for classification tasks based on the specified number of classes and case number. It sources external simulation scripts to generate training and testing datasets.
#'
#' @param n.class An integer indicating the number of classes for the simulation. Acceptable values are 2 or 3.
#' @param case.num An integer specifying the case number of the simulation to execute. This should correspond to the specific simulation scenarios available for the selected number of classes.
#'
#' @return A list containing:
#'   \item{train}{A data frame of the training dataset generated from the specified simulation case.}
#'   \item{test}{A data frame of the testing dataset generated from the specified simulation case.}
#'
#' @details
#' The function will source R scripts from predetermined paths based on the number of classes (2 or 3).
#' For each class scenario, the function can simulate up to four different cases, each generating distinct datasets.
#' It is important that the simulation scripts exist in the specified directories, and they should return a list containing `train.df` and `test.df`.
#' There are four different cases for simulation:
#'  (1) the groups are distinct across all time points, exhibiting similar average trends;
#'  (2) the groups diverge only within a specific time window, from the time point [5, 15], displaying different patterns;
#'  (3) the groups show separation over a fixed interval of 10 time points, though the exact time frame varies;
#'  (4) the groups are distinct during a random occurring period, with the length of the separation varying between 5 and 40 time points.
#'
#' @examples
#' # Simulate data for a 2-class problem, case 1
#' result <- simulation(n.class = 2, case.num = 1)
#'
#' # Access training and testing datasets
#' train_data <- result$train
#' test_data <- result$test
#'
#' @export
simulation <- function(n.class, case.num){
  if (n.class == 2){
    source("R/simulations/2classes/simulation2-case1.R")
    source("R/simulations/2classes/simulation2-case2.R")
    source("R/simulations/2classes/simulation2-case3.R")
    source("R/simulations/2classes/simulation2-case4.R")

    if (case.num == 1){
      train <- simulation2.case1()$train.df
      test <- simulation2.case1()$test.df
    }
    if (case.num == 2){
      train <- simulation2.case2()$train.df
      test <- simulation2.case2()$test.df
    }
    if (case.num == 3){
      train <- simulation2.case3()$train.df
      test <- simulation2.case3()$test.df
    }
    if (case.num == 4){
      train <- simulation2.case4()$train.df
      test <- simulation2.case4()$test.df
    }
  }

  if (n.class == 3){
    source("R/simulations/3classes/simulation3-case1.R")
    source("R/simulations/3classes/simulation3-case2.R")
    source("R/simulations/3classes/simulation3-case3.R")
    source("R/simulations/3classes/simulation3-case4.R")
    source("R/simulations/3classes/simulation3-case5.R")

    if (case.num == 1){
      train <- simulation3.case1()$train.df
      test <- simulation3.case1()$test.df
    }
    if (case.num == 2){
      train <- simulation3.case2()$train.df
      test <- simulation3.case2()$test.df
    }
    if (case.num == 3){
      train <- simulation3.case3()$train.df
      test <- simulation3.case3()$test.df
    }
    if (case.num == 4){
      train <- simulation3.case4()$train.df
      test <- simulation3.case4()$test.df
    }
    if (case.num == 5){
      train <- simulation3.case5()$train.df
      test <- simulation3.case5()$test.df
    }
  }
  result = list(train = as.data.frame(train), test = as.data.frame(test))
  return(result)
}

