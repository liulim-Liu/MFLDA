#' @title Spline-Based Estimation from Dataset
#'
#' @description This function takes a dataset of participant observations and uses spline basis functions to generate fitted values for each variable across unique time points. It constructs a linear model for each variable and predicts values using B-splines.
#'
#' @param olddat A data frame containing the dataset with participant observations. It must include columns for `id`, `group`, `time`, and additional variables to be predicted.
#'
#' @return A list containing:
#'   \item{pred.df}{A data frame with the predicted values for each variable across time points, along with participant IDs and group classifications.}
#'
#' @details
#' The function internally identifies unique time points from the dataset and generates B-spline basis matrices for these time points.
#' It uses linear models to predict values for each variable, storing the fitted values in a list.
#' The function also assigns group classifications based on participant IDs if the original data has irregular time points.
#'
#' @examples
#' # Assume `my_data` is a data frame with the necessary structure
#' prediction_results <- spline.prediction(my_data)
#'
#' # View the predicted values
#' head(prediction_results$pred.df)
#'
#' @importFrom splines bs
#' @export
spline.prediction <- function(olddat){
  dat <- olddat
  print("successfully import dataset")

  dat$group <- as.factor(dat$group)
  x = unique(dat$time)

  # parameters
  nk = 10 # number of knots 5~15
  nv = ncol(dat) - 3 # number of variables
  print(nv)
  id = unique(dat$id)
  nid = length(id)
  print(nid)
  n = nid
  n.class <- n/2

  B0 <- bs.generator(x, nk, q = 3)$Bx0
  B <- bs.generator(x, nk, q = 3)$Bx

  yhat.all <- vector(mode = 'list', length = nv)

  # main iteration
  print("Start iteration")
  for (j in 1:nv){ # for each variable
    yhat.j <- vector(mode = 'list', length = nid)
    dat.j = dat[, (j + 3)]

    for (i in 1:nid){ # for each participant
      dat.ij <- dat.j[dat$id == i] # get the value of each variable for each participants
      mfit.ij = lm(dat.ij ~ B)
      yhat.ij <- mfit.ij$fitted.values # fitted curve

      # save all the results
      yhat.j[[i]] <- yhat.ij
    }
    yhat.j <- unlist(yhat.j)
    yhat.all[[j]] <- yhat.j
  }

  # Prediction
  print("Start Prediction")
  pred.df <-  as.data.frame(do.call(cbind, yhat.all))
  time <- rep(1:40, n)
  id <- rep(1:n, each = 40)
  group <- ifelse(id <= 100, 1, ifelse(id > 200, 3, 2))
  group <- as.factor(group)
  pred.df <- cbind(id, time, group, pred.df)

  return(list(pred.df = pred.df))
}
