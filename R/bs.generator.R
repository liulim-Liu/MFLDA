######################################################################
# Input Arguments:
#	(1) x: Covariates. represent the time in the project
#	(2) N: Number of interior knots in generating spline matrix.
#       suggested to be 3,4,5 in our case of 40 time points
#	(3) q: Degree of polynomial spline. Default is 3. Cubic spline
#	(4) KnotsLocation: A character string naming the way for knots
#			locations. Default is "quantile". The only alternative is
#			"uniform".
#           "uniform" means the knots are uniformly distributed, such as
#           each month in a longitudinal study
#	(5) knots: An optional vector specifying the knots to be used in
#			constructing spline bases.
#
# Output Arguments:
#	(1) Bx0: Non-centered spline bases for covariates. (original phis) - use this as the phi results
#	(2) Bx: Centered spline bases for covariates. (phis centeralized)
#	(3) BxMean: Means of splinen bases substracted.
#	(4) knots: Knots used in constructing spline bases.


#' @title Generate B-spline Basis Matrices
#'
#' @description This function creates B-spline basis matrices for a given input vector or matrix, with options to customize the degree of the spline, the location of knots, and whether the knots are generated based on quantiles or uniform intervals.
#'
#' @param x A numeric vector or matrix representing the input data for which B-spline bases are to be generated. If a matrix, each column is treated as a separate variable.
#' @param N An integer specifying the number of internal knots to generate.
#' @param q An integer representing the degree of the spline. Default is 3 (cubic spline).
#' @param KnotsLocation A character string indicating the method for generating knots. Options are `"quantile"` for quantile-based knots or `"uniform"` for uniformly spaced knots. Default is `"quantile"`.
#' @param knots An optional vector of custom knots. If provided, this vector will be used directly, and `KnotsLocation` will be ignored.
#'
#' @return A list containing:
#'   \item{Bx0}{The raw B-spline basis matrix.}
#'   \item{Bx}{The B-spline basis matrix with centered columns.}
#'   \item{BxMean}{The column means of `Bx0`, used for centering.}
#'   \item{knots}{The knots used for generating the B-spline basis.}
#'
#' @details
#' Knots Generation: If `knots` is not provided, knots are generated based on the value of `KnotsLocation`. For quantile-based knots, they are placed at quantiles of the data in `x`, while for uniform spacing, knots are spread evenly across the range of `x`.
#' Matrix Handling: If `x` is a matrix, each column is treated as an independent variable, and the spline basis is constructed separately for each. This reduces memory usage by processing in blocks.
#'
#' @examples
#' # Generate a cubic B-spline basis for a vector with quantile-based knots
#' x <- rnorm(100)
#' bs.generator(x, N = 5, q = 3, KnotsLocation = "quantile")
#'
#' # Generate a cubic B-spline basis for a matrix with uniform knots
#' x_matrix <- matrix(rnorm(1000), ncol = 10)
#' bs.generator(x_matrix, N = 5, q = 3, KnotsLocation = "uniform")
#'
#' @importFrom splines bs
#' @export
######################################################################
bs.generator <- function(x, N, q = 3, KnotsLocation = "quantile", knots = NULL){
  library(splines)
	# Case I: x is a vector
	if (length(dim(x)) == 0) {
		# generate knots by supplied location
		if (is.null(knots)) {
			# knots generated in sample quantile
			if (KnotsLocation == "quantile") {
				knots1 = quantile(x, probs = seq(0, 1, by = 1/(N+1)))
			}
			# knots generated in uniform points
			if (KnotsLocation == "uniform") {
				knots1 = seq(min(x), max(x), by = (max(x)-min(x))/(N+1))
			}
		} else {
		  knots1 = knots
		}
		# generate constant / polynomial B-spline basis
		if (q == 0) {
			Bx0 = cbs(x, knots1[-c(1, N+2)], Boundary.knots = knots1[c(1, N+2)])
		} else {
		  Bx0 = bs(x, knots = knots1[-c(1, N+2)], degree = q, intercept = F,
					Boundary.knots = knots1[c(1, N+2)])
		}
		Knots = knots1
	} else { # Case II: x is a matrix
		d.x = ncol(x)
		block.size = floor(sqrt(d.x))
		nblock = ceiling(d.x/block.size)
		Bx0 = NULL
		Knots = NULL

		# This routine reduces the memory needed in constructing spline
		#	bases and improve the computational efficiency
		for(nj in 1:nblock) {
			Bxnj = NULL
			Knotsj = NULL
			if(nj < nblock) {
			  block.ind = (block.size*(nj-1)+1):(block.size*nj)
			}
			if(nj==nblock) {
			  block.ind = (block.size*(nj-1)+1):d.x
			}
			for (j in block.ind) {
				# generate knots by supplied location
				if (is.null(knots)) {
					# knots generated in sample quantile
					if (KnotsLocation == "quantile") {
					  knots1 = quantile(x[, j], probs = seq(0, 1, by = 1/(N+1)))
					}
					# knots generated in uniform points
					if (KnotsLocation == "uniform") {
					  knots1 = seq(min(x[, j]), max(x[, j]), by = (max(x[,j])-min(x[,j]))/(N+1))
					}
				} else {
				  knots1 = knots[(N+2)*(j-1)+1:(N+2)]
				}

				# generate constant / polynomial B-spline basis
				if (q == 0) {
					bx0 = cbs(x[, j], knots1[-c(1, N+2)], Boundary.knots = knots1[c(1, N+2)])
				} else {
				  bx0 = bs(x[, j], knots = knots1[-c(1, N+2)], degree = q,
							intercept = F, Boundary.knots = knots1[c(1, N+2)])
				}
				Bxnj = cbind(Bxnj, bx0)
				Knotsj = cbind(Knotsj, knots1)
			}
			Bx0 = cbind(Bx0, Bxnj)
			Knots = cbind(Knots, Knotsj)
		}
	}
	BxMean = colMeans(Bx0)
	Bx = sweep(Bx0, 2, BxMean, "-")
	list(Bx0 = Bx0, Bx = Bx, BxMean = BxMean, knots = Knots)
}

######################################################################
# This routine generates constant spline bases
# want no overlap within sub-intervals

#' @title Generate Constant B-spline Basis
#'
#' @description This function creates a constant spline basis matrix for a given vector of input values. It assigns binary indicators for which interval each value falls into based on specified knots.
#'
#' @param x A numeric vector for which the constant spline basis is to be generated.
#' @param knots A numeric vector of knots that define the intervals for the constant spline basis.
#' @param Boundary.knots A numeric vector of length two defining the boundary knots. Default is the range of `x`. These knots represent the outer boundaries for the spline basis.
#'
#' @return A matrix of binary indicators with rows corresponding to the input vector `x` and columns corresponding to the defined knots. Each entry is `1` if the corresponding value in `x` falls within the interval defined by the knots and `0` otherwise.
#'
#' @details
#' The function ensures that there is no overlap between the defined intervals by assigning values only to the corresponding interval where each input falls.
#' The first column of the output matrix corresponds to values less than the first knot, and the last column corresponds to values greater than or equal to the last knot.
#' @export
cbs <- function(x, knots, Boundary.knots = range(x)){
  n = length(x)
  N = length(knots)
  Ix = matrix(0L, nrow = n,ncol = N+1L)
  for (i in 1:n) {
    if (x[i] < knots[1L]) {
      Ix[i, 1] = 1L
    } else if (x[i] >= knots[N]) {
      Ix[i, N+1L] = 1L
    } else {
      kl = max(which(x[i] >= knots))
      Ix[i, kl+1] = 1L
    }
  }
  return(Ix)
}
