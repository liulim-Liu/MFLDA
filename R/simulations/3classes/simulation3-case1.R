#########################
## Simulation Function ##
## y_ij(t) = delta + eta0_j + eta1_j * t + eta2_j * t^2 + eta3_j * t^3 + eta4_j * t^4 + eta5_j * sin(eta6_j * t) + ep_ij

## for i : delta, ep_ij, eta (if the two curves are different)
## for j: eta0_j, eta1_j, eta2_j, eta3_j, eta4_j, eta5_j, eta6_j
#########################

simulation3.case1 <- function(){
  set.seed(031301)

  p <- rep(1:100, each = 1)    ## features
  t <- rep(1:40, each = 1)/10  ## time point [1, 1.01, ..., 3.99, 4]
  n1 <- 200                    ## class 1
  n2 <- 200                    ## class 2
  n3 <- 200                    ## class 3

  delta1 <- 0
  delta2 <- 500
  d <- c(10,15,20)         ## separation coefficient between two classes
  pos.neg <- c(1, -1)      ## class 1 is bigger or class 2 is bigger

  sigma <- c(100,200,300)

  ## ============ Find Eta ================== ##
  find_eta <- function(){
    x1 <- 0
    x6 <- 10


    ## random x2-5 and y1-6
    x2.5 <- runif(4, 0, 10)
    x <- c(x1, x2.5, x6)
    y <- runif(6, 50, 100)

    ## polynomial 4th degree for eta 0-4
    model <- lm(y ~ poly(x,4))
    eta0.4 <- c(model$coefficients)

    ## eta5
    f = function(t) {
      eta0.4[1] + eta0.4[2]*t + eta0.4[3]*t^2 + eta0.4[4]*t^3 + eta0.4[5]*t^4
    }

    # find the maximum of f(x) within the interval [1, 4]
    ans_min = optimize(f, interval = c(1,4), maximum = FALSE)
    ans_max = optimize(f, interval = c(1,4), maximum = TRUE)

    eta5 <- ans_max$objective - ans_min$objective

    ## eta6
    eta6 <- runif(1,0, 10)

    return(list(eta0 = eta0.4[1],
                eta1 = eta0.4[2],
                eta2 = eta0.4[3],
                eta3 = eta0.4[4],
                eta4 = eta0.4[5],
                eta5 = eta5,
                eta6 = eta6))
  }


  ############# Simulation ##################
  for (k in 1:1){
    features.all <- vector(mode='list', length=100)

    for (j in 1:100){
      result <- vector(mode='list', length=600)
      eta <- find_eta()
      delta.c <- sample(pos.neg, 1) * delta2 * sample(d,1)
      #print(delta.c)

      for (i in 1:600){     ## i from 1 to 200

        ## check whether the 10 var
        if (i > 200 & i <= 400 & j <= 10) {
          delta <- delta.c/10

          ## random the second class curve
          ## same curve as the previous or not
          if (i == 201){
            if (runif(1,0,1) >= 0.5){
              #delta.c <- delta.c/100
              eta <- find_eta()
            }
          }
        }else if (i > 400 & j <= 20 & j > 10){
          delta <- delta.c/10*2

          if (i == 401){
            if (runif(1,0,1) >= 0.5){
              eta <- find_eta()
            }
          }
        }else{
          delta <- delta1    ## class 1 OR no difference
        }

        ## epsilon
        sigma.c <- sample(sigma,1)
        ep <- rnorm(40, 1000, sigma.c)  ## larger between time point noise: rnorm(40)

        y_ij.t = delta + eta$eta0 + eta$eta1 * t + eta$eta2 * t^2 + eta$eta3 * t^3 + eta$eta4 * t^4 + eta$eta5 * sin(eta$eta6 * t) + ep
        result[[i]] <- y_ij.t
      }
      feature.j <- unlist(result)
      features.all[[j]] <- feature.j
    }


    ############# Data Frame ##################
    df <-  as.data.frame(do.call(cbind, features.all))
    time <- rep(1:40, 600)
    id <- rep(1:600, each = 40)
    group <- ifelse(id <= 200, 1, ifelse(id > 400, 3, 2))
    group <- as.factor(group)
    df <- cbind(id, time, group, df)

    #write.csv(df, paste("simulation_n300_p100.10_t40_c3_small.csv", sep = ""), row.names=FALSE)

  }

  #data <- read.csv("simulation_n300_p100.10_t40_c3_small.csv")
  #View(data)
  data <- df
  test.ind.1 <- sample(1:200, 100)
  test.ind.2 <- sample(201:400, 100)
  test.ind.3 <- sample(401:600, 100)
  test.ind <- c(test.ind.1, test.ind.2, test.ind.3)

  train.df <- subset(data, !(data$id %in% test.ind))
  test.df <- subset(data, (data$id %in% test.ind))

  id <- rep(1:300, each = 40)
  train.df$id <- id
  test.df$id <- id

  #write.csv(train.df, "train_simulation_n300_p100.10_t40_c3_small.csv", row.names=FALSE)
  #write.csv(test.df, "test_simulation_n300_p100.10_t40_c3_small.csv", row.names=FALSE)

  result = list(train.df = train.df, test.df = test.df)
  return(result)
}





