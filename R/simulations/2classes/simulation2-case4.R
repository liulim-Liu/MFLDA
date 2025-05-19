#########################
## Simulation Function ##
## y_ij(t) = delta + eta0_j + eta1_j * t + eta2_j * t^2 + eta3_j * t^3 + eta4_j * t^4 + eta5_j * sin(eta6_j * t) + ep_ij

## for i : delta, ep_ij, eta (if the two curves are different)
## for j: eta0_j, eta1_j, eta2_j, eta3_j, eta4_j, eta5_j, eta6_j
#########################

simulation2.case4 <- function(){
  set.seed(102373)

  p <- rep(1:100, each = 1)    ## features
  t <- rep(1:40, each = 1)/10  ## time point [1, 1.01, ..., 3.99, 4]
  n1 <- 200                    ## class 1
  n2 <- 200                    ## class 2

  delta1 <- 0
  delta2 <- 100
  d <- c(10,15,20)         ## separation coefficient between two classes
  pos.neg <- c(1, -1)     ## class 1 is bigger or class 2 is bigger

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
    features.all0 <- vector(mode='list', length=100)
    features.all1 <- vector(mode='list', length=100)

    for (j in 1:100){
      result0 <- vector(mode='list', length=200)
      result1 <- vector(mode='list', length=200)

      eta <- find_eta()
      eta1 <- eta

      delta.c <- sample(pos.neg, 1) * delta2 * sample(d,1)
      for (i in 1:400){     ## i from 1 to 200

        ## check whether the 10 var
        if (i > 200 & j <= 10) {
          delta <- delta.c

          ## random the second class curve
          ## same curve as the previous or not
          if (i == 201){
            eta1 <- find_eta()
            # if (runif(1,0,1) >= 0.5){
            #   eta <- find_eta()
            # }
          }

        }else{
          delta <- delta1    ## class 1 OR no difference
        }

        ## epsilon
        sigma.c <- sample(sigma,1)
        ep <- rnorm(40, 1000, sigma.c)  ## larger between time point noise: rnorm(40)

        y_ij.t.1 = delta + eta1$eta0 + eta1$eta1 * t + eta1$eta2 * t^2 + eta1$eta3 * t^3 + eta1$eta4 * t^4 + eta1$eta5 * sin(eta1$eta6 * t) + ep
        y_ij.t.0 = delta1 + eta$eta0 + eta$eta1 * t + eta$eta2 * t^2 + eta$eta3 * t^3 + eta$eta4 * t^4 + eta$eta5 * sin(eta$eta6 * t) + ep

         # changes in time period
        # y_ij.t.mean <- mean(y_ij.t)
        # if (i > 100 & j <= 10){
        #   for (t in 1:40){
        #     if (t %in% c(5:15)){
        #       y_ij.t[t] <- y_ij.t[t] * 0.5
        #     }
        #   }
        # }

        result1[[i]] <- y_ij.t.1
        result0[[i]] <- y_ij.t.0
      }
      feature.j.1 <- unlist(result1)
      features.all1[[j]] <- feature.j.1

      feature.j.0 <- unlist(result0)
      features.all0[[j]] <- feature.j.0
    }


    ############# Data Frame ##################
    df1 <-  as.data.frame(do.call(cbind, features.all1))
    time <- rep(1:40, 400)
    id <- rep(1:400, each = 40)
    group <- c(id > 200) + 1
    group <- as.factor(group)
    df1 <- cbind(id, time, group, df1)

    #write.csv(df1, paste("simulation1.csv", sep = ""), row.names=FALSE)

    df0 <-  as.data.frame(do.call(cbind, features.all0))
    time <- rep(1:40, 400)
    id <- rep(1:400, each = 40)
    group <- c(id > 200) + 1
    group <- as.factor(group)
    df0 <- cbind(id, time, group, df0)

    #write.csv(df0, paste("simulation0.csv", sep = ""), row.names=FALSE)

    ############# time period differences ##################

    #df0 <- read.csv("simulation0.csv")
    #df1 <- read.csv("simulation1.csv")

    ## replace time start:end from simulation 0 to simulation 1 by each variable
    for (j in 4:103){
      start <- sample(c(1:30), 1)
      start2 <- start+2
      end <- sample(c(start2:40),1)
      #print(start)
      #print(end)
      df0[df0$time %in% c(start:end),j] <- df1[df1$time %in% c(start:end),j]
    }

    #write.csv(df0, paste("simulation_n200_p100.10_t40_c2_tp7.3.csv", sep = ""), row.names=FALSE)

  }


  #data <- read.csv("simulation_n200_p100.10_t40_c2_tp7.3.csv")
  #View(data)
  data <- df0
  test.ind.1 <- sample(1:200, 100)
  test.ind.2 <- sample(201:400, 100)
  test.ind <- c(test.ind.1, test.ind.2)

  train.df <- subset(data, !(data$id %in% test.ind))
  test.df <- subset(data, (data$id %in% test.ind))

  id <- rep(1:200, each = 40)
  train.df$id <- id
  test.df$id <- id

  result = list(train.df = train.df, test.df = test.df)
  return(result)
}
