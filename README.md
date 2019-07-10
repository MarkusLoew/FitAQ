
# FitAQ
Yet another R package to analyse the relationship between photosynthetic assimilation rate and light. Fits a non-rectangular hyperbola of the form A ~ ((phi * Q + Amax - sqrt((phi * Q + Amax)^2 - 4 * theta * phi * Q * Amax)). A = assimilation rate, Q = PPFD, light level, phi = initial slope, theta = curvature, and Rd = respiration.

For documentation see the online help:

	help(package = "FitAQ")

Installation:

	devtools::install_github("MarkusLoew/FitAQ")


## FitAQ
Fits a non-rectangular hyperbola to assimilation rate vs light. Uses nls with the "port" algorithm by default. Other nls algorithms can be selected as well see ?nls for details. It is possible to restrict the curve fit paramters to stay within a specified range (e.g. the curvature theta is supposed to be between 0 and 1) using the "upper" and "lower" arguments of nls when the "port" algorithm is used. Custom start values for the four curve fit parameters can be provided if the buit-in estimate for the start values is failing. See example below.

## FitLCP
Calculates the light compensation point from a model provided by FitAQ.

## FitSat
Calculates the light level at which a specific ratio of the maximum assimilation rate is reached. By default the degree of saturation is 0.9, i.e. the function calculates the light level at which 90% of Amax is reached. Uses the model calculated by FitAQ.

## CalcAQ
This can be seen as the reverse of FitAQ: Given the four curve-fit parameters Amax, phi, Rd, and theta, the response curve is calculated.

### Examples
From the help file for FitAQ:

```r
 # A light response that can be fit with default settings
 p <- structure(list(Photo = c(-1.47735760973694, 0.607640418824763, 
 2.25569085178172, 4.57715647757407, 10.3609389085487, 14.9802351366038, 
 16.9201279124267, 19.2992412684557, 19.6156433699216, 19.5987435815687
 ), PARi = c(0, 25, 50, 100, 250, 500, 750, 1250, 1500, 2000)), class = "data.frame",
 row.names = c(NA, 10L))

 # do the AQ fit
 FitAQ(data = p, A = Photo, Q = PARi)

 # Example how to provide custom start values
 my.start <- list(Amax = 20, phi = 0.05, Rd = 5, theta = 0.9)
 FitAQ(data = p, A = Photo, Q = PARi, start = my.start)

 # Instead of just the coefficients, return the actual model
 model <- FitAQ(data = p, A = Photo, Q = PARi, start = my.start, provide.model = TRUE)
 # predict model results
 predict_range <- data.frame(Q = seq(0, 10000, length = 1000))
 my.line <- within(predict_range, A <- predict(model, newdata = predict_range))
 
 plot(Photo ~ PARi, data = p)
 lines(A ~ Q, data = my.line, col = "red")

 # Another plant where fit needs to be tweaked
 q <- structure(list(Photo = c(-5.1305665308579, -1.80737285839408, 
 1.29971850585751, 4.94585314505446, 13.9726643801131, 23.6419602864848, 
 28.9662664715206, 34.9085507079998, 37.3274268032741, 39.473492407329
 ), PARi = c(0, 25, 50, 100, 250, 500, 750, 1250, 1500, 2000)), class = "data.frame", 
 row.names = 1597:1606)

 # creating a figure with wide margins to visualise differences between fits
 plot(Photo ~ PARi, data = q, xlim = c(0, 4000), ylim = c(0, 50))
 
 FitAQ(q, A = Photo, Q = PARi) # results in negative theta
 
 model.3 <- FitAQ(q, A = Photo, Q = PARi, provide.model = TRUE)
 my.line.3 <- within(predict_range, A <- predict(model.3, newdata = predict_range))
 lines(A ~ Q, data = my.line.3, col = "orange")

 # re-fit with constraints and start values:
 model.4 <- FitAQ(q, A = Photo, Q = PARi, 
                  start = c(Amax = 50, phi = 0.09, Rd = 4, theta = 0.4),
                  lower = c(0, 0, 0, 0),
                  upper = c(56, 0.1, 10, 1),
                  provide.model = TRUE,
                  trace = TRUE)

 my.line.4 <- within(predict_range, A <- predict(model.4, newdata = predict_range))
 lines(A ~ Q, data = my.line.4, col = "green")
 ```
