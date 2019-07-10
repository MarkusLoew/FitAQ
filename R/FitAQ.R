#' Flexible function to calculate a non-rectangular hyperbola fit for a photosynthetic light response curve.
# 
#' @description Fit a non-rectangular hyperbola to the photosynthetic light response (assimilation rate as a function of light) following Lambers et al. 2008, Thornley 1976, Farquhar & Sharkey 1982: \code{(A ~ ((phi * Q + Amax - sqrt((phi * Q + Amax)^2 - 4 * theta * phi * Q * Amax)) / 2 * theta) - Rd)} Uses \code{nls} for the curve fit using the algorithm \code{"port"}. Starting values for the initial curve fit (Amax, phi, Rd, theta) can be provided as a named list. Otherwise, a simple estimation of starting parameters takes place.
# 
#' @param data Data frame with the light response gas exchange data. Needs to provide vectors for light level and assimilation rate
#' @param A Name of the vector in the data frame that holds the assimilation rate.
#' @param Q Name of the vector in the data frame that holds the light level information.
#' @param nlscontrol Optional object that controls internal settings of nls. See ?nls.control. Default is NULL.
#' @param start Optional named list that provides starting values for nls fit of Amax, phi, Rd, theta. (Note: Rd is a positive value)
#' @param provide.model Logical. If FALSE (the default) the model coefficients are returned. When TRUE, the model is returned.
#' @param ... Other options to the underlying functions.
#' @return Depending on the value of provide.model, the data frame with the coefficients of a non-rectangular hyperbola (Amax, phi, Rd, theta), is returned. When provide.model = TRUE, the curve fit model is returned.
#' @references 
#' Farquhar, Graham D., and Thomas D. Sharkey. 'Stomatal Conductance and Photosynthesis'. Annual Review of Plant Physiology 33, no. 1 (1982): 317–345.
#'
#' Lambers, Hans, F. Stuart Chapin III, and Thijs L. Pons. Plant Physiological Ecology. Springer Science & Business Media, 2008.
#'
#' Thornley, J. H. M. Mathematical Models in Plant Physiology: A Quantitative Approach to Problems in Plant and Crop Physiology. London; New York: Academic Press, 1976.
#' @examples
#' \dontrun{
#' # A light response that can be fit with default settings
#' p <- structure(list(Photo = c(-1.47735760973694, 0.607640418824763, 
#' 2.25569085178172, 4.57715647757407, 10.3609389085487, 14.9802351366038, 
#' 16.9201279124267, 19.2992412684557, 19.6156433699216, 19.5987435815687
#' ), PARi = c(0, 25, 50, 100, 250, 500, 750, 1250, 1500, 2000)), class = "data.frame",
#' row.names = c(NA, 10L))

#' # do the AQ fit
#' FitAQ(data = p, A = Photo, Q = PARi)

#' # Example how to provide custom start values
#' my.start <- list(Amax = 20, phi = 0.05, Rd = 5, theta = 0.9)
#' FitAQ(data = p, A = Photo, Q = PARi, start = my.start)
#'
#' # Instead of just the coefficients, return the actual model
#' model <- FitAQ(data = p, A = Photo, Q = PARi, start = my.start, provide.model = TRUE)
#' 
#' # predict model results over a wide x-axis range
#' predict_range <- data.frame(Q = seq(from = 0, to = 4000, by = 50))
#' my.line <- within(predict_range, A <- predict(model, newdata = predict_range))
#' 
#' plot(Photo ~ PARi, data = p)
#' lines(A ~ Q, data = my.line, col = "red")
#'
#' # Another plant where fit needs to be tweaked
#' q <- structure(list(Photo = c(-5.1305665308579, -1.80737285839408, 
#' 1.29971850585751, 4.94585314505446, 13.9726643801131, 23.6419602864848, 
#' 28.9662664715206, 34.9085507079998, 37.3274268032741, 39.473492407329
#' ), PARi = c(0, 25, 50, 100, 250, 500, 750, 1250, 1500, 2000)), class = "data.frame", 
#' row.names = 1597:1606)

#' # creating a figure with wide margins to visualise differences between fits
#' plot(Photo ~ PARi, data = q, xlim = c(0, 4000), ylim = c(0, 50))
#' FitAQ(q, A = Photo, Q = PARi) # results in negative theta
#' model.3 <- FitAQ(q, A = Photo, Q = PARi, provide.model = TRUE)
#' my.line.3 <- within(predict_range, A <- predict(model.3, newdata = predict_range))
#' lines(A ~ Q, data = my.line.3, col = "orange")

#' # re-fit with constraints and start values:
#' model.4 <- FitAQ(q, A = Photo, Q = PARi, 
#'                  start = c(Amax = 50, phi = 0.09, Rd = 4, theta = 0.4),
#'                  lower = c(0, 0, 0, 0),
#'                  upper = c(56, 0.1, 10, 1),
#'                  provide.model = TRUE,
#'                  trace = TRUE)
#'
#' my.line.4 <- within(predict_range, A <- predict(model.4, newdata = predict_range))
#' lines(A ~ Q, data = my.line.4, col = "green")
#' }
#' @export


FitAQ <- function(data, 
                  A             = A, 
                  Q             = Q, 
                  nlscontrol    = NULL,
                  start         = NULL, 
                  provide.model = FALSE, 
                  ...) {
  
  # match the arguments - arguments A and Q are mandatory
  arguments  <- as.list(match.call())
  A          <- eval(arguments$A, data)
  Q          <- eval(arguments$Q, data)
  nlscontrol <- eval(arguments$nlscontrol)
  start      <- eval(arguments$start)

  if (is.null(start)) {
     # some initial ballpark estimates to get the curve fit procedure going
     start.Amax    <- max(A) - min(A)
     # calculate a linear model for the lower end of the AQ response to get a first estimate of phi
     start.phi.est  <- stats::lm(A[Q <= 110] ~ Q[Q<=110])
     start.phi     <- stats::coef(start.phi.est)[2] # get slope from linear model
     start.Rd      <- -min(A)
     # to ensure start.Rd is a positive number
     start.Rd      <- ifelse(start.Rd < 0, abs(start.Rd), start.Rd) 
     start.theta   <- 0.9
        # create start list for nls
	start = list(Amax = start.Amax,
	             phi  = start.phi,
		     Rd   = start.Rd,
		    theta = start.theta)
  } else {
  	start = start
  }

 # non-rectangular hyperbola equation to fit light response. 
 # Following Lambers et al. (2008) "Plant Physiologcal Ecology", page 27, eq. 7, Farquhar & Sharkey 1982, Thornley 1976
 eq.lightresponse <-  A ~ ((phi * Q + Amax - sqrt((phi * Q + Amax)^2 - 4 * theta * phi * Q * Amax))) / (2 * theta) - Rd
 
 # try to fit the model
 model.light <- try(stats::nls(formula = eq.lightresponse,
                          start = start,
                        control = nlscontrol, 
                      algorithm = "port",
                          #lower = c(0, 0, 0, 0), 
                          #upper = c(120, 0.4, 10, 1),
                          #lower = c(-Inf, -Inf, -Inf, 0), 
                          #upper = c(Inf, Inf, Inf, 1),
                          ...))
 
 # parse model output
 if(inherits(model.light, "try-error")) {
   message("No fit found")
   # create dummy output for failed fits
   out <- data.frame(Amax  = NA,
                     phi   = NA,
                     Rd    = NA,
                     theta = NA)
   model <- NULL
 } else {
   # message("It fits!")
   # create output data frame
   out <- data.frame(Amax  = stats::coef(model.light)[1],
                     phi   = stats::coef(model.light)[2],
                     Rd    = stats::coef(model.light)[3],
                     theta = stats::coef(model.light)[4])
   model <- model.light
   }
                     
 # final output
 if (provide.model) {
   return(model)
 } else {
   return(out)
 }
}

