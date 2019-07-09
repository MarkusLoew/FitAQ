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
#' Farquhar, Graham D., and Thomas D. Sharkey. ‘Stomatal Conductance and Photosynthesis’. Annual Review of Plant Physiology 33, no. 1 (1982): 317–345.
#'
#' Lambers, Hans, F. Stuart Chapin III, and Thijs L. Pons. Plant Physiological Ecology. Springer Science & Business Media, 2008.
#'
#' Thornley, J. H. M. Mathematical Models in Plant Physiology : A Quantitative Approach to Problems in Plant and Crop Physiology. London; New York: Academic Press, 1976.
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
     print(start.phi)
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

