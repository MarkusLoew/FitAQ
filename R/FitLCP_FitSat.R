#' Calculate light compensation point and light saturation point from non-rectangular hyperbola model of a photosynthetic light reponse
#'
#' @param model Model object from either \code{FitAQ} or \code{nls}
#' @param coef.Amax Character string with the name of the curve fit coefficient representing maximum photosynthetic assimilation rate Amax. Default is "Amax".
#' @param coef.phi Character string with the name of the curve fit coefficient representing the initial slope phi. Default is "phi".
#' @param coef.theta Character string with the name of the curve fit coefficient representing the curvature of the hyperbola. Default is "theta".
#' @param coef.Rd Character string with the name of the curve fit coefficient representing the respiration rate Rd. Default is "Rd".
#' @return Light compensation point. Named numeric.
#' @export

FitLCP <- function(model, coef.Amax = "Amax", coef.phi = "phi", coef.theta = "theta", coef.Rd = "Rd") {
  
  Amax  <- stats::coef(model)[coef.Amax]
  phi   <- stats::coef(model)[coef.phi]
  theta <- stats::coef(model)[coef.theta]
  Rd    <- stats::coef(model)[coef.Rd]
  
  # non-rectangular hyperbola function to substitute the model coefficients
  x <- function(x) {
       (1 / (2 * theta)) * (phi * x + Amax - sqrt((phi * x + Amax)^2 - 4 * phi * theta * Amax * x)) - Rd}
  
  # search for zero (= light compensation point) in the interval of 0 to 100 PPFD
  out <- stats::uniroot(x, c(0, 100))$root
  names(out) <- "LCP"
 return(out)
}


#' Calculate the light level at which the photosynthetic assimilation is saturated
#'
#' Given a odel of a non-recangular curve fit of assimilation rate and light, this funciton calculates the light level at which a specific proportion of Amax is reached.
#'
#' @param model Model object created by either \code{FitAQ} or \code{nls}.
#' @param sat.fac Degree of saturation of Amax. Defaults to 0.9, i.e. 90\% of Amax.
#' @param range Numeric list with two entries: First entry is the lower end of a window for the light level at which the assimilation rate is saturated. Second entry is the upper end of the window.
#' @param coef.Amax Character string with the name of the curve fit coefficient representing maximum photosynthetic assimilation rate Amax. Default is "Amax".
#' @param coef.phi Character string with the name of the curve fit coefficient representing the initial slope phi. Default is "phi".
#' @param coef.theta Character string with the name of the curve fit coefficient representing the curvature of the hyperbola. Default is "theta".
#' @param coef.Rd Character string with the name of the curve fit coefficient representing the respiration rate Rd. Default is "Rd".
#' @return Light level at which the desired level of Amax is reached. Named numeric.
#' @export

FitSat <- function(model, sat.fac = 0.9, range = c(0, 2200), 
                   coef.Amax = "Amax", coef.phi = "phi", coef.theta = "theta", coef.Rd = "Rd") {
  Amax  <- stats::coef(model)[coef.Amax]
  phi   <- stats::coef(model)[coef.phi]
  theta <- stats::coef(model)[coef.theta]
  Rd    <- stats::coef(model)[coef.Rd]
  
  # non-rectangular hyperbola function to substitute the model coefficients including the level of saturation
  
  x <- function(x) {
        #(1 / (2 * theta)) * (phi * x + Amax - sqrt((phi * x + Amax)^2 - 4 * phi * theta * Amax * x)) - Rd - (0.75 * Amax)+ 0.75 * (Rd)}
        (1 / (2 * theta)) * (phi * x + Amax - sqrt((phi * x + Amax)^2 - 4 * phi * theta * Amax * x)) - Rd - (sat.fac * Amax)+ sat.fac * (Rd)}
  
  # search for sat.fac (= light saturation point) in the interval of 0 to 5000 PPFD
  #out <- uniroot(x, c(0, 5000))$root
  out <- stats::uniroot(x, range)$root
  names(out) <- "LSP"
  return(out)
}

