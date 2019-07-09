#' Calculate photosynthetic light response from Amax, phi, Rd and theta.
#'
#' This function calculates the light response of the net assimilation rate based on for parameters. Input paramters are maximum photosynthetic rate Amax, initial quantum efficiency phi, respiration rate Rd, and curvature theta. Using a non-rectangular hyperbola based on these four parameters, the net assimilation rate at a series of light levels is calcualted. This is the reverse of the function \code{FitAQ}.
#' @param Amax Maximum net assimilation rate
#' @param phi Quantum efficiency of the photosynthetic light response. I.e. initial slope at low light.
#' @param Rd Respiration rate (positive value!).
#' @param theta Curvature of the non-retangular hyperbola.
#' @param Q Vector of light levels (PPFD). By default \code{seq(from = 0, to = 2000, by = 50))} is used.
#' @seealso \code{\link[FitAQ:FitAQ]{FitAQ}}
#' @export
#' @examples
#' x <- CalcAQ(24, 0.07, 2.26, 0.7)
#' plot(A ~ Q, data = x)


CalcAQ <- function(Amax, phi, Rd, theta, Q = seq(from = 0, to = 2000, by = 50)) {
      A <- ((phi * Q + Amax - sqrt((phi * Q + 
             Amax)^2 - 4 * theta * phi * Q * Amax)))/(2 * theta) - 
             Rd
      out <- data.frame(Q = Q, A = A)
      return(out)
}

