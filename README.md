[![Build Status](https://travis-ci.org/MarkusLoew/FitAQ.svg?branch=master)](https://travis-ci.org/MarkusLoew/FitAQ)

# FitAQ
Yet another R package to analyse the relationship between photosynthetic assimilation rate and light. Fits a non-rectangular hyperbola of the form Photo ~ ((phi * PPFD + Amax - sqrt((phi * PPFD + Amax)^2 - 4 * theta * phi * PPFD * Amax)). Photo = assimilation rate, PPFD = light level, phi = initial slope, theta = curvature, and Rd = respiration.

For documentation see the online help:

	help(package = "FitAQ")

Installation:

	devtools::install_github("MarkusLoew/FitAQ")


## FitAQ
Fits a non-rectangular hyperbola to assimilation rate vs light. Uses \code{nls} with the "port" algorithm. It is possible to restrict the curve fit paramters to stay within a specified range (e.g. the curvature theta is supposed to be between 0 and 1).

## FitLCP
Calculates the light compensation point from a model provided by FitAQ.

## FitSat
Calculates the light level at which a specific ratio of the maximum assimilation rate is reached. By default the degree of saturation is 0.9, i.e. the function calculates the light leve at which 90% of Amax is reached. Uses the model calculated by FitAQ.

## CalcAQ
This can be seen as the reverse of FitAQ: Given the four curve-fit parameters Amax, phi, Rd, and theta, the response curve is calculated.

