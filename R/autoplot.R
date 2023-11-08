#' Autoplot
#'
#' @param object HB model
#' @param ... other argument
#' @return plot
#'
#' @examples
#' \dontrun{
#'  model1 <- hb_beta(
#'     y ~ x1 + x2 + x3,
#'     data = unemployment, seed = 1,
#'     burn.in = 2000, iter.mcmc = 10000, iter.update = 10,
#'  )
#' autoplot(model1)
#' }
#'
#' @export
autoplot.saehb <- function(object, ...){
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  graphics::par(mar = c(2, 2, 2, 2))
  coda::autocorr.plot(object$result_mcmc, col = "brown2", lwd = 2)
  plot(object$result_mcmc, col = "brown2", lwd = 2)
}

#' @inherit ggplot2::autoplot
#' @export
autoplot <- function(object, ...) {
  UseMethod("autoplot")
}
