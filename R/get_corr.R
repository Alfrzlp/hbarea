#' Get correlation between y and x variabl
#'
#' @param y a numeric vector
#' @param x a numeric matrix or data.frame
#' @param method a character string indicating which correlation coefficient is to be used for the test. One of "pearson", "kendall", or "spearman", can be abbreviated.
#' @param alpha significance level maximum
#' @param filter_cor filter minimum correlation
#' @param round integer indicating the number of decimal places
#'
#' @return data.frame contains two columns: estimation correlation and p-value of the t-test.
#' @export
#'
#' @examples
#' get_corr(unemployment$y, x = unemployment[-c(1:2)])
get_corr <- function(y, x, method = 'pearson', alpha = NULL, filter_cor = NULL, round = 3){
  res <- apply(
    X = x,
    MARGIN = 2,
    FUN = function(kolom){
      res <- stats::cor.test(y, kolom, method = method)
      return(c(res$estimate, pvalue = res$p.value))
    }
  )

  res <- as.data.frame(t(res))
  res <- round(dplyr::arrange(res, "pvalue"), round)

  if (!is.null(alpha)) {
    res <- stats::na.omit(res[res$pvalue <= alpha, ])
  }
  if (!is.null(filter_cor)) {
    res <- stats::na.omit(res[abs(res$cor) >= filter_cor, ])
  }
  return(res)
}
