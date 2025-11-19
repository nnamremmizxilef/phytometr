#' Plot host model diagnostics
#'
#' @param host_model Object returned by \code{define_host_model()}.
#'
#' @return A ggplot object with a simple textual summary.
#' @export
plot_host_model <- function(host_model) {

  if (is.null(host_model$type)) {
    stop("`host_model` must have a 'type' element (e.g. 'clone' or 'natural').")
  }

  if (host_model$type == "clone") {
    text <- "Clone system:\nHost variance = 0\n(no host-driven noise)"
  } else {
    text <- sprintf(
      "Natural population\nHost variance fraction = %.2f\nStructure = %s",
      host_model$var_host,
      host_model$pop_structure
    )
  }

  df <- data.frame(x = 1, y = 1, label = text)

  ggplot2::ggplot(df, ggplot2::aes(x, y, label = label)) +
    ggplot2::geom_text(size = 5) +
    ggplot2::theme_void()
}
