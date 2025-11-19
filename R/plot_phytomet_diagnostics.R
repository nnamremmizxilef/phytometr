#' Plot diagnostics for a holobiont design
#'
#' @param x Object of class \code{"phytomet_diagnostics"}
#' @param which Character, which plot(s) to show:
#'   \itemize{
#'     \item \code{"overview"} (default): 2x2 layout with env-space, PCA,
#'           component correlation and environment-holobiont correlation.
#'     \item \code{"env"}: environmental space plot.
#'     \item \code{"pca"}: PCA of holobiont components.
#'     \item \code{"components"}: holobiont component correlation heatmap.
#'     \item \code{"env_effects"}: environment-holobiont correlation heatmap.
#'     \item \code{"variance"}: variance decomposition barplot (if present).
#'   }
#' @param ... Not used, for compatibility.
#'
#' @return A ggplot or patchwork object (invisibly).
#' @export
plot.phytomet_diagnostics <- function(x,
                                      which = c("overview",
                                                "env",
                                                "pca",
                                                "components",
                                                "env_effects",
                                                "variance"),
                                      ...) {

  which <- match.arg(which)

  # convenience aliases
  p_env   <- x$plots$env_space
  p_pca   <- x$plots$pca
  p_comp  <- x$plots$comp_corr
  p_envhf <- x$plots$env_effects
  p_var   <- x$plots$variance  # may be NULL if you didn't store it

  # helper for "single" panels
  show_single <- function(p) {
    if (is.null(p)) {
      stop("Requested panel is not available in this diagnostics object.")
    }
    print(p)
    invisible(p)
  }

  if (which != "overview") {
    if (which == "env")         return(show_single(p_env))
    if (which == "pca")         return(show_single(p_pca))
    if (which == "components")  return(show_single(p_comp))
    if (which == "env_effects") return(show_single(p_envhf))
    if (which == "variance")    return(show_single(p_var))
  }

  # --- overview layout ----------------------------------------------
  # Try to use patchwork if available for a nice 2x2 grid
  have_patchwork <- requireNamespace("patchwork", quietly = TRUE)

  if (have_patchwork) {
    # if variance plot exists, put it under PCA
    if (!is.null(p_var)) {
      top_row    <- p_env | (p_pca / p_var)
    } else {
      top_row    <- p_env | p_pca
    }
    bottom_row <- p_comp | p_envhf

    p <- top_row / bottom_row

    # add a global title if you like
    p <- p + patchwork::plot_annotation(
      title = "phytometr design diagnostics",
      subtitle = "Environment, holobiont structure and environment-holobiont links"
    )

    print(p)
    return(invisible(p))
  }

  # Fallback: no patchwork installed, show a message and print sequentially
  message(
    "Package 'patchwork' is not installed; ",
    "showing diagnostic panels sequentially.\n",
    "Install patchwork for a combined overview: install.packages('patchwork')."
  )

  print(p_env)
  print(p_pca)
  if (!is.null(p_var)) print(p_var)
  print(p_comp)
  print(p_envhf)

  invisible(x)
}
