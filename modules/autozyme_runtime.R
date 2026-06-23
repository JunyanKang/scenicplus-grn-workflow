activate_autozyme_r <- function(families = c("seurat")) {
  status <- list()
  if (identical(Sys.getenv("AUTOZYME_DISABLED", unset = "0"), "1")) {
    message("[autozyme] disabled by AUTOZYME_DISABLED=1")
    status[["autozyme"]] <- "disabled_by_env"
    return(status)
  }
  if (!requireNamespace("autozyme", quietly = TRUE)) {
    message("[autozyme] R package not available")
    status[["autozyme"]] <- "not_available"
    return(status)
  }
  for (family in families) {
    ok <- tryCatch({
      autozyme::activate(family)
      TRUE
    }, error = function(e) {
      message("[autozyme] could not activate ", family, ": ", conditionMessage(e))
      FALSE
    })
    status[[family]] <- if (isTRUE(ok)) "active" else "not_active"
  }
  status
}
