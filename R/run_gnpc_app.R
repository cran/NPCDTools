#' Launch the GNPC Shiny app
#'
#' Launches the interactive demo shipped with the package (can be found in
#' \code{inst/shiny/GNPC_app}). The app demonstrates NPC, GNPC, and G-DINA
#' workflows and the ECPE real-data example.
#'
#' @param launch.browser Logical; open in a web browser? Defaults to \code{interactive()}.
#' @param host Host interface passed to \code{shiny::runApp}. Defaults to "127.0.0.1".
#' @param port Optional integer port. If \code{NULL}, Shiny picks a free port.
#' @param ...  Additional arguments forwarded to \code{shiny::runApp()}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   run_gnpc_app()
#' }
run_gnpc_app <- function(launch.browser = interactive(),
                         host = "127.0.0.1",
                         port = NULL,
                         ...) {
  # 1) Locate the app directory installed with the package
  app_dir <- system.file("shiny", "GNPC_app", package = utils::packageName())
  if (identical(app_dir, "") || !dir.exists(app_dir)) {
    stop("Could not find 'inst/shiny/GNPC_app' in the installed package. ",
         "Try reinstalling the package.", call. = FALSE)
  }
  
  # 2) Check runtime dependencies used by the app (keep CDM in Suggests if desired)
  req_pkgs <- c("shiny", "shinythemes", "GDINA", "CDM")
  missing <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing required package(s): ", paste(missing, collapse = ", "),
         ". Please install them and try again.", call. = FALSE)
  }
  
  # 3) Helpful banner: version + namespace path (useful when you have a local dev build)
  ver  <- as.character(utils::packageVersion(utils::packageName()))
  nsp  <- getNamespaceInfo(asNamespace(utils::packageName()), "path")
  msg  <- sprintf("\nLaunching GNPC app from %s\nPackage %s (%s)\n",
                  normalizePath(app_dir, winslash = "/"),
                  utils::packageName(), ver)
  message(msg)
  
  # 4) Run the app
  shiny::runApp(appDir = app_dir,
                launch.browser = launch.browser,
                host = host,
                port = port,
                display.mode = "normal",
                ...)
}
