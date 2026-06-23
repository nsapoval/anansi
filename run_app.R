#!/usr/bin/env Rscript
# Spin up the anansi Shiny app locally from this source clone (no install).
#
#   Rscript run_app.R            # default port 8100, opens your browser
#   Rscript run_app.R 8200       # custom port
#
# Stop the app with Ctrl-C.

args <- commandArgs(trailingOnly = TRUE)
port <- if (length(args) >= 1L && nzchar(args[[1]])) as.integer(args[[1]]) else 8100L
if (is.na(port)) stop("Port must be an integer, e.g. Rscript run_app.R 8200",
                      call. = FALSE)

if (!requireNamespace("shiny", quietly = TRUE)) {
  stop('The app needs shiny. Install it with: install.packages("shiny")',
       call. = FALSE)
}
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop('Loading the source clone needs devtools. ',
       'Install it with: install.packages("devtools")', call. = FALSE)
}

# Optional UI niceties; the app runs without them (text colour fields, no SVG).
optional <- c("bslib", "colourpicker", "svglite")
missing <- optional[!vapply(optional, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  message("Optional UI packages not installed (app still runs): ",
          paste(missing, collapse = ", "), "\n  install.packages(c(",
          paste(sprintf('"%s"', missing), collapse = ", "), "))")
}

# Load the package from this directory so local edits take effect immediately.
suppressMessages(devtools::load_all(".", quiet = TRUE))

message(sprintf("anansi app -> http://127.0.0.1:%d  (Ctrl-C to stop)", port))
run_anansi_app(port = port, host = "127.0.0.1", launch.browser = TRUE)
