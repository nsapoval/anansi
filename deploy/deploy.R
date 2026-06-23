#!/usr/bin/env Rscript
# deploy/deploy.R --------------------------------------------------------------
#
# Deploy the bundled anansi Shiny app (inst/shiny/app.R) to shinyapps.io.
#
# The app is a thin UI over the anansi package, so the build server installs
# `anansi` itself (from GitHub) plus its CRAN/Bioconductor dependency closure.
# rsconnect builds that dependency manifest from the *locally installed* library,
# so whatever environment runs this script must already have anansi (installed
# from GitHub, not local source) and its deps available -- the GitHub Actions
# workflow (.github/workflows/deploy-shinyapps.yaml) handles that for CI.
#
# Credentials come from the environment so the same script works in CI (repo
# secrets) and locally:
#   SHINYAPPS_NAME    shinyapps.io account/username
#   SHINYAPPS_TOKEN   token  (dashboard -> Account -> Tokens)
#   SHINYAPPS_SECRET  secret (revealed via "Show secret")
#
# Locally you can instead authenticate once interactively with
#   rsconnect::setAccountInfo(name = ..., token = ..., secret = ...)
# and then just run:  SHINYAPPS_NAME=<account> Rscript deploy/deploy.R

if (nzchar(Sys.getenv("SHINYAPPS_TOKEN"))) {
  rsconnect::setAccountInfo(
    name   = Sys.getenv("SHINYAPPS_NAME"),
    token  = Sys.getenv("SHINYAPPS_TOKEN"),
    secret = Sys.getenv("SHINYAPPS_SECRET"))
}

acct <- Sys.getenv("SHINYAPPS_NAME")

rsconnect::deployApp(
  appDir      = "inst/shiny",
  appName     = "anansi",
  appTitle    = "anansi — densinet",
  account     = if (nzchar(acct)) acct else NULL,
  forceUpdate = TRUE,
  logLevel    = "verbose")
