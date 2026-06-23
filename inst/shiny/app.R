# anansi Shiny web app -- launch with anansi::run_anansi_app().
#
# Upload a CSV (or paste/upload raw extended-Newick, or load the bundled
# example), tune the densinet() controls, preview the figure, and download it.
# The input -> densinet() mapping lives in anansi::build_densinet() so this file
# stays a thin UI shell. bslib / colourpicker / svglite are optional: the app
# degrades gracefully when they are not installed.

library(shiny)

has_bslib <- requireNamespace("bslib", quietly = TRUE)
has_cp    <- requireNamespace("colourpicker", quietly = TRUE)
has_svg   <- requireNamespace("svglite", quietly = TRUE)

# A colour input that uses colourpicker when available, else a plain text field.
colorInput <- function(id, label, value) {
  if (has_cp) colourpicker::colourInput(id, label, value)
  else        textInput(id, paste0(label, " (name or #hex)"), value)
}

# ---- UI pieces -------------------------------------------------------------

input_controls <- tagList(
  fileInput("file", "Upload networks CSV", accept = c(".csv")),
  tags$details(
    tags$summary("...or paste / upload extended-Newick"),
    textAreaInput("enewick_text", NULL,
                  placeholder = "One extended-Newick network per line",
                  height = "120px"),
    fileInput("enewick_file", "Upload .nwk / .txt (one network per line)",
              accept = c(".nwk", ".txt", ".tre", ".newick")),
    actionButton("load_enewick", "Parse pasted/uploaded networks")
  ),
  actionButton("load_example", "Load bundled example"),
  verbatimTextOutput("summary"),
  tags$hr()
)

key_controls <- tagList(
  selectInput("method", "Tip-order method",
              c("mode", "mds", "closest_leaf", "first")),
  radioButtons("reticulation_style", "Reticulation style",
               c("arrow", "hybrid"), inline = TRUE),
  radioButtons("mode", "Mode", c("cladogram", "phylogram"), inline = TRUE),
  radioButtons("layout", "Cloud layout", c("slanted", "rectangular"),
               inline = TRUE),
  sliderInput("consensus_p", "Consensus threshold (p)", 0, 1, 0.5, 0.05),
  numericInput("top_n", "Top-N networks (0 = all)", 0, min = 0, step = 1),
  selectInput("top_by", "Rank top-N by", choices = c("(auto)" = "")),
  selectizeInput("outgroup", "Outgroup taxa", choices = character(0),
                 multiple = TRUE),
  radioButtons("outgroup_position", "Outgroup position",
               c("top", "bottom"), inline = TRUE),
  selectizeInput("keep", "Keep taxa (subset; blank = all)",
                 choices = character(0), multiple = TRUE),
  textInput("title", "Title (blank = auto)", ""),
  checkboxInput("tip_labels", "Show tip labels", TRUE)
)

advanced_controls <- tagList(
  checkboxInput("snap_to_consensus", "Snap tip order to consensus", FALSE),
  checkboxInput("color_by_support", "Colour backbone by clade support", TRUE),
  checkboxInput("consensus_ret", "Draw consensus reticulations", TRUE),
  checkboxInput("consistent_only", "Drop divergent-taxa networks", TRUE),
  sliderInput("consensus_ret_min", "Min reticulation frequency", 0, 1, 0.1, 0.05),
  sliderInput("ret_edge_frac", "Reticulation edge anchor frac", 0, 1, 0.1, 0.05),
  numericInput("linewidth", "Cloud linewidth", 0.3, min = 0, step = 0.1),
  numericInput("consensus_linewidth", "Consensus linewidth", 0.7, min = 0,
               step = 0.1),
  numericInput("jitter", "Cloud jitter", 0, min = 0, step = 0.005),
  numericInput("tip_size", "Tip-label size", 3, min = 0, step = 0.5),
  numericInput("tip_offset", "Tip-label offset", 0.02, step = 0.01),
  numericInput("alpha", "Cloud alpha (blank = auto)", NA, min = 0, max = 1,
               step = 0.01),
  numericInput("ret_alpha", "Reticulation alpha (blank = auto)", NA, min = 0,
               max = 1, step = 0.01),
  selectInput("ret_linetype", "Reticulation line type",
              c("dashed", "dotted", "solid", "dotdash", "longdash")),
  colorInput("tree_color", "Backbone cloud colour", "steelblue"),
  colorInput("ret_color", "Reticulation cloud colour", "firebrick"),
  colorInput("consensus_color", "Consensus backbone colour", "black"),
  colorInput("consensus_ret_color", "Consensus reticulation colour", "darkred")
)

advanced_panel <- if (has_bslib) {
  bslib::accordion(
    open = FALSE,
    bslib::accordion_panel("Advanced options", advanced_controls))
} else {
  tags$details(tags$summary("Advanced options"), advanced_controls)
}

download_controls <- tagList(
  tags$hr(),
  selectInput("dl_format", "Download format",
              c("PNG", "PDF", if (has_svg) "SVG")),
  fluidRow(
    column(6, numericInput("dl_w", "Width (in)", 8, min = 1)),
    column(6, numericInput("dl_h", "Height (in)", 6, min = 1))),
  numericInput("dl_dpi", "DPI (raster)", 150, min = 36),
  downloadButton("download", "Download figure")
)

sidebar_content <- tagList(input_controls, key_controls, advanced_panel,
                           download_controls)

ui <- if (has_bslib) {
  bslib::page_sidebar(
    title = "anansi — densinet",
    sidebar = bslib::sidebar(width = 380, sidebar_content),
    plotOutput("plot", height = "680px"))
} else {
  fluidPage(
    titlePanel("anansi — densinet"),
    sidebarLayout(
      sidebarPanel(sidebar_content, width = 4),
      mainPanel(plotOutput("plot", height = "680px"), width = 8)))
}

# ---- Server ----------------------------------------------------------------

server <- function(input, output, session) {
  rv <- reactiveValues(netset = NULL)

  # Load a netset, surfacing parse errors / divergent-taxa warnings as
  # non-fatal notifications instead of crashing the app.
  set_netset <- function(expr) {
    tryCatch(
      withCallingHandlers(
        rv$netset <- expr,
        warning = function(w) {
          showNotification(conditionMessage(w), type = "warning", duration = 8)
          invokeRestart("muffleWarning")
        }),
      error = function(e)
        showNotification(conditionMessage(e), type = "error", duration = 10))
  }

  observeEvent(input$file, {
    req(input$file)
    set_netset(anansi::read_networks_csv(input$file$datapath))
  })

  # An uploaded .nwk/.txt fills the text box; the user then clicks "Parse".
  observeEvent(input$enewick_file, {
    req(input$enewick_file)
    txt <- tryCatch(paste(readLines(input$enewick_file$datapath, warn = FALSE),
                          collapse = "\n"),
                    error = function(e) "")
    updateTextAreaInput(session, "enewick_text", value = txt)
  })

  observeEvent(input$load_enewick, {
    txt <- input$enewick_text
    if (is.null(txt) || !nzchar(trimws(txt))) {
      showNotification("Paste or upload extended-Newick first.", type = "message")
      return()
    }
    set_netset(anansi::netset_from_enewick(txt))
  })

  observeEvent(input$load_example, {
    set_netset(anansi::read_networks_csv(anansi::anansi_example()))
  })

  # Populate taxa and numeric-metadata selectors whenever a set loads.
  observeEvent(rv$netset, {
    ns <- rv$netset
    req(ns)
    tx <- anansi::network_taxa(ns)
    updateSelectizeInput(session, "outgroup", choices = tx, selected = character(0),
                         server = TRUE)
    updateSelectizeInput(session, "keep", choices = tx, selected = character(0),
                         server = TRUE)
    num_cols <- if (!is.null(ns$meta))
      names(ns$meta)[vapply(ns$meta, is.numeric, logical(1))] else character(0)
    updateSelectInput(session, "top_by", choices = c("(auto)" = "", num_cols))
  })

  output$summary <- renderText({
    ns <- rv$netset
    if (is.null(ns)) return("No data loaded yet.")
    sprintf("%d networks | %d taxa | %d taxa-consistent",
            length(ns), length(ns$taxa), sum(ns$taxa_ok))
  })

  current_params <- reactive(list(
    method = input$method, reticulation_style = input$reticulation_style,
    mode = input$mode, layout = input$layout, consensus_p = input$consensus_p,
    top_n = input$top_n, top_by = input$top_by,
    outgroup = input$outgroup, outgroup_position = input$outgroup_position,
    keep = input$keep, title = input$title, tip_labels = input$tip_labels,
    snap_to_consensus = input$snap_to_consensus,
    color_by_support = input$color_by_support, consensus_ret = input$consensus_ret,
    consistent_only = input$consistent_only,
    consensus_ret_min = input$consensus_ret_min, ret_edge_frac = input$ret_edge_frac,
    linewidth = input$linewidth, consensus_linewidth = input$consensus_linewidth,
    jitter = input$jitter, tip_size = input$tip_size, tip_offset = input$tip_offset,
    alpha = input$alpha, ret_alpha = input$ret_alpha,
    ret_linetype = input$ret_linetype, tree_color = input$tree_color,
    ret_color = input$ret_color, consensus_color = input$consensus_color,
    consensus_ret_color = input$consensus_ret_color))

  plot_obj <- reactive({
    ns <- rv$netset
    validate(need(ns, "Upload a CSV, paste extended-Newick, or load the example to begin."))
    tryCatch(
      anansi::build_densinet(ns, current_params()),
      error = function(e) validate(need(FALSE, conditionMessage(e))))
  })

  output$plot <- renderPlot(plot_obj())

  output$download <- downloadHandler(
    filename = function() paste0("anansi_densinet.", tolower(input$dl_format)),
    content = function(file) {
      p <- plot_obj()
      dev <- switch(input$dl_format, PNG = "png", PDF = "pdf", SVG = "svg")
      ggplot2::ggsave(file, p, device = dev, width = input$dl_w,
                      height = input$dl_h, dpi = input$dl_dpi, bg = "white")
    })
}

shinyApp(ui, server)
