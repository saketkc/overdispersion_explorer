# source("imports.R")
suppressPackageStartupMessages({
  library(shiny)
  library(shinyjs)
  library(shinyWidgets)
  library(shinydashboard)
  library(httr)
  library(XML)
  library(dplyr)
  library(ggplot2)
  library(sctransform)
  library(patchwork)
  library(scales)
  library(sparseMatrixStats)
  library(reshape2)
  library(cowplot)
  library(ggrepel)
  library(tools)
  library(scattermore)
  library(Matrix)
})


source("helpers.R")
source("plots.R")
theme_set(theme_bw())

ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "GEO Overdispersion Explorer"),
  dashboardSidebar(
    width = 350,
    useShinyjs(),
    searchInput(
      inputId = "gse_id",
      label = "Enter GSE/GSM/GPL",
      value = "GSE132044", # "GSM4191941",#  "GSM3564450", #"GSM3564443",# "GSM3823986", # "GSE132044",
      btnSearch = icon("search"), btnReset = icon("remove"),
      width = "100%"
    ),
    uiOutput("file_options"),
    actionButton("go", "Run")
  ),
  dashboardBody(
    box(plotOutput("mean_var_plot", height = "800px", width = "800px"))
  )
)

server <- function(input, output) {
  GetCM <- eventReactive(input$go, {
    # ( (!is.null(input$mtx_file)) & (!is.null(input$feature_file)) & (!is.null(input$barcode_file)) )
    disable("go")
    geo_files_df <- FetchGEOFiles(input$gse_id)
    rownames(geo_files_df) <- geo_files_df$filename
    cm <- ReadMtx(
      mtx = geo_files_df[input$mtx_file, "url"],
      cells = geo_files_df[input$barcode_file, "url"],
      features = geo_files_df[input$feature_file, "url"],
      feature.column = as.numeric(input$feature_column),
      cell.column = as.numeric(input$cell_column)
    )
    cm <- as(object = cm, Class = "dgCMatrix")
    cm
  })

  output$file_options <- renderUI({
    geo_files_df <- FetchGEOFiles(input$gse_id)
    if (is.null(geo_files_df)) {
      return("Nothing found")
    }
    geo_filenames <- as.list(geo_files_df$filename)
    # names(geo_filenames) <- geo_files_df$url
    selected_mtx_file <- NULL
    selected_barcode_file <- NULL
    selected_feature_file <- NULL

    is_mtx <- grepl(".mtx|matrix|count", geo_filenames)
    if (sum(is_mtx) > 0) {
      selected_mtx_file <- geo_filenames[is_mtx][[1]]
    }
    is_gene <- grepl("gene|feature", geo_filenames)
    if (sum(is_gene) > 0) {
      selected_feature_file <- geo_filenames[is_gene][[1]]
    }
    is_barcode <- grepl("barcode|cell", geo_filenames)
    if (sum(is_barcode) > 0) {
      selected_barcode_file <- geo_filenames[is_barcode][[1]]
    }

    m <- pickerInput(
      inputId = "mtx_file",
      label = "Matrix file",
      choices = geo_filenames,
      selected = selected_mtx_file
    )
    f <- pickerInput(
      inputId = "feature_file",
      label = "Features file",
      choices = geo_filenames,
      selected = selected_feature_file
    )
    b <- pickerInput(
      inputId = "barcode_file",
      label = "Barcodes file",
      choices = geo_filenames,
      selected = selected_barcode_file
    )
    cc <- pickerInput(
      inputId = "cell_column",
      label = "Cell column",
      choices = c(1, 2)
    )
    fc <- pickerInput(
      inputId = "feature_column",
      label = "Feature column",
      choices = c(1, 2)
    )
    vm <-
      pickerInput(
        inputId = "vst_method",
        label = "SCT Method",
        choices = c("glmGamPoi2", "glmGamPoi", "offset-10", "offset-100", "offset-Inf")
      )
    list(m, f, b, cc, fc, vm)
  })
  nsteps <- 6
  output$mean_var_plot <- renderPlot({
    withProgress(message = "Status: ", value = 0, {
      incProgress(1 / nsteps, detail = "Fetching Data from GEO")
      out <- tryCatch(
        {
          cm <- GetCM()
        },
        error = function(err) {
          enable("go")
          message(err)
          return(NULL)
        }
      )
      incProgress(2 / nsteps, detail = "Plotting Mean Variance")
      p1 <- PlotMeanVar(cm)
      incProgress(3 / nsteps, detail = "Plotting Zeros Fraction")
      p2 <- PlotZerosFit(cm)
      incProgress(4 / nsteps, detail = "Running SCT")
      if (grepl("offset", input$vst_method)) {
        theta_given <- as.numeric(strsplit(input$vst_method, "-")[[1]][[2]])
        vst.out <- sctransform::vst(umi = cm, method = "offset", theta_given = theta_given, n_genes = 2000, n_cells = 5000, verbosity = 0)
      } else {
        vst.out <- sctransform::vst(umi = cm, method = input$vst_method, n_genes = 2000, n_cells = 5000, verbosity = 0)
      }
      gene_attr <- as.data.frame(vst.out$gene_attr)
      gene_attr$gene <- rownames(gene_attr)
      p3 <- residualVarPlot(gene_attr) + ggtitle("")
      p4 <- sctransform::plot_model_pars(vst.out, show_theta = TRUE)
      incProgress(5 / nsteps, detail = "Done!")
      enable("go")
      return((p1 | p2) / p3 / p4)
    })
  })
}

shinyApp(ui = ui, server = server)
