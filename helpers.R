#' Fetch supplemetary files from GEO
#' @importFrom XML htmlParse xpathSApply
#' @importFrom httr GET parse_url
#' @export
FetchGEOFiles <- function(geo, download.dir = getwd(), download.files = FALSE, ...) {
  geo <- trimws(toupper(geo))
  geo_type <- substr(geo, 1, 3)
  url.prefix <- "https://ftp.ncbi.nlm.nih.gov/geo/"
  if (geo_type == "GSE") {
    url.prefix <- paste0(url.prefix, "series/")
  } else if (geo_type == "GSM") {
    url.prefix <- paste0(url.prefix, "samples/")
  } else if (geotype == "GPL") {
    url.prefix <- paste0(url.prefix, "platform/")
  }


  # url.prefix <- "https://ftp.ncbi.nlm.nih.gov/geo/series/"

  geo_prefix <- paste0(substr(x = geo, start = 1, stop = nchar(geo) - 3), "nnn")
  url <- paste0(url.prefix, geo_prefix, "/", geo, "/", "suppl", "/")

  response <- GET(url = url)
  html_parsed <- htmlParse(file = response)
  links <- xpathSApply(doc = html_parsed, path = "//a/@href")
  suppl_files <- as.character(grep(pattern = "^G", x = links, value = TRUE))
  if (length(suppl_files) == 0) {
    return(NULL)
  }

  file.url <- paste0(url, suppl_files)
  file_list <- data.frame(filename = suppl_files, url = file.url)

  if (download.files) {
    names(file.url) <- suppl_files
    download_file <- function(url, filename, ...) {
      message(paste0("Downloading ", filename, " to ", download.dir))
      download.file(url = url, destfile = file.path(download.dir, filename), mode = "wb", ...)
      message("Done!")
    }
    lapply(seq_along(file.url), function(y, n, i) {
      download_file(y[[i]], n[[i]], ...)
    },
    y = file.url, n = names(file.url)
    )
  }

  return(file_list)
}


#' Load in data from remote or local mtx files
#' Adapted and inspired from Seurat's Read10X
#'
#' Enables easy loading of sparse data matrices
#'
#' @param mtx Name or remote URL of the mtx file
#' @param cells Name or remote URL of the cells/barcodes file
#' @param features Name or remote URL of the features/genes file
#' @param feature.column Specify which column of features files to use for feature/gene names; default is 2
#' @param cell.column Specify which column of cells file to use for cell names; default is 1
#' @param unique.features Make feature names unique (default TRUE)
#' @param strip.suffix Remove trailing "-1" if present in all cell barcodes.
#'
#' @return A sparse matrix containing the expression data.
#'
#' @importFrom Matrix readMM
#' @importFrom utils read.delim
#' @importFrom httr build_url parse_url
#' @importFrom tools file_ext
#'
#'
#' @export
#' @concept preprocessing
#'
#' @examples
#' \dontrun{
#' # For local files:
#'
#' expression_matrix <- ReadMtx(genes = "count_matrix.mtx.gz", features = "features.tsv.gz", cells = "barcodes.tsv.gz")
#' seurat_object <- CreateSeuratObject(counts = expression_matrix)
#'
#' # For remote files:
#'
#' expression_matrix <- ReadMtx(
#'   mtx = "http://localhost/matrix.mtx",
#'   cells = "http://localhost/barcodes.tsv",
#'   features = "http://localhost/genes.tsv"
#' )
#' seurat_object <- CreateSeuratObject(counts = data)
#' }
#'
ReadMtx <- function(mtx,
                    cells,
                    features,
                    cell.column = 1,
                    feature.column = 2,
                    unique.features = TRUE,
                    strip.suffix = FALSE) {
  mtx <- build_url(url = parse_url(url = mtx))
  cells <- build_url(url = parse_url(url = cells))
  features <- build_url(url = parse_url(url = features))
  all_files <- list("Expression matrix" = mtx, "Barcode" = cells, "Gene name" = features)

  check_file_exists <- function(filetype, filepath) {
    if (grepl(pattern = "^:///", x = filepath)) {
      filepath <- gsub(pattern = ":///", replacement = "", x = filepath)
      if (!file.exists(paths = filepath)) {
        stop(paste(filetype, "file missing. Expecting", filepath), call. = FALSE)
      }
    }
  }

  # check if all files exist
  lapply(seq_along(all_files), function(y, n, i) {
    check_file_exists(n[[i]], y[[i]])
  }, y = all_files, n = names(all_files))

  # convenience fucntion to read local or remote tab delimited files
  readTableUri <- function(uri) {
    if (grepl(pattern = "^:///", x = uri)) {
      uri <- gsub(pattern = ":///", replacement = "", x = uri)
      textcontent <- read.table(file = uri, header = FALSE, sep = "\t", row.names = NULL)
    } else {
      if (file_ext(uri) == "gz") {
        textcontent <- read.table(
          file = gzcon(url(uri), text = TRUE),
          header = FALSE, sep = "\t", row.names = NULL
        )
      } else {
        textcontent <- read.table(
          file = uri, header = FALSE,
          sep = "\t", row.names = NULL
        )
      }
    }
    return(textcontent)
  }

  # read barcodes
  cell.barcodes <- readTableUri(uri = cells)
  bcols <- ncol(x = cell.barcodes)
  if (bcols < cell.column) {
    stop(paste0(
      "cell.column was set to ", cell.column,
      " but ", cells, " only has ", bcols, " columns.",
      " Try setting the cell.column argument to a value <= to ", bcols, "."
    ))
  }
  cell.names <- cell.barcodes[, cell.column]

  if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
    cell.names <- as.vector(x = as.character(x = sapply(
      X = cell.names,
      FUN = ExtractField,
      field = 1,
      delim = "-"
    )))
  }

  # read features
  feature.names <- readTableUri(uri = features)
  fcols <- ncol(x = feature.names)
  if (fcols < feature.column) {
    stop(paste0(
      "feature.column was set to ", feature.column,
      " but ", features, " only has ", fcols, " column(s).",
      " Try setting the feature.column argument to a value <= to ", fcols, "."
    ))
  }
  if (any(is.na(x = feature.names[, feature.column]))) {
    na.features <- which(x = is.na(x = feature.names[, feature.column]))
    replacement.column <- ifelse(test = feature.column == 2, yes = 1, no = 2)
    if (replacement.column > fcols) {
      stop(
        paste0(
          "Some features names are NA in column ", feature.column,
          ". Try specifiying a different column."
        ),
        call. = FALSE,
        immediate. = TRUE
      )
    } else {
      warning(
        paste0(
          "Some features names are NA in column ", feature.column,
          ". Replacing NA names with ID from column ", replacement.column, "."
        ),
        call. = FALSE,
        immediate. = TRUE
      )
    }
    feature.names[na.features, feature.column] <- feature.names[na.features, replacement.column]
  }

  feature.names <- feature.names[, feature.column]
  if (unique.features) {
    feature.names <- make.unique(names = feature.names)
  }

  # read mtx
  if (grepl(pattern = "^:///", x = mtx)) {
    mtx <- gsub(pattern = ":///", replacement = "", x = mtx)
    data <- readMM(mtx)
  } else {
    if (file_ext(mtx) == "gz") {
      data <- readMM(gzcon(url(mtx)))
    } else {
      data <- readMM(mtx)
    }
  }

  colnames(x = data) <- cell.names
  rownames(x = data) <- feature.names

  return(data)
}
