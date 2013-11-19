unquote <- function(x) {
  gsub("\"(.*?)\"", "\\1", x)
}

parse.matlab.matrix <- function(s) {
  s <- gsub("\\[(.*?)\\]", "\\1", s, perl = TRUE)
  rows <- strsplit(s, "\\s*;\\s*", perl = TRUE)[[1]]
  mx <- strsplit(rows, "\\s*,\\s*", perl = TRUE)

  for (row in seq_along(mx)) {
    mx[[row]] <- unquote(mx[[row]])
  }

  mx
}

read.pccsv <- function(file, col.names) {
  ncolumn <- length(col.names)

  # Get file content
  lines <- readLines(file)
  # Trim lines
  lines <- str_trim(lines)

  # Remove comments
  lines <- lines[grep("^#", lines, perl = TRUE, invert = TRUE)]

  # Remove empty lines
  lines <- lines[grep("^\\s*$", lines, perl = TRUE, invert = TRUE)]

  N <- length(lines)

  mx <- matrix(NA, N, ncolumn)
  nn <- list()
  for (i in seq_along(lines)) {
    row <- as.numeric(str_split(lines[i], ",\\s*")[[1]])
    mx[i, ] <- row[1:ncolumn]
    nn[[i]] <- row[-(1:(ncolumn + 1))]
  }

  df <- as.data.frame(mx)
  df[[1]] <- as.integer(df[[1]])

  attr(df, "nn") <- nn

  names(df) <- col.names
  df
}

read.pce <- function(file) {
  # Get file content
  lines <- readLines(file)
  # Trim lines
  lines <- str_trim(lines)
  # Extract header
  lines <- lines[grep("^#", lines, perl = TRUE)]
  # Remove comments from header
  lines <- lines[grep("^##", lines, perl = TRUE, invert = TRUE)]
  # Remove leading sharpes
  lines <- sub("\\s*#+\\s*", "", lines, perl = TRUE)

  # Split lines by equality sign
  splitted.lines <- strsplit(lines, "\\s*=\\s*", perl = TRUE)

  # Extract and set names
  prop.names <- sapply(splitted.lines, function(el) el[1])
  props <- lapply(splitted.lines, function(el) paste0(el[-1], collapse = ""))
  prop.names <- gsub("_", ".", prop.names, fixed = TRUE)
  names(props) <- prop.names

  # Create empty list for result
  info <- list()

  info$name <- unquote(props$name)
  info$original.percent <- as.integer(unquote(props$original.percent))
  info$stage.percent <- as.integer(unquote(props$stage.percent))
  info$pointcloud.quality <- as.integer(unquote(props$pointcloud.quality))
  info$stage <- unquote(props$stage)
  info$phenotype.string <- unquote(props$phenotype.string)
  info$genotype.id <- as.integer(unquote(props$genotype.id))
  info$phenotype.number <- as.integer(unquote(props$phenotype.number))
  info$automated.quality <- if (is.null(props$automated.quality)) {
        props$pointcloud.quality
      } else {
        as.integer(unquote(props$automated.quality))
      }
  info$nuclear.count <- as.integer(unquote(props$nuclear.count))
  info$nuclear.stain <- unquote(props$nuclear.stain)


  info$col.names <- parse.matlab.matrix(props$column)[[1]]

  info$column.info <- parse.matlab.matrix(props$column.info)[[1]]

  info$column.info.bid <- parse.matlab.matrix(props$column.info.bid)
  info$column.info.bid <- lapply(info$column.info.bid, function(el) {
        names(el) <- c("dye", "stain", "gene", "gene.full.name", "data.type", "default.measurement")
        as.list(el)
      })
  genes <- sapply(info$column.info.bid, function(el) el$gene)
  dyes <- sapply(info$column.info.bid, function(el) el$dye)
  names(dyes) <- genes
  names(info$column.info.bid) <- genes
  info$genes <- genes
  info$dyes <- dyes

  info$channel.offset <- as.numeric(parse.matlab.matrix(props$channel.offset)[[1]])
  info$channel.gain <- as.numeric(parse.matlab.matrix(props$channel.gain)[[1]])
  names(info$channel.offset) <- names(info$channel.gain) <- dyes

  intensity.correction.raw <- parse.matlab.matrix(props$intensity.correction)
  info$intensity.correction <- sapply(intensity.correction.raw, function(el) el[[2]])
  names(info$intensity.correction) <- sapply(intensity.correction.raw, function(el) el[[1]])

  info$attenuation.offset <- as.numeric(parse.matlab.matrix(props$attenuation.offset)[[1]])
  names(info$attenuation.offset) <- parse.matlab.matrix(props$attenuation.correct)[[1]]

  # Store raw header in info
  info$header <- props

  pce <- read.pccsv(file, info$col.names)
  info$col.names <- NULL

  attr(pce, "info") <- info

  invisible(pce)
}

extract.gene.pce <- function(pce, gene,
                             measurement = c("default", "apical", "basal", "nuclear", "cellular")) {
  info <- attr(pce, "info")
  gene <- match.arg(gene, info$genes)

  other.genes <- info$genes[info$genes != gene][-1]
  other.dyes <- info$dyes[other.genes]

  column.info.bid <- info$column.info.bid[[gene]]

  dye <- column.info.bid$dye
  default.measurement <- column.info.bid$default.measurement
  measurement <- match.arg(measurement)
  if (identical(measurement, "default"))
    measurement <- default.measurement
  intensity.col <- if (measurement == "") dye else paste(dye, measurement, sep = "_")
  intensity.col <- gsub(" ", ".", intensity.col, fixed = TRUE)
  pce <- pce[, c("x", "y", "z", intensity.col)]
  names(pce) <- c("x3d", "y3d", "z3d", "values")

  # TODO Decide what to do with addittional info. Drop it or document it?
  info$genes <- NULL
  info$gene <- gene
  info$dyes <- NULL
  info$dye <- dye

  info$channel.offset <- info$channel.offset[dye]
  info$channel.gain <- info$channel.gain[dye]
  info$intensity.correction <- info$intensity.correction[dye]
  info$attenuation.offset <- info$attenuation.offset[dye]
  info$column.info.bid <- column.info.bid
  info$measurement <- measurement

  info$channel.offset <- info$channel.offset[dye]
  info$channel.gain <- info$channel.gain[dye]
  info$attenuation.offset <- info$attenuation.offset[dye]
  info$intensity.correction <- info$intensity.correction[dye]

  info$other.genes <- paste0(other.genes, collapse = ",")
  info$other.dyes <- paste0(other.dyes, collapse = ",")

  attr(pce, "info") <- info

  pce <- as.list(pce)

  class(pce) <- "embryo3d"
  invisible(pce)
}

read.pce.gene.old <- function(file, gene) {
  Zcol <- if (is.character(gene)) {
    gene <- match.arg(gene, c("eve", "Kr"))
    switch(gene, eve = 13, Kr = 17)
  } else {
    gene
  }

  ncol <- max(count.fields(file, sep = ",", comment.char = "#"))
  pce <- read.table(file, sep = ",", comment.char = "#", fill = TRUE,
      col.names = paste("V", seq_len(ncol), sep = ""))

  out <- pce[, c(2:4, Zcol)]
  names(out) <- c("x3d", "y3d", "z3d", "values")
  out <- as.list(out)

  class(out) <- "embryo3d"
  out
}

read.pce.gene <- function(file, gene, ...) {
  extract.gene.pce(read.pce(file), gene = gene, ...)
}
