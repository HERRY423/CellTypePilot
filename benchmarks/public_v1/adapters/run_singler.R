#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleR)
  library(SummarizedExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("usage: run_singler.R TRAIN_H5AD TEST_H5AD OUTPUT_CSV")
}
# Use zellkonverter's native R reader. The default Python reader lazily creates
# a hidden basilisk environment and is not an auditable dependency boundary for
# this benchmark release.
train <- readH5AD(args[[1]], use_hdf5 = FALSE, reader = "R")
test <- readH5AD(args[[2]], use_hdf5 = FALSE, reader = "R")
if (!("cell_type" %in% colnames(colData(train)))) {
  stop("fold-train H5AD lacks cell_type labels")
}

choose_assay <- function(x) {
  if ("logcounts" %in% assayNames(x)) return("logcounts")
  if ("X" %in% assayNames(x)) return("X")
  assayNames(x)[[1]]
}
train_assay <- choose_assay(train)
test_assay <- choose_assay(test)
assay(train, "logcounts") <- assay(train, train_assay)
assay(test, "logcounts") <- assay(test, test_assay)

shared <- intersect(rownames(train), rownames(test))
if (length(shared) < 50) stop("fewer than 50 shared genes")
train <- train[shared, , drop = FALSE]
test <- test[shared, , drop = FALSE]
prediction <- SingleR(
  test = test,
  ref = train,
  labels = as.character(colData(train)$cell_type),
  assay.type.test = "logcounts",
  assay.type.ref = "logcounts",
  de.method = "classic"
)
label <- as.character(prediction$pruned.labels)
label[is.na(label) | label == ""] <- "Unknown"
score <- apply(as.matrix(prediction$scores), 1, max, na.rm = TRUE)
confidence <- pmin(1, pmax(0, (score + 1) / 2))
output <- data.frame(
  cell_id = colnames(test),
  predicted_label = label,
  confidence = confidence,
  stringsAsFactors = FALSE
)
write.csv(output, args[[3]], row.names = FALSE, quote = TRUE)
