#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleR)
  library(SummarizedExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("usage: run_singler.R REFERENCE_H5AD QUERY_H5AD OUTPUT_CSV LABEL_KEY")
}
reference_path <- args[[1]]
query_path <- args[[2]]
output_path <- args[[3]]
label_key <- args[[4]]

reference <- readH5AD(reference_path, use_hdf5 = FALSE, reader = "R")
query <- readH5AD(query_path, use_hdf5 = FALSE, reader = "R")
if (!(label_key %in% colnames(colData(reference)))) {
  stop(sprintf("reference H5AD lacks label key: %s", label_key))
}

choose_assay <- function(x) {
  if ("logcounts" %in% assayNames(x)) return("logcounts")
  if ("X" %in% assayNames(x)) return("X")
  assayNames(x)[[1]]
}
reference_assay <- choose_assay(reference)
query_assay <- choose_assay(query)
assay(reference, "logcounts") <- assay(reference, reference_assay)
assay(query, "logcounts") <- assay(query, query_assay)

shared <- intersect(rownames(reference), rownames(query))
if (length(shared) < 50) stop("fewer than 50 shared genes")
reference <- reference[shared, , drop = FALSE]
query <- query[shared, , drop = FALSE]
prediction <- SingleR(
  test = query,
  ref = reference,
  labels = as.character(colData(reference)[[label_key]]),
  assay.type.test = "logcounts",
  assay.type.ref = "logcounts",
  de.method = "classic"
)
label <- as.character(prediction$pruned.labels)
label[is.na(label) | label == ""] <- "Unknown"
score <- apply(as.matrix(prediction$scores), 1, max, na.rm = TRUE)
output <- data.frame(
  cell_id = colnames(query),
  predicted_label = label,
  similarity = score,
  score_semantics = "maximum_spearman_similarity_with_pruning_not_probability",
  stringsAsFactors = FALSE
)
write.csv(output, output_path, row.names = FALSE, quote = TRUE)
