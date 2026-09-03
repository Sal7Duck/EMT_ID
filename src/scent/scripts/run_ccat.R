#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
    stop(
        "Usage: Rscript run_ccat.R expression.tsv ppi.rda output.tsv"
    )
}

expression_file <- args[1]
ppi_file <- args[2]
output_file <- args[3]

library(SCENT)

# ---------------------------------------------------------
# Load expression matrix
# Rows    = human Entrez IDs
# Columns = cell IDs
# ---------------------------------------------------------

exp.m <- as.matrix(
    read.delim(
        expression_file,
        header = TRUE,
        row.names = 1,
        check.names = FALSE
    )
)

# ---------------------------------------------------------
# Load official SCENT PPI
# ---------------------------------------------------------

load(ppi_file)

# SCENT's official object should be net13Jun12.m
if (!exists("net13Jun12.m")) {
    stop("net13Jun12.m was not found in the PPI RDA file.")
}

# ---------------------------------------------------------
# Run official SCENT CCAT
# ---------------------------------------------------------

ccat <- CompCCAT(
    exp.m = exp.m,
    ppiA.m = net13Jun12.m
)

# CompCCAT returns an unnamed vector
result <- data.frame(
    cell_id = colnames(exp.m),
    CCAT_score = as.numeric(ccat)
)

write.table(
    result,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat(
    "CCAT completed for",
    nrow(result),
    "cells\n"
)
