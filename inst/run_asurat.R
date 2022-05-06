#!/usr/bin/env Rscript

library(optparse)
library(ASURAT)

# Parse arguments from command line
options <- list(
  make_option(c("-a", "--args1"), action = "store", default = NA, type = "character", help="Character args"),
  make_option(c("-u", "--args2"), action = "store", default = NA, type = "numeric", help="Numeric args")
)

arguments <- parse_args(OptionParser(option_list = options))

run_asurat(arg1 = arguments$gff)
