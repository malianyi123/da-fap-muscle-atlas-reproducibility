#!/usr/bin/env Rscript

root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), "..", ".."), mustWork = TRUE)
out <- file.path(root, "environment", "sessionInfo.txt")

sink(out)
cat("R sessionInfo\n")
cat("=============\n\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")
print(sessionInfo())
sink()

