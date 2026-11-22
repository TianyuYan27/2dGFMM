src_dir <- local({
  ofile <- NULL
  for (i in seq_len(sys.nframe())) {
    candidate <- sys.frame(i)$ofile
    if (!is.null(candidate)) ofile <- candidate
  }
  if (is.null(ofile)) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) == 1) {
      ofile <- sub("^--file=", "", file_arg)
    }
  }
  dirname(normalizePath(ofile))
})

source(file.path(src_dir, "gfmm2d_gaussian.R"))
source(file.path(src_dir, "gfmm2d_binomial.R"))
source(file.path(src_dir, "utils.R"))
source(file.path(src_dir, "hypothesis_testing.R"))
