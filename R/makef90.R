fill_template <- function(template_file, values) {
  txt <- readLines(template_file, warn = FALSE)
  for (key in names(values)) {
    pattern <- paste0("\\{", key, "\\}")
    txt <- gsub(pattern, values[[key]], txt)
  }
  paste(txt, collapse = "\n")
}

#' @export
makef90 <- function(df, no_mainbm, runno = "run000", col_ID = "ID", col_serial = "ID", col_TIME = "TIME", col_BM = "CMT", cols_COVT = NULL, cols_COVY = NULL, plotmax = 40, plotmin = -15) {
  code_covt <- ""
  if (!is.null(cols_COVT)) {
    colnos_COVT <- which(names(df) %in% cols_COVT)
    if (length(cols_COVT) != length(colnos_COVT)) {
      stop("One or more specified column names do not exist in the dataset.")
    }
    numbm <- max(df$CMT)
    for (i in seq_along(colnos_COVT)) {
      code_covt <- paste0(code_covt, " + theta(", 3 * numbm + i, ") * datrec(", colnos_COVT[i], ")")
    }
  }

  code_covy <- ""
  if (!is.null(cols_COVY)) {
    colnos_COVY <- which(names(df) %in% cols_COVY)
    if (length(cols_COVY) != length(colnos_COVY)) {
      stop("One or more specified column names do not exist in the dataset.")
    }
    numbm <- max(df$CMT)
    numbm_thetas <- numbm * 3 + length(cols_COVT)
    for (i in seq_along(colnos_COVY)) {
      code_covy <- paste0(code_covy, " + theta(", numbm_thetas + numbm * (i - 1), "+bm) * datrec(", colnos_COVY[i], ")")
    }
  }

  output <- fill_template(system.file("template/pred_sreft.f90", package = "sreft"),
                          list(plotmax = plotmax,
                               plotmin = plotmin,
                               runno = runno,
                               no_mainbm = no_mainbm,
                               id_num = which(names(df) == col_ID),
                               serial_num = which(names(df) == col_serial),
                               time_num = which(names(df) == col_TIME),
                               bm_num = which(names(df) == col_BM),
                               meanx_num = paste(grep("MeanX", names(df)), collapse = ","),
                               meany_num = paste(grep("MeanY", names(df)), collapse = ","),
                               coun_num = paste(grep("Count", names(df)), collapse = ","),
                               covt_code = code_covt,
                               covy_code = code_covy
                          ))

  return(output)
}
