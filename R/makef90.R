#' @export
makef90 <- function(df, col_ID = "ID", col_serial = "ID", col_TIME = "TIME", col_BM = "CMT", cols_COVT = NULL, cols_COVY = NULL, plotmax = 40, plotmin = -15) {
  code_1st <- system.file("template/pred_Sreft_1st.f90", package = "sreft")
  code_1st <- readLines(code_1st)
  code_1st <- paste0(code_1st, collapse = "\n")

  code_2nd <- system.file("template/pred_Sreft_2nd.f90", package = "sreft")
  code_2nd <- readLines(code_2nd)
  code_2nd <- paste0(code_2nd, collapse = "\n")

  code_3rd <- system.file("template/pred_Sreft_3rd.f90", package = "sreft")
  code_3rd <- readLines(code_3rd)
  code_3rd <- paste0(code_3rd, collapse = "\n")

  code_range <- paste0("      real (kind=dpsize) :: plotmax = ", plotmax, "d0, plotmin = ", plotmin, "d0\n")

  code_prms <- paste0("      id = datrec(", which(names(df) == col_ID), ")\n",
                      "      serial = datrec(", which(names(df) == col_serial), ")\n",
                      "      time = datrec(", which(names(df) == col_TIME), ")\n",
                      "      bm = datrec(", which(names(df) == col_BM), ")\n",
                      "      meanx = datrec((/", paste(grep("MeanX", names(df)), collapse = ","), "/)\n",
                      "      meany = datrec((/", paste(grep("MeanY", names(df)), collapse = ","), "/))\n",
                      "      coun = datrec((/", paste(grep("Count", names(df)), collapse = ","), "/))\n",
                      "      a = theta(1:numbm) + eta(1:numbm)\n",
                      "      b = theta(1+numbm:2*numbm) + eta(1+numbm:2*numbm)\n",
                      "      c = theta(1+2*numbm:3*numbm) + eta(1+2*numbm:3*numbm)\n")

  code_covt <- "      covt = 1.0d0"
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
  code_covt <- paste0(code_covt, "\n")

  code_covy <- "      covy = 0.0d0"
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
  code_covy <- paste0(code_covy, "\n")

  output <- paste0(list(code_1st, code_range, code_2nd, code_prms, code_covt, code_covy, code_3rd), collapse = "\n")

  return(output)
}
