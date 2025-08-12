makef90 <- function(df, col_ID = "ID", col_serial = "ID", col_TIME = "TIME", col_BM = "CMT", cols_COVT = NULL, plotmax = 40, plotmin = -15) {
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

  code_prms <- paste0("      id = datrac(", which(names(df) == col_ID), ")\n",
                      "      serial = datrac(", which(names(df) == col_serial), ")\n",
                      "      time = datrac(", which(names(df) == col_TIME), ")\n",
                      "      bm = datrac(", which(names(df) == col_BM), ")\n",
                      "      meanx = datrac(/", paste(grep("MeanX", names(df)), collapse = ","), "/)\n",
                      "      meany = datrac(/", paste(grep("MeanY", names(df)), collapse = ","), "/)\n",
                      "      coun = datrac(/", paste(grep("Count", names(df)), collapse = ","), "/)\n",
                      "      a = theta(1:numbm) + eta(1:numbm)\n",
                      "      b = theta(1+numbm:2*numbm) + eta(1+numbm:2*numbm)\n",
                      "      c = theta(1+2*numbm:3*numbm) + eta(1+2*numbm:3*numbm)\n",
                      "      covt = 1.0d0\n",
                      "      covy = 0.0d0\n")

  # if (!is.null(cols_COVT)) {
  #   which(names(df) %in% cols_COVT)
  # }

  output <- paste0(list(code_1st, code_range, code_2nd, code_prms, code_3rd), collapse = "\n")

  return(output)
}
