check_lengths <- function(...) {
  lengths <- sapply(list(...), length)
  if (length(unique(lengths)) != 1) {
    stop("All vectors must have the same length.")
  }
}

getObs <- function(x, var_res) {
  a_ <- x[paste0("a", x["BM"])]
  b_ <- x[paste0("b", x["BM"])]
  c_ <- x[paste0("c", x["BM"])]
  t_ <- x["TIME2"]
  output <- a_ + b_ / c_ * (exp(c_ * t_) - 1) + rnorm(1, 0, sqrt(var_res[as.numeric(x["BM"])]))

}

#' @export
makeDemodata <- function(n_sub, prms_a, prms_b, prms_c, var_a, var_b, var_c, var_res, min_time, max_time, timepoint, isspread = TRUE) {
  check_lengths(prms_a, prms_b, prms_c, var_a, var_b, var_c)
  n_bm <- length(prms_a)

  df_indprms <- data.frame(ID = 1:n_sub,
                           offsetT = runif(n_sub, min_time, max_time))
  for (i in 1:n_bm) {
    df_indprms[[paste0("a", i)]] <- prms_a[i] + rnorm(n_sub, 0, sqrt(var_a[i]))
    df_indprms[[paste0("b", i)]] <- prms_b[i] * exp(rnorm(n_sub, 0, sqrt(var_b[i])))
    df_indprms[[paste0("c", i)]] <- prms_c[i] * exp(rnorm(n_sub, 0, sqrt(var_c[i])))
  }

  df <- expand.grid(ID = 1:n_sub,
                    TIME = timepoint,
                    BM = 1:n_bm)
  df <- df[order(df$ID, df$TIME), ]
  df <- merge(df, df_indprms, by = "ID", sort = FALSE)
  df$TIME2 <- df$TIME + df$offsetT
  df$DV <- apply(df, 1, function(row) getObs(row, var_res = var_res))

  if (isspread) {
    df$BM_str <- paste0("Biomarker", df$BM)

    wide_df <- unique(df[c("ID", "TIME", "TIME2", "offsetT")])
    wide_df <- wide_df[order(wide_df$ID, wide_df$TIME), ]

    for (bm in unique(df$BM_str)) {
      tmp <- df[df$BM_str == bm, c("ID", "TIME", "DV")]
      names(tmp)[3] <- bm
      wide_df <- merge(wide_df, tmp, by = c("ID", "TIME"), all.x = TRUE, sort = FALSE)
    }

    df <- wide_df
  }

  return(df)
}
