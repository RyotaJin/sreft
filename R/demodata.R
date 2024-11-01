check_lengths <- function(...) {
  lengths <- sapply(list(...), length)
  if (length(unique(lengths)) != 1) {
    stop("Error: All vectors must have the same length.")
  }
}

getObs <- function(x, sd_res) {
  a_ <- x[paste0("a", x["BM"])]
  b_ <- x[paste0("b", x["BM"])]
  c_ <- x[paste0("c", x["BM"])]
  t_ <- x["TIME2"]
  output <- a_ + b_ / c_ * (exp(c_ * t_) - 1) + rnorm(1, 0, sd_res[as.numeric(x["BM"])])

}

makeDemodata <- function(n_sub, prms_a, prms_b, prms_c, sd_a, sd_b, sd_c, sd_res, min_time, max_time, timepoint, isspread = TRUE) {
  check_lengths(prms_a, prms_b, prms_c, sd_a, sd_b, sd_c)
  n_bm <- length(prms_a)

  df_offsetT <- data.frame(ID = 1:n_sub,
                           offsetT = runif(n_sub, min_time, max_time))

  for (i in 1:n_bm) {
    df_offsetT <- df_offsetT %>%
      mutate(!!paste0("a", i) := prms_a[i] + rnorm(n_sub, 0, sd_a[i]),
             !!paste0("b", i) := prms_b[i] + rnorm(n_sub, 0, sd_b[i]),
             !!paste0("c", i) := prms_c[i] + rnorm(n_sub, 0, sd_c[i]))
  }

  df <- expand.grid(ID = 1:n_sub,
                    TIME = timepoint,
                    BM = 1:n_bm) %>%
    arrange(ID, TIME) %>%
    left_join(df_offsetT, by = "ID") %>%
    mutate(TIME2 = TIME + offsetT) %>%
    mutate(DV = apply(., 1, getObs, sd_res = sd_res))

  if (isspread) {
    df <- df %>%
      mutate(BM_str = paste0("Biomarker", BM)) %>%
      select(ID, TIME, TIME2, BM_str, DV, offsetT) %>%
      pivot_wider(names_from  = BM_str, values_from = DV) %>%
      arrange(ID, TIME)
  }

  return(df)
}
