#' Score Converter
#'
#' This function takes a vector of values and a maximum value as inputs, and performs a score conversion.
#'
#' @param x_vec A numeric vector of values to be converted.
#' @param max The maximum value for the conversion.
#'
#' @return A vector of converted scores.
#'
#' @details The function applies a score conversion formula to each value in the input vector. It adds 0.5 to each value, divides it by (max + 1), and performs a transformation to convert it to a score.
#'
#' @examples
#' x <- c(1, 2, 3)
#' max_value <- 5
#' converted_scores <- scoreConverter(x, max_value)
#'
scoreConverter <- function (x_vec, max) {
  output <- (x_vec + 0.5) / (max + 1)
  output <- output / (1 - output)
  return(output)
}


lognormConverter <- function (df, selected_bm, XMAX) {
  output <- df[selected_bm]
  min_values <- output %>% apply(2, min, na.rm = TRUE)
  under_zero_indices <- which(min_values <= 0)

  if (any(min_values <= 0 & XMAX[under_zero_indices] <= 0)) {
    biomarkers_with_issue <- selected_bm[min_values <= 0 & XMAX <= 0]
    message_ <- paste(biomarkers_with_issue, collapse = ", ") %>%
      paste("has values under 0, but appropriate XMAX was not supplied. Please provide XMAX or convert the values before passing them to this function.")
    stop(message_)
  } else {
    if (any(under_zero_indices)) {
      df[names(under_zero_indices)] <- mapply(scoreConverter, df[names(under_zero_indices)], XMAX[under_zero_indices])
    }
    output <- log(output)
    output <- apply(output, 2, scale)
  }
  return(output)
}


slope <- function (x, y) {
  if (length(x) != length(y)) {
    stop("x and y is not same length!")
  }

  x_ <- x[!is.na(x) & !is.na(y)]
  y_ <- y[!is.na(x) & !is.na(y)]

  if (length(x_) < 2 | sd(x_) == 0) {
    return(NA)
  }
  if (sd(y_) == 0) {
    return(0)
  }
  return(cor(y_, x_) * (sd(y_) / sd(x_)))
}


intercept <- function (x, y) {
  x_ <- x[!is.na(x) & !is.na(y)]
  y_ <- y[!is.na(x) & !is.na(y)]

  return(mean(y_) - slope(x_, y_) * mean(x_))
}


#' @title Create data frame for NONMEM
#'
#' @description This function takes a data frame and selected biomarkers and covariates as inputs, and generates a NONMEM sheet.
#'
#' @param df The input data frame."ID" and "TIME" columns should be prepared in advance.
#' @param selected_bm A character vector specifying the selected biomarkers.
#' @param selected_cov A character vector specifying the selected covariates.
#' @param XMAX A numeric vector specifying the maximum values for biomarkers. Default is NULL.
#' @param lognorm A logical value indicating whether to perform log-normalization. Default is TRUE.
#'
#' @return A data frame for NONMEM.
#'
#' @details The function filters the data frame, performs calculations on biomarkers, and creates the NONMEM sheet with calculated means, counts, and other relevant columns.
#'
#' @examples
#' df <- read.csv("data.csv")
#' selected_bm <- c("biomarker1", "biomarker2")
#' selected_cov <- c("covariate1", "covariate2")
#' sheet <- makeNonmemSheet(df, selected_bm, selected_cov)
#'
#' @export
makeNonmemSheet <- function (df, selected_bm, selected_cov, XMAX = NULL, lognorm = TRUE){
  output <- df
  if (lognorm) {
    output[selected_bm] <- lognormConverter(df, selected_bm, XMAX)
  }

  numbm <- length(selected_bm)

  meanxs <- lapply(selected_bm, function (col) {
    output %>%
      filter(!is.na(!!sym({{col}}))) %>%
      summarise(meanx = mean(TIME), .by = ID)
  }) %>%
    bind_rows(.id = "Biomarker") %>%
    pivot_wider(names_from  = Biomarker, values_from = meanx) %>%
    setNames(c("ID", paste0("MeanX", 1:numbm)))

  meanys <- output %>%
    summarise(across(all_of(selected_bm), \ (x) mean(x, na.rm = TRUE)), .by = ID) %>%
    setNames(c("ID", paste0("MeanY", 1:numbm)))

  counts <- output %>%
    summarise(across(all_of(selected_bm), \ (x) sum(!is.na(x))), .by = ID) %>%
    setNames(c("ID", paste0("Count", 1:numbm)))

  output <- output %>%
    select(ID, TIME, all_of(selected_bm), all_of(selected_cov)) %>%
    merge(meanxs) %>%
    merge(meanys) %>%
    merge(counts) %>%
    pivot_longer(cols = all_of(selected_bm), names_to = "CMT_name", values_to = "DV") %>%
    drop_na(DV) %>%
    mutate(CMT = match(CMT_name, (selected_bm))) %>%
    select(ID, TIME, CMT, CMT_name, DV, everything()) %>%
    arrange(ID, CMT, TIME)

  return(output)
}

#' Set Initial Parameters
#'
#' This function sets the initial parameters based on the given data frame and biomarker information.
#'
#' @param df The data frame containing the biomarker data.
#' @param selected_bm A character vector specifying the selected biomarkers.
#' @param definition_bm The biomarker used for definition.
#' @param definition_value The value used for the definition of the biomarker.
#' @param XMAX An optional parameter specifying the maximum value for log-normal conversion.
#' @param estimated_mean_offsetT The estimated mean offset.
#' @param lognorm Logical value indicating whether log-normal conversion should be applied.
#'
#' @return A data frame containing the calculated initial parameters.
#'
#' @examples
#' setInitialPrms(df, selected_bm, definition_bm, definition_value)
#'
#' @export
setInitialPrms <- function (df, selected_bm, definition_bm, definition_value, XMAX = NULL, estimated_mean_offsetT = 10, lognorm = TRUE) {
  df_ <- df
  if (lognorm) {
    df_[selected_bm] <- lognormConverter(df_, selected_bm, XMAX)
  }

  numbm <- length(selected_bm)

  slopes <- lapply(selected_bm, function (col) {
    df_ %>%
      summarise(slope = slope(TIME, !!sym({{col}})), .by = ID)
  }) %>%
    bind_rows(.id = "Biomarker") %>%
    mutate(Biomarker = as.numeric(Biomarker))

  meanys <- df_ %>%
    summarise(across(all_of(selected_bm), \ (x) mean(x, na.rm = TRUE)), .by = ID) %>%
    pivot_longer(cols = -ID, names_to = "Biomarker", values_to = "meany") %>%
    mutate(Biomarker = match(Biomarker, selected_bm))

  .Mean_sp  <- split(meanys$meany, meanys$Biomarker)
  .Slope_sp <- split(slopes$slope, slopes$Biomarker)

  output <- data.frame(Biomarker = selected_bm) %>%
    mutate(slope = mapply(slope, .Mean_sp, .Slope_sp),
           intercept = mapply(intercept, .Mean_sp, .Slope_sp),
           meanslope = mapply(mean, .Slope_sp, na.rm = TRUE),
           ave = apply(df_[selected_bm], 2, mean, na.rm = TRUE),
           sd = apply(df_[selected_bm], 2, sd, na.rm = TRUE),
           omega_α = 0.0001,
           omega_β = 0,
           omega_γ = 0.0001,
           sigma = 0.1)
  output[output$Biomarker == definition_bm, c("omega_β")] <- 0.0001
  output[output$Biomarker == definition_bm, c("omega_α", "omega_γ")] <- 0

  output[, "α"] <- output[, "meanslope"] / exp(output[, "slope"] * estimated_mean_offsetT)
  if (lognorm) {
    output[output$Biomarker == definition_bm, "α"] <- (log(definition_value) - output[1, "ave"]) / output[1, "sd"]
  } else {
    output[output$Biomarker == definition_bm, "α"] <- definition_value
  }
  output[, "β"] <- output[, "intercept"] + output[, "α"] * output[, "slope"]
  output[, "γ"] <- output[, "slope"]

  output[, -1] <- signif(output[, -1], digits = 5)
  return(output)
}


#' @export
makeControlStream <- function (init, df, no_definition_bm, runno = "", PROBLEM = "", DATA = "data.csv") {
  ctl <- list()
  ctl["problem"] <- paste0("$PROBLEM ", PROBLEM, "\n")
  ctl["input"] <- paste0("$INPUT ", paste(names(df), collapse = " "),"\n")
  ctl["data"] <- paste0("$DATA ", DATA, "\nIGNORE=@\n")
  ctl["subroutine"] <- paste0("$SUBROUTINE PRED=../../", runno, "_pred_sreft.f90\n")

  tmp_alpha <- init[, "α"]
  if (no_definition_bm >= length(tmp_alpha)) {
    tmp_alpha <- c(tmp_alpha, "FIXED")
  } else {
    tmp_alpha <- c(tmp_alpha[1:no_definition_bm], "FIXED", tmp_alpha[(no_definition_bm + 1):length(tmp_alpha)])
  }

  ctl["theta"] <- paste0("$THETA\n" ,
                         paste(tmp_alpha, collapse = " "), "\n",
                         paste(init[, "β"], collapse = " "), "\n",
                         paste(init[, "γ"], collapse = " "), "\n")
  ctl["omega"] <- paste0("$OMEGA\n",
                         paste(init[, "omega_α"] %>% replace(. == 0, "0 FIXED"), collapse = " "), "\n",
                         paste(init[, "omega_β"] %>% replace(. == 0, "0 FIXED"), collapse = " "), "\n",
                         paste(init[, "omega_γ"] %>% replace(. == 0, "0 FIXED"), collapse = " "), "\n")
  ctl["sigma"] <- paste0("$SIGMA ", paste(init[, "sigma"], collapse = " "), "\n")
  ctl["estimation"] <- "$ESTIMATION METHOD=1 MAXEVAL=99999 PRINT=1 NOABORT SIGDIGITS=2 NOTITLE=1\n"
  ctl["covariance"] <- "$COVARIANCE UNCONDITIONAL MATRIX=S\n"
  ctl["table"] <- paste0("$TABLE ", paste(names(df), collapse = " "), " PRED CIPRED CWRES\nNOPRINT FORMAT=,F10.5 FILE=",
                         runno, ".fit ONEHEADER NOTITLE NOAPPEND")

  return(ctl)
}

#' @export
evalModel <- function (x, prms) {
  n_bm <- length(prms) / 3
  .time  <- x["TIME"]
  .bm    <- x["BM"]
  .alpha <- prms[paste0("alpha", 1:n_bm)]
  .beta  <- prms[paste0("beta", 1:n_bm)]
  .gamma <- prms[paste0("gamma", 1:n_bm)]

  if (.bm == 1) {
    output <- .alpha[.bm] + .beta[.bm] / .gamma[.bm] * (exp(.gamma[.bm] * .time) - 1)
  }else{
    output <- .alpha[.bm] + .beta[.bm] / .gamma[.bm] * (exp(.gamma[.bm] * .time) - 1)
  }

  return(output %>% as.numeric())
}
