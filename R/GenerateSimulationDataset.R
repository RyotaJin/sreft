library(MASS)
library(dplyr)
library(tidyr)


slope <- function (x, y) {
  if (length(x) != length(y)) {
    stop("x and y is not same length!")
  }
  if (length(x) < 2 | sd(x) == 0) {
    return(NA)
  }
  if (sd(y) == 0) {
    return(0)
  }
  return(cor(y, x) * (sd(y) / sd(x)))
}

intercept <- function (x, y) {
  return(mean(y) - slope(x, y) * mean(x))
}

estimatecInitialParameter <- function (nmsheet, BASELINE, GTA = 10) {
  NUM_BM <- nmsheet$Biomarker %>% max()

  summary_eachbm <- nmsheet %>%
    dplyr::group_by(ID, Biomarker) %>%
    dplyr::summarise(Mean = mean(TIME),
                     Slope = slope(TIME, DV),
                     .groups = "drop") %>%
    dplyr::group_by(Biomarker) %>%
    dplyr::summarise(slopes = slope(Mean, Slope),
                     intercepts = intercept(Mean, Slope),
                     meanslope = mean(Slope),
                     .groups = "drop")

  output <- data.frame(Biomarker = paste0("Biomarker", 1:NUM_BM)) %>%
    dplyr::mutate(theta_alpha = c(BASELINE, summary_eachbm$meanslope[-1] / exp(summary_eachbm$slopes[-1] * GTA)),
                  theta_beta = summary_eachbm$intercepts + theta_alpha * summary_eachbm$slopes,
                  theta_gamma =  summary_eachbm$slopes,
                  omega_alpha = c(0, rep(0.0001, NUM_BM - 1)),
                  omega_beta = c(0.0001, rep(0, NUM_BM - 1)),
                  omega_gamma = c(0, rep(0.0001, NUM_BM - 1)),
                  sigma = 0.01)

  return(output)
}

makeControlStreamSimulationTemplate <- function (init, nmsheet) {
  NUM_BM <- nrow(init)

  ctl <- list()

  ctl["problem"] <- paste0("$PROBLEM simulation\n")
  ctl["input"] <- paste0("$INPUT ", paste(names(nmsheet), collapse = " "),"\n")
  ctl["data"] <- paste0("$DATA data.csv ignore=@\n")
  ctl["subroutine"] <- paste0("$SUBROUTINE PRED=pred_sreft.f90\n")
  ctl["theta"] <- paste0("$THETA\n" ,
                         paste(init[1, "theta_alpha"], "FIXED",
                               paste(init[2:NUM_BM, "theta_alpha"], collapse = " ")), "\n",
                         paste(init[, "theta_beta"], collapse = " "), "\n",
                         paste(init[, "theta_gamma"], collapse = " "), "\n")
  ctl["omega"] <- paste0("$OMEGA\n",
                         paste(init[, "omega_alpha"] %>% replace(. == 0, "0 FIXED"), collapse = " "), "\n",
                         paste(init[, "omega_beta"] %>% replace(. == 0, "0 FIXED"), collapse = " "), "\n",
                         paste(init[, "omega_gamma"] %>% replace(. == 0, "0 FIXED"), collapse = " "), "\n")
  ctl["sigma"] <- paste0("$SIGMA ", paste(init[, "sigma"], collapse = " "), "\n")
  ctl["estimation"] <-
    paste0("$ESTIMATION METHOD=1 MAXEVAL=99999 PRINT=1 NOABORT SIGDIGITS=2 FILE=iteration.csv NOTITLE=1\n")
  ctl["covariance"] <- "$COVARIANCE UNCONDITIONAL\n"
  ctl["table"] <-
    paste0("$TABLE ID TIME BM DV \nNOPRINT FORMAT=,F10.5 FILE=table.csv ONEHEADER NOTITLE")

  return(ctl)
}

makeSimulationNonmemData <- function (NUM_SUBJ, THETA, OMEGA, SIGMA, range_offsetT = c(-5, 20)) {
  simulation_data_spread <- makeSimulationDataSpread(NUM_SUBJ, THETA, OMEGA, SIGMA, range_offsetT)

  output <- convertSpreadDataToNonmemData(simulation_data_spread) %>%
    dplyr::select(ID, TIME, Biomarker, DV, everything()) %>%
    dplyr::arrange(ID, Biomarker, TIME)

  return(output)
}

convertSpreadDataToNonmemData <- function (df_input) {
  NUM_BM <- ncol(df_input) - 2
  output <- NULL
  df_summary <- df_input %>% distinct(ID)

  for (i in 1:NUM_BM) {
    temp <- df_input %>%
      dplyr::select(ID, TIME, DV = paste0("Biomarker", i)) %>%
      dplyr::mutate(Biomarker = i) %>%
      tidyr::drop_na()

    temp_summary <- temp %>%
      dplyr::group_by(ID) %>%
      dplyr::summarise(meanx = mean(TIME),
                       meany = mean(DV),
                       count = length(DV),
                       .groups = "drop") %>%
      setNames(c("ID", paste0(c("meanx", "meany", "count"), i)))

    output <- rbind(output, temp)
    df_summary <- merge(df_summary, temp_summary, by = "ID", all = TRUE)
  }

  output <- merge(output, df_summary, by = "ID") %>%
    select(ID, TIME, DV, Biomarker, paste0(rep(c("meanx", "meany", "count"), each = NUM_BM), 1:NUM_BM))

  return(output)
}


makeSimulationDataSpread <- function (NUM_SUBJ, THETA, OMEGA, SIGMA, TIME = 0:4, range_offsetT = c(-5, 20)) {
  NUM_BM <- length(THETA) / 3

  prms <- getParameterOfEachSubject(NUM_SUBJ, THETA, OMEGA, SIGMA, range_offsetT)

  output <- expand.grid(ID = 1:NUM_SUBJ,
                        Biomarker = 1:NUM_BM,
                        TIME = TIME) %>%
    merge(prms, by = "ID") %>%
    dplyr::mutate(TIME = TIME + offsetT) %>%
    dplyr::mutate(DV = apply(., 1, evalModelSimulation)) %>%
    dplyr::select(ID, TIME, DV, Biomarker) %>%
    tidyr::spread(key = Biomarker, value = DV) %>%
    setNames(c("ID", "TIME", paste0("Biomarker", 1:NUM_BM)))

  return(output)
}

getParameterOfEachSubject <- function (NUM_SUBJ, THETA, OMEGA, SIGMA, range_offsetT = c(-5, 20)) {
  NUM_BM <- length(THETA) / 3

  ETA <- MASS::mvrnorm(NUM_SUBJ, numeric(length(OMEGA)), diag(OMEGA))
  EPS <- MASS::mvrnorm(NUM_SUBJ, numeric(length(SIGMA)), diag(SIGMA))

  output <- (THETA + t(ETA)) %>%
    t() %>%
    cbind(EPS) %>%
    as.data.frame() %>%
    setNames(c(paste0(c("alpha", "beta", "gamma", "EPS") %>% rep(each = NUM_BM), 1:NUM_BM))) %>%
    dplyr::mutate(ID = 1:NUM_SUBJ,
                  offsetT = runif(nrow(.), min = range_offsetT[1], max = range_offsetT[2])) %>%
    dplyr::select(ID, offsetT, everything())

  return(output)
}

evalModelSimulation <- function (x) {
  prms_time  <- x["TIME"]
  prms_bm    <- x["Biomarker"]
  prms_alpha <- x[paste0("alpha", prms_bm)]
  prms_beta  <- x[paste0("beta", prms_bm)]
  prms_gamma <- x[paste0("gamma", prms_bm)]
  prms_eps <- x[paste0("EPS", prms_bm)]

  output <- prms_alpha + (prms_beta / prms_gamma) * (exp(prms_gamma * prms_time) - 1)
  output <- output + prms_eps

  return(output)
}


# User settings --------------------------------------------------------------------------------------


DIR_NAME <- "simulation"

NUM_SUBJ <- 3

THETA <- c(1.5, -2, 2,
           -0.20, 0.05, -0.10,
           -0.1, 0.1, 0.1)

OMEGA <- c(0, 0.1, 0.1,
           0.01, 0, 0,
           0, 0.01, 0.01)

SIGMA <- c(0.1, 0.1, 0.1)


# Parameter generation -------------------------------------------------------------------------------


df <- makeSimulationDataSpread(NUM_SUBJ, THETA, OMEGA, SIGMA)
prms <- getParameterOfEachSubject(NUM_SUBJ, THETA, OMEGA, SIGMA)

nmsheet <- makeSimulationNonmemData(NUM_SUBJ, THETA, OMEGA, SIGMA)

init <- estimatecInitialParameter(nmsheet, THETA[1])

ctl <- makeControlStreamSimulationTemplate(init, nmsheet)


# Observation generation -----------------------------------------------------------------------------



PROC_TIME <- getProcessTime()

ROOT <- "D:/Users/Ryota/OneDrive - 千葉大学/Project"
DIR_PATH <- paste0(ROOT, "notes/", PROC_TIME, "_", DIR_NAME, "/")
FILE_PATH <- paste0(DIR_PATH, PROC_TIME, "_")

dir.create(DIR_PATH)

fwrite(nmsheet, paste0(FILE_PATH, "data.csv"))
fwrite(init, paste0(FILE_PATH, "initialprms.csv"))
cat(paste(ctl, collapse = "\n"), file = paste0(FILE_PATH, "control.txt"))

DATA_TYPE <- "mean"
PRED_TYPE <- paste0("old_ana_", DATA_TYPE)
file.copy(from = paste0(ROOT, "src/pred/", PRED_TYPE, ".f90"), to = paste0(FILE_PATH, "pred_sreft.f90"))

