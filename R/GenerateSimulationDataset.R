# Set working directory ------------------------------------------------------------------------------


ROOT <- "D:/Users/Ryota/OneDrive - 千葉大学/Project/"
setwd(ROOT)
source("./src/r/00_Functions.R", encoding = "UTF-8")
source("./src/r/02_PreProcessFun.R", encoding = "UTF-8")
source("./src/r/03_PostProcessFun.R", encoding = "UTF-8")
library(MASS)

evalModelSimulation <- function (x) {
  .time  <- x["TIME"]
  .bm    <- x["BM"]
  .alpha <- x[paste0("alpha", .bm)]
  .beta  <- x[paste0("beta", .bm)]
  .gamma <- x[paste0("gamma", .bm)]
  .eps <- x[paste0("EPS", .bm)]
  
  output <- .alpha + (.beta / .gamma) * (exp(.gamma * .time) - 1)
  output <- output + .eps

  return(output)
}

makeNmsheetSimulation <- function (dataset) {
  dv <- dataset[, c("ID", "TIME", DV_NAME)]
  
  output <- switch(DATA_TYPE,
                   "mean" = makeNmsheetMeanSimulation(dv),
                   # "auc"  = makeNmsheetAUC(dv),
                   stop("Unexpected DATA_TYPE was specified.")) %>%
    mutate(serial = as.numeric(factor(.$ID))) %>%
    dplyr::select(ID, serial, TIME, DV, DVraw, BM, everything()) %>%
    arrange(ID, BM, TIME)
  
  return(output)
}

makeNmsheetMeanSimulation <- function (df) {
  for (i in 1:NUMBM) {
    tmp <- df %>% dplyr::select(ID, TIME, DV = DV_NAME[i], DVraw = DV_NAME[i]) %>% drop_na()
    
    tmp_x <- aggregate(TIME~ID, data = tmp, FUN = mean) %>% setNames(c("ID", paste0("meanx", i)))
    tmp_y <- aggregate(DV~ID, data = tmp, FUN = mean) %>% setNames(c("ID", paste0("meany", i)))
    tmp_c <- aggregate(DV~ID, data = tmp, FUN = length) %>% setNames(c("ID", paste0("count", i)))
    tmp$BM <- i
    if (i == 1) {
      meanx <- tmp_x
      meany <- tmp_y
      counts <- tmp_c
      output <- tmp
    } else {
      meanx <- merge(meanx, tmp_x, by = "ID", all = TRUE)
      meany <- merge(meany, tmp_y, by = "ID", all = TRUE)
      counts <- merge(counts, tmp_c, by = "ID", all = TRUE)
      output <- rbind(output, tmp)
    }
  }
  output <- list(output, meanx, meany, counts) %>% merge2(by = "ID") %>% signif(digits = 5)
  return(output)
}

setInitialPrms <- function (df) {
  output <- data.frame(Biomarker = DV_NAME)
  
  output <- calcInitialPrms(df, output)
  
  return(output)
}

calcInitialPrms <- function (df, output) {
  GTA <- 10
  
  .TIME_sp  <- split(df$TIME, list(df$ID, df$BM))
  .DV_sp    <- split(df$DV, list(df$ID, df$BM))
  
  .res <- data.frame(ID    = distinct(df, ID, BM) %>% arrange(BM, ID) %>% dplyr::select(ID),
                     BM    = distinct(df, ID, BM) %>% arrange(BM, ID) %>% dplyr::select(BM),
                     Mean  = mapply(mean, .TIME_sp),
                     Slope = mapply(slope, .TIME_sp, .DV_sp)) %>% drop_na()
  
  .Mean_sp  <- split(.res$Mean, .res$BM)
  .Slope_sp <- split(.res$Slope, .res$BM)
  
  .slopes <- mapply(slope, .Mean_sp, .Slope_sp)
  .intercepts <- mapply(intercept, .Mean_sp, .Slope_sp)
  .meanslope <- mapply(mean, .Slope_sp)
  
  output <- output %>%
    mutate(α = c(BASELINE, .meanslope[-1] / exp(.slopes[-1] * GTA)),
            β = .intercepts + α * .slopes,
            γ =  .slopes,
            omega_α = c(0, rep(0.0001, NUMBM - 1)),
            omega_β = c(0.0001, rep(0, NUMBM - 1)),
            omega_γ = c(0, rep(0.0001, NUMBM - 1)),
            sigma = 0.01)
  
  return(output)
}

makeControlStreamSimulation <- function (init, df) {
  ctl <- list()
  ctl["problem"] <- paste0("$PROBLEM simulation\n")
  ctl["input"] <- paste0("$INPUT ", paste(names(df), collapse = " "),"\n")
  ctl["data"] <- paste0("$DATA ", PROC_TIME, "_data.csv ignore=@\n")
  ctl["subroutine"] <- paste0("$SUBROUTINE PRED=./", PROC_TIME, "_pred_sreft.f90\n")
  ctl["theta"] <- paste0("$THETA\n" ,
                         paste(init[1, "α"], "FIXED", paste(init[2:NUMBM, "α"], collapse = " ")), "\n",
                         paste(init[, "β"], collapse = " "), "\n",
                         paste(init[, "γ"], collapse = " "), "\n")
  ctl["omega"] <- paste0("$OMEGA\n",
                         paste(init[, "omega_α"] %>% replace(. == 0, "0 FIXED"), collapse = " "), "\n",
                         paste(init[, "omega_β"] %>% replace(. == 0, "0 FIXED"), collapse = " "), "\n",
                         paste(init[, "omega_γ"] %>% replace(. == 0, "0 FIXED"), collapse = " "), "\n")
  ctl["sigma"] <- paste0("$SIGMA ", paste(init[, "sigma"], collapse = " "), "\n")
  ctl["estimation"] <- paste0("$ESTIMATION METHOD=1 MAXEVAL=99999 PRINT=1 NOABORT SIGDIGITS=2 FILE=",
                              PROC_TIME, "_iteration.csv NOTITLE=1\n")
  ctl["covariance"] <- "$COVARIANCE UNCONDITIONAL MATRIX=S\n"
  ctl["table"] <- paste0("$TABLE ID TIME DV DVraw BM PRED CIPRED CWRES\nNOPRINT FORMAT=,F10.5 FILE=",
                         PROC_TIME, "_table.csv ONEHEADER NOTITLE NOAPPEND")
  return(ctl)
}


# User settings --------------------------------------------------------------------------------------


DIR_NAME <- "simulation"

NUMSJ <- 3

THETA <- c(1.5, -2, 2,
           -0.20, 0.05, -0.10,
           -0.1, 0.1, 0.1)

OMEGA <- c(0, 0.1, 0.1,
           0.01, 0, 0,
           0, 0.01, 0.01)

SIGMA <- c(0.1, 0.1, 0.1)


# Parameter generation -------------------------------------------------------------------------------


NUMBM <- len(THETA) / 3
range_offsetT <- c(-5, 20)

ETA <- mvrnorm(NUMSJ, numeric(len(OMEGA)), diag(OMEGA))
EPS <- mvrnorm(NUMSJ, numeric(len(SIGMA)), diag(SIGMA))

prms <- (THETA + t(ETA)) %>%
  t() %>%
  cbind(EPS) %>%
  as.data.frame() %>%
  setNames(c(paste0(c("alpha", "beta", "gamma", "EPS") %>% rep(each = NUMBM), 1:NUMBM))) %>% 
  mutate(ID = 1:NUMSJ,
       offsetT = runif(nrow(.), min = range_offsetT[1], max = range_offsetT[2])) %>%
  dplyr::select(ID, offsetT, everything())


# Observation generation -----------------------------------------------------------------------------


DV_NAME <- paste0("BM", 1:NUMBM)
DATA_TYPE <- "mean"
PRED_TYPE <- paste0("old_ana_", DATA_TYPE)
BASELINE <- THETA[1]
PROC_TIME <- getProcessTime()

nmsheet <- expand.grid(ID = 1:NUMSJ,
                       BM = 1:NUMBM,
                       TIME = 0:4) %>%
  merge(prms, by = "ID") %>% 
  mutate(TIME = TIME + offsetT) %>% 
  mutate(DV = apply(., 1, evalModelSimulation)) %>% 
  dplyr::select(ID, TIME, DV, BM) %>%
  spread(key = BM, value = DV) %>%
  setNames(c("ID", "TIME", paste0("BM", 1:NUMBM))) %>%
  makeNmsheetSimulation()

init <- setInitialPrms(nmsheet)
ctl <- makeControlStreamSimulation(init, nmsheet)

DIR_PATH <- paste0(ROOT, "notes/", PROC_TIME, "_", DIR_NAME, "/")
FILE_PATH <- paste0(DIR_PATH, PROC_TIME, "_")

dir.create(DIR_PATH)

fwrite(nmsheet, paste0(FILE_PATH, "data.csv"))
fwrite(init, paste0(FILE_PATH, "initialprms.csv"))
cat(paste(ctl, collapse = "\n"), file = paste0(FILE_PATH, "control.txt"))
file.copy(from = paste0(ROOT, "src/pred/", PRED_TYPE, ".f90"), to = paste0(FILE_PATH, "pred_sreft.f90"))



# nmrun (optional) -----------------------------------------------------------------------------------


setwd(DIR_PATH)
system(paste("nmrun", PROC_TIME, "1"))
setwd(ROOT)

renameAfterNmrun(PROC_TIME, DIR_PATH)
