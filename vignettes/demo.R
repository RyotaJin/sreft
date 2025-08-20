library(sreft)
library(dplyr)
library(ggplot2)

# setwd("")
# dir.create("demo")

set.seed(42)
df <- makeDemodata(500,
                   c(-2, 2, 1.5),
                   c(0.05, -0.05, -0.2),
                   c(0.1, 0.1, -0.1),
                   sqrt(c(0.125, 0.18, 0)),
                   sqrt(c(0, 0, 0.001)),
                   sqrt(c(0.0013, 0.0005, 0)),
                   c(0.2, 0.2, 0.2),
                   -5,
                   20,
                   0:4,
                   FALSE)

ggplot(df, aes(x = TIME, y = DV, group = ID, color = factor(BM))) +
  geom_point() +
  geom_line() +
  facet_wrap(.~BM, ncol = 3, scales = "free") +
  guides(color = "none")

ggplot(df, aes(x = TIME2, y = DV, group = ID, color = factor(BM))) +
  geom_point() +
  geom_line() +
  facet_wrap(.~BM, ncol = 3, scales = "free") +
  guides(color = "none")


set.seed(42)
df_ <- makeDemodata(500,
                    c(-2, 2, 1.5),
                    c(0.05, -0.05, -0.2),
                    c(0.1, 0.1, -0.1),
                    sqrt(c(0.125, 0.18, 0)),
                    sqrt(c(0, 0, 0.001)),
                    sqrt(c(0.0013, 0.0005, 0)),
                    c(0.2, 0.2, 0.2),
                    -5,
                    20,
                    0:4)

nmsheet <- makeNonmemSheet(df_, paste0("Biomarker", 1:3), "offsetT", lognorm = FALSE) %>%
  select(-CMT_name)

write.csv(nmsheet, "data.csv", row.names = FALSE)


inits <- setInitialPrms(df_, paste0("Biomarker", 1:3), "Biomarker3", 1.5, estimated_mean_offsetT = 10, lognorm = FALSE)
runno <- "run001"
ctl <- makeControlStream(inits,
                         nmsheet,
                         runno,
                         3)
cat(ctl, file = paste0(runno, ".mod"))

f90 <- makef90(nmsheet, runno)
cat(f90, file = paste0(runno, ".f90"))


# system("execute run001.mod)


df_offsetT_pred <- read.csv("offsetT.csv") %>%
  rename(offsetT_pred = offsetT)

df_fit <- read.csv(paste0(runno, ".fit")) %>%
  left_join(df_offsetT_pred) %>%
  mutate(TIME = TIME + offsetT_pred)

# Load Final estimate
prms <- read.table(paste0(runno, ".ext"), skip = 1) %>%
  filter(ITERATION == -1000000000) %>%
  select(matches("THETA")) %>%
  setNames(paste0(rep(c("alpha", "beta", "gamma"), each = (ncol(.) / 3)), 1:(ncol(.) / 3)))

# prediction
df_pred <- expand.grid(TIME = seq(min(df_fit$TIME), max(df_fit$TIME), len = 100),
                       CMT = 1:max(df_fit$CMT),
                       KEEP.OUT.ATTRS = FALSE) %>%
  mutate(pred = apply(., 1, evalModel, prms))

ggplot(df_fit, aes(x = TIME, y = DV)) +
  geom_point(size = 0.7, alpha = 0.7, shape = 16) +
  geom_line(linewidth = 0.3, alpha = 0.7, aes(group = ID)) +
  geom_line(data = df_pred, aes(x = TIME, y = pred, group = CMT), size = 1.5, colour = "#3366FF") +
  facet_wrap(~CMT, scales = "free_y") +
  labs(x = "Time (year)", y = "DV")

# offsetT
df_fit %>%
  distinct(ID, offsetT, offsetT_pred) %>%
  ggplot(aes(x = offsetT, y = offsetT_pred)) +
  geom_abline(slope = 1) +
  geom_point() +
  stat_smooth()
