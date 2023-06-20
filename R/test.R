# df <- read.csv("C:/Users/ryota/OneDrive - 千葉大学/Project_ml/data/SOLVD/data_solvd_bbc0cea45017694b6792961a26ade527730acb92.csv") %>%
#   dplyr::rename(TIME = vTIME_y)
# nmsheet <- makeNonmemSheet(df = df,
#                            selected_bm = c("WBC", "BUN", "Na", "SBP", "DBP"),
#                            selected_cov = c("ALLD"),
#                            lognorm = TRUE)
# int_prms <- setInitialPrms(df = df,
#                selected_bm = c("WBC", "BUN", "Na", "SBP", "DBP"),
#                definition_bm = "Na",
#                definition_value = 130)
