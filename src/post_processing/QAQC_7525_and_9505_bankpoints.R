## purpose: remove data above and below threshold percentiles and summarize stats by reach
# options(echo=TRUE)

library(foreign)
library(dplyr)
library(pastecs)

### summarize channel metrics from 1D calculation
### 1) filter metrics >95th/<5th and >75th/<25th
### 2) group metrics by LINKNO and calculate summary stats
### 3) export stats by LINKNO to dbf

# get args from python
myArgs <- commandArgs(trailingOnly = TRUE)
myArgs <- as.list(myArgs)
shp <- as.character(myArgs[1])

# read the dbf
raw_Metrics <- read.dbf(shp)

# if WBDFlag doesn't exist
if(!("WBDFlag" %in% colnames(raw_Metrics)))
{
  raw_Metrics["WBDFlag"] <- NA
}

raw_Metrics[is.na(raw_Metrics)] <- 0

# subset
Metrics <- subset(raw_Metrics, NHDFlag == 0 & WBDFlag == 0, linkno:NHDFlag)

# replace all metrics <= 0 w/NAs
Metrics$bankht_1d[Metrics$bankht_1d <= 0] = NA
Metrics$bank_elev[Metrics$bank_elev <= 0] = NA
Metrics$bnk_ang_1[Metrics$bnk_ang_1 <= 0] = NA
Metrics$bnk_ang_2[Metrics$bnk_ang_2 <= 0] = NA
Metrics$chan_area[Metrics$chan_area <= 0] = NA
Metrics$chan_width[Metrics$chan_width <= 0] = NA
Metrics$obank_rat[Metrics$obank_rat <= 0] = NA
Metrics$area_ratio[Metrics$area_ratio <= 0] = NA

# add new columns for additional filtering to be calculated, while keeping original fields
# 95th/5th
Metrics$BH955 <- Metrics$bankht_1d
Metrics$BEL955 <- Metrics$bank_elev
Metrics$BA1955 <- Metrics$bnk_ang_1
Metrics$BA2955 <- Metrics$bnk_ang_2
Metrics$BFA955 <- Metrics$chan_area
Metrics$CW955 <- Metrics$chan_width
Metrics$OBR955 <- Metrics$obank_rat
Metrics$AR955 <- Metrics$area_ratio

# 75th/25th
Metrics$BH7525 <- Metrics$bankht_1d
Metrics$BEL7525 <- Metrics$bank_elev
Metrics$BA17525 <- Metrics$bnk_ang_1
Metrics$BA27525 <- Metrics$bnk_ang_2
Metrics$BFA7525 <- Metrics$chan_area
Metrics$CW7525 <- Metrics$chan_width
Metrics$OBR7525 <- Metrics$obank_rat
Metrics$AR7525 <- Metrics$area_ratio


# group metrics by each reach to filter out outliers
linkno_list = unique(Metrics$linkno)
data_list = list()
for (i in seq_along(linkno_list)) {
    dat <- filter(Metrics, linkno == linkno_list[i])

    # replace all metrics >95th and <5th percentiles w/ NAs
    dat$BH955[dat$BH955 >= quantile(dat$BH955, c(0.95), na.rm = T) | dat$BH955 <= quantile(dat$BH955, c(0.05), na.rm = T)] = NA
    dat$BEL955[dat$BEL955 >= quantile(dat$BEL955, c(0.95), na.rm = T) | dat$BEL955 <= quantile(dat$BEL955, c(0.05), na.rm = T)] = NA
    dat$BA1955[dat$BA1955 >= quantile(dat$BA1955, c(0.95), na.rm = T) | dat$BA1955 <= quantile(dat$BA1955, c(0.05), na.rm = T)] = NA
    dat$BA2955[dat$BA2955 >= quantile(dat$BA2955, c(0.95), na.rm = T) | dat$BA2955 <= quantile(dat$BA2955, c(0.05), na.rm = T)] = NA
    dat$BFA955[dat$BFA955 >= quantile(dat$BFA955, c(0.95), na.rm = T) | dat$BFA955 <= quantile(dat$BFA955, c(0.05), na.rm = T)] = NA
    dat$CW955[dat$CW955 >= quantile(dat$CW955, c(0.95), na.rm = T) | dat$CW955 <= quantile(dat$CW955, c(0.05), na.rm = T)] = NA
    dat$OBR955[dat$OBR955 >= quantile(dat$OBR955, c(0.95), na.rm = T) | dat$OBR955 <= quantile(dat$OBR955, c(0.05), na.rm = T)] = NA
    dat$AR955[dat$AR955 >= quantile(dat$AR955, c(0.95), na.rm = T) | dat$AR955 <= quantile(dat$AR955, c(0.05), na.rm = T)] = NA

    # replace all dat >75th and <25th percentiles w/ NAs
    dat$BH7525[dat$BH7525 >= quantile(dat$BH7525, c(0.75), na.rm = T) | dat$BH7525 <= quantile(dat$BH7525, c(0.25), na.rm = T)] = NA
    dat$BEL7525[dat$BEL7525 >= quantile(dat$BEL7525, c(0.75), na.rm = T) | dat$BEL7525 <= quantile(dat$BEL7525, c(0.25), na.rm = T)] = NA
    dat$BA17525[dat$BA17525 >= quantile(dat$BA17525, c(0.75), na.rm = T) | dat$BA17525 <= quantile(dat$BA17525, c(0.25), na.rm = T)] = NA
    dat$BA27525[dat$BA27525 >= quantile(dat$BA27525, c(0.75), na.rm = T) | dat$BA27525 <= quantile(dat$BA27525, c(0.25), na.rm = T)] = NA
    dat$BFA7525[dat$BFA7525 >= quantile(dat$BFA7525, c(0.75), na.rm = T) | dat$BFA7525 <= quantile(dat$BFA7525, c(0.25), na.rm = T)] = NA
    dat$CW7525[dat$CW7525 >= quantile(dat$CW7525, c(0.75), na.rm = T) | dat$CW7525 <= quantile(dat$CW7525, c(0.25), na.rm = T)] = NA
    dat$OBR7525[dat$OBR7525 >= quantile(dat$OBR7525, c(0.75), na.rm = T) | dat$OBR7525 <= quantile(dat$OBR7525, c(0.25), na.rm = T)] = NA
    dat$AR7525[dat$AR7525 >= quantile(dat$AR7525, c(0.75), na.rm = T) | dat$AR7525 <= quantile(dat$AR7525, c(0.25), na.rm = T)] = NA

    data_list[[i]] <- dat
    }

tmp_metrics <- bind_rows(data_list)

#summarize each reach/linkno by mean, sd, median, IQR, and mean absolute deviation
Metrics_reach <- tmp_metrics %>%
    group_by(linkno) %>%
    summarise(mean(bankht_1d, na.rm = T), sd(bankht_1d, na.rm = T), median(bankht_1d, na.rm = T), IQR(bankht_1d, na.rm = T),
              median(BH955, na.rm = T), mean(BH955, na.rm = T), IQR(BH955, na.rm = T), mad(BH955, na.rm = T), sd(BH955, na.rm = T),
              median(BH7525, na.rm = T), mean(BH7525, na.rm = T), IQR(BH7525, na.rm = T), mad(BH7525, na.rm = T), sd(BH7525, na.rm = T),

              mean(bank_elev, na.rm = T), sd(bank_elev, na.rm = T), median(bank_elev, na.rm = T), IQR(bank_elev, na.rm = T),
              median(BEL955, na.rm = T), mean(BEL955, na.rm = T), IQR(BEL955, na.rm = T), mad(BEL955, na.rm = T), sd(BEL955, na.rm = T),
              median(BEL7525, na.rm = T), mean(BEL7525, na.rm = T), IQR(BEL7525, na.rm = T), mad(BEL7525, na.rm = T), sd(BEL7525, na.rm = T),

              mean(bnk_ang_1, na.rm = T), sd(bnk_ang_1, na.rm = T), median(bnk_ang_1, na.rm = T), IQR(bnk_ang_1, na.rm = T),
              median(BA1955, na.rm = T), mean(BA1955, na.rm = T), IQR(BA1955, na.rm = T), mad(BA1955, na.rm = T), sd(BA1955, na.rm = T),
              median(BA17525, na.rm = T), mean(BA17525, na.rm = T), IQR(BA17525, na.rm = T), mad(BA17525, na.rm = T), sd(BA17525, na.rm = T),

              mean(bnk_ang_2, na.rm = T), sd(bnk_ang_2, na.rm = T),median(bnk_ang_2, na.rm = T), IQR(bnk_ang_2, na.rm = T),
              median(BA2955, na.rm = T), mean(BA2955, na.rm = T), IQR(BA2955, na.rm = T), mad(BA2955, na.rm = T), sd(BA2955, na.rm = T),
              median(BA27525, na.rm = T), mean(BA27525, na.rm = T), IQR(BA27525, na.rm = T), mad(BA27525, na.rm = T), sd(BA27525, na.rm = T),

              mean(chan_area, na.rm = T), sd(chan_area, na.rm = T), median(chan_area, na.rm = T), IQR(chan_area, na.rm = T),
              median(BFA955, na.rm = T), mean(BFA955, na.rm = T), IQR(BFA955, na.rm = T), mad(BFA955, na.rm = T), sd(BFA955, na.rm = T),
              median(BFA7525, na.rm = T), mean(BFA7525, na.rm = T), IQR(BFA7525, na.rm = T), mad(BFA7525, na.rm = T), sd(BFA7525, na.rm = T),

              mean(chan_width, na.rm = T), sd(chan_width, na.rm = T), median(chan_width, na.rm = T), IQR(chan_width, na.rm = T),
              median(CW955, na.rm = T), mean(CW955, na.rm = T), IQR(CW955, na.rm = T), mad(CW955, na.rm = T), sd(CW955, na.rm = T),
              median(CW7525, na.rm = T), mean(CW7525, na.rm = T), IQR(CW7525, na.rm = T), mad(CW7525, na.rm = T), sd(CW7525, na.rm = T),

              mean(obank_rat, na.rm = T), sd(obank_rat, na.rm = T), median(obank_rat, na.rm = T), IQR(obank_rat, na.rm = T),
              median(OBR955, na.rm = T), mean(OBR955, na.rm = T), IQR(OBR955, na.rm = T), mad(OBR955, na.rm = T), sd(OBR955, na.rm = T),
              median(OBR7525, na.rm = T), mean(OBR7525, na.rm = T), IQR(OBR7525, na.rm = T), mad(OBR7525, na.rm = T), sd(OBR7525, na.rm = T),

              mean(area_ratio, na.rm = T), sd(area_ratio, na.rm = T), median(area_ratio, na.rm = T), IQR(area_ratio, na.rm = T),
              median(AR955, na.rm = T), mean(AR955, na.rm = T), IQR(AR955, na.rm = T), mad(AR955, na.rm = T), sd(AR955, na.rm = T),
              median(AR7525, na.rm = T), mean(AR7525, na.rm = T), IQR(AR7525, na.rm = T), mad(AR7525, na.rm = T), sd(AR7525, na.rm = T))
# rename column headings
col_headings <- c("LINKNO",
                "BHmean", "BHsd", "BHmed", "BHIQR", "BH955med", "BH955mean", "BH955IQR", "BH955mad", "BH955sd", "BH7525med", "BH7525mean", "BH7525IQR", "BH7525mad", "BH7525sd",
                "BELmean", "BELsd", "BELmed", "BELIQR", "BEL955med", "BEL955mean", "BEL955IQR", "BEL955mad", "BEL955sd", "BEL7525med", "BEL7525mean", "BEL7525IQR", "BEL7525mad", "BEL7525sd",
                "BA1mean", "BA1sd", "BA1med", "BA1IQR", "BA1955med", "BA1955mean", "BA1955IQR", "BA1955mad", "BA1955sd", "BA17525med", "BA17525mean", "BA17525IQR", "BA17525mad", "BA17525sd",
                "BA2mean", "BA2sd", "BA2med", "BA2IQR", "BA2955med", "BA2955mean", "BA2955IQR", "BA2955mad", "BA2955sd", "BA27525med", "BA27525mean", "BA27525IQR", "BA27525mad", "BA27525sd",
                "CNAmean", "CNAsd", "CNAmed", "CNAIQR", "CNA955med", "CNA955mean", "CNA955IQR", "CNA955mad", "CNA955sd", "CNA7525med", "CNA7525mean", "CNA7525IQR", "CNA7525mad", "CNA7525sd",
                "CWmean", "CWsd", "CWmed", "CWIQR", "CW955med", "CW955mean", "CW955IQR", "CW955mad", "CW955sd", "CW7525med", "CW7525mean", "CW7525IQR", "CW7525mad", "CW7525sd",
                "OBRmean", "OBRsd", "OBRmed", "OBRIQR", "OBR955med", "OBR955mean", "OBR955IQR", "OBR955mad", "OBR955sd", "OBR7525med", "OBR7525mean", "OBR7525IQR", "OBR7525mad", "OBR7525sd",
                "ARmean", "ARsd", "ARmed", "ARIQR", "AR955med", "AR955mean", "AR955IQR", "AR955mad", "AR955sd", "AR7525med", "AR7525mean", "AR7525IQR", "AR7525mad", "AR7525sd")

names(Metrics_reach) <- col_headings

Metrics_reach$BHCOV <- Metrics_reach$BHsd/Metrics_reach$BHmean*100
Metrics_reach$BH955COV <- Metrics_reach$BH955sd/Metrics_reach$BH955mean*100
Metrics_reach$BH7525COV <- Metrics_reach$BH7525sd/Metrics_reach$BH7525mean*100

Metrics_reach$BELCOV <- Metrics_reach$BELsd/Metrics_reach$BELmean*100
Metrics_reach$BEL955COV <- Metrics_reach$BEL955sd/Metrics_reach$BEL955mean*100
Metrics_reach$BEL7525COV <- Metrics_reach$BEL7525sd/Metrics_reach$BEL7525mean*100

Metrics_reach$BA1COV <- Metrics_reach$BA1sd/Metrics_reach$BA1mean*100
Metrics_reach$BA1955COV <- Metrics_reach$BA1955sd/Metrics_reach$BA1955mean*100
Metrics_reach$BA17525COV <- Metrics_reach$BA17525sd/Metrics_reach$BA17525mean*100

Metrics_reach$BA2COV <- Metrics_reach$BA2sd/Metrics_reach$BA2mean*100
Metrics_reach$BA2955COV <- Metrics_reach$BA2955sd/Metrics_reach$BA2955mean*100
Metrics_reach$BA27525COV <- Metrics_reach$BA27525sd/Metrics_reach$BA27525mean*100

Metrics_reach$CNACOV <- Metrics_reach$CNAsd/Metrics_reach$CNAmean*100
Metrics_reach$CNA955COV <- Metrics_reach$CNA955sd/Metrics_reach$CNA955mean*100
Metrics_reach$CNA7525COV <- Metrics_reach$CNA7525sd/Metrics_reach$CNA7525mean*100

Metrics_reach$CWCOV <- Metrics_reach$CWsd/Metrics_reach$CWmean*100
Metrics_reach$CW955COV <- Metrics_reach$CW955sd/Metrics_reach$CW955mean*100
Metrics_reach$CW7525COV <- Metrics_reach$CW7525sd/Metrics_reach$CW7525mean*100

Metrics_reach$OBRCOV <- Metrics_reach$OBRsd/Metrics_reach$OBRmean*100
Metrics_reach$OBR955COV <- Metrics_reach$OBR955sd/Metrics_reach$OBR955mean*100
Metrics_reach$OBR7525COV <- Metrics_reach$OBR7525sd/Metrics_reach$OBR7525mean*100

Metrics_reach$ARCOV <- Metrics_reach$ARsd/Metrics_reach$ARmean*100
Metrics_reach$AR955COV <- Metrics_reach$AR955sd/Metrics_reach$AR955mean*100
Metrics_reach$AR7525COV <- Metrics_reach$AR7525sd/Metrics_reach$AR7525mean*100

#export back to .dbf
Metrics_reach <- as.data.frame(Metrics_reach)

# write out output
file_name = basename(tools::file_path_sans_ext(shp)) # e.g. breach_bankpts_flagged
output_file = file.path(dirname(shp), paste0(file_name, '_all_stats.dbf'))
write.dbf(Metrics_reach, output_file)