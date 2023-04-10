## purpose: remove data above and below threshold percentiles and summarize stats by reach
# options(echo=TRUE)

library(foreign)
library(dplyr)
library(pastecs)

### summarize floodplain metrics from 1D calculation
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
Metrics$fpwid_1d[Metrics$fpwid_1d <= 0] = NA
Metrics$mind_1d[Metrics$mind_1d <= 0] = NA
Metrics$maxd_1d[Metrics$maxd_1d <= 0] = NA
Metrics$rngd_1d[Metrics$rngd_1d <= 0] = NA
Metrics$meand_1d[Metrics$meand_1d <= 0] = NA
Metrics$stdd_1d[Metrics$stdd_1d <= 0] = NA
Metrics$sumd_1d[Metrics$sumd_1d <= 0] = NA
Metrics$mine_1d[Metrics$mine_1d <= 0] = NA
Metrics$maxe_1d[Metrics$maxe_1d <= 0] = NA
Metrics$rnge_1d[Metrics$rnge_1d <= 0] = NA
Metrics$meane_1d[Metrics$meane_1d <= 0] = NA
Metrics$stde_1d[Metrics$stde_1d <= 0] = NA
Metrics$sume_1d[Metrics$sume_1d <= 0] = NA

# add new columns for additional filtering to be calculated, while keeping original fields too
# 95th/5th
Metrics$TW955 <- Metrics$fpwid_1d
Metrics$MIND955 <- Metrics$mind_1d
Metrics$MAXD955 <- Metrics$maxd_1d
Metrics$RNGD955 <- Metrics$rngd_1d
Metrics$MEAND955 <- Metrics$meand_1d
Metrics$STDD955 <- Metrics$stdd_1d
Metrics$SUMD955 <- Metrics$sumd_1d
Metrics$MINE955 <- Metrics$mine_1d
Metrics$MAXE955 <- Metrics$maxe_1d
Metrics$RNGE955 <- Metrics$rnge_1d
Metrics$MEANE955 <- Metrics$meane_1d
Metrics$STDE955 <- Metrics$stde_1d
Metrics$SUME955 <- Metrics$sume_1d

# 75th/25th
Metrics$TW7525 <- Metrics$fpwid_1d
Metrics$MIND7525 <- Metrics$mind_1d
Metrics$MAXD7525 <- Metrics$maxd_1d
Metrics$RNGD7525 <- Metrics$rngd_1d
Metrics$MEAND7525 <- Metrics$meand_1d
Metrics$STDD7525 <- Metrics$stdd_1d
Metrics$SUMD7525 <- Metrics$sumd_1d
Metrics$MINE7525 <- Metrics$mine_1d
Metrics$MAXE7525 <- Metrics$maxe_1d
Metrics$RNGE7525 <- Metrics$rnge_1d
Metrics$MEANE7525 <- Metrics$meane_1d
Metrics$STDE7525 <- Metrics$stde_1d
Metrics$SUME7525 <- Metrics$sume_1d

# group metrics by each reach to filter out outliers
linkno_list = unique(Metrics$linkno)
data_list = list()
for (i in seq_along(linkno_list)) {
    dat <- filter(Metrics, linkno == linkno_list[i])

    # replace all metrics >95th and <5th percentiles w/ NAs
    dat$TW955[dat$TW955 >= quantile(dat$TW955, c(0.95), na.rm = T) | dat$TW955 <= quantile(dat$TW955, c(0.05), na.rm = T)] = NA
    dat$MIND955[dat$MIND955 >= quantile(dat$MIND955, c(0.95), na.rm = T) | dat$MIND955 <= quantile(dat$MIND955, c(0.05), na.rm = T)] = NA
    dat$MAXD955[dat$MAXD955 >= quantile(dat$MAXD955, c(0.95), na.rm = T) | dat$MAXD955 <= quantile(dat$MAXD955, c(0.05), na.rm = T)] = NA
    dat$RNGD955[dat$RNGD955 >= quantile(dat$RNGD955, c(0.95), na.rm = T) | dat$RNGD955 <= quantile(dat$RNGD955, c(0.05), na.rm = T)] = NA
    dat$MEAND955[dat$MEAND955 >= quantile(dat$MEAND955, c(0.95), na.rm = T) | dat$MEAND955 <= quantile(dat$MEAND955, c(0.05), na.rm = T)] = NA
    dat$STDD955[dat$STDD955 >= quantile(dat$STDD955, c(0.95), na.rm = T) | dat$STDD955 <= quantile(dat$STDD955, c(0.05), na.rm = T)] = NA
    dat$SUMD955[dat$SUMD955 >= quantile(dat$SUMD955, c(0.95), na.rm = T) | dat$SUMD955 <= quantile(dat$SUMD955, c(0.05), na.rm = T)] = NA
    dat$MINE955[dat$MINE955 >= quantile(dat$MINE955, c(0.95), na.rm = T) | dat$MINE955 <= quantile(dat$MINE955, c(0.05), na.rm = T)] = NA
    dat$MAXE955[dat$MAXE955 >= quantile(dat$MAXE955, c(0.95), na.rm = T) | dat$MAXE955 <= quantile(dat$MAXE955, c(0.05), na.rm = T)] = NA
    dat$RNGE955[dat$RNGE955 >= quantile(dat$RNGE955, c(0.95), na.rm = T) | dat$RNGE955 <= quantile(dat$RNGE955, c(0.05), na.rm = T)] = NA
    dat$MEANE955[dat$MEANE955 >= quantile(dat$MEANE955, c(0.95), na.rm = T) | dat$MEANE955 <= quantile(dat$MEANE955, c(0.05), na.rm = T)] = NA
    dat$STDE955[dat$STDE955 >= quantile(dat$STDE955, c(0.95), na.rm = T) | dat$STDE955 <= quantile(dat$STDE955, c(0.05), na.rm = T)] = NA
    dat$SUME955[dat$SUME955 >= quantile(dat$SUME955, c(0.95), na.rm = T) | dat$SUME955 <= quantile(dat$SUME955, c(0.05), na.rm = T)] = NA

    # replace all metrics >75th and <25th percentiles w/ NAs
    dat$TW7525[dat$TW7525 >= quantile(dat$TW7525, c(0.75), na.rm = T) | dat$TW7525 <= quantile(dat$TW7525, c(0.25), na.rm = T)] = NA
    dat$MIND7525[dat$MIND7525 >= quantile(dat$MIND7525, c(0.75), na.rm = T) | dat$MIND7525 <= quantile(dat$MIND7525, c(0.25), na.rm = T)] = NA
    dat$MAXD7525[dat$MAXD7525 >= quantile(dat$MAXD7525, c(0.75), na.rm = T) | dat$MAXD7525 <= quantile(dat$MAXD7525, c(0.25), na.rm = T)] = NA
    dat$RNGD7525[dat$RNGD7525 >= quantile(dat$RNGD7525, c(0.75), na.rm = T) | dat$RNGD7525 <= quantile(dat$RNGD7525, c(0.25), na.rm = T)] = NA
    dat$MEAND7525[dat$MEAND7525 >= quantile(dat$MEAND7525, c(0.75), na.rm = T) | dat$MEAND7525 <= quantile(dat$MEAND7525, c(0.25), na.rm = T)] = NA
    dat$STDD7525[dat$STDD7525 >= quantile(dat$STDD7525, c(0.75), na.rm = T) | dat$STDD7525 <= quantile(dat$STDD7525, c(0.25), na.rm = T)] = NA
    dat$SUMD7525[dat$SUMD7525 >= quantile(dat$SUMD7525, c(0.75), na.rm = T) | dat$SUMD7525 <= quantile(dat$SUMD7525, c(0.25), na.rm = T)] = NA
    dat$MINE7525[dat$MINE7525 >= quantile(dat$MINE7525, c(0.75), na.rm = T) | dat$MINE7525 <= quantile(dat$MINE7525, c(0.25), na.rm = T)] = NA
    dat$MAXE7525[dat$MAXE7525 >= quantile(dat$MAXE7525, c(0.75), na.rm = T) | dat$MAXE7525 <= quantile(dat$MAXE7525, c(0.25), na.rm = T)] = NA
    dat$RNGE7525[dat$RNGE7525 >= quantile(dat$RNGE7525, c(0.75), na.rm = T) | dat$RNGE7525 <= quantile(dat$RNGE7525, c(0.25), na.rm = T)] = NA
    dat$MEANE7525[dat$MEANE7525 >= quantile(dat$MEANE7525, c(0.75), na.rm = T) | dat$MEANE7525 <= quantile(dat$MEANE7525, c(0.25), na.rm = T)] = NA
    dat$STDE7525[dat$STDE7525 >= quantile(dat$STDE7525, c(0.75), na.rm = T) | dat$STDE7525 <= quantile(dat$STDE7525, c(0.25), na.rm = T)] = NA
    dat$SUME7525[dat$SUME7525 >= quantile(dat$SUME7525, c(0.75), na.rm = T) | dat$SUME7525 <= quantile(dat$SUME7525, c(0.25), na.rm = T)] = NA

    data_list[[i]] <- dat
    }

tmp_metrics <- bind_rows(data_list)

#summarize each reach/linkno by mean, sd, median, IQR, and mean absolute deviation
Metrics_reach <- tmp_metrics %>%
    group_by(linkno) %>%
    summarise(mean(fpwid_1d, na.rm = T), sd(fpwid_1d, na.rm = T), median(fpwid_1d, na.rm = T), IQR(fpwid_1d, na.rm = T),
            median(TW955, na.rm = T), mean(TW955, na.rm = T), IQR(TW955, na.rm = T), mad(TW955, na.rm = T), sd(TW955, na.rm = T),
            median(TW7525, na.rm = T), mean(TW7525, na.rm = T), IQR(TW7525, na.rm = T), mad(TW7525, na.rm = T), sd(TW7525, na.rm = T),

            mean(mind_1d, na.rm = T), sd(mind_1d, na.rm = T), median(mind_1d, na.rm = T), IQR(mind_1d, na.rm = T),
            median(MIND955, na.rm = T), mean(MIND955, na.rm = T), IQR(MIND955, na.rm = T), mad(MIND955, na.rm = T), sd(MIND955, na.rm = T),
            median(MIND7525, na.rm = T), mean(MIND7525, na.rm = T), IQR(MIND7525, na.rm = T), mad(MIND7525, na.rm = T), sd(MIND7525, na.rm = T),

            mean(maxd_1d, na.rm = T), sd(maxd_1d, na.rm = T), median(maxd_1d, na.rm = T), IQR(maxd_1d, na.rm = T),
            median(MAXD955, na.rm = T), mean(MAXD955, na.rm = T), IQR(MAXD955, na.rm = T), mad(MAXD955, na.rm = T), sd(MAXD955, na.rm = T),
            median(MAXD7525, na.rm = T), mean(MAXD7525, na.rm = T), IQR(MAXD7525, na.rm = T), mad(MAXD7525, na.rm = T), sd(MAXD7525, na.rm = T),

            mean(rngd_1d, na.rm = T), sd(rngd_1d, na.rm = T),median(rngd_1d, na.rm = T), IQR(rngd_1d, na.rm = T),
            median(RNGD955, na.rm = T), mean(RNGD955, na.rm = T), IQR(RNGD955, na.rm = T), mad(RNGD955, na.rm = T), sd(RNGD955, na.rm = T),
            median(RNGD7525, na.rm = T), mean(RNGD7525, na.rm = T), IQR(RNGD7525, na.rm = T), mad(RNGD7525, na.rm = T), sd(RNGD7525, na.rm = T),

            mean(meand_1d, na.rm = T), sd(meand_1d, na.rm = T), median(meand_1d, na.rm = T), IQR(meand_1d, na.rm = T),
            median(MEAND955, na.rm = T), mean(MEAND955, na.rm = T), IQR(MEAND955, na.rm = T), mad(MEAND955, na.rm = T), sd(MEAND955, na.rm = T),
            median(MEAND7525, na.rm = T), mean(MEAND7525, na.rm = T), IQR(MEAND7525, na.rm = T), mad(MEAND7525, na.rm = T), sd(MEAND7525, na.rm = T),

            mean(stdd_1d, na.rm = T), sd(stdd_1d, na.rm = T), median(stdd_1d, na.rm = T), IQR(stdd_1d, na.rm = T),
            median(STDD955, na.rm = T), mean(STDD955, na.rm = T), IQR(STDD955, na.rm = T), mad(STDD955, na.rm = T), sd(STDD955, na.rm = T),
            median(STDD7525, na.rm = T), mean(STDD7525, na.rm = T), IQR(STDD7525, na.rm = T), mad(STDD7525, na.rm = T), sd(STDD7525, na.rm = T),

            mean(sumd_1d, na.rm = T), sd(sumd_1d, na.rm = T), median(sumd_1d, na.rm = T), IQR(sumd_1d, na.rm = T),
            median(SUMD955, na.rm = T), mean(SUMD955, na.rm = T), IQR(SUMD955, na.rm = T), mad(SUMD955, na.rm = T), sd(SUMD955, na.rm = T),
            median(SUMD7525, na.rm = T), mean(SUMD7525, na.rm = T), IQR(SUMD7525, na.rm = T), mad(SUMD7525, na.rm = T), sd(SUMD7525, na.rm = T),

            mean(mine_1d, na.rm = T), sd(mine_1d, na.rm = T), median(mine_1d, na.rm = T), IQR(mine_1d, na.rm = T),
            median(MINE955, na.rm = T), mean(MINE955, na.rm = T), IQR(MINE955, na.rm = T), mad(MINE955, na.rm = T), sd(MINE955, na.rm = T),
            median(MINE7525, na.rm = T), mean(MINE7525, na.rm = T), IQR(MINE7525, na.rm = T), mad(MINE7525, na.rm = T), sd(MINE7525, na.rm = T),

            mean(maxe_1d, na.rm = T), sd(maxe_1d, na.rm = T), median(maxe_1d, na.rm = T), IQR(maxe_1d, na.rm = T),
            median(MAXE955, na.rm = T), mean(MAXE955, na.rm = T), IQR(MAXE955, na.rm = T), mad(MAXE955, na.rm = T), sd(MAXE955, na.rm = T),
            median(MAXE7525, na.rm = T), mean(MAXE7525, na.rm = T), IQR(MAXE7525, na.rm = T), mad(MAXE7525, na.rm = T), sd(MAXE7525, na.rm = T),

            mean(rnge_1d, na.rm = T), sd(rnge_1d, na.rm = T), median(rnge_1d, na.rm = T), IQR(rnge_1d, na.rm = T),
            median(RNGE955, na.rm = T), mean(RNGE955, na.rm = T), IQR(RNGE955, na.rm = T), mad(RNGE955, na.rm = T), sd(RNGE955, na.rm = T),
            median(RNGE7525, na.rm = T), mean(RNGE7525, na.rm = T), IQR(RNGE7525, na.rm = T), mad(RNGE7525, na.rm = T), sd(RNGE7525, na.rm = T),

            mean(meane_1d, na.rm = T), sd(meane_1d, na.rm = T), median(meane_1d, na.rm = T), IQR(meane_1d, na.rm = T),
            median(MEANE955, na.rm = T), mean(MEANE955, na.rm = T), IQR(MEANE955, na.rm = T), mad(MEANE955, na.rm = T), sd(MEANE955, na.rm = T),
            median(MEANE7525, na.rm = T), mean(MEANE7525, na.rm = T), IQR(MEANE7525, na.rm = T), mad(MEANE7525, na.rm = T), sd(MEANE7525, na.rm = T),

            mean(stde_1d, na.rm = T), sd(stde_1d, na.rm = T), median(stde_1d, na.rm = T), IQR(stde_1d, na.rm = T),
            median(STDE955, na.rm = T), mean(STDE955, na.rm = T), IQR(STDE955, na.rm = T), mad(STDE955, na.rm = T), sd(STDE955, na.rm = T),
            median(STDE7525, na.rm = T), mean(STDE7525, na.rm = T), IQR(STDE7525, na.rm = T), mad(STDE7525, na.rm = T), sd(STDE7525, na.rm = T),

            mean(sume_1d, na.rm = T), sd(sume_1d, na.rm = T), median(sume_1d, na.rm = T), IQR(sume_1d, na.rm = T),
            median(SUME955, na.rm = T), mean(SUME955, na.rm = T), IQR(SUME955, na.rm = T), mad(SUME955, na.rm = T), sd(SUME955, na.rm = T),
            median(SUME7525, na.rm = T), mean(SUME7525, na.rm = T), IQR(SUME7525, na.rm = T), mad(SUME7525, na.rm = T), sd(SUME7525, na.rm = T))

  # rename column headings
  col_headings <- c("LINKNO",
                    "FWmean", "FWsd", "FWmed", "FWIQR", "FW955med", "FW955mean", "FW955IQR", "FW955mad", "FW955sd", "FW7525med", "FW7525mean", "FW7525IQR", "FW7525mad", "FW7525sd",
                    "MINDmean", "MINDsd", "MINDmed", "MINDIQR", "MIND955med", "MIND955mean", "MIND955IQR", "MIND955mad", "MIND955sd", "MIND75med", "MIND75mean", "MIND75IQR", "MIND75mad", "MIND75sd",
                    "MAXDmean", "MAXDsd", "MAXDmed", "MAXDIQR", "MAXD955med", "MAXD955mean", "MAXD955IQR", "MAXD955mad", "MAXD955sd", "MAXD75med", "MAXD75mean", "MAXD75IQR", "MAXD75mad", "MAXD75sd",
                    "RNGDmean", "RNGDsd", "RNGDmed", "RNGDIQR", "RNGD955med", "RNGD955mean", "RNGD955IQR", "RNGD955mad", "RNGD955sd", "RNGD75med", "RNGD75mean", "RNGD75IQR", "RNGD75mad", "RNGD75sd",
                    "MEANDmean", "MEANDsd", "MEANDmed", "MEANDIQR", "MEAND95med", "MEAND95mean", "MEAND95IQR", "MEAND95mad", "MEAND95sd", "MEAND75med", "MEAND75mean", "MEAND75IQR", "MEAND75mad", "MEAND75sd",
                    "STDDmean", "STDDsd", "STDDmed", "STDDIQR", "STDD955med", "STDD955mean", "STDD955IQR", "STDD955mad", "STDD955sd", "STDD75med", "STDD75mean", "STDD75IQR", "STDD75mad", "STDD75sd",
                    "SUMDmean", "SUMDsd", "SUMDmed", "SUMDIQR", "SUMD955med", "SUMD955mean", "SUMD955IQR", "SUMD955mad", "SUMD955sd", "SUMD75med", "SUMD75mean", "SUMD75IQR", "SUMD75mad", "SUMD75sd",
                    "MINEmean", "MINEsd", "MINEmed", "MINEIQR", "MINE955med", "MINE955mean", "MINE955IQR", "MINE955mad", "MINE955sd", "MINE75med", "MINE75mean", "MINE75IQR", "MINE75mad", "MINE75sd",
                    "MAXEmean", "MAXEsd", "MAXEmed", "MAXEIQR", "MAXE955med", "MAXE955mean", "MAXE955IQR", "MAXE955mad", "MAXE955sd", "MAXE75med", "MAXE75mean", "MAXE75IQR", "MAXE75mad", "MAXE75sd",
                    "RNGEmean", "RNGEsd", "RNGEmed", "RNGEIQR", "RNGE955med", "RNGE955mean", "RNGE955IQR", "RNGE955mad", "RNGE955sd", "RNGE75med", "RNGE75mean", "RNGE75IQR", "RNGE75mad", "RNGE75sd",
                    "MEANEmean", "MEANEsd", "MEANEmed", "MEANEIQR", "MEANE95med", "MEANE95mean", "MEANE95IQR", "MEANE95mad", "MEANE95sd", "MEANE75med", "MEANE75mean", "MEANE75IQR", "MEANE75mad", "MEANE75sd",
                    "STDEmean", "STDEsd", "STDEmed", "STDEIQR", "STDE955med", "STDE955mean", "STDE955IQR", "STDE955mad", "STDE955sd", "STDE75med", "STDE75mean", "STDE75IQR", "STDE75mad", "STDE75sd",
                    "SUMEmean", "SUMEsd", "SUMEmed", "SUMEIQR", "SUME955med", "SUME955mean", "SUME955IQR", "SUME955mad", "SUME955sd", "SUME75med", "SUME75mean", "SUME75IQR", "SUME75mad", "SUME75sd")

names(Metrics_reach) <- col_headings

Metrics_reach$FWCOV <- Metrics_reach$FWsd/Metrics_reach$FWmean*100
Metrics_reach$FW955COV <- Metrics_reach$FW955sd/Metrics_reach$FW955mean*100
Metrics_reach$FW7525COV <- Metrics_reach$FW7525sd/Metrics_reach$FW7525mean*100

Metrics_reach$MINDCOV <- Metrics_reach$MINDsd/Metrics_reach$MINDmean*100
Metrics_reach$MIND955COV <- Metrics_reach$MIND955sd/Metrics_reach$MIND955mean*100
Metrics_reach$MIND75COV <- Metrics_reach$MIND75sd/Metrics_reach$MIND75mean*100

Metrics_reach$MAXDCOV <- Metrics_reach$MAXDsd/Metrics_reach$MAXDmean*100
Metrics_reach$MAXD955COV <- Metrics_reach$MAXD955sd/Metrics_reach$MAXD955mean*100
Metrics_reach$MAXD75COV <- Metrics_reach$MAXD75sd/Metrics_reach$MAXD75mean*100

Metrics_reach$RNGDCOV <- Metrics_reach$RNGDsd/Metrics_reach$RNGDmean*100
Metrics_reach$RNGD955COV <- Metrics_reach$RNGD955sd/Metrics_reach$RNGD955mean*100
Metrics_reach$RNGD75COV <- Metrics_reach$RNGD75sd/Metrics_reach$RNGD75mean*100

Metrics_reach$MEANDCOV <- Metrics_reach$MEANDsd/Metrics_reach$MEANDmean*100
Metrics_reach$MEAND95COV <- Metrics_reach$MEAND95sd/Metrics_reach$MEAND95mean*100
Metrics_reach$MEAND75COV <- Metrics_reach$MEAND75sd/Metrics_reach$MEAND75mean*100

Metrics_reach$STDDCOV <- Metrics_reach$STDDsd/Metrics_reach$STDDmean*100
Metrics_reach$STDD955COV <- Metrics_reach$STDD955sd/Metrics_reach$STDD955mean*100
Metrics_reach$STDD75COV <- Metrics_reach$STDD75sd/Metrics_reach$STDD75mean*100

Metrics_reach$SUMDCOV <- Metrics_reach$SUMDsd/Metrics_reach$SUMDmean*100
Metrics_reach$SUMD955COV <- Metrics_reach$SUMD955sd/Metrics_reach$SUMD955mean*100
Metrics_reach$SUMD75COV <- Metrics_reach$SUMD75sd/Metrics_reach$SUMD75mean*100

Metrics_reach$MINECOV <- Metrics_reach$MINEsd/Metrics_reach$MINEmean*100
Metrics_reach$MINE955COV <- Metrics_reach$MINE955sd/Metrics_reach$MINE955mean*100
Metrics_reach$MINE75COV <- Metrics_reach$MINE75sd/Metrics_reach$MINE75mean*100

Metrics_reach$MAXECOV <- Metrics_reach$MAXEsd/Metrics_reach$MAXEmean*100
Metrics_reach$MAXE955COV <- Metrics_reach$MAXE955sd/Metrics_reach$MAXE955mean*100
Metrics_reach$MAXE75COV <- Metrics_reach$MAXE75sd/Metrics_reach$MAXE75mean*100

Metrics_reach$RNGECOV <- Metrics_reach$RNGEsd/Metrics_reach$RNGEmean*100
Metrics_reach$RNGE955COV <- Metrics_reach$RNGE955sd/Metrics_reach$RNGE955mean*100
Metrics_reach$RNGE75COV <- Metrics_reach$RNGE75sd/Metrics_reach$RNGE75mean*100

Metrics_reach$MEANECOV <- Metrics_reach$MEANEsd/Metrics_reach$MEANEmean*100
Metrics_reach$MEANE95COV <- Metrics_reach$MEANE95sd/Metrics_reach$MEANE95mean*100
Metrics_reach$MEANE75COV <- Metrics_reach$MEANE75sd/Metrics_reach$MEANE75mean*100

Metrics_reach$STDECOV <- Metrics_reach$STDEsd/Metrics_reach$STDEmean*100
Metrics_reach$STDE955COV <- Metrics_reach$STDE955sd/Metrics_reach$STDE955mean*100
Metrics_reach$STDE75COV <- Metrics_reach$STDE75sd/Metrics_reach$STDE75mean*100

Metrics_reach$SUMECOV <- Metrics_reach$SUMEsd/Metrics_reach$SUMEmean*100
Metrics_reach$SUME955COV <- Metrics_reach$SUME955sd/Metrics_reach$SUME955mean*100
Metrics_reach$SUME75COV <- Metrics_reach$SUME75sd/Metrics_reach$SUME75mean*100

#export back to .dbf
Metrics_reach <- as.data.frame(Metrics_reach)

# write out output
file_name = basename(tools::file_path_sans_ext(shp)) # e.g. breach_bankpts_flagged
output_file = file.path(dirname(shp), paste0(file_name, '_all_stats.dbf'))
write.dbf(Metrics_reach, output_file)