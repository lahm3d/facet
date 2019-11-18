## purpose: remove data above and below threshold percentiles and summarize stats by reach
# options(echo=TRUE)

library(foreign)
library(dplyr)
library(pastecs)

### summarize channel metrics from bank pixel method and floodplain metrics from 2D calculation
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
Metrics$chnwid_px[Metrics$chnwid_px <= 0] = NA
Metrics$chnwid1_px[Metrics$chnwid1_px <= 0] = NA
Metrics$chnwid2_px[Metrics$chnwid2_px <= 0] = NA
Metrics$dist_sl[Metrics$dist_sl <= 0] = NA
Metrics$dist[Metrics$dist <= 0] = NA
Metrics$sinuosity[Metrics$sinuosity <= 0] = NA
Metrics$strmorder[Metrics$strmorder <= 0] = NA
Metrics$fpwid_2dc[Metrics$fpwid_2dc <= 0] = NA
Metrics$rngd_2dc[Metrics$rngd_2dc <= 0] = NA
Metrics$mind_2dc[Metrics$mind_2dc <= 0] = NA
Metrics$maxd_2dc[Metrics$maxd_2dc <= 0] = NA
Metrics$stdd_2dc[Metrics$stdd_2dc <= 0] = NA
Metrics$rug_2dc[Metrics$rug_2dc <= 0] = NA
Metrics$bankht_2dh[Metrics$bankht_2dh <= 0] = NA
Metrics$chnshp_2dh[Metrics$chnshp_2dh <= 0] = NA
Metrics$chnwid_2dh[Metrics$chnwid_2dh <= 0] = NA
Metrics$mind_2dh[Metrics$mind_2dh <= 0] = NA
Metrics$maxd_2dh[Metrics$maxd_2dh <= 0] = NA
Metrics$stdd_2dh[Metrics$stdd_2dh <= 0] = NA
Metrics$fpwid_2dh[Metrics$fpwid_2dh <= 0] = NA
Metrics$rug_2dh[Metrics$rug_2dh <= 0] = NA
Metrics$rngd_2dh[Metrics$rngd_2dh <= 0] = NA
Metrics$mine_2dh[Metrics$mine_2dh <= 0] = NA
Metrics$maxe_2dh[Metrics$maxe_2dh <= 0] = NA
Metrics$stde_2dh[Metrics$stde_2dh <= 0] = NA
Metrics$rnge_2dh[Metrics$rnge_2dh <= 0] = NA

# add new columns for additional filtering to be calculated, while keeping original fields too
# 95th/5th
Metrics$CWTPX_9505 <- Metrics$chnwid_px
Metrics$CW1PX_9505 <- Metrics$chnwid1_px
Metrics$CW2PX_9505 <- Metrics$chnwid2_px
Metrics$FPW2_9505 <- Metrics$fpwid_2dc
Metrics$FPR2_9505 <- Metrics$rngd_2dc
Metrics$FPMI2_9505 <- Metrics$mind_2dc
Metrics$FPMA2_9505 <- Metrics$maxd_2dc
Metrics$FPST2_9505 <- Metrics$stdd_2dc
Metrics$FPRU2_9505 <- Metrics$rug_2dc
Metrics$BKHT_9505 <- Metrics$bankht_2dh
Metrics$CHSHP_9505 <- Metrics$chnshp_2dh
Metrics$CWPX3_9505 <- Metrics$chnwid_2dh
Metrics$FPMI3_9505 <- Metrics$mind_2dh
Metrics$FPMA3_9505 <- Metrics$maxd_2dh
Metrics$FPST3_9505 <- Metrics$stdd_2dh
Metrics$FPW3_9505 <- Metrics$fpwid_2dh
Metrics$FPRU3_9505 <- Metrics$rug_2dh
Metrics$FPRN3_9505 <- Metrics$rngd_2dh
Metrics$FMI3E_9505 <- Metrics$mine_2dh
Metrics$FMA3E_9505 <- Metrics$maxe_2dh
Metrics$FST3E_9505 <- Metrics$stde_2dh
Metrics$FRN3E_9505 <- Metrics$rnge_2dh

# 75th/25th
Metrics$CWTPX_7525 <- Metrics$chnwid_px
Metrics$CW1PX_7525 <- Metrics$chnwid1_px
Metrics$CW2PX_7525 <- Metrics$chnwid2_px
Metrics$FPW2_7525 <- Metrics$fpwid_2dc
Metrics$FPR2_7525 <- Metrics$rngd_2dc
Metrics$FPMI2_7525 <- Metrics$mind_2dc
Metrics$FPMA2_7525 <- Metrics$maxd_2dc
Metrics$FPST2_7525 <- Metrics$stdd_2dc
Metrics$FPRU2_7525 <- Metrics$rug_2dc
Metrics$BKHT_7525 <- Metrics$bankht_2dh
Metrics$CHSHP_7525 <- Metrics$chnshp_2dh
Metrics$CWPX3_7525 <- Metrics$chnwid_2dh
Metrics$FPMI3_7525 <- Metrics$mind_2dh
Metrics$FPMA3_7525 <- Metrics$maxd_2dh
Metrics$FPST3_7525 <- Metrics$stdd_2dh
Metrics$FPW3_7525 <- Metrics$fpwid_2dh
Metrics$FPRU3_7525 <- Metrics$rug_2dh
Metrics$FPRN3_7525 <- Metrics$rngd_2dh
Metrics$FMI3E_7525 <- Metrics$mine_2dh
Metrics$FMA3E_7525 <- Metrics$maxe_2dh
Metrics$FST3E_7525 <- Metrics$stde_2dh
Metrics$FRN3E_7525 <- Metrics$rnge_2dh

# group metrics by each reach to filter out outliers
linkno_list = unique(Metrics$linkno)
data_list = list()
for (i in seq_along(linkno_list)) {
    dat <- filter(Metrics, linkno == linkno_list[i])

    # replace all metrics >95th and <5th percentiles w/ NAs
    dat$CWTPX_9505[dat$CWTPX_9505 >= quantile(dat$CWTPX_9505, c(0.95), na.rm = T) | dat$CWTPX_9505 <= quantile(dat$CWTPX_9505, c(0.05), na.rm = T)] = NA
    dat$CW1PX_9505[dat$CW1PX_9505 >= quantile(dat$CW1PX_9505, c(0.95), na.rm = T) | dat$CW1PX_9505 <= quantile(dat$CW1PX_9505, c(0.05), na.rm = T)] = NA
    dat$CW2PX_9505[dat$CW2PX_9505 >= quantile(dat$CW2PX_9505, c(0.95), na.rm = T) | dat$CW2PX_9505 <= quantile(dat$CW2PX_9505, c(0.05), na.rm = T)] = NA
    dat$FPW2_9505[dat$FPW2_9505 >= quantile(dat$FPW2_9505, c(0.95), na.rm = T) | dat$FPW2_9505 <= quantile(dat$FPW2_9505, c(0.05), na.rm = T)] = NA
    dat$FPR2_9505[dat$FPR2_9505 >= quantile(dat$FPR2_9505, c(0.95), na.rm = T) | dat$FPR2_9505 <= quantile(dat$FPR2_9505, c(0.05), na.rm = T)] = NA
    dat$FPMI2_9505[dat$FPMI2_9505 >= quantile(dat$FPMI2_9505, c(0.95), na.rm = T) | dat$FPMI2_9505 <= quantile(dat$FPMI2_9505, c(0.05), na.rm = T)] = NA
    dat$FPMA2_9505[dat$FPMA2_9505 >= quantile(dat$FPMA2_9505, c(0.95), na.rm = T) | dat$FPMA2_9505 <= quantile(dat$FPMA2_9505, c(0.05), na.rm = T)] = NA
    dat$FPST2_9505[dat$FPST2_9505 >= quantile(dat$FPST2_9505, c(0.95), na.rm = T) | dat$FPST2_9505 <= quantile(dat$FPST2_9505, c(0.05), na.rm = T)] = NA
    dat$FPRU2_9505[dat$FPRU2_9505 >= quantile(dat$FPRU2_9505, c(0.95), na.rm = T) | dat$FPRU2_9505 <= quantile(dat$FPRU2_9505, c(0.05), na.rm = T)] = NA
    dat$BKHT_9505[dat$BKHT_9505 >= quantile(dat$BKHT_9505, c(0.95), na.rm = T) | dat$BKHT_9505 <= quantile(dat$BKHT_9505, c(0.05), na.rm = T)] = NA
    dat$CHSHP_9505[dat$CHSHP_9505 >= quantile(dat$CHSHP_9505, c(0.95), na.rm = T) | dat$CHSHP_9505 <= quantile(dat$CHSHP_9505, c(0.05), na.rm = T)] = NA
    dat$CWPX3_9505[dat$CWPX3_9505 >= quantile(dat$CWPX3_9505, c(0.95), na.rm = T) | dat$CWPX3_9505 <= quantile(dat$CWPX3_9505, c(0.05), na.rm = T)] = NA
    dat$FPMI3_9505[dat$FPMI3_9505 >= quantile(dat$FPMI3_9505, c(0.95), na.rm = T) | dat$FPMI3_9505 <= quantile(dat$FPMI3_9505, c(0.05), na.rm = T)] = NA
    dat$FPMA3_9505[dat$FPMA3_9505 >= quantile(dat$FPMA3_9505, c(0.95), na.rm = T) | dat$FPMA3_9505 <= quantile(dat$FPMA3_9505, c(0.05), na.rm = T)] = NA
    dat$FPST3_9505[dat$FPST3_9505 >= quantile(dat$FPST3_9505, c(0.95), na.rm = T) | dat$FPST3_9505 <= quantile(dat$FPST3_9505, c(0.05), na.rm = T)] = NA
    dat$FPW3_9505[dat$FPW3_9505 >= quantile(dat$FPW3_9505, c(0.95), na.rm = T) | dat$FPW3_9505 <= quantile(dat$FPW3_9505, c(0.05), na.rm = T)] = NA
    dat$FPRU3_9505[dat$FPRU3_9505 >= quantile(dat$FPRU3_9505, c(0.95), na.rm = T) | dat$FPRU3_9505 <= quantile(dat$FPRU3_9505, c(0.05), na.rm = T)] = NA
    dat$FPRN3_9505[dat$FPRN3_9505 >= quantile(dat$FPRN3_9505, c(0.95), na.rm = T) | dat$FPRN3_9505 <= quantile(dat$FPRN3_9505, c(0.05), na.rm = T)] = NA
    dat$FMI3E_9505[dat$FMI3E_9505 >= quantile(dat$FMI3E_9505, c(0.95), na.rm = T) | dat$FMI3E_9505 <= quantile(dat$FMI3E_9505, c(0.05), na.rm = T)] = NA
    dat$FMA3E_9505[dat$FMA3E_9505 >= quantile(dat$FMA3E_9505, c(0.95), na.rm = T) | dat$FMA3E_9505 <= quantile(dat$FMA3E_9505, c(0.05), na.rm = T)] = NA
    dat$FST3E_9505[dat$FST3E_9505 >= quantile(dat$FST3E_9505, c(0.95), na.rm = T) | dat$FST3E_9505 <= quantile(dat$FST3E_9505, c(0.05), na.rm = T)] = NA
    dat$FRN3E_9505[dat$FRN3E_9505 >= quantile(dat$FRN3E_9505, c(0.95), na.rm = T) | dat$FRN3E_9505 <= quantile(dat$FRN3E_9505, c(0.05), na.rm = T)] = NA

    # replace all dat >75th and <25th percentiles w/ NAs
    dat$CWTPX_7525[dat$CWTPX_7525 >= quantile(dat$CWTPX_7525, c(0.75), na.rm = T) | dat$CWTPX_7525 <= quantile(dat$CWTPX_7525, c(0.25), na.rm = T)] = NA
    dat$CW1PX_7525[dat$CW1PX_7525 >= quantile(dat$CW1PX_7525, c(0.75), na.rm = T) | dat$CW1PX_7525 <= quantile(dat$CW1PX_7525, c(0.25), na.rm = T)] = NA
    dat$CW2PX_7525[dat$CW2PX_7525 >= quantile(dat$CW2PX_7525, c(0.75), na.rm = T) | dat$CW2PX_7525 <= quantile(dat$CW2PX_7525, c(0.25), na.rm = T)] = NA
    dat$FPW2_7525[dat$FPW2_7525 >= quantile(dat$FPW2_7525, c(0.75), na.rm = T) | dat$FPW2_7525 <= quantile(dat$FPW2_7525, c(0.25), na.rm = T)] = NA
    dat$FPR2_7525[dat$FPR2_7525 >= quantile(dat$FPR2_7525, c(0.75), na.rm = T) | dat$FPR2_7525 <= quantile(dat$FPR2_7525, c(0.25), na.rm = T)] = NA
    dat$FPMI2_7525[dat$FPMI2_7525 >= quantile(dat$FPMI2_7525, c(0.75), na.rm = T) | dat$FPMI2_7525 <= quantile(dat$FPMI2_7525, c(0.25), na.rm = T)] = NA
    dat$FPMA2_7525[dat$FPMA2_7525 >= quantile(dat$FPMA2_7525, c(0.75), na.rm = T) | dat$FPMA2_7525 <= quantile(dat$FPMA2_7525, c(0.25), na.rm = T)] = NA
    dat$FPST2_7525[dat$FPST2_7525 >= quantile(dat$FPST2_7525, c(0.75), na.rm = T) | dat$FPST2_7525 <= quantile(dat$FPST2_7525, c(0.25), na.rm = T)] = NA
    dat$FPRU2_7525[dat$FPRU2_7525 >= quantile(dat$FPRU2_7525, c(0.75), na.rm = T) | dat$FPRU2_7525 <= quantile(dat$FPRU2_7525, c(0.25), na.rm = T)] = NA
    dat$BKHT_7525[dat$BKHT_7525 >= quantile(dat$BKHT_7525, c(0.75), na.rm = T) | dat$BKHT_7525 <= quantile(dat$BKHT_7525, c(0.25), na.rm = T)] = NA
    dat$CHSHP_7525[dat$CHSHP_7525 >= quantile(dat$CHSHP_7525, c(0.75), na.rm = T) | dat$CHSHP_7525 <= quantile(dat$CHSHP_7525, c(0.25), na.rm = T)] = NA
    dat$CWPX3_7525[dat$CWPX3_7525 >= quantile(dat$CWPX3_7525, c(0.75), na.rm = T) | dat$CWPX3_7525 <= quantile(dat$CWPX3_7525, c(0.25), na.rm = T)] = NA
    dat$FPMI3_7525[dat$FPMI3_7525 >= quantile(dat$FPMI3_7525, c(0.75), na.rm = T) | dat$FPMI3_7525 <= quantile(dat$FPMI3_7525, c(0.25), na.rm = T)] = NA
    dat$FPMA3_7525[dat$FPMA3_7525 >= quantile(dat$FPMA3_7525, c(0.75), na.rm = T) | dat$FPMA3_7525 <= quantile(dat$FPMA3_7525, c(0.25), na.rm = T)] = NA
    dat$FPST3_7525[dat$FPST3_7525 >= quantile(dat$FPST3_7525, c(0.75), na.rm = T) | dat$FPST3_7525 <= quantile(dat$FPST3_7525, c(0.25), na.rm = T)] = NA
    dat$FPW3_7525[dat$FPW3_7525 >= quantile(dat$FPW3_7525, c(0.75), na.rm = T) | dat$FPW3_7525 <= quantile(dat$FPW3_7525, c(0.25), na.rm = T)] = NA
    dat$FPRU3_7525[dat$FPRU3_7525 >= quantile(dat$FPRU3_7525, c(0.75), na.rm = T) | dat$FPRU3_7525 <= quantile(dat$FPRU3_7525, c(0.25), na.rm = T)] = NA
    dat$FPRN3_7525[dat$FPRN3_7525 >= quantile(dat$FPRN3_7525, c(0.75), na.rm = T) | dat$FPRN3_7525 <= quantile(dat$FPRN3_7525, c(0.25), na.rm = T)] = NA
    dat$FMI3E_7525[dat$FMI3E_7525 >= quantile(dat$FMI3E_7525, c(0.75), na.rm = T) | dat$FMI3E_7525 <= quantile(dat$FMI3E_7525, c(0.25), na.rm = T)] = NA
    dat$FMA3E_7525[dat$FMA3E_7525 >= quantile(dat$FMA3E_7525, c(0.75), na.rm = T) | dat$FMA3E_7525 <= quantile(dat$FMA3E_7525, c(0.25), na.rm = T)] = NA
    dat$FST3E_7525[dat$FST3E_7525 >= quantile(dat$FST3E_7525, c(0.75), na.rm = T) | dat$FST3E_7525 <= quantile(dat$FST3E_7525, c(0.25), na.rm = T)] = NA
    dat$FRN3E_7525[dat$FRN3E_7525 >= quantile(dat$FRN3E_7525, c(0.75), na.rm = T) | dat$FRN3E_7525 <= quantile(dat$FRN3E_7525, c(0.25), na.rm = T)] = NA

    data_list[[i]] <- dat

    }

  tmp_metrics <- bind_rows(data_list)

#summarize each reach/linkno by median, mean, IQR, sd, and mean absolute deviation
Metrics_reach <- tmp_metrics %>%
    group_by(linkno) %>%
    summarise(mean(chnwid_px, na.rm = T), sd(chnwid_px, na.rm = T), median(chnwid_px, na.rm = T), IQR(chnwid_px, na.rm = T),
            median(CWTPX_7525, na.rm = T), mean(CWTPX_7525, na.rm = T), IQR(CWTPX_7525, na.rm = T), mad(CWTPX_7525, na.rm = T), sd(CWTPX_7525, na.rm = T),
            median(CWTPX_9505, na.rm = T), mean(CWTPX_9505, na.rm = T), IQR(CWTPX_9505, na.rm = T), mad(CWTPX_9505, na.rm = T), sd(CWTPX_9505, na.rm = T),

            mean(chnwid1_px, na.rm = T), sd(chnwid1_px, na.rm = T), median(chnwid1_px, na.rm = T), IQR(chnwid1_px, na.rm = T),
            median(CW1PX_7525, na.rm = T), mean(CW1PX_7525, na.rm = T), IQR(CW1PX_7525, na.rm = T), mad(CW1PX_7525, na.rm = T), sd(CW1PX_7525, na.rm = T),
            median(CW1PX_9505, na.rm = T), mean(CW1PX_9505, na.rm = T), IQR(CW1PX_9505, na.rm = T), mad(CW1PX_9505, na.rm = T), sd(CW1PX_9505, na.rm = T),

            mean(chnwid2_px, na.rm = T), sd(chnwid2_px, na.rm = T), median(chnwid2_px, na.rm = T), IQR(chnwid2_px, na.rm = T),
            median(CW2PX_7525, na.rm = T), mean(CW2PX_7525, na.rm = T), IQR(CW2PX_7525, na.rm = T), mad(CW2PX_7525, na.rm = T), sd(CW2PX_7525, na.rm = T),
            median(CW2PX_9505, na.rm = T), mean(CW2PX_9505, na.rm = T), IQR(CW2PX_9505, na.rm = T), mad(CW2PX_9505, na.rm = T), sd(CW2PX_9505, na.rm = T),

            mean(dist_sl, na.rm = T), sd(dist_sl, na.rm = T),median(dist_sl, na.rm = T), IQR(dist_sl, na.rm = T),

            mean(dist, na.rm = T), sd(dist, na.rm = T), median(dist, na.rm = T), IQR(dist, na.rm = T),

            mean(sinuosity, na.rm = T), sd(sinuosity, na.rm = T), median(sinuosity, na.rm = T), IQR(sinuosity, na.rm = T),

            mean(strmorder, na.rm = T),

            mean(fpwid_2dc, na.rm = T), sd(fpwid_2dc, na.rm = T), median(fpwid_2dc, na.rm = T), IQR(fpwid_2dc, na.rm = T),
            median(FPW2_7525, na.rm = T), mean(FPW2_7525, na.rm = T), IQR(FPW2_7525, na.rm = T), mad(FPW2_7525, na.rm = T), sd(FPW2_7525, na.rm = T),
            median(FPW2_9505, na.rm = T), mean(FPW2_9505, na.rm = T), IQR(FPW2_9505, na.rm = T), mad(FPW2_9505, na.rm = T), sd(FPW2_9505, na.rm = T),

            mean(rngd_2dc, na.rm = T), sd(rngd_2dc, na.rm = T), median(rngd_2dc, na.rm = T), IQR(rngd_2dc, na.rm = T),
            median(FPR2_7525, na.rm = T), mean(FPR2_7525, na.rm = T), IQR(FPR2_7525, na.rm = T), mad(FPR2_7525, na.rm = T), sd(FPR2_7525, na.rm = T),
            median(FPR2_9505, na.rm = T), mean(FPR2_9505, na.rm = T), IQR(FPR2_9505, na.rm = T), mad(FPR2_9505, na.rm = T), sd(FPR2_9505, na.rm = T),

            mean(mind_2dc, na.rm = T), sd(mind_2dc, na.rm = T), median(mind_2dc, na.rm = T), IQR(mind_2dc, na.rm = T),
            median(FPMI2_7525, na.rm = T), mean(FPMI2_7525, na.rm = T), IQR(FPMI2_7525, na.rm = T), mad(FPMI2_7525, na.rm = T), sd(FPMI2_7525, na.rm = T),
            median(FPMI2_9505, na.rm = T), mean(FPMI2_9505, na.rm = T), IQR(FPMI2_9505, na.rm = T), mad(FPMI2_9505, na.rm = T), sd(FPMI2_9505, na.rm = T),

            mean(maxd_2dc, na.rm = T), sd(maxd_2dc, na.rm = T), median(maxd_2dc, na.rm = T), IQR(maxd_2dc, na.rm = T),
            median(FPMA2_7525, na.rm = T), mean(FPMA2_7525, na.rm = T), IQR(FPMA2_7525, na.rm = T), mad(FPMA2_7525, na.rm = T), sd(FPMA2_7525, na.rm = T),
            median(FPMA2_9505, na.rm = T), mean(FPMA2_9505, na.rm = T), IQR(FPMA2_9505, na.rm = T), mad(FPMA2_9505, na.rm = T), sd(FPMA2_9505, na.rm = T),

            mean(stdd_2dc, na.rm = T), sd(stdd_2dc, na.rm = T), median(stdd_2dc, na.rm = T), IQR(stdd_2dc, na.rm = T),
            median(FPST2_7525, na.rm = T), mean(FPST2_7525, na.rm = T), IQR(FPST2_7525, na.rm = T), mad(FPST2_7525, na.rm = T), sd(FPST2_7525, na.rm = T),
            median(FPST2_9505, na.rm = T), mean(FPST2_9505, na.rm = T), IQR(FPST2_9505, na.rm = T), mad(FPST2_9505, na.rm = T), sd(FPST2_9505, na.rm = T),

            mean(rug_2dc, na.rm = T), sd(rug_2dc, na.rm = T), median(rug_2dc, na.rm = T), IQR(rug_2dc, na.rm = T),
            median(FPRU2_7525, na.rm = T), mean(FPRU2_7525, na.rm = T), IQR(FPRU2_7525, na.rm = T), mad(FPRU2_7525, na.rm = T), sd(FPRU2_7525, na.rm = T),
            median(FPRU2_9505, na.rm = T), mean(FPRU2_9505, na.rm = T), IQR(FPRU2_9505, na.rm = T), mad(FPRU2_9505, na.rm = T), sd(FPRU2_9505, na.rm = T),

            mean(bankht_2dh, na.rm = T), sd(bankht_2dh, na.rm = T), median(bankht_2dh, na.rm = T), IQR(bankht_2dh, na.rm = T),
            median(BKHT_7525, na.rm = T), mean(BKHT_7525, na.rm = T), IQR(BKHT_7525, na.rm = T), mad(BKHT_7525, na.rm = T), sd(BKHT_7525, na.rm = T),
            median(BKHT_9505, na.rm = T), mean(BKHT_9505, na.rm = T), IQR(BKHT_9505, na.rm = T), mad(BKHT_9505, na.rm = T), sd(BKHT_9505, na.rm = T),

            mean(chnshp_2dh, na.rm = T), sd(chnshp_2dh, na.rm = T), median(chnshp_2dh, na.rm = T), IQR(chnshp_2dh, na.rm = T),
            median(CHSHP_7525, na.rm = T), mean(CHSHP_7525, na.rm = T), IQR(CHSHP_7525, na.rm = T), mad(CHSHP_7525, na.rm = T), sd(CHSHP_7525, na.rm = T),
            median(CHSHP_9505, na.rm = T), mean(CHSHP_9505, na.rm = T), IQR(CHSHP_9505, na.rm = T), mad(CHSHP_9505, na.rm = T), sd(CHSHP_9505, na.rm = T),

            mean(chnwid_2dh, na.rm = T), sd(chnwid_2dh, na.rm = T), median(chnwid_2dh, na.rm = T), IQR(chnwid_2dh, na.rm = T),
            median(CWPX3_7525, na.rm = T), mean(CWPX3_7525, na.rm = T), IQR(CWPX3_7525, na.rm = T), mad(CWPX3_7525, na.rm = T), sd(CWPX3_7525, na.rm = T),
            median(CWPX3_9505, na.rm = T), mean(CWPX3_9505, na.rm = T), IQR(CWPX3_9505, na.rm = T), mad(CWPX3_9505, na.rm = T), sd(CWPX3_9505, na.rm = T),

            mean(mind_2dh, na.rm = T), sd(mind_2dh, na.rm = T), median(mind_2dh, na.rm = T), IQR(mind_2dh, na.rm = T),
            median(FPMI3_7525, na.rm = T), mean(FPMI3_7525, na.rm = T), IQR(FPMI3_7525, na.rm = T), mad(FPMI3_7525, na.rm = T), sd(FPMI3_7525, na.rm = T),
            median(FPMI3_9505, na.rm = T), mean(FPMI3_9505, na.rm = T), IQR(FPMI3_9505, na.rm = T), mad(FPMI3_9505, na.rm = T), sd(FPMI3_9505, na.rm = T),

            mean(maxd_2dh, na.rm = T), sd(maxd_2dh, na.rm = T), median(maxd_2dh, na.rm = T), IQR(maxd_2dh, na.rm = T),
            median(FPMA3_7525, na.rm = T), mean(FPMA3_7525, na.rm = T), IQR(FPMA3_7525, na.rm = T), mad(FPMA3_7525, na.rm = T), sd(FPMA3_7525, na.rm = T),
            median(FPMA3_9505, na.rm = T), mean(FPMA3_9505, na.rm = T), IQR(FPMA3_9505, na.rm = T), mad(FPMA3_9505, na.rm = T), sd(FPMA3_9505, na.rm = T),

            mean(stdd_2dh, na.rm = T), sd(stdd_2dh, na.rm = T), median(stdd_2dh, na.rm = T), IQR(stdd_2dh, na.rm = T),
            median(FPST3_7525, na.rm = T), mean(FPST3_7525, na.rm = T), IQR(FPST3_7525, na.rm = T), mad(FPST3_7525, na.rm = T), sd(FPST3_7525, na.rm = T),
            median(FPST3_9505, na.rm = T), mean(FPST3_9505, na.rm = T), IQR(FPST3_9505, na.rm = T), mad(FPST3_9505, na.rm = T), sd(FPST3_9505, na.rm = T),

            mean(fpwid_2dh, na.rm = T), sd(fpwid_2dh, na.rm = T), median(fpwid_2dh, na.rm = T), IQR(fpwid_2dh, na.rm = T),
            median(FPW3_7525, na.rm = T), mean(FPW3_7525, na.rm = T), IQR(FPW3_7525, na.rm = T), mad(FPW3_7525, na.rm = T), sd(FPW3_7525, na.rm = T),
            median(FPW3_9505, na.rm = T), mean(FPW3_9505, na.rm = T), IQR(FPW3_9505, na.rm = T), mad(FPW3_9505, na.rm = T), sd(FPW3_9505, na.rm = T),

            mean(rug_2dh, na.rm = T), sd(rug_2dh, na.rm = T), median(rug_2dh, na.rm = T), IQR(rug_2dh, na.rm = T),
            median(FPRU3_7525, na.rm = T), mean(FPRU3_7525, na.rm = T), IQR(FPRU3_7525, na.rm = T), mad(FPRU3_7525, na.rm = T), sd(FPRU3_7525, na.rm = T),
            median(FPRU3_9505, na.rm = T), mean(FPRU3_9505, na.rm = T), IQR(FPRU3_9505, na.rm = T), mad(FPRU3_9505, na.rm = T), sd(FPRU3_9505, na.rm = T),

            mean(rngd_2dh, na.rm = T), sd(rngd_2dh, na.rm = T), median(rngd_2dh, na.rm = T), IQR(rngd_2dh, na.rm = T),
            median(FPRN3_7525, na.rm = T), mean(FPRN3_7525, na.rm = T), IQR(FPRN3_7525, na.rm = T), mad(FPRN3_7525, na.rm = T), sd(FPRN3_7525, na.rm = T),
            median(FPRN3_9505, na.rm = T), mean(FPRN3_9505, na.rm = T), IQR(FPRN3_9505, na.rm = T), mad(FPRN3_9505, na.rm = T), sd(FPRN3_9505, na.rm = T),

            mean(mine_2dh, na.rm = T), sd(mine_2dh, na.rm = T), median(mine_2dh, na.rm = T), IQR(mine_2dh, na.rm = T),
            median(FMI3E_7525, na.rm = T), mean(FMI3E_7525, na.rm = T), IQR(FMI3E_7525, na.rm = T), mad(FMI3E_7525, na.rm = T), sd(FMI3E_7525, na.rm = T),
            median(FMI3E_9505, na.rm = T), mean(FMI3E_9505, na.rm = T), IQR(FMI3E_9505, na.rm = T), mad(FMI3E_9505, na.rm = T), sd(FMI3E_9505, na.rm = T),

            mean(maxe_2dh, na.rm = T), sd(maxe_2dh, na.rm = T), median(maxe_2dh, na.rm = T), IQR(maxe_2dh, na.rm = T),
            median(FMA3E_7525, na.rm = T), mean(FMA3E_7525, na.rm = T), IQR(FMA3E_7525, na.rm = T), mad(FMA3E_7525, na.rm = T), sd(FMA3E_7525, na.rm = T),
            median(FMA3E_9505, na.rm = T), mean(FMA3E_9505, na.rm = T), IQR(FMA3E_9505, na.rm = T), mad(FMA3E_9505, na.rm = T), sd(FMA3E_9505, na.rm = T),

            mean(stde_2dh, na.rm = T), sd(stde_2dh, na.rm = T), median(stde_2dh, na.rm = T), IQR(stde_2dh, na.rm = T),
            median(FST3E_7525, na.rm = T), mean(FST3E_7525, na.rm = T), IQR(FST3E_7525, na.rm = T), mad(FST3E_7525, na.rm = T), sd(FST3E_7525, na.rm = T),
            median(FST3E_9505, na.rm = T), mean(FST3E_9505, na.rm = T), IQR(FST3E_9505, na.rm = T), mad(FST3E_9505, na.rm = T), sd(FST3E_9505, na.rm = T),

            mean(rnge_2dh, na.rm = T), sd(rnge_2dh, na.rm = T), median(rnge_2dh, na.rm = T), IQR(rnge_2dh, na.rm = T),
            median(FRN3E_7525, na.rm = T), mean(FRN3E_7525, na.rm = T), IQR(FRN3E_7525, na.rm = T), mad(FRN3E_7525, na.rm = T), sd(FRN3E_7525, na.rm = T),
            median(FRN3E_9505, na.rm = T), mean(FRN3E_9505, na.rm = T), IQR(FRN3E_9505, na.rm = T), mad(FRN3E_9505, na.rm = T), sd(FRN3E_9505, na.rm = T))

# rename column headings
col_headings <- c("LINKNO",
                "CWTPXmean", "CWTPXsd", "CWTPXmed", "CWTPXIQR",
                "CWTPX75med", "CWTPX75mean", "CWTPX75IQR", "CWTPX75mad", "CWTPX75sd",
                "CWTPX95med", "CWTPX95mean", "CWTPX95IQR", "CWTPX95mad", "CWTPX95sd",

                "CW1PXmean", "CW1PXsd", "CW1PXmed", "CW1PXIQR",
                "CW1PX75med", "CW1PX75mean", "CW1PX75IQR", "CW1PX75mad", "CW1PX75sd",
                "CW1PX95med", "CW1PX95mean", "CW1PX95IQR", "CW1PX95mad", "CW1PX95sd",

                "CW2PXmean", "CW2PXsd", "CW2PXmed", "CW2PXIQR",
                "CW2PX75med", "CW2PX75mean", "CW2PX75IQR", "CW2PX75mad", "CW2PX75sd",
                "CW2PX95med", "CW2PX95mean", "CW2PX95IQR", "CW2PX95mad", "CW2PX95sd",

                "distslmean", "distslsd", "distslmed", "distslIQR",
                "distmean", "distsd", "distmed", "distIQR",
                "sinumean", "sinusd", "sinumed", "sinuIQR",
                "ordmean",

                "FPW2mean", "FPW2sd", "FPW2med", "FPW2IQR",
                "FPW275med", "FPW275mean", "FPW275IQR", "FPW275mad", "FPW275sd",
                "FPW295med", "FPW295mean", "FPW295IQR", "FPW295mad", "FPW295sd",

                "FPR2mean", "FPR2sd", "FPR2med", "FPR2IQR",
                "FPR275med", "FPR275mean", "FPR275IQR", "FPR275mad", "FPR275sd",
                "FPR295med", "FPR295mean", "FPR295IQR", "FPR295mad", "FPR295sd",

                "FPMI2mean", "FPMI2sd", "FPMI2med", "FPMI2IQR",
                "FPMI275med", "FPMI275mean", "FPMI275IQR", "FPMI275mad", "FPMI275sd",
                "FPMI295med", "FPMI295mean", "FPMI295IQR", "FPMI295mad", "FPMI295sd",

                "FPMA2mean", "FPMA2sd", "FPMA2med", "FPMA2IQR",
                "FPMA275med", "FPMA275mean", "FPMA275IQR", "FPMA275mad", "FPMA275sd",
                "FPMA295med", "FPMA295mean", "FPMA295IQR", "FPMA295mad", "FPMA295sd",

                "FPST2mean", "FPST2sd", "FPST2med", "FPST2IQR",
                "FPST275med", "FPST275mean", "FPST275IQR", "FPST275mad", "FPST275sd",
                "FPST295med", "FPST295mean", "FPST295IQR", "FPST295mad", "FPST295sd",

                "FPRU2mean", "FPRU2sd", "FPRU2med", "FPRU2IQR",
                "FPRU275med", "FPRU275mean", "FPRU275IQR", "FPRU275mad", "FPRU275sd",
                "FPRU295med", "FPRU295mean", "FPRU295IQR", "FPRU295mad", "FPRU295sd",

                "BKHTmean", "BKHTsd", "BKHTmed", "BKHTIQR",
                "BKHT75med", "BKHT75mean", "BKHT75IQR", "BKHT75mad", "BKHT75sd",
                "BKHT95med", "BKHT95mean", "BKHT95IQR", "BKHT95mad", "BKHT95sd",

                "CHSHPmean", "CHSHPsd", "CHSHPmed", "CHSHPIQR",
                "CHSHP75med", "CHSHP75mean", "CHSHP75IQR", "CHSHP75mad", "CHSHP75sd",
                "CHSHP95med", "CHSHP95mean", "CHSHP95IQR", "CHSHP95mad", "CHSHP95sd",

                "CWPX3mean", "CWPX3sd", "CWPX3med", "CWPX3IQR",
                "CWPX375med", "CWPX375mean", "CWPX375IQR", "CWPX375mad", "CWPX375sd",
                "CWPX395med", "CWPX395mean", "CWPX395IQR", "CWPX395mad", "CWPX395sd",

                "FPMI3mean", "FPMI3sd", "FPMI3med", "FPMI3IQR",
                "FPMI375med", "FPMI375mean", "FPMI375IQR", "FPMI375mad", "FPMI375sd",
                "FPMI395med", "FPMI395mean", "FPMI395IQR", "FPMI395mad", "FPMI395sd",

                "FPMA3mean", "FPMA3sd", "FPMA3med", "FPMA3IQR",
                "FPMA375med", "FPMA375mean", "FPMA375IQR", "FPMA375mad", "FPMA375sd",
                "FPMA395med", "FPMA395mean", "FPMA395IQR", "FPMA395mad", "FPMA395sd",

                "FPST3mean", "FPST3sd", "FPST3med", "FPST3IQR",
                "FPST375med", "FPST375mean", "FPST375IQR", "FPST375mad", "FPST375sd",
                "FPST395med", "FPST395mean", "FPST395IQR", "FPST395mad", "FPST395sd",

                "FPW3mean", "FPW3sd", "FPW3med", "FPW3IQR",
                "FPW375med", "FPW375mean", "FPW375IQR", "FPW375mad", "FPW375sd",
                "FPW395med", "FPW395mean", "FPW395IQR", "FPW395mad", "FPW395sd",

                "FPRU3mean", "FPRU3sd", "FPRU3med", "FPRU3IQR",
                "FPRU375med", "FPRU375mean", "FPRU375IQR", "FPRU375mad", "FPRU375sd",
                "FPRU395med", "FPRU395mean", "FPRU395IQR", "FPRU395mad", "FPRU395sd",

                "FPRN3mean", "FPRN3sd", "FPRN3med", "FPRN3IQR",
                "FPRN375med", "FPRN375mean", "FPRN375IQR", "FPRN375mad", "FPRN375sd",
                "FPRN395med", "FPRN395mean", "FPRN395IQR", "FPRN395mad", "FPRN395sd",

                "FMI3Emean", "FMI3Esd", "FMI3Emed", "FMI3EIQR",
                "FMI3E75med", "FMI3E75mean", "FMI3E75IQR", "FMI3E75mad", "FMI3E75sd",
                "FMI3E95med", "FMI3E95mean", "FMI3E95IQR", "FMI3E95mad", "FMI3E95sd",

                "FMA3Emean", "FMA3Esd", "FMA3Emed", "FMA3EIQR",
                "FMA3E75med", "FMA3E75mean", "FMA3E75IQR", "FMA3E75mad", "FMA3E75sd",
                "FMA3E95med", "FMA3E95mean", "FMA3E95IQR", "FMA3E95mad", "FMA3E95sd",

                "FST3Emean", "FST3Esd", "FST3Emed", "FST3EIQR",
                "FST3E75med", "FST3E75mean", "FST3E75IQR", "FST3E75mad", "FST3E75sd",
                "FST3E95med", "FST3E95mean", "FST3E95IQR", "FST3E95mad", "FST3E95sd",

                "FRN3Emean", "FRN3Esd", "FRN3Emed", "FRN3EIQR",
                "FRN3E75med", "FRN3E75mean", "FRN3E75IQR", "FRN3E75mad", "FRN3E75sd",
                "FRN3E95med", "FRN3E95mean", "FRN3E95IQR", "FRN3E95mad", "FRN3E95sd")

names(Metrics_reach) <- col_headings

Metrics_reach$CWTPXCOV <- Metrics_reach$CWTPXsd/Metrics_reach$CWTPXmean*100
Metrics_reach$CWTPX75COV <- Metrics_reach$CWTPX75sd/Metrics_reach$CWTPX75mean*100
Metrics_reach$CWTPX95COV <- Metrics_reach$CWTPX95sd/Metrics_reach$CWTPX95mean*100

Metrics_reach$CW1PXCOV <- Metrics_reach$CW1PXsd/Metrics_reach$CW1PXmean*100
Metrics_reach$CW1PX75COV <- Metrics_reach$CW1PX75sd/Metrics_reach$CW1PX75mean*100
Metrics_reach$CW1PX95COV <- Metrics_reach$CW1PX95sd/Metrics_reach$CW1PX95mean*100

Metrics_reach$CW2PXCOV <- Metrics_reach$CW2PXsd/Metrics_reach$CW2PXmean*100
Metrics_reach$CW2PX75COV <- Metrics_reach$CW2PX75sd/Metrics_reach$CW2PX75mean*100
Metrics_reach$CW2PX95COV <- Metrics_reach$CW2PX95sd/Metrics_reach$CW2PX95mean*100

Metrics_reach$FPW2COV <- Metrics_reach$FPW2sd/Metrics_reach$FPW2mean*100
Metrics_reach$FPW275COV <- Metrics_reach$FPW275sd/Metrics_reach$FPW275mean*100
Metrics_reach$FPW295COV <- Metrics_reach$FPW295sd/Metrics_reach$FPW295mean*100

Metrics_reach$FPR2COV <- Metrics_reach$FPR2sd/Metrics_reach$FPR2mean*100
Metrics_reach$FPR275COV <- Metrics_reach$FPR275sd/Metrics_reach$FPR275mean*100
Metrics_reach$FPR295COV <- Metrics_reach$FPR295sd/Metrics_reach$FPR295mean*100

Metrics_reach$FPMI2COV <- Metrics_reach$FPMI2sd/Metrics_reach$FPMI2mean*100
Metrics_reach$FPMI275COV <- Metrics_reach$FPMI275sd/Metrics_reach$FPMI275mean*100
Metrics_reach$FPMI295COV <- Metrics_reach$FPMI295sd/Metrics_reach$FPMI295mean*100

Metrics_reach$FPMA2COV <- Metrics_reach$FPMA2sd/Metrics_reach$FPMA2mean*100
Metrics_reach$FPMA275COV <- Metrics_reach$FPMA275sd/Metrics_reach$FPMA275mean*100
Metrics_reach$FPMA295COV <- Metrics_reach$FPMA295sd/Metrics_reach$FPMA295mean*100

Metrics_reach$FPST2COV <- Metrics_reach$FPST2sd/Metrics_reach$FPST2mean*100
Metrics_reach$FPST275COV <- Metrics_reach$FPST275sd/Metrics_reach$FPST275mean*100
Metrics_reach$FPST295COV <- Metrics_reach$FPST295sd/Metrics_reach$FPST295mean*100

Metrics_reach$FPRU2COV <- Metrics_reach$FPRU2sd/Metrics_reach$FPRU2mean*100
Metrics_reach$FPRU275COV <- Metrics_reach$FPRU275sd/Metrics_reach$FPRU275mean*100
Metrics_reach$FPRU295COV <- Metrics_reach$FPRU295sd/Metrics_reach$FPRU295mean*100

Metrics_reach$BKHTCOV <- Metrics_reach$BKHTsd/Metrics_reach$BKHTmean*100
Metrics_reach$BKHT75COV <- Metrics_reach$BKHT75sd/Metrics_reach$BKHT75mean*100
Metrics_reach$BKHT95COV <- Metrics_reach$BKHT95sd/Metrics_reach$BKHT95mean*100

Metrics_reach$CHSHPCOV <- Metrics_reach$CHSHPsd/Metrics_reach$CHSHPmean*100
Metrics_reach$CHSHP75COV <- Metrics_reach$CHSHP75sd/Metrics_reach$CHSHP75mean*100
Metrics_reach$CHSHP95COV <- Metrics_reach$CHSHP95sd/Metrics_reach$CHSHP95mean*100

Metrics_reach$CWPX3COV <- Metrics_reach$CWPX3sd/Metrics_reach$CWPX3mean*100
Metrics_reach$CWPX375COV <- Metrics_reach$CWPX375sd/Metrics_reach$CWPX375mean*100
Metrics_reach$CWPX395COV <- Metrics_reach$CWPX395sd/Metrics_reach$CWPX395mean*100

Metrics_reach$FPMI3COV <- Metrics_reach$FPMI3sd/Metrics_reach$FPMI3mean*100
Metrics_reach$FPMI375COV <- Metrics_reach$FPMI375sd/Metrics_reach$FPMI375mean*100
Metrics_reach$FPMI395COV <- Metrics_reach$FPMI395sd/Metrics_reach$FPMI395mean*100

Metrics_reach$FPMA3COV <- Metrics_reach$FPMA3sd/Metrics_reach$FPMA3mean*100
Metrics_reach$FPMA375COV <- Metrics_reach$FPMA375sd/Metrics_reach$FPMA375mean*100
Metrics_reach$FPMA395COV <- Metrics_reach$FPMA395sd/Metrics_reach$FPMA395mean*100

Metrics_reach$FPST3COV <- Metrics_reach$FPST3sd/Metrics_reach$FPST3mean*100
Metrics_reach$FPST375COV <- Metrics_reach$FPST375sd/Metrics_reach$FPST375mean*100
Metrics_reach$FPST395COV <- Metrics_reach$FPST395sd/Metrics_reach$FPST395mean*100

Metrics_reach$FPW3COV <- Metrics_reach$FPW3sd/Metrics_reach$FPW3mean*100
Metrics_reach$FPW375COV <- Metrics_reach$FPW375sd/Metrics_reach$FPW375mean*100
Metrics_reach$FPW395COV <- Metrics_reach$FPW395sd/Metrics_reach$FPW395mean*100

Metrics_reach$FPRU3COV <- Metrics_reach$FPRU3sd/Metrics_reach$FPRU3mean*100
Metrics_reach$FPRU375COV <- Metrics_reach$FPRU375sd/Metrics_reach$FPRU375mean*100
Metrics_reach$FPRU395COV <- Metrics_reach$FPRU395sd/Metrics_reach$FPRU395mean*100

Metrics_reach$FPRN3COV <- Metrics_reach$FPRN3sd/Metrics_reach$FPRN3mean*100
Metrics_reach$FPRN375COV <- Metrics_reach$FPRN375sd/Metrics_reach$FPRN375mean*100
Metrics_reach$FPRN395COV <- Metrics_reach$FPRN395sd/Metrics_reach$FPRN395mean*100

Metrics_reach$FMI3ECOV <- Metrics_reach$FMI3Esd/Metrics_reach$FMI3Emean*100
Metrics_reach$FMI3E75COV <- Metrics_reach$FMI3E75sd/Metrics_reach$FMI3E75mean*100
Metrics_reach$FMI3E95COV <- Metrics_reach$FMI3E95sd/Metrics_reach$FMI3E95mean*100

Metrics_reach$FMA3ECOV <- Metrics_reach$FMA3Esd/Metrics_reach$FMA3Emean*100
Metrics_reach$FMA3E75COV <- Metrics_reach$FMA3E75sd/Metrics_reach$FMA3E75mean*100
Metrics_reach$FMA3E95COV <- Metrics_reach$FMA3E95sd/Metrics_reach$FMA3E95mean*100

Metrics_reach$FST3ECOV <- Metrics_reach$FST3Esd/Metrics_reach$FST3Emean*100
Metrics_reach$FST3E75COV <- Metrics_reach$FST3E75sd/Metrics_reach$FST3E75mean*100
Metrics_reach$FST3E95COV <- Metrics_reach$FST3E95sd/Metrics_reach$FST3E95mean*100

Metrics_reach$FRN3ECOV <- Metrics_reach$FRN3Esd/Metrics_reach$FRN3Emean*100
Metrics_reach$FRN3E75COV <- Metrics_reach$FRN3E75sd/Metrics_reach$FRN3E75mean*100
Metrics_reach$FRN3E95COV <- Metrics_reach$FRN3E95sd/Metrics_reach$FRN3E95mean*100


#export back to .dbf
Metrics_reach <- as.data.frame(Metrics_reach)

# write out output
file_name = basename(tools::file_path_sans_ext(shp)) # e.g. breach_bankpts_flagged
output_file = file.path(dirname(shp), paste0(file_name, '_all_stats.dbf'))
write.dbf(Metrics_reach, output_file)