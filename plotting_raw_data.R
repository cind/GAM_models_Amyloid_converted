library(ggplot2)
library(dplyr)
library(patchwork)

######################################################################################################################################
root = "gam_modelling_data"

centiloid_plot_data <- read.csv(file.path(root, "centiloid_gam_data.csv"))
tau_plot_data <- read.csv(file.path(root, "tau_gam_data.csv"))
ptau_plot_data <- read.csv(file.path(root, "ptau_gam_data.csv"))
abeta_plot_data <- read.csv(file.path(root, "abeta_gam_data.csv"))
npi_plot_data <- read.csv(file.path(root, "npi_gam_data.csv"))
mpacctrailsb_plot_data <- read.csv(file.path(root, "mpacctrailsb_gam_data.csv"))
cdrsb_plot_data <- read.csv(file.path(root, "cdrsb_gam_data.csv"))
mmse_plot_data <- read.csv(file.path(root, "mmse_gam_data.csv"))
fdg_plot_data <- read.csv(file.path(root, "fdg_gam_data.csv"))
meta_roi_plot_data <- read.csv(file.path(root, "meta_roi_gam_data_no_adni1.csv"))
hippocampal_volume_plot_data <- read.csv(file.path(root, "hippocampal_volume_gam_data_no_adni1.csv"))
adas13_plot_data <- read.csv(file.path(root, "adas13_gam_data.csv"))
ecog_s_plot_data <- read.csv(file.path(root, "ecog_s_gam_data.csv"))
ecog_p_plot_data <- read.csv(file.path(root, "ecog_p_gam_data.csv"))

#####################################################################################################
#Getting CDGLOBAL into datasets where it is missing
#####################################################################################################
#CDGLOBAL
cdglobal_data <- read.csv("~/Data/CDR_24Apr2024.csv")
cdglobal_data$CDR.DATE <- cdglobal_data$USERDATE
cdglobal_data <- cdglobal_data[cdglobal_data$CDGLOBAL>=0,]
cdglobal_data <- cdglobal_data %>%
  dplyr::filter(RID %in% unique(centiloid_plot_data$RID)) %>%
  dplyr::select(RID, CDR.DATE, CDGLOBAL)

#centiloid
centiloid_plot_data <- merge(cdglobal_data, centiloid_plot_data, all = TRUE) %>%
  dplyr::mutate(time_diff = abs(difftime(EXAMDATE_pet, as.Date(CDR.DATE)))) %>%
  dplyr::group_by(RID, EXAMDATE_pet) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()

#fdg
fdg_plot_data <- merge(cdglobal_data, fdg_plot_data, all = TRUE) %>%
  dplyr::mutate(time_diff = abs(difftime(EXAMDATE, as.Date(CDR.DATE)))) %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()

#tau
tau_plot_data <- merge(cdglobal_data, tau_plot_data, all = TRUE) %>%
  dplyr::mutate(time_diff = abs(difftime(EXAMDATE, as.Date(CDR.DATE)))) %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()

#ptau
ptau_plot_data <- merge(cdglobal_data, ptau_plot_data, all = TRUE) %>%
  dplyr::mutate(time_diff = abs(difftime(EXAMDATE, as.Date(CDR.DATE)))) %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()

#abeta
abeta_plot_data <- merge(cdglobal_data, abeta_plot_data, all = TRUE) %>%
  dplyr::mutate(time_diff = abs(difftime(EXAMDATE, as.Date(CDR.DATE)))) %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()

#meta-roi cortical thickness
meta_roi_plot_data <- merge(cdglobal_data, meta_roi_plot_data, all = TRUE) %>%
  dplyr::mutate(time_diff = abs(difftime(EXAMDATE, as.Date(CDR.DATE)))) %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()

#hippocampal volume
hippocampal_volume_plot_data <- merge(cdglobal_data, hippocampal_volume_plot_data, all = TRUE) %>%
  dplyr::mutate(time_diff = abs(difftime(EXAMDATE, as.Date(CDR.DATE)))) %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()

#####################################################################################################
# Centiloid
#####################################################################################################

meas_fn_centiloid <- "PET Aβ Centiloid"

#plotting raw data
raw_centiloid_plot <- ggplot(data = centiloid_plot_data, aes(x = adjusted_new_time, y = Centiloid, group = RID, color = as.factor(CDGLOBAL))) +
  theme_classic() + 
  ylab(unique(meas_fn_centiloid)) + 
  scale_x_continuous(name = "Years from Aβ+", breaks = seq(-10, 10, by = 2)) + #limits=c(-8, 8)) +
  theme(legend.position="none") + 
  geom_point(size = .2) +
  geom_line(size = .2) +
  scale_color_manual(values=c("darkblue", "forestgreen", "darkorange", "red"))

#####################################################################################################
# fdg
#####################################################################################################

meas_fn_fdg <- "FDG PET meta-ROI"

#plotting raw data
raw_fdg_plot <- ggplot(data = fdg_plot_data, aes(x = adjusted_new_time, y = adjusted_Meta_ROI, group = RID, color = as.factor(CDGLOBAL))) +
  theme_classic() + 
  ylab(unique(meas_fn_fdg)) + 
  scale_x_continuous(name = "Years from Aβ+", breaks = seq(-10, 10, by = 2)) + #limits=c(-8, 8)) +
  # theme(legend.position="none") +
  geom_point(size = .2) +
  geom_line(size = .2) +
  scale_color_manual(values=c("darkblue", "forestgreen", "darkorange", "red"))
raw_fdg_plot$labels$colour <- "CDR"


#####################################################################################################
# csf tau
#####################################################################################################

meas_fn_tau <- "CSF total tau"

#plotting raw data
raw_tau_plot <- ggplot(data = tau_plot_data, aes(x = adjusted_new_time, y = TAU, group = RID, color = as.factor(CDGLOBAL))) +
  theme_classic() + 
  ylab(unique(meas_fn_tau)) + 
  scale_x_continuous(name = "Years from Aβ+", breaks = seq(-10, 10, by = 2)) + #limits=c(-8, 8)) +
  theme(legend.position="none") +
  geom_point(size = .2) +
  geom_line(size = .2) +
  scale_color_manual(values=c("darkblue", "forestgreen", "darkorange", "red"))

#####################################################################################################
# csf ptau
#####################################################################################################

meas_fn_ptau <- "CSF p-tau181"

#plotting raw data
raw_ptau_plot <- ggplot(data = ptau_plot_data, aes(x = adjusted_new_time, y = PTAU, group = RID, color = as.factor(CDGLOBAL))) +
  theme_classic() + 
  ylab(unique(meas_fn_ptau)) + 
  scale_x_continuous(name = "Years from Aβ+", breaks = seq(-10, 10, by = 2)) + #limits=c(-8, 8)) +
  theme(legend.position="none") + 
  geom_point(size = .2) +
  geom_line(size = .2) +
  scale_color_manual(values=c("darkblue", "forestgreen", "darkorange", "red"))

#####################################################################################################
# csf Amyloid
#####################################################################################################

meas_fn_abeta <- "CSF Aβ42"

#plotting raw data
raw_abeta_plot <- ggplot(data = abeta_plot_data, aes(x = adjusted_new_time, y = ABETA, group = RID, color = as.factor(CDGLOBAL))) +
  theme_classic() + 
  ylab(unique(meas_fn_abeta)) + 
  scale_x_continuous(name = "Years from Aβ+", breaks = seq(-10, 10, by = 2)) + #limits=c(-8, 8)) +
  theme(legend.position="none") + 
  geom_point(size = .2) +
  geom_line(size = .2) +
  scale_color_manual(values=c("darkblue", "forestgreen", "darkorange", "red"))

#####################################################################################################
# meta-roi
#####################################################################################################

meas_fn_meta_roi <- "meta-ROI cortical thickness"

#plotting raw data
raw_meta_roi_plot <- ggplot(data = meta_roi_plot_data, aes(x = adjusted_new_time, y = meta_ROI, group = RID, color = as.factor(CDGLOBAL))) +
  theme_classic() + 
  ylab(unique(meas_fn_meta_roi)) + 
  scale_x_continuous(name = "Years from Aβ+", breaks = seq(-10, 10, by = 2)) + #limits=c(-8, 8)) +
  theme(legend.position="none") + 
  geom_point(size = .2) +
  geom_line(size = .2) +
  scale_color_manual(values=c("darkblue", "forestgreen", "darkorange", "red"))

#####################################################################################################
# hippocampal volume
#####################################################################################################

meas_fn_hippocampal_volume <- "Hippocampal Volume"

#plotting raw data
raw_hippocampal_volume_plot <- ggplot(data = hippocampal_volume_plot_data, aes(x = adjusted_new_time, y = hippocampal_volume, group = RID, color = as.factor(CDGLOBAL))) +
  theme_classic() + 
  ylab(unique(meas_fn_hippocampal_volume)) + 
  scale_x_continuous(name = "Years from Aβ+", breaks = seq(-10, 10, by = 2)) + #limits=c(-8, 8)) +
  theme(legend.position="none") + 
  geom_point(size = .2) +
  geom_line(size = .2) +
  scale_color_manual(values=c("darkblue", "forestgreen", "darkorange", "red"))

#####################################################################################################
# Ecog - Subject
#####################################################################################################

meas_fn_ecog_s <- "ECog - Subject"

#plotting raw data
raw_ecog_s_plot <- ggplot(data = ecog_s_plot_data, aes(x = adjusted_new_time, y = EcogGlobal, group = RID, color = as.factor(CDGLOBAL))) +
  theme_classic() + 
  ylab(unique(meas_fn_ecog_s)) + 
  scale_x_continuous(name = "Years from Aβ+", breaks = seq(-10, 10, by = 2)) + #limits=c(-8, 8)) +
  # theme(legend.position="none") + 
  geom_point(size = .2) +
  geom_line(size = .2) +
  scale_color_manual(values=c("darkblue", "forestgreen", "darkorange", "red"))
raw_ecog_s_plot$labels$colour <- "CDR"

#####################################################################################################
# Ecog - Study Partner
#####################################################################################################

meas_fn_ecog_p <- "ECog - Study Partner"

#plotting raw data
raw_ecog_p_plot <- ggplot(data = ecog_p_plot_data, aes(x = adjusted_new_time, y = EcogGlobal, group = RID, color = as.factor(CDGLOBAL))) +
  theme_classic() + 
  ylab(unique(meas_fn_ecog_p)) + 
  scale_x_continuous(name = "Years from Aβ+", breaks = seq(-10, 10, by = 2)) + #limits=c(-8, 8)) +
  theme(legend.position="none") + 
  geom_point(size = .2) +
  geom_line(size = .2) +
  scale_color_manual(values=c("darkblue", "forestgreen", "darkorange", "red"))

#####################################################################################################
# ADAS13
#####################################################################################################

meas_fn_adas13 <- "ADAS-Cog13"

#plotting raw data
raw_adas13_plot <- ggplot(data = adas13_plot_data, aes(x = adjusted_new_time, y = ADAS13, group = RID, color = as.factor(CDGLOBAL))) +
  theme_classic() + 
  ylab(unique(meas_fn_adas13)) + 
  scale_x_continuous(name = "Years from Aβ+", breaks = seq(-10, 10, by = 2)) + #limits=c(-8, 8)) +
  theme(legend.position="none") + 
  geom_point(size = .2) +
  geom_line(size = .2) +
  scale_color_manual(values=c("darkblue", "forestgreen", "darkorange", "red"))

#####################################################################################################
# MMSE
#####################################################################################################

meas_fn_mmse <- "MMSE"

#plotting raw data
raw_mmse_plot <- ggplot(data = mmse_plot_data, aes(x = adjusted_new_time, y = MMSE, group = RID, color = as.factor(CDGLOBAL))) +
  theme_classic() + 
  ylab(unique(meas_fn_mmse)) + 
  scale_x_continuous(name = "Years from Aβ+", breaks = seq(-10, 10, by = 2)) + #limits=c(-8, 8)) +
  theme(legend.position="none") + 
  geom_point(size = .2) +
  geom_line(size = .2) +
  scale_color_manual(values=c("darkblue", "forestgreen", "darkorange", "red"))

#####################################################################################################
# mPACCtrailsB
#####################################################################################################

meas_fn_mpacctrailsb <- "PACC"

#plotting raw data
raw_mpacctrailsb_plot <- ggplot(data = mpacctrailsb_plot_data, aes(x = adjusted_new_time, y = mPACCtrailsB, group = RID, color = as.factor(CDGLOBAL))) +
  theme_classic() + 
  ylab(unique(meas_fn_mpacctrailsb)) + 
  scale_x_continuous(name = "Years from Aβ+", breaks = seq(-10, 10, by = 2)) + #limits=c(-8, 8)) +
  theme(legend.position="none") + 
  geom_point(size = .2) +
  geom_line(size = .2) +
  scale_color_manual(values=c("darkblue", "forestgreen", "darkorange", "red"))

#####################################################################################################
# CDRSB
#####################################################################################################

meas_fn_cdrsb <- "CDR-SB"

#plotting raw data
raw_cdrsb_plot <- ggplot(data = cdrsb_plot_data, aes(x = adjusted_new_time, y = CDRSB, group = RID, color = as.factor(CDGLOBAL))) +
  theme_classic() + 
  ylab(unique(meas_fn_cdrsb)) + 
  scale_x_continuous(name = "Years from Aβ+", breaks = seq(-10, 10, by = 2)) + #limits=c(-8, 8)) +
  theme(legend.position="none") + 
  geom_point(size = .2) +
  geom_line(size = .2) +
  scale_color_manual(values=c("darkblue", "forestgreen", "darkorange", "red"))

#####################################################################################################
# plotting
#####################################################################################################

#raw plots
csf_centiloid_raw <- raw_centiloid_plot + raw_abeta_plot + raw_ptau_plot + raw_tau_plot + plot_layout(nrow = 2)

scans_raw <- raw_hippocampal_volume_plot + raw_meta_roi_plot + raw_fdg_plot + plot_layout(nrow = 1)

fig3 <- (csf_centiloid_raw) / (scans_raw) + plot_annotation(tag_levels = 'A') &  theme(plot.tag = element_text(face = 'bold'))

ragg::agg_tiff("~/plots/S1_plot.jpg", width = 9, height = 8, units = "in", res = 300)
fig3
dev.off()

fig4 <- raw_cdrsb_plot + raw_adas13_plot + raw_mmse_plot + raw_mpacctrailsb_plot + raw_ecog_p_plot + raw_ecog_s_plot + plot_layout(nrow = 2) + plot_annotation(tag_levels = 'A') &  theme(plot.tag = element_text(face = 'bold'))

ragg::agg_tiff("~/plots/S2_plot.jpg", width = 9, height = 5, units = "in", res = 300)
fig4
dev.off()
