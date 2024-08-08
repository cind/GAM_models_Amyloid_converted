library(dplyr)

######################################################################################################################################
# Script used to get information for table 1
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

#in order to get predicted dates of amyloid onset
lm_data <- read.csv(file.path(root, "lm_data.csv")

########################################
# pulling diagnosis in
adni_diagnoses <- readr::read_delim("~/Data/DXSUM_PDXCONV_23Jul2024.csv") %>% # change to your file name and location here
  dplyr::select(RID,DIAGNOSIS,EXAMDATE,VISCODE2) %>%
  dplyr::rename(DX.DATE = EXAMDATE) %>%
  dplyr::mutate(DX = case_when(
    (DIAGNOSIS == 1) ~ "CU",
    (DIAGNOSIS == 2) ~ "MCI",
    (DIAGNOSIS == 3) ~ "Dementia")
  ) %>%
  dplyr::filter(!is.na(DX))

diagnoses_all <- adni_diagnoses

diagnoses <- adni_diagnoses %>%
  dplyr::rename(diags = DX) %>%
  dplyr::filter(RID %in% centiloid_plot_data$RID,
                !(is.na(DX.DATE))) %>%
  dplyr::mutate(RID = as.character(RID)) %>%
  dplyr::select(RID, DX.DATE, diags)

adni_diagnoses <- adni_diagnoses %>%
  dplyr::select(RID,DX,VISCODE2)

diagnoses_bl <- diagnoses %>%
  dplyr::group_by(RID) %>%
  dplyr::arrange(as.Date(DX.DATE)) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

diagnoses_at_0 <- diagnoses %>%
  dplyr::left_join(lm_data %>%
                     mutate(RID = as.character(RID))) %>%
  dplyr::mutate(date_diff = abs(as.Date(predicted_date) - as.Date(DX.DATE))) %>%
  dplyr::filter(!is.na(DX.DATE),
                RID %in% unique(centiloid_plot_data$RID))

diagnoses_at_0 <- diagnoses_at_0 %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(date_diff == min(date_diff)) %>%
  dplyr::filter(row_number() == 1) %>% #4485 has the same diagnosis an equal number of days before and after
  dplyr::ungroup() %>%
  dplyr::rename(diag_closest_0 = diags,
                diag_date_closest_0 = DX.DATE) %>%
  dplyr::select(RID, diag_date_closest_0, diag_closest_0)

################################################################################
# Demographics
################################################################################

# getting APOE type
apoeres <- read.csv("~/Data/APOERES.csv") %>% #OR Local: "C:\\Documents\\paper data longitudinal phases\\APOERES.csv"
  dplyr::select(RID, APGEN1, APGEN2) %>%
  dplyr::mutate(apoe = paste(paste("E", APGEN1, sep = ""), paste("E", APGEN2, sep = ""), sep = "/"),
                RID = as.character(RID)) %>%
  dplyr::distinct()

# getting information to calculate age
dem <- read.csv("~/Data/demographics.csv") %>%
  dplyr::mutate(day = 1,
                birthdate = as.character(paste(day, PTDOB, sep = "/")))

dem <- dem %>%
  dplyr::mutate(birthdate = as.Date(birthdate, format = "%d/%m/%Y")) %>%
  dplyr::select(RID, PTGENDER, PTEDUCAT, birthdate, PTRACCAT, PTETHCAT) %>%
  dplyr::distinct() %>%
  dplyr::filter(!is.na(birthdate))

dem <- merge(dem, apoeres, all = TRUE) %>%
  dplyr::mutate(apoe = case_when(apoe == "E2/E3" | apoe == "E2/E2" ~ "E2",
                                 apoe == "E3/E3" ~ "E3",
                                 apoe == "E3/E4" | apoe == "E4/E4" ~ "E4"),
                RID = as.character(RID))

dem_all <- dem

dem_all <- dem_all %>%
  dplyr::group_by(RID) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

dem <- dem %>%
  dplyr::filter(RID %in% unique(centiloid_plot_data$RID))

############################# demographics mean, sd #############################
# education mean, sd
education_table <- merge(dem, diagnoses_at_0) %>%
  dplyr::group_by(diag_closest_0) %>%
  dplyr::summarise(mean=mean(PTEDUCAT), sd=sd(PTEDUCAT), n = n(), range = paste(min(PTEDUCAT), max(PTEDUCAT)))

# age mean, sd
age_table <- merge(dem, diagnoses_at_0) %>%
  dplyr::left_join(lm_data %>%
                     dplyr::mutate(RID = as.character(RID))) %>%
  dplyr::mutate(age = lubridate::time_length(difftime(predicted_date, birthdate), "years")) %>%
  dplyr::group_by(diag_closest_0) %>%
  dplyr::summarise(mean=mean(age), sd=sd(age), n = n(), range = paste(min(age), max(age)))

# apoe percentages
merge(dem, diagnoses_at_0) %>%
  dplyr::group_by(diag_closest_0, paste(APGEN1, APGEN2)) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))

# gender percentages
merge(dem, diagnoses_at_0) %>%
  dplyr::group_by(diag_closest_0, PTGENDER) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))

# race percentages
merge(dem, diagnoses_at_0) %>%
  dplyr::group_by(diag_closest_0, PTRACCAT) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))

# ethnicity percentages
merge(dem, diagnoses_at_0) %>%
  dplyr::group_by(diag_closest_0, PTETHCAT) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))
  
################################################################################
# Centiloid
################################################################################

amyloid_pet <- readr::read_delim("~/Data/UCBERKELEY_AMY_6MM_31Jul2024.csv") [,1:16]
amyloid_pet <- amyloid_pet %>% dplyr::rename(suvr_summary=SUMMARY_SUVR,
                                             Centiloid=CENTILOIDS,
                                             AmyloidPosPET=AMYLOID_STATUS,
                                             EXAMDATE_pet = SCANDATE) %>%
  dplyr::mutate(EXAMDATE_pet = as.Date(EXAMDATE_pet),
                RID = as.character(RID)) %>%
  dplyr::select(RID, EXAMDATE_pet, AmyloidPosPET, Centiloid, suvr_summary)

amyloid_longitudinal <- amyloid_pet %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(n()>1) %>%
  dplyr::ungroup()

change_in_amyloid <- amyloid_longitudinal %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(all(c(0, 1) %in% AmyloidPosPET)) %>%
  dplyr::distinct()

# getting important dates for all ids
last_a_negative_date_pet <- amyloid_pet %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(AmyloidPosPET == 0) %>%
  dplyr::mutate(last_a_neg_date_pet = max(EXAMDATE_pet),
                AmyNeg_Centiloid = Centiloid) %>%
  dplyr::select(RID, Centiloid, EXAMDATE_pet, last_a_neg_date_pet) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()
last_a_negative_date_pet <- last_a_negative_date_pet %>%
  dplyr::filter(EXAMDATE_pet == last_a_neg_date_pet,
                RID %in% unique(centiloid_plot_data$RID)) %>%
  dplyr::mutate(AmyNeg_Centiloid = Centiloid) %>%
  dplyr::select(-Centiloid, -EXAMDATE_pet)

# PET - getting first A+ dates by RID
first_a_positive_date_pet <- amyloid_pet %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(AmyloidPosPET == 1) %>%
  dplyr::mutate(first_a_pos_date_pet = min(EXAMDATE_pet)) %>%
  dplyr::select(RID, Centiloid, EXAMDATE_pet, first_a_pos_date_pet) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()
first_a_positive_date_pet <- first_a_positive_date_pet %>%
  dplyr::filter(EXAMDATE_pet == first_a_pos_date_pet,
                RID %in% unique(centiloid_plot_data$RID)) %>%
  dplyr::mutate(AmyPos_Centiloid = Centiloid) %>%
  dplyr::select(-Centiloid, -EXAMDATE_pet)

########################################
# first getting diagnoses closest to dates

last_a_negative_date_pet <- merge(diagnoses, last_a_negative_date_pet, all = TRUE) %>%
  dplyr::mutate(time_diff = abs(difftime(last_a_neg_date_pet, DX.DATE))) %>%
  dplyr::group_by(RID, last_a_neg_date_pet) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup() %>%
  dplyr::rename(closest_diag = diags, 
                closest_diag_date = DX.DATE)

first_a_positive_date_pet <- merge(diagnoses, first_a_positive_date_pet, all = TRUE) %>%
  dplyr::mutate(time_diff = abs(difftime(first_a_pos_date_pet, DX.DATE))) %>%
  dplyr::group_by(RID, first_a_pos_date_pet) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup() %>%
  dplyr::rename(closest_diag = diags,
                closest_diag_date = DX.DATE)

length(unique(last_a_negative_date_pet$RID))
length(unique(first_a_positive_date_pet$RID))

########################################
# getting mean, sd, range of Centiloid at last visit and first visit

table_info_centiloid_before <- merge(diagnoses_at_0, last_a_negative_date_pet) %>%
    group_by(diag_closest_0) %>%
  dplyr::summarise( mean=mean(AmyNeg_Centiloid), sd=sd(AmyNeg_Centiloid), n = n(), range = paste(min(AmyNeg_Centiloid), max(AmyNeg_Centiloid)))

table_info_centiloid_after <- merge(diagnoses_at_0, first_a_positive_date_pet) %>%
  group_by(diag_closest_0) %>%
  dplyr::summarise(mean=mean(AmyPos_Centiloid), sd=sd(AmyPos_Centiloid), n = n(), range = paste(min(AmyPos_Centiloid), max(AmyPos_Centiloid)))

########################################
# getting time between estimated amyloid onset and diagnosis (reviewer's comment)
dates_between <- merge(lm_data, diagnoses_at_0) %>%
  dplyr::mutate(time_diff_years = abs(lubridate::time_length(difftime(predicted_date, diag_date_closest_0), "years")),
                time_diff_days = abs(lubridate::time_length(difftime(predicted_date, diag_date_closest_0), "days")))

dates_between %>%
  dplyr::summarise(mean=mean(time_diff_years), sd=sd(time_diff_years), range = paste(min(time_diff_years), max(time_diff_years)))
dates_between %>%
  dplyr::summarise(mean=mean(time_diff_days), sd=sd(time_diff_days), range = paste(min(time_diff_days), max(time_diff_days)))

################################################################################
# Centiloid
################################################################################

# getting centiloid mean years before onset/after onset (reviewer's comment)                    
centiloid_years_table_data <- centiloid_plot_data %>%
  dplyr::group_by(RID) %>%
  dplyr::mutate(years_pre = abs(min(adjusted_new_time)),
                count = n(),
                years_post = max(adjusted_new_time)) %>%
  dplyr::ungroup() %>%
  dplyr::select(RID, years_pre, years_post, count) %>%
  dplyr::distinct()

centiloid_years_table <- merge(centiloid_years_table_data, diagnoses_at_0) %>%
  dplyr::group_by(diag_closest_0) %>%
  dplyr::summarise(mean_years_pre=mean(years_pre), sd_years_pre=sd(years_pre), range_years_pre = paste(min(years_pre), max(years_pre)), mean_years_post=mean(years_post), sd_years_post=sd(years_post), range_years_post = paste(min(years_post), max(years_post)), mean_count=mean(count), sd_count=sd(count), range_count = paste(min(count), max(count)), n = n())

################################################################################
# CSF ABETA
################################################################################

# getting mean, sd, range of ABETA                    
abeta_table_data <- abeta_plot_data %>%
  dplyr::mutate(time_from_conversion = abs(adjusted_new_time)) %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(time_from_conversion == min(time_from_conversion)) %>%
  dplyr::ungroup()

abeta_table <- merge(abeta_table_data, diagnoses_at_0) %>%
  dplyr::group_by(diag_closest_0) %>%
  dplyr::summarise(mean=mean(ABETA), sd=sd(ABETA), n = n(), range = paste(min(ABETA), max(ABETA)))

################################################################################
# CSF TAU
################################################################################

# getting mean, sd, range of TAU                                        
tau_table_data <- tau_plot_data %>%
  dplyr::mutate(time_from_conversion = abs(adjusted_new_time)) %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(time_from_conversion == min(time_from_conversion)) %>%
  dplyr::ungroup()

tau_table <- merge(tau_table_data, diagnoses_at_0) %>%
  dplyr::group_by(diag_closest_0) %>%
  dplyr::summarise(mean=mean(TAU), sd=sd(TAU), n = n(), range = paste(min(TAU), max(TAU)))

################################################################################
# CSF PTAU
################################################################################

# getting mean, sd, range of PTAU                    
ptau_table_data <- ptau_plot_data %>%
  dplyr::mutate(time_from_conversion = abs(adjusted_new_time)) %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(time_from_conversion == min(time_from_conversion)) %>%
  dplyr::ungroup()

ptau_table <- merge(ptau_table_data, diagnoses_at_0) %>%
  dplyr::group_by(diag_closest_0) %>%
  dplyr::summarise(mean=mean(PTAU), sd=sd(PTAU), n = n(), range = paste(min(PTAU), max(PTAU)))

################################################################################
# Hippocampal Volume
################################################################################

# getting mean, sd, range of hippocampal volume                    
hippocampal_volume_table_data <- hippocampal_plot_data %>%
  dplyr::mutate(time_from_conversion = abs(adjusted_new_time)) %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(time_from_conversion == min(time_from_conversion)) %>%
  dplyr::ungroup()

hippocampal_volume_table <- merge(hippocampal_volume_table_data, diagnoses_at_0) %>%
  dplyr::group_by(diag_closest_0) %>%
  dplyr::summarise(mean=mean(hippocampal_volume), sd=sd(hippocampal_volume), n = n(), range = paste(min(hippocampal_volume), max(hippocampal_volume)))

# getting number of years/visits before/after AB onset
hippocampal_volume_years_table_data <- hippocampal_plot_data %>%
  dplyr::group_by(RID) %>%
  dplyr::mutate(years_pre = abs(min(adjusted_new_time)),
                count = n(),
                years_post = max(adjusted_new_time)) %>%
  dplyr::ungroup() %>%
  dplyr::select(RID, years_pre, years_post, count) %>%
  dplyr::distinct()

hippocampal_volume_years_table <- merge(hippocampal_volume_years_table_data, diagnoses_at_0) %>%
  dplyr::group_by(diag_closest_0) %>%
  dplyr::summarise(mean_years_pre=mean(years_pre), sd_years_pre=sd(years_pre), range_years_pre = paste(min(years_pre), max(years_pre)), mean_years_post=mean(years_post), sd_years_post=sd(years_post), range_years_post = paste(min(years_post), max(years_post)), mean_count=mean(count), sd_count=sd(count), range_count = paste(min(count), max(count)), n = n())

################################################################################
# Meta-ROI Cortical Thickness
################################################################################

# getting mean, sd, range of meta-roi cortical thickness                    
meta_roi_table_data <- meta_roi_plot_data %>%
  dplyr::mutate(time_from_conversion = abs(adjusted_new_time)) %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(time_from_conversion == min(time_from_conversion)) %>%
  dplyr::ungroup()

meta_roi_table <- merge(meta_roi_table_data, diagnoses_at_0) %>%
  dplyr::group_by(diag_closest_0) %>%
  dplyr::summarise(mean=mean(meta_ROI), sd=sd(meta_ROI), n = n(), range = paste(min(meta_ROI), max(meta_ROI)))

################################################################################
# FDG PET
################################################################################

# getting mean, sd, range of FDG PET data                    
fdg_table_data <- fdg_plot_data %>%
  dplyr::mutate(time_from_conversion = abs(adjusted_new_time)) %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(time_from_conversion == min(time_from_conversion)) %>%
  dplyr::ungroup()

fdg_table <- merge(fdg_table_data, diagnoses_at_0) %>%
  dplyr::group_by(diag_closest_0) %>%
  dplyr::summarise(mean=mean(adjusted_Meta_ROI), sd=sd(adjusted_Meta_ROI), n = n(), range = paste(min(adjusted_Meta_ROI), max(adjusted_Meta_ROI)))

################################################################################
# mPACCtrailsB
################################################################################

# getting mean, sd, range of PACC                    
mpacctrailsb_table_data <- mpacctrailsb_plot_data %>%
  dplyr::mutate(time_from_conversion = abs(adjusted_new_time)) %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(time_from_conversion == min(time_from_conversion)) %>%
  dplyr::ungroup()

mpacctrailsb_table <- merge(mpacctrailsb_table_data, diagnoses_at_0) %>%
  dplyr::group_by(diag_closest_0) %>%
  dplyr::summarise(mean=mean(mPACCtrailsB), sd=sd(mPACCtrailsB), n = n(), range = paste(min(mPACCtrailsB), max(mPACCtrailsB)))

# getting number of years/visits before/after AB onset
mpacctrailsb_years_table_data <- mpacctrailsb_plot_data %>%
  dplyr::group_by(RID) %>%
  dplyr::mutate(years_pre = abs(min(adjusted_new_time)),
                years_post = max(adjusted_new_time),
                count_pre = sum(adjusted_new_time > 0),
                count_post = sum(adjusted_new_time < 0),
                count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::select(RID, years_pre, years_post, count_pre, count_post, count) %>%
  dplyr::distinct()

mpacctrailsb_years_table <- merge(mpacctrailsb_years_table_data, diagnoses_at_0) %>%
  dplyr::group_by(diag_closest_0) %>%
  dplyr::summarise(mean_years_pre=mean(years_pre), sd_years_pre=sd(years_pre), range_years_pre = paste(min(years_pre), max(years_pre)), mean_years_post=mean(years_post), sd_years_post=sd(years_post), range_years_post = paste(min(years_post), max(years_post)), mean_count_pre=mean(count_pre), sd_count_pre=sd(count_pre), range_count_pre = paste(min(count_pre), max(count_pre)), mean_count_post=mean(count_post), sd_count_post=sd(count_post), range_count_post = paste(min(count_post), max(count_post)), n = n())

################################################################################
# adas13
################################################################################

# getting mean, sd, range of ADAS13                    
adas13_table_data <- adas13_plot_data %>%
  dplyr::mutate(time_from_conversion = abs(adjusted_new_time)) %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(time_from_conversion == min(time_from_conversion)) %>%
  dplyr::ungroup()

adas13_table <- merge(adas13_table_data, diagnoses_at_0) %>%
  dplyr::group_by(diag_closest_0) %>%
  dplyr::summarise(mean=mean(ADAS13), sd=sd(ADAS13), n = n(), range = paste(min(ADAS13), max(ADAS13)))

# getting number of years/visits before/after AB onset
adas13_years_table_data <- adas13_plot_data %>%
  dplyr::group_by(RID) %>%
  dplyr::mutate(years_pre = abs(min(adjusted_new_time)),
                years_post = max(adjusted_new_time),
                count_pre = sum(adjusted_new_time > 0),
                count_post = sum(adjusted_new_time < 0),
                count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::select(RID, years_pre, years_post, count_pre, count_post, count) %>%
  dplyr::distinct()

adas13_years_table <- merge(adas13_years_table_data, diagnoses_at_0) %>%
  dplyr::group_by(diag_closest_0) %>%
  dplyr::summarise(mean_years_pre=mean(years_pre), sd_years_pre=sd(years_pre), range_years_pre = paste(min(years_pre), max(years_pre)), mean_years_post=mean(years_post), sd_years_post=sd(years_post), range_years_post = paste(min(years_post), max(years_post)), mean_count_pre=mean(count_pre), sd_count_pre=sd(count_pre), range_count_pre = paste(min(count_pre), max(count_pre)), mean_count_post=mean(count_post), sd_count_post=sd(count_post), range_count_post = paste(min(count_post), max(count_post)), n = n())

################################################################################
# cdrsb - this used for table 1
################################################################################

# getting mean, sd, range of CDRSB                    
cdrsb_table_data <- cdrsb_plot_data %>%
  dplyr::mutate(time_from_conversion = abs(adjusted_new_time)) %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(time_from_conversion == min(time_from_conversion)) %>%
  dplyr::ungroup()

cdrsb_table <- merge(cdrsb_table_data, diagnoses_at_0) %>%
  dplyr::group_by(diag_closest_0) %>%
  dplyr::summarise(mean=mean(CDRSB), sd=sd(CDRSB), n = n(), range = paste(min(CDRSB), max(CDRSB)))

# getting number of years/visits before/after AB onset (ended up using info from CDRSB for Table 1)
cdrsb_years_table_data <- cdrsb_plot_data %>%
  dplyr::group_by(RID) %>%
  dplyr::mutate(years_pre = abs(min(adjusted_new_time)),
                years_post = max(adjusted_new_time),
                count_pre = sum(adjusted_new_time > 0),
                count_post = sum(adjusted_new_time < 0),
                count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::select(RID, years_pre, years_post, count_pre, count_post, count) %>%
  dplyr::distinct()

cdrsb_years_table <- merge(cdrsb_years_table_data, diagnoses_at_0) %>%
  dplyr::group_by(diag_closest_0) %>%
  dplyr::summarise(mean_years_pre=mean(years_pre), sd_years_pre=sd(years_pre), range_years_pre = paste(min(years_pre), max(years_pre)), mean_years_post=mean(years_post), sd_years_post=sd(years_post), range_years_post = paste(min(years_post), max(years_post)), mean_count_pre=mean(count_pre), sd_count_pre=sd(count_pre), range_count_pre = paste(min(count_pre), max(count_pre)), mean_count_post=mean(count_post), sd_count_post=sd(count_post), range_count_post = paste(min(count_post), max(count_post)), n = n())


################################################################################
# mmse
################################################################################

# getting mean, sd, range of MMSE                   
mmse_table_data <- mmse_plot_data %>%
  dplyr::mutate(time_from_conversion = abs(adjusted_new_time)) %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(time_from_conversion == min(time_from_conversion)) %>%
  dplyr::ungroup()

mmse_table <- merge(mmse_table_data, diagnoses_at_0) %>%
  dplyr::group_by(diag_closest_0) %>%
  dplyr::summarise(mean=mean(MMSE), sd=sd(MMSE), n = n(), range = paste(min(MMSE), max(MMSE)))

# getting number of years/visits before/after AB onset
mmse_years_table_data <- mmse_plot_data %>%
  dplyr::group_by(RID) %>%
  dplyr::mutate(years_pre = abs(min(adjusted_new_time)),
                years_post = max(adjusted_new_time),
                count_pre = sum(adjusted_new_time > 0),
                count_post = sum(adjusted_new_time < 0),
                count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::select(RID, years_pre, years_post, count_pre, count_post, count) %>%
  dplyr::distinct()

mmse_years_table <- merge(mmse_years_table_data, diagnoses_at_0) %>%
  dplyr::group_by(diag_closest_0) %>%
  dplyr::summarise(mean_years_pre=mean(years_pre), sd_years_pre=sd(years_pre), range_years_pre = paste(min(years_pre), max(years_pre)), mean_years_post=mean(years_post), sd_years_post=sd(years_post), range_years_post = paste(min(years_post), max(years_post)), mean_count_pre=mean(count_pre), sd_count_pre=sd(count_pre), range_count_pre = paste(min(count_pre), max(count_pre)), mean_count_post=mean(count_post), sd_count_post=sd(count_post), range_count_post = paste(min(count_post), max(count_post)), n = n())

################################################################################
# ECog - Study Partner
################################################################################

# getting mean, sd, range of ECog - Study Partner                    
ecog_p_table_data <- ecog_p_plot_data %>%
  dplyr::mutate(time_from_conversion = abs(adjusted_new_time)) %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(time_from_conversion == min(time_from_conversion)) %>%
  dplyr::ungroup()

ecog_p_table <- merge(ecog_p_table_data, diagnoses_at_0) %>%
  dplyr::group_by(diag_closest_0) %>%
  dplyr::summarise(mean=mean(EcogGlobal), sd=sd(EcogGlobal), n = n(), range = paste(min(EcogGlobal), max(EcogGlobal)))

# getting number of years/visits before/after AB onset
ecog_p_years_table_data <- ecog_p_plot_data %>%
  dplyr::group_by(RID) %>%
  dplyr::mutate(years_pre = abs(min(adjusted_new_time)),
                count = n(),
                years_post = max(adjusted_new_time)) %>%
  dplyr::ungroup() %>%
  dplyr::select(RID, years_pre, years_post, count) %>%
  dplyr::distinct()

ecog_p_years_table <- merge(ecog_p_years_table_data, diagnoses_at_0) %>%
  dplyr::group_by(diag_closest_0) %>%
  dplyr::summarise(mean_years_pre=mean(years_pre), sd_years_pre=sd(years_pre), range_years_pre = paste(min(years_pre), max(years_pre)), mean_years_post=mean(years_post), sd_years_post=sd(years_post), range_years_post = paste(min(years_post), max(years_post)), mean_count=mean(count), sd_count=sd(count), range_count = paste(min(count), max(count)), n = n())

################################################################################
# ECog - Subject
################################################################################

# getting mean, sd, range of ECog - Subject                   
ecog_s_table_data <- ecog_s_plot_data %>%
  dplyr::mutate(time_from_conversion = abs(adjusted_new_time)) %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(time_from_conversion == min(time_from_conversion)) %>%
  dplyr::ungroup()

ecog_s_table <- merge(ecog_s_table_data, diagnoses_at_0) %>%
  dplyr::group_by(diag_closest_0) %>%
  dplyr::summarise(mean=mean(EcogGlobal), sd=sd(EcogGlobal), n = n(), range = paste(min(EcogGlobal), max(EcogGlobal)))

################################################################################
#percentages for paper population description
################################################################################

# getting percentages of diagnoses per gender
dem %>%
  dplyr::group_by(PTGENDER) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))

merge(dem, diagnoses_at_0) %>%
  dplyr::group_by(diag_closest_0, PTGENDER) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(freq = n / sum(n))

# getting percentages of diagnoses per apoe
dem %>%
  dplyr::group_by(paste(APGEN1, APGEN2)) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))

merge(dem, diagnoses_at_0) %>%
  dplyr::group_by(diag_closest_0, paste(APGEN1, APGEN2)) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(freq = n / sum(n)) %>%
  dplyr::select(-n)

# getting percentages of diagnoses per race
dem %>%
  dplyr::mutate(race = case_when(PTRACCAT == "5" ~ "white",
                                 TRUE ~ "black/mixed")) %>%
  dplyr::group_by(race) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))

merge(dem, diagnoses_at_0) %>%
  dplyr::mutate(race = case_when(PTRACCAT == "5" ~ "white",
                                 TRUE ~ "black/mixed")) %>%
  dplyr::group_by(diag_closest_0, race) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(freq = n / sum(n)) %>%
  dplyr::select(-n)

# getting percentages of diagnoses per ethnicity
dem %>%
  dplyr::group_by(PTETHCAT) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))

merge(dem, diagnoses_at_0) %>%
  dplyr::group_by(diag_closest_0, PTETHCAT) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(freq = n / sum(n)) %>%
  dplyr::select(-n)

##############                    
# getting percentages of data availability

# MRI
length(unique(hippocampal_volume_plot_data$RID))/length(unique(centiloid_plot_data$RID))

# CSF
length(unique(ptau_plot_data$RID))/length(unique(centiloid_plot_data$RID))

# FDG PET
length(unique(fdg_plot_data$RID))/length(unique(centiloid_plot_data$RID))
