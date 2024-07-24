#can also use SUVR 1.42 as a cutoff for amyloid
library(dplyr)
library(ggplot2)
library(segmented)
library(ggrepel)

source("~/Projects/ComGamPackage/ComGamFunctionHelpers.R")
source("~/Projects/ComGamPackage/ComGamHarmFunction.R")

#####################################################################################
#getting the participant ID's that have changed from amyloid negative to amyloid positive
#####################################################################################
amyloid_pet <- readr::read_delim("~/Data/UCBERKELEY_AMY_6MM_30Sep2023.csv") [,1:14]
amyloid_pet <- amyloid_pet %>% dplyr::rename(suvr_summary=SUMMARY_SUVR,
                                             Centiloid=CENTILOIDS,
                                             AmyloidPosPET=AMYLOID_STATUS,
                                             EXAMDATE_pet = SCANDATE) %>%
  dplyr::mutate(EXAMDATE_pet = as.Date(EXAMDATE_pet),
                RID = as.character(RID)) %>%
  dplyr::select(RID, EXAMDATE_pet, AmyloidPosPET, Centiloid, suvr_summary)

change_in_amyloid <- amyloid_pet %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(all(c(0, 1) %in% AmyloidPosPET)) %>%
  dplyr::distinct()

#getting important dates for all ids
last_a_negative_date_pet <- amyloid_pet %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(AmyloidPosPET == 0) %>%
  dplyr::mutate(last_a_neg_date_pet = max(EXAMDATE_pet),
                AmyNeg_Centiloid = Centiloid) %>%
  dplyr::select(RID, Centiloid, EXAMDATE_pet, last_a_neg_date_pet) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()
last_a_negative_date_pet <- last_a_negative_date_pet %>%
  dplyr::filter(EXAMDATE_pet == last_a_neg_date_pet) %>%
  dplyr::mutate(AmyNeg_Centiloid = Centiloid) %>%
  dplyr::select(-Centiloid, -EXAMDATE_pet)

#PET - getting first A+ dates by RID
first_a_positive_date_pet <- amyloid_pet %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(AmyloidPosPET == 1) %>%
  dplyr::mutate(first_a_pos_date_pet = min(EXAMDATE_pet)) %>%
  dplyr::select(RID, Centiloid, EXAMDATE_pet, first_a_pos_date_pet) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()
first_a_positive_date_pet <- first_a_positive_date_pet %>%
  dplyr::filter(EXAMDATE_pet == first_a_pos_date_pet) %>%
  dplyr::mutate(AmyPos_Centiloid = Centiloid) %>%
  dplyr::select(-Centiloid, -EXAMDATE_pet)

change_in_amyloid <- merge(change_in_amyloid, last_a_negative_date_pet, by = "RID")
change_in_amyloid <- merge(change_in_amyloid, first_a_positive_date_pet, by = "RID")

change_in_amyloid_filtered <- change_in_amyloid %>%
  dplyr::filter(EXAMDATE_pet == last_a_neg_date_pet | EXAMDATE_pet == first_a_pos_date_pet)

#####################################################################################
# building a date to frame new x-axis
#####################################################################################

#getting weight-adjusted predicted date
change_in_amyloid_filtered$predicted_date <- as.Date(change_in_amyloid_filtered$last_a_neg_date_pet) + as.numeric(as.Date(change_in_amyloid_filtered$first_a_pos_date_pet) - as.Date(change_in_amyloid_filtered$last_a_neg_date_pet)) * (20 - change_in_amyloid_filtered$AmyNeg_Centiloid) / (change_in_amyloid_filtered$AmyPos_Centiloid -change_in_amyloid_filtered$AmyNeg_Centiloid)

lm_data <- change_in_amyloid_filtered %>%
  dplyr::select(RID, predicted_date) %>%
  dplyr::distinct()

#####################################################################################
# looking at change in amyloid status
#####################################################################################
centiloid_plot_data <- merge(change_in_amyloid, lm_data, all.x = TRUE) %>%
  dplyr::mutate(adjusted_new_time = lubridate::time_length(difftime(EXAMDATE_pet, predicted_date), "years"))

#getting rid of the RID's that have weird Centiloid patterns
wonky_RIDs <- centiloid_plot_data %>%
  filter(adjusted_new_time < 0 & Centiloid > 20)

centiloid_plot_data <- centiloid_plot_data %>%
  dplyr::filter(!(RID %in% unique(wonky_RIDs$RID) | RID == "1261" | RID == "6234"))

#############################################################################
# getting CDGLOBAL, APOE, education, and birth year into a dataset
#############################################################################
#CDGLOBAL
cdglobal_data <- read.csv("~/Data/CDR_24Apr2024.csv")
cdglobal_data$CDR.DATE <- cdglobal_data$USERDATE
cdglobal_data <- cdglobal_data[cdglobal_data$CDGLOBAL>=0,]
cdglobal_data <- cdglobal_data %>%
  dplyr::filter(RID %in% unique(centiloid_plot_data$RID)) %>%
  dplyr::select(RID, CDR.DATE, CDGLOBAL)

#getting APOE type
apoeres <- read.csv("~/Data/APOERES.csv") %>% #OR Local: "C:\\Documents\\paper data longitudinal phases\\APOERES.csv"
  dplyr::select(RID, APGEN1, APGEN2) %>%
  dplyr::mutate(apoe = paste(paste("E", APGEN1, sep = ""), paste("E", APGEN2, sep = ""), sep = "/"),
                RID = as.character(RID)) %>%
  dplyr::distinct()

#getting information to calculate age
dem <- read.csv("~/Data/demographics.csv") %>%
  dplyr::mutate(day = 1,
                birthdate = as.character(paste(day, PTDOB, sep = "/")))

dem <- dem %>%
  dplyr::mutate(birthdate = as.Date(birthdate, format = "%d/%m/%Y")) %>%
  dplyr::select(RID, PTGENDER, PTEDUCAT, birthdate) %>%
  dplyr::distinct() %>%
  dplyr::filter(!is.na(birthdate))

dem <- merge(dem, apoeres, all = TRUE) %>%
  dplyr::mutate(apoe = case_when(apoe == "E2/E3" | apoe == "E2/E2" ~ "E2",
                                 apoe == "E3/E3" ~ "E3",
                                 apoe == "E3/E4" | apoe == "E4/E4" ~ "E4"),
                RID = as.character(RID))

dem_all <- dem

dem <- dem %>%
  dplyr::filter(RID %in% unique(centiloid_plot_data$RID))

#adding diagnoses as CN, MCI, or AD
adni_diagnoses <-
  readr::read_delim("~/Data/DXSUM_PDXCONV_23Jul2024.csv") %>% # change to your file name and location here 
  dplyr::select(RID,DIAGNOSIS,EXAMDATE,VISCODE2) %>%
  dplyr::rename(DX.DATE = EXAMDATE) %>%
  dplyr::mutate(DX = case_when(
    (DIAGNOSIS == 1) ~ "CU",
    (DIAGNOSIS == 2) ~ "MCI",
    (DIAGNOSIS == 3) ~ "Dementia")
  ) %>%
  # dplyr::select(RID,DX,VISCODE2) %>%
  dplyr::filter(!is.na(DX))

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
  dplyr::left_join(lm_data) %>%
  dplyr::mutate(date_diff = abs(as.Date(predicted_date) - as.Date(DX.DATE))) %>%
  dplyr::filter(!is.na(DX.DATE))

diagnoses_at_0 <- diagnoses_at_0 %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(date_diff == min(date_diff)) %>%
  dplyr::ungroup()

#############################################################################
#continuing with Centiloid
#############################################################################

centiloid_plot_data <- as.data.frame(centiloid_plot_data) %>%
  dplyr::left_join(dem, by = "RID") %>%
  dplyr::mutate(age = lubridate::time_length(difftime(EXAMDATE_pet, birthdate), "years"))
centiloid_plot_data <- centiloid_plot_data %>%
  dplyr::mutate(age = round(age, digits = 1))
centiloid_plot_data$RID <- as.numeric(centiloid_plot_data$RID)
centiloid_plot_data <- merge(diagnoses, centiloid_plot_data, all = TRUE) %>%
  dplyr::mutate(time_diff = abs(difftime(EXAMDATE_pet, DX.DATE))) %>%
  dplyr::group_by(RID, EXAMDATE_pet) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup()

centiloid_plot_data <- centiloid_plot_data %>%
  dplyr::select(RID, diags, Centiloid, PTGENDER, PTEDUCAT, apoe, age, adjusted_new_time)

write.csv(centiloid_plot_data, "C:\\Work Folder\\paper data longitudinal phases\\centiloid_lcmm_data.csv")

#####################################################################################
# looking at change in fdg PET Meta-ROI
#####################################################################################
fdg_pet_mean <- read.csv("~/Data/UCBERKELEYFDG_8mm_02_17_23_01Feb2024.csv") %>%
  dplyr::filter(!(ROINAME == "Top50PonsVermis")) %>%
  dplyr::mutate(EXAMDATE = as.Date(EXAMDATE),
                RID = as.character(RID)) %>%
  dplyr::rename(Meta_ROI_fdg_mean = MEAN) %>%
  dplyr::select(-ROINAME)

fdg_pons <- read.csv("~/Data/UCBERKELEYFDG_8mm_02_17_23_01Feb2024.csv") %>%
  dplyr::filter(ROINAME == "Top50PonsVermis") %>%
  dplyr::mutate(mean_pons = MEAN) %>%
  dplyr::select(RID, VISCODE, VISCODE2, EXAMDATE, mean_pons, MAX, STDEV, TOTVOX, update_stamp)

fdg_pet <- merge(fdg_pet_mean, fdg_pons, by = c("RID", "VISCODE", "VISCODE2", "EXAMDATE"), all = TRUE) %>%
  dplyr::select(RID, EXAMDATE, Meta_ROI_fdg_mean, mean_pons) %>%
  dplyr::mutate(adjusted_Meta_ROI = Meta_ROI_fdg_mean/mean_pons)

fdg_plot_data <- fdg_pet %>%
  dplyr::filter(RID %in% unique(centiloid_plot_data$RID)) %>%
  dplyr::left_join(lm_data) %>%
  dplyr::mutate(adjusted_new_time = lubridate::time_length(difftime(EXAMDATE, predicted_date), "years"))

fdg_plot_data <- as.data.frame(fdg_plot_data) %>%
  dplyr::left_join(dem, by = "RID") %>%
  dplyr::mutate(age = lubridate::time_length(difftime(EXAMDATE, birthdate), "years"))
fdg_plot_data <- fdg_plot_data %>%
  dplyr::mutate(age = round(age, digits = 1))
fdg_plot_data$RID <- as.numeric(fdg_plot_data$RID)

fdg_plot_data <- merge(diagnoses, fdg_plot_data, all = TRUE) %>%
  dplyr::mutate(time_diff = abs(difftime(EXAMDATE, DX.DATE))) %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup()

fdg_plot_data <- fdg_plot_data %>%
  dplyr::select(RID, diags, adjusted_Meta_ROI, PTGENDER, PTEDUCAT, apoe, age, adjusted_new_time)

write.csv(fdg_plot_data, "("~/Projects/GAM_models_Amyloid_converted/gam_modelling_data/fdg_gam_data.csv")
