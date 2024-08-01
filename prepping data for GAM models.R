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
  dplyr::filter(!(RID %in% unique(wonky_RIDs$RID) | RID == "1261" | RID == "6234" | RID == "6454")) #6454 is amyloid positive based on SUVR, not centiloid. the others shift from positive to negative

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
apoeres <- read.csv("~/Data/APOERES.csv") %>% #OR Local: "~/Data/APOERES.csv"
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

write.csv(centiloid_plot_data, "~/Projects/GAM_models_Amyloid_converted/gam_modelling_data/centiloid_lcmm_data.csv")

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

write.csv(fdg_plot_data, "~/Projects/GAM_models_Amyloid_converted/gam_modelling_data/fdg_gam_data.csv")

#############################################################################
# MRI: Meta-ROI
#############################################################################

#############################
# Loading in data
#############################
# ADNIGO/2 data
adni2_3 <- read.csv("~/Data/freesurfer_5.1_ADNI2GO_protocol_3T.csv") %>%
  dplyr::mutate(STUDY = "ADNI2",
                Field_Strength = "3T")
adni2_3 <- adni2_3 %>%
  dplyr::filter(OVERALLQC == "Pass")

#getting demographics and closest diagnosis into the data
adni2_3 <- merge(adni2_3, dem_all, all.x = TRUE) %>% 
  dplyr::mutate(age = round(lubridate::time_length(difftime(EXAMDATE, birthdate), "years"), digits = 1))
adni2_3$RID <- as.factor(adni2_3$RID)
adni_diagnoses$RID <- as.factor(adni_diagnoses$RID)

adni2_3 <- merge(adni2_3, diagnoses_all %>% 
                        dplyr::select(-VISCODE2) %>%
                        dplyr::filter(!is.na(DX.DATE)), 
                      by = "RID") %>%
  dplyr::mutate(date_diff = abs(lubridate::time_length(difftime(EXAMDATE, DX.DATE), "years"))) %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::mutate(min_match_date = min(date_diff)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(min_match_date == date_diff) %>%
  dplyr::distinct() %>%
  dplyr::rename(diags = DX)

adni2_3 <- adni2_3 %>%
  dplyr::rename(ICV = ST10CV)

# ADNI3 data
adni3_3 <- read.csv("~/Data/freesurfer_6.0_ADNI3_protocol_3T.csv") %>%
  dplyr::mutate(STUDY = "ADNI3",
                Field_Strength = "3T")
freesurfer_data <- read.csv("~/Data/adni3_freesurfer.csv") %>% # received this file from Mark - files on LONI were not updated to include recent scans
  dplyr::rename(EXAMDATE = SCANDATE) %>%
  dplyr::mutate(Field_Strength = "3T")

adni_3 <- merge(adni3_3, freesurfer_data, all = TRUE)

adni_3 <- adni_3 %>%
  dplyr::filter(OVERALLQC == "Pass")

#getting demographics and closest diagnosis into the data
adni_3 <- merge(adni_3, dem_all, all.x = TRUE) %>% 
  dplyr::mutate(age = round(lubridate::time_length(difftime(EXAMDATE, birthdate), "years"), digits = 1))
adni_3$RID <- as.factor(adni_3$RID)
adni_diagnoses$RID <- as.factor(adni_diagnoses$RID)

adni_3 <- merge(adni_3, diagnoses_all %>% 
                   dplyr::select(-VISCODE2) %>%
                   dplyr::filter(!is.na(DX.DATE)), 
                 by = "RID") %>%
  dplyr::mutate(date_diff = abs(lubridate::time_length(difftime(EXAMDATE, DX.DATE), "years"))) %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::mutate(min_match_date = min(date_diff)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(min_match_date == date_diff) %>%
  dplyr::distinct() %>%
  dplyr::rename(diags = DX)

adni_3 <- adni_3 %>%
  dplyr::rename(ICV = ST10CV)
  
#############################
# ICV-adjustments: getting CU, amyloid-negative scans
#############################
#getting amyloid negative RID's that were consistently negative for every scan they had
amyloid_negs <- amyloid_pet %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(AmyloidPosPET == 0) %>%
  dplyr::mutate(first_a_neg_date_pet = min(EXAMDATE_pet),
                RID = as.character(RID)) %>%
  dplyr::select(RID, first_a_neg_date_pet) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()

#getting earliest, amyloid-negative scan per RID
adni2_3_temp <- adni2_3 %>%
  dplyr::filter(diags == "CU")
adni2_3_temp <- adni2_3_temp %>%
  dplyr::left_join(amyloid_negs) %>%
  dplyr::mutate(diff_time = abs(lubridate::time_length(difftime(EXAMDATE, first_a_neg_date_pet), "years"))) %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(diff_time == min(diff_time)) %>%
  dplyr::ungroup()
adni2_3_temp <- adni2_3_temp %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup()

adni_3_temp <- adni_3 %>%
  dplyr::filter(diags == "CU")
adni_3_temp <- adni_3_temp %>%
  dplyr::left_join(amyloid_negs) %>%
  dplyr::mutate(diff_time = abs(lubridate::time_length(difftime(EXAMDATE, first_a_neg_date_pet), "years"))) %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(diff_time == min(diff_time)) %>%
  dplyr::ungroup()
adni_3_temp <- adni_3_temp %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup()

#############################
# ICV-adjustments: setting up regression
#############################

#ADNIGO/2
# creating list of variables that need to be ICV-adjusted
volumes_adni2_3 <- adni2_3_temp %>%
  dplyr::select(ICV | contains("CV") | contains("SV"))

variables_adni2_3 <- c()

for (roi_num in 1:length(names(volumes_adni2_3))){
  current_roi <- names(volumes_adni2_3)[roi_num]
  variables_adni2_3 <- append(variables_adni2_3, current_roi)
}

for (i in variables_adni2_3[-1]){
    fit <- lm(adni2_3_temp[[i]] ~ ICV, data = adni2_3_temp, na.action=na.exclude)
    adni2_3[paste0(i,'_adjusted')] <- adni2_3[paste0(i)] - (adni2_3$ICV - mean(adni2_3_temp$ICV, na.rm = T)) * fit$coefficients[2]
}

adni2_3 <- adni2_3 %>%
  dplyr::select(-names(volumes_adni2_3))

#ADNI3
# creating list of variables that need to be ICV-adjusted
volumes_adni_3 <- adni_3_temp %>%
  dplyr::select(ICV | contains("CV") | contains("SV"))

variables_adni_3 <- c()

for (roi_num in 1:length(names(volumes_adni_3))){
  current_roi <- names(volumes_adni_3)[roi_num]
  variables_adni_3 <- append(variables_adni_3, current_roi)
}

for (i in variables_adni_3[-1]){
  fit <- lm(adni_3_temp[[i]] ~ ICV, data = adni_3_temp, na.action=na.exclude)
  adni_3[paste0(i,'_adjusted')] <- adni_3[paste0(i)] - (adni_3$ICV - mean(adni_3_temp$ICV, na.rm = T)) * fit$coefficients[2]
}

adni_3 <- adni_3 %>%
  dplyr::select(-names(volumes_adni_3))

#############################
# merging variables together
#############################

mri_data <- merge(adni2_3, adni_3, all = TRUE)

#############################
# ComBat Harmonization
#############################
# getting mri data split into features and covariates
mri_data <- mri_data %>%
  tidyr::drop_na(ST125SV_adjusted | ST66SV_adjusted | ST147SV_adjusted | ST148SV_adjusted | ST149SV_adjusted | ST150SV_adjusted | ST151SV_adjusted | ST152SV_adjusted | ST153SV_adjusted | ST154SV_adjusted | ST155SV_adjusted) #these variables had NA values, and combat wont work with NA

all_features <- mri_data[, c(20:225, 241:358, 363:378)] 
all_features <- all_features %>%
  dplyr::select(-contains("TS")) # has NA values and we do not need the variables

all_features <- all_features %>%
  dplyr::select_if(~mean(is.na(.)) < 0.2)
colSums(is.na(all_features))

all_covariates <- mri_data[, c(1:19, 226:240, 359:362, 379:381)] 

# now creating dataset for harmonization
all_data_combined_CN <- cbind(all_features, all_covariates) %>%
  dplyr::filter(diags == "CU")

all_features_CN <- all_data_combined_CN[, 1:254] 
all_covariates_CN <- all_data_combined_CN[, c(275, 277, 283)]  
extra_covariates_CN <- all_data_combined_CN[, c(255:273, 276, 278:282, 286, 289:293, 295)]
all_covariates_CN$STUDY <- as.factor(all_covariates_CN$STUDY)
all_covariates_CN$age <- as.numeric(all_covariates_CN$age)
all_covariates_CN$PTGENDER <- as.factor(all_covariates_CN$PTGENDER)

CN_data_harmonized <- ComGamHarm(feature.data = all_features_CN,
                                 covar.data = all_covariates_CN, 
                                 eb                = TRUE,
                                 parametric        = TRUE,
                                 smooth.terms      = c("age"),  # c("Examdate_Age"),
                                 k.val             = 5,
                                 verbose           = TRUE,
                                 model.diagnostics = FALSE)

CN_harmonized_data <- as.data.frame(t(CN_data_harmonized$harm.results))

CN_harmonized_data <- cbind(extra_covariates_CN, all_covariates_CN, CN_harmonized_data)

all_data_combined_MCI_AD <- cbind(all_features, all_covariates) %>%
  dplyr::filter(!diags == "CN")

all_features_MCI_AD <- all_data_combined_MCI_AD[, 1:254] 
all_covariates_MCI_AD <- all_data_combined_MCI_AD[, c(275, 277, 283)] 
extra_covariates_MCI_AD <- all_data_combined_MCI_AD[, c(255:273, 276, 278:282, 286, 289:293, 295)]
all_covariates_MCI_AD$STUDY <- as.factor(all_covariates_MCI_AD$STUDY)
all_covariates_MCI_AD$age <- as.numeric(all_covariates_MCI_AD$age)
all_covariates_MCI_AD$PTGENDER <- as.factor(all_covariates_MCI_AD$PTGENDER)

MCI_AD_harmonized <- t(ApplyHarm(feature.data   = all_features_MCI_AD,
                                 covariate.data = all_covariates_MCI_AD,
                                 comgam.out     = CN_data_harmonized))

MCI_AD_harmonized <- as.data.frame(MCI_AD_harmonized)

MCI_AD_harmonized <- cbind(extra_covariates_MCI_AD, all_covariates_MCI_AD, MCI_AD_harmonized)

harmonized_data_from_CN <- rbind(CN_harmonized_data, MCI_AD_harmonized) #%>%

#############################
# Age Adjustment
#############################
mri_plot_data <- harmonized_data_from_CN %>%
  dplyr::distinct()

#getting amyloid-negative CN cases
mri_plot_data_temp <- mri_plot_data %>%
  dplyr::filter(diags == "CU")
mri_plot_data_temp <- mri_plot_data_temp %>%
  dplyr::left_join(amyloid_negs) %>%
  dplyr::mutate(diff_time = abs(lubridate::time_length(difftime(EXAMDATE, first_a_neg_date_pet), "years"))) %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(diff_time == min(diff_time)) %>%
  dplyr::ungroup()
mri_plot_data_temp <- mri_plot_data_temp %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup()

# creating list of variables that need to be ICV-adjusted
volumes_mri_plot_data <- mri_plot_data %>%
  dplyr::select(contains("CV") | contains("SV") | contains("TA"), -STATUS)

variables_mri_plot_data <- c()

for (roi_num in 1:length(names(volumes_mri_plot_data))){
  current_roi <- names(volumes_mri_plot_data)[roi_num]
  variables_mri_plot_data <- append(variables_mri_plot_data, current_roi)
}

# regression for age-adjustment
for (i in variables_mri_plot_data){
  fit <- lm(mri_plot_data_temp[[i]] ~ age + poly(age, 2, raw = TRUE)[,"2"], data = mri_plot_data_temp, na.action=na.exclude)
  mri_plot_data[i] <- mri_plot_data[paste0(i)] - (mri_plot_data$age - mean(mri_plot_data_temp$age, na.rm = T)) * fit$coefficients[2]
}

#############################
# Filtering to relevant RID's
#############################

#getting lm_data into the dataset
mri_plot_data <- mri_plot_data %>%
  dplyr::filter(RID %in% unique(centiloid_plot_data$RID)) %>%
  dplyr::left_join(lm_data)

#############################
# Calculating Meta-ROI and Hippocampal Volume
#############################
mri_plot_data <- mri_plot_data %>%
  dplyr::mutate(meta_ROI_left = (ST24TA + ST32TA + ST40TA + ST26TA)/4,  # meta_ROI is mean cortical thickness in the following individual ROIs: entorhinal, inferior temporal, middle temporal, and fusiform.
                meta_ROI_right = (ST83TA + ST91TA + ST99TA + ST85TA)/4,
                right_hippocampal_volume = ST88SV_adjusted,
                left_hippocampal_volume = ST29SV_adjusted)
mri_plot_data <- mri_plot_data %>%
  dplyr::mutate(adjusted_new_time = lubridate::time_length(difftime(EXAMDATE, predicted_date), "years"),
                hippocampal_volume = (right_hippocampal_volume + left_hippocampal_volume)/2,
                meta_ROI = (meta_ROI_right + meta_ROI_left)/2)

#putting together Meta-ROI data
meta_roi_plot_data <- mri_plot_data %>%
  dplyr::filter(!is.na(meta_ROI))
meta_roi_plot_data <- meta_roi_plot_data %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup()
meta_roi_plot_data <- meta_roi_plot_data %>%
  dplyr::select(RID, diags, meta_ROI, PTGENDER, PTEDUCAT, apoe, age, adjusted_new_time)

write.csv(meta_roi_plot_data, "~/Projects/GAM_models_Amyloid_converted/gam_modelling_data/meta_roi_gam_data_no_adni1.csv")

#putting together Hippocampal Volume data
hippocampal_volume_plot_data <- mri_plot_data %>%
  dplyr::filter(!is.na(hippocampal_volume))
hippocampal_volume_plot_data <- hippocampal_volume_plot_data %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup()
# hippocampal_volume_plot_data <- hippocampal_volume_plot_data %>%
#   dplyr::mutate(hippocampal_volume = scale(hippocampal_volume))
hippocampal_volume_plot_data <- hippocampal_volume_plot_data %>%
  dplyr::select(RID, diags, hippocampal_volume, PTGENDER, PTEDUCAT, apoe, age, adjusted_new_time)

write.csv(hippocampal_volume_plot_data, "~/Projects/GAM_models_Amyloid_converted/gam_modelling_data/hippocampal_volume_gam_data_no_adni1.csv")

#####################################################################################
# ADAS13
#####################################################################################

#Using Processing Pipeline from CIND
adas_1_score <- readr::read_delim("~/Data/ADASSCORES_23Jul2024.csv") %>% # edit to file name and directory for ADAS Scores, ADNI 1
  dplyr::rename(VISCODE2 = VISCODE)
adas_2_3_go <- readr::read_delim("~/Data/ADAS_ADNIGO23_23Jul2024.csv") %>% # edit to file name and directory for ADAS Scores, ADNI GO/2/3
  dplyr::mutate(ADAS13=TOTAL13,EXAMDATE=VISDATE)

## recodes missing data for total ADAS scores in ADNI1
adas_1_score$ADAS13 <- ifelse(adas_1_score$TOTALMOD==-4,NA,adas_1_score$TOTALMOD)

## joins two ADAS tables
adas_scores <- dplyr::bind_rows(adas_2_3_go,adas_1_score) %>%
  dplyr::select(RID,VISCODE2,EXAMDATE,ADAS13) %>%
  dplyr::group_by(RID,VISCODE2) %>%
  dplyr::mutate(across(.cols=dplyr::where(is.numeric),.fns=~mean(.,na.rm=TRUE))) %>%
  dplyr::distinct_at(vars(RID),.keep_all=TRUE)

adas13_plot_data <- adas_scores %>%
  dplyr::mutate(EXAMDATE = as.Date(EXAMDATE),
                RID = as.character(RID))

adas13_plot_data <- adas13_plot_data %>%
  dplyr::filter(RID %in% unique(centiloid_plot_data$RID)) %>%
  dplyr::left_join(lm_data) %>%
  dplyr::mutate(adjusted_new_time = lubridate::time_length(difftime(EXAMDATE, predicted_date), "years"))

#getting demographics
adas13_plot_data <- as.data.frame(adas13_plot_data) %>%
  dplyr::left_join(dem, by = "RID") %>%
  dplyr::mutate(age = lubridate::time_length(difftime(EXAMDATE, birthdate), "years"))
adas13_plot_data <- adas13_plot_data %>%
  dplyr::mutate(age = round(age, digits = 1))
adas13_plot_data$RID <- as.numeric(adas13_plot_data$RID)
adas13_plot_data <- merge(diagnoses, adas13_plot_data, all = TRUE) %>%
  dplyr::mutate(time_diff = abs(difftime(EXAMDATE, DX.DATE))) %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup()
adas13_plot_data <- merge(cdglobal_data, adas13_plot_data, all = TRUE) %>%
  dplyr::mutate(time_diff = abs(difftime(EXAMDATE, as.Date(CDR.DATE)))) %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()

adas13_plot_data <- adas13_plot_data %>%
  dplyr::filter(!is.na(ADAS13))
adas13_plot_data <- adas13_plot_data %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup()

adas13_plot_data <- adas13_plot_data %>%
  dplyr::select(RID, CDGLOBAL, diags, ADAS13, PTGENDER, PTEDUCAT, apoe, age, adjusted_new_time)

write.csv(adas13_plot_data, "~/Projects/GAM_models_Amyloid_converted/gam_modelling_data/adas13_gam_data.csv")

########################################################################################################
#MMSE
########################################################################################################
mmse <- readr::read_delim("~/Data/MMSE_18Jun2024.csv") %>% # edit to your source file name + location 
  dplyr::mutate(MMSE=MMSCORE,
                EXAMDATE = as.Date(VISDATE),
                RID = as.character(RID)) %>%
  dplyr::mutate(VISCODE2=case_when(
    VISCODE2 == "sc" ~ "bl",
    TRUE ~ VISCODE2
  )) %>%
  dplyr::group_by(RID,VISCODE2) %>%
  dplyr::mutate(across(.cols=dplyr::where(is.numeric),.fns=~mean(.,na.rm=TRUE))) %>%
  dplyr::distinct_at(vars(RID),.keep_all=TRUE)

mmse_plot_data <- mmse %>%
  dplyr::filter(RID %in% unique(centiloid_plot_data$RID)) %>%
  dplyr::left_join(lm_data) %>%
  dplyr::mutate(adjusted_new_time = lubridate::time_length(difftime(EXAMDATE, predicted_date), "years"))

mmse_plot_data <- mmse_plot_data %>%
  dplyr::filter(!is.na(MMSE))
mmse_plot_data <- mmse_plot_data %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup()

#getting demographics
mmse_plot_data <- as.data.frame(mmse_plot_data) %>%
  dplyr::left_join(dem, by = "RID") %>%
  dplyr::mutate(age = lubridate::time_length(difftime(EXAMDATE, birthdate), "years"))
mmse_plot_data <- mmse_plot_data %>%
  dplyr::mutate(age = round(age, digits = 1))
mmse_plot_data <- merge(diagnoses, mmse_plot_data, all = TRUE) %>%
  dplyr::mutate(time_diff = abs(difftime(EXAMDATE, DX.DATE))) %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup()
mmse_plot_data <- merge(cdglobal_data, mmse_plot_data, all = TRUE) %>%
  dplyr::mutate(time_diff = abs(difftime(EXAMDATE, as.Date(CDR.DATE)))) %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()

mmse_plot_data <- mmse_plot_data %>%
  dplyr::select(RID, CDGLOBAL, diags, MMSE, PTGENDER, PTEDUCAT, apoe, age, adjusted_new_time)

write.csv(mmse_plot_data, "~/Projects/GAM_models_Amyloid_converted/gam_modelling_data/mmse_gam_data.csv")

#########################################################################################################
# CDRSB
#########################################################################################################
cdrsb_plot_data <- read.csv("~/Data/CDR_24Apr2024.csv")
cdrsb_plot_data$CDR.DATE <- cdrsb_plot_data$USERDATE
cdrsb_plot_data <- cdrsb_plot_data[cdrsb_plot_data$CDGLOBAL>=0,]
cdrsb_plot_data$CDRSB <- rowSums(cdrsb_plot_data[,c('CDMEMORY','CDORIENT','CDJUDGE','CDCOMMUN','CDHOME','CDCARE')])
cdrsb_plot_data$RID <- as.character(cdrsb_plot_data$RID)

#getting lm_data in
cdrsb_plot_data <- cdrsb_plot_data %>%
  dplyr::filter(RID %in% unique(centiloid_plot_data$RID)) %>%
  dplyr::left_join(lm_data) %>%
  dplyr::mutate(adjusted_new_time = lubridate::time_length(difftime(CDR.DATE, predicted_date), "years"))

#getting demographics
cdrsb_plot_data <- as.data.frame(cdrsb_plot_data) %>%
  dplyr::left_join(dem, by = "RID") %>%
  dplyr::mutate(age = lubridate::time_length(difftime(CDR.DATE, birthdate), "years"))
cdrsb_plot_data <- cdrsb_plot_data %>%
  dplyr::mutate(age = round(age, digits = 1))
cdrsb_plot_data$RID <- as.numeric(cdrsb_plot_data$RID)

cdrsb_plot_data <- merge(diagnoses, cdrsb_plot_data, all = TRUE) %>%
  dplyr::mutate(time_diff = abs(difftime(CDR.DATE, DX.DATE))) %>%
  dplyr::group_by(RID, CDR.DATE) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup()

cdrsb_plot_data <- cdrsb_plot_data %>%
  dplyr::filter(!is.na(CDRSB))

cdrsb_plot_data <- cdrsb_plot_data %>%
  dplyr::select(RID, CDGLOBAL, diags, CDRSB, PTGENDER, PTEDUCAT, apoe, age, adjusted_new_time)

write.csv(cdrsb_plot_data, "~/Projects/GAM_models_Amyloid_converted/gam_modelling_data/cdrsb_gam_data.csv")

##########################################################################################################
# mPACCtrailsB
##########################################################################################################
# Using Processing Pipeline from CIND
neurobat_source <- readr::read_delim("~/Data/NEUROBAT_23Jul2024.csv") # change to name and location of Neuropsychological Battery

debugged_pacc_fn<-function (dd, keepComponents = FALSE)
{
  require(dplyr)
  require(tidyr)
  holdNames <- colnames(dd)
  dd <- mutate(dd, log.TRABSCOR = log(TRABSCOR + 1))
  dd$log.TRABSCOR <- ifelse(is.finite(dd$log.TRABSCOR),dd$log.TRABSCOR,NA)
  bl.summary <- dplyr::filter(dd, VISCODE == "bl") %>% dplyr::select(DX.bl,
                                                                     ADASQ4, LDELTOTAL, DIGITSCOR, log.TRABSCOR, MMSE) %>%
    tidyr::gather(VARIABLE, SCORE, -DX.bl) %>% dplyr::filter(!is.na(SCORE)) %>%
    dplyr::group_by(DX.bl, VARIABLE) %>% dplyr::summarize(N = n(),
                                                          mean = mean(SCORE), sd = sd(SCORE))
  zscore <- function(x, var) {
    (x - dplyr::filter(bl.summary, DX.bl == "CN" & VARIABLE == var)$mean)/dplyr::filter(bl.summary,
                                                                                        DX.bl == "CN" & VARIABLE == var)$sd
  }
  dd <- dplyr::mutate(dd, ADASQ4.z = -zscore(ADASQ4, "ADASQ4"), LDELTOTAL.z = zscore(LDELTOTAL,
                                                                                     "LDELTOTAL"), DIGITSCOR.z = zscore(DIGITSCOR, "DIGITSCOR"),
                      log.TRABSCOR.z = -zscore(log.TRABSCOR, "log.TRABSCOR"),
                      MMSE.z = zscore(MMSE, "MMSE"))
  corTest <- cor(dplyr::select(dd, ADASQ4.z, LDELTOTAL.z, DIGITSCOR.z,
                               log.TRABSCOR.z, MMSE.z), use = "pairwise.complete.obs")
  if (any(corTest < 0))
    stop("Some PACC z scores are negatively correlated!")
  compscore <- function(x, n.components = 4, n.missing = 2) {
    ifelse(sum(is.na(x)) > n.missing, NA, mean(x, na.rm = TRUE)) *
      n.components
  }
  dd$mPACCdigit <- apply(dd[, c("ADASQ4.z", "LDELTOTAL.z",
                                "DIGITSCOR.z", "MMSE.z")], 1, compscore)
  dd$mPACCtrailsB <- apply(dd[, c("ADASQ4.z", "LDELTOTAL.z",
                                  "log.TRABSCOR.z", "MMSE.z")], 1, compscore)
  if (!keepComponents)
    dd <- dd[, c(holdNames, "mPACCdigit", "mPACCtrailsB")]
  dd <- as.data.frame(dd) %>%
    dplyr::distinct_at(vars(RID,VISCODE,mPACCdigit,mPACCtrailsB))
}

## LDEL and DIGITSCOR/TRABSCOR are recorded at separate baseline/screening visits - below code aligns to one visit
neurobat_non_bl <- neurobat_source %>%
  dplyr::filter(!(VISCODE2 %in% c("bl","sc","f"))) %>%
  dplyr::select(RID,VISCODE2,LDELTOTAL,DIGITSCOR,TRABSCOR)

neurobat_ldel_bl <- neurobat_source %>%
  dplyr::filter(VISCODE2 %in% c("bl","sc","f")) %>%
  dplyr::select(RID,VISCODE2,LDELTOTAL) %>%
  dplyr::mutate(VISCODE2="bl") %>%
  dplyr::arrange(RID)

neurobat_scores_bl <- neurobat_source %>%
  dplyr::filter(VISCODE2 %in% c("bl","sc","f")) %>%
  dplyr::select(RID,VISCODE2,DIGITSCOR,TRABSCOR) %>%
  dplyr::mutate(VISCODE2="bl") %>%
  dplyr::arrange(RID)

## selects correct merged rows
neurobat_bl <- dplyr::full_join(neurobat_ldel_bl,neurobat_scores_bl,by=c("RID","VISCODE2")) %>%
  dplyr::mutate(total_na = is.na(DIGITSCOR)+is.na(TRABSCOR)+is.na(LDELTOTAL)) %>%
  dplyr::arrange(total_na) %>%
  dplyr::distinct_at(vars(RID,VISCODE2),.keep_all=TRUE)

## averages records for subjects with multiple batteries evaluated
neurobat <- dplyr::bind_rows(neurobat_bl,neurobat_non_bl) %>%
  dplyr::group_by(RID,VISCODE2) %>%
  dplyr::mutate(across(.cols=dplyr::where(is.numeric),.fns=~mean(.,na.rm=TRUE))) %>%
  dplyr::distinct_at(vars(RID),.keep_all=TRUE)

adas_q4_for_pacc_1 <- adas_1_score %>%
  dplyr::select(RID,VISCODE2,Q4) %>%
  dplyr::rename(ADASQ4=Q4)

adas_q4_for_pacc_2_3_go <- adas_2_3_go %>%
  dplyr::select(RID,VISCODE2,Q4SCORE) %>%
  dplyr::rename(ADASQ4=Q4SCORE)

adas_q4_for_pacc <- dplyr::bind_rows(adas_q4_for_pacc_1,adas_q4_for_pacc_2_3_go) %>%
  dplyr::group_by(RID,VISCODE2) %>%
  dplyr::mutate(across(.cols=dplyr::where(is.numeric),.fns=~mean(.,na.rm=TRUE))) %>%
  dplyr::distinct_at(vars(RID),.keep_all=TRUE)

pacc_df <- dplyr::left_join(neurobat,
                            adas_q4_for_pacc,
                            by=c("RID","VISCODE2")) %>%
  dplyr::left_join(.,
                   mmse %>%
                     dplyr::mutate(RID = as.numeric(RID)) %>%
                     dplyr::select(RID,VISCODE2,MMSE),
                   by=c("RID","VISCODE2")) %>%
  dplyr::left_join(.,
                   adni_diagnoses %>%
                     dplyr::rename(DX.bl=DX) %>%
                     dplyr::mutate(RID = as.numeric(RID)) %>%
                     dplyr::filter(VISCODE2=="bl") %>%
                     dplyr::select(RID,DX.bl),
                   by=c("RID")) %>%
  dplyr::mutate(log.TRABSCOR=log(TRABSCOR+1)) %>%
  dplyr::select(RID,VISCODE2,DX.bl,ADASQ4, LDELTOTAL, DIGITSCOR,TRABSCOR,MMSE) %>%
  tidyr::drop_na(DX.bl) %>%
  dplyr::mutate(DX.bl =
                  case_when(DX.bl == "CU" ~ "CN",
                            TRUE ~ DX.bl)) %>%
  dplyr::rename(VISCODE=VISCODE2)

pacc <- debugged_pacc_fn(dd=data.frame(pacc_df)) %>%
  dplyr::group_by(RID,VISCODE) %>%
  dplyr::mutate(across(.cols=dplyr::where(is.numeric),.fns=~mean(.,na.rm=TRUE))) %>%
  dplyr::distinct_at(vars(RID),.keep_all=TRUE)

#getting dates
registry <- read.csv("~/Data/REGISTRY_23Jul2024.csv") %>%
  dplyr::select(RID, VISCODE2, EXAMDATE)

mpacctrailsb_plot_data <- pacc %>%
  dplyr::filter(RID %in% unique(centiloid_plot_data$RID))
mpacctrailsb_plot_data <- merge(mpacctrailsb_plot_data, registry, by.x = c("RID", "VISCODE"), by.y = c("RID", "VISCODE2"), all.x = TRUE)
mpacctrailsb_plot_data <- mpacctrailsb_plot_data %>%
  dplyr::filter(!(EXAMDATE == "" | is.na(EXAMDATE)),
                !(is.na(mPACCtrailsB))) %>%
  dplyr::mutate(mPACCtrailsB = as.numeric(mPACCtrailsB),
                RID = as.character(RID))
#getting lm_data in
mpacctrailsb_plot_data <- mpacctrailsb_plot_data %>%
  dplyr::left_join(lm_data) %>%
  dplyr::mutate(adjusted_new_time = lubridate::time_length(difftime(EXAMDATE, predicted_date), "years"))

#getting demographics
mpacctrailsb_plot_data <- as.data.frame(mpacctrailsb_plot_data) %>%
  dplyr::left_join(dem, by = "RID") %>%
  dplyr::mutate(age = lubridate::time_length(difftime(EXAMDATE, birthdate), "years"))
mpacctrailsb_plot_data <- mpacctrailsb_plot_data %>%
  dplyr::mutate(age = round(age, digits = 1))
mpacctrailsb_plot_data$RID <- as.numeric(mpacctrailsb_plot_data$RID)

mpacctrailsb_plot_data <- merge(diagnoses, mpacctrailsb_plot_data, all = TRUE) %>%
  dplyr::mutate(time_diff = abs(difftime(EXAMDATE, DX.DATE))) %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup()
mpacctrailsb_plot_data <- merge(cdglobal_data, mpacctrailsb_plot_data, all = TRUE) %>%
  dplyr::mutate(time_diff = abs(difftime(EXAMDATE, as.Date(CDR.DATE)))) %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()

mpacctrailsb_plot_data <- mpacctrailsb_plot_data %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup()
# mpacctrailsb_plot_data <- mpacctrailsb_plot_data %>%
#   dplyr::mutate(mPACCtrailsB = scale(mPACCtrailsB))
mpacctrailsb_plot_data <- mpacctrailsb_plot_data %>%
  dplyr::select(RID, CDGLOBAL, diags, mPACCtrailsB, PTGENDER, PTEDUCAT, apoe, age, adjusted_new_time)

write.csv(mpacctrailsb_plot_data, "~/Projects/GAM_models_Amyloid_converted/gam_modelling_data/mpacctrailsb_gam_data.csv")

###########################################################################################################
# TAU
###########################################################################################################
## CSF data load-in
upenn_csf_master <- readr::read_delim("~/Data/UPENNBIOMK_MASTER_FINAL.csv")
upenn_csf_master <- upenn_csf_master %>% 
  dplyr::filter(BATCH %in% c("UPENNBIOMK9","UPENNBIOMK10",
                             "UPENNBIOMK12","UPENNBIOMK13"))
upenn_csf_master <- upenn_csf_master %>% 
  dplyr::arrange(RUNDATE) %>% 
  dplyr::distinct_at(vars(RID,VISCODE2),.keep_all=TRUE)

## upenn_merged_csf_biomarkers <- dplyr::bind_rows(upenn_mk_12,upenn_mk_10,upenn_mk_9)
upenn_merged_csf_biomarkers <- upenn_csf_master %>% 
  dplyr::mutate(ABETA=ABETA42)
upenn_merged_csf_biomarkers <- upenn_merged_csf_biomarkers %>% 
  dplyr::select(-VISCODE) %>% dplyr::rename(VISCODE=VISCODE2)
upenn_merged_csf_biomarkers <- upenn_merged_csf_biomarkers %>% 
  dplyr::mutate(ABETA=round(ABETA),
                ptau_pos=PTAU>21.8,
                AmyloidPosCSF=ABETA<980,
                ptau_ab_ratio=(as.numeric(PTAU)/as.numeric(ABETA)),ad_path_pos=ptau_ab_ratio>0.025)
upenn_merged_csf_biomarkers <- upenn_merged_csf_biomarkers %>% 
  dplyr::arrange(RUNDATE) %>% 
  dplyr::distinct_at(vars(RID,VISCODE),.keep_all=TRUE)

## split ABETA values here; values <200 or >1700 are usable for status assignment but not for numerical analysis
## don't need to split CSF PTAU - all values are within technical limits
upenn_merged_csf_biomarkers <- upenn_merged_csf_biomarkers %>% 
  dplyr::mutate(ABETA_for_status=ABETA,
                ABETA_for_biomarkers=case_when(ABETA>=200 & ABETA<=1700 ~ ABETA),
                RID = as.character(RID),
                VISCODE_csf = VISCODE)
#getting lm_data in
csf_plot_data <- upenn_merged_csf_biomarkers %>%
  dplyr::filter(RID %in% unique(centiloid_plot_data$RID)) %>%
  dplyr::left_join(lm_data) %>%
  dplyr::mutate(adjusted_new_time = lubridate::time_length(difftime(EXAMDATE, predicted_date), "years"),
                TAU = as.numeric(TAU),
                PTAU = as.numeric(PTAU))

#getting demographics
csf_plot_data <- as.data.frame(csf_plot_data) %>%
  dplyr::left_join(dem, by = "RID") %>%
  dplyr::mutate(age = lubridate::time_length(difftime(EXAMDATE, birthdate), "years"))
csf_plot_data <- csf_plot_data %>%
  dplyr::mutate(age = round(age, digits = 1))
csf_plot_data$RID <- as.numeric(csf_plot_data$RID)

csf_plot_data <- merge(diagnoses, csf_plot_data, all = TRUE) %>%
  dplyr::mutate(time_diff = abs(difftime(EXAMDATE, DX.DATE))) %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup()
csf_plot_data <- merge(cdglobal_data, csf_plot_data, all = TRUE) %>%
  dplyr::mutate(time_diff = abs(difftime(EXAMDATE, as.Date(CDR.DATE)))) %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()

tau_plot_data <- csf_plot_data %>%
  dplyr::filter(!is.na(TAU)) #%>%
  # dplyr::mutate(TAU = scale(TAU))
tau_plot_data <- tau_plot_data %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup()

tau_plot_data <- tau_plot_data %>%
  dplyr::select(RID, CDGLOBAL, diags, TAU, PTGENDER, PTEDUCAT, apoe, age, adjusted_new_time)

write.csv(tau_plot_data, "~/Projects/GAM_models_Amyloid_converted/gam_modelling_data/tau_gam_data.csv")

###########################################################################################################
# PTAU
###########################################################################################################
ptau_plot_data <- csf_plot_data %>%
  dplyr::filter(!is.na(PTAU))
ptau_plot_data <- ptau_plot_data %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup()

ptau_plot_data <- ptau_plot_data %>%
  dplyr::select(RID, CDGLOBAL, diags, PTAU, PTGENDER, PTEDUCAT, apoe, age, adjusted_new_time)

write.csv(ptau_plot_data, "~/Projects/GAM_models_Amyloid_converted/gam_modelling_data/ptau_gam_data.csv")

###########################################################################################################
# ABETA
###########################################################################################################
abeta_plot_data <- csf_plot_data %>%
  dplyr::filter(!is.na(ABETA))
abeta_plot_data <- abeta_plot_data %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup()

abeta_plot_data <- abeta_plot_data %>%
  dplyr::select(RID, CDGLOBAL, diags, ABETA, PTGENDER, PTEDUCAT, apoe, age, adjusted_new_time)

write.csv(abeta_plot_data, "~/Projects/GAM_models_Amyloid_converted/gam_modelling_data/abeta_gam_data.csv")

############################################################################################################
#E-COG: Subject
############################################################################################################
d_s <- read.csv("~/Data/ECOGPT_08Apr2024.csv")
d_s[d_s == "9"] <- NA
# d_s <- d_s[, colMeans(is.na(d_s)) <= .3]
ds_test <- d_s %>%
  dplyr::filter(!is.na(EcogPtTotal))
ds_test <- ds_test %>%
  dplyr::filter(RID %in% unique(centiloid_plot_data$RID))

total_s <- c('MEMORY1', 'MEMORY2', 'MEMORY3', 'MEMORY4', 'MEMORY5', 'MEMORY6', 'MEMORY7', 'MEMORY8',
           'LANG1', 'LANG2', 'LANG3', 'LANG4', 'LANG5', 'LANG6', 'LANG7', 'LANG8', 'LANG9',
           'VISSPAT1', 'VISSPAT2', 'VISSPAT3', 'VISSPAT4',
           'VISSPAT6', 'VISSPAT7', 'VISSPAT8',
           'PLAN1', 'PLAN2', 'PLAN3', 'PLAN4', 'PLAN5',
           'ORGAN1', 'ORGAN2', 'ORGAN3', 'ORGAN4', 'ORGAN5', 'ORGAN6',
           'DIVATT1', 'DIVATT2', 'DIVATT3', 'DIVATT4')
total.coef_s <- c(0.89,0.83,0.83,0.87,0.82,0.90,0.83,0.90,0.73,0.78,0.72,0.81,0.75,
                0.63,0.55,0.68,0.65,0.82,0.84,0.91,0.94,0.89,0.90,0.89,0.95,0.86,
                0.91,0.80,0.87,0.74,0.91,0.91,0.90,0.91,0.85,0.85,0.88,0.84,0.85)
mem_s <- c('MEMORY1', 'MEMORY2', 'MEMORY3', 'MEMORY4', 'MEMORY5', 'MEMORY6', 'MEMORY7', 'MEMORY8')
mem.coef_s <- c(0.26,0.45,0.47,0.21,0.47,0.23,0.41,0.17)
lang_s <-  c('LANG1', 'LANG2', 'LANG3', 'LANG4', 'LANG5', 'LANG6', 'LANG7', 'LANG8', 'LANG9')
lang.coef_s <- c(0.55,0.45,0.55,0.45,0.44,0.63,0.71,0.62,0.67)
visspat_s <- c('VISSPAT1', 'VISSPAT2', 'VISSPAT3', 'VISSPAT4',
             'VISSPAT6', 'VISSPAT7', 'VISSPAT8')
visspat.coef_s <- c(0.56,0.51,0.25,0.23,0.35,0.37,0.33)
plan_s <- c('PLAN1', 'PLAN2', 'PLAN3', 'PLAN4', 'PLAN5')
plan.coef_s <- c(0.20,0.27,0.29,0.58,0.42)
organ_s <- c('ORGAN1', 'ORGAN2', 'ORGAN3', 'ORGAN4', 'ORGAN5', 'ORGAN6')
organ.coef_s <- c(0.46,0.35,0.38,0.36,0.23,0.38)
divatt_s <- c('DIVATT1', 'DIVATT2', 'DIVATT3', 'DIVATT4')
divatt.coef_s <- c(0.42,0.37,0.42,0.44)

d_s$EcogGlobal <- rowMeans(d_s[,total_s]*total.coef_s,na.rm=T)

ecog_s_plot_data <- d_s %>%
  dplyr::select(RID, VISDATE, EcogGlobal) %>%
  dplyr::filter(!is.na(EcogGlobal)) %>%
  dplyr::mutate(RID = as.character(RID),
                EXAMDATE = as.Date(VISDATE))

#getting lm_data in
ecog_s_plot_data <- ecog_s_plot_data %>%
  dplyr::filter(RID %in% unique(centiloid_plot_data$RID)) %>%
  dplyr::left_join(lm_data) %>%
  dplyr::mutate(adjusted_new_time = lubridate::time_length(difftime(EXAMDATE, predicted_date), "years"),
                EcogGlobal = as.numeric(EcogGlobal))

#getting demographics
ecog_s_plot_data <- as.data.frame(ecog_s_plot_data) %>%
  dplyr::left_join(dem, by = "RID") %>%
  dplyr::mutate(age = lubridate::time_length(difftime(EXAMDATE, birthdate), "years"))
ecog_s_plot_data <- ecog_s_plot_data %>%
  dplyr::mutate(age = round(age, digits = 1))
ecog_s_plot_data$RID <- as.numeric(ecog_s_plot_data$RID)

ecog_s_plot_data <- merge(diagnoses, ecog_s_plot_data, all = TRUE) %>%
  dplyr::mutate(time_diff = abs(difftime(VISDATE, DX.DATE))) %>%
  dplyr::group_by(RID, VISDATE) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup() #621 cases if global measures, 246 if EcogTotal measures
ecog_s_plot_data <- merge(cdglobal_data, ecog_s_plot_data, all = TRUE) %>%
  dplyr::mutate(time_diff = abs(difftime(VISDATE, as.Date(CDR.DATE)))) %>%
  dplyr::group_by(RID, VISDATE) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()

ecog_s_plot_data <- ecog_s_plot_data %>%
  dplyr::select(RID, CDGLOBAL, diags, EcogGlobal, PTGENDER, PTEDUCAT, apoe, age, adjusted_new_time)

write.csv(ecog_s_plot_data, "~/Projects/GAM_models_Amyloid_converted/gam_modelling_data/ecog_s_gam_data.csv")

############################################################################################################
#E-COG: Partner
############################################################################################################
d_p <- read.csv("~/Data/ECOGSP_08Apr2024.csv")
d_p[d_p == "9"] <- NA
total_p <- c('MEMORY1', 'MEMORY2', 'MEMORY3', 'MEMORY4', 'MEMORY5', 'MEMORY6', 'MEMORY7', 'MEMORY8',
             'LANG1', 'LANG2', 'LANG3', 'LANG4', 'LANG5', 'LANG6', 'LANG7', 'LANG8', 'LANG9',
             'VISSPAT1', 'VISSPAT2', 'VISSPAT3', 'VISSPAT4',
             'VISSPAT6', 'VISSPAT7', 'VISSPAT8',
             'PLAN1', 'PLAN2', 'PLAN3', 'PLAN4', 'PLAN5',
             'ORGAN1', 'ORGAN2', 'ORGAN3', 'ORGAN4', 'ORGAN5', 'ORGAN6',
             'DIVATT1', 'DIVATT2', 'DIVATT3', 'DIVATT4')
total.coef_p <- c(0.89,0.83,0.83,0.87,0.82,0.90,0.83,0.90,0.73,0.78,0.72,0.81,0.75,
                  0.63,0.55,0.68,0.65,0.82,0.84,0.91,0.94,0.89,0.90,0.89,0.95,0.86,
                  0.91,0.80,0.87,0.74,0.91,0.91,0.90,0.91,0.85,0.85,0.88,0.84,0.85)
mem_p <- c('MEMORY1', 'MEMORY2', 'MEMORY3', 'MEMORY4', 'MEMORY5', 'MEMORY6', 'MEMORY7', 'MEMORY8')
mem.coef_p <- c(0.26,0.45,0.47,0.21,0.47,0.23,0.41,0.17)
lang_p <-  c('LANG1', 'LANG2', 'LANG3', 'LANG4', 'LANG5', 'LANG6', 'LANG7', 'LANG8', 'LANG9')
lang.coef_p <- c(0.55,0.45,0.55,0.45,0.44,0.63,0.71,0.62,0.67)
visspat_p <- c('VISSPAT1', 'VISSPAT2', 'VISSPAT3', 'VISSPAT4',
               'VISSPAT6', 'VISSPAT7', 'VISSPAT8')
visspat.coef_p <- c(0.56,0.51,0.25,0.23,0.35,0.37,0.33)
plan_p <- c('PLAN1', 'PLAN2', 'PLAN3', 'PLAN4', 'PLAN5')
plan.coef_p <- c(0.20,0.27,0.29,0.58,0.42)
organ_p <- c('ORGAN1', 'ORGAN2', 'ORGAN3', 'ORGAN4', 'ORGAN5', 'ORGAN6')
organ.coef_p <- c(0.46,0.35,0.38,0.36,0.23,0.38)
divatt_p <- c('DIVATT1', 'DIVATT2', 'DIVATT3', 'DIVATT4')
divatt.coef_p <- c(0.42,0.37,0.42,0.44)

d_p$EcogGlobal <- rowMeans(d_p[,total_p]*total.coef_p,na.rm=T)

ecog_p_plot_data <- d_p %>%
  dplyr::select(RID, VISDATE, EcogGlobal) %>%
  dplyr::filter(!is.na(EcogGlobal)) %>%
  dplyr::mutate(RID = as.character(RID),
                EXAMDATE = as.Date(VISDATE))

#getting lm_data in
ecog_p_plot_data <- ecog_p_plot_data %>%
  dplyr::filter(RID %in% unique(centiloid_plot_data$RID)) %>%
  dplyr::left_join(lm_data) %>%
  dplyr::mutate(adjusted_new_time = lubridate::time_length(difftime(EXAMDATE, predicted_date), "years"),
                EcogGlobal = as.numeric(EcogGlobal))

#getting demographics
ecog_p_plot_data <- as.data.frame(ecog_p_plot_data) %>%
  dplyr::left_join(dem, by = "RID") %>%
  dplyr::mutate(age = lubridate::time_length(difftime(EXAMDATE, birthdate), "years"))
ecog_p_plot_data <- ecog_p_plot_data %>%
  dplyr::mutate(age = round(age, digits = 1))
ecog_p_plot_data$RID <- as.numeric(ecog_p_plot_data$RID)

ecog_p_plot_data <- merge(diagnoses, ecog_p_plot_data, all = TRUE) %>%
  dplyr::mutate(time_diff = abs(difftime(VISDATE, DX.DATE))) %>%
  dplyr::group_by(RID, VISDATE) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup()
ecog_p_plot_data <- merge(cdglobal_data, ecog_p_plot_data, all = TRUE) %>%
  dplyr::mutate(time_diff = abs(difftime(VISDATE, as.Date(CDR.DATE)))) %>%
  dplyr::group_by(RID, VISDATE) %>%
  dplyr::filter(time_diff == min(time_diff)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()

ecog_p_plot_data <- ecog_p_plot_data %>%
  dplyr::select(RID, CDGLOBAL, diags, EcogGlobal, PTGENDER, PTEDUCAT, apoe, age, adjusted_new_time)

write.csv(ecog_p_plot_data, "~/Projects/GAM_models_Amyloid_converted/gam_modelling_data/ecog_p_gam_data.csv")
