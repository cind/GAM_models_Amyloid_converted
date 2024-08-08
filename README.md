# GAM_models_Amyloid_converted_SAA
############################################################################################
# Repo used for Paper: "Timing of Alzheimer’s Disease Biomarker Progressions: A Two-Decade 
# Observational Study from the Alzheimer’s Disease Neuroimaging Initiative (ADNI)
# By: Schaap, Thropp, & Tosun
#
# Project using sub-cohort of ADNI participants who converted from AB PET- to AB PET+
# using Centiloid = 20 as a threshold for AB status. We used GAM models to estimate 
# trajectories in several different biological and cognitive outcome measures to detect
# the point in time that slopes of the trajectories begin to significantly differ from
# a slope of 0. Below is a short description of the different files/folders in the repo.
############################################################################################

gam_helpers: contains functions needed to run gam_models

folder "gam_modelling_data": contains processed data ready for modelling

prepping data for GAM models: contains processing steps for data that is in the folder "gam_modelling_data"

table_info: contains steps taken to get information in Table 1 as well as answering some of the reviewer's comments

plotting_raw_data: script to create the supplementary plots of raw data for subjects

gam_models: script to run GAM models as well as plot creation for plot 2/3 in revised manuscript
