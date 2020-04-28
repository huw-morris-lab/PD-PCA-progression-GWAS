#------------GWAS OF PCA PROGRESSION IN PROBAND AND OXFORD AND PPMI COMBINED------------#
#Using PPMI MDS-UPDRS3 and H&Y data assessed in the OFF medication state
#Using PPMI data from baseline, only using ANNUAL visits for all measures
#Progression scores created in combined dataset
#Non-PD cases have been removed
#Normalising percentage scores within each cohort before mixed effects model
#This is to account for differences in scales between cohorts e.g. semantic fluency, H&Y stage
#Also creating progression scores for each scale separately
#Using same method of extracting random slope from mixed effects model but without PCA
#Including an export for sensitivity analyses excluding potential non-PD cases based on extreme 5%
#It will create a separate export file with the extreme 5% cases removed

#---Load packages---####
library(ggplot2)
library(readstata13)
library(gridExtra)
library(gtable)
library(grid)
library(reshape2)
library(pROC)
library(ROCR)
library(readxl)
library(survival)
library(survminer)
library(lubridate)
library(lme4)
library(car)
library(MASS)
library(factoextra)
library(data.table)
library(tidyverse)
library(ggpubr)
library(Hmisc)
library(corrplot)
library(MuMIn)

#---Load functions---####

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

## Norms the data within specified groups in a data frame; it normalizes each
## subject (identified by idvar) so that they have the same mean, within each group
## specified by betweenvars.
##   data: a data frame.
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   na.rm: a boolean that indicates whether to ignore NA's
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
  library(plyr)
  
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
                         .fun = function(xx, col, na.rm) {
                           c(subjMean = mean(xx[,col], na.rm=na.rm))
                         },
                         measurevar,
                         na.rm
  )
  
  # Put the subject means with original data
  data <- merge(data, data.subjMean)
  
  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "_norm", sep="")
  data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
    mean(data[,measurevar], na.rm=na.rm)
  
  # Remove this subject mean column
  data$subjMean <- NULL
  
  return(data)
}


## Summarizes data, handling within-subjects variables by removing inter-subject variability.
## It will still work if there are no within-S variables.
## Gives count, un-normed mean, normed mean (with same between-group mean),
##   standard deviation, standard error of the mean, and confidence interval.
## If there are within-subject variables, calculate adjusted values using method from Morey (2008).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   withinvars: a vector containing names of columns that are within-subjects variables
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  
  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
                       FUN=is.factor, FUN.VALUE=logical(1))
  
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  
  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  
  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
  
  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")
  
  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                                  FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  
  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  
  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}

#Normalise function
normFunc <- function(x){
  (x-mean(x, na.rm = T))/sd(x, na.rm = T)
}

#---TABLE 1 all cohorts: Summary stats f and demographics at baseline---####

#Count of N in each cohort
PD_only %>% 
  group_by(cohort) %>% 
  summarise(count = n())

#Count of gender in each cohort
PD_only %>% 
  group_by(cohort, gender) %>% 
  summarise(count = n())


#Summary stats for age at onset and age at entry at baseline
PD_only %>% 
  group_by(cohort) %>% 
  dplyr::summarise(AAO_mean = mean(age_onset_imput, na.rm = TRUE),
                   AAO_sd = sd(age_onset_imput, na.rm = TRUE))

#Mean disease duration at onset
PD_only %>% 
  group_by(cohort) %>% 
  dplyr::summarise(disdur_onset_mean = mean(BL_disease_duration_onset, na.rm = TRUE),
                   disdur_onset_sd = sd(BL_disease_duration_onset, na.rm = TRUE))

#Mean and SD MDS-UPDRS Part III at each timepoint
PD_only %>% 
  group_by(cohort) %>% 
  dplyr::summarise(UPDRSIII_mean = mean(BL_UPDRS_III_total, na.rm = TRUE),
                   UPDRSIII_sd = sd(BL_UPDRS_III_total, na.rm = TRUE))

#Mean and SD MDS-UPDRS Part II at each timepoint
PD_only %>% 
  group_by(cohort) %>% 
  dplyr::summarise(UPDRSII_mean = mean(BL_UPDRS_II_total, na.rm = TRUE),
                   UPDRSII_sd = sd(BL_UPDRS_II_total, na.rm = TRUE))

#Mean H&Y stage at baseline
PD_only %>% 
  group_by(cohort) %>% 
  dplyr::summarise(HY_mean = mean(BL_hoehn_and_yahr_stage, na.rm = TRUE),
                   HY_sd =sd(BL_hoehn_and_yahr_stage, na.rm = TRUE))

#Create new variable for Hoehn and Yahr stage grouped - 0 to 1.5 vs 2 to 2.5 vs. 3+
PD_only <- PD_only %>% 
  mutate(BL_hoehn_and_yahr_grouped = ifelse(BL_hoehn_and_yahr_stage == 0 | BL_hoehn_and_yahr_stage == 1 | BL_hoehn_and_yahr_stage == 1.5, "0 to 1.5",
                                            ifelse(BL_hoehn_and_yahr_stage == 2 | BL_hoehn_and_yahr_stage == 2.5, "2 to 2.5",
                                                   ifelse(BL_hoehn_and_yahr_stage == 3 | BL_hoehn_and_yahr_stage == 4 | BL_hoehn_and_yahr_stage == 5, "3+", NA))))

#H&Y stage proportions at baseline
PD_only %>% 
  group_by(cohort, BL_hoehn_and_yahr_grouped) %>% 
  summarise(count = n())

#Mean MoCA at baseline
PD_only %>%
  group_by(cohort) %>% 
  dplyr::summarise(MOCA_mean = mean(BL_MOCA_total, na.rm = TRUE),
                   MOCA_sd = sd(BL_MOCA_total, na.rm = TRUE))

#Mean MoCA at SC for PPMI
PD_only %>% 
  group_by(cohort) %>% 
  summarise(SC_MOCA_mean = mean(SC_MOCA_total, na.rm = TRUE),
            SC_MOCA_sd = sd(SC_MOCA_total, na.rm = TRUE))

#Semantic fluency scores
PD_only %>% 
  group_by(cohort) %>% 
  dplyr::summarise(SF_mean = mean(BL_seman_flu_score, na.rm = TRUE),
                   SF_sd =sd(BL_seman_flu_score, na.rm = TRUE))

#MDS-UPDRS Part 1.1 mean
PD_only %>% 
  group_by(cohort) %>% 
  summarise(updrs1_1_mean = mean(BL_UPDRS_I_1recode, na.rm = TRUE),
            updrs1_1_sd = sd(BL_UPDRS_I_1recode, na.rm = TRUE))

#---Summary stats on follow-up---####

#Total number of visits completed
#Calculate the number of visits per patient
disease_durations <- PD_only %>% 
  select(BL_hoehn_and_yahr_stage,
         BL_disease_duration_onset,
         V04_disease_duration_onset,
         V05_disease_duration_onset,
         V06_disease_duration_onset,
         V08_disease_duration_onset,
         V10_disease_duration_onset,
         V11_disease_duration_onset,
         V12_disease_duration_onset,
         V13_disease_duration_onset) %>% 
  mutate(BL_visit = ifelse(is.na(BL_disease_duration_onset), NA,
                           ifelse(!is.na(BL_disease_duration_onset), 1, "check"))) %>% 
  select(contains("disease_duration_onset"), -BL_disease_duration_onset, BL_visit) %>% 
  mutate(number_visits = rowSums(!is.na(.)))
#Note that this uses disease duration onset variable as an indication of whether the visit was completed
#There are some individuals who do not have disease duration data but have other data
#However we do not analyse the other clinical data if the time (disease duration) data is not available

#Calculate total number of visits across all patients
disease_durations %>% 
  summarise(total_visits = sum(number_visits))

#Calculate total number of visits for each cohort
PD_only_visits %>% 
  group_by(PD_only$cohort) %>% 
  summarise(total_visits = sum(number_visits))

#---Calculate the mean and median follow-up time in years---####

#Select all disease duration variables except baseline
#This is so we can find the last nonmissing disease duration value (the latest/most recent follow-up visit)
followup <- PD_only %>% 
  select(ID,
         V04_disease_duration_onset,
         V05_disease_duration_onset,
         V06_disease_duration_onset,
         V08_disease_duration_onset,
         V10_disease_duration_onset,
         V11_disease_duration_onset,
         V12_disease_duration_onset,
         V13_disease_duration_onset) 

#Find the disease duration onset value for the most recent follow-up visit
lastdiseaseduration <- followup %>% 
  gather(key = Key, value = Value, -ID) %>% 
  group_by(ID) %>% 
  dplyr::filter(!is.na(Value)) %>%
  slice(n()) %>% 
  select(-Key)

#Merge with main dataset
#Calculate time from baseline to the most recent follow-up visit for each patient
followup_summary <- PD_only %>% 
  select(ID, cohort, BL_disease_duration_onset) %>% 
  left_join(lastdiseaseduration, by = "ID") %>% 
  mutate(followup_time = Value - BL_disease_duration_onset)

followup_summary %>%
  summarise(mean_followup = mean(followup_time, na.rm = TRUE),
            median_followup = median(followup_time, na.rm = TRUE),
            sd_followup = sd(followup_time, na.rm = TRUE))

#---Identify how many patients have completed the 72 month visit---####

followup %>%
  summarise(count = sum(is.na(followup$V13_disease_duration_onset)))

########## NORMALISE VARIABLES TO POPULATION BASELINE MEAN AND SD WITHIN COHORT ########## 
#---Transform to percentages and normalise MDS-UPDRS Part III scores---####

#First transform UPDRSIII scores to a percentage of the maximum score
UPDRS3 <- PD_only %>% 
  select(contains("UPDRS_III_total"))

#For loop - convert UPDRS3 score to percentage (divide by 132 which is the total score)
for (nm in names(UPDRS3)){
  PD_only <- PD_only %>% 
    mutate(perc = PD_only[[nm]]/132 * 100)
  visitid <- word(nm, start = 1, sep = "\\_")
  names(PD_only)[names(PD_only) == "perc"] <- paste(visitid, "UPDRS_III_total_perc", sep = "_")
} 

#---Normalise variables to population baseline (BASELINE visit) mean and SD within each cohort---####

#Clear old variables
rm(results)
rm(nm)
rm(visitid, varname)

#First list the names of the columns containing visitdate
#Using the disease duration columns as this covers all the visit dates
datecols <- PD_only %>% 
  select(contains("disease_duration_onset"))
colnames(datecols)

#Not every assessment is done at every timepoint, so only use the visit IDs that match the measurevar
#select_df is a subset of the main dataframe that just has the assessment of interest at each timepoint
#Normalise using the baseline population mean and SD WITHIN EACH COHORT
normalise_popbaseline_selectvisits <- function(data = NULL, measurevar, select_df) {
  
  #Get names of all the visit date variables
  datecol_names <- names(datecols) 
  
  #Create a new variable with just the visit ID
  visit_ids <- substr(datecol_names,1,nchar(datecol_names)-23) 
  
  #From the dataframe with just the assessment of interest, take the visit IDs
  #This is because not all assessments are done at every visit
  names <- names(select_df)
  select_visitid <- word(names, start = 1, sep = "\\_")
  new_visit_ids <- visit_ids[visit_ids %in% select_visitid]
  
  #Create results dataframe for new normalised variables
  results <- as.data.frame(matrix(ncol = 1+length(new_visit_ids), nrow = nrow(data)))
  results[1] <- data$ID #Put study IDs into first column
  
  names(results) <- c("ID", new_visit_ids)
  
  #Create variable names for each visit
  for (visitno in new_visit_ids) {
    
    #Create variable name
    varname <- paste(visitno, "_", measurevar, sep = "")
    
    #Create baseline variable name
    bl_varname <- paste("BL_", measurevar, sep = "")
    
    #Select data from dataframe using variable name
    df <- data[, which(colnames(data) == varname[1] | colnames(data) == "ID")]
    
    #Rename variable name
    df <- df %>% 
      rename(original_var = paste(visitno, "_", measurevar, sep = ""))
    
    #Select baseline variable to calculate population baseline mean  and SD within each cohort
    bl_values <- data %>% 
      select(ID, cohort, bl_varname) %>% 
      rename(bl_values = bl_varname) %>% 
      group_by(cohort) %>% 
      mutate(bl_mean = mean(bl_values, na.rm = TRUE),
             bl_sd = sd(bl_values, na.rm = TRUE)) %>% 
      select(ID, cohort, bl_mean, bl_sd)
    
    
    #Merge with df of values (cycling through each visit)
    df_merged <- df %>%
      inner_join(bl_values, by = "ID")
    
    #Standardise variable according to population baseline mean and SD
    df_new <- df_merged %>% 
      mutate(newvar = (original_var - bl_mean)/bl_sd)
    
    #Put new scaled variable into results dataframe
    results[names(results) == visitno] <- df_new$newvar
    #Rename results column name
    names(results)[names(results) == visitno] <- paste(visitno, "_", measurevar, "_scaled", sep = "")
  }
  merged <- data %>% 
    left_join(results, by = "ID")
  return(merged)
}

#Use function to normalise UPDRSIII percentage scores to population baseline mean and SD
PD_only <- normalise_popbaseline_selectvisits(PD_only, "UPDRS_III_total_perc", UPDRS3)

#Check means and SDs
PD_only %>% 
  group_by(cohort) %>% 
  summarise(BL_mean = mean(BL_UPDRS_III_total_perc_scaled, na.rm = TRUE),
            BL_sd = sd(BL_UPDRS_III_total_perc_scaled, na.rm = TRUE),
            V01_mean = mean(V01_UPDRS_III_total_perc_scaled, na.rm = TRUE),
            V01_sd = sd(V01_UPDRS_III_total_perc_scaled, na.rm = TRUE),
            V04_mean = mean(V04_UPDRS_III_total_perc_scaled, na.rm = TRUE),
            V04_sd = sd(V04_UPDRS_III_total_perc_scaled, na.rm = TRUE),
            V05_mean = mean(V05_UPDRS_III_total_perc_scaled, na.rm = TRUE),
            V05_sd = sd(V05_UPDRS_III_total_perc_scaled, na.rm = TRUE),
            V08_mean = mean(V08_UPDRS_III_total_perc_scaled, na.rm = TRUE),
            V08_sd = sd(V08_UPDRS_III_total_perc_scaled, na.rm = TRUE),
            V11_mean = mean(V11_UPDRS_III_total_perc_scaled, na.rm = TRUE),
            V11_sd = sd(V11_UPDRS_III_total_perc_scaled, na.rm = TRUE),
            V13_mean = mean(V13_UPDRS_III_total_perc_scaled, na.rm = TRUE),
            V13_sd = sd(V13_UPDRS_III_total_perc_scaled, na.rm = TRUE))

#---Transform to percentages and normalise MDS-UPDRS Part II scores---####

#First transform UPDRSIII scores to a percentage of the maximum score
UPDRS2 <- PD_only %>% 
  select(contains("UPDRS_II_total"))

#For loop - convert UPDRS3 score to percentage (divide by 52 which is the total score)
for (nm in names(UPDRS2)){
  PD_only <- PD_only %>% 
    mutate(perc = PD_only[[nm]]/52 * 100)
  visitid <- word(nm, start = 1, sep = "\\_")
  names(PD_only)[names(PD_only) == "perc"] <- paste(visitid, "UPDRS_II_total_perc", sep = "_")
} 


#Use function to normalise UPDRSIII percentage scores to population baseline mean and SD
PD_only <- normalise_popbaseline_selectvisits(PD_only, "UPDRS_II_total_perc", UPDRS2)

#Check means and SDs
PD_only %>% 
  group_by(cohort) %>% 
  summarise(BL_mean = mean(BL_UPDRS_II_total_perc_scaled, na.rm = TRUE),
            BL_sd = sd(BL_UPDRS_II_total_perc_scaled, na.rm = TRUE),
            V01_mean = mean(V01_UPDRS_II_total_perc_scaled, na.rm = TRUE),
            V01_sd = sd(V01_UPDRS_II_total_perc_scaled, na.rm = TRUE),
            V04_mean = mean(V04_UPDRS_II_total_perc_scaled, na.rm = TRUE),
            V04_sd = sd(V04_UPDRS_II_total_perc_scaled, na.rm = TRUE),
            V05_mean = mean(V05_UPDRS_II_total_perc_scaled, na.rm = TRUE),
            V05_sd = sd(V05_UPDRS_II_total_perc_scaled, na.rm = TRUE),
            V08_mean = mean(V08_UPDRS_II_total_perc_scaled, na.rm = TRUE),
            V08_sd = sd(V08_UPDRS_II_total_perc_scaled, na.rm = TRUE),
            V11_mean = mean(V11_UPDRS_II_total_perc_scaled, na.rm = TRUE),
            V11_sd = sd(V11_UPDRS_II_total_perc_scaled, na.rm = TRUE),
            V13_mean = mean(V13_UPDRS_II_total_perc_scaled, na.rm = TRUE),
            V13_sd = sd(V13_UPDRS_II_total_perc_scaled, na.rm = TRUE))


#---Transform to percentages and normalise Hoehn and Yahr stage scores---####

#First transform Hoehn and Yahr scores to a percentage of the maximum score
hoehn_yahr <- PD_only %>% 
  select(contains("hoehn_and_yahr_stage"))

#For loop - convert UPDRS3 score to percentage (divide by 5 which is the total score)
for (nm in names(hoehn_yahr)){
  PD_only <- PD_only %>% 
    mutate(perc = PD_only[[nm]]/5 * 100)
  visitid <- word(nm, start = 1, sep = "\\_")
  names(PD_only)[names(PD_only) == "perc"] <- paste(visitid, "hoehn_and_yahr_stage_perc", sep = "_")
} 

#Use function to normalise Hoehn and Yahr scores
PD_only <- normalise_popbaseline_selectvisits(data = PD_only, "hoehn_and_yahr_stage_perc", hoehn_yahr)

#Summarise means and SDs of scaled variables
PD_only %>%  
  group_by(cohort) %>% 
  summarise(BL_mean = mean(BL_hoehn_and_yahr_stage_perc_scaled, na.rm = TRUE),
            BL_sd = sd(BL_hoehn_and_yahr_stage_perc_scaled, na.rm = TRUE),
            V1_mean = mean(V01_hoehn_and_yahr_stage_perc_scaled, na.rm = TRUE),
            V1_sd = sd(V01_hoehn_and_yahr_stage_perc_scaled, na.rm = TRUE),
            V05_mean = mean(V05_hoehn_and_yahr_stage_perc_scaled, na.rm = TRUE),
            V05_sd = sd(V05_hoehn_and_yahr_stage_perc_scaled, na.rm = TRUE),
            V08_mean = mean(V08_hoehn_and_yahr_stage_perc_scaled, na.rm = TRUE),
            V08_sd = sd(V08_hoehn_and_yahr_stage_perc_scaled, na.rm = TRUE),
            V9_mean = mean(V09_hoehn_and_yahr_stage_perc_scaled, na.rm = TRUE),
            V9_sd = sd(V09_hoehn_and_yahr_stage_perc_scaled, na.rm = TRUE),
            V11_mean = mean(V11_hoehn_and_yahr_stage_perc_scaled, na.rm = TRUE),
            V11_sd = sd(V11_hoehn_and_yahr_stage_perc_scaled, na.rm = TRUE))


#---Merge MoCA baseline BL and screening SC scores---####
#In PPMI the MOCA is only done at the screening (SC) visit and not at baseline
#Need to merge this with PROBAND+Oxford BL in order to normalise all the scores to population baseline mean and SD
#This should be ok because in long format it will relate to the time from PD diagnosis
#Just need to remember to use the SC_disease_duration_onset for PPMI patients when transforming to long format

#Create new variable for the BL_MOCA_total which is either BL in PROBAND+Oxford, or SC in PPMI
PD_only <- PD_only %>% 
  mutate(BL_MOCA_total_new = ifelse(is.na(BL_MOCA_total), SC_MOCA_total, BL_MOCA_total)) %>% 
  select(-BL_MOCA_total, -SC_MOCA_total) %>% 
  rename(BL_MOCA_total = BL_MOCA_total_new)

#---Invert MoCA scores, transform to percentages and normalise---####

#Select MoCA total score variables
MOCA <- PD_only %>% 
  select(contains("MOCA_total"))

#For loop - invert MoCA scores out of 30 and convert to percentage
for (nm in names(MOCA)){
  PD_only <- PD_only %>% 
    mutate(perc = (30 - PD_only[[nm]])/30 * 100)
  visitid <- word(nm, start = 1, sep = "\\_")
  names(PD_only)[names(PD_only) == "perc"] <- paste(visitid, "MOCA_total_inverse_perc", sep = "_")
} 

#Normalise inverted MoCA total scores
PD_only <- normalise_popbaseline_selectvisits(data = PD_only, "MOCA_total_inverse_perc", MOCA)


PD_only %>% 
  group_by(cohort) %>% 
  summarise(BL_mean = mean(BL_MOCA_total_inverse_perc_scaled, na.rm = TRUE),
            BL_sd = sd(BL_MOCA_total_inverse_perc_scaled, na.rm = TRUE),
            V03_mean = mean(V03_MOCA_total_inverse_perc_scaled, na.rm = TRUE),
            V03_sd = sd(V03_MOCA_total_inverse_perc_scaled, na.rm = TRUE),
            V05_mean = mean(V05_MOCA_total_inverse_perc_scaled, na.rm = TRUE),
            V05_sd = sd(V05_MOCA_total_inverse_perc_scaled, na.rm = TRUE),
            V06_mean = mean(V06_MOCA_total_inverse_perc_scaled, na.rm = TRUE),
            V06_sd = sd(V06_MOCA_total_inverse_perc_scaled, na.rm = TRUE),
            V08_mean = mean(V08_MOCA_total_inverse_perc_scaled, na.rm = TRUE),
            V08_sd = sd(V08_MOCA_total_inverse_perc_scaled, na.rm = TRUE))


#---Invert semantic fluency scores, transform to percentages and normalise---####
#Invert the scoring so that higher values indicate worse semantic fluency
#Use the highest score at baseline as the total maximum

semantic_fluency <- PD_only %>% 
  select(contains("seman_flu_score"))

#Calculate maximum semantic fluency score at baseline within each cohort
max_cohort <- PD_only %>% 
  group_by(cohort) %>% 
  mutate(max_BL_seman_flu_cohort = max(BL_seman_flu_score, na.rm =TRUE)) %>% 
  select(ID, cohort, max_BL_seman_flu_cohort)

#Merge with main dataset
PD_only <-PD_only %>% 
  select(-cohort) %>% 
  left_join(max_cohort, by = "ID")

#For loop - invert semantic fluency scores out of maximum score at baseline within each cohort and convert to percentage
#So that higher scores indicate worse fluency
for (nm in names(semantic_fluency)){
  PD_only <- PD_only %>% 
    mutate(perc = (max_BL_seman_flu_cohort - PD_only[[nm]])/max_BL_seman_flu_cohort * 100)
  
  visitid <- word(nm, start = 1, sep = "\\_")
  names(PD_only)[names(PD_only) == "perc"] <- paste(visitid, "seman_flu_score_inverse_perc", sep = "_")
} 

#Normalise to the population baseline mean and SD
PD_only <- normalise_popbaseline_selectvisits(PD_only, "seman_flu_score_inverse_perc", semantic_fluency)

#Check means and SDs of normalised variable
PD_only %>% 
  group_by(cohort) %>% 
  summarise(BL_mean = mean(BL_seman_flu_score_inverse_perc_scaled, na.rm = TRUE),
            BL_sd = sd(BL_seman_flu_score_inverse_perc_scaled, na.rm = TRUE),
            V01_mean = mean(V01_seman_flu_score_inverse_perc_scaled, na.rm = TRUE),
            V01_sd = sd(V01_seman_flu_score_inverse_perc_scaled, na.rm = TRUE),
            V4_mean = mean(V04_seman_flu_score_inverse_perc_scaled, na.rm = TRUE),
            V4_sd = sd(V04_seman_flu_score_inverse_perc_scaled, na.rm = TRUE),
            V8_mean = mean(V08_seman_flu_score_inverse_perc_scaled, na.rm = TRUE),
            V8_sd = sd(V08_seman_flu_score_inverse_perc_scaled, na.rm = TRUE),
            V10_mean = mean(V10_seman_flu_score_inverse_perc_scaled, na.rm = TRUE),
            V10_sd = sd(V10_seman_flu_score_inverse_perc_scaled, na.rm = TRUE))

#Check means of raw scores
PD_only %>% 
  group_by(cohort) %>% 
  summarise(BL_mean = mean(BL_seman_flu_score, na.rm = TRUE),
            V4_mean = mean(V04_seman_flu_score, na.rm = TRUE),
            V8_mean = mean(V08_seman_flu_score, na.rm = TRUE),
            V10_mean = mean(V10_seman_flu_score, na.rm = TRUE))

#---Transform to percentage and normalise MDS-UPDRS 1.1 scores---####

#Select MDS-UPDRS 1.1
updrs1_1 <- PD_only %>% 
  select(contains("UPDRS_I_1recode"))

#Transform to percentages
#For loop - convert UPDRS1.1 score to percentage and then normalise to baseline mean and SD
for (nm in names(updrs1_1)){
  PD_only <- PD_only %>% 
    mutate(perc = PD_only[[nm]]/4 * 100)
  visitid <- word(nm, start = 1, sep = "\\_")
  names(PD_only)[names(PD_only) == "perc"] <- paste(visitid, "UPDRS_I_1recode_perc", sep = "_")
} 

#Normalise to population baseline mean and SD
PD_only <- normalise_popbaseline_selectvisits(data = PD_only, "UPDRS_I_1recode_perc", updrs1_1)

#Summarise means and SDs of scaled variables
PD_only %>%  
  group_by(cohort) %>% 
  summarise(BL_mean = mean(BL_UPDRS_I_1recode_perc_scaled, na.rm = TRUE),
            BL_sd = sd(BL_UPDRS_I_1recode_perc_scaled, na.rm = TRUE),
            V1_mean = mean(V01_UPDRS_I_1recode_perc_scaled, na.rm = TRUE),
            V1_sd = sd(V01_UPDRS_I_1recode_perc_scaled, na.rm = TRUE),
            V4_mean = mean(V04_UPDRS_I_1recode_perc_scaled, na.rm = TRUE),
            V4_sd = sd(V04_UPDRS_I_1recode_perc_scaled, na.rm = TRUE),
            V5_mean = mean(V05_UPDRS_I_1recode_perc_scaled, na.rm = TRUE),
            V5_sd = sd(V05_UPDRS_I_1recode_perc_scaled, na.rm = TRUE),
            V7_mean = mean(V07_UPDRS_I_1recode_perc_scaled, na.rm = TRUE),
            V7_sd = sd(V07_UPDRS_I_1recode_perc_scaled, na.rm = TRUE),
            V9_mean = mean(V09_UPDRS_I_1recode_perc_scaled, na.rm = TRUE),
            V9_sd = sd(V09_UPDRS_I_1recode_perc_scaled, na.rm = TRUE),
            V10_mean = mean(V10_UPDRS_I_1recode_perc_scaled, na.rm = TRUE),
            V10_sd = sd(V10_UPDRS_I_1recode_perc_scaled, na.rm = TRUE))

########## EXTRACT RANDOM SLOPES FOR EACH PROGRESSION MEASURE ########## 
#---First remove data for the PPMI non-annual visits---####
#As only the annual visits are done in the technically defined OFF state for motor assessments

#Delete all the data for PPMI at V05 and V11
PPMI_annualvisits <- PD_only %>% 
  filter(cohort == "PPMI") %>% 
  mutate(V05_UPDRS_III_total_perc_scaled = NA,
         V11_UPDRS_III_total_perc_scaled = NA,
         V05_UPDRS_II_total_perc_scaled = NA,
         V11_UPDRS_II_total_perc_scaled = NA,
         V05_hoehn_and_yahr_stage_perc_scaled = NA,
         V11_hoehn_and_yahr_stage_perc_scaled = NA,
         V05_MOCA_total_inverse_perc_scaled = NA,
         V11_MOCA_total_inverse_perc_scaled = NA,
         V05_seman_flu_score_inverse_perc_scaled = NA,
         V11_seman_flu_score_inverse_perc_scaled = NA,
         V05_UPDRS_I_1recode_perc_scaled = NA,
         V11_UPDRS_I_1recode_perc_scaled = NA)

#Join back with PROBAND+Oxford
PD_only_PRO_Ox <- PD_only %>% 
  filter(cohort == "PROBAND" | cohort == "Oxford")

PD_only <- rbind(PD_only_PRO_Ox, PPMI_annualvisits)

##---MDS-UPDRS III---####
#---Convert data to long format for MDS-UPDRSIII---####

#Time variable is time from disease onset
dis_dur_vars <- PD_only %>% 
  select(ends_with("disease_duration_onset"))

#Convert to dataframe (this is needed for reshape)
PD_only_df <- as.data.frame(PD_only)

#Function to list the variable of interest names and MATCHING disease_duration_onset variables
#Again because not every assessment is done at every timepoint
#Measurevar is the final percentage transformed and scaled variable in quotes
#EXCLUDING SCREENING VISIT also - BL should be the first visit
#Only for ANNUAL visits for the UPDRS3 because only these are done in the practically defined OFF state
list_vars_for_reshape_annualvisits <- function(data = NULL, measurevar) {
  
  #Get names of all the visit date variables
  disdur_names <- names(dis_dur_vars) 
  
  #Create a new variable with just the visit ID
  visit_ids <- substr(disdur_names,1,nchar(disdur_names)-23) 
  
  #Remove the unscheduled visits (U visits) and other random visit IDs from list
  #ALSO EXCLUDING VISIT 14, 15 and 16 - only up to visit 13
  visit_ids <- visit_ids[visit_ids %in% c("BL", "V04", "V05", "V06", "V08", "V10", "V11", "V12", "V13")]
  
  #From the dataframe with just the assessment of interest, take the visit IDs
  #This is because not all assessments are done at every visit
  select_df <- PD_only %>% 
    select(contains(paste(measurevar)))
  
  names <- names(select_df)
  select_visitid <- word(names, start = 1, sep = "\\_")
  new_visit_ids <- visit_ids[visit_ids %in% select_visitid]
  
  #Now write the final vectors which will be used for reshape
  disduration <- c("disease_duration_onset")
  start_names <- paste(new_visit_ids, disduration, sep = "_")
  value_names <- paste(new_visit_ids, measurevar, sep = "_")
  return(list("start_names" = start_names, "value_names" = value_names))
}


#Get variable names for UPDRSIII - only ANNUAL visits
UPDRSIII_names <- list_vars_for_reshape_annualvisits(PD_only_longitudinal, "UPDRS_III_total_perc_scaled")

#Reshape to long format
#Gather both UPDRSIII values and time from BL values
#Using normalised variables
PD_only_UPDRSIII_df_long <- reshape(PD_only_df, idvar="ID", direction="long", 
                                    varying = list(UPDRSIII_names$start_names, UPDRSIII_names$value_names),
                                    v.names = c("time_from_onset", "UPDRSIII_total_perc_scaled"))

#Arrange by ID and time from onset
PD_only_UPDRSIII_df_long <- PD_only_UPDRSIII_df_long %>% 
  arrange(ID, time_from_onset)

check <- PD_only %>% 
  filter(ID == "3001") %>% 
  select(contains("disease_duration_onset"), contains("UPDRS_III_total_perc_scaled"))

#Make back into a tibble
PD_only_UPDRSIII_df_long <- as_tibble(PD_only_UPDRSIII_df_long) 

PD_only_UPDRSIII_df_long %>% 
  dplyr::select(ID, time_from_onset, time, UPDRSIII_total_perc_scaled)

#---Get random slopes for UPDRSIII regressed off age at onset, gender and time from onset---####

#Fit mixed effects model for motor progression in MDS-UPDRS III
#Assuming correlated slope and intercept
#Including cohort and cohort x time interaction as covariates
lmer_UPDRSIII_aao_gender <-  lmer(UPDRSIII_total_perc_scaled ~ age_onset_imput + age_onset_imput*time_from_onset 
                                  + as.factor(gender) + as.factor(gender)*time_from_onset 
                                  + as.factor(cohort) + as.factor(cohort)*time_from_onset
                                  + (time_from_onset | ID),
                                  data = PD_only_UPDRSIII_df_long,
                                  REML = TRUE)

summary(lmer_UPDRSIII_aao_gender)
Anova(lmer_UPDRSIII_aao_gender)

#Take random slope for each individual
UPDRSIII_randomslope <- ranef(lmer_UPDRSIII_aao_gender)$ID

#Convert to table and convert ID to separate column
UPDRSIII_randomslope_tbl <- tibble::rownames_to_column(UPDRSIII_randomslope, "ID")

#Rename random intercept column
UPDRSIII_randomslope_tbl <- UPDRSIII_randomslope_tbl %>% 
  dplyr::rename(randomslope_UPDRSIIIxtime_from_onset = time_from_onset)

#Merge with main dataset
PD_only <- PD_only %>% 
  left_join(UPDRSIII_randomslope_tbl, by = "ID")

#---Histogram of UPDRSIII random slopes---####

ggplot(data = PD_only, mapping = aes(randomslope_UPDRSIIIxtime_from_onset)) +
  geom_histogram(colour = "black", fill = "white") +
  theme_bw() +
  ggsave("plots/PROBAND_Oxford_PPMI_combined/histogram_randomslope_UPDRSIIIxtime_from_onset.png")

##---MOCA---####
#---Convert data to long format for MOCA---####
#The MoCA has been done at BL for PROBAND+Oxford, but at SC for PPMI
#Need to take the appropriate disease_duration_onset for each cohort

PD_only <- PD_only %>% 
  mutate(MOCA_newBL_disease_duration_onset = ifelse(cohort == "PROBAND" | cohort == "Oxford", BL_disease_duration_onset,
                                                    ifelse(cohort == "PPMI", SC_disease_duration_onset, NA)))


#Function to list the variable of interest names and disease duration onset variables
#Only for ANNUAL visits
#Excluding the baseline and screening visits because the MOCA is not done at BL
#Will add in new BL variable for MOCA when we transform to long format
list_vars_for_reshape_annualvisits_SC <- function(data = NULL, measurevar) {
  
  #Get names of all the visit date variables
  disdur_names <- names(dis_dur_vars) 
  
  #Create a new variable with just the visit ID
  visit_ids <- substr(disdur_names,1,nchar(disdur_names)-23) 
  
  #Remove the unscheduled visits (U visits) and other random visit IDs from list
  #ALSO EXCLUDING VISIT 14, 15 and 16 - only up to visit 13
  visit_ids <- visit_ids[visit_ids %in%  c("V04", "V05", "V06", "V08", "V10", "V11", "V12", "V13")]
  
  #From the dataframe with just the assessment of interest, take the visit IDs
  #This is because not all assessments are done at every visit
  select_df <- PD_only %>% 
    select(contains(paste(measurevar)))
  
  names <- names(select_df)
  select_visitid <- word(names, start = 1, sep = "\\_")
  new_visit_ids <- visit_ids[visit_ids %in% select_visitid]
  
  #Now write the final vectors which will be used for reshape
  disduration <- c("disease_duration_onset")
  start_names <- paste(new_visit_ids, disduration, sep = "_")
  value_names <- paste(new_visit_ids, measurevar, sep = "_")
  return(list("start_names" = start_names, "value_names" = value_names))
}

MOCA_names <- list_vars_for_reshape_annualvisits_SC(PD_only, "MOCA_total_inverse_perc_scaled")

#Convert main dataset to dataframe (because we created a new variable above that was not included previously)
PD_only_df <- as.data.frame(PD_only)
#View(colnames(PD_only_df))

#Reshape to long format
#Gather both MOCA values and time from BL values
#Using normalised variables
PD_only_MOCA_df_long <- reshape(PD_only_df, idvar="ID", direction="long", 
                                #Adding extra variables into the varying list, so that we use the new BL disease duration
                                #Also use the new BL MOCA score that combines the BL and SC visits from PROBAND and PPMI
                                varying = list(c("MOCA_newBL_disease_duration_onset", MOCA_names$start_names), 
                                               c("BL_MOCA_total_inverse_perc_scaled", MOCA_names$value_names)),
                                v.names = c("time_from_onset", "MOCA_total_inverse_perc_scaled"))

#Arrange by ID
PD_only_MOCA_df_long <- PD_only_MOCA_df_long %>% 
  arrange(ID, time_from_onset)

#Make back into a tibble
PD_only_MOCA_df_long <- as_tibble(PD_only_MOCA_df_long) 

PD_only_MOCA_df_long %>% 
  dplyr::select(ID, time_from_onset, time, MOCA_total_inverse_perc_scaled)


#---Get random slopes for MOCA adjusted regressed off age at onset, gender, education, and time from onset---####

#Fit mixed effects model for cognitive progression in MOCA
#Using raw education rata (age left school and higher education)
#Using inverse MoCA score (number of incorrect responses)
lmer_MOCA_aao_gender_educ2 <-  lmer(MOCA_total_inverse_perc_scaled ~ age_onset_imput + age_onset_imput*time_from_onset 
                                    + as.factor(gender) + as.factor(gender)*time_from_onset 
                                    + ed_years_before_further_ed + ed_years_before_further_ed*time_from_onset 
                                    + higher_education + higher_education*time_from_onset 
                                    + as.factor(cohort) + as.factor(cohort)* time_from_onset 
                                    + (time_from_onset | ID),
                                    data = PD_only_MOCA_df_long,
                                    REML = TRUE)

summary(lmer_MOCA_aao_gender_educ2)
Anova(lmer_MOCA_aao_gender_educ2)

#Take random slope for each individual
MOCA_randomslope <- ranef(lmer_MOCA_aao_gender_educ2)$ID

#Convert to table and convert ID to separate column
MOCA_randomslope_tbl <- tibble::rownames_to_column(MOCA_randomslope, "ID")

#Rename random intercept column
MOCA_randomslope_tbl <- MOCA_randomslope_tbl %>% 
  dplyr::rename(randomslope_MOCAxtime_from_onset = time_from_onset)

#Merge with main dataset
PD_only <- PD_only %>% 
  left_join(MOCA_randomslope_tbl, by = "ID")

#---Histogram of MOCA random slopes---####

ggplot(data = PD_only, mapping = aes(randomslope_MOCAxtime_from_onset)) +
  geom_histogram(colour = "black", fill = "white") +
  theme_bw() +
  ggsave("plots/PROBAND_Oxford_PPMI_combined/histogram_randomslope_MOCAxtime_from_onset.png")


#---Test normality of random slope distributions---####

#UPDRSIII random slope distribution
shapiro.test(PD_only$randomslope_UPDRSIIIxtime_from_onset)

#QQ plot of UPDRSIII random slope
ggqqplot(PD_only$randomslope_UPDRSIIIxtime_from_onset, ylab = "UPDRSIII random slope")

#MOCA random slope distribution
shapiro.test(PD_only$randomslope_MOCAxtime_from_onset)

#QQ plot of MOCA random slope
ggqqplot(PD_only$randomslope_MOCAxtime_from_onset, ylab = "MOCA random slope")

#---Plot random slopes by cohort---####

#Plot density plot by cohort
ggplot(data = PD_only, mapping = aes(x = randomslope_MOCAxtime_from_onset, colour = cohort)) +
  geom_density()

#Plot boxplots
ggplot(PD_only, aes(x = cohort, y = randomslope_MOCAxtime_from_onset, fill = cohort)) +
  geom_boxplot()

#Plot histogram
ggplot(PD_only, aes(x=randomslope_MOCAxtime_from_onset, fill=cohort)) +
  geom_histogram(alpha=.5, position = "identity")

##---MDS-UPDRS II---####
#---Convert data to long format for MDS-UPDRSII---####

#Get variable names
UPDRSII_names <- list_vars_for_reshape_annualvisits(PD_only, "UPDRS_II_total_perc_scaled")

#Reshape to long format
#Gather both UPDRSIII values and time from BL values
#Using normalised variables
PD_only_UPDRSII_df_long <- reshape(PD_only_df, idvar="ID", direction="long", 
                                   varying = list(UPDRSII_names$start_names, UPDRSII_names$value_names),
                                   v.names = c("time_from_onset", "UPDRSII_total_perc_scaled"))

#Arrange by ID
PD_only_UPDRSII_df_long <- PD_only_UPDRSII_df_long %>% 
  arrange(ID, time_from_onset)

#Make back into a tibble
PD_only_UPDRSII_df_long <- as_tibble(PD_only_UPDRSII_df_long) 

PD_only_UPDRSII_df_long %>% 
  dplyr::select(ID, time_from_onset, time, UPDRSII_total_perc_scaled)

#---Get random slopes for UPDRSII regressed off age at onset, gender and time from onset---####

#Fit mixed effects model for motor progression in MDS-UPDRS III
#Assuming correlated slope and intercept
lmer_UPDRSII_aao_gender <-  lmer(UPDRSII_total_perc_scaled ~ age_onset_imput + age_onset_imput*time_from_onset 
                                 + as.factor(gender) + as.factor(gender)*time_from_onset 
                                 + as.factor(cohort) + as.factor(cohort)* time_from_onset
                                 + (time_from_onset | ID),
                                 data = PD_only_UPDRSII_df_long,
                                 REML = TRUE)

summary(lmer_UPDRSII_aao_gender)
Anova(lmer_UPDRSII_aao_gender)

#Take random slope for each individual
UPDRSII_randomslope <- ranef(lmer_UPDRSII_aao_gender)$ID

#Convert to table and convert ID to separate column
UPDRSII_randomslope_tbl <- tibble::rownames_to_column(UPDRSII_randomslope, "ID")

#Rename random intercept column
UPDRSII_randomslope_tbl <- UPDRSII_randomslope_tbl %>% 
  dplyr::rename(randomslope_UPDRSIIxtime_from_onset = time_from_onset)

#Merge with main dataset
PD_only <- PD_only %>% 
  left_join(UPDRSII_randomslope_tbl, by = "ID")

#---Histogram of UPDRSII random slopes---####

ggplot(data = PD_only, mapping = aes(randomslope_UPDRSIIxtime_from_onset)) +
  geom_histogram(colour = "black", fill = "white") +
  theme_bw() +
  ggsave("plots/PROBAND_Oxford_PPMI_combined/histogram_randomslope_UPDRSIIxtime_from_onset.png")

##---HOEHN AND YAHR STAGE---####
#---Convert data to long format for Hoehn and Yahr stage---####

#Get variable names
hoehnyahr_names <- list_vars_for_reshape_annualvisits(PD_only, "hoehn_and_yahr_stage_perc_scaled")

#Reshape to long format
#Gather both UPDRSIII values and time from PD onset values
#Using normalised variables
PD_only_hy_df_long <- reshape(PD_only_df, idvar="ID", direction="long", 
                              varying = list(hoehnyahr_names$start_names, hoehnyahr_names$value_names),
                              v.names = c("time_from_onset", "hoehn_and_yahr_stage_perc_scaled"))

#Arrange by ID
PD_only_hy_df_long <- PD_only_hy_df_long %>% 
  arrange(ID, time_from_onset)

#Make back into a tibble
PD_only_hy_df_long <- as_tibble(PD_only_hy_df_long) 

PD_only_hy_df_long %>% 
  dplyr::select(ID, time_from_onset, hoehn_and_yahr_stage_perc_scaled)

#---Get random slopes for Hoehn and Yahr regressed off age at onset, gender and time from onset---####

#Fit mixed effects model for motor progression in H&Y stage
#Assuming correlated slope and intercept
lmer_HY_aao_gender <-  lmer(hoehn_and_yahr_stage_perc_scaled ~ age_onset_imput + age_onset_imput*time_from_onset 
                            + as.factor(gender) + as.factor(gender)*time_from_onset 
                            + as.factor(cohort) + as.factor(cohort)* time_from_onset
                            + (time_from_onset | ID),
                            data = PD_only_hy_df_long,
                            REML = TRUE)

summary(lmer_HY_aao_gender)
Anova(lmer_HY_aao_gender)

#Take random slope for each individual
HY_randomslope <- ranef(lmer_HY_aao_gender)$ID

#Convert to table and convert ID to separate column
HY_randomslope_tbl <- tibble::rownames_to_column(HY_randomslope, "ID")

#Rename random intercept column
HY_randomslope_tbl <- HY_randomslope_tbl %>% 
  dplyr::rename(randomslope_HYxtime_from_onset = time_from_onset)

#Merge with main dataset
PD_only <- PD_only %>% 
  left_join(HY_randomslope_tbl, by = "ID")

#---Histogram of HY random slopes---####

ggplot(data = PD_only, mapping = aes(randomslope_HYxtime_from_onset)) +
  geom_histogram(colour = "black", fill = "white") +
  theme_bw() +
  ggsave("plots/PROBAND_Oxford_PPMI_combined/histogram_randomslope_HYxtime_from_onset.png")

##---SEMANTIC FLUENCY---####
#---Convert data to long format for fluency---####

#Get variable names
semanflu_names <- list_vars_for_reshape_annualvisits(PD_only, "seman_flu_score_inverse_perc_scaled")

#Reshape to long format
#Gather both UPDRSIII values and time from PD onset values
#Using normalised variables
PD_only_semanflu_df_long <- reshape(PD_only_df, idvar="ID", direction="long", 
                                    varying = list(semanflu_names$start_names, semanflu_names$value_names),
                                    v.names = c("time_from_onset", "seman_flu_score_inverse_perc_scaled"))

#Arrange by ID
PD_only_semanflu_df_long <- PD_only_semanflu_df_long %>% 
  arrange(ID, time_from_onset)

#Make back into a tibble
PD_only_semanflu_df_long <- as_tibble(PD_only_semanflu_df_long) 

PD_only_semanflu_df_long %>% 
  dplyr::select(ID, time_from_onset, seman_flu_score_inverse_perc_scaled)

#---Get random slopes for fluency regressed off age at onset, gender and time from onset---####

#Fit mixed effects model for progression in fluency
#Assuming correlated slope and intercept
#Using 2 raw education variables
lmer_fluency_aao_gender2 <-  lmer(seman_flu_score_inverse_perc_scaled ~ age_onset_imput + age_onset_imput*time_from_onset 
                                  + as.factor(gender) + as.factor(gender)*time_from_onset 
                                  + ed_years_before_further_ed + ed_years_before_further_ed*time_from_onset 
                                  + higher_education + higher_education*time_from_onset 
                                  + as.factor(cohort) + as.factor(cohort)* time_from_onset 
                                  + (time_from_onset | ID),
                                  data = PD_only_semanflu_df_long,
                                  REML = TRUE)

summary(lmer_fluency_aao_gender2)
Anova(lmer_fluency_aao_gender2)

#Take random slope for each individual
fluency_randomslope <- ranef(lmer_fluency_aao_gender2)$ID

#Convert to table and convert ID to separate column
fluency_randomslope_tbl <- tibble::rownames_to_column(fluency_randomslope, "ID")

#Rename random intercept column
fluency_randomslope_tbl <- fluency_randomslope_tbl %>% 
  dplyr::rename(randomslope_fluencyxtime_from_onset = time_from_onset)

#Merge with main dataset
PD_only <- PD_only %>% 
  left_join(fluency_randomslope_tbl, by = "ID")

#---Histogram of fluency random slopes---####

ggplot(data = PD_only, mapping = aes(randomslope_fluencyxtime_from_onset)) +
  geom_histogram(colour = "black", fill = "white") +
  theme_bw() +
  ggsave("plots/PROBAND_Oxford_PPMI_combined/histogram_randomslope_fluencyxtime_from_onset.png")

#---Test for normality---####

#Test whether random slopes for semantic fluency are normally distributed
shapiro.test(PD_only$randomslope_fluencyxtime_from_onset)

#QQ plot of fluency random slopes
ggqqplot(PD_only$randomslope_fluencyxtime_from_onset, ylab = "randomslope_fluency")


#---Plot random slope scores by cohort---####

#Summarise means and medians and SDs of semantic fluency random slope scores by cohort
PD_only %>% 
  group_by(cohort) %>% 
  summarise(mean = mean(randomslope_fluencyxtime_from_onset, na.rm = TRUE),
            median = median(randomslope_fluencyxtime_from_onset, na.rm = TRUE),
            sd = sd(randomslope_fluencyxtime_from_onset, na.rm = TRUE))

#Plot density plot by cohort
ggplot(data = PD_only, mapping = aes(x = randomslope_fluencyxtime_from_onset, colour = cohort)) +
  geom_density()

#Plot histogram
ggplot(PD_only, aes(x=randomslope_fluencyxtime_from_onset, fill=cohort)) +
  geom_histogram(alpha=.5, position = "identity")

#Plot boxplots
ggplot(PD_only, aes(x = cohort, y = randomslope_fluencyxtime_from_onset, fill = cohort)) +
  geom_boxplot()

##---MDS-UPDRS 1.1 MEMORY---####
#---Convert data to long format for MDS-UPDRS 1.1---####

#Get variable names
UPDRS1_1_names <- list_vars_for_reshape_annualvisits(PD_only, "UPDRS_I_1recode_perc_scaled")

#Reshape to long format
#Gather both UPDRSIII values and time from PD onset values
#Using normalised variables
PD_only_UPDRSI1_df_long <- reshape(PD_only_df, idvar="ID", direction="long", 
                                   varying = list(UPDRS1_1_names$start_names, UPDRS1_1_names$value_names),
                                   v.names = c("time_from_onset", "UPDRS_I_1_perc_scaled"))

#Arrange by ID
PD_only_UPDRSI1_df_long <- PD_only_UPDRSI1_df_long %>% 
  arrange(ID, time_from_onset)

#Make back into a tibble
PD_only_UPDRSI1_df_long <- as_tibble(PD_only_UPDRSI1_df_long) 

PD_only_UPDRSI1_df_long %>% 
  dplyr::select(ID, time_from_onset, UPDRS_I_1_perc_scaled)

PD_only %>% 
  select(BL_UPDRS_I_1recode_perc_scaled)

#---Get random slopes for MDS-UPDRS1.1 regressed off age at onset, gender and time from onset---####

#Fit mixed effects model for cognitive progression in MDS-UPDRS 1.1
#Including 2 raw education variables as covars
lmer_UPDRSI_1_aao_gender2 <-  lmer(UPDRS_I_1_perc_scaled ~ age_onset_imput + age_onset_imput*time_from_onset 
                                   + as.factor(gender) + as.factor(gender)*time_from_onset 
                                   + ed_years_before_further_ed + ed_years_before_further_ed*time_from_onset 
                                   + higher_education + higher_education*time_from_onset 
                                   + as.factor(cohort) + as.factor(cohort)* time_from_onset
                                   + (time_from_onset | ID),
                                   data = PD_only_UPDRSI1_df_long,
                                   REML = TRUE)

summary(lmer_UPDRSI_1_aao_gender2)
Anova(lmer_UPDRSI_1_aao_gender2)

#Take random slope for each individual
UPDRSI_1_randomslope <- ranef(lmer_UPDRSI_1_aao_gender2)$ID

#Convert to table and convert ID to separate column
UPDRSI_1_randomslope_tbl <- tibble::rownames_to_column(UPDRSI_1_randomslope, "ID")

#Rename random intercept column
UPDRSI_1_randomslope_tbl <- UPDRSI_1_randomslope_tbl %>% 
  dplyr::rename(randomslope_UPDRSI1xtime_from_onset = time_from_onset)

#Merge with main dataset
PD_only <- PD_only %>% 
  left_join(UPDRSI_1_randomslope_tbl, by = "ID")

#---Histogram of MDS-UPDRS 1.1 random slopes---####

ggplot(data = PD_only, mapping = aes(randomslope_UPDRSI1xtime_from_onset)) +
  geom_histogram(colour = "black", fill = "white") +
  theme_bw() +
  ggsave("plots/PROBAND_Oxford_PPMI_combined/histogram_randomslope_UPDRSI1xtime_from_onset2.png")

########## MOTOR VS. COGNITIVE PROGRESSION ########## 
#---Correlation between MoCA and UPDRSIII random slopes---####

#Pearson correlation between UPDRSIII and MOCA random slopes
#Use Pearson correlation as this is what Mike suggested - Spearman only for ordinal data?
UPDRSIII_MOCA_randomslope_spearman_cor <- cor.test(PD_only$randomslope_UPDRSIIIxtime_from_onset, PD_only$randomslope_MOCAxtime_from_onset,
                                                   use = "complete.obs",
                                                   method = "pearson")

UPDRSIII_MOCA_randomslope_spearman_cor

#---Scatterplot of MOCA vs UPDRSIII random slopes---####

#Scatterplot of MOCA random slope vs. UPDRSIII random slope
#Line is linear model
ggplot(data = PD_only, mapping = aes(x = randomslope_MOCAxtime_from_onset, y = randomslope_UPDRSIIIxtime_from_onset)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()

#Scatterplot with Loess smooth line
ggplot(data = PD_only, mapping = aes(x = randomslope_MOCAxtime_from_onset, y = randomslope_UPDRSIIIxtime_from_onset)) +
  geom_point() +
  geom_smooth() +
  theme_bw()

#---PCA of motor progression measures---####
#Includes: MDS-UPDRS III, MDS-UPDRS II, Hoehn and Yahr stage

#Select random slope variables for motor progression measures
PD_only_motorprogression <- PD_only %>% 
  select(randomslope_UPDRSIIIxtime_from_onset, randomslope_UPDRSIIxtime_from_onset, randomslope_HYxtime_from_onset)

rownames(PD_only_motorprogression) <- PD_only$ID

#Run PCA on random slopes of motor measures
motor_progression.pca <- prcomp(~ ., data = PD_only_motorprogression, center = TRUE, scale. = TRUE)

summary(motor_progression.pca)

#Create scree plot
fviz_eig(motor_progression.pca)
screeplot(motor_progression.pca, type = "lines")

# Variability of each principal component: pr.var (these are the eigenvalues)
#sdev is the square root of the eigenvalues
motor_progression.var <- (motor_progression.pca$sdev)^2

# Variance explained by each principal component: pve
motor_progression.pve <- motor_progression.var / sum(motor_progression.var)
motor_progression.pve

#Plot proportion of variance explained
png(filename = "plots/PROBAND_Oxford_PPMI_combined/scree_motor_progression_varexplained.png")
plot(motor_progression.pve, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     main = "Motor progression PCA: Variance explained",
     ylim = c(0, 1), type = "b")
dev.off()

#Plot the eigenvalues
png(filename = "plots/PROBAND_Oxford_PPMI_combined/scree_motor_progression_eigenvals.png")
plot(motor_progression.var, xlab = "Principal Component",
     ylab = "Eigenvalue",
     main = "Motor progression PCA: Scree plot",
     type = "b")
dev.off()

#Get individual results for PCA
motor_progression_ind.coord <- motor_progression.pca$x

#Histogram plot of PC1
ggplot(data = as.data.frame(motor_progression_ind.coord), mapping = aes(x = PC1)) +
  geom_histogram(colour = "black", fill = "white") +
  theme_bw()

#Convert to tibble and include individual IDs as a column
motor_progression_ind.coord_tbl <- tibble::rownames_to_column(as.data.frame(motor_progression_ind.coord), "ID")

#Rename columns
motor_progression_ind.coord_tbl <- motor_progression_ind.coord_tbl %>% 
  rename(motorprogression_PC1 = "PC1",
         motorprogression_PC2 = "PC2",
         motorprogression_PC3 = "PC3")

#Merge with main dataframe
PD_only <- PD_only %>% 
  left_join(motor_progression_ind.coord_tbl, by = "ID")

#---Correlogram of motor progression measures---####

cor_motor <- cor(PD_only_motorprogression, method = "pearson", use = "pairwise.complete.obs")

png(file = "plots/PROBAND_Oxford_PPMI_combined/correlogram_motorprog.png")
corrplot(cor_motor, method= "circle", type = "upper",
         addCoef.col = "black")
dev.off()

#---PCA of cognitive progression measures---####
#Includes: MOCA total score (unadjusted), phonemic fluency, MDS-UPDRS 1.1

PD_only_cogprogression <- PD_only %>% 
  select(randomslope_MOCAxtime_from_onset, randomslope_fluencyxtime_from_onset, randomslope_UPDRSI1xtime_from_onset)

rownames(PD_only_cogprogression) <- PD_only$ID

#Run PCA on random slopes of motor measures
cog_progression.pca <- prcomp(~ ., data = PD_only_cogprogression, center = TRUE, scale. = TRUE)

summary(cog_progression.pca)

#Create scree plot
fviz_eig(cog_progression.pca)

# Variability of each principal component: pr.var (these are the eigenvalues)
#sdev is the square root of the eigenvalues
cog_progression.var <- (cog_progression.pca$sdev)^2

# Variance explained by each principal component: pve
cog_progression.pve <- cog_progression.var / sum(cog_progression.var)
cog_progression.pve

#Plot proportion of variance explained
png(filename = "plots/PROBAND_Oxford_PPMI_combined/scree_cog_progression_varexplained.png")
plot(cog_progression.pve, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     main = "Cognitive progression PCA: Variance explained",
     ylim = c(0, 1), type = "b")
dev.off()

#Plot the eigenvalues
png(filename = "plots/PROBAND_Oxford_PPMI_combined/scree_cog_progression_eigenvals.png")
plot(cog_progression.var, xlab = "Principal Component",
     ylab = "Eigenvalue",
     main = "Cognitive progression PCA: Scree plot",
     type = "b")
dev.off()

#Get individual results for PCA
cog_progression_ind.coord <- cog_progression.pca$x

#Histogram plot of PC1
ggplot(data = as.data.frame(cog_progression_ind.coord), mapping = aes(x = PC1)) +
  geom_histogram(colour = "black", fill = "white") +
  theme_bw()

#Convert to tibble and include individual IDs as a column
cog_progression_ind.coord_tbl <- tibble::rownames_to_column(as.data.frame(cog_progression_ind.coord), "ID")

#Rename columns
cog_progression_ind.coord_tbl <- cog_progression_ind.coord_tbl %>% 
  rename(cogprogression_PC1 = "PC1",
         cogprogression_PC2 = "PC2",
         cogprogression_PC3 = "PC3")

#Merge with main dataframe
PD_only <- PD_only %>% 
  left_join(cog_progression_ind.coord_tbl, by = "ID")

#---Correlogram of cognitive progression measures---####

#Correlation matrix
cor_cognitive <- cor(PD_only_cogprogression, method = "pearson", use = "pairwise.complete.obs")

#Plot correlogram
png(file = "plots/PROBAND_Oxford_PPMI_combined/correlogram_cogprog.png")
corrplot(cor_cognitive, method= "circle", type = "upper",
         addCoef.col = "black")
dev.off()

#---Correlation between motor and cognitive PCA first principal components---####

#Pearson correlation between first PCs for motor and cognitive progression
cor_motor_cognitive <- cor.test(PD_only$motorprogression_PC1, PD_only$cogprogression_PC1, method = "pearson")

cor_motor_cognitive

#---Correlation matrix between all principal components for motor and cognitive progression---####

PD_only_progression_PCs <- PD_only %>% 
  select(motorprogression_PC1, motorprogression_PC2, motorprogression_PC3,
         cogprogression_PC1, cogprogression_PC2, cogprogression_PC3)

progression_cor_matrix <- rcorr(as.matrix(PD_only_progression_PCs), type = "pearson")$r

progression_cor_matrix_pval <- rcorr(as.matrix(PD_only_progression_PCs), type = "pearson")$P

rsquared <- progression_cor_matrix*progression_cor_matrix

#---Scatterplot of motor and cognitive progression first principal components---####

#Scatterplot of motor progression PC1 vs. cognitive progression PC1
ggplot(data = PD_only, mapping = aes(x = motorprogression_PC1, y = cogprogression_PC1)) +
  geom_point(alpha = 1/4) +
  geom_smooth(method = "lm") +
  theme_bw() + 
  ggsave("plots/PROBAND_Oxford_PPMI_combined/scatter_motorPC1_cogPC1.png")


#---Correlation between all principal components for motor progression and study site---####

PD_only_motorprogression_PCs_cohort <- PD_only %>% 
  select(ID, cohort, motorprogression_PC1, motorprogression_PC2, motorprogression_PC3)

PD_only_motorprogression_PCs_cohort <- as.data.frame(PD_only_motorprogression_PCs_cohort)

PD_only_motorprogression_PCs_cohort_melt <- reshape2::melt(PD_only_motorprogression_PCs_cohort, measure.vars = 3:5)

PD_only_motorprogression_PCs_cohort_melt <- PD_only_motorprogression_PCs_cohort_melt %>% 
  arrange(ID)

#Plot boxplot of motor principal components vs. cohort
ggplot(PD_only_motorprogression_PCs_cohort_melt, aes(x = variable, y = value, fill = cohort)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#---Correlation between all principal components for cognitive progression and study site---####

PD_only_cogprogression_PCs_cohort <- PD_only %>% 
  select(ID, cohort, cogprogression_PC1, cogprogression_PC2, cogprogression_PC3)

PD_only_cogprogression_PCs_cohort <- as.data.frame(PD_only_cogprogression_PCs_cohort)

PD_only_cogprogression_PCs_cohort_melt <- reshape2::melt(PD_only_cogprogression_PCs_cohort, measure.vars = 3:5)

PD_only_cogprogression_PCs_cohort_melt <- PD_only_cogprogression_PCs_cohort_melt %>% 
  arrange(ID)

#Plot boxplot of motor principal components vs. cohort
ggplot(PD_only_cogprogression_PCs_cohort_melt, aes(x = variable, y = value, fill = cohort)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#---Association between composite progression scores and cohorts using glm---####

#Get the names of the motor and cognitive progression principal components
motor_cog_progression_names <- names(PD_only %>% 
                                       select(contains("motorprogression_PC"), contains("cogprogression_PC")))

#Make results dataframe
results_glm_progression_cohort <- as.data.frame(matrix(ncol = 6, nrow = 2))
names(results_glm_progression_cohort) <- motor_cog_progression_names
rownames(results_glm_progression_cohort) <- c("Chiqs", "Pval")

#For loop to run glm of each all progression principal component vs. cohort
#Output chisq and p value results into results data frame
for (nm in motor_cog_progression_names) {
  glm <- glm(PD_only[[nm]] ~ cohort, data = PD_only)
  results_glm_progression_cohort[names(results_glm_progression_cohort) == nm][1,] <- Anova(glm)[1]
  results_glm_progression_cohort[names(results_glm_progression_cohort) == nm][2,] <- Anova(glm)[3]
}


########## COMBINED PROGRESSION SCORE ########## 
#---PCA of all progression measures combined---####
#Includes: MDS-UPDRS III, MDS-UPDRS II, Hoehn and Yahr stage

#Select random slope variables for motor progression measures
PD_only_allprogression <- PD_only %>% 
  select(randomslope_UPDRSIIIxtime_from_onset, randomslope_UPDRSIIxtime_from_onset, randomslope_HYxtime_from_onset,
         randomslope_MOCAxtime_from_onset, randomslope_fluencyxtime_from_onset, randomslope_UPDRSI1xtime_from_onset)

rownames(PD_only_allprogression) <- PD_only$ID

#Run PCA on random slopes of motor measures
all_progression.pca <- prcomp(~ ., data = PD_only_allprogression, center = TRUE, scale. = TRUE)

summary(all_progression.pca)

#Create scree plot
fviz_eig(all_progression.pca)

# Variability of each principal component: pr.var (these are the eigenvalues)
#sdev is the square root of the eigenvalues
all_progression.var <- (all_progression.pca$sdev)^2

# Variance explained by each principal component: pve
all_progression.pve <- all_progression.var / sum(all_progression.var)
all_progression.pve

#Plot proportion of variance explained
png(filename = "plots/PROBAND_Oxford_PPMI_combined/scree_all_progression_varexplained.png")
plot(all_progression.pve, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     main = "Combined progression PCA: Variance explained",
     ylim = c(0, 1), type = "b")
dev.off()

#Plot the eigenvalues
png(filename = "plots/PROBAND_Oxford_PPMI_combined/scree_all_progression_eigenvals.png")
plot(all_progression.pve, xlab = "Principal Component",
     ylab = "Eigenvalue",
     main = "Combined progression PCA: Scree plot",
     type = "b")
dev.off()

#Get individual results for PCA
all_progression_ind.coord <- all_progression.pca$x

all_progression_ind.coord <- as.data.frame(all_progression_ind.coord)

#Histogram plot of PC1
ggplot(data = all_progression_ind.coord, mapping = aes(x = PC1)) +
  geom_histogram(colour = "black", fill = "white") +
  theme_bw()

#Convert to tibble and include individual IDs as a column
all_progression_ind.coord_tbl <- tibble::rownames_to_column(as.data.frame(all_progression_ind.coord), "ID")

#Rename columns
all_progression_ind.coord_tbl <- all_progression_ind.coord_tbl %>% 
  rename(allprogression_PC1 = "PC1",
         allprogression_PC2 = "PC2",
         allprogression_PC3 = "PC3",
         allprogression_PC4 = "PC4",
         allprogression_PC5 = "PC5",
         allprogression_PC6 = "PC6")

#Merge with main dataframe
PD_only <- PD_only %>% 
  left_join(all_progression_ind.coord_tbl, by = "ID")

PD_only %>% 
  select(allprogression_PC1)

#---Distribution of PC1 for combined progression---####

#Density plot with normal distribution of PC1
ggplot(data = PD_only, aes(x = allprogression_PC1)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", bins = 20) +
  theme_bw() +
  stat_function(fun = dnorm, 
                color = "red",
                args = list(mean = mean(PD_only$allprogression_PC1, na.rm = TRUE), sd = sd(PD_only$allprogression_PC1, na.rm = TRUE)))


#Regular histogram of PC1 distribution
ggplot(data = PD_only, aes(x = allprogression_PC1)) +
  geom_histogram(colour = "black", fill = "white", bins = 20) +
  theme_bw() +
  ggsave("plots/PROBAND_Oxford_PPMI_combined/histogram_allprogression_PC1.png")

#PC1 test for normality
shapiro.test(PD_only$allprogression_PC1)

#QQ plot of PC1
ggqqplot(PD_only$allprogression_PC1, ylab = "PC1")

#---Correlation between combined progression PCs and raw measures---####

#Select combined progression PCs and random slope variables for all the original measures
PD_only_allprogression_PCs_rawvars <- PD_only %>% 
  select(randomslope_UPDRSIIIxtime_from_onset, randomslope_UPDRSIIxtime_from_onset, randomslope_HYxtime_from_onset,
         randomslope_MOCAxtime_from_onset, randomslope_fluencyxtime_from_onset, randomslope_UPDRSI1xtime_from_onset,
         allprogression_PC1, allprogression_PC2, allprogression_PC3, allprogression_PC4, allprogression_PC5, allprogression_PC6)

allprogression_cor_matrix <- rcorr(as.matrix(PD_only_allprogression_PCs_rawvars), type = "pearson")$r

write.table(allprogression_cor_matrix, "outputs/allprogression_cor_matrix.txt", quote = F, row.names = T)

allprogression_cor_matrix_pval <- rcorr(as.matrix(PD_only_allprogression_PCs_rawvars), type = "pearson")$P

#---Correlation between all principal components for motor progression and study site---####

PD_only_allprogression_PCs_cohort <- PD_only %>% 
  select(ID, cohort, allprogression_PC1, allprogression_PC2, allprogression_PC3,
         allprogression_PC4, allprogression_PC5, allprogression_PC6)

PD_only_allprogression_PCs_cohort <- as.data.frame(PD_only_allprogression_PCs_cohort)

PD_only_allprogression_PCs_cohort_melt <- reshape2::melt(PD_only_allprogression_PCs_cohort, measure.vars = 3:8)

PD_only_allprogression_PCs_cohort_melt <- PD_only_allprogression_PCs_cohort_melt %>% 
  arrange(ID)

#Plot boxplot of motor principal components vs. cohort
ggplot(PD_only_allprogression_PCs_cohort_melt, aes(x = variable, y = value, fill = cohort)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#---Association between composite progression scores and cohorts using glm---####

#Get the names of the allprogression principal components
allprogression_names <- names(PD_only %>% 
                                select(contains("allprogression_PC")))

#Make results dataframe
results_glm_allprogression_cohort <- as.data.frame(matrix(ncol = 6, nrow = 2))
names(results_glm_allprogression_cohort) <- allprogression_names
rownames(results_glm_allprogression_cohort) <- c("Chiqs", "Pval")

#For loop to run glm of each all progression principal component vs. cohort
#Output chisq and p value results into results data frame
for (nm in allprogression_names) {
  glm <- glm(PD_only[[nm]] ~ cohort, data = PD_only)
  results_glm_allprogression_cohort[names(results_glm_allprogression_cohort) == nm][1,] <- Anova(glm)[1]
  results_glm_allprogression_cohort[names(results_glm_allprogression_cohort) == nm][2,] <- Anova(glm)[3]
}

#---Association between composite progression scores and raw scales using mixed effects models---####
#Mike Nalls suggested getting the r2 of real outcomes per cohort relating to the progression score

#Convert data to df
PD_only_df <- as.data.frame(PD_only)

#For loop to get r2 for each outcome (raw values rather than percentage scaled) for each cohort
#Create results dataframe
results_r2 <- as.data.frame(matrix(ncol = 6, nrow = 3))
names(results_r2) <- c("UPDRS_III_total", "UPDRS_II_total", "hoehn_and_yahr_stage", "MOCA_total", "seman_flu_score", "UPDRS_I_1recode")
rownames(results_r2) <- c("PROBAND", "Oxford", "PPMI")

#For loop
for (nm in names(results_r2)) {
  if (nm == "MOCA_total") {
    
    #For the MOCA, this leaves out the BL visit as the MOCA is actually done at screening in PPMI
    names <- list_vars_for_reshape_annualvisits_SC(data = PD_only, nm)
    
    #Reshape to long format
    #Using raw variables rather than percentage scaled variables
    #For the MOCA, the 
    df_long <- reshape(PD_only_df, idvar="ID", direction="long", 
                       varying = list(c("MOCA_newBL_disease_duration_onset", names$start_names),
                                      c("BL_MOCA_total", names$value_names)),
                       v.names = c("time_from_onset", nm))
    
  } else {
    names <- list_vars_for_reshape_annualvisits(data = PD_only, nm)
    
    #Reshape to long format
    #Gather both UPDRSIII values and time from BL values
    #Using raw variables rather than percentage scaled variables
    df_long <- reshape(PD_only_df, idvar="ID", direction="long", 
                       varying = list(names$start_names, names$value_names),
                       v.names = c("time_from_onset", nm))
    
  }
  
  #Arrange by ID and time from onset
  df_long <- as_tibble(df_long %>% 
                         arrange(ID, time_from_onset) %>% 
                         select(ID, cohort, time_from_onset, all_of(nm), allprogression_PC1))
  
  #For each cohort at a time
  for (cohort_nm in rownames(results_r2)) {
    
    #Filter for just one cohort at a time
    df_cohort <- df_long %>% 
      filter(cohort == cohort_nm)
    
    #Remove rows with missing values of original scale
    df_cohort <- df_cohort %>% 
      filter(!is.na(df_cohort[nm]))
    
    #Fit mixed effects model of raw variable against allprogression PC1
    #Assuming correlated slope and intercept
    lmer <-  lmer(df_cohort[[nm]] ~ allprogression_PC1 + (time_from_onset | ID),
                  data = df_cohort,
                  REML = TRUE)
    
    results_r2[rownames(results_r2) %in% cohort_nm, names(results_r2) == nm] <- r.squaredGLMM(lmer)[1]
  }
}

#Export r-squared results
write.csv(results_r2, "outputs/r2_allprogressionPC1_vs_rawscales.csv", quote = F,
          row.names = T, col.names = T)


########## EXPORT COMPOSITE PROGRESSION SCORES ########## 
#---Merge with PROBAND imputed data fam file---####

#Read in fam file from final imputed dataset
imputed_fam <- read.table("../genetic_data/allchromosomes.filtered.converted_R2_0.7.plink.hardcall_0.4999.fam")

imputed_fam <- imputed_fam %>% 
  separate(V1, c("FID", "IID", "pos")) %>% 
  mutate(IID_new = paste(IID, pos, sep = "_"))

#Merge with clinical dataset
#This is to get the correct IID for the PROBAND samples which are FID_IID
PD_only <- PD_only %>% 
  left_join(imputed_fam, by = c("IID" = "IID_new"))

#---Export phenotype, cohort, batch number and gender data---####

#Count how many indivduals have composite progression score
PD_only %>% 
  filter(!is.na(allprogression_PC1)) %>% 
  summarise(count = n())

#Count how many indivduals have motor progression score
PD_only %>% 
  filter(!is.na(motorprogression_PC1)) %>% 
  summarise(count = n())

#Count how many indivduals have progression score
PD_only %>% 
  filter(!is.na(cogprogression_PC1)) %>% 
  summarise(count = n())

#Select variables to export. Remove individuals who do not have IID (are not in genotype files)
#The Oxford genetic IID is just ID_ID
to_export <- PD_only %>% 
  select(ID, IID, V2, cohort, allprogression_PC1, gender, new_batch_number, array,
         motorprogression_PC1, cogprogression_PC1) %>% 
  filter(!is.na(IID)) %>% #Remove individuals missing IID
  filter(!is.na(new_batch_number)) %>%  #Remove individuals missing batch number
  filter(!is.na(V2) & cohort == "PROBAND" | 
           is.na(V2) & cohort == "Oxford" |
           is.na(V2) & cohort == "PPMI") %>%  #Remove only PROBAND individuals missing V2 (new IID)
  mutate(IID_final = ifelse(!is.na(V2), paste(V2),
                            ifelse(is.na(V2) & cohort == "Oxford", IID,
                                   ifelse(is.na(V2) & cohort == "PPMI", paste("PPMISI", IID, sep =""), NA)))) %>% 
  #Create final IID - for Oxford this is the IID, for PPMI it has a prefix of PPMISI
  select(IID_final, gender, cohort, new_batch_number, array, allprogression_PC1, motorprogression_PC1, cogprogression_PC1) %>% 
  mutate(FID = IID_final) %>% 
  rename(IID = IID_final) %>% 
  mutate(batch = paste("batch", new_batch_number, sep = "")) %>% #Create final batch variable
  mutate(cohort_PROBAND = ifelse(cohort == "PROBAND", 1, 0),
         cohort_Oxford = ifelse(cohort == "Oxford", 1, 0),
         cohort_PPMI = ifelse(cohort == "PPMI", 1, 0),
         array1 = ifelse(array == 1, 1, 0), #Create binary variables for each array
         array2 = ifelse(array == 2, 1, 0),
         array3 = ifelse(array == 3, 1, 0),
         array4 = ifelse(array == 4, 1, 0),
         array5 = ifelse(array == 5, 1, 0),
         array6 = ifelse(array == "WGS", 1, 0)) %>% 
  select(FID, IID, gender, cohort_PROBAND, cohort_Oxford, cohort_PPMI,
         array1, array2, array3, array4, array5, array6,
         allprogression_PC1, motorprogression_PC1, cogprogression_PC1) 

#Count how many individuals have progression score as well as imputed genetic data after QC
to_export %>% 
  filter(!is.na(allprogression_PC1)) %>% 
  summarise(count = n())

#Write as output
write.table(to_export, "outputs/PROBAND_OXFORD_PPMI_progression_perc_scaled_annualvisits_OFF_forGWAS_normfluency.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

#---Export IDs and cohort only---####

to_export_cohort <- PD_only %>% 
  select(ID, IID, V2, cohort, allprogression_PC1, gender, new_batch_number, array,
         motorprogression_PC1, cogprogression_PC1) %>% 
  filter(!is.na(IID)) %>% #Remove individuals missing IID
  filter(!is.na(new_batch_number)) %>%  #Remove individuals missing batch number
  filter(!is.na(V2) & cohort == "PROBAND" | 
           is.na(V2) & cohort == "Oxford" |
           is.na(V2) & cohort == "PPMI") %>%  #Remove only PROBAND individuals missing V2 (new IID)
  mutate(IID_final = ifelse(!is.na(V2), paste(V2),
                            ifelse(is.na(V2) & cohort == "Oxford", IID,
                                   ifelse(is.na(V2) & cohort == "PPMI", paste("PPMISI", IID, sep =""), NA)))) %>% 
  #Create final IID - for Oxford this is the IID, for PPMI it has a prefix of PPMISI
  select(IID_final, gender, cohort, new_batch_number, array, allprogression_PC1, motorprogression_PC1, cogprogression_PC1) %>% 
  mutate(FID = IID_final) %>% 
  rename(IID = IID_final) %>% 
  mutate(batch = paste("batch", new_batch_number, sep = "")) %>% #Create final batch variable
  mutate(cohort_PROBAND = ifelse(cohort == "PROBAND", 1, 0),
         cohort_Oxford = ifelse(cohort == "Oxford", 1, 0),
         cohort_PPMI = ifelse(cohort == "PPMI", 1, 0),
         array1 = ifelse(array == 1, 1, 0), #Create binary variables for each array
         array2 = ifelse(array == 2, 1, 0),
         array3 = ifelse(array == 3, 1, 0),
         array4 = ifelse(array == 4, 1, 0),
         array5 = ifelse(array == 5, 1, 0),
         array6 = ifelse(array == "WGS", 1, 0)) %>% 
  select(FID, IID, cohort) 

to_export_cohort %>% 
  group_by(cohort) %>% 
  summarise(count = n())

write.table(to_export_cohort, "outputs/IDs_cohorts.txt", quote = F, row.names = F, col.names = T, sep = "\t")

#---Export phenotype, cohort, batch number and gender data excluding extreme 5% of progressors---####

#Select top and bottom 5% of COMPOSITE progression scores (excluding people who are missing this score)
to_export_allprogression_top5pc <- to_export %>% 
  filter(!is.na(allprogression_PC1)) %>% 
  top_frac(0.05, allprogression_PC1)

to_export_allprogression_bot5pc <- to_export %>% 
  filter(!is.na(allprogression_PC1)) %>% 
  top_frac(-0.05, allprogression_PC1)

#Combine top and bottom 5% for composite progression
to_export_allprogression_extreme5pc <- rbind(to_export_allprogression_top5pc, to_export_allprogression_bot5pc)

to_export_allprogression_extreme5pc <- to_export_allprogression_extreme5pc %>% 
  select(IID) %>% 
  mutate(allprogression_extremes = "remove")


#Select top and bottom 5% of MOTOR progression scores (excluding people who are missing this score)
to_export_motorprogression_top5pc <- to_export %>% 
  filter(!is.na(motorprogression_PC1)) %>% 
  top_frac(0.05, motorprogression_PC1)

to_export_motorprogression_bot5pc <- to_export %>% 
  filter(!is.na(motorprogression_PC1)) %>% 
  top_frac(-0.05, motorprogression_PC1)

to_export_motorprogression_extreme5pc <- rbind(to_export_motorprogression_top5pc, to_export_motorprogression_bot5pc)

to_export_motorprogression_extreme5pc <- to_export_motorprogression_extreme5pc %>% 
  select(IID) %>% 
  mutate(motorprogression_extremes = "remove")

#Select top and bottom 5% of COGNITIVE progression scores (excluding people who are missing this score)
to_export_cogprogression_top5pc <- to_export %>% 
  filter(!is.na(cogprogression_PC1)) %>% 
  top_frac(0.05, cogprogression_PC1)

to_export_cogprogression_bot5pc <- to_export %>% 
  filter(!is.na(cogprogression_PC1)) %>% 
  top_frac(-0.05, cogprogression_PC1)

to_export_cogprogression_extreme5pc <- rbind(to_export_cogprogression_top5pc, to_export_cogprogression_bot5pc)

to_export_cogprogression_extreme5pc <- to_export_cogprogression_extreme5pc %>% 
  select(IID) %>% 
  mutate(cogprogression_extremes = "remove")



#Merge with main export dataset
#This is so we can keep just one single export file with the extreme progressors removed in each progression score
to_export_extreme5pc_removed <- to_export %>% 
  left_join(to_export_allprogression_extreme5pc, by = "IID") %>% 
  left_join(to_export_motorprogression_extreme5pc, by = "IID") %>% 
  left_join(to_export_cogprogression_extreme5pc, by = "IID") %>% 
  mutate(allprogression_PC1_extremesremoved = ifelse(is.na(allprogression_extremes), allprogression_PC1,
                                                     ifelse(allprogression_extremes == "remove", NA, "CHECK"))) %>% 
  mutate(motorprogression_PC1_extremesremoved = ifelse(is.na(motorprogression_extremes), motorprogression_PC1,
                                                       ifelse(motorprogression_extremes == "remove", NA, "CHECK"))) %>%
  mutate(cogprogression_PC1_extremesremoved = ifelse(is.na(cogprogression_extremes), cogprogression_PC1,
                                                     ifelse(cogprogression_extremes == "remove", NA, "CHECK"))) %>% 
  select(-allprogression_PC1, -motorprogression_PC1, -cogprogression_PC1, 
         -allprogression_extremes, -motorprogression_extremes, -cogprogression_extremes)

#Write as output
write.table(to_export_extreme5pc_removed, "outputs/PROBAND_OXFORD_PPMI_progression_perc_scaled_annualvisits_OFF_forGWAS_normfluency_extreme5pc.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")


########## EXPORT SINGLE SCALE PROGRESSION SCORES ########## 
#---Export random slopes for each scale---####
to_export_singlescales <- PD_only %>% 
  select(ID, IID, V2, cohort, allprogression_PC1, gender, new_batch_number, array,
         randomslope_UPDRSIIIxtime_from_onset,
         randomslope_UPDRSIIxtime_from_onset,
         randomslope_HYxtime_from_onset,
         randomslope_MOCAxtime_from_onset,
         randomslope_fluencyxtime_from_onset,
         randomslope_UPDRSI1xtime_from_onset) %>% 
  filter(!is.na(IID)) %>% #Remove individuals missing IID
  filter(!is.na(new_batch_number)) %>%  #Remove individuals missing batch number
  filter(!is.na(V2) & cohort == "PROBAND" | 
           is.na(V2) & cohort == "Oxford" |
           is.na(V2) & cohort == "PPMI") %>%  #Remove only PROBAND individuals missing V2 (new IID)
  mutate(IID_final = ifelse(!is.na(V2), paste(V2),
                            ifelse(is.na(V2) & cohort == "Oxford", IID,
                                   ifelse(is.na(V2) & cohort == "PPMI", paste("PPMISI", IID, sep =""), NA)))) %>% 
  #Create final IID - for Oxford this is the IID, for PPMI it has a prefix of PPMISI
  select(IID_final, gender, cohort, new_batch_number, array,
         randomslope_UPDRSIIIxtime_from_onset,
         randomslope_UPDRSIIxtime_from_onset,
         randomslope_HYxtime_from_onset,
         randomslope_MOCAxtime_from_onset,
         randomslope_fluencyxtime_from_onset,
         randomslope_UPDRSI1xtime_from_onset) %>% 
  mutate(FID = IID_final) %>% 
  rename(IID = IID_final) %>% 
  mutate(batch = paste("batch", new_batch_number, sep = "")) %>% #Create final batch variable
  mutate(cohort_PROBAND = ifelse(cohort == "PROBAND", 1, 0),
         cohort_Oxford = ifelse(cohort == "Oxford", 1, 0),
         cohort_PPMI = ifelse(cohort == "PPMI", 1, 0),
         array1 = ifelse(array == 1, 1, 0), #Create binary variables for each array
         array2 = ifelse(array == 2, 1, 0),
         array3 = ifelse(array == 3, 1, 0),
         array4 = ifelse(array == 4, 1, 0),
         array5 = ifelse(array == 5, 1, 0),
         array6 = ifelse(array == "WGS", 1, 0)) %>% 
  select(FID, IID, gender, cohort_PROBAND, cohort_Oxford, cohort_PPMI,
         array1, array2, array3, array4, array5, array6,
         randomslope_UPDRSIIIxtime_from_onset,
         randomslope_UPDRSIIxtime_from_onset,
         randomslope_HYxtime_from_onset,
         randomslope_MOCAxtime_from_onset,
         randomslope_fluencyxtime_from_onset,
         randomslope_UPDRSI1xtime_from_onset) 

#Write as output
write.table(to_export_singlescales, "outputs/PROBAND_OXFORD_PPMI_progression_perc_scaled_annualvisits_OFF_sepscales_normfluency_forGWAS.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")


