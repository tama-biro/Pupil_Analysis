##################################
#### Preprocessing Pupil Data ####
##################################

library(seewave)

# Load data (this can be found at https://osf.io/yz9e6/)
sample_level <- read.csv('decision.csv')
data_summary <- read.csv('eye.csv')

# StimulusName on data_summary
data_summary$StimulusName <- data_summary$Gain * 100 + data_summary$Loss

# Make new column with pupil left to work on
sample_level$PupilLeft_Work <- sample_level$PupilLeft

# Count number of missing items (-1) and extremely small (<1)
sum(sample_level$PupilLeft_Work == -1) # 425827 NAs
sum(sample_level$PupilLeft_Work > -1 & sample_level$PupilLeft_Work < 1) # 2674

# Make -1 (missing) and <1 (extremely small) values NA
sample_level$PupilLeft_Work[sample_level$PupilLeft_Work < 1] <- NA

# Filter through data and remove NAs
# Interpolate data from 100ms on either side
while(sum(is.na(sample_level$PupilLeft_Work)) > 1) {
  for(i in which(is.na(sample_level$PupilLeft_Work))){
    
    # Find first relevant data point for subject
    if(i <= 12) {
      first <- 1
    } else if(sample_level$Run[i] != sample_level$Run[i - 12]) {
      first <- which(sample_level$Subject == sample_level$Subject[i] &
                       sample_level$Run == sample_level$Run[i])[1]
    } else {
      first <- i - 12
    }
    # Find last relevant data point for subject
    if(i + 12 > nrow(sample_level)) {
      last <- nrow(sample_level)
    } else if(sample_level$Run[i] != sample_level$Run[i + 12]) {
      last <- which(sample_level$Subject == sample_level$Subject[i] &
                      sample_level$Run == sample_level$Run[i])[length(which(sample_level$Subject == sample_level$Subject[i] &
                                                                              sample_level$Run == sample_level$Run[i]))]
    } else {
      last <- i + 12
    }
    
    # Interpolate
    sample_level$PupilLeft_Work[i] <- mean(sample_level$PupilLeft_Work[first:last], na.rm = TRUE)
    
    print(i)
  }
}


#### Filtering, baseline correction and z-scoring ####

# Butterwoth filtering 0.02 Hz to 4 Hz, 3rd order
pdataBW <- cbind(sample_level[1,], data.frame(Pupil_Filtered = 1))

for(i in unique(sample_level$Subject*100 + sample_level$Run)) {
  a <- sample_level[(sample_level$Subject*100 + sample_level$Run) == i, ]
  TestBW <- as.data.frame(bwfilter(wave = a$PupilLeft_Work, f = 120,
                                   n = 3, from = 0.02, to = 4,
                                   bandpass = TRUE,
                                   output = "matrix"))
  
  pdataTemp <- cbind(a, TestBW$V1)
  colnames(pdataTemp) <- colnames(pdataBW)
  
  pdataBW <- rbind(pdataBW, pdataTemp)
}

pdataBW <- pdataBW[2:nrow(pdataBW),]

colnames(pdataBW)[length(pdataBW)] <- 'Pupil_Filtered'

MeanX <- aggregate(Pupil_Filtered ~ Subject + Run, data = pdataBW, FUN = mean, na.action = na.omit)
SDX <- aggregate(Pupil_Filtered ~ Subject + Run, data = pdataBW, FUN = sd, na.action = na.omit)
colnames(MeanX) <- c('Subject', 'Run', 'Mean')
colnames(SDX) <- c('Subject', 'Run', 'SD')

butterS <- data.frame(Subject = pdataBW$Subject, Run = pdataBW$Run)

for(j in unique(butterS$Subject*100 + butterS$Run)) {
  butterS$Mean[(butterS$Subject*100 + butterS$Run) == j] <- MeanX$Mean[(MeanX$Subject*100 + MeanX$Run) == j]
  butterS$SD[(butterS$Subject*100 + butterS$Run) == j] <- SDX$SD[(SDX$Subject*100 + SDX$Run) == j]
}

pdataBW$Pupil_Z <- (pdataBW$Pupil_Filtered - butterS$Mean)/butterS$SD

# Move decision, trial inclusion, and subject inclusion to pdataBW
for(i in 1:nrow(data_summary)) {
  
  pdataBW$trial_inclusion[pdataBW$Subject == data_summary$Subject[i] & 
                            pdataBW$Run == data_summary$Block[i] &
                            pdataBW$StimulusName == data_summary$StimulusName[i]] <-
    data_summary$TrialInclusion[i]
  
  pdataBW$subject_inclusion[pdataBW$Subject == data_summary$Subject[i] & 
                              pdataBW$Run == data_summary$Block[i] &
                              pdataBW$StimulusName == data_summary$StimulusName[i]] <-
    data_summary$SubjectInclusion[i]
  
  pdataBW$Decision[pdataBW$Subject == data_summary$Subject[i] & 
                     pdataBW$Run == data_summary$Block[i] &
                     pdataBW$StimulusName == data_summary$StimulusName[i]] <-
    data_summary$Decision[i]
  
  print(i)
}

write.csv(pdataBW, 'pdataBW_processed.csv', row.names = FALSE)


# Make NAs in all data_summary columns we're about to fill
data_summary[,59:89] <- NA

# Find data for each trial and move to data_summary and baseline correct
for(i in which(pdataBW$StimulusName > 100)) {
  if(pdataBW$StimulusName[i + 1] != pdataBW$StimulusName[i]) {
    
    # Save all rows of pdataBW that are in the trial
    k <- which(pdataBW$StimulusName == pdataBW$StimulusName[i] &
                 pdataBW$Run == pdataBW$Run[i] &
                 pdataBW$Subject == pdataBW$Subject[i])
    
    # Baseline pupil
    data_summary[data_summary$Subject == pdataBW$Subject[i] &
                  data_summary$Block == pdataBW$Run[i] &
                  data_summary$StimulusName == pdataBW$StimulusName[i], 59] <-
                      mean(pdataBW$Pupil_Z[(k[1] - 60):(k[1] - 1)])
    
    # Take the RT for the gamble and divide by 120
    # This will give the first bin with data we pick out
    firstBin <- round((data_summary$RT[data_summary$Subject == pdataBW$Subject[i] &
                                  data_summary$StimulusName == pdataBW$StimulusName[i] &
                                  data_summary$Block == pdataBW$Run[i]])/100, 0)
    
    # We only want data 1500ms before bin, so make sure firstBin is 15 or less
    if(firstBin > 15) {
      firstBin <- 15
    }
    
    # Enter baseline corrected data into every column
    for(j in -firstBin:14) {
      data_summary[data_summary$Subject == pdataBW$Subject[i] &
                     data_summary$Block == pdataBW$Run[i] &
                     data_summary$StimulusName == pdataBW$StimulusName[i], j + 75] <-
        mean(pdataBW$Pupil_Z[(i + 12*j + 1):(i + 12*(j + 1))]) - 
        mean(pdataBW$Pupil_Z[(k[1] - 60):(k[1] - 1)])
    }
    print(i)
  }
}


colnames(data_summary)[59:89] <- c('Baseline', 
                                   'Neg15', 'Neg14', 'Neg13', 'Neg12', 'Neg11',
                                   'Neg10', 'Neg9', 'Neg8', 'Neg7', 'Neg6',
                                   'Neg5', 'Neg4', 'Neg3', 'Neg2', 'Neg1',
                                   'Pos1', 'Pos2', 'Pos3', 'Pos4', 'Pos5',
                                   'Pos6', 'Pos7', 'Pos8', 'Pos9', 'Pos10',
                                   'Pos11', 'Pos12', 'Pos13', 'Pos14', 'Pos15')



write.csv(data_summary[,c(1:57, 59:ncol(data_summary))], 
          'data_summary_pupil.csv', row.names = FALSE)


