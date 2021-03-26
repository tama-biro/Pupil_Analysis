
library(ggplot2)

#### From Trial 360 data points onward ####

# Run through pupil signal and find point of gamble onset
# When found, pick out 60 values before (500ms) and 360 values after (3000ms)
# Subtract baseline (mean 500ms to 0ms before gamble onset) from all values
# Split by decision (accept/reject)

RejectLong <- data.frame(Numbers = 1, Pupil = 1)
AcceptLong <- data.frame(Numbers = 1, Pupil = 1)

for(i in which(pdataBW$StimulusName > 100 &
               pdataBW$trial_inclusion == 1 &
               pdataBW$subject_inclusion == 1)) {
  
  # Reject
  if(pdataBW$StimulusName[i] != pdataBW$StimulusName[i - 1] & 
     pdataBW$Decision[i] == 0) {
    
    x <- data.frame(Numbers = 1:420, 
                    Pupil = (pdataBW$Pupil_Z[(i - 60):(i + 359)] - 
                               mean(pdataBW$Pupil_Z[(i - 60):(i - 1)])))
    RejectLong <- rbind(RejectLong, x)
    print(i)
  }
  # Accept
  if(pdataBW$StimulusName[i] != pdataBW$StimulusName[i - 1] & 
     pdataBW$Decision[i] == 1) {
    
    x <- data.frame(Numbers = 1:420, 
                    Pupil = (pdataBW$Pupil_Z[(i - 60):(i + 359)] - 
                               mean(pdataBW$Pupil_Z[(i - 60):(i - 1)])))
    AcceptLong <- rbind(AcceptLong, x)
    print(i)
  }
}

# Remove first row of data frames
RejectLong <- RejectLong[2:nrow(RejectLong), ]
AcceptLong <- AcceptLong[2:nrow(AcceptLong), ]

# Getting standard error for plotting
# Data for each data point is isolated and SE for that point is calculated
# Output is a data frame with one value (SE) for each data point

RFOLSE <- data.frame(Numbers = 1, SE = 2)
AFOLSE <- data.frame(Numbers = 1, SE = 2)

for (i in 1:420) {
  a <- RejectLong[RejectLong$Numbers == i, ]
  RFOLSE <- rbind(RFOLSE, c(i, (sd(a$Pupil)/sqrt(length(a$Pupil)))))
  
  a <- AcceptLong[AcceptLong$Numbers == i, ]
  AFOLSE <- rbind(AFOLSE, c(i, (sd(a$Pupil)/sqrt(length(a$Pupil)))))
  
  print(i)
}

# Remove first row
RFOLSE <- RFOLSE[2:nrow(RFOLSE), ]
AFOLSE <- AFOLSE[2:nrow(AFOLSE), ]

# Make aggregated data frames
# Output is mean pupil size for each data point
RejectLongAgg <- aggregate(Pupil ~ Numbers, data = RejectLong, FUN = mean)
AcceptLongAgg <- aggregate(Pupil ~ Numbers, data = AcceptLong, FUN = mean)

# Make one big data frame for all data points and their SEs
LongPupilPlot <- cbind(RejectLongAgg, RFOLSE$SE,
                       AcceptLongAgg$Pupil, AFOLSE$SE)

# Rename columns
colnames(LongPupilPlot) <- c("Numbers", "Reject", "RejectSE",
                             "Accept", "AcceptSE")

# Smooth out plot
for(i in 7:414) {
  for(j in 2:5) {
    LongPupilPlot[i, j] <- mean(LongPupilPlot[(i - 6):(i + 6), j])
  }
}

# Rename columns
colnames(LongPupilPlot) <- c("Numbers", "Reject", "RejectSE",
                             "Accept", "AcceptSE")

# Plot all values with horizontal lines at mean RT for accept/reject
ggplot(LongPupilPlot, aes(x = Numbers)) +
  scale_x_continuous(breaks = c(0, 60, 120, 180, 240, 300, 360, 420), 
                     labels = c("-500", "0", "500", "1000", "1500",
                                "2000", "2500", "3000"),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(-.5, .5), expand = c(0, 0)) +
  labs(x = "Time from FixOff/ms", y = "Pupil Size (arbitrary unit)") +
  geom_vline(xintercept = 60) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_ribbon(aes(ymin = Reject - RejectSE, ymax = Reject + RejectSE),
              fill = '#DAE3F3') +
  geom_ribbon(aes(ymin = Accept - AcceptSE, ymax = Accept + AcceptSE),
              fill = '#FFCCCC') +
  geom_line(aes(y = Reject), col = '#4472C4', size = 1) +
  geom_line(aes(y = Accept), col = '#C00000', size = 1) +
  geom_vline(xintercept = 221, col = '#4472C4') +
  geom_vline(xintercept = 251, col = '#C00000')


ggsave("Long Split by Decision.png")


#### Decision Plot split by outcome ####

# Data frames for Decision plot
RejectDec <- data.frame(Numbers = 1, Pupil = 1)
AcceptDec <- data.frame(Numbers = 1, Pupil = 1)

# Runs through data frame and finds point of decision
# 120 data points before and 180 points (1/1.5s) after decision is collected
# Data also filtered by decision type (accept vs. reject)
# Trials that are not 1 second long are not used (removing approx. 1/4 of trial)


for(i in which(pdataBW$StimulusName > 100 &
               pdataBW$trial_inclusion == 1 &
               pdataBW$subject_inclusion == 1)) {
  
  # Reject
  if(pdataBW$StimulusName[i] == pdataBW$StimulusName[i - 119] & 
     pdataBW$StimulusName[i] != pdataBW$StimulusName[i + 1] & 
     pdataBW$Decision[i] == 0) {
    
    k <- which(pdataBW$StimulusName == pdataBW$StimulusName[i] &
                 pdataBW$Run == pdataBW$Run[i] &
                 pdataBW$Subject == pdataBW$Subject[i])
    
    x <- data.frame(Numbers = 1:300, 
                    Pupil = (pdataBW$Pupil_Z[(i - 119):(i + 180)] -
                               mean(pdataBW$Pupil_Z[(k[1] - 60):(k[length(k)] - 1)])))
    RejectDec <- rbind(RejectDec, x)
    print(i)
  }
  
  # Accept
  if(pdataBW$StimulusName[i] == pdataBW$StimulusName[i - 119] & 
     pdataBW$StimulusName[i] != pdataBW$StimulusName[i + 1] & 
     pdataBW$Decision[i] == 1) {
    
    k <- which(pdataBW$StimulusName == pdataBW$StimulusName[i] &
                 pdataBW$Run == pdataBW$Run[i] &
                 pdataBW$Subject == pdataBW$Subject[i])
    
    x <- data.frame(Numbers = 1:300, 
                    Pupil = (pdataBW$Pupil_Z[(i - 119):(i + 180)] -
                               mean(pdataBW$Pupil_Z[(k[1] - 60):(k[length(k)] - 1)])))
    AcceptDec <- rbind(AcceptDec, x)
    print(i)
  }
}

# Remove first row of data frames
RejectDec <- RejectDec[2:nrow(RejectDec), ]
AcceptDec <- AcceptDec[2:nrow(AcceptDec), ]

# Getting standard error for plotting
# Data for each data point is isolated and SE for that point is calculated
# Output is a data frame with one value (SE) for each data point

RDecSE <- data.frame(Numbers = 1, SE = 2)
ADecSE <- data.frame(Numbers = 1, SE = 2)

for(i in 1:300) {
  
  a <- RejectDec[RejectDec$Numbers == i, ]
  RDecSE <- rbind(RDecSE, c(i, (sd(a$Pupil)/sqrt(length(a$Pupil)))))
  
  a <- AcceptDec[AcceptDec$Numbers == i, ]
  ADecSE <- rbind(ADecSE, c(i, (sd(a$Pupil)/sqrt(length(a$Pupil)))))
  
  print(i)
}

# Remove first row
RDecSE <- RDecSE[2:nrow(RDecSE), ]
ADecSE <- ADecSE[2:nrow(ADecSE), ]

# Make aggregated data frames
# Output is mean pupil size for each data point
RejectDecAgg <- aggregate(Pupil ~ Numbers, data = RejectDec, FUN = mean)
AcceptDecAgg <- aggregate(Pupil ~ Numbers, data = AcceptDec, FUN = mean)

# Make one big data frame for all data points and their SEs
DecisionPlot <- cbind(RejectDecAgg, RDecSE$SE, 
                      AcceptDecAgg$Pupil,ADecSE$SE)

# Rename columns
colnames(DecisionPlot) <- c("Numbers", "Reject", "RejectSE",
                            "Accept", "AcceptSE")

# Smooth out plot
for(i in 7:294) {
  for(j in 2:5) {
    DecisionPlot[i, j] <- mean(DecisionPlot[(i - 6):(i + 6), j])
  }
}

# Plot all values
ggplot(DecisionPlot, aes(x = Numbers)) +
  scale_x_continuous(breaks = c(0, 60, 120, 180, 240, 300), 
                     labels = c('-1000', '-500', '0', '500', '1000', '1500'),
                     expand = c(0, 1)) +
  scale_y_continuous(limits = c(-.25, .25), expand = c(0, 0)) +
  labs(x = "Time from Decision/ms", y = "Pupil Size (arbitrary unit)") +
  geom_vline(xintercept = 120) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_ribbon(aes(ymin = Reject - RejectSE, ymax = Reject + RejectSE),
              fill = '#DAE3F3') +
  geom_ribbon(aes(ymin = Accept - AcceptSE, ymax = Accept + AcceptSE),
              fill = '#FFCCCC') +
  geom_line(aes(y = Reject), col = '#4472C4', size = 1) +
  geom_line(aes(y = Accept), col = '#C00000', size = 1)

ggsave('Decision Split by Decision.png')



