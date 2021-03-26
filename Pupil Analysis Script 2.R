
library(lme4)
library(ggplot2)

# Remove missing trials
data_analyze <- data_summary[data_summary$TrialInclusion == 1 &
                               data_summary$SubjectInclusion == 1, ]

#### T-tests ####
# Running t-tests for difference P(A) - P(R)
t_values <- data.frame(Bin = c(seq(-1450, -50, 100), seq(50, 1450, 100)))

for(j in 60:89) {
    
  t <- t.test(data_analyze[data_analyze$Decision == 1, j],
              data_analyze[data_analyze$Decision == 0, j])
    
  t_values$t_statistic[j - 59] <- t$statistic
  t_values$p_value[j - 59] <- t$p.value
    
}

ggplot(t_values, aes(x = Bin, y = t_statistic)) +
  geom_bar(stat = 'identity', col = 'black', fill = 'white') +
  geom_hline(yintercept = 1.96, col = 'red') +
  labs(x = 'Bin (100ms)', y = 't-statistic') +
  theme_minimal()

ggsave('t_values_pupil_bar.png')

write.csv(t_values, 't_values_pupil.csv')

# T-test & Cohen's D for whole time window at decision
data_analyze$Pupil_Dec <- apply(data_analyze[, 71:89], 1, mean, na.rm = TRUE)

t.test(data_analyze$Pupil_Dec[data_analyze$Decision == 1],
       data_analyze$Pupil_Dec[data_analyze$Decision == 0])

(mean(data_analyze$Pupil_Dec[data_analyze$Decision == 1], na.rm = T) - 
    mean(data_analyze$Pupil_Dec[data_analyze$Decision == 0], na.rm = T))/
  sqrt(
    (sd(data_analyze$Pupil_Dec[data_analyze$Decision == 1], na.rm = T)^2 +
       sd(data_analyze$Pupil_Dec[data_analyze$Decision == 0], na.rm = T)^2)/
      2
  )

# Correlation pupil difference & accceptance rate
PupilDiff <- data.frame(Subject = unique(data_analyze$Subject),
                        AccRate = rep(0, 62),
                        PupilDiff = rep(0, 62))

for(i in unique(data_analyze$Subject)) {
  
  PupilDiff$AccRate[PupilDiff$Subject == i] <- mean(data_analyze$Decision[data_analyze$Subject == i], na.rm = T)
  PupilDiff$PupilDiff[PupilDiff$Subject == i] <-
    mean(data_analyze$Pupil_Dec[data_analyze$Subject == i & data_analyze$Decision == 1], na.rm = T) -
    mean(data_analyze$Pupil_Dec[data_analyze$Subject == i & data_analyze$Decision == 0], na.rm = T)
  
}

cor.test(PupilDiff$PupilDiff, PupilDiff$AccRate)

for(i in 1:62) {
  data_analyze$AccRate[data_analyze$Subject == i] <- PupilDiff$AccRate[PupilDiff$Subject == i]
}

#### Mixed Effects Modeling ####

b_values <- data.frame(Bin = c(seq(-1450, -50, 100), seq(50, 1450, 100)))

for(i in 60:89) {
  
  mixed <- glmer(Decision ~ data_analyze[, i] + Gain + Loss +
                   (1 | Subject), data = data_analyze,
                 family = binomial(link = 'logit'),
                 control = glmerControl(optCtrl = list(maxfun = 2e5)))
  
  beta_phasic <- as.data.frame(summary(mixed)[10])[2,1]
  p <- as.data.frame(summary(mixed)[10])[2,4]
  
  b_values$beta[i - 59] <- beta_phasic
  b_values$p_value[i - 59] <- p
  
}


# LMM models
mixed1 <- lmer(Pupil_Dec ~ Gain + Loss + Decision +
                 (1 + Gain + Loss + Decision | Subject), 
               data = data_analyze, REML = F,
               control = lmerControl(optCtrl = list(maxfun = 2e5)))

mixed2 <- lmer(Pupil_Dec ~ Gain + Loss + Decision + AccRate +
                 (1 + Gain + Loss + Decision | Subject), 
               data = data_analyze, REML = F,
               control = lmerControl(optCtrl = list(maxfun = 2e5)))

mixed3 <- lmer(Pupil_Dec ~ Gain + Loss + Decision*AccRate +
                 (1 + Gain + Loss + Decision | Subject), 
               data = data_analyze, REML = F,
               control = lmerControl(optCtrl = list(maxfun = 2e5)))


AIC(mixed1, mixed2, mixed3)
anova(mixed1, mixed2, mixed3)

summary(mixed3)

