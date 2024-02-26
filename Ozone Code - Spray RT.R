########################
#### Ozone Project #####
########################

# Sahaana laptop
setwd("~/All School Stuff/All Graduate School/Research/Non-Thesis Paper Writing/Adam Papers/Ozone")

# import required packages
require(ggplot2) # for plotting
library(tidyverse)
library(gridExtra)
require(gridExtra)
library(dplyr)
require(car) # to get p-values for glm
require(emmeans); require(multcomp) # for post-hoc analysis
require(multcompView)
library(multcomp)
require(gdata) # contains keep()
require(ggpubr) # contains ggarrange()
library(rstatix)
require(lme4)
require(lmerTest) # for mixed-models and GLMMs
library(kableExtra)

# import data
ozone <- read.csv(file="Ozone.csv", sep=",", header=T)

# format data
ozone <- within(ozone, {Ozone..ppm. <- as.factor(ozone$Ozone..ppm.);
Pathogen <- as.factor(Pathogen);
Produce <- as.factor(Produce);
Treatment <- as.factor(Treatment);
Time..s. <- as.factor(Time..s.);
Prep <- as.factor(Prep);
LogReductionTenth <- as.numeric(LogReductionTenth);
Log.Reduction <- as.numeric(Log.Reduction)})

###################################################################
###################################################################
############# Spray Treatment - Room Temperature ##################
###################################################################
###################################################################

# keep only spray, room temperature data
ozonesprayRT <-ozone[ozone$Prep == "Room Temp", ]

# keep only Lm and Salmonella data
ozonesprayRTLS <-ozonesprayRT[ozonesprayRT$Pathogen %in% c("L. monocytogenes", "Salmonella"),]

# visualize data
ggplot(data=ozonesprayRTLS, aes(x=Produce, color=Ozone..ppm., y=Log.Reduction))+
  geom_boxplot()+
  facet_wrap(~Pathogen)

# raw data plot
SprayRT_raw <- ggplot(ozonesprayRTLS, aes(x = Produce, y = Log.Reduction, shape =Ozone..ppm.))+ geom_point(size = 5) + facet_grid(~Pathogen) + 
  labs(shape = "Ozone (ppm):") +
  xlab("Produce") +
  ylab("Log Reduction") +
  ylim(0, 6.5) +
  scale_shape_manual(values = c(21,22)) + 
  theme(plot.title = element_text(hjust = 0.5, face = "plain"),
        panel.background = element_rect(fill="white", color="black"),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(face = "plain", size = 20),
        axis.title = element_text(size=20),
        axis.text = element_text(size=16, color="black"),
        panel.grid.minor.y = element_line(color="grey45", linetype="solid", size=0.25),
        legend.title = element_text(size=20),
        legend.text = element_text(size = 19),
        axis.text.y = element_text(color= "black"),
        legend.key=element_blank(),
        axis.text.x = element_text(color= "black")) 

ggsave(filename = "SprayRT_Raw.jpg", width=12, height=8, unit="in", plot = SprayRT_raw)

# --------------- Linear Model, Assuming Normality and Homogeneity of Variances ------------
# optimum linear model

SprayRT_LM <- lm(Log.Reduction ~ Produce*Ozone..ppm.*Pathogen , data = ozonesprayRTLS)

# display anova table in the console
anova(SprayRT_LM) 
summary(SprayRT_LM)

plot(SprayRT_LM, which = 1) 
plot(SprayRT_LM, which = 2) ## Q-Q Plot 

## Assumption of linear model is not validated. Need to consider GLM-Poisson

----------------------------------------------------------------------------
## GLM - Poisson
  
SprayRT_Poisson <- glm(LogReductionTenth ~ Ozone..ppm. * Produce * Pathogen, poisson(link = "log"), data=ozonesprayRTLS)
summary(SprayRT_Poisson)
anova(SprayRT_Poisson)
str(ozonesprayRTLS)

# Residual residual deviance is higher than the residual degrees of freedom 
# GLM- Quasipoisson

SprayRT_Poisson <- glm(LogReductionTenth ~ Ozone..ppm. * Produce * Pathogen
                       -Ozone..ppm.:Produce:Pathogen
                       -Ozone..ppm.:Produce
                       -Ozone..ppm.:Pathogen, quasipoisson(link = "log"), data=ozonesprayRTLS)
summary(SprayRT_Poisson)
anova(SprayRT_Poisson)

# no interactions hence removed everything above

SprayRT_Poisson.res <- emmeans(SprayRT_Poisson, ~ Ozone..ppm.+Produce+Pathogen, level=0.95)
SprayRT_Poisson.res <- cld(SprayRT_Poisson.res, alpha = 0.05, Letters = letters, decreasing = TRUE)
SprayRT_Poisson.res <- data.frame(SprayRT_Poisson.res)
SprayRT_Poisson.res

# Convert log values of estimate, asump.LCL, and asymp.UCL. Divide by 100
SprayRT_Poisson.res$estimate <- exp(SprayRT_Poisson.res$emmean)/100
SprayRT_Poisson.res$LCL <- exp(SprayRT_Poisson.res$asymp.LCL)/100
SprayRT_Poisson.res$UCL <- exp(SprayRT_Poisson.res$asymp.UCL)/100

# remove excess spaces in .group
SprayRT_Poisson.res$.group <- gsub(" ", "", SprayRT_Poisson.res$.group)
SprayRT_Poisson.res

# plot results
SprayRT_Poisson.res <- ggplot(data = SprayRT_Poisson.res, aes(x=Produce, y=estimate, fill=Ozone..ppm.))+
  facet_wrap(~Pathogen)+
  geom_bar(stat="identity", position = position_dodge(width=0.7), color="black", width=0.7)+
  geom_errorbar(aes(ymin = LCL, ymax = UCL), position = position_dodge(width=0.7), width=0.2)+
  geom_text(aes(y = UCL+0.2, label=.group), position = position_dodge(width=0.7), size=5)+
  xlab("Produce") + ylab("Log Reduction")+
  #ggtitle("Spray Treatment-Room Temperature Prep")+
  theme(plot.title = element_text(size = 20, face = "bold", hjust=0.5))+
  scale_fill_manual(values=c("grey", "black"), name = "Ozone (ppm)") +
  theme(axis.text.x = element_text(size = 18)) +   
  theme(axis.text.y = element_text(size = 18)) +
  theme(axis.text.y = element_text(colour = "black")) +
  theme(axis.text.x = element_text(colour = "black")) +
  theme(axis.title = element_text(size = 20)) +
  theme(strip.text = element_text(face = "plain", size = 18)) +
  theme(strip.background = element_rect(fill=NA))+
  theme(legend.title = element_text(size = 20)) +
  theme(legend.text = element_text(size = 20)) +
  theme(panel.background = element_rect(fill=NA, color="grey15"),
        panel.grid.major.x =  element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color="grey15",linetype="dashed", size=0.1),
        panel.grid.minor.y = element_line(color="grey45",linetype="dotted", size=0.05))
SprayRT_Poisson.res

ggsave(filename = "Ozone Spray-RT.jpg", width=12, height=6, unit="in", plot = SprayRT_Poisson.res, dpi=100)