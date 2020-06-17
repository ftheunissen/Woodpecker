getwd()
setwd("../8_Playbacks/")

library(ade4)
library(lme4)
library(ggplot2)
library(car)
library(RVAideMemoire)
library(lsmeans)
library(reshape)
library(optimx)


rm(list=ls())#empty the workspace
dat = read.csv("PB_Analysis_Manip1.csv", header=T)

numdat = dat[10:15]
numdat = scale(numdat, center =T, scale = T)
numdat = as.data.frame(numdat)

dat.pca<-dudi.pca(numdat,center = F, scale = F, scannf = T) #from examining PC values, n = 2 components to retain
summary(dat.pca)
write.csv(dat.pca$c1[1:2], 'Playback_PCA_loadings.csv')

pc1<-dat.pca$li[,1]
pc2 <- dat.pca$li[,2]


dat$Group=as.factor(dat$Group)
dat$StudyDay = scale(dat$StudyDay, center = T, scale = T)
dat$StudyTime = scale(dat$StudyTime, center = T, scale = T)


StudyDay <- dat$StudyDay
StudyTime <- dat$StudyTime
Group <- dat$Group
Treatment <- dat$Treatment
SpeciesType <- dat$TypeBroadcast
SpeciesID <- dat$SpeciesBroadcast
Subject <- dat$Test.subject
Order <- dat$Order


mFull1 <- lmer(pc1~Treatment*Group + (1|StudyDay) + (1|StudyTime) + (1|Order) + (1|Subject))
plotresid(mFull1)
shapiro.test(residuals(mFull1))
Anova(mFull1) #package (car)
lsmFull1 <- lsmeans(mFull1, ~Treatment|Group, adjust="tukey")
lsmFull1
contrast(lsmFull1)


mFull2 <- lmer(pc2~Treatment*Group + (1|StudyDay) + (1|StudyTime) + (1|Order) + (1|Subject))

plotresid(mFull2)
shapiro.test(residuals(mFull2))
Anova(mFull2) #package (car)
lsmFull2 <- lsmeans(mFull2, ~Treatment|Group, adjust="tukey")
lsmFull2
contrast(lsmFull2)



#preparing data for plotting: 
SpeciesID <- factor(SpeciesID, levels = c('Dmajor','Dhyperythrus', 'Dsyriacus', 'Dminor', 'Pcanus'),ordered = TRUE)
#Groupordered <- factor(Group,c('1', '2', '3', '4'), ordered = T)
Groupordered <- factor(Group,c('2','1', '4','3'), ordered = T) # ordered following phylogenetic distance with D.Major

# possibility to write a PDF for the plot:
#pdf ("Manip_INTERSPECIES_Response_to_PB.pdf")
ggplot(dat, aes(factor(Groupordered), pc1)) + 
  geom_boxplot(aes(fill = SpeciesID), outlier.shape = NA) +
  theme_bw() +
  ylab("Behavioral response (Dim.1)") +
  theme(axis.line = element_line(colour = "black", size=0.7),
        axis.title.x=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.ticks = element_line(size=0.7),
        axis.text = element_text(size=20),
        axis.text.x=element_text(angle=45, size=15, vjust=0.5, face ="bold.italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
#  scale_x_discrete(labels=c("Group1","Group2","Group3","Group4")) +
  scale_x_discrete(labels=c("n = 6","n = 6","n = 6","n = 6")) +
  scale_fill_manual(values=c("Grey100", "Cyan1", "Cyan3", "grey40", "grey75"))
#dev.off()

# possibility to write a PDF for the plot:
#pdf ("Manip_INTERSPECIES_Response_to_PB_PC2.pdf")
ggplot(dat, aes(factor(Groupordered), pc2)) + 
  geom_boxplot(aes(fill = SpeciesID), outlier.shape = NA) +
  theme_bw() +
  ylab("Behavioral response (Dim.2)") +
  theme(axis.line = element_line(colour = "black", size=0.7),
        axis.title.x=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.ticks = element_line(size=0.7),
        axis.text = element_text(size=20),
        axis.text.x=element_text(angle=45, size=15, vjust=0.5, face ="bold.italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  #  scale_x_discrete(labels=c("Group1","Group2","Group3","Group4")) +
  scale_x_discrete(labels=c("n = 6","n = 6","n = 6","n = 6")) +
  scale_fill_manual(values=c("Grey100", "Cyan1", "Cyan3", "grey40", "grey75"))
#dev.off()





# now on to analyzing the other experimental set using resynthesis and for which PC values are computed from the scores
# obtained in the first experiment

# read in the csv file
dat2 = read.csv("PB_Analysis_Manip2.csv", header=T)

# center-scale the behavioral data before applying the linear equation of the two significant PCs obtained from the PCA in the first experiment
numdat2 = dat2[9:14]
numdat2 = scale(numdat2, center =T, scale = T)
numdat2 = as.data.frame(numdat2)

#prepare to store the newly computed PC values
dat2$PC1 <- 0
dat2$PC2 <- 0

# load PC scores obtained from PCA in first experiment
loadingMatrix <- t(dat.pca$c1[1:2])

for (experiment in 1:length(dat2[,1])) {
  dat2$PC1[experiment] <- sum(numdat2[experiment, ]*loadingMatrix[1,])
  dat2$PC2[experiment] <- sum(numdat2[experiment, ]*loadingMatrix[2,])
}

PC1_exp2 <- dat2$PC1
PC2_exp2 <- dat2$PC2

# similar to experiment 1, prepare data for stats and plots:
dat2$Group=as.factor(dat2$Group)
dat2$StudyDay = scale(dat2$StudyDay, center = T, scale = T)
dat2$StudyTime = scale(dat2$StudyTime, center = T, scale = T)

StudyDay_exp2 <- dat2$StudyDay
StudyTime_exp2 <- dat2$StudyTime
Group_exp2 <- dat2$Group
Treatment_exp2 <- dat2$Treatment
SpeciesType_exp2 <- dat2$TypeBroadcast
Subject_exp2 <- dat2$Test.subject
Order_exp2 <- dat2$Order

mFull1_exp2 <- lmer(PC1_exp2~Treatment_exp2*Group_exp2 + (1|StudyDay_exp2) + (1|StudyTime_exp2) + (1|Order_exp2) + (1|Subject_exp2))

plotresid(mFull1_exp2)
shapiro.test(residuals(mFull1_exp2))
Anova(mFull1_exp2) #package (car)
lsmFull1_exp2 <- lsmeans(mFull1_exp2, ~Treatment_exp2|Group_exp2, adjust="tukey")
lsmFull1_exp2
contrast(lsmFull1_exp2)


mFull2_exp2 <- lmer(PC2_exp2~Treatment_exp2*Group_exp2 + (1|StudyDay_exp2) + (1|StudyTime_exp2) + (1|Order_exp2) + (1|Subject_exp2))

plotresid(mFull2_exp2)
shapiro.test(residuals(mFull2_exp2))
Anova(mFull2_exp2) #package (car)
lsmFull2_exp2 <- lsmeans(mFull2_exp2, ~Treatment_exp2|Group_exp2, adjust="tukey")
lsmFull2_exp2
contrast(lsmFull2_exp2)

#preparing data for plotting: 
SpeciesType_exp2 <- factor(SpeciesType_exp2, levels = c('Major','No_Accel', 'No_Amp', 'No_Spec', 'Flat'),ordered = TRUE)
Groupordered_exp2 <- factor(Group_exp2,c('1','3','4','2'), ordered = T)


# possibility to write a PDF for the plotting:
#pdf ("Manip_Resynthesis_Response_to_PB.pdf")
ggplot(dat2, aes(factor(Groupordered_exp2), PC1_exp2)) + 
  geom_boxplot(aes(fill = SpeciesType_exp2), outlier.shape = NA) +
  theme_bw() +
  ylab("Behavioral response (Dim.1)") +
  theme(axis.line = element_line(colour = "black", size=0.7),
        axis.title.x=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.ticks = element_line(size=0.7),
        axis.text = element_text(size=20),
        axis.text.x=element_text(angle=45, size=15, vjust=0.5, face ="bold.italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  #  scale_x_discrete(labels=c("Group1","Group2","Group3","Group4")) +
  scale_x_discrete(labels=c("n = 6","n = 6","n = 6","n = 6")) +
  scale_fill_manual(values=c("Grey100", "burlywood1", "sandybrown", "lightsalmon3", "darkorange4")) 
#dev.off()

# possibility to write a PDF for the plotting:
#pdf ("Manip_Resynthesis_Response_to_PB_PC2.pdf")
ggplot(dat2, aes(factor(Groupordered_exp2), PC2_exp2)) + 
  geom_boxplot(aes(fill = SpeciesType_exp2), outlier.shape = NA) +
  theme_bw() +
  ylab("Behavioral response (Dim.2)") +
  theme(axis.line = element_line(colour = "black", size=0.7),
        axis.title.x=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.ticks = element_line(size=0.7),
        axis.text = element_text(size=20),
        axis.text.x=element_text(angle=45, size=15, vjust=0.5, face ="bold.italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  #  scale_x_discrete(labels=c("Group1","Group2","Group3","Group4")) +
  scale_x_discrete(labels=c("n = 6","n = 6","n = 6","n = 6")) +
  scale_fill_manual(values=c("Grey100", "burlywood1", "sandybrown", "lightsalmon3", "darkorange4")) 
#dev.off()
