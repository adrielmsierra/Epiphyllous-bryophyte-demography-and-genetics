###################################################################################
#               Estimates of epiphyllous population size dynamics
#         Change of epiphyllous population occupancy patterns in fifteen years
#
#                               Adriel M. Sierra
#
#
#                                   2021
#
###################################################################################

rm(list=ls())
####Directory
setwd("/Users/adrielsierra/Documents/MS publication/Epiphyllous_fragment_pop_gen/manus_analisis")

#Library
library(ggpubr)
library(Rmisc)
library(ggplot2)
library(car)
library(lme4)
library(Matrix)
library(hrbrthemes)
library(RColorBrewer)
library(scales)
library(raincloudplots)
library(dplyr)

#Data
demo<-read.csv("demographic dataII.csv", sep=";")
str(demo)

#Variables as factors
demo$Year<-as.factor(demo$Year)
demo$Reserve<-as.factor(demo$Reserve)
demo$Fragsize<-as.factor(demo$Fragsize)
demo$frag<-as.factor(demo$frag)

###### Cloud plots


##Radula flaccida
rf_small_1x1 <- data_1x1(
  array_1 = demo$Rflaccida.popsize[1:7],
  array_2 = demo$Rflaccida.popsize[21:27],
  jit_distance = .09,
  jit_seed = 321)
head (rf_small_1x1)
str(rf_small_1x1)
rf_small_1x1$year=as.factor(c(rep("2000", 7), rep("2016", 7)))
wilcox.test(y_axis ~ year, data=rf_small_1x1) 

raincloud_rf_small <- raincloud_1x1_repmes(
  data = rf_small_1x1,
  colors = (c('#d95f02', '#d95f02')),
  fills = (c('#d95f02', '#d95f02')),
  line_color = 'gray',
  line_alpha = .3,
  size = 1,
  alpha = .6,
  align_clouds = FALSE) +
  scale_x_continuous(breaks=c(1,2), labels=c("2000", "2016"), limits=c(0, 3)) +
  ylim(0,1000)+
  xlab("Year") + 
  ylab("Estimated number of colonies per 1-ha") +
  theme_classic()+
  ggtitle('W = 2, p-value = 0.002')+
  theme(plot.title = element_text(hjust = 0.5))

raincloud_rf_small

rf_large_1x1 <- data_1x1(
  array_1 = demo$Rflaccida.popsize[8:19],
  array_2 = demo$Rflaccida.popsize[30:41],
  jit_distance = .09,
  jit_seed = 321)
head (rf_large_1x1)
str(rf_large_1x1)
rf_large_1x1$year=as.factor(c(rep("2000", 12), rep("2016", 12)))
wilcox.test(y_axis ~ year, data=rf_large_1x1) 

raincloud_rf_large <- raincloud_1x1_repmes(
  data = rf_large_1x1,
  colors = (c('#1b9e77', '#1b9e77')),
  fills = (c('#1b9e77', '#1b9e77')),
  line_color = 'gray',
  line_alpha = .3,
  size = 1,
  alpha = .6,
  align_clouds = FALSE) +
  
  scale_x_continuous(breaks=c(1,2), labels=c("2000", "2016"), limits=c(0, 3)) +
  xlab("Year") + 
  ylab("Estimated number of colonies per 1-ha") +
  theme_classic()+
  ggtitle('W = 47, p-value = 0.16')+
  theme(plot.title = element_text(hjust = 0.5))

raincloud_rf_large

### Cololejeunea surinamensis

cs_small_1x1 <- data_1x1(
  array_1 = demo$Csurinamensis.popsize[1:7],
  array_2 = demo$Csurinamensis.popsize[21:27],
  jit_distance = .09,
  jit_seed = 321)
head (cs_small_1x1)
str(cs_small_1x1)
cs_small_1x1$year=as.factor(c(rep("2000", 7), rep("2016", 7)))
wilcox.test(y_axis ~ year, data=cs_small_1x1) 

raincloud_cs_small <- raincloud_1x1_repmes(
  data = cs_small_1x1,
  colors = (c('#d95f02', '#d95f02')),
  fills = (c('#d95f02', '#d95f02')),
  line_color = 'gray',
  line_alpha = .3,
  size = 1.5,
  alpha = .6,
  align_clouds = FALSE) +
  ylim(0, 850)+
  scale_x_continuous(breaks=c(1,2), labels=c("2000", "2016"), limits=c(0, 3)) +
  xlab("Year") + 
  ylab("Estimated number of colonies per 1-ha") +
  theme_classic()+
  ggtitle('W = 0, p-value = 0.0005')+
  theme(plot.title = element_text(hjust = 0.5))

raincloud_cs_small

cs_large_1x1 <- data_1x1(
  array_1 = demo$Csurinamensis.popsize[8:19],
  array_2 = demo$Csurinamensis.popsize[30:41],
  jit_distance = .09,
  jit_seed = 321)
head (cs_large_1x1)
str(cs_large_1x1)
cs_large_1x1$year=as.factor(c(rep("2000", 12), rep("2016", 12)))
wilcox.test(y_axis ~ year, data=cs_large_1x1) 

raincloud_cs_large <- raincloud_1x1_repmes(
  data = cs_large_1x1,
  colors = (c('#1b9e77', '#1b9e77')),
  fills = (c('#1b9e77', '#1b9e77')),
  line_color = 'gray',
  line_alpha = .3,
  size = 1.5,
  alpha = .6,
  align_clouds = FALSE) +
  ylim(0, 850)+
  scale_x_continuous(breaks=c(1,2), labels=c("2000", "2016"), limits=c(0, 3)) +
  xlab("Year") + 
  ylab("Estimated number of colonies per 1-ha") +
  theme_classic()+
  ggtitle('W = 33, p-value = 0.024')+
  theme(plot.title = element_text(hjust = 0.5))


raincloud_cs_large


#Radula flaccida 

### Group mean
# by size class in the year 2000 and 2016
mura <- ddply(demo, .(Fragsize, Year), summarise, grp.mean=mean(Rflaccida.popsize))

#Histogram plot
histra<-ggplot(demo, aes(x=Rflaccida.popsize, color=Fragsize, fill=Fragsize)) +
  geom_histogram(aes(y=..density..), bins = 20, position="identity", alpha=0.5)+
  geom_density(alpha=0.6)+
  geom_vline(data=mura, aes(xintercept=grp.mean, color=Fragsize),
             linetype=c("solid", "dashed","solid", "dashed"))+
  scale_color_brewer(palette="Dark2", breaks=c("Small", "Large"))+
  scale_fill_brewer(palette="Dark2", breaks=c("Small", "Large"))+
  labs(x="Population size", y = "Density")+
  theme_classic()
histra 

# specific fragment size for the year 2000
mura20 <- ddply(demo[c(1:20),], .(frag, Year), summarise, grp.mean=mean(Rflaccida.popsize))
mura20$Fragsize = c("Small","Small","Large","Large")

#Histogram plot
histra20<-ggplot(demo[c(1:20),], aes(x=Rflaccida.popsize, color=Fragsize, fill=Fragsize)) +
  geom_histogram(aes(y=..density..), bins = 14, position="identity", alpha=0.5)+
  geom_density(alpha=0.6)+
  geom_vline(data=mura20, aes(xintercept=grp.mean, color=Fragsize),
             linetype=c("solid", "dashed","solid", "dashed"))+
  scale_color_brewer(palette="Dark2", breaks=c("Small", "Large"))+
  scale_fill_brewer(palette="Dark2", breaks=c("Small", "Large"))+
  labs(x=" ", y = "Density")+
  theme_classic()
histra20 

# non-parametric Wilcoxon rank sum exact test
# 1 ha and 10 ha vs 100 ha and Continuous forest
wilcox.test(Rflaccida.popsize ~ Fragsize, data=demo[c(1:20),])
# 10 ha vs 100 ha and Continuous forest
wilcox.test(Rflaccida.popsize ~ Fragsize, data=demo[c(4:20),])

# specific fragment size for the year 2016
mura16 <- ddply(demo[c(21:44),], .(frag, Year), summarise, grp.mean=mean(Rflaccida.popsize))
mura16$Fragsize = c("Small","Small","Large","Large")

#Histogram plot
histra16<-ggplot(demo[c(21:44),], aes(x=Rflaccida.popsize, color=Fragsize, fill=Fragsize)) +
  geom_histogram(aes(y=..density..), bins = 14, position="identity", alpha=0.5)+
  geom_density(alpha=0.6)+
  geom_vline(data=mura16, aes(xintercept=grp.mean, color=Fragsize),
             linetype=c("solid", "dashed","solid", "dashed"))+
  scale_color_brewer(palette="Dark2", breaks=c("Small", "Large"))+
  scale_fill_brewer(palette="Dark2", breaks=c("Small", "Large"))+
  labs(x="Estimated number of colonies per 1-ha", y = "Density")+
  theme_classic()
histra16 
                       
# non-parametric Wilcoxon rank sum exact test
# 1 ha and 10 ha vs 100 ha and Continuous forest
wilcox.test(Rflaccida.popsize ~ Fragsize, data=demo[c(21:44),])
# 10 ha vs 100 ha and Continuous forest
wilcox.test(Rflaccida.popsize ~ Fragsize, data=demo[c(24:28,30:44),])


################################################################################
#Cololejeunea surinamensis

### Group mean
# by size class in the year 2000 and 2016
mucs <- ddply(demo, .(Fragsize, Year), summarise, grp.mean=mean(Csurinamensis.popsize))

# specific fragment size for the year 2000
mucs20 <- ddply(demo[c(1:20),], .(frag, Year), summarise, grp.mean=mean(Csurinamensis.popsize))
mucs20$Fragsize = c("Small","Small","Large","Large")
#Histogram plot
histcs20<-ggplot(demo[c(1:20),], aes(x=Csurinamensis.popsize, color=Fragsize, fill=Fragsize))+
  geom_histogram(aes(y=..density..), bins = 14, position="identity", alpha=0.5)+
  geom_density(alpha=0.6)+
  geom_vline(data=mucs20, aes(xintercept=grp.mean, color=Fragsize),
             linetype=c("solid", "dashed","solid", "dashed"))+
  scale_color_brewer(palette="Dark2", breaks=c("Small", "Large"))+
  scale_fill_brewer(palette="Dark2", breaks=c("Small", "Large"))+
  labs(x=" ", y = "Density")+
  theme_classic()
histcs20

#### non-parametric Wilcoxon rank sum exact test
# 1 ha and 10 ha vs 100 ha and Continuous forest
wilcox.test(Csurinamensis.popsize ~ Fragsize, data=demo[c(1:20),])
# 10 ha vs 100 ha and Continuous forest
wilcox.test(Csurinamensis.popsize ~ Fragsize, data=demo[c(4:20),])


# specific fragment size for the year 2016
mucs16 <- ddply(demo[c(21:44),], .(frag, Year), summarise, grp.mean=mean(Csurinamensis.popsize))
mucs16$Fragsize = c("Small","Small","Large","Large")

#Histogram plot
histcs16<-ggplot(demo[c(21:44),], aes(x=Csurinamensis.popsize, color=Fragsize, fill=Fragsize))+
  geom_histogram(aes(y=..density..), bins = 14, position="identity", alpha=0.5)+
  geom_density(alpha=0.6)+
  geom_vline(data=mucs16, aes(xintercept=grp.mean, color=Fragsize),
             linetype=c("solid", "dashed","solid", "dashed"))+
  scale_color_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  labs(x="Estimated number of colonies per 1-ha", y = "Density")+
  theme_classic()
histcs16


#### non-parametric Wilcoxon rank sum exact test
# 1 ha and 10 ha vs 100 ha and Continuous forest
wilcox.test(Csurinamensis.popsize ~ Fragsize, data=demo[c(21:44),])
# 10 ha vs 100 ha and Continuous forest
wilcox.test(Csurinamensis.popsize ~ Fragsize, data=demo[c(24:28,30:44),])


#Final figures
# 
# cloud plots
figure.cloud <- ggarrange(raincloud_rf_small, raincloud_rf_large, raincloud_cs_small, raincloud_cs_large,
                          labels = c("Radula flaccida", " " , "Cololejeunea surinamensis", " "), hjust= 0, vjust = 0.85,
                          font.label = list(size = 9, color = "black", family = NULL),
                          ncol = 2, nrow = 2)
figure.cloud
# histograms
figure.histog <- ggarrange(histra20, histra16, histcs20, histcs16,
                          labels = c("2000", "2016" , "2000", "2016"), hjust= 0, vjust = 0.85,
                          font.label = list(size = 9, color = "black", family = NULL),
                          ncol = 1, nrow = 4)
figure.histog

# Compile Final figure
figure1 <- ggarrange(figure.cloud, figure.histog,
                     ncol = 2, nrow = 1)
figure1
