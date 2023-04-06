###################################################################################
#                             Genetic summary statistic
#         Population genetic summary statistics and its relationship to patch size
#
#                               Adriel M. Sierra
#
#
#                                   2022
#
###################################################################################


#########genetic diversity relationship
#library
library(ggpubr)
library(Rmisc)
library(ggplot2)
library(car)
library(lme4)
library(Matrix)
library(hrbrthemes)
library(RColorBrewer)
library(scales)
library(pegas)
library(poppr)
library(hierfstat)
library(cowplot)

# Read genind objects for both species
#Radula flaccida
radula.pop.105 = readRDS(file = "Radula_105R15_genind.rds")
radula.pop.70 = readRDS(file = "Radula_70R20_genind.rds")
radula.popmap.105 = read.table("Radula_105_popmap_com.txt", header=T)
radula.popmap.70 = read.table("Radula_71_popmap_com.txt", header=T)

#Cololejeunea surinamensis
cololej.pop.107 = readRDS(file = "Cololejeunea_107R15_genind.rds")
cololej.pop.80 = readRDS(file = "Cololejeunea_80R20_genind.rds")
cololej.popmap.108 = read.table("Cololej_108_popmap_com.txt", header=T)
cololej.popmap.71 = read.table("Cololej_71_popmap_com.txt", header=T)

# basic genetic statistics
basic_radula.105 = basic.stats(radula.pop.105, diploid = FALSE)
basic_radula.105
basic_radula.70 = basic.stats(radula.pop.70, diploid = FALSE)
basic_radula.70

basic_cololej.107 = basic.stats(cololej.pop.107, diploid = FALSE)
basic_cololej.107
basic_cololej.71 = basic.stats(cololej.pop.71, diploid = FALSE)
basic_cololej.71

# Mean  gene diversity per site
Hs_Radula.105 = apply(basic_radula.105$Hs, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 2)
Hs_Radula.105

Hs_Radula.70 = apply(basic_radula.70$Hs, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 2)
Hs_Radula.70

Hs_cololej.107 = apply(basic_cololej.107$Hs, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 2)
Hs_cololej.107

Hs_cololej.71 = apply(basic_cololej.71$Hs, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 2)
Hs_cololej.71

#Print the number of multilocus genotypes
mlg(radula.pop.105)
mlg(radula.pop.70)
mlg(cololej.pop.107)
mlg(cololej.pop.71)

##Genotype accumulation curve
gac <- genotype_curve(cololej.pop, sample = 1000, quiet = TRUE)
library("ggplot2")
p <- last_plot()
p + geom_smooth()

p + geom_smooth()+ 
  labs(x="Number of loci sampled", y ="Multilocus genotype (MLG)")+ 
  theme_bw() + 
  theme(text = element_text(size = 12)) +
  scale_x_continuous(breaks=seq(0, 3100, 200), limits=c(0, 3100))+
  theme_classic()

##Calculating genotypic diversity
poppr(radula.pop.105)
poppr(radula.pop.70)
poppr(cololej.pop.107)
poppr(cololej.pop.71)

#private alleles by pop
private_alleles(radula.pop.105) %>% apply(MARGIN = 1, FUN = sum)
private_alleles(radula.pop.70) %>% apply(MARGIN = 1, FUN = sum)
private_alleles(cololej.pop.107) %>% apply(MARGIN = 1, FUN = sum)
private_alleles(cololej.pop.71) %>% apply(MARGIN = 1, FUN = sum)

#convert data.genind to a "hierfstat" object
radula.hf.105 <- genind2hierfstat(radula.pop.105)
dim(radula.hf.105)

radula.hf.70 <- genind2hierfstat(radula.pop.70)
dim(radula.hf.70)

cololej.hf.107 <- genind2hierfstat(cololej.pop.107)
dim(cololej.hf.107)

cololej.hf.71 <- genind2hierfstat(cololej.pop.71)
dim(cololej.hf.71)

##allelic richness

allelic.richness(radula.hf.105, min.n = 107, diploid=FALSE)$Ar %>%
  apply(MARGIN = 2, FUN = sum, na.rm=TRUE)

allelic.richness(radula.hf.70, min.n = 107, diploid=FALSE)$Ar %>%
  apply(MARGIN = 2, FUN = sum, na.rm=TRUE)

allelic.richness(cololej.pop.107, min.n = 107, diploid=FALSE)$Ar %>%
  apply(MARGIN = 2, FUN = sum, na.rm=TRUE)

allelic.richness(cololej.hf.71, min.n = 107, diploid=FALSE)$Ar %>%
  apply(MARGIN = 2, FUN = sum, na.rm=TRUE)


# Genetic diversity parameters were compiled in a csv
gen<-read.csv("Gen div Stats.csv", header=T, sep=",")
str(gen)
gen

# factors variables
gen$Species = as.factor(gen$Species)
gen$Population.ID = as.factor(gen$Population.ID)
gen$Size = as.factor(gen$Size)
gen$Size.class = as.factor(gen$Size.class)
gen$Population.size.2000 = as.integer(gen$Population.size.2000)
gen$change = as.integer(gen$change)

gen$Size.class <- factor(gen$Size.class, levels=c("Small", "Large"))

# Plot loop
pd <- position_dodge(0.1) 

plot_for_loop <- function(df, Size.class, y_var) {
  
  ggplot(df, aes(x = Size, y = .data[[y_var]], color=Species, linetype= Species, group=Species)) + 
    geom_point(position = pd, size=3)+ labs(x="Forest fragment size", y =y_var) + 
    #geom_line() +
    geom_smooth(method=lm , color="Black", fill="lightgrey", se=TRUE)+
    # colouring
    scale_color_brewer(palette="Dark2")+
    theme_classic(base_size = 12)
}

plot_list1 <- colnames(gen)[-c(1:9,20:29)] %>% 
  map( ~ plot_for_loop(gen, colnames(gen)[1], .x))

plot_list2 <- colnames(gen)[-c(1:19, 22:23)] %>% 
  map( ~ plot_for_loop(gen, colnames(gen)[1], .x))

# view all plots individually (not shown)
plot_list1
plot_list2

# Combine all plots

Figure_S6 = plot_grid(plotlist = plot_list1,
          ncol = 2)

Figure_S7 = plot_grid(plotlist = plot_list2,
          ncol = 2)

# statistical regressions
###Radula flaccida
summary(lm(Variant.sites.15~Size.class,data=gen[c(1:11),]))
summary(lm(Variant.sites.20~Size.class,data=gen[c(1:11),]))

summary(lm(Polymorphic.sites.15~Size.class,data=gen[c(1:11),]))
summary(lm(Polymorphic.sites.20~Size.class,data=gen[c(1:11),]))

summary(lm(Allelic.richness.15~Size.class,data=gen[c(1:11),]))
summary(lm(Allelic.richness.20~Size.class,data=gen[c(1:11),]))

summary(lm(Private.alleles.15~Size.class,data=gen[c(1:11),]))
summary(lm(Private.alleles.20~Size.class,data=gen[c(1:11),]))

summary(lm(Pi.πT.15~Size.class,data=gen[c(1:11),]))
summary(lm(Pi.πT.20~Size.class,data=gen[c(1:11),]))

summary(lm(MLG.15~Size.class,data=gen[c(1:11),]))
summary(lm(MLG.20~Size.class,data=gen[c(1:11),]))

summary(lm(H.15~Size.class,data=gen[c(1:11),]))
summary(lm(H.20~Size.class,data=gen[c(1:11),]))

summary(lm(G.15~Size.class,data=gen[c(1:11),]))
summary(lm(G.20~Size.class,data=gen[c(1:11),]))

summary(lm(lambda.15~Size.class,data=gen[c(1:11),]))
summary(lm(lambda.20~Size.class,data=gen[c(1:11),]))

###Cololejeunea surinamensis

summary(lm(Variant.sites.15~Size.class,data=gen[c(12:22),]))
summary(lm(Variant.sites.20~Size.class,data=gen[c(12:22),]))

summary(lm(Polymorphic.sites.15~Size.class,data=gen[c(12:22),]))
summary(lm(Polymorphic.sites.20~Size,data=gen[c(12:22),]))

summary(lm(Allelic.richness.15~Size.class,data=gen[c(12:22),]))
summary(lm(Allelic.richness.20~Size.class,data=gen[c(12:22),]))

summary(lm(Private.alleles.15~Size.class,data=gen[c(12:22),]))
summary(lm(Private.alleles.20~Size.class,data=gen[c(12:22),]))

summary(lm(Pi.πT.15~Size.class,data=gen[c(12:22),]))
summary(lm(Pi.πT.20~Size.class,data=gen[c(12:22),]))

summary(lm(MLG.15~Size.class,data=gen[c(12:22),]))
summary(lm(MLG.20~Size.class,data=gen[c(12:22),]))

summary(lm(H.15~Size.class,data=gen[c(12:22),]))
summary(lm(H.20~Size.class,data=gen[c(12:22),]))

summary(lm(G.15~Size.class,data=gen[c(12:22),]))
summary(lm(G.20~Size.class,data=gen[c(12:22),]))

summary(lm(lambda.15~Size.class,data=gen[c(12:22),]))
summary(lm(lambda.20~Size.class,data=gen[c(12:22),]))

