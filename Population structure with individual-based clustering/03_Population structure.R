###################################################################################
#     Population structure: individual-based clustering and pairwise differentiation
#
#                               Adriel M. Sierra
#
#
#                                   2022
#
###################################################################################


#library
library(ggpubr)
library(Rmisc)
library(ggplot2)
library(car)
library(Matrix)
library(hrbrthemes)
library(RColorBrewer)
library(scales)
library(pegas)
library(poppr)
library(hierfstat)
library(cowplot)
library(ape)
library(adegenet)

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

#--------------#
#
# Perform DAPC
#
#--------------#

# Perform cross validation to find the optimal number of PCs to retain in DAPC
set.seed(123)
# Replace missing data with zero
radx <- tab(radula.pop, NA.method = "zero")

str(radx)
radcrossval <- xvalDapc(radx, radula.pop$pop, result = c("groupMean","overall"), xval.plot = TRUE)

# Number of PCs with best stats (lower score = better)
radcrossval$`Root Mean Squared Error by Number of PCs of PCA`
radcrossval$`Number of PCs Achieving Highest Mean Success`
radcrossval$`Number of PCs Achieving Lowest MSE`
numPCs <- as.numeric(radcrossval$`Number of PCs Achieving Lowest MSE`)

# Run a DAPC using population IDs as priors
rad.dapc1 = dapc(radula.pop, radula.pop$pop, n.pca = numPCs, n.da = 5)
rad.dapc1

# Analyse how much percent of genetic variance is explained by each axis
radpercent <- rad.dapc1$eig/sum(rad.dapc1$eig)*100
radpercent
is.numeric(radpercent)
rad.percent.plot<-barplot(radpercent, ylab = "Percent of genetic variance explained by eigenvectors", 
                          names.arg = round(radpercent, 2))

# Create a dataframe containing individual coordinates
rad.ind_coords <- as.data.frame(rad.dapc1$ind.coord)

#--------------#
#
# Visualise results
#
#--------------#

# Rename columns of dataframe
colnames(rad.ind_coords) = c("Axis1","Axis2","Axis3")
# Add a column containing individuals
rad.ind_coords$Ind = indNames(radula.pop)

# Add a column with the population IDs
rad.ind_coords$Pop = radula.pop$pop


# Adding labels of class size to the dataframe 
addlabel = function(x){
  # If population label x is present function will output y
  if(x=="col1"|x=="pa1"|x=="dim1") y = "1 ha"
  if(x=="col10"|x=="pa10"|x=="dim10") y = "10 ha"
  if(x=="dim100") y = "100 ha"
  if(x=="Cabofrio"|x=="flor"|x=="km41"|x=="dimcont") y = "Continuous"
  return(y)
}
rad.ind_coords$FragmentSize = sapply(rad.ind_coords$Pop, addlabel)  
head(rad.ind_coords)

# Calculate centroid (average) position for each population
centroid = aggregate(cbind(Axis1, Axis2, Axis3) ~ Pop,
                     data = rad.ind_coords,
                     FUN = mean)
head(centroid)
#for size class data

# Add centroid coordinates to ind_coords dataframe
rad.ind_coords<-rad.ind_coords[-c(4:5)]
#rad.ind_coords = left_join.(rad.ind_coords, centroid, by = "Pop", suffix = c("",".cen"))
head(rad.ind_coords)

rad.ind_coords = full_join.(rad.ind_coords, centroid, by = "Pop", suffix = c("",".cen"))
head(rad.ind_coords)
rad.ind_coords$SizeClass = radula.popmap[,4]

# Add region labels to centroid dataframe
centroid = centroid[-5]
centroid$FragmentSize = sapply(centroid$Pop, addlabel)  
centroid

# Define colour palette
show_col(brewer.pal(11, "Paired"))
cols = brewer.pal(nPop(radula.pop), "Paired") ; cols

##few colors
show_col(brewer.pal(8, "Dark2"))
cols = brewer.pal(nPop(radula.pop), "Dark2") ; cols

# Custom x and y labels
xlab = paste("Axis 1 (", format(round(radpercent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(radpercent[2], 1), nsmall=1)," %)", sep="")

# Custom ggplot2 theme
ggtheme = theme(legend.title = element_blank(),
                axis.text.y = element_text(colour="black", size=14),
                axis.text.x = element_text(colour="black", size=14),
                axis.title = element_text(colour="black", size=14),
                legend.position = "right",
                legend.text = element_text(size=15),
                legend.key = element_rect(fill = NA),
                legend.key.size = unit(0.7, "cm"),
                legend.box.spacing = unit(0, "cm"),
                panel.border = element_rect(colour="black", fill=NA, size=1),
                panel.background = element_blank(),
                # title centered
                plot.title = element_text(hjust=0.5, size=25) 
)

# Scatter plot axis 1 vs. 2 [colour by population]
ggplot(data = rad.ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spider segments
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Pop), show.legend = FALSE)+
  # points
  geom_point(aes(fill = Pop), shape = 21, size = 2, show.legend = TRUE)+
  # centroids
  #geom_label(data = centroid, aes(label = Pop, fill = Pop), size = 5, show.legend = FALSE)+
  # stat ellipse
  stat_ellipse(aes(fill = Pop, colour = Pop, show.legend = TRUE), level=0.95, size=0.75, show.legend = FALSE)+
  # colouring
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  # custom labels
  labs(x = xlab, y = ylab)+
  # custom theme
  ggtheme

#Export plot

#################
#
#compute fst
#
#################

radula_fst = genet.dist(radula.pop, diploid=FALSE, method = "WC84")
#radula_fst = genet.dist(radula.pop, diploid=FALSE, method = "Fst")
#radula_fst = genet.dist(radula.pop2, diploid=FALSE, method = "Nei87")


radula_fst[radula_fst < 0] = 0
radula_fst = radula_fst %>% round(digits = 4)

radula_fst_boot = boot.ppfst(radula.pop)
str(radula_fst_boot)

radula_fst_boot$vc.per.loc
radula_fst_boot$ll[is.na(radula_fst_boot$ll)] = 00
radula_fst_boot$ul[is.na(radula_fst_boot$ul)] = 00



radula_fst %>% round(digits = 5)
radula_fst_boot %>% round(digits = 4)

#Visualize results

# Convert dist object to data.frame
fst.matrix = as.matrix(radula_fst)
ind = which( upper.tri(fst.matrix), arr.ind = TRUE)
fst.df = data.frame(Site1 = dimnames(fst.matrix)[[2]][ind[,2]],
                    Site2 = dimnames(fst.matrix)[[1]][ind[,1]],
                    Fst = fst.matrix[ ind ] %>% round(digits = 3))
fst.df$class = rep("Fst", 55)

fst.matrix_boot = as.matrix(radula_fst_boot$ul)
ind = which( upper.tri(fst.matrix_boot), arr.ind = TRUE)
fst.df_boot = data.frame(Site1 = dimnames(fst.matrix_boot)[[2]][ind[,2]],
                         Site2 = dimnames(fst.matrix_boot)[[1]][ind[,1]],
                         Fst = fst.matrix_boot[ ind ] %>% round(digits = 3))

fst.df_boot$class = rep("upperlim", 55)
fst.df_boot_up = fst.df_boot
fst.df_boot_low = fst.df_boot

fst.df.com = data.frame(rbind(fst.df,fst.df_boot_low, fst.df_boot_up))

# Convert minus values to zero
fst.df$Fst[fst.df$Fst < 0] = 0

# Print data.frame summary
fst.df %>% str

# Fst italic label
fst.label = expression(italic("F")[ST])

# Extract middle Fst value for gradient argument
mid = max(fst.df$Fst) / 2

# Plot heatmap
ggplot(data = fst.df, aes(x = Site1, y = Site2, fill = Fst))+
  geom_tile(colour = "black")+
  geom_text(aes(label = Fst), color="black", size = 4)+
  scale_fill_gradient2(low = "blue", mid = "pink", high = "red", midpoint = mid, name = fst.label, limits = c(0, max(fst.df$Fst)), breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), position = "right")+
  theme(axis.text = element_text(colour = "black", size = 12),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12)
  )


###
str(fst.df.com)
fst.df.com$Sitescomp=as.factor(paste(fst.df.com$Site1,fst.df.com$Site2))

head(fst.df.com)
addcomparison = function(x){
  if(x=="col10 col1"|x=="dim1 col1"|x=="dim1 col10"|x=="dim10 col1"|x=="dim10 col10"|x=="dim10 dim1"|x=="pa1 col1"|x=="pa1 col10"|x=="pa1 dim1"|x=="pa1 dim10"|x=="pa10 col1"|x=="pa10 col10"|x=="pa10 dim1"|x=="pa10 dim10"|x=="pa10 pa1") y = "Small-Small"
  if(x=="col1 Cabofrio"|x=="col10 Cabofrio"|x=="dim1 Cabofrio"|x=="dim10 Cabofrio"|x=="dim100 col1"|x=="dim100 col10"|x=="dim100 dim1"|x=="dim100 dim10"|x=="dimcont col1"|x=="dimcont col10"|x=="dimcont dim1"|x=="dimcont dim10"|x=="flor col1"|x=="flor col10"|x=="flor dim1"|x=="flor dim10"|x=="km41 col1"|x=="km41 col10"|x=="km41 dim1"|x=="km41 dim10"|x=="pa1 Cabofrio"|x=="pa1 dim100"|x=="pa1 dimcont"|x=="pa1 flor"|x=="pa1 km41"|x=="pa10 Cabofrio"|x=="pa10 dim100"|x=="pa10 dimcont"|x=="pa10 flor"|x=="pa10 km41") y = "Small-Large"
  if(x=="dim100 Cabofrio"|x=="dimcont Cabofrio"|x=="dimcont dim100"|x=="flor Cabofrio"| x=="flor dim100"|x=="flor dimcont"|x=="km41 Cabofrio"|x=="km41 dim100"|x=="km41 dimcont"|x=="km41 flor") y = "Large-Large"
  return(y)
}

str(fst.df.com)
fst.df.com$comp = as.factor(sapply(c(fst.df.com$Sitescomp), addcomparison))  

levels(fst.df.com$comp)
levels(fst.df.com$comp) = c("Small-Small", "Small-Large","Large-Large")
head(fst.df.com)


str(fst.df.com)
fst.df.com$class<-as.factor(fst.df.com$class)
uplow = rep("Upper_lower_limits", 110)
manovel = rep("Mean value", 55)
fst.df.com$values = as.factor(c(manovel, uplow))

#######

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
# sample size
library(dplyr)

fst.df.com.boot = fst.df.com[c(56:165),]
fst.df.com.mean = fst.df.com[c(1:55),]
sample_size = fst.df.com %>% group_by(comp) %>% summarize(num=n())

#install.packages("see")
library(see)
fst.df.com.mean = fst.df.com.mean %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(comp, "\n", "n=", num))

fst.df.com.mean$myaxis = as.factor(c(rep("Small-Small\nn=10",10), rep("Small-Large\nn=30",30), rep("Large-Large\nn=15",15)))
levels(fst.df.com.mean$myaxis) = c("Small-Small\nn=10", "Small-Large\nn=30", "Large-Large\nn=15")

write.table(fst.df.com, "cololej_80R20_fstdf_oct22.csv", sep = ",")

plotfst.boot<- ggplot(fst.df.com.boot, aes(x=myaxis, y=Fst)) +
  #geom_violin(trim=FALSE, width=0.4, fill="#D95F02")+
  geom_violindot(trim=FALSE, width=0.6, fill="#D95F02", size_dots=0.20) +
  #geom_boxplot(width=0.1)+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  labs(x="Forest fragment size comparison", y =fst.label)+
  stat_summary(fun.data = data_summary, size = 1.5)+
  theme_classic()
plotfst.boot

plotfst<-ggplot(fst.df.com.mean, aes(x=myaxis, y=Fst)) +
  #geom_violin(trim=FALSE, width=0.6, fill="#1B9E77")+
  geom_violindot(trim=FALSE, width=0.6, fill="#1B9E77", size_dots=0.20) +
  #geom_boxplot(width=0.4, fill="darkgrey")+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  labs(x=" ", y = fst.label)+
  scale_color_brewer(palette="Dark2")+
  stat_summary(fun.data = data_summary, size = 1.5)+
  theme_classic()
plotfst

figure.fst <- ggarrange(plotfst, plotfst.boot,
                        labels = c("A", "B"), hjust= 0, vjust = 0.85,
                        font.label = list(size = 9, color = "black", family = NULL),
                        ncol = 1, nrow = 2)
##AMOVA
radula.pop
strata(radula.pop) <- strata(radula.pop)[,4:2]
radula.pop.bis <- as.genclone(radula.pop)
radula.pop.bis

amova.result <- poppr.amova(radula.pop.bis, ~sizeClass/size/subpop, missing = 'zero')

amova.test <- randtest(amova.result, nrepet = 999) # Test for significance
plot(amova.test)
############
#
#
#Cololejeunea surinamensis pop strut
#
#
############
# Perform cross validation to find the optimal number of PCs to retain in DAPC
set.seed(124)
# Replace missing data with zero
colx <- tab(cololej.pop, NA.method = "zero")
#--------------#
#
# Perform DAPC
#
#--------------#

# Perform cross validation to find the optimal number of PCs to retain in DAPC
set.seed(123)
colx <- tab(cololej.pop, NA.method = "mean")
crossvalcol <- xvalDapc(colx, cololej.pop$pop, result = "groupMean", xval.plot = TRUE)

# Number of PCs with best stats (lower score = better)
crossvalcol$`Root Mean Squared Error by Number of PCs of PCA`
crossvalcol$`Number of PCs Achieving Highest Mean Success`
crossvalcol$`Number of PCs Achieving Lowest MSE`
numPCs = as.numeric(crossvalcol$`Number of PCs Achieving Lowest MSE`)

# Run a DAPC using population IDs as priors
col.dapc1 = dapc(cololej.pop, cololej.pop$pop, n.pca = numPCs, n.da = 5)
col.dapc1

# Analyse how much percent of genetic variance is explained by each axis
col.percent = col.dapc1$eig/sum(col.dapc1$eig)*100
col.percent
barplot(col.percent, ylab = "Percent of genetic variance explained by eigenvectors", 
        names.arg = round(col.percent, 2))

# Create a dataframe containing individual coordinates
col.ind_coords = as.data.frame(col.dapc1$ind.coord)


#--------------#
#
# Perform PCA
#
#--------------#

# Replace missing data with the mean allele frequencies
#col.x = tab(cololej.pop, NA.method = "zero")
#col.x

# Perform PCA
#col.pca1 = dudi.pca(col.x, scannf = FALSE, scale = FALSE, nf = 3)
#col.pca1

# Analyse how much percent of genetic variance is explained by each axis
#col.percent = col.pca1$eig/sum(col.pca1$eig)*100
#col.percent
#barplot(col.percent, ylab = "Percent of genetic variance explained by eigenvectors",
#        names.arg = round(col.percent, 2))

# Create a dataframe containing individual coordinates
#col.ind_coords = as.data.frame(col.pca1$li)


#--------------#
#
# Visualise results
#
#--------------#

# Rename columns of dataframe
colnames(col.ind_coords) = c("Axis1","Axis2","Axis3")

# Add a column containing individuals
col.ind_coords$Ind = indNames(cololej.pop)

# Add a column with the population IDs
col.ind_coords$Pop = cololej.pop$pop
head(col.ind_coords)

# Conditional function that adds other labels to dataframe (optional)
addlabel = function(x){
  # If population label x is present function will output y
  if(x=="col1"|x=="pa1"|x=="dim1") y = "1 ha"
  if(x=="col10"|x=="pa10"|x=="dim10") y = "10 ha"
  if(x=="dim100") y = "100 ha"
  if(x=="Cabofrio"|x=="flor"|x=="km41"|x=="dimcont") y = "Continuous"
  return(y)
}
col.ind_coords$FragmentSize = sapply(cololej.pop$pop, addlabel)  
head(col.ind_coords)




# Reorder labels for plotting (optional)
# ind_coords$Region = factor(ind_coords$Region, levels = c("Region2","Region1"))

# Calculate centroid (average) position for each population
col.centroid = aggregate(cbind(Axis1, Axis2, Axis3) ~ FragmentSize,
                         data = col.ind_coords,
                         FUN = mean)

# Add centroid coordinates to ind_coords dataframe
#col.ind_coords = left_join.(col.ind_coords, col.centroid, by = "Pop")
head(col.ind_coords)
#head(col.ind_coords[-(4:7)])
#col.ind_coords<-col.ind_coords[-(4:7)]
col.ind_coords<-col.ind_coords[-c(4:5)]
head(col.ind_coords)
col.ind_coords = full_join.(col.ind_coords, col.centroid, by = "FragmentSize", suffix = c("",".cen"))
head(col.ind_coords)
col.centroid$FragmentSize = sapply(col.centroid$Pop, addlabel)  

# Add region labels to centroid dataframe
#col.centroid$Region = sapply(col.centroid$Pop, addlabel)  

col.centroid

# Define colour palette
show_col(brewer.pal(8, "Dark2"))
cols = brewer.pal(nPop(cololej.pop), "Dark2") ; cols

# Custom x and y labels
xlab = paste("Axis 1 (", format(round(col.percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(col.percent[2], 1), nsmall=1)," %)", sep="")

# Custom ggplot2 theme
ggtheme = theme(legend.title = element_blank(),
                axis.text.y = element_text(colour="black", size=14),
                axis.text.x = element_text(colour="black", size=14),
                axis.title = element_text(colour="black", size=14),
                legend.position = "bottom",
                legend.text = element_text(size=15),
                legend.key = element_rect(fill = NA),
                legend.key.size = unit(0.7, "cm"),
                legend.box.spacing = unit(0, "cm"),
                panel.border = element_rect(colour="black", fill=NA, size=1),
                panel.background = element_blank(),
                # title centered
                plot.title = element_text(hjust=0.5, size=25) 
)

# Scatter plot axis 1 vs. 2 [colour by population]
ggplot(data = col.ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spider segments
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = FragmentSize), show.legend = FALSE)+
  # points
  geom_point(aes(fill = FragmentSize), shape = 21, size = 2, show.legend = TRUE)+
  # centroids
  #geom_label(data = centroid, aes(label = Pop, fill = Pop), size = 5, show.legend = FALSE)+
  # stat ellipse
  stat_ellipse(aes(fill = FragmentSize, colour = FragmentSize), level=0.95, size=0.7, show.legend = FALSE)+
  # colouring
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  # custom labels
  labs(x = xlab, y = ylab)+
  # custom theme
  ggtheme

# Scatter plot axis 1 vs. 2 [colour by region]
p1 = ggplot(data = col.ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spider segments
  
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = FragmentSize), show.legend = FALSE)+
  # points
  geom_point(aes(fill = FragmentSize, colour = FragmentSize), shape = 21, size = 2, show.legend = TRUE)+
  # centroids
  #geom_label(data = centroid, aes(label = Pop, fill = FragmentSize), size = 5, show.legend = FALSE)+
  # stat ellipse
  stat_ellipse(aes(fill = FragmentSize, colour = FragmentSize, show.legend = TRUE),  level=0.95, size=0.7, show.legend = FALSE)+
  # colouring
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  # custom labels
  labs(x = xlab, y = ylab)+
  # custom theme
  ggtheme

# Add density curves to y and x axis
library(cowplot)
xdens <- 
  axis_canvas(p1, axis = "x") + 
  geom_density(data = col.ind_coords, aes(x = Axis1, fill = FragmentSize, colour = FragmentSize), alpha = 0.5)+
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)
ydens <-
  axis_canvas(p1, axis = "y", coord_flip = TRUE) + 
  geom_density(data = col.ind_coords, aes(x = Axis2, fill = FragmentSize, colour = FragmentSize), alpha = 0.5) +
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  coord_flip()
p1 %>%
  insert_xaxis_grob(xdens, grid::unit(1.2, "in"), position = "top") %>%
  insert_yaxis_grob(ydens, grid::unit(1.2, "in"), position = "right") %>%
  ggdraw()
# Export plot

###########
#
#Compute fst
#
###########
cololej_fst = genet.dist(cololej.pop, diploid=FALSE, method = "WC84")
#cololej_fst = genet.dist(cololej.pop, diploid=FALSE, method = "Fst")
#cololej_fst = genet.dist(cololej.pop, diploid=FALSE, method = "Nei87")

cololej_fst %>% round(digits = 4)
cololej_fst[cololej_fst < 0] = 0
cololej_fst = cololej_fst %>% round(digits = 4)
#Visualize results
# Convert dist object to data.frame
fst.matrix = as.matrix(cololej_fst)
ind = which( upper.tri(fst.matrix), arr.ind = TRUE)
fst.df = data.frame(Site1 = dimnames(fst.matrix)[[2]][ind[,2]],
                    Site2 = dimnames(fst.matrix)[[1]][ind[,1]],
                    Fst = fst.matrix[ ind ] %>% round(digits = 3))

# Convert minus values to zero
fst.df$Fst[fst.df$Fst < 0] = 0

# Print data.frame summary
fst.df %>% str

# Fst italic label
fst.label = expression(italic("F")[ST])

# Extract middle Fst value for gradient argument
mid = max(fst.df$Fst) / 2

# Plot heatmap
ggplot(data = fst.df, aes(x = Site1, y = Site2, fill = Fst))+
  geom_tile(colour = "black")+
  geom_text(aes(label = Fst), color="black", size = 3)+
  scale_fill_gradient2(low = "blue", mid = "pink", high = "red", midpoint = mid, name = fst.label, limits = c(0, max(fst.df$Fst)), breaks = c(0, 0.05, 0.10, 0.15))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), position = "right")+
  theme(axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 10)
  )

cololej_fst_boot = boot.ppfst(cololej.pop)
str(cololej_fst_boot)

cololej_fst_boot$vc.per.loc
cololej_fst_boot$ll[is.na(cololej_fst_boot$ll)] = 00
cololej_fst_boot$ul[is.na(cololej_fst_boot$ul)] = 00



cololej_fst %>% round(digits = 5)
cololej_fst_boot %>% round(digits = 4)
#Visualize results
# Convert dist object to data.frame
fst.matrix = as.matrix(cololej_fst)
ind = which( upper.tri(fst.matrix), arr.ind = TRUE)
fst.df = data.frame(Site1 = dimnames(fst.matrix)[[2]][ind[,2]],
                    Site2 = dimnames(fst.matrix)[[1]][ind[,1]],
                    Fst = fst.matrix[ ind ] %>% round(digits = 3))
fst.df$class = rep("Fst", 55)

fst.matrix_boot = as.matrix(cololej_fst_boot$ul)
ind = which( upper.tri(fst.matrix_boot), arr.ind = TRUE)
fst.df_boot = data.frame(Site1 = dimnames(fst.matrix_boot)[[2]][ind[,2]],
                         Site2 = dimnames(fst.matrix_boot)[[1]][ind[,1]],
                         Fst = fst.matrix_boot[ ind ] %>% round(digits = 3))

fst.df_boot$class = rep("upperlim", 55)
fst.df_boot_up = fst.df_boot
fst.df_boot_low = fst.df_boot

fst.df.com = data.frame(rbind(fst.df,fst.df_boot_low, fst.df_boot_up))
# Convert minus values to zero
fst.df$Fst[fst.df$Fst < 0] = 0

# Print data.frame summary
fst.df %>% str

# Fst italic label
fst.label = expression(italic("F")[ST])

# Extract middle Fst value for gradient argument
mid = max(fst.df$Fst) / 2

# Plot heatmap
ggplot(data = fst.df, aes(x = Site1, y = Site2, fill = Fst))+
  geom_tile(colour = "black")+
  geom_text(aes(label = Fst), color="black", size = 4)+
  scale_fill_gradient2(low = "blue", mid = "pink", high = "red", midpoint = mid, name = fst.label, limits = c(0, max(fst.df$Fst)), breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), position = "right")+
  theme(axis.text = element_text(colour = "black", size = 12),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12)
  )


###
str(fst.df.com)


fst.df.com$Sitescomp=as.factor(paste(fst.df.com$Site1,fst.df.com$Site2))

head(fst.df.com)
#View(fst.df.com)
addcomparison = function(x){
  if(x=="col10 col1"|x=="dim1 col1"|x=="dim1 col10"|x=="dim10 col1"|x=="dim10 col10"|x=="dim10 dim1"|x=="pa1 col1"|x=="pa1 col10"|x=="pa1 dim1"|x=="pa1 dim10"|x=="pa10 col1"|x=="pa10 col10"|x=="pa10 dim1"|x=="pa10 dim10"|x=="pa10 pa1") y = "Small-Small"
  if(x=="col1 Cabofrio"|x=="col10 Cabofrio"|x=="dim1 Cabofrio"|x=="dim10 Cabofrio"|x=="dim100 col1"|x=="dim100 col10"|x=="dim100 dim1"|x=="dim100 dim10"|x=="dimcont col1"|x=="dimcont col10"|x=="dimcont dim1"|x=="dimcont dim10"|x=="flor col1"|x=="flor col10"|x=="flor dim1"|x=="flor dim10"|x=="km41 col1"|x=="km41 col10"|x=="km41 dim1"|x=="km41 dim10"|x=="pa1 Cabofrio"|x=="pa1 dim100"|x=="pa1 dimcont"|x=="pa1 flor"|x=="pa1 km41"|x=="pa10 Cabofrio"|x=="pa10 dim100"|x=="pa10 dimcont"|x=="pa10 flor"|x=="pa10 km41") y = "Small-Large"
  if(x=="dim100 Cabofrio"|x=="dimcont Cabofrio"|x=="dimcont dim100"|x=="flor Cabofrio"| x=="flor dim100"|x=="flor dimcont"|x=="km41 Cabofrio"|x=="km41 dim100"|x=="km41 dimcont"|x=="km41 flor") y = "Large-Large"
  return(y)
}

str(fst.df.com)
fst.df.com$comp = as.factor(sapply(c(fst.df.com$Sitescomp), addcomparison))  

levels(fst.df.com$comp)
levels(fst.df.com$comp) = c("Small-Small", "Small-Large","Large-Large")
head(fst.df.com)
#View(fst.df.com)

str(fst.df.com)
fst.df.com$class<-as.factor(fst.df.com$class)
uplow = rep("Upper_lower_limits", 110)
manovel = rep("Mean value", 55)
fst.df.com$values = as.factor(c(manovel, uplow))
#######

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
# sample size
library(dplyr)

fst.df.com.boot = fst.df.com[c(56:165),]
fst.df.com.mean = fst.df.com[c(1:55),]


sample_size = fst.df.com.mean %>% group_by(comp) %>% summarize(num=n())
#install.packages("see")
library(see)
fst.df.com.mean = fst.df.com.mean %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(comp, "\n", "n=", num))

fst.df.com.mean$myaxis = as.factor(c(rep("Small-Small\nn=10",10), rep("Small-Large\nn=30",30), rep("Large-Large\nn=15",15)))
levels(fst.df.com.mean$myaxis) = c("Small-Small\nn=10", "Small-Large\nn=30", "Large-Large\nn=15")

aggregate(fst.df.com.mean$Fst, list(fst.df.com.mean$comp), FUN=mean) 
aggregate(fst.df.com.boot$Fst, list(fst.df.com.boot$comp), FUN=mean) 


write.table(fst.df.com, "cololej_80R20_fstdf_oct22.csv", sep = ",")

write.table(fst.df.com, "cololej_105R15_fstdf_oct22.csv", sep = ",")

plotfst.boot<- ggplot(fst.df.com.boot, aes(x=myaxis, y=Fst)) +
  #geom_violin(trim=FALSE, width=0.4, fill="#D95F02")+
  geom_violindot(trim=FALSE, width=0.6, fill="#D95F02", size_dots=0.20) +
  #geom_boxplot(width=0.1)+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  labs(x="Forest fragment size comparison", y =fst.label)+
  stat_summary(fun.data = data_summary, size = 1.5)+
  theme_classic()
plotfst.boot

plotfst<-ggplot(fst.df.com.mean, aes(x=myaxis, y=Fst)) +
  #geom_violin(trim=FALSE, width=0.6, fill="#1B9E77")+
  geom_violindot(trim=FALSE, width=0.6, fill="#1B9E77", size_dots=0.20) +
  #geom_boxplot(width=0.4, fill="darkgrey")+
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+
  labs(x=" ", y = fst.label)+
  scale_color_brewer(palette="Dark2")+
  stat_summary(fun.data = data_summary, size = 1.5)+
  theme_classic()
plotfst

figure.fst <- ggarrange(plotfst, plotfst.boot,
                        labels = c("A", "B"), hjust= 0, vjust = 0.85,
                        font.label = list(size = 9, color = "black", family = NULL),
                        ncol = 1, nrow = 2)


##AMOVA
cololej.pop
strata(cololej.pop) <- strata(cololej.pop)[,4:2]
cololej.pop.bis <- as.genclone(cololej.pop)
cololej.pop.bis

amova.result <- poppr.amova(cololej.pop.bis, ~sizeClass/size/subpop, missing = 'zero')

amova.test <- randtest(amova.result, nrepet = 999) # Test for significance
plot(amova.test)