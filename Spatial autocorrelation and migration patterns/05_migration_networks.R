###################################################################################
#                Spatial autocorrelation (Isolation by Distance)
#                                 Mantel test
#
#                               Adriel M. Sierra
#
#
#                                   2022
#
###################################################################################

rm(list=ls())
####Directory
setwd("/Users/adrielsierra/Documents/MS publication/Epiphyllous_fragment_pop_gen/manus_analisis")

# Read genind objects for both species
#Radula flaccida
radula.pop.105 = readRDS(file = "Radula_105R15_genind.rds")
radula.popmap.105 = read.table("Radula_105_popmap_com.txt", header=T)

#Cololejeunea surinamensis
cololej.pop.107 = readRDS(file = "Cololejeunea_107R15_genind.rds")
cololej.popmap.108 = read.table("Cololej_108_popmap_com.txt", header=T)

# install devtools
#install.packages("devtools")
# load devtools
library(devtools)
# load diveRsity
library(diveRsity)
library(graph4lg)
# install the dev version of diveRsity
#install_github("diveRsity", "kkeenan02")
#install.packages("remotes")
#remotes::install_github("spflanagan/gwscaR")
#install.packages("remotes")
#remotes::install_github("romunov/zvau")
library(zvau)
library(gwscaR)
library(adegenet)
library(dartR)

####Radula migration
radula.gl = gi2gl(radula.pop, parallel = FALSE, verbose = NULL)

radula.genepop = gl2genepop(radula.gl, outfile="radula.gen", outpath=getwd())

diffCalc(radula.genepop,  fst = TRUE, pairwise = TRUE,  bs_locus = TRUE, 
         bs_pairwise = TRUE, boots = 1000,  ci_type = "individuals", alpha = 0.05, para = TRUE)

reserves = c("Cabo frio",	"Colosso 1 ha",	"Colosso 10 ha",	"Dimona 1 ha",	"Dimona 10 ha",	"Dimona 100 ha",	"Dimona continuous",	"Florestal continuous",	"Km 41 continuous",	"Porto Alegre 1 ha", "Porto Alegre 10 ha")        
###calculate migration patterns with divMigrate
rad.migr.gst <- divMigrate("./radula.gen", plot_network = TRUE, 
                           boots = 1000, para = TRUE, stat = "gst")

rad.migr.gst <- divMigrate("./radula.gen", plot_network = TRUE, 
                           boots = 1000, para = TRUE, stat = "gst", filter_threshold = 0.35)


rad.migr.gst.df = as.matrix(rad.migr.gst$gRelMig)

rad.migr.gst.df.tsh = as.matrix(rad.migr.gst$gRelMigSig)

colnames(rad.migr.gst.df.tsh) = reserves
rownames(rad.migr.gst.df.tsh) = reserves

#library(RColorBrewer)
#display.brewer.all()
#brewer.pal(5, "Dark2")
#edit network graph
qgraph::qgraph(rad.migr.gst.df.tsh, nodeNames = colnames(rad.migr.gst.df.tsh), 
               color = c("#E7298A", "#1B9E77", "#D95F02", 
                         "#1B9E77", "#D95F02", "#7570B3","#E7298A",
                         "#E7298A","#E7298A","#1B9E77", "#D95F02"), 
               posCol = "#66A61E",
               legend = TRUE, edge.labels = TRUE, 
               curve = 3, width=5, height=5, mar = c(5,2,4,5))


###Nm
rad.migr.Nm <- divMigrate("./radula.gen", plot_network = TRUE, 
                          boots = 1000, para = TRUE, stat = "Nm", filter_threshold = 0.35)

rad.migr.Nm.df = as.matrix(rad.migr.Nm$nmRelMig)
View(rad.migr.Nm.df)
rad.migr.Nm.df.tsh = as.matrix(rad.migr.Nm$nmRelMigSig)

colnames(rad.migr.Nm.df) = reserves
rownames(rad.migr.Nm.df) = reserves

#library(RColorBrewer)
#display.brewer.all()
#brewer.pal(5, "Dark2")
#edit network graph
qgraph::qgraph(rad.migr.Nm.df.tsh, nodeNames = colnames(rad.migr.Nm.df.tsh), 
               color = c("#E7298A", "#1B9E77", "#D95F02", 
                         "#1B9E77", "#D95F02", "#7570B3","#E7298A",
                         "#E7298A","#E7298A","#1B9E77", "#D95F02"), 
               posCol = "#66A61E",
               legend = TRUE, edge.labels = TRUE, 
               curve = 3, width=5, height=5, mar = c(5,2,4,5))

####Cololejeunea migration
cololej.gl = gi2gl(cololej.pop, parallel = FALSE, verbose = NULL)

cololej.genepop = gl2genepop(cololej.gl, outfile="cololej.gen", outpath=getwd())

diffCalc(cololej.genepop,  fst = TRUE, pairwise = TRUE,  bs_locus = TRUE, 
         bs_pairwise = TRUE, boots = 1000,  ci_type = "individuals", alpha = 0.05, para = TRUE)





reserves = c("Cabo frio",	"Colosso 1 ha",	"Colosso 10 ha",	"Dimona 1 ha",	"Dimona 10 ha",	"Dimona 100 ha",	"Dimona continuous",	"Florestal continuous",	"Km 41 continuous",	"Porto Alegre 1 ha", "Porto Alegre 10 ha")        

###calculate migration patterns with divMigrate
col.migr.gst <- divMigrate("./cololej.gen", plot_network = TRUE, 
                           boots = 1000, para = TRUE, stat = "gst")

col.migr.gst <- divMigrate("./cololej.gen", plot_network = TRUE, 
                  boots = 1000, para = TRUE, stat = "gst", filter_threshold = 0)

View(col.migr.gst)
col.migr.gst.df = as.matrix(col.migr.gst$gRelMig)
View(col.migr.gst.df)
col.migr.gst.df.tsh = as.matrix(col.migr.gst$gRelMigSig)

colnames(col.migr.gst.df) = reserves
rownames(col.migr.gst.df) = reserves

#library(RColorBrewer)
#display.brewer.all()
#brewer.pal(5, "Dark2")
#edit network graph
qgraph::qgraph(col.migr.gst.df, nodeNames = colnames(col.migr.gst.df), 
               color = c("#E7298A", "#1B9E77", "#D95F02", 
                         "#1B9E77", "#D95F02", "#7570B3","#E7298A",
                        "#E7298A","#E7298A","#1B9E77", "#D95F02"), 
               posCol = "#66A61E",
               legend = TRUE, edge.labels = TRUE, 
               curve = 3, width=5, height=5, mar = c(5,2,4,5))


###Nm
col.migr.Nm <- divMigrate("./cololej.gen", plot_network = TRUE, 
                      boots = 1000, para = TRUE, stat = "Nm", filter_threshold = 0.35)

col.migr.Nm.df = as.matrix(col.migr.Nm$nmRelMig)
View(col.migr.Nm.df)
col.migr.Nm.df.tsh = as.matrix(col.migr.Nm$nmRelMigSig)

colnames(col.migr.Nm.df.tsh) = reserves
rownames(col.migr.Nm.df.tsh) = reserves

#library(RColorBrewer)
#display.brewer.all()
#brewer.pal(5, "Dark2")
#edit network graph
qgraph::qgraph(col.migr.Nm.df.tsh, nodeNames = colnames(col.migr.Nm.df.tsh), 
               color = c("#E7298A", "#1B9E77", "#D95F02", 
                         "#1B9E77", "#D95F02", "#7570B3","#E7298A",
                         "#E7298A","#E7298A","#1B9E77", "#D95F02"), 
               posCol = "#66A61E",
               legend = TRUE, edge.labels = TRUE, 
               curve = 3, width=5, height=5, mar = c(5,2,4,5))

