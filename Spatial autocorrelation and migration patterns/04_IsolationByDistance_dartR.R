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
radula.pop.70 = readRDS(file = "Radula_70R20_genind.rds")
radula.popmap.105 = read.table("Radula_105_popmap_com.txt", header=T)
radula.popmap.70 = read.table("Radula_71_popmap_com.txt", header=T)

#Cololejeunea surinamensis
cololej.pop.107 = readRDS(file = "Cololejeunea_107R15_genind.rds")
cololej.pop.80 = readRDS(file = "Cololejeunea_80R20_genind.rds")
cololej.popmap.108 = read.table("Cololej_108_popmap_com.txt", header=T)
cololej.popmap.71 = read.table("Cololej_71_popmap_com.txt", header=T)

#Libraries
library("spaa")
library(dartR)
library(ggplot2)
library("ggpubr")
library("adegenet")
library("hierfstat")
library(diveRsity)
library(pegas)
library(poppr)
library(hierfstat)

### Radula 107 dataset
#transform objects
radula.pop.gl = gi2gl(radula.pop.105)
latlong = radula.popmap.105[,6:7]
radula.pop.gl@other$latlong <- data.frame(radula.popmap.105[,6:7])
radula.pop.gl@other$size <- data.frame(radula.popmap.105[,5])

ibdall <- gl.ibd(radula.pop.gl)
str(ibdall)
plot(ibdall$Dgeo, ibdall$Dgen)

col.ibdmatrix = as.matrix(ibdall$Dgen)
ind = which( upper.tri(col.ibdmatrix), arr.ind = TRUE)
col.ibddf = data.frame(Site1 = dimnames(col.ibdmatrix)[[2]][ind[,2]],
                    Site2 = dimnames(col.ibdmatrix)[[1]][ind[,1]],
                    radula.ibd = col.ibdmatrix[ ind ] %>% round(digits = 3))

dgeo.matrix = as.matrix(ibdall$Dgeo)
ind = which( upper.tri(dgeo.matrix), arr.ind = TRUE)
dgeo.df = data.frame(Site1 = dimnames(dgeo.matrix)[[2]][ind[,2]],
                     Site2 = dimnames(dgeo.matrix)[[1]][ind[,1]],
                     radula.dgeo = dgeo.matrix[ ind ] %>% round(digits = 3))
dgeo.df
col.ibddf$geodist = dgeo.df$radula.dgeo


col.ibddf$Sitescomp =as.factor(paste(col.ibddf$Site1,col.ibddf$Site2))

head(col.ibddf)
#View(fst.df.com)
addcomparison = function(x){
  if(x=="col10 col1"|x=="dim1 col1"|x=="dim1 col10"|x=="dim10 col1"|x=="dim10 col10"|x=="dim10 dim1"|x=="pa1 col1"|x=="pa1 col10"|x=="pa1 dim1"|x=="pa1 dim10"|x=="pa10 col1"|x=="pa10 col10"|x=="pa10 dim1"|x=="pa10 dim10"|x=="pa10 pa1") y = "Small-Small"
  if(x=="col1 Cabofrio"|x=="col10 Cabofrio"|x=="dim1 Cabofrio"|x=="dim10 Cabofrio"|x=="dim100 col1"|x=="dim100 col10"|x=="dim100 dim1"|x=="dim100 dim10"|x=="dimcont col1"|x=="dimcont col10"|x=="dimcont dim1"|x=="dimcont dim10"|x=="flor col1"|x=="flor col10"|x=="flor dim1"|x=="flor dim10"|x=="km41 col1"|x=="km41 col10"|x=="km41 dim1"|x=="km41 dim10"|x=="pa1 Cabofrio"|x=="pa1 dim100"|x=="pa1 dimcont"|x=="pa1 flor"|x=="pa1 km41"|x=="pa10 Cabofrio"|x=="pa10 dim100"|x=="pa10 dimcont"|x=="pa10 flor"|x=="pa10 km41") y ="Small-Large"
  if(x=="dim100 Cabofrio"|x=="dimcont Cabofrio"|x=="dimcont dim100"|x=="flor Cabofrio"| x=="flor dim100"|x=="flor dimcont"|x=="km41 Cabofrio"|x=="km41 dim100"|x=="km41 dimcont"|x=="km41 flor") y = "Large-Large"
  return(y)
}

str(col.ibddf)
col.ibddf$classes = as.factor(sapply(c(col.ibddf$Sitescomp), addcomparison))  

levels(col.ibddf$classes)
#levels(col.ibddf$classes) = c("Small-Small", "Small-Large","Large-Large")
head(col.ibddf)

radsmall.lar = col.ibddf$classes == "Small-Large"
table(radsmall.lar)
radula.col.ibdsmall<-col.ibddf[radsmall.lar,]
radula.col.ibdsmall
library(vegan)
#mantel(radula.col.ibdsmall$radula.ibd, radula.col.ibdsmall$geodist, permutations = 999)



matgen = dist2list(ibdall$Dgen)
matgeo = dist2list(ibdall$Dgeo)
matgen$geo = matgeo$value

matgen = matgen[c(2:5,10:12,17:20,23,28:31,34,39:42,45,50:53,57:60,65:66,68:71,76:77,79:82,87:88,90:93,98:100,105:108,111,116:119),]

sl.gen = list2dist(matgen[,c(1:3)])
sl.geo = list2dist(matgen[,c(1:2,4)])
ibdsmall.large = gl.ibd(Dgen = sl.gen, Dgeo = sl.geo)
mantel(sl.gen, sl.geo, permutations = 999, na.rm = TRUE)
#
ibdsmall.large.matrix = as.matrix(ibdsmall.large$Dgen)
ind = which( upper.tri(ibdsmall.large.matrix), arr.ind = TRUE)
ibdsmall.large.df = data.frame(Site1 = dimnames(ibdsmall.large.matrix)[[2]][ind[,2]],
                                  Site2 = dimnames(ibdsmall.large.matrix)[[1]][ind[,1]],
                                  radula.col.ibdsmall = ibdsmall.large.matrix[ ind ] %>% round(digits = 3))

small.largedgeo.matrix = as.matrix(ibdsmall.large$Dgeo)
ind = which( upper.tri(small.largedgeo.matrix), arr.ind = TRUE)
small.largedgeo.df = data.frame(Site1 = dimnames(small.largedgeo.matrix)[[2]][ind[,2]],
                                   Site2 = dimnames(small.largedgeo.matrix)[[1]][ind[,1]],
                                   radula.smalldgeo = small.largedgeo.matrix[ ind ] %>% round(digits = 3))
small.largedgeo.df
ibdsmall.large.df$geodist = small.largedgeo.df$radula.smalldgeo
#
isoldistsmall.large<- ggplot(ibdsmall.large.df, aes(x=geodist, y=radula.col.ibdsmall)) +
  #geom_errorbar(aes(ymin=VarPi-VarPiStdErr, ymax=VarPi+VarPiStdErr), colour="black", width=.1, position=pd) +
  geom_point(size=3)+ 
  labs(x="log (Geographic distance)", y ="Fst/1-Fst")+
  geom_smooth(method=lm , color="Black", fill="lightgrey", se=TRUE)+
  theme_classic()+
  theme(legend.justification=c(),
        legend.position=c()) 
isoldistsmall.large


isoldist<- ggplot(col.ibddf, aes(x=geodist, y=radula.ibd, fill=classes)) +
  #geom_errorbar(aes(ymin=VarPi-VarPiStdErr, ymax=VarPi+VarPiStdErr), colour="black", width=.1, position=pd) +
  #geom_point()+
  geom_point(aes(shape = classes, col=classes), size=2.5)+ 
  labs(x="log (Geographic distance)", y ="Fst/1-Fst")+
  scale_color_brewer(palette="Dark2")+
  geom_smooth(aes(col=classes), method=lm, fill="lightgrey", se=TRUE, level=0.95)+
  theme_classic()+
  theme(legend.justification=c(),
        legend.position=c()) 
isoldist
###
### subset small
radsmall = radula.pop.gl@other$size == "Small"
table(radsmall)
radula.pop.gl.small<-radula.pop.gl[radsmall,]
radula.pop.gl.small@other$size
#
ibdsmall = gl.ibd(radula.pop.gl.small)
#
ibdsmall.matrix = as.matrix(ibdsmall$Dgen)
ind = which( upper.tri(ibdsmall.matrix), arr.ind = TRUE)
ibdsmall.df = data.frame(Site1 = dimnames(ibdsmall.matrix)[[2]][ind[,2]],
                         Site2 = dimnames(ibdsmall.matrix)[[1]][ind[,1]],
                         radula.col.ibdsmall = ibdsmall.matrix[ ind ] %>% round(digits = 3))

smalldgeo.matrix = as.matrix(ibdsmall$Dgeo)
ind = which( upper.tri(smalldgeo.matrix), arr.ind = TRUE)
smalldgeo.df = data.frame(Site1 = dimnames(smalldgeo.matrix)[[2]][ind[,2]],
                          Site2 = dimnames(smalldgeo.matrix)[[1]][ind[,1]],
                          radula.smalldgeo = smalldgeo.matrix[ ind ] %>% round(digits = 3))
smalldgeo.df
ibdsmall.df$geodist = smalldgeo.df$radula.smalldgeo

ibdsmall.df$sizeplot = as.factor(c("1 ha vs 10 ha", "1 ha vs 1 ha", "1 ha vs 10 ha", "1 ha vs 10 ha", "10 ha vs 10 ha", "1 ha vs 10 ha", "1 ha vs 1 ha", "1 ha vs 10 ha", "1 ha vs 1 ha", "1 ha vs 10 ha", "1 ha vs 10 ha", "10 ha vs 10 ha", "1 ha vs 10 ha", "10 ha vs 10 ha", "1 ha vs 10 ha"))
str(ibdsmall.df)
isoldistsmall<- ggplot(ibdsmall.df, aes(x=geodist, y=radula.col.ibdsmall)) +
  #geom_errorbar(aes(ymin=VarPi-VarPiStdErr, ymax=VarPi+VarPiStdErr), colour="black", width=.1, position=pd) +
  geom_point( size=3)+ 
  labs(x="log (Geographic distance)", y ="Fst/1-Fst")+
  geom_smooth(method=lm , color="Black", fill="lightgrey", se=TRUE)+
  theme_classic()+
  theme(legend.justification=c(),
        legend.position=c()) 
isoldistsmall

###subset large
radlarge = radula.pop.gl@other$size == "Large"
table(radlarge)
radula.pop.gl.large<-radula.pop.gl[radlarge,]
radula.pop.gl.large@other$size
ibdlarge = gl.ibd(radula.pop.gl.large)

#
ibdlarge.matrix = as.matrix(ibdlarge$Dgen)
ind = which( upper.tri(ibdlarge.matrix), arr.ind = TRUE)
ibdlarge.df = data.frame(Site1 = dimnames(ibdlarge.matrix)[[2]][ind[,2]],
                         Site2 = dimnames(ibdlarge.matrix)[[1]][ind[,1]],
                         radulalarge.ibd = ibdlarge.matrix[ ind ] %>% round(digits = 3))

ibdgeolarge.matrix = as.matrix(ibdlarge$Dgeo)
ind = which( upper.tri(ibdgeolarge.matrix), arr.ind = TRUE)
ibdgeolarge.df = data.frame(Site1 = dimnames(ibdgeolarge.matrix)[[2]][ind[,2]],
                          Site2 = dimnames(ibdgeolarge.matrix)[[1]][ind[,1]],
                          radulalarge.dgeo = ibdgeolarge.matrix[ ind ] %>% round(digits = 3))
ibdgeolarge.df
ibdlarge.df$geodist = ibdgeolarge.df$radulalarge.dgeo


isoldistlarge<- ggplot(ibdlarge.df, aes(x=geodist, y=radulalarge.ibd)) +
  #geom_errorbar(aes(ymin=VarPi-VarPiStdErr, ymax=VarPi+VarPiStdErr), colour="black", width=.1, position=pd) +
  geom_point(size=3)+ labs(x="log (Geographic distance)", y ="Fst/1-Fst")+
  geom_smooth(method=lm , color="Black", fill="lightgrey", se=TRUE)+
  theme_classic()+
  theme(legend.justification=c(),
        legend.position=c()) 
isoldistlarge

### Radula 70 dataset
#
#


#transform objects
radula.pop.gl.70 = gi2gl(radula.pop.70)
latlong = radula.popmap.70[,6:7]
radula.pop.gl.70@other$latlong <- data.frame(radula.popmap.70[,6:7])
radula.pop.gl.70@other$size <- data.frame(radula.popmap.70[,5])

ibdall.70 <- gl.ibd(radula.pop.gl.70)
str(ibdall.70)
plot(ibdall.70$Dgeo, ibdall.70$Dgen)


ibd.matrix.70 = as.matrix(ibdall.70$Dgen)
ind = which( upper.tri(ibd.matrix.70), arr.ind = TRUE)
ibd.df.70 = data.frame(Site1 = dimnames(ibd.matrix.70)[[2]][ind[,2]],
                    Site2 = dimnames(ibd.matrix.70)[[1]][ind[,1]],
                    radula.ibd = ibd.matrix.70[ ind ] %>% round(digits = 3))

dgeo.matrix.70 = as.matrix(ibdall.70$Dgeo)
ind = which( upper.tri(dgeo.matrix), arr.ind = TRUE)
dgeo.df.70 = data.frame(Site1 = dimnames(dgeo.matrix.70)[[2]][ind[,2]],
                     Site2 = dimnames(dgeo.matrix.70)[[1]][ind[,1]],
                     radula.dgeo = dgeo.matrix.70[ ind ] %>% round(digits = 3))
dgeo.df.70
ibd.df.70$geodist = dgeo.df.70$radula.dgeo


ibd.df.70$Sitescomp =as.factor(paste(ibd.df.70$Site1,ibd.df.70$Site2))

head(ibd.df.70)
#View(fst.df.com)
addcomparison = function(x){
  if(x=="col10 col1"|x=="dim1 col1"|x=="dim1 col10"|x=="dim10 col1"|x=="dim10 col10"|x=="dim10 dim1"|x=="pa1 col1"|x=="pa1 col10"|x=="pa1 dim1"|x=="pa1 dim10"|x=="pa10 col1"|x=="pa10 col10"|x=="pa10 dim1"|x=="pa10 dim10"|x=="pa10 pa1") y = "Small-Small"
  if(x=="col1 Cabofrio"|x=="col10 Cabofrio"|x=="dim1 Cabofrio"|x=="dim10 Cabofrio"|x=="dim100 col1"|x=="dim100 col10"|x=="dim100 dim1"|x=="dim100 dim10"|x=="dimcont col1"|x=="dimcont col10"|x=="dimcont dim1"|x=="dimcont dim10"|x=="flor col1"|x=="flor col10"|x=="flor dim1"|x=="flor dim10"|x=="km41 col1"|x=="km41 col10"|x=="km41 dim1"|x=="km41 dim10"|x=="pa1 Cabofrio"|x=="pa1 dim100"|x=="pa1 dimcont"|x=="pa1 flor"|x=="pa1 km41"|x=="pa10 Cabofrio"|x=="pa10 dim100"|x=="pa10 dimcont"|x=="pa10 flor"|x=="pa10 km41") y ="Small-Large"
  if(x=="dim100 Cabofrio"|x=="dimcont Cabofrio"|x=="dimcont dim100"|x=="flor Cabofrio"| x=="flor dim100"|x=="flor dimcont"|x=="km41 Cabofrio"|x=="km41 dim100"|x=="km41 dimcont"|x=="km41 flor") y = "Large-Large"
  return(y)
}

str(ibd.df.70)
ibd.df.70$classes = as.factor(sapply(c(ibd.df.70$Sitescomp), addcomparison))  

levels(ibd.df.70$classes)
#levels(ibd.df$classes) = c("Small-Small", "Small-Large","Large-Large")
head(ibd.df.70)



matgen.70 = dist2list(ibdall.70$Dgen)
matgeo.70 = dist2list(ibdall.70$Dgeo)
matgen.70$geo = matgeo.70$value

matgen.70 = matgen.70[c(2:5,10:12,17:20,23,28:31,34,39:42,45,50:53,57:60,65:66,68:71,76:77,79:82,87:88,90:93,98:100,105:108,111,116:119),]

sl.gen.70 = list2dist(matgen.70[,c(1:3)])
sl.geo.70 = list2dist(matgen.70[,c(1:2,4)])
ibdsmall.large.70 = gl.ibd(Dgen = sl.gen.70, Dgeo = sl.geo.70)
mantel(sl.gen.70, sl.geo.70, permutations = 999, na.rm = TRUE)

#
ibdsmall.large.matrix.70 = as.matrix(ibdsmall.large.70$Dgen)
ind = which( upper.tri(ibdsmall.large.matrix.70), arr.ind = TRUE)
ibdsmall.large.df.70 = data.frame(Site1 = dimnames(ibdsmall.large.matrix.70)[[2]][ind[,2]],
                            Site2 = dimnames(ibdsmall.large.matrix.70)[[1]][ind[,1]],
                            radula.ibd.small = ibdsmall.large.matrix.70[ ind ] %>% round(digits = 3))

small.largedgeo.matrix.70 = as.matrix(ibdsmall.large.70$Dgeo)
ind = which( upper.tri(small.largedgeo.matrix.70), arr.ind = TRUE)
small.largedgeo.df.70 = data.frame(Site1 = dimnames(small.largedgeo.matrix.70)[[2]][ind[,2]],
                             Site2 = dimnames(small.largedgeo.matrix.70)[[1]][ind[,1]],
                             radula.smalldgeo = small.largedgeo.matrix.70[ ind ] %>% round(digits = 3))
small.largedgeo.df.70
ibdsmall.large.df.70$geodist = small.largedgeo.df.70$radula.smalldgeo
#
isoldistsmall.large.70<- ggplot(ibdsmall.large.df.70, aes(x=geodist, y=radula.ibd.small)) +
  #geom_errorbar(aes(ymin=VarPi-VarPiStdErr, ymax=VarPi+VarPiStdErr), colour="black", width=.1, position=pd) +
  geom_point(size=3)+ 
  labs(x="log (Geographic distance)", y ="Fst/1-Fst")+
  geom_smooth(method=lm , color="Black", fill="lightgrey", se=TRUE)+
  theme_classic()+
  theme(legend.justification=c(),
        legend.position=c()) 
isoldistsmall.large.70


#all
isoldist.70<- ggplot(ibd.df.70, aes(x=geodist, y=radula.ibd, fill=classes)) +
  #geom_errorbar(aes(ymin=VarPi-VarPiStdErr, ymax=VarPi+VarPiStdErr), colour="black", width=.1, position=pd) +
  #geom_point()+
  geom_point(aes(shape = classes, col=classes), size=2.5)+ 
  labs(x="log (Geographic distance)", y ="Fst/1-Fst")+
  scale_color_brewer(palette="Dark2")+
  geom_smooth(aes(col=classes), method=lm, fill="lightgrey", se=TRUE, level=0.95)+
  theme_classic()+
  theme(legend.justification=c(),
        legend.position=c()) 
isoldist.70
###
### subset small
radsmall.70 = radula.pop.gl.70@other$size == "Small"
table(radsmall.70)
radula.pop.gl.small.70<-radula.pop.gl.70[radsmall.70,]
radula.pop.gl.small.70@other$size
#
ibdsmall.70 = gl.ibd(radula.pop.gl.small.70)
#
ibdsmall.matrix.70 = as.matrix(ibdsmall.70$Dgen)
ind = which( upper.tri(ibdsmall.matrix.70), arr.ind = TRUE)
ibdsmall.df.70 = data.frame(Site1 = dimnames(ibdsmall.matrix.70)[[2]][ind[,2]],
                         Site2 = dimnames(ibdsmall.matrix.70)[[1]][ind[,1]],
                         radula.ibd.small = ibdsmall.matrix.70[ ind ] %>% round(digits = 3))

smalldgeo.matrix.70 = as.matrix(ibdsmall.70$Dgeo)
ind = which( upper.tri(smalldgeo.matrix), arr.ind = TRUE)
smalldgeo.df.70 = data.frame(Site1 = dimnames(smalldgeo.matrix.70)[[2]][ind[,2]],
                          Site2 = dimnames(smalldgeo.matrix.70)[[1]][ind[,1]],
                          radula.smalldgeo = smalldgeo.matrix.70[ ind ] %>% round(digits = 3))
smalldgeo.df.70
ibdsmall.df.70$geodist = smalldgeo.df.70$radula.smalldgeo

ibdsmall.df.70$sizeplot = as.factor(c("1 ha vs 10 ha", "1 ha vs 1 ha", "1 ha vs 10 ha", "1 ha vs 10 ha", "10 ha vs 10 ha", "1 ha vs 10 ha", "1 ha vs 1 ha", "1 ha vs 10 ha", "1 ha vs 1 ha", "1 ha vs 10 ha", "1 ha vs 10 ha", "10 ha vs 10 ha", "1 ha vs 10 ha", "10 ha vs 10 ha", "1 ha vs 10 ha"))
str(ibdsmall.df.70)
isoldistsmall.70<- ggplot(ibdsmall.df.70, aes(x=geodist, y=radula.ibd.small)) +
  #geom_errorbar(aes(ymin=VarPi-VarPiStdErr, ymax=VarPi+VarPiStdErr), colour="black", width=.1, position=pd) +
  geom_point(size=3)+ 
  labs(x="log (Geographic distance)", y ="Fst/1-Fst")+
  geom_smooth(method=lm , color="Black", fill="lightgrey", se=TRUE)+
  theme_classic()+
  theme(legend.justification=c(),
        legend.position=c()) 
isoldistsmall.70



###subset large
radlarge.70 = radula.pop.gl.70@other$size == "Large"
table(radlarge.70)
radula.pop.gl.large.70<-radula.pop.gl.70[radlarge.70,]
radula.pop.gl.large.70@other$size
ibdlarge.70 = gl.ibd(radula.pop.gl.large.70)

#
ibdlarge.matrix.70 = as.matrix(ibdlarge.70$Dgen)
ind = which( upper.tri(ibdlarge.matrix.70), arr.ind = TRUE)
ibdlarge.df.70 = data.frame(Site1 = dimnames(ibdlarge.matrix.70)[[2]][ind[,2]],
                         Site2 = dimnames(ibdlarge.matrix.70)[[1]][ind[,1]],
                         radulalarge.ibd = ibdlarge.matrix.70[ ind ] %>% round(digits = 3))

ibdgeolarge.matrix.70 = as.matrix(ibdlarge.70$Dgeo)
ind = which( upper.tri(ibdgeolarge.matrix.70), arr.ind = TRUE)
ibdgeolarge.df.70 = data.frame(Site1 = dimnames(ibdgeolarge.matrix.70)[[2]][ind[,2]],
                          Site2 = dimnames(ibdgeolarge.matrix.70)[[1]][ind[,1]],
                          radulalarge.dgeo = ibdgeolarge.matrix.70[ ind ] %>% round(digits = 3))
ibdgeolarge.df.70
ibdlarge.df.70$geodist = ibdgeolarge.df.70$radulalarge.dgeo


isoldistlarge.70<- ggplot(ibdlarge.df.70, aes(x=geodist, y=radulalarge.ibd)) +
  #geom_errorbar(aes(ymin=VarPi-VarPiStdErr, ymax=VarPi+VarPiStdErr), colour="black", width=.1, position=pd) +
  geom_point(size=3)+ labs(x="log (Geographic distance)", y ="Fst/1-Fst")+
  geom_smooth(method=lm , color="Black", fill="lightgrey", se=TRUE)+
  theme_classic()+
  theme(legend.justification=c(),
        legend.position=c()) 
#plot list
isoldistlarge
isoldistsmall
isoldistsmall.large

isoldistlarge.70
isoldistsmall.70
isoldistsmall.large.70


### Cololejeunea 107 dataset
cololej.pop.107 = readRDS(file = "Cololejeunea_107R15_genind.rds")
cololej.pop.80 = readRDS(file = "Cololejeunea_80R20_genind.rds")
cololej.popmap.108 = read.table("Cololej_108_popmap_com.txt", header=T)
cololej.popmap.71
#transform objects
col.pop.gl = gi2gl(cololej.pop.107)
latlong = cololej.popmap.108[,6:7]
col.pop.gl@other$latlong <- data.frame(cololej.popmap.108[,6:7])
col.pop.gl@other$size <- data.frame(cololej.popmap.108[,5])

col.ibdall <- gl.ibd(col.pop.gl)
str(col.ibdall)
plot(col.ibdall$Dgeo, col.ibdall$Dgen)

col.ibdmatrix = as.matrix(col.ibdall$Dgen)
ind = which( upper.tri(col.ibdmatrix), arr.ind = TRUE)
col.ibddf = data.frame(Site1 = dimnames(col.ibdmatrix)[[2]][ind[,2]],
                    Site2 = dimnames(col.ibdmatrix)[[1]][ind[,1]],
                    col.ibd = col.ibdmatrix[ ind ] %>% round(digits = 3))

col.dgeo.matrix = as.matrix(col.ibdall$Dgeo)
ind = which( upper.tri(dgeo.matrix), arr.ind = TRUE)
col.dgeo.df = data.frame(Site1 = dimnames(dgeo.matrix)[[2]][ind[,2]],
                     Site2 = dimnames(dgeo.matrix)[[1]][ind[,1]],
                     col.dgeo = dgeo.matrix[ ind ] %>% round(digits = 3))
col.dgeo.df
col.ibddf$geodist = col.dgeo.df$col.dgeo


col.ibddf$Sitescomp =as.factor(paste(col.ibddf$Site1,col.ibddf$Site2))

head(col.ibddf)
#View(fst.df.com)
addcomparison = function(x){
  if(x=="col10 col1"|x=="dim1 col1"|x=="dim1 col10"|x=="dim10 col1"|x=="dim10 col10"|x=="dim10 dim1"|x=="pa1 col1"|x=="pa1 col10"|x=="pa1 dim1"|x=="pa1 dim10"|x=="pa10 col1"|x=="pa10 col10"|x=="pa10 dim1"|x=="pa10 dim10"|x=="pa10 pa1") y = "Small-Small"
  if(x=="col1 Cabofrio"|x=="col10 Cabofrio"|x=="dim1 Cabofrio"|x=="dim10 Cabofrio"|x=="dim100 col1"|x=="dim100 col10"|x=="dim100 dim1"|x=="dim100 dim10"|x=="dimcont col1"|x=="dimcont col10"|x=="dimcont dim1"|x=="dimcont dim10"|x=="flor col1"|x=="flor col10"|x=="flor dim1"|x=="flor dim10"|x=="km41 col1"|x=="km41 col10"|x=="km41 dim1"|x=="km41 dim10"|x=="pa1 Cabofrio"|x=="pa1 dim100"|x=="pa1 dimcont"|x=="pa1 flor"|x=="pa1 km41"|x=="pa10 Cabofrio"|x=="pa10 dim100"|x=="pa10 dimcont"|x=="pa10 flor"|x=="pa10 km41") y ="Small-Large"
  if(x=="dim100 Cabofrio"|x=="dimcont Cabofrio"|x=="dimcont dim100"|x=="flor Cabofrio"| x=="flor dim100"|x=="flor dimcont"|x=="km41 Cabofrio"|x=="km41 dim100"|x=="km41 dimcont"|x=="km41 flor") y = "Large-Large"
  return(y)
}

str(col.ibddf)
col.ibddf$classes = as.factor(sapply(c(col.ibddf$Sitescomp), addcomparison))  

levels(col.ibddf$classes)
#levels(col.ibddf$classes) = c("Small-Small", "Small-Large","Large-Large")
head(col.ibddf)

colsmall.lar = col.ibddf$classes == "Small-Large"
table(colsmall.lar)
col.col.ibdsmall<-col.ibddf[colsmall.lar,]
col.col.ibdsmall
#library(vegan)
#mantel(col.col.ibdsmall$col.ibd, col.col.ibdsmall$geodist, permutations = 999)



col.matgen = dist2list(col.ibdall$Dgen)
col.matgeo = dist2list(col.ibdall$Dgeo)
col.matgen$geo = col.matgeo$value

col.matgen = col.matgen[c(2:5,10:12,17:20,23,28:31,34,39:42,45,50:53,57:60,65:66,68:71,76:77,79:82,87:88,90:93,98:100,105:108,111,116:119),]

col.sl.gen = list2dist(col.matgen[,c(1:3)])
col.sl.geo = list2dist(col.matgen[,c(1:2,4)])
col.ibdsmall.large = gl.ibd(Dgen = col.sl.gen, Dgeo = col.sl.geo)
mantel(col.sl.gen, col.sl.geo, permutations = 999, na.rm = TRUE)
#
col.ibdsmall.large.matrix = as.matrix(col.ibdsmall.large$Dgen)
ind = which( upper.tri(col.ibdsmall.large.matrix), arr.ind = TRUE)
col.ibdsmall.large.df = data.frame(Site1 = dimnames(col.ibdsmall.large.matrix)[[2]][ind[,2]],
                               Site2 = dimnames(col.ibdsmall.large.matrix)[[1]][ind[,1]],
                               col.col.ibdsmall = col.ibdsmall.large.matrix[ ind ] %>% round(digits = 3))

col.small.largedgeo.matrix = as.matrix(col.ibdsmall.large$Dgeo)
ind = which( upper.tri(col.small.largedgeo.matrix), arr.ind = TRUE)
col.small.largedgeo.df = data.frame(Site1 = dimnames(col.small.largedgeo.matrix)[[2]][ind[,2]],
                                Site2 = dimnames(col.small.largedgeo.matrix)[[1]][ind[,1]],
                                col.smalldgeo = col.small.largedgeo.matrix[ ind ] %>% round(digits = 3))
col.small.largedgeo.df
col.ibdsmall.large.df$geodist = col.small.largedgeo.df$col.smalldgeo
#
col.isoldistsmall.large<- ggplot(col.ibdsmall.large.df, aes(x=geodist, y=col.col.ibdsmall)) +
  #geom_errorbar(aes(ymin=VarPi-VarPiStdErr, ymax=VarPi+VarPiStdErr), colour="black", width=.1, position=pd) +
  geom_point(size=3)+ 
  labs(x="log (Geographic distance)", y ="Fst/1-Fst")+
  geom_smooth(method=lm , color="Black", fill="lightgrey", se=TRUE)+
  theme_classic()+
  theme(legend.justification=c(),
        legend.position=c()) 
col.isoldistsmall.large


col.isoldist<- ggplot(col.ibddf, aes(x=geodist, y=col.ibd, fill=classes)) +
  #geom_errorbar(aes(ymin=VarPi-VarPiStdErr, ymax=VarPi+VarPiStdErr), colour="black", width=.1, position=pd) +
  #geom_point()+
  geom_point(aes(shape = classes, col=classes), size=2.5)+ 
  labs(x="log (Geographic distance)", y ="Fst/1-Fst")+
  scale_color_brewer(palette="Dark2")+
  geom_smooth(aes(col=classes), method=lm, fill="lightgrey", se=TRUE, level=0.95)+
  theme_classic()+
  theme(legend.justification=c(),
        legend.position=c()) 
col.isoldist
###
### subset small
colsmall = col.pop.gl@other$size == "Small"
table(colsmall)
col.pop.gl.small<-col.pop.gl[colsmall,]
col.pop.gl.small@other$size
#
col.ibdsmall = gl.ibd(col.pop.gl.small)
#
col.ibdsmall.matrix = as.matrix(col.ibdsmall$Dgen)
ind = which( upper.tri(ibdsmall.matrix), arr.ind = TRUE)
col.ibdsmall.df = data.frame(Site1 = dimnames(col.ibdsmall.matrix)[[2]][ind[,2]],
                         Site2 = dimnames(col.ibdsmall.matrix)[[1]][ind[,1]],
                         col.col.ibdsmall = col.ibdsmall.matrix[ ind ] %>% round(digits = 3))

col.smalldgeo.matrix = as.matrix(col.ibdsmall$Dgeo)
ind = which( upper.tri(col.smalldgeo.matrix), arr.ind = TRUE)
col.smalldgeo.df = data.frame(Site1 = dimnames(col.smalldgeo.matrix)[[2]][ind[,2]],
                          Site2 = dimnames(col.smalldgeo.matrix)[[1]][ind[,1]],
                          col.smalldgeo = col.smalldgeo.matrix[ ind ] %>% round(digits = 3))
col.smalldgeo.df
col.ibdsmall.df$geodist = col.smalldgeo.df$col.smalldgeo

col.ibdsmall.df$sizeplot = as.factor(c("1 ha vs 10 ha", "1 ha vs 1 ha", "1 ha vs 10 ha", "1 ha vs 10 ha", "10 ha vs 10 ha", "1 ha vs 10 ha", "1 ha vs 1 ha", "1 ha vs 10 ha", "1 ha vs 1 ha", "1 ha vs 10 ha", "1 ha vs 10 ha", "10 ha vs 10 ha", "1 ha vs 10 ha", "10 ha vs 10 ha", "1 ha vs 10 ha"))
str(ibdsmall.df)
col.isoldistsmall<- ggplot(col.ibdsmall.df, aes(x=geodist, y=col.col.ibdsmall)) +
  #geom_errorbar(aes(ymin=VarPi-VarPiStdErr, ymax=VarPi+VarPiStdErr), colour="black", width=.1, position=pd) +
  geom_point(size=3)+ 
  labs(x="log (Geographic distance)", y ="Fst/1-Fst")+
  geom_smooth(method=lm , color="Black", fill="lightgrey", se=TRUE)+
  theme_classic()+
  theme(legend.justification=c(),
        legend.position=c()) 
col.isoldistsmall

###subset large
collarge = col.pop.gl@other$size == "Large"
table(collarge)
col.pop.gl.large<-col.pop.gl[collarge,]
col.pop.gl.large@other$size
col.ibdlarge = gl.ibd(col.pop.gl.large)

#
col.ibdlarge.matrix = as.matrix(col.ibdlarge$Dgen)
ind = which( upper.tri(col.ibdlarge.matrix), arr.ind = TRUE)
col.ibdlarge.df = data.frame(Site1 = dimnames(col.ibdlarge.matrix)[[2]][ind[,2]],
                         Site2 = dimnames(col.ibdlarge.matrix)[[1]][ind[,1]],
                         collarge.ibd = col.ibdlarge.matrix[ ind ] %>% round(digits = 3))

col.dgeolarge.matrix = as.matrix(col.ibdlarge$Dgeo)
ind = which( upper.tri(col.dgeolarge.matrix), arr.ind = TRUE)
col.ibdgeolarge.df = data.frame(Site1 = dimnames(col.dgeolarge.matrix)[[2]][ind[,2]],
                            Site2 = dimnames(col.dgeolarge.matrix)[[1]][ind[,1]],
                            collarge.dgeo = col.dgeolarge.matrix[ ind ] %>% round(digits = 3))
col.ibdgeolarge.df
col.ibdlarge.df$geodist = col.ibdgeolarge.df$collarge.dgeo


col.isoldistlarge<- ggplot(col.ibdlarge.df, aes(x=geodist, y=collarge.ibd)) +
  #geom_errorbar(aes(ymin=VarPi-VarPiStdErr, ymax=VarPi+VarPiStdErr), colour="black", width=.1, position=pd) +
  geom_point(size=3)+ labs(x="log (Geographic distance)", y ="Fst/1-Fst")+
  geom_smooth(method=lm , color="Black", fill="lightgrey", se=TRUE)+
  theme_classic()+
  theme(legend.justification=c(),
        legend.position=c()) 
col.isoldistlarge

### col 70 dataset
#
#


#transform objects
col.pop.gl.70 = gi2gl(cololej.pop.80)
latlong = cololej.popmap.71[,6:7]
col.pop.gl.70@other$latlong <- data.frame(cololej.popmap.71[,6:7])
col.pop.gl.70@other$size <- data.frame(cololej.popmap.71[,5])

col.ibdall.70 <- gl.ibd(col.pop.gl.70)
str(col.ibdall.70)
plot(col.ibdall.70$Dgeo, col.ibdall.70$Dgen)


col.ibdmatrix.70 = as.matrix(col.ibdall.70$Dgen)
ind = which( upper.tri(col.ibdmatrix.70), arr.ind = TRUE)
col.ibddf.70 = data.frame(Site1 = dimnames(col.ibdmatrix.70)[[2]][ind[,2]],
                       Site2 = dimnames(col.ibdmatrix.70)[[1]][ind[,1]],
                       col.ibd = col.ibdmatrix.70[ ind ] %>% round(digits = 3))

col.dgeo.matrix.70 = as.matrix(col.ibdall.70$Dgeo)
ind = which( upper.tri(dgeo.matrix), arr.ind = TRUE)
col.dgeo.df.70 = data.frame(Site1 = dimnames(col.dgeo.matrix.70)[[2]][ind[,2]],
                        Site2 = dimnames(col.dgeo.matrix.70)[[1]][ind[,1]],
                        col.dgeo = col.dgeo.matrix.70[ ind ] %>% round(digits = 3))
col.dgeo.df.70
col.ibddf.70$geodist = col.dgeo.df.70$col.dgeo


col.ibddf.70$Sitescomp =as.factor(paste(col.ibddf.70$Site1,col.ibddf.70$Site2))

head(col.ibddf.70)
#View(fst.df.com)
addcomparison = function(x){
  if(x=="col10 col1"|x=="dim1 col1"|x=="dim1 col10"|x=="dim10 col1"|x=="dim10 col10"|x=="dim10 dim1"|x=="pa1 col1"|x=="pa1 col10"|x=="pa1 dim1"|x=="pa1 dim10"|x=="pa10 col1"|x=="pa10 col10"|x=="pa10 dim1"|x=="pa10 dim10"|x=="pa10 pa1") y = "Small-Small"
  if(x=="col1 Cabofrio"|x=="col10 Cabofrio"|x=="dim1 Cabofrio"|x=="dim10 Cabofrio"|x=="dim100 col1"|x=="dim100 col10"|x=="dim100 dim1"|x=="dim100 dim10"|x=="dimcont col1"|x=="dimcont col10"|x=="dimcont dim1"|x=="dimcont dim10"|x=="flor col1"|x=="flor col10"|x=="flor dim1"|x=="flor dim10"|x=="km41 col1"|x=="km41 col10"|x=="km41 dim1"|x=="km41 dim10"|x=="pa1 Cabofrio"|x=="pa1 dim100"|x=="pa1 dimcont"|x=="pa1 flor"|x=="pa1 km41"|x=="pa10 Cabofrio"|x=="pa10 dim100"|x=="pa10 dimcont"|x=="pa10 flor"|x=="pa10 km41") y ="Small-Large"
  if(x=="dim100 Cabofrio"|x=="dimcont Cabofrio"|x=="dimcont dim100"|x=="flor Cabofrio"| x=="flor dim100"|x=="flor dimcont"|x=="km41 Cabofrio"|x=="km41 dim100"|x=="km41 dimcont"|x=="km41 flor") y = "Large-Large"
  return(y)
}

str(col.ibddf.70)
col.ibddf.70$classes = as.factor(sapply(c(col.ibddf.70$Sitescomp), addcomparison))  

levels(col.ibddf.70$classes)
#levels(col.ibddf$classes) = c("Small-Small", "Small-Large","Large-Large")
head(col.ibddf.70)



col.matgen.70 = dist2list(col.ibdall.70$Dgen)
col.matgeo.70 = dist2list(col.ibdall.70$Dgeo)
col.matgen.70$geo = col.matgeo.70$value

col.matgen.70 = col.matgen.70[c(2:5,10:12,17:20,23,28:31,34,39:42,45,50:53,57:60,65:66,68:71,76:77,79:82,87:88,90:93,98:100,105:108,111,116:119),]

col.sl.gen.70 = list2dist(col.matgen.70[,c(1:3)])
col.sl.geo.70 = list2dist(col.matgen.70[,c(1:2,4)])
col.ibdsmall.large.70 = gl.ibd(Dgen = col.sl.gen.70, Dgeo = col.sl.geo.70)
mantel(col.sl.gen.70, col.sl.geo.70, permutations = 999, na.rm = TRUE)

#
col.ibdsmall.large.matrix.70 = as.matrix(col.ibdsmall.large.70$Dgen)
ind = which( upper.tri(col.ibdsmall.large.matrix.70), arr.ind = TRUE)
col.ibdsmall.large.df.70 = data.frame(Site1 = dimnames(col.ibdsmall.large.matrix.70)[[2]][ind[,2]],
                                  Site2 = dimnames(col.ibdsmall.large.matrix.70)[[1]][ind[,1]],
                                  col.col.ibdsmall = col.ibdsmall.large.matrix.70[ ind ] %>% round(digits = 3))

col.small.largedgeo.matrix.70 = as.matrix(col.ibdsmall.large.70$Dgeo)
ind = which( upper.tri(col.small.largedgeo.matrix.70), arr.ind = TRUE)
col.small.largedgeo.df.70 = data.frame(Site1 = dimnames(col.small.largedgeo.matrix.70)[[2]][ind[,2]],
                                   Site2 = dimnames(col.small.largedgeo.matrix.70)[[1]][ind[,1]],
                                   col.smalldgeo = col.small.largedgeo.matrix.70[ ind ] %>% round(digits = 3))
col.small.largedgeo.df.70
col.ibdsmall.large.df.70$geodist = col.small.largedgeo.df.70$col.smalldgeo
#
col.isoldistsmall.large.70<- ggplot(col.ibdsmall.large.df.70, aes(x=geodist, y=col.col.ibdsmall)) +
  #geom_errorbar(aes(ymin=VarPi-VarPiStdErr, ymax=VarPi+VarPiStdErr), colour="black", width=.1, position=pd) +
  geom_point(size=3)+ 
  labs(x="log (Geographic distance)", y ="Fst/1-Fst")+
  geom_smooth(method=lm , color="Black", fill="lightgrey", se=TRUE)+
  theme_classic()+
  theme(legend.justification=c(),
        legend.position=c()) 
col.isoldistsmall.large.70


#all
col.isoldist.70<- ggplot(col.ibddf.70, aes(x=geodist, y=col.ibd, fill=classes)) +
  #geom_errorbar(aes(ymin=VarPi-VarPiStdErr, ymax=VarPi+VarPiStdErr), colour="black", width=.1, position=pd) +
  #geom_point()+
  geom_point(aes(shape = classes, col=classes), size=2.5)+ 
  labs(x="log (Geographic distance)", y ="Fst/1-Fst")+
  scale_color_brewer(palette="Dark2")+
  geom_smooth(aes(col=classes), method=lm, fill="lightgrey", se=TRUE, level=0.95)+
  theme_classic()+
  theme(legend.justification=c(),
        legend.position=c()) 
col.isoldist.70
###
### subset small
colsmall.70 = col.pop.gl.70@other$size == "Small"
table(colsmall.70)
col.pop.gl.small.70<-col.pop.gl.70[colsmall.70,]
col.pop.gl.small.70@other$size
#
col.ibdsmall.70 = gl.ibd(col.pop.gl.small.70)
#
col.ibdsmall.matrix.70 = as.matrix(col.ibdsmall.70$Dgen)
ind = which( upper.tri(col.ibdsmall.matrix.70), arr.ind = TRUE)
col.ibdsmall.df.70 = data.frame(Site1 = dimnames(col.ibdsmall.matrix.70)[[2]][ind[,2]],
                            Site2 = dimnames(col.ibdsmall.matrix.70)[[1]][ind[,1]],
                            col.col.ibdsmall = col.ibdsmall.matrix.70[ ind ] %>% round(digits = 3))

col.smalldgeo.matrix.70 = as.matrix(col.ibdsmall.70$Dgeo)
ind = which( upper.tri(col.smalldgeo.matrix), arr.ind = TRUE)
col.smalldgeo.df.70 = data.frame(Site1 = dimnames(col.smalldgeo.matrix.70)[[2]][ind[,2]],
                             Site2 = dimnames(col.smalldgeo.matrix.70)[[1]][ind[,1]],
                             col.smalldgeo = col.smalldgeo.matrix.70[ ind ] %>% round(digits = 3))
col.smalldgeo.df.70
col.ibdsmall.df.70$geodist = col.smalldgeo.df.70$col.smalldgeo

col.ibdsmall.df.70$sizeplot = as.factor(c("1 ha vs 10 ha", "1 ha vs 1 ha", "1 ha vs 10 ha", "1 ha vs 10 ha", "10 ha vs 10 ha", "1 ha vs 10 ha", "1 ha vs 1 ha", "1 ha vs 10 ha", "1 ha vs 1 ha", "1 ha vs 10 ha", "1 ha vs 10 ha", "10 ha vs 10 ha", "1 ha vs 10 ha", "10 ha vs 10 ha", "1 ha vs 10 ha"))
str(ibdsmall.df.70)
col.isoldistsmall.70<- ggplot(col.ibdsmall.df.70, aes(x=geodist, y=col.col.ibdsmall)) +
  #geom_errorbar(aes(ymin=VarPi-VarPiStdErr, ymax=VarPi+VarPiStdErr), colour="black", width=.1, position=pd) +
  geom_point(size=3)+ 
  labs(x="log (Geographic distance)", y ="Fst/1-Fst")+
  geom_smooth(method=lm , color="Black", fill="lightgrey", se=TRUE)+
  theme_classic()+
  theme(legend.justification=c(),
        legend.position=c()) 
col.isoldistsmall.70



###subset large
collarge.70 = col.pop.gl.70@other$size == "Large"
table(collarge.70)
col.pop.gl.large.70<-col.pop.gl.70[collarge.70,]
col.pop.gl.large.70@other$size
col.ibdlarge.70 = gl.ibd(col.pop.gl.large.70)

#
col.ibdlarge.matrix.70 = as.matrix(col.ibdlarge.70$Dgen)
ind = which( upper.tri(col.ibdlarge.matrix.70), arr.ind = TRUE)
col.ibdlarge.df.70 = data.frame(Site1 = dimnames(col.ibdlarge.matrix.70)[[2]][ind[,2]],
                            Site2 = dimnames(col.ibdlarge.matrix.70)[[1]][ind[,1]],
                            collarge.ibd = col.ibdlarge.matrix.70[ ind ] %>% round(digits = 3))

col.dgeolarge.matrix.70 = as.matrix(col.ibdlarge.70$Dgeo)
ind = which( upper.tri(col.dgeolarge.matrix.70), arr.ind = TRUE)
col.ibdgeolarge.df.70 = data.frame(Site1 = dimnames(col.dgeolarge.matrix.70)[[2]][ind[,2]],
                             Site2 = dimnames(col.dgeolarge.matrix.70)[[1]][ind[,1]],
                             collarge.dgeo = col.dgeolarge.matrix.70[ ind ] %>% round(digits = 3))
col.ibdgeolarge.df.70
col.ibdlarge.df.70$geodist = col.ibdgeolarge.df.70$collarge.dgeo


col.isoldistlarge.70<- ggplot(col.ibdlarge.df.70, aes(x=geodist, y=collarge.ibd)) +
  #geom_errorbar(aes(ymin=VarPi-VarPiStdErr, ymax=VarPi+VarPiStdErr), colour="black", width=.1, position=pd) +
  geom_point(size=3)+ labs(x="log (Geographic distance)", y ="Fst/1-Fst")+
  geom_smooth(method=lm , color="Black", fill="lightgrey", se=TRUE)+
  theme_classic()+
  theme(legend.justification=c(),
        legend.position=c()) 
col.isoldistlarge.70

# Separate plot list
isoldistlarge
isoldistsmall
isoldistsmall.large

isoldistlarge.70
isoldistsmall.70
isoldistsmall.large.70

col.isoldistlarge
isoldistsmall
isoldistsmall.large

col.isoldistlarge.70
col.isoldistsmall.70
col.isoldistsmall.large.70

####puting together IBD
figureibd <- ggarrange(isoldistsmall, isoldistlarge, isoldistsmall.large,
                       isoldistsmall.70, isoldistlarge.70,isoldistsmall.large.70,
                       col.isoldistsmall, col.isoldistlarge, col.isoldistsmall.large,
                       col.isoldistsmall.70, col.isoldistlarge.70,col.isoldistsmall.large.70,
                       labels = c("A", "B", "C", "D", "E", "F",
                                 "G", "H", "I", "J", "K", "L"), vjust = 0.5,
                       ncol = 3, nrow = 4)
figureibd
