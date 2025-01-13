
###
# Script for the seminar on Community Phylogenetics
# Cutting vEdge Tools (IAVS EcoInformatics)
# 31/01/2025
###

# Load the necessary packages
library(ape)
library(phytools)
library(picante)
library(ggplot2)
library(phyloregion)
library(reshape2)

###
# Table of contents:
#
# A. Generate a random phylogenetic tree
# B. Generate a random set of sites/communities
# C. Calculate metrics of alpha diversity
# C.1. Faith's PD
# C.2. Mean Pairwise Distance (MPD)
# C.3. Mean Nearest Taxon Distance (MNTD)
# C.4. Community Evolutionary Distinctiviness (CED)
# D. Calculate metrics of beta diversity
# D.1. PhyloSor
# D.2. Dpw/Dnn
#
###

##
# A. Generate a random phylogenetic tree:
set.seed(134) #for reproducibility
tree <- rcoal(10) #10 species

# Plot the tree:
tree$tip.label <- paste0("sp", 1:10) #assign species names
plot(tree)

##
# B. Generate a random set of sites/communities
com<-data.frame(
  sp1=c(1,1,0,1,1,0,1,0,0,0),
  sp2=c(0,0,0,1,1,0,1,1,0,0),
  sp3=c(0,0,0,1,1,0,1,0,0,0),
  sp4=c(1,1,0,1,0,1,0,1,0,0),
  sp5=c(0,1,0,1,0,1,0,0,0,0),
  sp6=c(0,1,1,0,0,1,1,1,0,0),
  sp7=c(0,1,1,0,0,1,1,1,0,0),
  sp8=c(0,0,1,0,1,0,0,0,1,0),
  sp9=c(0,1,0,0,1,0,1,0,1,1),
  sp10=c(1,0,1,0,1,0,0,0,1,1)
)
rownames(com)<-paste0("site", 1:10) #set sites names
t(com) #transpose

##
# C. Calculate metrics of alpha diversity

##
# C.1. Faith's PD
res <- pd(samp = com,
          tree = tree,
          include.root = FALSE)
res #visualize the result

# Test for correlation with species richness and plot the relationship:
cor.test(res$PD, res$SR)
plot(res$SR, 
     res$PD, 
     xlab = "Species Richness", ylab = "Faith's PD", 
     pch = 16)

# Get ses.PD (i.e, standardized effect size of Faith's PD):
res <- ses.pd(samp = com,
              tree = tree,
              null.model = "independentswap",
              runs = 999)
res #visualize the result

# Test for correlation with species richness and plot the relationship:
cor.test(res$pd.obs.z, res$ntaxa)
plot(res$ntaxa, 
     res$pd.obs.z, 
     xlab = "Species Richness", ylab = "Faith's PD (SES)", 
     pch = 16)

##
# C.2. Mean Pairwise Distance (MPD)
res_mpd <- ses.mpd(samp = com,
               dis = cophenetic(tree),
               null.model = "independentswap",
               abundance.weighted = FALSE, 
               runs = 999)
res_mpd #visualize the result

# Test for correlation with ses.PD and plot the relationship:
cor.test(res$pd.obs.z, res_mpd$mpd.obs.z)
plot(res$pd.obs.z, 
     res_mpd$mpd.obs.z, 
     xlab = "Faith's PD (SES)", ylab = "MPD (SES)", 
     pch = 16)

# C.3. Mean Nearest Taxon Distance (MNTD)
res_mntd <- ses.mntd(samp = com,
                   dis = cophenetic(tree),
                   null.model = "independentswap",
                   abundance.weighted = FALSE, 
                   runs = 999)
res_mntd #visualize the result

# Test for correlation with ses.MPD and plot the relationship:
cor.test(res_mpd$mpd.obs.z, res_mntd$mntd.obs.z)
plot(res_mpd$mpd.obs.z, 
     res_mntd$mntd.obs.z, 
     xlab = "MPD (SES)", ylab = "MNTD (SES)", 
     pch = 16)

# Test for correlation with ses.PD and plot the relationship:
cor.test(res$pd.obs.z, res_mntd$mntd.obs.z)
plot(res$pd.obs.z, 
     res_mntd$mntd.obs.z, 
     xlab = " Faith's PD (SES)", ylab = "MNTD (SES)", 
     pch = 16)

##
# C.4. Community Evolutionary Distinctiveness (CED)
ev<-as.data.frame(evol_distinct(tree, type="fair.proportion"))
ev$species<-rownames(ev) #add a column with species names
names(ev)[1]<-"ED" #change name of column containing ED values

# Convert community data from wide to long format:
long_com<-com #create a copty
long_com$sites<-rownames(long_com) #create a column with rownames
long_com<-melt(long_com, id.vars = c("sites"), variable.name = "species") #convert into long format
long_com<-subset(long_com, value >0) #subset rows with present species

#Merge tables:
out<-merge(long_com, ev, by="species", all.x=T)

#Get mean ED:
ed<-aggregate(out$ED, list(out$sites), mean)
names(ed)<-c("sites", "CED")
ed

##
# D. Calculate metrics of beta diversity

##
# D.1. PhyloSor
t<-round(vegdist(com, method="bray"), 2) #1st, taxonomic (Sorensen), for comparison
p<-1-round(phylosor(com, tree), 2) #2n, phylogenetic (PhyloSor)

# Mantel test between the two matrices:
mantel(t, p, method = "spearman", permutations = 999, na.rm = TRUE)

# Run NMDS for taxonomic beta diversity:
set.seed(123) #for reproducibility
vare.mds <- metaMDS(t) #nmds
data.scores <- as.data.frame(scores(vare.mds))  #using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  #create a column of site names, from the rownames of data.scores

#Plot:
ggplot() + 
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2),size=3) + # add the point markers
  geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site),size=6,vjust=0) +  # add the site labels
  theme_bw()

# Run NMDS for phylogenetic beta diversity:
set.seed(123) #for reproducibility
vare.mds <- metaMDS(p) #nmds
data.scores <- as.data.frame(scores(vare.mds))  #using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  #create a column of site names, from the rownames of data.scores

#Plot:
ggplot() + 
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2),size=3) + # add the point markers
  geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site),size=6,vjust=0) +  # add the site labels
  theme_bw()

##
# D.2. Dpw/Dnn
dpw<-comdist(com, cophenetic(tree))
dnn<-comdistnt(com, cophenetic(tree))

# Test for correlations between matrices:
mantel(p, dpw, method = "spearman", permutations = 999, na.rm = TRUE)
mantel(p, dnn, method = "spearman", permutations = 999, na.rm = TRUE)

# Run NMDS for Dpw:
set.seed(123) #for reproducibility
vare.mds <- metaMDS(dpw) #nmds
data.scores <- as.data.frame(scores(vare.mds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  #create a column of site names, from the rownames of data.scores

#Plot:
ggplot() + 
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2),size=3) + # add the point markers
  geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site),size=6,vjust=0) +  # add the site labels
  theme_bw()

# Run NMDS for Dnn:
set.seed(123) #for reproducibility
vare.mds <- metaMDS(dnn) #nmds
data.scores <- as.data.frame(scores(vare.mds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  #create a column of site names, from the rownames of data.scores

#Plot:
ggplot() + 
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2),size=3) + # add the point markers
  geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site),size=6,vjust=0) +  # add the site labels
  theme_bw()

# Clean-up environment:
rm(list=ls())
