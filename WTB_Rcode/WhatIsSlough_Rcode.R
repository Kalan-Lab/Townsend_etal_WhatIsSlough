# What is Slough - R code for figures 

# needed packages
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggforce)
library(RColorBrewer)
library(pheatmap) # for heatmaps 
library(EnhancedVolcano) # for volcanoplots 
library(limma)
library(qvalue)
library(DEqMS) 
library(matrixStats)
library(phyloseq)
library(microbiome)
library(qiime2R)
library(xlsx)
library(decontam)
library(vegan)
library(DESeq2)
library(ANCOMBC) 
library(mixOmics) 
library(circlize)
library(factoextra) # for k-means cluistering 
library(cluster) # for kmeans clustering
library(reshape2)

# working directory
setwd(dir = "/Users/liztown/Documents/KalanLab/Papers/WTB/WTB_Rcode")

# FIGURE 1 
### DEfining the top enriched proteins in Slough 
# only proteins with an FDR q-value < 1E-4 are included 
GOproteomicsEnrichedTables <- read.csv("./GOproteomicsEnrichedTables.csv")
GOpro <- as.data.frame(GOproteomicsEnrichedTables)

GoProcess <- subset.data.frame(GOpro, GOpro$ProcessFunctionComponent == "Process") 
TopGoProc <-subset.data.frame(GoProcess, GoProcess$Log_FDRqValue <= -4) 
ProcesLST <- sort((unique(TopGoProc$Classification)))
GoFun <- subset.data.frame(GOpro, GOpro$ProcessFunctionComponent == "Function") 
FunLst <-sort(unique(GoFun$Classification))
GoCom <- subset.data.frame(GOpro, GOpro$ProcessFunctionComponent == "Component") 
ComLst <-sort(unique(GoCom$Classification))

GOpro2 <- rbind(GoFun, GoCom, TopGoProc)
GOpro2$Classification <- factor(GOpro2$Classification, levels = c("Blood or Clot","Collagen","Cytoskeleton","Extracellular Lipoprotein",
                                                                  "Extracellular Matrix or Space","Intracellular Organelle Component ", "Membrane",
                                                                  "Platelet","Skin Specific ", # component N=9
                                                                  "Carbohydrate Derivate Binding","Enzyme Regulator Activity","Extracellular Matrix Structure",
                                                                  "Heparin Binding", "Ion Binding","Structural Molecule Activity","Sulfur Compound Binding", # Function = 7
                                                                  "Biosynthetic Process","Blood Clotting","Cell Activation","Cell Death",
                                                                  "Cell Localization","Cell Motility","Complement Cascade","Extracellular Matrix", 
                                                                  "Homeostasis","Immune Function","Inter-cellular Transport","Metabolic Process",
                                                                  "Multicellular Organismal Process","Response to Bacteria","Response to Stimulus","Skin Barrier Formation",
                                                                  "Wound Healing")) # Process = 17) )



PCFcolors <- c(colorRampPalette(c("#412535", "#913C44","#C8443F","#DF7350","#EE925B"))(9), # reds
               colorRampPalette(c("#F5A261","#F6D27A"))(7), # Yellows
               colorRampPalette(c("#E9CE7B","#8AB17D", "#2A9D8F","#053061", "#31273D"))(17)) # BLues


ggplot(GOpro2, aes(x=Log_FDRqValue, y=reorder(Description, -Log_FDRqValue), color = Classification, fill = Classification))+
  geom_point(size = 2.75)+ 
  geom_col(width = .02)+
  scale_color_manual(values = PCFcolors)+
  scale_fill_manual(values = PCFcolors)+
  facet_col(~ProcessFunctionComponent, scales = "free_y", space = "free")+
  theme_classic() +
  scale_y_discrete(position = "right")+
  theme(strip.text = element_blank(), axis.text=element_text(size=5), plot.title = element_text(hjust = 0.5)) #legend.position = "none")


# Figure 2A: Protein HEATMAP 
# # upload ProteinMatrix.csv
ProteinMatrix<-read.csv("./ProteinMatrix.csv")
ProteinMatrix <- as.data.frame(ProteinMatrix)
PM <- ProteinMatrix[,-1]
PM1 <- subset.data.frame(PM, select = c(WTB.001.AVE, WTB.002, WTB.003, WTB.004, WTB.005, WTB.006, WTB.007.AVE,WTB.008, WTB.009, WTB.010.AVE))
Proteins <- as.list(ProteinMatrix$Proteins)

PMatrix <- data.matrix(PM1)
rownames(PMatrix) <- Proteins
WTBmatrix <- PMatrix

WTBcolMetadataAVE <- read.csv("./WTBcolMetadataAVE.csv")
cant<- data.frame(WTBcolMetadataAVE)
colnames(cant) <- c("WTB",	"Subject",	"BodySite","Etiology",	"WndAge(yr)",	"Outcome")
#cant<-cant[!duplicated(cant), ]
WTBn <- cant$WTB
rownames(cant) <- WTBn
head(cant)
cant_row <- cant[-1]
head(cant_row)

pheatmap(WTBmatrix, scale="row",  annotation_col = cant_row, color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100), fontsize_row = 1 )



# Figure 2 B-D VOLCANO PLOTS 
# https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
# USE DEqMS package - for the differential proptein expression analysis 
WTBsimpleProteinMatrix <- read.csv("./WTBsimpleProteinMatrix.csv")
WTBsimpleProteinM <-as.data.frame(WTBsimpleProteinMatrix)# a paired down csv with the Assession numbers in col 1 and the samples in the remainter of the columns (no other metadata)
WTBsimpleProteinM[WTBsimpleProteinM == 0] <- 1   
dat.prot <- as.matrix(WTBsimpleProteinM[,-1]) 
rownames(dat.prot) = WTBsimpleProteinM$Accession 

dat.log = log2(dat.prot) # log transform data
dat.log = na.omit(dat.log) #remove rows with NAs
dat.log = equalMedianNormalization(dat.log) # median centring the data 
dat.log
boxplot(dat.log,las=2)

cond = as.factor(c("Superficial","Deep","Worse","Ongoing","Ongoing",
                   "Healed","Healed","Healed","Superficial","Med", "Deep", "Worse",
                   "Ongoing", "Ongoing","Superficial", "Deep", "Worse", 
                   "PooledControl", "CellLineControl"))
design = model.matrix(~0+cond) # generating the design matrix. 0 means no intercept for the linear model
colnames(design) = gsub("cond","",colnames(design))
x <- c("Worse-Healed","Ongoing-Healed","Worse-Ongoing",
       "Superficial-Deep","Healed-CellLineControl",
       "Worse-CellLineControl","PooledControl-CellLineControl")
contrast =  makeContrasts(contrasts=x,levels=design)
fit1 <- lmFit(dat.log, design)
fit2 <- contrasts.fit(fit1,contrasts = contrast)
fit3 <- eBayes(fit2)

psm.count.table = data.frame(count = rowMins(dat.prot), row.names =  WTBsimpleProteinM$Accession)
fit3$count = psm.count.table[rownames(fit3$coefficients),"count"]
fit4 = spectraCounteBayes(fit3)
head(fit4)

WorseHealed.results = outputResult(fit4,coef_col = 1) 
OngoingHealed.results = outputResult(fit4,coef_col = 2) 
WorseOngoing.results = outputResult(fit4,coef_col = 3) 

geneDescriptionLinker <- read.csv("./geneDescriptionLinker.csv")
linker <- as.data.frame(geneDescriptionLinker) # the CSV file with the gene (assession numbers) and the associated descriptors and functions
WorseHealed <- merge(WorseHealed.results, linker, by = 'gene')
OngoingHealed <- merge (OngoingHealed.results, linker, by = 'gene')
WorseOngoing <- merge (WorseOngoing.results, linker, by = 'gene')

EnhancedVolcano(WorseHealed, x = "logFC", y = "P.Value", lab = WorseHealed$Description,
                pCutoff = 1e-2, FCcutoff = 1 ,
                ylim = c(0,5),
                xlab= "log2 Fold change - Worse / Healed",
                pointSize = 3.0, labSize = 2,
                #drawConnectors = TRUE, widthConnectors = 0.1,
                legendPosition = 'right',legendLabSize = 10, legendIconSize = 3.0,
                title = 'Worse / Healed')
EnhancedVolcano(OngoingHealed, x = "logFC", y = "P.Value", lab = OngoingHealed$Description,
                pCutoff = 1e-2, FCcutoff = 1 ,
                ylim = c(0,5),
                xlab= "log2 Fold change - Ongoing / Healed",
                pointSize = 3.0, labSize = 2,
                #drawConnectors = TRUE, widthConnectors = 0.1,
                legendPosition = 'right',legendLabSize = 10, legendIconSize = 3.0,
                title = 'Ongoing / Healed')
EnhancedVolcano(WorseOngoing, x = "logFC", y = "P.Value", lab = WorseOngoing$Description,
                pCutoff = 1e-2, FCcutoff = 1 ,
                ylim = c(0,5),
                xlab= "log2 Fold change - Worse / Ongoing",
                pointSize = 3.0, labSize = 2,
                #drawConnectors = TRUE, widthConnectors = 0.1,
                legendPosition = 'right',legendLabSize = 10, legendIconSize = 3.0,
                title = 'Worse / Ongoing')

# Supplemental Figure 3 Volcano plots separating by wound age > 1 year, < 1 year 
WTBsimpleProteinM <-as.data.frame(WTBsimpleProteinMatrix)# a paired down csv with the Assession numbers in col 1 and the samples in the remainter of the columns (no other metadata)
WTBsimpleProteinM[WTBsimpleProteinM == 0] <- 1   # making 0 NA
dat.prot <- as.matrix(WTBsimpleProteinM[,-1]) 
rownames(dat.prot) = WTBsimpleProteinM$Accession # adding the row names to the matrix

dat.log = log2(dat.prot) # log transform data
dat.log = na.omit(dat.log) #remove rows with NAs
dat.log = equalMedianNormalization(dat.log) # median centring the data 
dat.log
cond = as.factor(c("Superficial","Deep","Old","Old","Old",
                   "Young","Young","Young","Superficial","Med", "Deep", "Old",
                   "Old", "Young","Superficial", "Deep", "Old", 
                   "PooledControl", "CellLineControl"))
design = model.matrix(~0+cond) # generating the design matrix. 0 means no intercept for the linear model
colnames(design) = gsub("cond","",colnames(design))
x <- c("Old-Young")
contrast =  makeContrasts(contrasts=x,levels=design)
fit1 <- lmFit(dat.log, design)
fit2 <- contrasts.fit(fit1,contrasts = contrast)
fit3 <- eBayes(fit2)
psm.count.table = data.frame(count = rowMins(dat.prot), row.names =  WTBsimpleProteinM$Accession)
fit3$count = psm.count.table[rownames(fit3$coefficients),"count"]
fit4 = spectraCounteBayes(fit3)
head(fit4)

OldYoung.results = outputResult(fit4,coef_col = 1) 

geneDescriptionLinker <- read.csv("./geneDescriptionLinker.csv")
linker <- as.data.frame(geneDescriptionLinker) # the CSV file with the gene (assession numbers) and the associated descriptors and functions

OldYoung <- merge(OldYoung.results, linker, by = 'gene')

EnhancedVolcano(OldYoung, x = "logFC", y = "P.Value", lab = OldYoung$Description,
                pCutoff = 1e-2, FCcutoff = 1 ,
                ylim = c(0,5),
                xlab= "log2 Fold change - Old wounds / Young Wounds",
                pointSize = 3.0, labSize = 2,
                #drawConnectors = TRUE, widthConnectors = 0.1,
                legendPosition = 'right',legendLabSize = 10, legendIconSize = 3.0,
                title = 'Old wounds > 1yr / Young Wounds < 1 yr')


## Figure 3A 
qPCRvCFU <- read.csv("./qPCRvCFU.csv")
PC <-as.data.frame(qPCRvCFU)
ggplot(data = PC, aes(x=log10(CFU), y = log10(NumberBacteriaPerSample)))+
  geom_point(size = 6)+
  stat_smooth(method = lm, level = 0.90, color = "black")+
  #geom_hline(yintercept=c(log10(300)), size = 1, linetype="dotted") + 
  geom_point(size = 7, aes(color = Subject))+
  theme_bw()+
  scale_y_continuous(breaks = c(2,4,6,8))+
  scale_x_continuous(breaks = c(2,4,6,8))+
  ggtitle("Number of Bacteria as determined by CFU vs. qPCR") +
  ylab("Log10(Bacteria) by qPCR")+
  xlab("Log10(Bacteria) by CFU")+
  scale_color_manual(values = c("#053061", "#175169","#21777C", "#2A9D8F", "#8AB17D","#F5D279","#EFB366","#F4A261","#DE7350","#C8443F" ))
 
# Figure 3B 
ps.WTB <- qza_to_phyloseq(features = "./table-WTB-dada2.qza",
                          taxonomy = "./taxonomy-WTB.qza",
                          metadata = "./WTB16Smanufest.txt",
                          tree = "./WTB-rooted-tree.qza")

ps.Decontam <- subset_taxa(ps.WTB, Family !="mitochondria")
ps.Decontam <- subset_taxa(ps.Decontam, Class !="Chloroplast")
ps.Decontam <- subset_taxa(ps.Decontam, Phylum !="Cyanobacteria")
ps.Decontam <- subset_taxa(ps.Decontam, Phylum !=" ")
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Methylobacterium") # Highly previlent contaminants
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Sphingomonas") # Highly previlent contaminants
ps.DecontamWTB <- subset_samples(ps.Decontam, Subject != "NegEX") # removing the negative control
ps.DecontamWTB# THE CLEANED FILE 
sample_data(ps.DecontamWTB)$Etiology <- data.frame(sample_data(ps.DecontamWTB)) %>% mutate(Etiology)
relAll <- transform(ps.DecontamWTB, "compositional")
relAllout <- psmelt(relAll)

Other <- as.data.frame(WTBRelativeOUT)
Other <- Other %>% mutate(Genus = ifelse(Abundance < 0.01 , "Other", Genus)) # replacing the genus catagory with "Other" if the abundance is less than 1% 
Other <- Other %>% mutate(Phylum = ifelse(Abundance < 0.01 , "Other", Phylum)) 
Other$Outcome <- factor(Other$Outcome, levels = c("Healed", "Ongoing", "Worse"))
Other$EdgeCenter<- factor(Other$EddeCenter,  levels = c("Center", "Edge"))
Other$ECoutcome <- paste(Other$EdgeCenter, Other$Outcome) #making column with both outcome and EC variables for easy faceting 
Other$ECbody <- paste(Other$EdgeCenter, Other$BodySite)
Other$Phylum <- factor(Other$Phylum, levels = c("Other", "Actinobacteria", "Firmicutes", "Bacteroidetes", "Proteobacteria") )
Other$Genus<- factor(Other$Genus, levels = c("Other",
                                             "Actinobaculum","Actinomyces", "Brachybacterium", "Corynebacterium", "Dermabacter","Micrococcus","Varibaculum", #Actinobacteria / actinomycetota (7)
                                             "Alloiococcus","Anaerobacillus","Anaerococcus","Dialister","Enterococcus", "Facklamia", "Finegoldia","Helcococcus","Parvimonas","Peptococcus","Peptoniphilus", "Peptostreptococcus","Staphylococcus", "Streptococcus", # Firmicutes/ Bacillota (14)
                                             "Bacteroides", "Porphyromonas", "Prevotella", #Bacteroidetes / Bacteroidota
                                             "Campylobacter", # Campylobacterota 
                                             "Acinetobacter", "Eikenella","Enhydrobacter","Neisseria","Pseudomonas", # Proteobacteria /pseudomonadota (5)
                                             "Desulfovibrio")) # Thermodesulfobacteriota
BigMicro2 <-c("#C1C1C1", #Gray
              colorRampPalette(c("#053061","#2A9D8F"))(7),
              colorRampPalette(c("#8AB17D", "#F6D27A","#F5A261"))(14),
              colorRampPalette(c("#DF7350","#C8443F","#913C44"))(10))

ggplot(Other, aes(x= Subject, y =Abundance, fill = Genus))+
  geom_bar(stat="identity", position="fill") +
  facet_wrap("EdgeCenter", nrow = 2)+
  ylab("Relative Abundance")+
  theme_light()+
  scale_fill_manual(values = BigMicro2)+
  theme(axis.text.x = element_text(size = 7)) +
  ggtitle("Genus present > 0.1% of reads in a samplee")

# Figure 3 C-D 
# Running a PCA on all the MICROBIOME samples individually. (not averaged by subject)
# dataset used here is filtered so that OTU with less than 1% are re-classified as "Other"
MicroAllSamples.df <- read.csv("./WTBRelativeAbundanceWithOTHERotuPIVOT.csv")
MicroAllSamples.df <- as.data.frame(MicroAllSamples.df) 
rownames(MicroAllSamples.df) <- MicroAllSamples.df$Sample
MicroAllSamples.df<-MicroAllSamples.df[,-1]

All.micro.clin.df<- read.csv("./WTBMicroColMetadata.csv")
All.micro.clin.df <- as.data.frame(All.micro.clin.df)
Ogenus <-OTUgenus

Micro.pca <- pca(MicroAllSamples.df) # running the PCA
plotIndiv(Micro.pca, 
          group = All.micro.clin.df$Etiology, # group assigns a color
          pch = as.factor(All.micro.clin.df$BodySite), # Pch assigns shape
          cex = 6,
          title = 'WTB: Microbiome, PCA comp 1 - 2',
          legend.title = 'Etiology', legend.title.pch = 'BodySite',
          legend = TRUE)
coordinates <- plotVar(Micro.pca, plot = FALSE) # to save the variables

PCAcord<-Micro.pca$variates[["X"]]
names <- rownames(PCAcord)
PCAcord1<- cbind(names,PCAcord)
colnames(PCAcord1)[colnames(PCAcord1) == "names"] = "Sample"
PCAcord2<-merge(PCAcord1, All.micro.clin.df, by = "Sample")
PCAcord2$PC1 <- as.numeric(as.character(PCAcord2$PC1))
PCAcord2$PC2 <- as.numeric(as.character(PCAcord2$PC2))

PCAcord2$BodySite <- factor(PCAcord2$BodySite, levels = c("Sacrum & coccyx", "Lower Leg", "Shin", "Ankle"))
ggplot(PCAcord2, aes(x=PC1, y = PC2, color = as.factor(Etiology), fill= as.factor(Etiology), shape = as.factor(BodySite)))+
  geom_jitter(size = 6) +
  scale_shape_manual(values = c(21,22,23,24))+
  stat_ellipse(type = "t", level = 0.75) +
  xlab("PC1: 50% expl. var")+
  ylab("PC2: 29% expl. var")+
  theme_bw()+
  scale_color_manual(values = c( "#053061","#166176","#2C9E8F", "#89B17D"))+ # Blue to Green 
  scale_fill_manual(values = c( "#053061","#166176","#2C9E8F", "#89B17D"))+ # Blue to Green 
  ggtitle("WTB: Microbiome PCA plot")

ggplot(PCAcord2, aes(x=PC1, y = PC2, color = as.factor(BodySite), fill= as.factor(BodySite), shape = as.factor(Etiology)))+
  geom_jitter(size = 6) +
  scale_shape_manual(values = c(21,22,23,24))+
  stat_ellipse(type = "t", level = 0.75) +
  xlab("PC1: 50% expl. var")+
  ylab("PC2: 29% expl. var")+
  theme_bw()+
  scale_color_manual(values = c("#F6D27A", "#F4A261","#DF7350", "#BB4241"))+ # yellows to reds
  scale_fill_manual(values = c("#F6D27A", "#F4A261","#DF7350", "#BB4241"))+ 
  ggtitle("WTB: Microbiome PCA plot")

# Figure 3E 
plotVar(Micro.pca, cutoff = 0.5, cex = 5 )
Driv_OTU <- c("X836ac518d86fdb0dd587761ae39b8845", "X6a418787996565e7641dbbf39b7d3e18", "X18af7b7f2b61429936fcd63a453fcefd", "X5f92ec0b039c9127715e033a44803b16")
Driving_OTU <- as.data.frame(subset.data.frame(coordinates, names %in% Driv_OTU))
colnames(Driving_OTU)[colnames(Driving_OTU) == "x"] = "X1"
colnames(Driving_OTU)[colnames(Driving_OTU) == "y"] = "Y1"

ggplot(Driving_OTU, aes(x =0, y = 0, xend =X1, yend =Y1, label = names))+
  geom_segment()+
  #geom_point(aes(x=X1, y= Y1))+
  geom_text(aes(x=X1, y = Y1, label=as.character(names)))+
  theme_bw()+
  ggtitle("WTB: Microbiome PCA plot DRIVING OTU")


# Supplemental table 6 PERMENOVA on the Micro.PCA data above. 
sampledf <- data.frame(All.micro.clin.df) 
PCAcord_adonis <- as.data.frame(PCAcord1[,-1])
PCAcord_adonis$PC1 <- as.numeric(as.character(PCAcord_adonis$PC1))
PCAcord_adonis$PC2 <- as.numeric(as.character(PCAcord_adonis$PC2))

adonis2(PCAcord_adonis ~ BodySite, data = sampledf, method = "eu", permutations = 9999 )
adonis2(PCAcord_adonis ~ Etiology, data = sampledf, method = "eu", permutations = 9999 )
adonis2(PCAcord_adonis ~ Outcome, data = sampledf, method = "eu", permutations = 9999 )
adonis2(PCAcord_adonis ~ Subject, data = sampledf, method = "eu", permutations = 9999 )
adonis2(PCAcord_adonis ~ EdgeCenter, data = sampledf, method = "eu", permutations = 9999 )


# Figure 4
# PREPING Datasets for DIABLO 
#   Kmeans clustering --- https://www.statology.org/k-means-clustering-in-r/#:~:text=K%2Dmeans%20clustering%20is%20a,quite%20different%20from%20each%20other.
WTBsimpleProteinM <- read.csv("./WTBRelativeAbundanceWithOTHERotuPIVOT.csv")
WTBsimpleProteinM <-as.data.frame(WTBsimpleAVEProteinMatrix)# a paired down csv with the Assession numbers in col 1 and the samples in the remainter of the columns (no other metadata)
WTBsimpleProteinM[WTBsimpleProteinM == 0] <- 1 # make 0 into 1 (as a pseudo count)
dat.prot <- as.matrix(WTBsimpleProteinM[,-1]) 
WTB.protein.scaled <- scale(dat.prot) # scaling each of the rows to have a center of 0 and std of 1 (fvizclust does not work on log transformed data)

gap_stat <- clusGap(WTB.protein.scaled,FUN = kmeans, nstart = 25, K.max = 50, B = 50) # calculating the gap stat
fviz_gap_stat(gap_stat) #ignore the warnings ^ above. Optimal # clusters looks like its 23

set.seed(1)
WTB.prot.km <- kmeans(WTB.protein.scaled, centers = 23, nstart = 25)
WTB.prot.km
fviz_cluster(WTB.prot.km , data = WTB.protein.scaled, labelsize = 1)
WTB.prot.Kclusters <- aggregate(WTB.protein.scaled, by=list(cluster=WTB.prot.km$cluster), mean)
final_WTB_prot_clustered_data <- cbind(WTB.protein.scaled, cluster =WTB.prot.km$cluster)

WTBproteins_annotated<- read.csv("./WTBproteins_annotated.csv")
anot <- as.data.frame(WTBproteins_annotated)
WTB.pro.scaled.A <-tibble::rownames_to_column(as.data.frame(WTB.protein.scaled), "Accession")
WTB.pro.scaled.Annotated <- merge(anot,WTB.pro.scaled.A, by = "Accession", all = F )

# Figure 5: Pls-da 
WTBwndAssessment_numbers<- read.csv("./WTBwndAssessment_numbers.csv")
WndAsmt.df<- as.data.frame(WTBwndAssessment_numbers)
rownames(WndAsmt.df) <- WndAsmt.df$Subject
WndAsmt.df <- WndAsmt.df[,-1]
Wnd_List <- colnames(WndAsmt.df)

WTBRelativeAbundanceOTUAveragedBySubjectFILTERED<- read.csv("./WTBRelativeAbundanceWithOTHERaveragedBySubject.csv")
MicroAveOTUSamples.df <- as.data.frame(WTBRelativeAbundanceOTUAveragedBySubjectFILTERED) # the microbiome dataset from the PCA above, with <1% OTU incorporated to other. and all Edge and Center wound samples were averaged by subject
rownames(MicroAveOTUSamples.df) <- MicroAveOTUSamples.df$Sample
MicroAveOTUSamples.df<-MicroAveOTUSamples.df[,-1]
Micro_list<-colnames(MicroAveOTUSamples.df)

WTBproteinsKclusterMeans<-read.csv("./WTBproteinsKclusterMeans.csv")
WTB.pClust.df <- as.data.frame(WTBproteinsKclusterMeans)
rownames(WTB.pClust.df) <-WTB.pClust.df$Subject
WTB.pClust.df<-WTB.pClust.df[-1]

X <- list(Proteins = WTB.pClust.df,
          Microbes = MicroAveOTUSamples.df,
          Wound = WndAsmt.df) 
Y <- as.factor(WTB.clin.df$Outcome) 

WTBResult.diablo <- block.splsda(X, Y)
coordinates <- plotVar(WTBResult.diablo, plot = FALSE) # to save the variables
plotVar(WTBResult.diablo, var.names = c(FALSE,F,F),
        cutoff = 0.5,
        legend=TRUE, pch=c(16,16,16),
        title = "Diablo variable plot - cutoff 0.5")

plotIndiv(WTBResult.diablo,
          ind.names = F, 
          legend=TRUE, cex=6,
          col =c("#8AB17D","#F6D27A","#F5A261"),
          xlim = c(-5,5), ylim = c(-5,5),
          title = 'WTB with DIABLO')

plotArrow(WTBResult.diablo, 
          ind.names = F,
          group = WTB.clin.df$Outcome,
          legend=TRUE, #cex=c(5,5,5,5),
          pch.size = 6, arrow.size = 1.5,
          col.per.group =c("#8AB17D","#F6D27A","#F5A261"),
          #xlim = c(-4,4), ylim = c(-5,5),
          title = 'WTB with DIABLO')
plotArrow(WTBResult.diablo, 
          ind.names = T,
          ind.names.size = 4,
          group = WTB.clin.df$Outcome,
          legend=TRUE, #cex=c(5,5,5,5),
          pch.size = 5.5, arrow.size = 1.5,
          col.per.group =c("#8AB17D","#F6D27A","#F5A261"),
          #xlim = c(-4,4), ylim = c(-5,5),
          title = 'WTB with DIABLO')

plotLoadings(WTBResult.diablo, comp = 1, contrib = 'max', method = 'median')
plotLoadings(WTBResult.diablo, comp = 2, contrib = 'max', method = 'median')

cord1_Load_Pro<-plotLoadings(WTBResult.diablo, comp = 1, block = "Proteins", contrib = 'max', method = 'median')
cord1_Load_Micro<-plotLoadings(WTBResult.diablo, comp = 1, block = "Microbes", contrib = 'max', method = 'median')
cord1_Load_Wound<-plotLoadings(WTBResult.diablo, comp = 1, block = "Wound", contrib = 'max', method = 'median')
Cord_1_Loadings_dataframe <- rbind(cord1_Load_Pro, cord1_Load_Micro,cord1_Load_Wound)

cord2_Load_Pro<-plotLoadings(WTBResult.diablo, comp = 2, block = "Proteins", contrib = 'max', method = 'median')
cord2_Load_Micro<-plotLoadings(WTBResult.diablo, comp = 2, block = "Microbes", contrib = 'max', method = 'median')
cord2_Load_Wound<-plotLoadings(WTBResult.diablo, comp = 2, block = "Wound", contrib = 'max', method = 'median')
Cord_2_Loadings_dataframe <- rbind(cord2_Load_Pro, cord2_Load_Micro, cord2_Load_Wound)

Cord_1_Loadings_dataframe <- tibble::rownames_to_column(Cord_1_Loadings_dataframe, "Variable")
TiesX3 <- filter(Cord_1_Loadings_dataframe, Contrib.Healed == "TRUE" & Contrib.Ongoing == "TRUE" & Contrib.Worse == "TRUE") # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_1_Loadings_dataframe_filt<-Cord_1_Loadings_dataframe[ !Cord_1_Loadings_dataframe$Variable %in% Tie3Var,]

Healed_1 <- filter(Cord_1_Loadings_dataframe_filt, Contrib.Healed == "TRUE")
Healed_1$GroupContrib[Healed_1$GroupContrib == "tie"] <- "Healed"
Ongoing_1 <- filter(Cord_1_Loadings_dataframe_filt, Contrib.Ongoing == "TRUE")
Ongoing_1$GroupContrib[Ongoing_1$GroupContrib == "tie"] <- "Ongoing"
Worse_1 <- filter(Cord_1_Loadings_dataframe_filt, Contrib.Worse == "TRUE")
Worse_1$GroupContrib[Worse_1$GroupContrib == "tie"] <- "Worse"

Cord_1_Filtered_Loadings.DF <- rbind(Healed_1,Ongoing_1,Worse_1)
Cord_1_Filtered_Loadings.DF <- Cord_1_Filtered_Loadings.DF %>% mutate(ProMicWou= case_when(Variable %in% Protein_list ~ 'Protein Cluster', Variable %in% Micro_list ~ 'Micro OTU', Variable %in% Wnd_List ~ 'Wound Assessment'))

ggplot(Cord_1_Filtered_Loadings.DF, aes(x = importance, y = reorder(Variable, abs(importance)), color = ProMicWou, fill = ProMicWou))+
  geom_point(aes(size = 1.5*abs(importance)))+
  geom_col(width = .04)+
  xlim(-.5,0.5)+
  scale_color_manual(values = c("#F0975D", "#0F4C6D", "#C8443F"))+
  scale_fill_manual(values = c("#F0975D", "#0F4C6D", "#C8443F"))+
  facet_wrap(~factor(GroupContrib, levels = c("Worse", "Ongoing", "Healed")), scales = "free_y")+
  theme_bw()

# Figure 6 and Supplemental Figure 
# go terms identified through the GO panther software. due to different numbers and variety of proteins in each cluster the significant of the biologic process enrichment varied. 
# so for each cluster GO biologic process enrichment was ranked by the LOWEST corrected P-value (FDR), followed by lowest uncorrected P-value if there were < 25 significant terms. 
# (1 = most enriched biologic process, 25 = least enriched of the top 25)
TopClusterGoTerms <- read.csv("./TopClusterGoTerms.csv")
TopGoTerms <- as.data.frame(TopClusterGoTerms)
TopGoTerms$Cluster <- factor(TopGoTerms$Cluster, levels = c("Cluster-1", "Cluster-2", "Cluster-3","Cluster-4","Cluster-5","Cluster-6","Cluster-7","Cluster-8","Cluster-9","Cluster-10","Cluster-11","Cluster-12","Cluster-13","Cluster-14","Cluster-15","Cluster-16","Cluster-17","Cluster-18","Cluster-19","Cluster-20","Cluster-21","Cluster-22", "Cluster-23"))
table(TopGoTerms["Classification"]) # 25 variables 

PCFcolors <- c(colorRampPalette(c("#412535", "#913C44","#C8443F","#DF7350","#EE925B"))(10), # reds
               colorRampPalette(c("#F5A261","#F6D27A"))(5), # Yellows
               colorRampPalette(c("#E9CE7B","#8AB17D", "#2A9D8F","#053061", "#31273D"))(10)) # BLues
ggplot(TopGoTerms, aes(x= Cluster, y = reorder(GO.biological.process,-ArbitraryNumberRank), color = Classification, fill = Classification)) +
  geom_point(size = 2.5)+
  facet_wrap(factor(Cluster) ~., scales = "free", ncol = 6) + # rows = 6, col = 4)+
  scale_color_manual(values = PCFcolors)+
  scale_fill_manual(values = PCFcolors)+
  theme_bw()+
  scale_x_discrete(expand = c(0.001,0))+
  theme(strip.background = element_blank(),strip.text.x = element_blank(),axis.text.y=element_text(size=2), axis.text.x=element_text(size=5))

KeyClu <- c("Cluster-6", "Cluster-21", "Cluster-11","Cluster-7", "Cluster-1", "Cluster-12", "Cluster-22", "Cluster-19", "Cluster-5" )
KeyClusters<-filter(TopGoTerms, Cluster %in% KeyClu)
PCFcolors2 <- c(colorRampPalette(c("#412535", "#913C44","#C8443F","#DF7350","#EE925B"))(7), # reds
                colorRampPalette(c("#F5A261","#F6D27A"))(3), # Yellows
                colorRampPalette(c("#E9CE7B","#8AB17D", "#2A9D8F","#053061", "#31273D"))(9)) # BLues
ggplot(KeyClusters, aes(x= Cluster, y = reorder(`GO biological process`,-ArbitraryNumberRank), color = Classification, fill = Classification)) +
  geom_point(size = 2.5)+
  facet_wrap(factor(Cluster) ~., scales = "free", ncol = 3) + # rows = 6, col = 4)+
  scale_color_manual(values = PCFcolors2)+
  scale_fill_manual(values = PCFcolors2)+
  theme_bw()+
  scale_x_discrete(expand = c(0.001,0))+
  theme(strip.background = element_blank(),strip.text.x = element_blank(),axis.text.y=element_text(size=2), axis.text.x=element_text(size=5))





