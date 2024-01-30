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
library(rrvgo)
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
library(speedyseq)

# working directory
setwd(dir = "/Users/liztown/Documents/KalanLab/Papers/WTB/WTB_Files_for_Rcode")

# FIGURE 1 (PCoA scatterplots) and supplemental Figure 2 (corresponding Treemap plots)
### visualizing the most significantly  enriched proteins in Slough 
GOproteomicsEnrichedTables <- read.csv("./GOproteomicsEnrichedTables.csv")
GOpro <- as.data.frame(GOproteomicsEnrichedTables)

GoProcess <- subset.data.frame(GOpro, GOpro$ProcessFunctionComponent == "Process") 
GoFun <- subset.data.frame(GOpro, GOpro$ProcessFunctionComponent == "Function") 
GoCom <- subset.data.frame(GOpro, GOpro$ProcessFunctionComponent == "Component") 

BP_simMatrix <- calculateSimMatrix(GoProcess$Goterm,
                                   orgdb="org.Hs.eg.db",
                                   ont="BP",
                                   method="Rel")
BP_scores <- setNames(-log10(GoProcess$FDR.q.value), GoProcess$Goterm)
BP_reducedTerms <- reduceSimMatrix(BP_simMatrix, BP_scores,
                                   threshold=0.7,
                                   orgdb="org.Hs.eg.db")

scatterPlot(BP_simMatrix, BP_reducedTerms)+  # Fig 1A
  scale_color_manual(values = c(colorRampPalette(c("#021327", "#053061","#2A9D8F","#8AB17D","#F6D27A",
                                                   "#F5A261","#DF7350","#C8443F", "#7E2825"))(24)))+
  theme(legend.position = "none")

pdf(file="./BP_tree.pdf") # fig S2A
treemapPlot(BP_reducedTerms)
dev.off()

CC_simMatrix <- calculateSimMatrix(GoCom$Goterm,
                                   orgdb="org.Hs.eg.db",
                                   ont="CC",
                                   method="Rel")
CC_scores <- setNames(-log10(GoCom$FDR.q.value), GoCom$Goterm)
CC_reducedTerms <- reduceSimMatrix(CC_simMatrix, CC_scores,
                                   threshold=0.7,
                                   orgdb="org.Hs.eg.db")

scatterPlot(CC_simMatrix, CC_reducedTerms)+ # Fig 1B
  scale_color_manual(values = c(colorRampPalette(c("#021327", "#053061","#2A9D8F","#8AB17D","#F2BF40",
                                                   "#F5A261","#DF7350","#C8443F", "#7E2825"))(14)))+
  theme(legend.position = "none")

pdf(file="./CC_tree.pdf") # Fig S2B
treemapPlot(CC_reducedTerms) 
dev.off()

Fun_simMatrix <- calculateSimMatrix(GoFun$Goterm,
                                    orgdb="org.Hs.eg.db",
                                    ont="MF",
                                    method="Rel")
MF_scores <- setNames(-log10(GoFun$FDR.q.value), GoFun$Goterm)
MF_reducedTerms <- reduceSimMatrix(Fun_simMatrix, MF_scores,
                                   threshold=0.7,
                                   orgdb="org.Hs.eg.db")

scatterPlot(Fun_simMatrix, MF_reducedTerms)+ # Fig 1C
  scale_color_manual(values = c(colorRampPalette(c( "#053061","#2A9D8F","#8AB17D",
                                                    "#F5A261","#C8443F"))(5)))+
  theme(legend.position = "none")

pdf(file="./MF_tree.pdf") # Fig S2C
treemapPlot(MF_reducedTerms)
dev.off()


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
## additional corresponding plots for Fig 2-B-D with BH adjusted p-values 
EnhancedVolcano(WorseHealed, x = "logFC", y = "adj.P.Val", lab = WorseHealed$Description,
                pCutoff = 0.05, FCcutoff = 1 ,
                ylim = c(0,2),
                xlab= "log2 Fold change - Worse / Healed",
                pointSize = 3.0, labSize = 2,
                #drawConnectors = TRUE, widthConnectors = 0.1,
                legendPosition = 'right',legendLabSize = 10, legendIconSize = 3.0,
                title = 'Worse / Healed')
EnhancedVolcano(OngoingHealed, x = "logFC", y = "adj.P.Val", lab = OngoingHealed$Description,
                pCutoff = 0.05, FCcutoff = 1 ,
                ylim = c(0,2),
                xlab= "log2 Fold change - Ongoing / Healed",
                pointSize = 3.0, labSize = 2,
                #drawConnectors = TRUE, widthConnectors = 0.1,
                legendPosition = 'right',legendLabSize = 10, legendIconSize = 3.0,
                title = 'Ongoing / Healed')
EnhancedVolcano(WorseOngoing, x = "logFC", y = "adj.P.Val", lab = WorseOngoing$Description,
                pCutoff = 0.05, FCcutoff = 1 ,
                ylim = c(0,2),
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

# Figure 3 C-D: Bray Curtis Beta diversity NMDS ordinations 
min(sample_sums(ps.DecontamWTB))# minimum sample read is 7
median(sample_sums(ps.DecontamWTB)) # 40021
max(sample_sums(ps.DecontamWTB)) #77764
table(sample_sums(ps.DecontamWTB))

rarecurve(t(otu_table(ps.DecontamWTB)), step=100) 
rarecurve(t(otu_table(ps.DecontamWTB)), step=100, ylim =c(0,60), xlim=c(0,1500))  ## considering using a read cut off of 5000 for beta diversity metrics 

ps.rarefiedWTB = rarefy_even_depth(ps.DecontamWTB, rngseed=1, sample.size=840, replace=F)

GP = ps.rarefiedWTB
GP.ord <- ordinate(GP, "NMDS",  "bray") # For Bray Curtis of the rareified dataset
plot_ordination(GP, GP.ord, type="samples", color="Etiology2", shape = "BodySite") + 
  geom_point(size=7) +
  scale_shape_manual(values = c(19,17, 15,23))+
  theme_bw()+
  scale_color_manual(values = c( "#053061","#166176","#2C9E8F", "#89B17D"))+ # Blue to Green 
  scale_fill_manual(values = c( "#053061","#166176","#2C9E8F", "#89B17D"))+ # Blue to Green 
  ggtitle("Bray Curtis")

plot_ordination(GP, GP.ord, type="samples", color="Etiology2", shape = "EdgeCenter") + 
  geom_point(size=6) +
  scale_shape_manual(values = c(19,17, 15,23))+
  theme_bw()+
  scale_color_manual(values = c( "#053061","#166176","#2C9E8F", "#89B17D"))+ # Blue to Green 
  scale_fill_manual(values = c( "#053061","#166176","#2C9E8F", "#89B17D"))+ # Blue to Green 
  ggtitle("Bray Curtis- Edge center plot")

plot_ordination(GP, GP.ord, type="samples", color="BodySite", shape = "Etiology2") + 
  geom_point(size=6) + 
  theme_bw()+
  scale_shape_manual(values = c(19,17, 15,23))+
  scale_color_manual(values = c("#F6D27A", "#F4A261","#DF7350", "#BB4241"))+ # yellows to reds
  scale_fill_manual(values = c("#F6D27A", "#F4A261","#DF7350", "#BB4241"))+ 
  ggtitle("Bray Curtis")

## Supplemental Table 6
sampledf <- data.frame(sample_data(ps.rarefiedWTB)) # viable dry Surg samples
WTB_bray <- phyloseq::distance(ps.rarefiedWTB, method = "bray")
Bray_clust <-hclust(WTB_bray)
adonis2(WTB_bray ~ BodySite, by= "margin", data = sampledf, permutations = 9999) 
adonis2(WTB_bray ~ Etiology2, by= "margin", data = sampledf, permutations = 9999)
adonis2(WTB_bray ~ EdgeCenter, by= "margin", data = sampledf, permutations = 9999)
adonis2(WTB_bray ~ Outcome, by= "margin", data = sampledf, permutations = 9999)

WTB.p = c(0.002, 0.001,0.9488, 0.1787) # correcting p-values
p.adjust(WTB.p, method="BH")


# Figure 3E: Vector plot of taxa that corelate with points positions on hte NMDS plot
GP = ps.rarefiedWTB
GP.ord <- ordinate(GP, "NMDS",  "bray")
BrayDF <- as.data.frame(GP.ord$points)
BrayPositions <- subset.data.frame(BrayDF, select = c("MDS1", "MDS2")) %>% as_tibble(rownames = "Sample")

WTB_OUT <-as.data.frame(psmelt(ps.rarefiedWTB)) # this is the CSV with the full metadata and the count of the abundance for each otu in each sample

OTUtaxa <- subset.data.frame(WTB_OUT , select = c("OTU", "Phylum", "Genus"))

Brayshaired<- merge(WTB_OUT , BrayPositions, by = "Sample") # merging the two dataframes to get corelations 
cor_x <- Brayshaired %>%
  nest(data = -OTU) %>%
  mutate(cor_x = purrr::map(data, ~ cor.test(.x$Abundance, .x$MDS1, method = "spearman",exact = FALSE) %>% tidy())) %>%
  dplyr::select(OTU, cor_x) %>%
  unnest(cor_x) %>%
  dplyr::select(OTU, estimate, p.value)

cor_y  <- Brayshaired %>%
  nest(data = -OTU) %>%
  mutate(cor_y = purrr::map(data, ~ cor.test(.x$Abundance, .x$MDS2, method = "spearman",exact = FALSE) %>% tidy())) %>%
  dplyr::select(OTU, cor_y) %>%
  unnest(cor_y) %>%
  dplyr::select(OTU, estimate, p.value)

correlations <- inner_join(cor_x, cor_y, by = "OTU") 

filtered <- correlations %>%
  filter(p.value.x < 0.05 | p.value.y < 0.05)

corTaxa <- merge(filtered, OTUtaxa, by = "OTU") 
corTaxa <- unique(corTaxa)
corTaxa # THIS IS THE Data fram with the list of the significant OTU
write.csv(corTaxa, "./BrayWTB_all_CorTaxa.csv")

OTUtaxaAnnotations <- subset.data.frame(corTaxa, select = c("OTU", "Phylum", "Genus"))
write.csv(OTUtaxaAnnotations, "./Bray_WTB_ALL_CorTaxaAnnotations.csv")

high_corr <- filter(corTaxa, abs(estimate.x >0.2) | abs(estimate.y>0.2))
high_corr%>% 
  ggplot(aes(x=0, xend=estimate.x, y=0, yend =estimate.y)) +
  geom_segment() 

BrayPositions %>%
  ggplot(aes(x=MDS1, y = MDS2)) +
  geom_point()+
  geom_segment(data = high_corr, aes(x=0, xend=estimate.x, y=0, yend =estimate.y))+
  geom_text_repel(data = high_corr, aes(x=estimate.x, y =estimate.y, label = OTU), min.segment.length = 0.01, max.overlaps = 15, size = 2) +
  theme_bw()

BrayPositions %>%
  ggplot(aes(x=MDS1, y = MDS2)) +
  geom_point()+
  geom_segment(data = high_corr, aes(x=0, xend=estimate.x, y=0, yend =estimate.y, color = Genus), size= 1.1)+
  geom_text_repel(data = high_corr, aes(x=estimate.x, y =estimate.y, label = Genus), min.segment.length = 0.01, max.overlaps = 15, size = 3) +
  theme_bw()

BrayPositions %>%
  ggplot(aes(x =NMDS1, y = NMDS2)) +
  #geom_point()+
  geom_segment(data = high_corr, aes(x=0, xend=estimate.x, y=0, yend =estimate.y), size= 1.1)+
  geom_text_repel(data = high_corr, aes(x=estimate.x, y =estimate.y, label = Genus), min.segment.length = 0.01, max.overlaps = 15, size = 3) +
  scale_y_continuous(limits = c(-1,1))+
  scale_x_continuous(limits = c(-1,1))+
  theme_bw()


# Supplemental Figure 4 
ps.wisc <- qza_to_phyloseq(
  features="QZA-Wisc-table.qza",
  taxonomy="QZA-Wisc-taxonomy.qza",
  metadata = "QZA-Wisc-metadata.txt",
  tree = "QZA-Wisc-rooted_tree.qza")
ps.decont <- ps.wisc
ps.decont <- prune_samples(sample_sums(ps.decont) > 1000, ps.decont)
ps.decont <- subset_taxa(ps.decont, is.na(Phylum) == F)
ps.decont <- subset_taxa(ps.decont, Family != 'mitochondria')
tax_table(ps.decont)["829ee0e126800d181ab4a09065174a4d", "Genus"] <- "Klebsiella/Enterobacter" # manually blasted and renamed abundant but "unclassified" ASVs at the genus level
tax_table(ps.decont)["5bfb87b61ec26b3347297c36d287f034", "Genus"] <- "Escherichia"
tax_table(ps.decont)["dd4f1e0c4ea9768f0d6371bdb5eb9ab7", "Genus"] <- "Escherichia"
tax_table(ps.decont)["430c3059353af8d4cf475ae776c5fb31", "Genus"] <- "Klebsiella/Enterobacter"
tax_table(ps.decont)["a8fead32776fba662549a2e627961359", "Genus"] <- "Klebsiella/Enterobacter"
tax_table(ps.decont)["93f7cd0fcb9c489339bcd538865d2b68", "Genus"] <- "Klebsiella/Enterobacter"
tax_table(ps.decont)["1557d8aa883f653f076f0c185bd523a8", "Genus"] <- "Alcaligenes"
tax_table(ps.decont)["01acdfc6b905737920c30d90f0db7852", "Genus"] <- "Alcaligenes"
sample_data(ps.decont)$Patient_ID <- as.character(sample_data(ps.decont)$Patient_ID)
ps.decont <- merge_samples(ps.decont, "Patient_ID")
ps.comp <- transform_sample_counts(ps.decont, function(x) x / sum(x))
df.ps = psmelt(ps.comp)
df.ps.long <- pivot_longer(df.ps, c(Kingdom:Species), names_to = 'level', values_to = 'taxon') %>%
  mutate(Patient_ID = fct_relevel(as.character(Patient_ID),
                                  "4", "6", "9",
                                  "8", "5", "3", "7",
                                  "10", "1", "2"))
tax_Gen_Fam_Phylum <- df.ps %>%
  select(Genus, Family, Phylum) %>%
  unique()
tax_Gen_Fam_Phylum[nrow(tax_Gen_Fam_Phylum)+1,] <- c("Other (<1%)", "Other (<1%)", "Other (<1%)")
tax_Gen_Fam_Phylum[nrow(tax_Gen_Fam_Phylum)+1,] <- c("Unclassified", "Other (<1%)", "Other (<1%)")
taxon_rel_abund <- df.ps.long %>%
  replace_na(list(taxon = "Unclassified", level = "Unclassified", Abundance = 0)) %>%
  filter(level=="Genus") %>%
  #filter(Patient_ID == 1) %>%
  group_by(Patient_ID, taxon) %>%
  summarize(rel_abund = sum(Abundance), .groups="drop") %>%
  group_by(Patient_ID, taxon) %>%
  summarize(mean_rel_abund = mean(rel_abund)*100, .groups="drop")
taxon_pool <- taxon_rel_abund %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 5,
            mean = mean(mean_rel_abund),
            .groups="drop")
phylum_color2 = c("#C8443F","#175169","#8AB17D", "#EFB366")
pWisc<-
  inner_join(taxon_rel_abund, taxon_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(Patient_ID, taxon) %>%
  summarize(mean_rel_abund = mean(mean_rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  inner_join(., tax_Gen_Fam_Phylum, by=c('taxon'='Genus')) %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc=F),
         Phylum = fct_relevel(Phylum, "Proteobacteria", "Actinobacteria", "Firmicutes", "Bacteroidetes", "None")
  ) %>%
  ggplot(aes(x=Patient_ID, y=taxon, size=mean_rel_abund, #alpha=mean_rel_abund,
             fill=Phylum)) +
  geom_point(shape=21, color="black") +
  scale_size(range=c(-1,10), breaks=c(0,1,10,20,40,60,80,100), name="Relative\nabundance (%)") +
  #scale_size_area(max_size = 10, n.breaks=6, name="Relative\nabundance (%)") +
  scale_alpha(n.breaks=6, name="Relative\nabundance (%)", range = c(0.3,1)) +
  theme_bw(base_size=10) +
  scale_fill_manual(values = phylum_color2,
                    labels=c("Proteobacteria"="Pseudomonadota",
                             "Firmicutes"="Bacillota",
                             "Actinobacteria"="Actinomycetota",
                             "Bacteroidetes"="Bacteroidota")) +
  labs(x="Wisconsin Patient ID",
       y="Bacterial Genera") +
  guides(fill = guide_legend(override.aes = list(size = 7), order=1)) +
  theme(axis.text.y = element_text(face = "italic"))
pWisc # Pannel A 

ps.AUS<-qza_to_phyloseq( 
  features="./AUS-table2.qza",
  taxonomy="./AUS-taxonomy.qza",
  metadata = "./AUS-sample_metadata.txt",
  tree = "./AUs-rooted-tree.qza") # 16S samples from australian cohort
ps.decont <- subset_samples(ps.AUS,sample == "Patient")
ps.decont <- subset_taxa(ps.decont, is.na(Phylum) == F)
ps.decont <- subset_taxa(ps.decont, Family != 'mitochondria')
t <- as.data.frame(tax_table(ps.decont))
t_no_Genus <- t %>% filter(is.na(Genus) == T)
t_Family_as_Genus <- t %>%
  mutate(Genus = ifelse(is.na(Genus)==T, Family, Genus)) %>%
  mutate(Genus = ifelse(Genus == "1-68", "1-68 [Tissierellaceae]", Genus)) %>%
  mutate(Genus = ifelse(Genus == "ph2", "ph2 [Tissierellaceae]", Genus))
t_Family_as_Genus <- tax_table(as.matrix(t_Family_as_Genus))
tax_table(ps.decont) <- t_Family_as_Genus
sample_data(ps.decont)$subject <- as.character(sample_data(ps.decont)$subject)

ps.comp <- transform_sample_counts(ps.decont, function(x) x / sum(x))
df.ps = psmelt(ps.comp)
df.ps <- as.data.frame(df.ps)
df.ps.long <- pivot_longer(df.ps, c(Kingdom:Species), names_to = 'level', values_to = 'taxon') %>%
    mutate(subject = fct_relevel(as.character(subject),
                               "13", "15", "1", "2",
                               "18", "6", "11", "4", "14",
                               "16", "3", "7", "17",
                               "12", "5", "10"))
tax_Gen_Fam_Phylum <- df.ps %>% dplyr::select(Genus, Family, Phylum) %>% unique()
tax_Gen_Fam_Phylum[nrow(tax_Gen_Fam_Phylum)+1,] <- c("Other (<1%)", "Other (<1%)", "Other (<1%)")
tax_Gen_Fam_Phylum[nrow(tax_Gen_Fam_Phylum)+1,] <- c("Unclassified", "Other (<1%)", "Other (<1%)")

taxon_rel_abund <- df.ps.long %>%
  replace_na(list(taxon = "Unclassified", level = "Unclassified", Abundance = 0)) %>%
  filter(level=="Genus") %>%
  group_by(subject, taxon) %>%
  summarize(rel_abund = sum(Abundance), .groups="drop") %>%
  group_by(subject, taxon) %>%
  summarize(mean_rel_abund = mean(rel_abund)*100, .groups="drop")
taxon_pool <- taxon_rel_abund %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 5,
            mean = mean(mean_rel_abund),
            .groups="drop")

phylum_color2 = c("#C8443F","#175169","#8AB17D", "#EFB366")
pAus<-
  inner_join(taxon_rel_abund, taxon_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(subject, taxon) %>%
  summarize(mean_rel_abund = mean(mean_rel_abund),
            mean = min(mean),
            .groups="drop") %>%
  inner_join(., tax_Gen_Fam_Phylum, by=c('taxon'='Genus')) %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc=F),
        ) %>%
  ggplot(aes(x=subject, y=taxon, size=mean_rel_abund, #alpha=mean_rel_abund,
             fill=Phylum)) +
  geom_point(shape=21,color="black") +
  scale_size(range=c(-1,10), breaks=c(0,1,10,20,40,60,80,100)) +
  theme_bw(base_size=10) +
  scale_fill_manual(values = phylum_color2,
                    labels=c("Proteobacteria"="Pseudomonadota",
                             "Firmicutes"="Bacillota",
                             "Actinobacteria"="Actinomycetota",
                             "Bacteroidetes"="Bacteroidota"),
                    limits=c("Proteobacteria",
                             "Actinobacteria",
                             "Firmicutes"
                    )) +
  labs(x="Australia Subject ID",
       y="Bacterial Genera") +
  guides(fill = guide_legend(override.aes = list(size = 7)))+
  theme(axis.text.y = element_text(face = "italic"),
        axis.text = element_text(color="black"),
        legend.position = "none"
  )
pAus # Pannel B 

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


# Figure 6 and Supplemental Figure 7 
TopHOWclusters <- read.csv("./TopClustersInHOW.csv")
GO.HOW <- as.data.frame(TopHOWclusters)

GO.Healed <- subset.data.frame(GO.HOW, GO.HOW$Enriched.in.HOW == "Healed") 
GO.Ongoing <- subset.data.frame(GO.HOW, GO.HOW$Enriched.in.HOW == "Ongoing") 
GO.Worse <- subset.data.frame(GO.HOW , GO.HOW$Enriched.in.HOW == "Worse") 

Heal_BP_simMatrix <- calculateSimMatrix(GO.Healed$Goterm,
                                        orgdb="org.Hs.eg.db",
                                        ont="BP",
                                        method="Rel")
Heal_BP_scores <- setNames(-log10(GO.Healed$FDR), GO.Healed$Goterm)
Heal_BP_reducedTerms <- reduceSimMatrix(Heal_BP_simMatrix, Heal_BP_scores,
                                        threshold=0.75,
                                        orgdb="org.Hs.eg.db")

scatterPlot(Heal_BP_simMatrix, Heal_BP_reducedTerms)+ #Fig 6A
  scale_color_manual(values = c(colorRampPalette(c("#021327", "#053061","#2A9D8F","#8AB17D","#D29A0F",
                                                   "#F5A261","#DF7350","#C8443F", "#7E2825"))(20)))+
  
  theme(legend.position = "none")
pdf(file="./Heal_Top3Clusters_BP_tree.pdf") # Fig S7A
treemapPlot(Heal_BP_reducedTerms)
dev.off()

On_BP_simMatrix <- calculateSimMatrix(GO.Ongoing$Goterm,
                                      orgdb="org.Hs.eg.db",
                                      ont="BP",
                                      method="Rel")
On_BP_scores <- setNames(-log10(GO.Ongoing$rawP.value), GO.Ongoing$Goterm) 
On_BP_reducedTerms <- reduceSimMatrix(On_BP_simMatrix, On_BP_scores,
                                      threshold=0.9,
                                      orgdb="org.Hs.eg.db")

scatterPlot(On_BP_simMatrix, On_BP_reducedTerms)+ # Fig 6B
  scale_color_manual(values = c(colorRampPalette(c("#021327", "#053061","#2A9D8F","#8AB17D","#D29A0F",
                                                   "#F5A261","#DF7350","#C8443F", "#7E2825"))(11)))+
  theme(legend.position = "none")
pdf(file="./Ongoing_Top3Clusters_BP_tree.pdf") # Fig S7B 
treemapPlot(On_BP_reducedTerms)
dev.off()

Worse_BP_simMatrix <- calculateSimMatrix(GO.Worse$Goterm,
                                         orgdb="org.Hs.eg.db",
                                         ont="BP",
                                         method="Rel")
Worse_BP_scores <- setNames(-log10(GO.Worse$FDR), GO.Worse$Goterm)
Worse_BP_reducedTerms <- reduceSimMatrix(Worse_BP_simMatrix, Worse_BP_scores,
                                         threshold=0.7,
                                         orgdb="org.Hs.eg.db")

scatterPlot(Worse_BP_simMatrix, Worse_BP_reducedTerms)+ # Fig 6C
  scale_color_manual(values = c(colorRampPalette(c("#021327", "#053061","#2A9D8F","#8AB17D","#D29A0F",
                                                   "#F5A261","#DF7350","#C8443F", "#7E2825"))(19)))+
  theme(legend.position = "none")

pdf(file="./Worse_Top3Clusters_BP_tree.pdf") # Fig S7C
treemapPlot(Worse_BP_reducedTerms)
dev.off()


# supplemental Figure 8 
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
ggplot(KeyClusters, aes(x= Cluster, y = reorder(GO.biological.process,-ArbitraryNumberRank), color = Classification, fill = Classification)) +
  geom_point(size = 2.5)+
  facet_wrap(factor(Cluster) ~., scales = "free", ncol = 3) + # rows = 6, col = 4)+
  scale_color_manual(values = PCFcolors2)+
  scale_fill_manual(values = PCFcolors2)+
  theme_bw()+
  scale_x_discrete(expand = c(0.001,0))+
  theme(strip.background = element_blank(),strip.text.x = element_blank(),axis.text.y=element_text(size=2), axis.text.x=element_text(size=5))

KeyClu2 <- c("Cluster-6", "Cluster-21", "Cluster-11", "Cluster-22", "Cluster-19", "Cluster-5" )
KeyClusters2<-filter(TopGoTerms, Cluster %in% KeyClu2)
PCFcolors3 <- c(colorRampPalette(c("#412535", "#913C44","#C8443F","#DF7350","#EE925B"))(6), # reds
                colorRampPalette(c("#F5A261","#F6D27A"))(3), # Yellows
                colorRampPalette(c("#E9CE7B","#8AB17D", "#2A9D8F","#053061", "#31273D"))(8))

Cl_order <- c("Blood Clotting", "Wound Healing", "Homeostasis", "Immune Function", "Response to Stimulus",
              "Multicellular Organismal Process", "Biologic Process",
              "Biosynthetic Process", 
              "Intercellular Transport",
              "Intracellular Process",
              "Metabolic Process",
              "Molecular Function",
              "Cell Adheasion", "Cell Death", "Cell Motility"
              "Gene Expression", "Developmental Process",
              
              )
ggplot(KeyClusters2, aes(x=Classification,  fill = Classification))+
  geom_bar(position = "dodge", width =0.8)+
  facet_grid(. ~factor(Cluster) ,scales = "free_x",  space = "free_x")+ 
 # scale_color_manual(values = PCFcolors3)+
  scale_fill_manual(values = PCFcolors3)+
  theme_light()+
  scale_x_discrete(expand = c(0.001,0))+
  theme(axis.text.y=element_text(size=10), axis.text.x=element_text(size=1))





