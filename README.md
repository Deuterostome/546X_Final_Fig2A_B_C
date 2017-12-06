library(dplyr)
f2 = read.table("S_Table_2.txt", header = TRUE)
f2_df = tbl_df(f2)


# for ra1@1mm 
f2_ra1_1 = filter(f2_df, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
head(f2_ra1_1$gene_id)
str(f2_ra1_1$gene_id)
ra1_1 = as.character(f2_ra1_1$gene_id)
str(ra1_1)

# for ra2@1mm 
f2_ra2_1 = filter(f2_df, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
ra2_1 = as.character(f2_ra2_1$gene_id)

# for ra3@1mm 
f2_ra3_1 = filter(f2_df, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
ra3_1 = as.character(f2_ra3_1$gene_id)

# Venn Diagram for 1mm
library(VennDiagram)
venn.diagram(
  x = list(ra1_1 , ra2_1, ra3_1),
  category.names = c("ra1@1mm", "ra2@1mm", "ra3@1mm"),
  filename = 'F2_A_venn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 ,
  width = 480 ,
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  lty = 'blank',
  fill = c('yellow', 'purple', 'green'),
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

# Getting TF from Supplementary_Table6 with Linux
grep "EntrezGene" S_Table6.csv > TF.csv
grep "transcription" S_Table6.csv >> TF.csv

# Transcription Factor @1mm
library(dplyr)
f6_TF = read.csv("TF.csv", header = TRUE)
f6_TF_df = tbl_df(f6_TF)
TF = data.frame(f6_TF_df$maize.gene.id)

TF_1 = arrange(TF, f6_TF_df.maize.gene.id)
f2_df_1 = arrange(f2_df, gene_id)
f2_TF = merge(f2_df_1, TF_1, by.x = "gene_id", by.y = "f6_TF_df.maize.gene.id")


# for ra1@1mm 
f2_TF_ra1_1 = filter(f2_TF, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
ra1_TF_1 = as.character(f2_TF_ra1_1$gene_id)

# for ra2@1mm 
f2_TF_ra2_1 = filter(f2_TF, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
ra2_TF_1 = as.character(f2_TF_ra2_1$gene_id)

# for ra3@1mm 
f2_TF_ra3_1 = filter(f2_TF, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
ra3_TF_1 = as.character(f2_TF_ra3_1$gene_id)

#Venn Diagram for 1mm
library(VennDiagram)
venn.diagram(
  x = list(ra1_TF_1 , ra3_TF_1, ra2_TF_1),
  category.names = c("ra1_TF@1", "ra3_TF@1", "ra2_TF@1"),
  filename = 'F2_A_TF_1_venn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 ,
  width = 480 ,
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  lty = 'blank',
  fill = c('yellow', 'purple', 'green'),
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

############################################
# Transcription Factor@2mm
library(dplyr)
f6_TF = read.csv("TF.csv", header = TRUE)
f6_TF_df = tbl_df(f6_TF)
TF = data.frame(f6_TF_df$maize.gene.id)

TF_1 = arrange(TF, f6_TF_df.maize.gene.id)
f2_df_1 = arrange(f2_df, gene_id)
f2_TF = merge(f2_df_1, TF_1, by.x = "gene_id", by.y = "f6_TF_df.maize.gene.id")


# for ra1@2mm 
f2_TF_ra1_2 = filter(f2_TF, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
ra1_TF_2 = as.character(f2_TF_ra1_2$gene_id)

# for ra2@2mm 
f2_TF_ra2_2 = filter(f2_TF, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
ra2_TF_2 = as.character(f2_TF_ra2_2$gene_id)

# for ra3@2mm 
f2_TF_ra3_2 = filter(f2_TF, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")
ra3_TF_2 = as.character(f2_TF_ra3_2$gene_id)

#Venn Diagram for 2mm
library(VennDiagram)
venn.diagram(
  x = list(ra1_TF_2 , ra3_TF_2, ra2_TF_2),
  category.names = c("ra1_TF@2", "ra3_TF@2", "ra2_TF@2"),
  filename = 'F2_A_TF_2_venn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 ,
  width = 480 ,
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  lty = 'blank',
  fill = c('yellow', 'purple', 'green'),
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

###################################
# GO enrichment
source("https://bioconductor.org/biocLite.R")
biocLite("GO.db")
biocLite("biomaRt")
biocLite("topGO")
biocLite("GOSemSim")
library("biomaRt")
mart <- useMart(biomart = "plants_mart", host="plants.ensembl.org", dataset="zmays_eg_gene")
univ.geneID <- getBM(attributes=c("ensembl_gene_id", "entrezgene"), mart = mart) # 40481
univ.geneID
univ.geneID4
## remove genes with no corresponding Entrez Gene ID
univ.geneID2 <- univ.geneID[!is.na(univ.geneID[,2]),] # 14142 
## remove duplicated Entrez Gene ID
univ.geneID3 <- univ.geneID2[ !duplicated(univ.geneID2[,2]),] # 13630
##Get GO terms
univ.geneID4<-getBM(attributes=c("entrezgene","goslim_goa_accession","name_1006","namespace_1003","go_linkage_type"),filters='entrezgene',mart=mart,values=univ.geneID3$entrezgene)
# Start from here
head(univ.geneID4)
write.csv(univ.geneID4, file = "univ_geneID4.csv")
univ.geneID4[,2]


## Code evidence for genes without GO terms as NA
univ.geneID5<-univ.geneID4[-which(univ.geneID4[,2]==""),]
###Make dataframe for GOStats
goframeData <- data.frame(go_id = univ.geneID4$goslim_goa_accession, Evidence = univ.geneID4$go_linkage_type, gene_id = univ.geneID4$entrezgene,stringsAsFactors=F)




## Make my geneID
f2_ra1_1 = filter(f2_df, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
f2_ra1_1_cut = f2_ra1_1[, 1:9]
colnames(f2_ra1_1_cut)[1] <- "ensembl_gene_id"
f2_ra1_1_cut_A = arrange (f2_ra1_1_cut, ensembl_gene_id)
univ.geneID3_A = arrange (univ.geneID3, ensembl_gene_id)
f2_ra1_1_GO <- merge(f2_ra1_1_cut_A, univ.geneID3_A, by ="ensembl_gene_id") # 620
f2_ra1_1_GO

## GO enrichment
library("GOstats")
library("GOSemSim")
##Prepare GO to gene mappings
goFrame=GOFrame(goframeData,organism="Zea mays")
goAllFrame=GOAllFrame(goFrame)


library(GSEABase)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

params <- GSEAGOHyperGParams(name="Domestication Zea mays GO", geneSetCollection=gsc, geneIds = f2_ra1_1_GO[,10], universeGeneIds = univ.geneID4$entrezgene, ontology = "BP", pvalueCutoff = 0.05, conditional = TRUE, testDirection = "over")



f2_ra1_1_GO


my.geneID = read.csv("S2.csv", header=T,stringsAsFactors = F)
my.geneID <- my.geneID[,1:5]
colnames(my.geneID)[1] <- "ensembl_gene_id"
my.geneID_A = arrange (my.geneID, ensembl_gene_id)

my.geneID2 <- merge(my.geneID_A, univ.geneID3_A, by ="ensembl_gene_id") # 620
my.geneID2
my.geneID_A

my.geneID_A
univ.geneID3_A
## remove duplicated Entrez Gene ID
my.geneID3 <- my.geneID2[ !duplicated(my.geneID2$entrezgene),] # 620
