# Necessary packages and library
library(dplyr)
library(plyr)

# Read the supplementary table 2
f2 = read.table("S_Table_2.txt", header = TRUE)
f2_df = tbl_df(f2)

## Fig 2A
# creating the vector of genes_id of ra1@1mm 
f2_ra1_1 = filter(f2_df, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
head(f2_ra1_1$gene_id)
str(f2_ra1_1$gene_id)
ra1_1 = as.character(f2_ra1_1$gene_id)
str(ra1_1)
# creating the vector of genes_id of ra2@1mm 
f2_ra2_1 = filter(f2_df, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
ra2_1 = as.character(f2_ra2_1$gene_id)
# creating the vector of genes_id of ra3@1mm 
f2_ra3_1 = filter(f2_df, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
ra3_1 = as.character(f2_ra3_1$gene_id)

#Venn Diagram for 1mm
library(VennDiagram)
venn.diagram(
  x = list(ra1_1 , ra3_1, ra2_1),
  category.names = c("ra1@1mm", "ra3@1mm", "ra2@1mm"),
  filename = 'F2_A_1_venn.png',
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

# creating the vector of genes_id of ra1@2mm 
f2_ra1_2 = filter(f2_df, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
head(f2_ra1_2$gene_id)
ra1_2 = as.character(f2_ra1_2$gene_id)
# creating the vector of genes_id of ra2@1mm 
f2_ra2_2 = filter(f2_df, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
ra2_2 = as.character(f2_ra2_2$gene_id)
# creating the vector of genes_id of ra3@1mm 
f2_ra3_2 = filter(f2_df, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")
ra3_2 = as.character(f2_ra3_2$gene_id)

#Venn Diagram for 2mm
library(VennDiagram)
venn.diagram(
  x = list(ra1_2 , ra3_2, ra2_2),
  category.names = c("ra1@2mm", "ra3@2mm", "ra2@2mm"),
  filename = 'F2_A_2_venn.png',
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

##########################################
## Fig 2B
# Picking up Transcription Factors from Supplementary table 2
library(dplyr)
f6_TF = read.csv("TF.csv", header = TRUE) #grepping the transcription factor on Linux, creating the TF.csv files
f6_TF_df = tbl_df(f6_TF)
TF = data.frame(f6_TF_df$maize.gene.id)
TF_1 = arrange(TF, f6_TF_df.maize.gene.id)
f2_df_1 = arrange(f2_df, gene_id)
f2_TF = merge(f2_df_1, TF_1, by.x = "gene_id", by.y = "f6_TF_df.maize.gene.id")

# ra1_TF@1mm 
f2_TF_ra1_1 = filter(f2_TF, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
ra1_TF_1 = as.character(f2_TF_ra1_1$gene_id)
# for ra2_TF@1mm 
f2_TF_ra2_1 = filter(f2_TF, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
ra2_TF_1 = as.character(f2_TF_ra2_1$gene_id)
# for ra3_TF@1mm 
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

#Transcription Factor@2mm
library(dplyr)
f6_TF = read.csv("TF.csv", header = TRUE)
f6_TF_df = tbl_df(f6_TF)
TF = data.frame(f6_TF_df$maize.gene.id)
TF_1 = arrange(TF, f6_TF_df.maize.gene.id)
f2_df_1 = arrange(f2_df, gene_id)
f2_TF = merge(f2_df_1, TF_1, by.x = "gene_id", by.y = "f6_TF_df.maize.gene.id")

# for ra1_TF@1mm 
f2_TF_ra1_2 = filter(f2_TF, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
ra1_TF_2 = as.character(f2_TF_ra1_2$gene_id)

# for ra2_TF@1mm 
f2_TF_ra2_2 = filter(f2_TF, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
ra2_TF_2 = as.character(f2_TF_ra2_2$gene_id)

# for ra3_TF@1mm 
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

#################################################
## Under each GO category, picking up the overlapping genes expression between ra1, ra2, and ra3 strains (1mm and 2mm)
f2_df_1 = arrange(f2_df, gene_id)
# nucleosome
nucleosome = read.csv("nucleosome.csv", header = TRUE)
nucleosome = tbl_df(nucleosome)
nucleosome_df = data.frame(nucleosome$maize.gene.id)
nucleosome_1 = arrange(nucleosome_df, nucleosome.maize.gene.id)
f2_nucleosome = merge(f2_df_1, nucleosome_1, by.x = "gene_id", by.y = "nucleosome.maize.gene.id")
# for ra1@1mm 
nucleosome_GO_ra1_1 = filter(f2_nucleosome, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
nucleosome_GO_ra1_1 = arrange (nucleosome_GO_ra1_1, gene_id)
# for ra2@1mm 
nucleosome_GO_ra2_1 = filter(f2_nucleosome, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
nucleosome_GO_ra2_1 = arrange (nucleosome_GO_ra2_1, gene_id)
# for ra3@1mm 
nucleosome_GO_ra3_1 = filter(f2_nucleosome, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
nucleosome_GO_ra3_1 = arrange (nucleosome_GO_ra3_1, gene_id)
# for ra1@2mm 
nucleosome_GO_ra1_2 = filter(f2_nucleosome, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
nucleosome_GO_ra1_2 = arrange (nucleosome_GO_ra1_2, gene_id)
# for ra2@2mm 
nucleosome_GO_ra2_2 = filter(f2_nucleosome, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
nucleosome_GO_ra2_2 = arrange (nucleosome_GO_ra2_2, gene_id)
# for ra3@2mm 
nucleosome_GO_ra3_2 = filter(f2_nucleosome, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")
nucleosome_GO_ra3_2 = arrange (nucleosome_GO_ra3_2, gene_id)
# ra1+ra2@1mm
nucleosome_ra1_ra2_1mm = merge (nucleosome_GO_ra1_1, nucleosome_GO_ra2_1, by = "gene_id")
# ra1+ra3@1mm
nucleosome_ra1_ra3_1mm = merge (nucleosome_GO_ra1_1, nucleosome_GO_ra3_1, by = "gene_id")
# ra1+ra2@2mm
nucleosome_ra1_ra2_2mm = merge (nucleosome_GO_ra1_2, nucleosome_GO_ra2_2, by = "gene_id")
# ra1+ra3@2mm
nucleosome_ra1_ra3_2mm = merge (nucleosome_GO_ra1_2, nucleosome_GO_ra3_2, by = "gene_id")
# chromatin
chromatin = read.csv("chromatin.csv", header = TRUE)
chromatin = tbl_df(chromatin)
chromatin_df = data.frame(chromatin$maize.gene.id)
chromatin_1 = arrange(chromatin_df, chromatin.maize.gene.id)
f2_chromatin = merge(f2_df_1, chromatin_1, by.x = "gene_id", by.y = "chromatin.maize.gene.id")
# for ra1@1mm 
chromatin_GO_ra1_1 = filter(f2_chromatin, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
chromatin_GO_ra1_1 = arrange (chromatin_GO_ra1_1, gene_id)
# for ra2@1mm 
chromatin_GO_ra2_1 = filter(f2_chromatin, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
chromatin_GO_ra2_1 = arrange (chromatin_GO_ra2_1, gene_id)
# for ra3@1mm 
chromatin_GO_ra3_1 = filter(f2_chromatin, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
chromatin_GO_ra3_1 = arrange (chromatin_GO_ra3_1, gene_id)
# for ra1@2mm 
chromatin_GO_ra1_2 = filter(f2_chromatin, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
chromatin_GO_ra1_2 = arrange (chromatin_GO_ra1_2, gene_id)
# for ra2@2mm 
chromatin_GO_ra2_2 = filter(f2_chromatin, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
chromatin_GO_ra2_2 = arrange (chromatin_GO_ra2_2, gene_id)
# for ra3@2mm 
chromatin_GO_ra3_2 = filter(f2_chromatin, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")
chromatin_GO_ra3_2 = arrange (chromatin_GO_ra3_2, gene_id)
# ra1+ra2@1mm
chromatin_ra1_ra2_1mm = merge (chromatin_GO_ra1_1, chromatin_GO_ra2_1, by = "gene_id")
# ra1+ra3@1mm
chromatin_ra1_ra3_1mm = merge (chromatin_GO_ra1_1, chromatin_GO_ra3_1, by = "gene_id")
# ra1+ra2@2mm
chromatin_ra1_ra2_2mm = merge (chromatin_GO_ra1_2, chromatin_GO_ra2_2, by = "gene_id")
# ra1+ra3@2mm
chromatin_ra1_ra3_2mm = merge (chromatin_GO_ra1_2, chromatin_GO_ra3_2, by = "gene_id")
# vesicle
vesicle = read.csv("vesicle.csv", header = TRUE)
vesicle = tbl_df(vesicle)
vesicle_df = data.frame(vesicle$maize.gene.id)
vesicle_1 = arrange(vesicle_df, vesicle.maize.gene.id)
f2_vesicle = merge(f2_df_1, vesicle_1, by.x = "gene_id", by.y = "vesicle.maize.gene.id")
# for ra1@1mm 
vesicle_GO_ra1_1 = filter(f2_vesicle, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
vesicle_GO_ra1_1 = arrange (vesicle_GO_ra1_1, gene_id)
# for ra2@1mm 
vesicle_GO_ra2_1 = filter(f2_vesicle, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
vesicle_GO_ra2_1 = arrange (vesicle_GO_ra2_1, gene_id)
# for ra3@1mm 
vesicle_GO_ra3_1 = filter(f2_vesicle, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
vesicle_GO_ra3_1 = arrange (vesicle_GO_ra3_1, gene_id)
# for ra1@2mm 
vesicle_GO_ra1_2 = filter(f2_vesicle, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
vesicle_GO_ra1_2 = arrange (vesicle_GO_ra1_2, gene_id)
# for ra2@2mm 
vesicle_GO_ra2_2 = filter(f2_vesicle, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
vesicle_GO_ra2_2 = arrange (vesicle_GO_ra2_2, gene_id)
# for ra3@2mm 
vesicle_GO_ra3_2 = filter(f2_vesicle, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")
vesicle_GO_ra3_2 = arrange (vesicle_GO_ra3_2, gene_id)
# ra1+ra2@1mm
vesicle_ra1_ra2_1mm = merge (vesicle_GO_ra1_1, vesicle_GO_ra2_1, by = "gene_id")
# ra1+ra3@1mm
vesicle_ra1_ra3_1mm = merge (vesicle_GO_ra1_1, vesicle_GO_ra3_1, by = "gene_id")
# ra1+ra2@2mm
vesicle_ra1_ra2_2mm = merge (vesicle_GO_ra1_2, vesicle_GO_ra2_2, by = "gene_id")
# ra1+ra3@2mm
vesicle_ra1_ra3_2mm = merge (vesicle_GO_ra1_2, vesicle_GO_ra3_2, by = "gene_id")
# temperature
temperature = read.csv("temperature.csv", header = TRUE)
temperature = tbl_df(temperature)
temperature_df = data.frame(temperature$maize.gene.id)
temperature_1 = arrange(temperature_df, temperature.maize.gene.id)
f2_temperature = merge(f2_df_1, temperature_1, by.x = "gene_id", by.y = "temperature.maize.gene.id")
# for ra1@1mm 
temperature_GO_ra1_1 = filter(f2_temperature, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
temperature_GO_ra1_1 = arrange (temperature_GO_ra1_1, gene_id)
# for ra2@1mm 
temperature_GO_ra2_1 = filter(f2_temperature, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
temperature_GO_ra2_1 = arrange (temperature_GO_ra2_1, gene_id)
# for ra3@1mm 
temperature_GO_ra3_1 = filter(f2_temperature, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
temperature_GO_ra3_1 = arrange (temperature_GO_ra3_1, gene_id)
# for ra1@2mm 
temperature_GO_ra1_2 = filter(f2_temperature, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
temperature_GO_ra1_2 = arrange (temperature_GO_ra1_2, gene_id)
# for ra2@2mm 
temperature_GO_ra2_2 = filter(f2_temperature, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
temperature_GO_ra2_2 = arrange (temperature_GO_ra2_2, gene_id)
# for ra3@2mm 
temperature_GO_ra3_2 = filter(f2_temperature, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")
temperature_GO_ra3_2 = arrange (temperature_GO_ra3_2, gene_id)
# ra1+ra2@1mm
temperature_ra1_ra2_1mm = merge (temperature_GO_ra1_1, temperature_GO_ra2_1, by = "gene_id")
# ra1+ra3@1mm
temperature_ra1_ra3_1mm = merge (temperature_GO_ra1_1, temperature_GO_ra3_1, by = "gene_id")
# ra1+ra2@2mm
temperature_ra1_ra2_2mm = merge (temperature_GO_ra1_2, temperature_GO_ra2_2, by = "gene_id")
# ra1+ra3@2mm
temperature_ra1_ra3_2mm = merge (temperature_GO_ra1_2, temperature_GO_ra3_2, by = "gene_id")
# stress
stress = read.csv("stress.csv", header = TRUE)
stress = tbl_df(stress)
stress_df = data.frame(stress$maize.gene.id)
stress_1 = arrange(stress_df, stress.maize.gene.id)
f2_stress = merge(f2_df_1, stress_1, by.x = "gene_id", by.y = "stress.maize.gene.id")
# for ra1@1mm 
stress_GO_ra1_1 = filter(f2_stress, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
stress_GO_ra1_1 = arrange (stress_GO_ra1_1, gene_id)
# for ra2@1mm 
stress_GO_ra2_1 = filter(f2_stress, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
stress_GO_ra2_1 = arrange (stress_GO_ra2_1, gene_id)
# for ra3@1mm 
stress_GO_ra3_1 = filter(f2_stress, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
stress_GO_ra3_1 = arrange (stress_GO_ra3_1, gene_id)
# for ra1@2mm 
stress_GO_ra1_2 = filter(f2_stress, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
stress_GO_ra1_2 = arrange (stress_GO_ra1_2, gene_id)
# for ra2@2mm 
stress_GO_ra2_2 = filter(f2_stress, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
stress_GO_ra2_2 = arrange (stress_GO_ra2_2, gene_id)
# for ra3@2mm 
stress_GO_ra3_2 = filter(f2_stress, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")
stress_GO_ra3_2 = arrange (stress_GO_ra3_2, gene_id)
# ra1+ra2@1mm
stress_ra1_ra2_1mm = merge (stress_GO_ra1_1, stress_GO_ra2_1, by = "gene_id")
# ra1+ra3@1mm
stress_ra1_ra3_1mm = merge (stress_GO_ra1_1, stress_GO_ra3_1, by = "gene_id")
# ra1+ra2@2mm
stress_ra1_ra2_2mm = merge (stress_GO_ra1_2, stress_GO_ra2_2, by = "gene_id")
# ra1+ra3@2mm
stress_ra1_ra3_2mm = merge (stress_GO_ra1_2, stress_GO_ra3_2, by = "gene_id")
# transcription
transcription = read.csv("transcription.csv", header = TRUE)
transcription = tbl_df(transcription)
transcription_df = data.frame(transcription$maize.gene.id)
transcription_1 = arrange(transcription_df, transcription.maize.gene.id)
f2_transcription = merge(f2_df_1, transcription_1, by.x = "gene_id", by.y = "transcription.maize.gene.id")
# for ra1@1mm 
transcription_GO_ra1_1 = filter(f2_transcription, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
transcription_GO_ra1_1 = arrange (transcription_GO_ra1_1, gene_id)
# for ra2@1mm 
transcription_GO_ra2_1 = filter(f2_transcription, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
transcription_GO_ra2_1 = arrange (transcription_GO_ra2_1, gene_id)
# for ra3@1mm 
transcription_GO_ra3_1 = filter(f2_transcription, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
transcription_GO_ra3_1 = arrange (transcription_GO_ra3_1, gene_id)
# for ra1@2mm 
transcription_GO_ra1_2 = filter(f2_transcription, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
transcription_GO_ra1_2 = arrange (transcription_GO_ra1_2, gene_id)
# for ra2@2mm 
transcription_GO_ra2_2 = filter(f2_transcription, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
transcription_GO_ra2_2 = arrange (transcription_GO_ra2_2, gene_id)
# for ra3@2mm 
transcription_GO_ra3_2 = filter(f2_transcription, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")
transcription_GO_ra3_2 = arrange (transcription_GO_ra3_2, gene_id)
# ra1+ra2@1mm
transcription_ra1_ra2_1mm = merge (transcription_GO_ra1_1, transcription_GO_ra2_1, by = "gene_id")
# ra1+ra3@1mm
transcription_ra1_ra3_1mm = merge (transcription_GO_ra1_1, transcription_GO_ra3_1, by = "gene_id")
# ra1+ra2@2mm
transcription_ra1_ra2_2mm = merge (transcription_GO_ra1_2, transcription_GO_ra2_2, by = "gene_id")
# ra1+ra3@2mm
transcription_ra1_ra3_2mm = merge (transcription_GO_ra1_2, transcription_GO_ra3_2, by = "gene_id")
# nitrogen
nitrogen = read.csv("nitrogen.csv", header = TRUE)
nitrogen = tbl_df(nitrogen)
nitrogen_df = data.frame(nitrogen$maize.gene.id)
nitrogen_1 = arrange(nitrogen_df, nitrogen.maize.gene.id)
f2_nitrogen = merge(f2_df_1, nitrogen_1, by.x = "gene_id", by.y = "nitrogen.maize.gene.id")
# for ra1@1mm 
nitrogen_GO_ra1_1 = filter(f2_nitrogen, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
nitrogen_GO_ra1_1 = arrange (nitrogen_GO_ra1_1, gene_id)
# for ra2@1mm 
nitrogen_GO_ra2_1 = filter(f2_nitrogen, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
nitrogen_GO_ra2_1 = arrange (nitrogen_GO_ra2_1, gene_id)
# for ra3@1mm 
nitrogen_GO_ra3_1 = filter(f2_nitrogen, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
nitrogen_GO_ra3_1 = arrange (nitrogen_GO_ra3_1, gene_id)
# for ra1@2mm 
nitrogen_GO_ra1_2 = filter(f2_nitrogen, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
nitrogen_GO_ra1_2 = arrange (nitrogen_GO_ra1_2, gene_id)
# for ra2@2mm 
nitrogen_GO_ra2_2 = filter(f2_nitrogen, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
nitrogen_GO_ra2_2 = arrange (nitrogen_GO_ra2_2, gene_id)
# for ra3@2mm 
nitrogen_GO_ra3_2 = filter(f2_nitrogen, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")
nitrogen_GO_ra3_2 = arrange (nitrogen_GO_ra3_2, gene_id)
# ra1+ra2@1mm
nitrogen_ra1_ra2_1mm = merge (nitrogen_GO_ra1_1, nitrogen_GO_ra2_1, by = "gene_id")
# ra1+ra3@1mm
nitrogen_ra1_ra3_1mm = merge (nitrogen_GO_ra1_1, nitrogen_GO_ra3_1, by = "gene_id")
# ra1+ra2@2mm
nitrogen_ra1_ra2_2mm = merge (nitrogen_GO_ra1_2, nitrogen_GO_ra2_2, by = "gene_id")
# ra1+ra3@2mm
nitrogen_ra1_ra3_2mm = merge (nitrogen_GO_ra1_2, nitrogen_GO_ra3_2, by = "gene_id")
# RNA_synthetase
RNA_synthetase = read.csv("RNA_synthetase.csv", header = TRUE)
RNA_synthetase = tbl_df(RNA_synthetase)
RNA_synthetase_df = data.frame(RNA_synthetase$maize.gene.id)
RNA_synthetase_1 = arrange(RNA_synthetase_df, RNA_synthetase.maize.gene.id)
f2_RNA_synthetase = merge(f2_df_1, RNA_synthetase_1, by.x = "gene_id", by.y = "RNA_synthetase.maize.gene.id")
# for ra1@1mm 
RNA_synthetase_GO_ra1_1 = filter(f2_RNA_synthetase, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
RNA_synthetase_GO_ra1_1 = arrange (RNA_synthetase_GO_ra1_1, gene_id)
# for ra2@1mm 
RNA_synthetase_GO_ra2_1 = filter(f2_RNA_synthetase, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
RNA_synthetase_GO_ra2_1 = arrange (RNA_synthetase_GO_ra2_1, gene_id)
# for ra3@1mm 
RNA_synthetase_GO_ra3_1 = filter(f2_RNA_synthetase, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
RNA_synthetase_GO_ra3_1 = arrange (RNA_synthetase_GO_ra3_1, gene_id)
# for ra1@2mm 
RNA_synthetase_GO_ra1_2 = filter(f2_RNA_synthetase, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
RNA_synthetase_GO_ra1_2 = arrange (RNA_synthetase_GO_ra1_2, gene_id)
# for ra2@2mm 
RNA_synthetase_GO_ra2_2 = filter(f2_RNA_synthetase, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
RNA_synthetase_GO_ra2_2 = arrange (RNA_synthetase_GO_ra2_2, gene_id)
# for ra3@2mm 
RNA_synthetase_GO_ra3_2 = filter(f2_RNA_synthetase, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")
RNA_synthetase_GO_ra3_2 = arrange (RNA_synthetase_GO_ra3_2, gene_id)
# ra1+ra2@1mm
RNA_synthetase_ra1_ra2_1mm = merge (RNA_synthetase_GO_ra1_1, RNA_synthetase_GO_ra2_1, by = "gene_id")
# ra1+ra3@1mm
RNA_synthetase_ra1_ra3_1mm = merge (RNA_synthetase_GO_ra1_1, RNA_synthetase_GO_ra3_1, by = "gene_id")
# ra1+ra2@2mm
RNA_synthetase_ra1_ra2_2mm = merge (RNA_synthetase_GO_ra1_2, RNA_synthetase_GO_ra2_2, by = "gene_id")
# ra1+ra3@2mm
RNA_synthetase_ra1_ra3_2mm = merge (RNA_synthetase_GO_ra1_2, RNA_synthetase_GO_ra3_2, by = "gene_id")
# G_protein
G_protein = read.csv("G_protein.csv", header = TRUE)
G_protein = tbl_df(G_protein)
G_protein_df = data.frame(G_protein$maize.gene.id)
G_protein_1 = arrange(G_protein_df, G_protein.maize.gene.id)
f2_G_protein = merge(f2_df_1, G_protein_1, by.x = "gene_id", by.y = "G_protein.maize.gene.id")
# for ra1@1mm 
G_protein_GO_ra1_1 = filter(f2_G_protein, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
G_protein_GO_ra1_1 = arrange (G_protein_GO_ra1_1, gene_id)
# for ra2@1mm 
G_protein_GO_ra2_1 = filter(f2_G_protein, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
G_protein_GO_ra2_1 = arrange (G_protein_GO_ra2_1, gene_id)
# for ra3@1mm 
G_protein_GO_ra3_1 = filter(f2_G_protein, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
G_protein_GO_ra3_1 = arrange (G_protein_GO_ra3_1, gene_id)
# for ra1@2mm 
G_protein_GO_ra1_2 = filter(f2_G_protein, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
G_protein_GO_ra1_2 = arrange (G_protein_GO_ra1_2, gene_id)
# for ra2@2mm 
G_protein_GO_ra2_2 = filter(f2_G_protein, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
G_protein_GO_ra2_2 = arrange (G_protein_GO_ra2_2, gene_id)
# for ra3@2mm 
G_protein_GO_ra3_2 = filter(f2_G_protein, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")
G_protein_GO_ra3_2 = arrange (G_protein_GO_ra3_2, gene_id)
# ra1+ra2@1mm
G_protein_ra1_ra2_1mm = merge (G_protein_GO_ra1_1, G_protein_GO_ra2_1, by = "gene_id")
# ra1+ra3@1mm
G_protein_ra1_ra3_1mm = merge (G_protein_GO_ra1_1, G_protein_GO_ra3_1, by = "gene_id")
# ra1+ra2@2mm
G_protein_ra1_ra2_2mm = merge (G_protein_GO_ra1_2, G_protein_GO_ra2_2, by = "gene_id")
# ra1+ra3@2mm
G_protein_ra1_ra3_2mm = merge (G_protein_GO_ra1_2, G_protein_GO_ra3_2, by = "gene_id")
# light
light = read.csv("light.csv", header = TRUE)
light = tbl_df(light)
light_df = data.frame(light$maize.gene.id)
light_1 = arrange(light_df, light.maize.gene.id)
f2_light = merge(f2_df_1, light_1, by.x = "gene_id", by.y = "light.maize.gene.id")
# for ra1@1mm 
light_GO_ra1_1 = filter(f2_light, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
light_GO_ra1_1 = arrange (light_GO_ra1_1, gene_id)
# for ra2@1mm 
light_GO_ra2_1 = filter(f2_light, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
light_GO_ra2_1 = arrange (light_GO_ra2_1, gene_id)
# for ra3@1mm 
light_GO_ra3_1 = filter(f2_light, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
light_GO_ra3_1 = arrange (light_GO_ra3_1, gene_id)
# for ra1@2mm 
light_GO_ra1_2 = filter(f2_light, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
light_GO_ra1_2 = arrange (light_GO_ra1_2, gene_id)
# for ra2@2mm 
light_GO_ra2_2 = filter(f2_light, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
light_GO_ra2_2 = arrange (light_GO_ra2_2, gene_id)
# for ra3@2mm 
light_GO_ra3_2 = filter(f2_light, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")
light_GO_ra3_2 = arrange (light_GO_ra3_2, gene_id)
# ra1+ra2@1mm
light_ra1_ra2_1mm = merge (light_GO_ra1_1, light_GO_ra2_1, by = "gene_id")
# ra1+ra3@1mm
light_ra1_ra3_1mm = merge (light_GO_ra1_1, light_GO_ra3_1, by = "gene_id")
# ra1+ra2@2mm
light_ra1_ra2_2mm = merge (light_GO_ra1_2, light_GO_ra2_2, by = "gene_id")
# ra1+ra3@2mm
light_ra1_ra3_2mm = merge (light_GO_ra1_2, light_GO_ra3_2, by = "gene_id")
# redox
redox = read.csv("redox.csv", header = TRUE)
redox = tbl_df(redox)
redox_df = data.frame(redox$maize.gene.id)
redox_1 = arrange(redox_df, redox.maize.gene.id)
f2_redox = merge(f2_df_1, redox_1, by.x = "gene_id", by.y = "redox.maize.gene.id")
# for ra1@1mm 
redox_GO_ra1_1 = filter(f2_redox, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
redox_GO_ra1_1 = arrange (redox_GO_ra1_1, gene_id)
# for ra2@1mm 
redox_GO_ra2_1 = filter(f2_redox, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
redox_GO_ra2_1 = arrange (redox_GO_ra2_1, gene_id)
# for ra3@1mm 
redox_GO_ra3_1 = filter(f2_redox, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
redox_GO_ra3_1 = arrange (redox_GO_ra3_1, gene_id)
# for ra1@2mm 
redox_GO_ra1_2 = filter(f2_redox, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
redox_GO_ra1_2 = arrange (redox_GO_ra1_2, gene_id)
# for ra2@2mm 
redox_GO_ra2_2 = filter(f2_redox, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
redox_GO_ra2_2 = arrange (redox_GO_ra2_2, gene_id)
# for ra3@2mm 
redox_GO_ra3_2 = filter(f2_redox, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")
redox_GO_ra3_2 = arrange (redox_GO_ra3_2, gene_id)
# ra1+ra2@1mm
redox_ra1_ra2_1mm = merge (redox_GO_ra1_1, redox_GO_ra2_1, by = "gene_id")
# ra1+ra3@1mm
redox_ra1_ra3_1mm = merge (redox_GO_ra1_1, redox_GO_ra3_1, by = "gene_id")
# ra1+ra2@2mm
redox_ra1_ra2_2mm = merge (redox_GO_ra1_2, redox_GO_ra2_2, by = "gene_id")
# ra1+ra3@2mm
redox_ra1_ra3_2mm = merge (redox_GO_ra1_2, redox_GO_ra3_2, by = "gene_id")
# monosaccharide
monosaccharide = read.csv("monosaccharide.csv", header = TRUE)
monosaccharide = tbl_df(monosaccharide)
monosaccharide_df = data.frame(monosaccharide$maize.gene.id)
monosaccharide_1 = arrange(monosaccharide_df, monosaccharide.maize.gene.id)
f2_monosaccharide = merge(f2_df_1, monosaccharide_1, by.x = "gene_id", by.y = "monosaccharide.maize.gene.id")
# for ra1@1mm 
monosaccharide_GO_ra1_1 = filter(f2_monosaccharide, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
monosaccharide_GO_ra1_1 = arrange (monosaccharide_GO_ra1_1, gene_id)
# for ra2@1mm 
monosaccharide_GO_ra2_1 = filter(f2_monosaccharide, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
monosaccharide_GO_ra2_1 = arrange (monosaccharide_GO_ra2_1, gene_id)
# for ra3@1mm 
monosaccharide_GO_ra3_1 = filter(f2_monosaccharide, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
monosaccharide_GO_ra3_1 = arrange (monosaccharide_GO_ra3_1, gene_id)
# for ra1@2mm 
monosaccharide_GO_ra1_2 = filter(f2_monosaccharide, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
monosaccharide_GO_ra1_2 = arrange (monosaccharide_GO_ra1_2, gene_id)
# for ra2@2mm 
monosaccharide_GO_ra2_2 = filter(f2_monosaccharide, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
monosaccharide_GO_ra2_2 = arrange (monosaccharide_GO_ra2_2, gene_id)
# for ra3@2mm 
monosaccharide_GO_ra3_2 = filter(f2_monosaccharide, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")
monosaccharide_GO_ra3_2 = arrange (monosaccharide_GO_ra3_2, gene_id)
# ra1+ra2@1mm
monosaccharide_ra1_ra2_1mm = merge (monosaccharide_GO_ra1_1, monosaccharide_GO_ra2_1, by = "gene_id")
# ra1+ra3@1mm
monosaccharide_ra1_ra3_1mm = merge (monosaccharide_GO_ra1_1, monosaccharide_GO_ra3_1, by = "gene_id")
# ra1+ra2@2mm
monosaccharide_ra1_ra2_2mm = merge (monosaccharide_GO_ra1_2, monosaccharide_GO_ra2_2, by = "gene_id")
# ra1+ra3@2mm
monosaccharide_ra1_ra3_2mm = merge (monosaccharide_GO_ra1_2, monosaccharide_GO_ra3_2, by = "gene_id")
# cell_wall
cell_wall = read.csv("cell_wall.csv", header = TRUE)
cell_wall = tbl_df(cell_wall)
cell_wall_df = data.frame(cell_wall$maize.gene.id)
cell_wall_1 = arrange(cell_wall_df, cell_wall.maize.gene.id)
f2_cell_wall = merge(f2_df_1, cell_wall_1, by.x = "gene_id", by.y = "cell_wall.maize.gene.id")
# for ra1@1mm 
cell_wall_GO_ra1_1 = filter(f2_cell_wall, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
cell_wall_GO_ra1_1 = arrange (cell_wall_GO_ra1_1, gene_id)
# for ra2@1mm 
cell_wall_GO_ra2_1 = filter(f2_cell_wall, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
cell_wall_GO_ra2_1 = arrange (cell_wall_GO_ra2_1, gene_id)
# for ra3@1mm 
cell_wall_GO_ra3_1 = filter(f2_cell_wall, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
cell_wall_GO_ra3_1 = arrange (cell_wall_GO_ra3_1, gene_id)
# for ra1@2mm 
cell_wall_GO_ra1_2 = filter(f2_cell_wall, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
cell_wall_GO_ra1_2 = arrange (cell_wall_GO_ra1_2, gene_id)
# for ra2@2mm 
cell_wall_GO_ra2_2 = filter(f2_cell_wall, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
cell_wall_GO_ra2_2 = arrange (cell_wall_GO_ra2_2, gene_id)
# for ra3@2mm 
cell_wall_GO_ra3_2 = filter(f2_cell_wall, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")
cell_wall_GO_ra3_2 = arrange (cell_wall_GO_ra3_2, gene_id)
# ra1+ra2@1mm
cell_wall_ra1_ra2_1mm = merge (cell_wall_GO_ra1_1, cell_wall_GO_ra2_1, by = "gene_id")
# ra1+ra3@1mm
cell_wall_ra1_ra3_1mm = merge (cell_wall_GO_ra1_1, cell_wall_GO_ra3_1, by = "gene_id")
# ra1+ra2@2mm
cell_wall_ra1_ra2_2mm = merge (cell_wall_GO_ra1_2, cell_wall_GO_ra2_2, by = "gene_id")
# ra1+ra3@2mm
cell_wall_ra1_ra3_2mm = merge (cell_wall_GO_ra1_2, cell_wall_GO_ra3_2, by = "gene_id")
# sucrose
sucrose = read.csv("sucrose.csv", header = TRUE)
sucrose = tbl_df(sucrose)
sucrose_df = data.frame(sucrose$maize.gene.id)
sucrose_1 = arrange(sucrose_df, sucrose.maize.gene.id)
f2_sucrose = merge(f2_df_1, sucrose_1, by.x = "gene_id", by.y = "sucrose.maize.gene.id")
# for ra1@1mm 
sucrose_GO_ra1_1 = filter(f2_sucrose, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
sucrose_GO_ra1_1 = arrange (sucrose_GO_ra1_1, gene_id)
# for ra2@1mm 
sucrose_GO_ra2_1 = filter(f2_sucrose, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
sucrose_GO_ra2_1 = arrange (sucrose_GO_ra2_1, gene_id)
# for ra3@1mm 
sucrose_GO_ra3_1 = filter(f2_sucrose, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
sucrose_GO_ra3_1 = arrange (sucrose_GO_ra3_1, gene_id)
# for ra1@2mm 
sucrose_GO_ra1_2 = filter(f2_sucrose, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
sucrose_GO_ra1_2 = arrange (sucrose_GO_ra1_2, gene_id)
# for ra2@2mm 
sucrose_GO_ra2_2 = filter(f2_sucrose, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
sucrose_GO_ra2_2 = arrange (sucrose_GO_ra2_2, gene_id)
# for ra3@2mm 
sucrose_GO_ra3_2 = filter(f2_sucrose, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")
sucrose_GO_ra3_2 = arrange (sucrose_GO_ra3_2, gene_id)
# ra1+ra2@1mm
sucrose_ra1_ra2_1mm = merge (sucrose_GO_ra1_1, sucrose_GO_ra2_1, by = "gene_id")
# ra1+ra3@1mm
sucrose_ra1_ra3_1mm = merge (sucrose_GO_ra1_1, sucrose_GO_ra3_1, by = "gene_id")
# ra1+ra2@2mm
sucrose_ra1_ra2_2mm = merge (sucrose_GO_ra1_2, sucrose_GO_ra2_2, by = "gene_id")
# ra1+ra3@2mm
sucrose_ra1_ra3_2mm = merge (sucrose_GO_ra1_2, sucrose_GO_ra3_2, by = "gene_id")
# trehalose
trehalose = read.csv("trehalose.csv", header = TRUE)
trehalose = tbl_df(trehalose)
trehalose_df = data.frame(trehalose$maize.gene.id)
trehalose_1 = arrange(trehalose_df, trehalose.maize.gene.id)
f2_trehalose = merge(f2_df_1, trehalose_1, by.x = "gene_id", by.y = "trehalose.maize.gene.id")
# for ra1@1mm 
trehalose_GO_ra1_1 = filter(f2_trehalose, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
trehalose_GO_ra1_1 = arrange (trehalose_GO_ra1_1, gene_id)
# for ra2@1mm 
trehalose_GO_ra2_1 = filter(f2_trehalose, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
trehalose_GO_ra2_1 = arrange (trehalose_GO_ra2_1, gene_id)
# for ra3@1mm 
trehalose_GO_ra3_1 = filter(f2_trehalose, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
trehalose_GO_ra3_1 = arrange (trehalose_GO_ra3_1, gene_id)
# for ra1@2mm 
trehalose_GO_ra1_2 = filter(f2_trehalose, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
trehalose_GO_ra1_2 = arrange (trehalose_GO_ra1_2, gene_id)
# for ra2@2mm 
trehalose_GO_ra2_2 = filter(f2_trehalose, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
trehalose_GO_ra2_2 = arrange (trehalose_GO_ra2_2, gene_id)
# for ra3@2mm 
trehalose_GO_ra3_2 = filter(f2_trehalose, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")
trehalose_GO_ra3_2 = arrange (trehalose_GO_ra3_2, gene_id)
# ra1+ra2@1mm
trehalose_ra1_ra2_1mm = merge (trehalose_GO_ra1_1, trehalose_GO_ra2_1, by = "gene_id")
# ra1+ra3@1mm
trehalose_ra1_ra3_1mm = merge (trehalose_GO_ra1_1, trehalose_GO_ra3_1, by = "gene_id")
# ra1+ra2@2mm
trehalose_ra1_ra2_2mm = merge (trehalose_GO_ra1_2, trehalose_GO_ra2_2, by = "gene_id")
# ra1+ra3@2mm
trehalose_ra1_ra3_2mm = merge (trehalose_GO_ra1_2, trehalose_GO_ra3_2, by = "gene_id")
# intracellular_protein
intracellular_protein = read.csv("intracellular_protein.csv", header = TRUE)
intracellular_protein = tbl_df(intracellular_protein)
intracellular_protein_df = data.frame(intracellular_protein$maize.gene.id)
intracellular_protein_1 = arrange(intracellular_protein_df, intracellular_protein.maize.gene.id)
f2_intracellular_protein = merge(f2_df_1, intracellular_protein_1, by.x = "gene_id", by.y = "intracellular_protein.maize.gene.id")
# for ra1@1mm 
intracellular_protein_GO_ra1_1 = filter(f2_intracellular_protein, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
intracellular_protein_GO_ra1_1 = arrange (intracellular_protein_GO_ra1_1, gene_id)
# for ra2@1mm 
intracellular_protein_GO_ra2_1 = filter(f2_intracellular_protein, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
intracellular_protein_GO_ra2_1 = arrange (intracellular_protein_GO_ra2_1, gene_id)
# for ra3@1mm 
intracellular_protein_GO_ra3_1 = filter(f2_intracellular_protein, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
intracellular_protein_GO_ra3_1 = arrange (intracellular_protein_GO_ra3_1, gene_id)
# for ra1@2mm 
intracellular_protein_GO_ra1_2 = filter(f2_intracellular_protein, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
intracellular_protein_GO_ra1_2 = arrange (intracellular_protein_GO_ra1_2, gene_id)
# for ra2@2mm 
intracellular_protein_GO_ra2_2 = filter(f2_intracellular_protein, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
intracellular_protein_GO_ra2_2 = arrange (intracellular_protein_GO_ra2_2, gene_id)
# for ra3@2mm 
intracellular_protein_GO_ra3_2 = filter(f2_intracellular_protein, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")
intracellular_protein_GO_ra3_2 = arrange (intracellular_protein_GO_ra3_2, gene_id)
# ra1+ra2@1mm
intracellular_protein_ra1_ra2_1mm = merge (intracellular_protein_GO_ra1_1, intracellular_protein_GO_ra2_1, by = "gene_id")
# ra1+ra3@1mm
intracellular_protein_ra1_ra3_1mm = merge (intracellular_protein_GO_ra1_1, intracellular_protein_GO_ra3_1, by = "gene_id")
# ra1+ra2@2mm
intracellular_protein_ra1_ra2_2mm = merge (intracellular_protein_GO_ra1_2, intracellular_protein_GO_ra2_2, by = "gene_id")
# ra1+ra3@2mm
intracellular_protein_ra1_ra3_2mm = merge (intracellular_protein_GO_ra1_2, intracellular_protein_GO_ra3_2, by = "gene_id")
# small_GTPase
small_GTPase = read.csv("small_GTPase.csv", header = TRUE)
small_GTPase = tbl_df(small_GTPase)
small_GTPase_df = data.frame(small_GTPase$maize.gene.id)
small_GTPase_1 = arrange(small_GTPase_df, small_GTPase.maize.gene.id)
f2_small_GTPase = merge(f2_df_1, small_GTPase_1, by.x = "gene_id", by.y = "small_GTPase.maize.gene.id")
# for ra1@1mm 
small_GTPase_GO_ra1_1 = filter(f2_small_GTPase, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
small_GTPase_GO_ra1_1 = arrange (small_GTPase_GO_ra1_1, gene_id)
# for ra2@1mm 
small_GTPase_GO_ra2_1 = filter(f2_small_GTPase, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
small_GTPase_GO_ra2_1 = arrange (small_GTPase_GO_ra2_1, gene_id)
# for ra3@1mm 
small_GTPase_GO_ra3_1 = filter(f2_small_GTPase, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
small_GTPase_GO_ra3_1 = arrange (small_GTPase_GO_ra3_1, gene_id)
# for ra1@2mm 
small_GTPase_GO_ra1_2 = filter(f2_small_GTPase, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
small_GTPase_GO_ra1_2 = arrange (small_GTPase_GO_ra1_2, gene_id)
# for ra2@2mm 
small_GTPase_GO_ra2_2 = filter(f2_small_GTPase, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
small_GTPase_GO_ra2_2 = arrange (small_GTPase_GO_ra2_2, gene_id)
# for ra3@2mm 
small_GTPase_GO_ra3_2 = filter(f2_small_GTPase, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")
small_GTPase_GO_ra3_2 = arrange (small_GTPase_GO_ra3_2, gene_id)
# ra1+ra2@1mm
small_GTPase_ra1_ra2_1mm = merge (small_GTPase_GO_ra1_1, small_GTPase_GO_ra2_1, by = "gene_id")
# ra1+ra3@1mm
small_GTPase_ra1_ra3_1mm = merge (small_GTPase_GO_ra1_1, small_GTPase_GO_ra3_1, by = "gene_id")
# ra1+ra2@2mm
small_GTPase_ra1_ra2_2mm = merge (small_GTPase_GO_ra1_2, small_GTPase_GO_ra2_2, by = "gene_id")
# ra1+ra3@2mm
small_GTPase_ra1_ra3_2mm = merge (small_GTPase_GO_ra1_2, small_GTPase_GO_ra3_2, by = "gene_id")
# intracellular_signaling
intracellular_signaling = read.csv("intracellular_signaling.csv", header = TRUE)
intracellular_signaling = tbl_df(intracellular_signaling)
intracellular_signaling_df = data.frame(intracellular_signaling$maize.gene.id)
intracellular_signaling_1 = arrange(intracellular_signaling_df, intracellular_signaling.maize.gene.id)
f2_intracellular_signaling = merge(f2_df_1, intracellular_signaling_1, by.x = "gene_id", by.y = "intracellular_signaling.maize.gene.id")
# for ra1@1mm 
intracellular_signaling_GO_ra1_1 = filter(f2_intracellular_signaling, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
intracellular_signaling_GO_ra1_1 = arrange (intracellular_signaling_GO_ra1_1, gene_id)
# for ra2@1mm 
intracellular_signaling_GO_ra2_1 = filter(f2_intracellular_signaling, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
intracellular_signaling_GO_ra2_1 = arrange (intracellular_signaling_GO_ra2_1, gene_id)
# for ra3@1mm 
intracellular_signaling_GO_ra3_1 = filter(f2_intracellular_signaling, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
intracellular_signaling_GO_ra3_1 = arrange (intracellular_signaling_GO_ra3_1, gene_id)
# for ra1@2mm 
intracellular_signaling_GO_ra1_2 = filter(f2_intracellular_signaling, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
intracellular_signaling_GO_ra1_2 = arrange (intracellular_signaling_GO_ra1_2, gene_id)
# for ra2@2mm 
intracellular_signaling_GO_ra2_2 = filter(f2_intracellular_signaling, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
intracellular_signaling_GO_ra2_2 = arrange (intracellular_signaling_GO_ra2_2, gene_id)
# for ra3@2mm 
intracellular_signaling_GO_ra3_2 = filter(f2_intracellular_signaling, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")
intracellular_signaling_GO_ra3_2 = arrange (intracellular_signaling_GO_ra3_2, gene_id)
# ra1+ra2@1mm
intracellular_signaling_ra1_ra2_1mm = merge (intracellular_signaling_GO_ra1_1, intracellular_signaling_GO_ra2_1, by = "gene_id")
# ra1+ra3@1mm
intracellular_signaling_ra1_ra3_1mm = merge (intracellular_signaling_GO_ra1_1, intracellular_signaling_GO_ra3_1, by = "gene_id")
# ra1+ra2@2mm
intracellular_signaling_ra1_ra2_2mm = merge (intracellular_signaling_GO_ra1_2, intracellular_signaling_GO_ra2_2, by = "gene_id")
# ra1+ra3@2mm
intracellular_signaling_ra1_ra3_2mm = merge (intracellular_signaling_GO_ra1_2, intracellular_signaling_GO_ra3_2, by = "gene_id")
# proteasome
proteasome = read.csv("proteasome.csv", header = TRUE)
proteasome = tbl_df(proteasome)
proteasome_df = data.frame(proteasome$maize.gene.id)
proteasome_1 = arrange(proteasome_df, proteasome.maize.gene.id)
f2_proteasome = merge(f2_df_1, proteasome_1, by.x = "gene_id", by.y = "proteasome.maize.gene.id")
# for ra1@1mm 
proteasome_GO_ra1_1 = filter(f2_proteasome, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
proteasome_GO_ra1_1 = arrange (proteasome_GO_ra1_1, gene_id)
# for ra2@1mm 
proteasome_GO_ra2_1 = filter(f2_proteasome, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
proteasome_GO_ra2_1 = arrange (proteasome_GO_ra2_1, gene_id)
# for ra3@1mm 
proteasome_GO_ra3_1 = filter(f2_proteasome, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
proteasome_GO_ra3_1 = arrange (proteasome_GO_ra3_1, gene_id)
# for ra1@2mm 
proteasome_GO_ra1_2 = filter(f2_proteasome, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
proteasome_GO_ra1_2 = arrange (proteasome_GO_ra1_2, gene_id)
# for ra2@2mm 
proteasome_GO_ra2_2 = filter(f2_proteasome, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
proteasome_GO_ra2_2 = arrange (proteasome_GO_ra2_2, gene_id)
# for ra3@2mm 
proteasome_GO_ra3_2 = filter(f2_proteasome, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")
proteasome_GO_ra3_2 = arrange (proteasome_GO_ra3_2, gene_id)
# ra1+ra2@1mm
proteasome_ra1_ra2_1mm = merge (proteasome_GO_ra1_1, proteasome_GO_ra2_1, by = "gene_id")
# ra1+ra3@1mm
proteasome_ra1_ra3_1mm = merge (proteasome_GO_ra1_1, proteasome_GO_ra3_1, by = "gene_id")
# ra1+ra2@2mm
proteasome_ra1_ra2_2mm = merge (proteasome_GO_ra1_2, proteasome_GO_ra2_2, by = "gene_id")
# ra1+ra3@2mm
proteasome_ra1_ra3_2mm = merge (proteasome_GO_ra1_2, proteasome_GO_ra3_2, by = "gene_id")
# nucleocytoplasmic
nucleocytoplasmic = read.csv("nucleocytoplasmic.csv", header = TRUE)
nucleocytoplasmic = tbl_df(nucleocytoplasmic)
nucleocytoplasmic_df = data.frame(nucleocytoplasmic$maize.gene.id)
nucleocytoplasmic_1 = arrange(nucleocytoplasmic_df, nucleocytoplasmic.maize.gene.id)
f2_nucleocytoplasmic = merge(f2_df_1, nucleocytoplasmic_1, by.x = "gene_id", by.y = "nucleocytoplasmic.maize.gene.id")
# for ra1@1mm 
nucleocytoplasmic_GO_ra1_1 = filter(f2_nucleocytoplasmic, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
nucleocytoplasmic_GO_ra1_1 = arrange (nucleocytoplasmic_GO_ra1_1, gene_id)
# for ra2@1mm 
nucleocytoplasmic_GO_ra2_1 = filter(f2_nucleocytoplasmic, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
nucleocytoplasmic_GO_ra2_1 = arrange (nucleocytoplasmic_GO_ra2_1, gene_id)
# for ra3@1mm 
nucleocytoplasmic_GO_ra3_1 = filter(f2_nucleocytoplasmic, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
nucleocytoplasmic_GO_ra3_1 = arrange (nucleocytoplasmic_GO_ra3_1, gene_id)
# for ra1@2mm 
nucleocytoplasmic_GO_ra1_2 = filter(f2_nucleocytoplasmic, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
nucleocytoplasmic_GO_ra1_2 = arrange (nucleocytoplasmic_GO_ra1_2, gene_id)
# for ra2@2mm 
nucleocytoplasmic_GO_ra2_2 = filter(f2_nucleocytoplasmic, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
nucleocytoplasmic_GO_ra2_2 = arrange (nucleocytoplasmic_GO_ra2_2, gene_id)
# for ra3@2mm 
nucleocytoplasmic_GO_ra3_2 = filter(f2_nucleocytoplasmic, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")
nucleocytoplasmic_GO_ra3_2 = arrange (nucleocytoplasmic_GO_ra3_2, gene_id)
# ra1+ra2@1mm
nucleocytoplasmic_ra1_ra2_1mm = merge (nucleocytoplasmic_GO_ra1_1, nucleocytoplasmic_GO_ra2_1, by = "gene_id")
# ra1+ra3@1mm
nucleocytoplasmic_ra1_ra3_1mm = merge (nucleocytoplasmic_GO_ra1_1, nucleocytoplasmic_GO_ra3_1, by = "gene_id")
# ra1+ra2@2mm
nucleocytoplasmic_ra1_ra2_2mm = merge (nucleocytoplasmic_GO_ra1_2, nucleocytoplasmic_GO_ra2_2, by = "gene_id")
# ra1+ra3@2mm
nucleocytoplasmic_ra1_ra3_2mm = merge (nucleocytoplasmic_GO_ra1_2, nucleocytoplasmic_GO_ra3_2, by = "gene_id")
# mitochondrial_transport
mitochondrial_transport = read.csv("mitochondrial_transport.csv", header = TRUE)
mitochondrial_transport = tbl_df(mitochondrial_transport)
mitochondrial_transport_df = data.frame(mitochondrial_transport$maize.gene.id)
mitochondrial_transport_1 = arrange(mitochondrial_transport_df, mitochondrial_transport.maize.gene.id)
f2_mitochondrial_transport = merge(f2_df_1, mitochondrial_transport_1, by.x = "gene_id", by.y = "mitochondrial_transport.maize.gene.id")
# for ra1@1mm 
mitochondrial_transport_GO_ra1_1 = filter(f2_mitochondrial_transport, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
mitochondrial_transport_GO_ra1_1 = arrange (mitochondrial_transport_GO_ra1_1, gene_id)
# for ra2@1mm 
mitochondrial_transport_GO_ra2_1 = filter(f2_mitochondrial_transport, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
mitochondrial_transport_GO_ra2_1 = arrange (mitochondrial_transport_GO_ra2_1, gene_id)
# for ra3@1mm 
mitochondrial_transport_GO_ra3_1 = filter(f2_mitochondrial_transport, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
mitochondrial_transport_GO_ra3_1 = arrange (mitochondrial_transport_GO_ra3_1, gene_id)
# for ra1@2mm 
mitochondrial_transport_GO_ra1_2 = filter(f2_mitochondrial_transport, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
mitochondrial_transport_GO_ra1_2 = arrange (mitochondrial_transport_GO_ra1_2, gene_id)
# for ra2@2mm 
mitochondrial_transport_GO_ra2_2 = filter(f2_mitochondrial_transport, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
mitochondrial_transport_GO_ra2_2 = arrange (mitochondrial_transport_GO_ra2_2, gene_id)
# for ra3@2mm 
mitochondrial_transport_GO_ra3_2 = filter(f2_mitochondrial_transport, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")
mitochondrial_transport_GO_ra3_2 = arrange (mitochondrial_transport_GO_ra3_2, gene_id)
# ra1+ra2@1mm
mitochondrial_transport_ra1_ra2_1mm = merge (mitochondrial_transport_GO_ra1_1, mitochondrial_transport_GO_ra2_1, by = "gene_id")
# ra1+ra3@1mm
mitochondrial_transport_ra1_ra3_1mm = merge (mitochondrial_transport_GO_ra1_1, mitochondrial_transport_GO_ra3_1, by = "gene_id")
# ra1+ra2@2mm
mitochondrial_transport_ra1_ra2_2mm = merge (mitochondrial_transport_GO_ra1_2, mitochondrial_transport_GO_ra2_2, by = "gene_id")
# ra1+ra3@2mm
mitochondrial_transport_ra1_ra3_2mm = merge (mitochondrial_transport_GO_ra1_2, mitochondrial_transport_GO_ra3_2, by = "gene_id")
# ion_trans
ion_trans = read.csv("ion_trans.csv", header = TRUE)
ion_trans = tbl_df(ion_trans)
ion_trans_df = data.frame(ion_trans$maize.gene.id)
ion_trans_1 = arrange(ion_trans_df, ion_trans.maize.gene.id)
f2_ion_trans = merge(f2_df_1, ion_trans_1, by.x = "gene_id", by.y = "ion_trans.maize.gene.id")
# for ra1@1mm 
ion_trans_GO_ra1_1 = filter(f2_ion_trans, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
ion_trans_GO_ra1_1 = arrange (ion_trans_GO_ra1_1, gene_id)
# for ra2@1mm 
ion_trans_GO_ra2_1 = filter(f2_ion_trans, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
ion_trans_GO_ra2_1 = arrange (ion_trans_GO_ra2_1, gene_id)
# for ra3@1mm 
ion_trans_GO_ra3_1 = filter(f2_ion_trans, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
ion_trans_GO_ra3_1 = arrange (ion_trans_GO_ra3_1, gene_id)
# for ra1@2mm 
ion_trans_GO_ra1_2 = filter(f2_ion_trans, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
ion_trans_GO_ra1_2 = arrange (ion_trans_GO_ra1_2, gene_id)
# for ra2@2mm 
ion_trans_GO_ra2_2 = filter(f2_ion_trans, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
ion_trans_GO_ra2_2 = arrange (ion_trans_GO_ra2_2, gene_id)
# for ra3@2mm 
ion_trans_GO_ra3_2 = filter(f2_ion_trans, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")
ion_trans_GO_ra3_2 = arrange (ion_trans_GO_ra3_2, gene_id)
# ra1+ra2@1mm
ion_trans_ra1_ra2_1mm = merge (ion_trans_GO_ra1_1, ion_trans_GO_ra2_1, by = "gene_id")
# ra1+ra3@1mm
ion_trans_ra1_ra3_1mm = merge (ion_trans_GO_ra1_1, ion_trans_GO_ra3_1, by = "gene_id")
# ra1+ra2@2mm
ion_trans_ra1_ra2_2mm = merge (ion_trans_GO_ra1_2, ion_trans_GO_ra2_2, by = "gene_id")
# ra1+ra3@2mm
ion_trans_ra1_ra3_2mm = merge (ion_trans_GO_ra1_2, ion_trans_GO_ra3_2, by = "gene_id")
# ATP_synthesis_coupled
ATP_synthesis_coupled = read.csv("ATP_synthesis_coupled.csv", header = TRUE)
ATP_synthesis_coupled = tbl_df(ATP_synthesis_coupled)
ATP_synthesis_coupled_df = data.frame(ATP_synthesis_coupled$maize.gene.id)
ATP_synthesis_coupled_1 = arrange(ATP_synthesis_coupled_df, ATP_synthesis_coupled.maize.gene.id)
f2_ATP_synthesis_coupled = merge(f2_df_1, ATP_synthesis_coupled_1, by.x = "gene_id", by.y = "ATP_synthesis_coupled.maize.gene.id")
# for ra1@1mm 
ATP_synthesis_coupled_GO_ra1_1 = filter(f2_ATP_synthesis_coupled, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
ATP_synthesis_coupled_GO_ra1_1 = arrange (ATP_synthesis_coupled_GO_ra1_1, gene_id)
# for ra2@1mm 
ATP_synthesis_coupled_GO_ra2_1 = filter(f2_ATP_synthesis_coupled, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
ATP_synthesis_coupled_GO_ra2_1 = arrange (ATP_synthesis_coupled_GO_ra2_1, gene_id)
# for ra3@1mm 
ATP_synthesis_coupled_GO_ra3_1 = filter(f2_ATP_synthesis_coupled, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
ATP_synthesis_coupled_GO_ra3_1 = arrange (ATP_synthesis_coupled_GO_ra3_1, gene_id)
# for ra1@2mm 
ATP_synthesis_coupled_GO_ra1_2 = filter(f2_ATP_synthesis_coupled, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
ATP_synthesis_coupled_GO_ra1_2 = arrange (ATP_synthesis_coupled_GO_ra1_2, gene_id)
# for ra2@2mm 
ATP_synthesis_coupled_GO_ra2_2 = filter(f2_ATP_synthesis_coupled, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
ATP_synthesis_coupled_GO_ra2_2 = arrange (ATP_synthesis_coupled_GO_ra2_2, gene_id)
# for ra3@2mm 
ATP_synthesis_coupled_GO_ra3_2 = filter(f2_ATP_synthesis_coupled, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")
ATP_synthesis_coupled_GO_ra3_2 = arrange (ATP_synthesis_coupled_GO_ra3_2, gene_id)
# ra1+ra2@1mm
ATP_synthesis_coupled_ra1_ra2_1mm = merge (ATP_synthesis_coupled_GO_ra1_1, ATP_synthesis_coupled_GO_ra2_1, by = "gene_id")
# ra1+ra3@1mm
ATP_synthesis_coupled_ra1_ra3_1mm = merge (ATP_synthesis_coupled_GO_ra1_1, ATP_synthesis_coupled_GO_ra3_1, by = "gene_id")
# ra1+ra2@2mm
ATP_synthesis_coupled_ra1_ra2_2mm = merge (ATP_synthesis_coupled_GO_ra1_2, ATP_synthesis_coupled_GO_ra2_2, by = "gene_id")
# ra1+ra3@2mm
ATP_synthesis_coupled_ra1_ra3_2mm = merge (ATP_synthesis_coupled_GO_ra1_2, ATP_synthesis_coupled_GO_ra3_2, by = "gene_id")

# ra1+ra2+ra3@1mm
nucleosome_ra1_ra2_ra3_1mm = merge(nucleosome_ra1_ra2_1mm, nucleosome_ra1_ra3_1mm, by = "gene_id")
chromatin_ra1_ra2_ra3_1mm = merge(chromatin_ra1_ra2_1mm, chromatin_ra1_ra3_1mm, by = "gene_id")
vesicle_ra1_ra2_ra3_1mm = merge(vesicle_ra1_ra2_1mm, vesicle_ra1_ra3_1mm, by = "gene_id")
temperature_ra1_ra2_ra3_1mm = merge(temperature_ra1_ra2_1mm, temperature_ra1_ra3_1mm, by = "gene_id")
stress_ra1_ra2_ra3_1mm = merge(stress_ra1_ra2_1mm, stress_ra1_ra3_1mm, by = "gene_id")
transcription_ra1_ra2_ra3_1mm = merge(transcription_ra1_ra2_1mm, transcription_ra1_ra3_1mm, by = "gene_id")
nitrogen_ra1_ra2_ra3_1mm = merge(nitrogen_ra1_ra2_1mm, nitrogen_ra1_ra3_1mm, by = "gene_id")
RNA_synthetase_ra1_ra2_ra3_1mm = merge(RNA_synthetase_ra1_ra2_1mm, RNA_synthetase_ra1_ra3_1mm, by = "gene_id")
G_protein_ra1_ra2_ra3_1mm = merge(G_protein_ra1_ra2_1mm, G_protein_ra1_ra3_1mm, by = "gene_id")
light_ra1_ra2_ra3_1mm = merge(light_ra1_ra2_1mm, light_ra1_ra3_1mm, by = "gene_id")
redox_ra1_ra2_ra3_1mm = merge(redox_ra1_ra2_1mm, redox_ra1_ra3_1mm, by = "gene_id")
monosaccharide_ra1_ra2_ra3_1mm = merge(monosaccharide_ra1_ra2_1mm, monosaccharide_ra1_ra3_1mm, by = "gene_id")
cell_wall_ra1_ra2_ra3_1mm = merge(cell_wall_ra1_ra2_1mm, cell_wall_ra1_ra3_1mm, by = "gene_id")
sucrose_ra1_ra2_ra3_1mm = merge(sucrose_ra1_ra2_1mm, sucrose_ra1_ra3_1mm, by = "gene_id")
trehalose_ra1_ra2_ra3_1mm = merge(trehalose_ra1_ra2_1mm, trehalose_ra1_ra3_1mm, by = "gene_id")
intracellular_protein_ra1_ra2_ra3_1mm = merge(intracellular_protein_ra1_ra2_1mm, intracellular_protein_ra1_ra3_1mm, by = "gene_id")
small_GTPase_ra1_ra2_ra3_1mm = merge(small_GTPase_ra1_ra2_1mm, small_GTPase_ra1_ra3_1mm, by = "gene_id")
intracellular_signaling_ra1_ra2_ra3_1mm = merge(intracellular_signaling_ra1_ra2_1mm, intracellular_signaling_ra1_ra3_1mm, by = "gene_id")
proteasome_ra1_ra2_ra3_1mm = merge(proteasome_ra1_ra2_1mm, proteasome_ra1_ra3_1mm, by = "gene_id")
nucleocytoplasmic_ra1_ra2_ra3_1mm = merge(nucleocytoplasmic_ra1_ra2_1mm, nucleocytoplasmic_ra1_ra3_1mm, by = "gene_id")
mitochondrial_transport_ra1_ra2_ra3_1mm = merge(mitochondrial_transport_ra1_ra2_1mm, mitochondrial_transport_ra1_ra3_1mm, by = "gene_id")
ion_trans_ra1_ra2_ra3_1mm = merge(ion_trans_ra1_ra2_1mm, ion_trans_ra1_ra3_1mm, by = "gene_id")
ATP_synthesis_coupled_ra1_ra2_ra3_1mm = merge(ATP_synthesis_coupled_ra1_ra2_1mm, ATP_synthesis_coupled_ra1_ra3_1mm, by = "gene_id")

# ra1+ra2+ra3@2mm
nucleosome_ra1_ra2_ra3_2mm = merge(nucleosome_ra1_ra2_2mm, nucleosome_ra1_ra3_2mm, by = "gene_id")
chromatin_ra1_ra2_ra3_2mm = merge(chromatin_ra1_ra2_2mm, chromatin_ra1_ra3_2mm, by = "gene_id")
vesicle_ra1_ra2_ra3_2mm = merge(vesicle_ra1_ra2_2mm, vesicle_ra1_ra3_2mm, by = "gene_id")
temperature_ra1_ra2_ra3_2mm = merge(temperature_ra1_ra2_2mm, temperature_ra1_ra3_2mm, by = "gene_id")
stress_ra1_ra2_ra3_2mm = merge(stress_ra1_ra2_2mm, stress_ra1_ra3_2mm, by = "gene_id")
transcription_ra1_ra2_ra3_2mm = merge(transcription_ra1_ra2_2mm, transcription_ra1_ra3_2mm, by = "gene_id")
nitrogen_ra1_ra2_ra3_2mm = merge(nitrogen_ra1_ra2_2mm, nitrogen_ra1_ra3_2mm, by = "gene_id")
RNA_synthetase_ra1_ra2_ra3_2mm = merge(RNA_synthetase_ra1_ra2_2mm, RNA_synthetase_ra1_ra3_2mm, by = "gene_id")
G_protein_ra1_ra2_ra3_2mm = merge(G_protein_ra1_ra2_2mm, G_protein_ra1_ra3_2mm, by = "gene_id")
light_ra1_ra2_ra3_2mm = merge(light_ra1_ra2_2mm, light_ra1_ra3_2mm, by = "gene_id")
redox_ra1_ra2_ra3_2mm = merge(redox_ra1_ra2_2mm, redox_ra1_ra3_2mm, by = "gene_id")
monosaccharide_ra1_ra2_ra3_2mm = merge(monosaccharide_ra1_ra2_2mm, monosaccharide_ra1_ra3_2mm, by = "gene_id")
cell_wall_ra1_ra2_ra3_2mm = merge(cell_wall_ra1_ra2_2mm, cell_wall_ra1_ra3_2mm, by = "gene_id")
sucrose_ra1_ra2_ra3_2mm = merge(sucrose_ra1_ra2_2mm, sucrose_ra1_ra3_2mm, by = "gene_id")
trehalose_ra1_ra2_ra3_2mm = merge(trehalose_ra1_ra2_2mm, trehalose_ra1_ra3_2mm, by = "gene_id")
intracellular_protein_ra1_ra2_ra3_2mm = merge(intracellular_protein_ra1_ra2_2mm, intracellular_protein_ra1_ra3_2mm, by = "gene_id")
small_GTPase_ra1_ra2_ra3_2mm = merge(small_GTPase_ra1_ra2_2mm, small_GTPase_ra1_ra3_2mm, by = "gene_id")
intracellular_signaling_ra1_ra2_ra3_2mm = merge(intracellular_signaling_ra1_ra2_2mm, intracellular_signaling_ra1_ra3_2mm, by = "gene_id")
proteasome_ra1_ra2_ra3_2mm = merge(proteasome_ra1_ra2_2mm, proteasome_ra1_ra3_2mm, by = "gene_id")
nucleocytoplasmic_ra1_ra2_ra3_2mm = merge(nucleocytoplasmic_ra1_ra2_2mm, nucleocytoplasmic_ra1_ra3_2mm, by = "gene_id")
mitochondrial_transport_ra1_ra2_ra3_2mm = merge(mitochondrial_transport_ra1_ra2_2mm, mitochondrial_transport_ra1_ra3_2mm, by = "gene_id")
ion_trans_ra1_ra2_ra3_2mm = merge(ion_trans_ra1_ra2_2mm, ion_trans_ra1_ra3_2mm, by = "gene_id")
ATP_synthesis_coupled_ra1_ra2_ra3_2mm = merge(ATP_synthesis_coupled_ra1_ra2_2mm, ATP_synthesis_coupled_ra1_ra3_2mm, by = "gene_id")

###########################################################
## Calculating the p value for each overlapping region
# column 1
n1_1_1c = (mean(nucleosome_ra1_ra2_ra3_1mm$qval.x.x) + mean(nucleosome_ra1_ra2_ra3_1mm$qval.y.x)+mean(nucleosome_ra1_ra2_ra3_1mm$qval.y.y))/3
n1_2_1c = (mean(nucleosome_ra1_ra2_ra3_2mm$qval.x.x) + mean(nucleosome_ra1_ra2_ra3_2mm$qval.y.x)+mean(nucleosome_ra1_ra2_ra3_2mm$qval.y.y))/3
c2_1_1c = (mean(chromatin_ra1_ra2_ra3_1mm$qval.x.x) + mean(chromatin_ra1_ra2_ra3_1mm$qval.y.x)+mean(chromatin_ra1_ra2_ra3_1mm$qval.y.y))/3
c2_2_1c = (mean(chromatin_ra1_ra2_ra3_2mm$qval.x.x) + mean(chromatin_ra1_ra2_ra3_2mm$qval.y.x)+mean(chromatin_ra1_ra2_ra3_2mm$qval.y.y))/3
v3_1_1c = (mean(vesicle_ra1_ra2_ra3_1mm$qval.x.x) + mean(vesicle_ra1_ra2_ra3_1mm$qval.y.x)+mean(vesicle_ra1_ra2_ra3_1mm$qval.y.y))/3
v3_2_1c = (mean(vesicle_ra1_ra2_ra3_2mm$qval.x.x) + mean(vesicle_ra1_ra2_ra3_2mm$qval.y.x)+mean(vesicle_ra1_ra2_ra3_2mm$qval.y.y))/3
t4_1_1c = (mean(temperature_ra1_ra2_ra3_1mm$qval.x.x) + mean(temperature_ra1_ra2_ra3_1mm$qval.y.x)+mean(temperature_ra1_ra2_ra3_1mm$qval.y.y))/3
t4_2_1c = (mean(temperature_ra1_ra2_ra3_2mm$qval.x.x) + mean(temperature_ra1_ra2_ra3_2mm$qval.y.x)+mean(temperature_ra1_ra2_ra3_2mm$qval.y.y))/3
s5_1_1c = (mean(stress_ra1_ra2_ra3_1mm$qval.x.x) + mean(stress_ra1_ra2_ra3_1mm$qval.y.x)+mean(stress_ra1_ra2_ra3_1mm$qval.y.y))/3
s5_2_1c = (mean(stress_ra1_ra2_ra3_2mm$qval.x.x) + mean(stress_ra1_ra2_ra3_2mm$qval.y.x)+mean(stress_ra1_ra2_ra3_2mm$qval.y.y))/3
t6_1_1c = (mean(transcription_ra1_ra2_ra3_1mm$qval.x.x) + mean(transcription_ra1_ra2_ra3_1mm$qval.y.x)+mean(transcription_ra1_ra2_ra3_1mm$qval.y.y))/3
t6_2_1c = (mean(transcription_ra1_ra2_ra3_2mm$qval.x.x) + mean(transcription_ra1_ra2_ra3_2mm$qval.y.x)+mean(transcription_ra1_ra2_ra3_2mm$qval.y.y))/3
n7_1_1c = (mean(nitrogen_ra1_ra2_ra3_1mm$qval.x.x) + mean(nitrogen_ra1_ra2_ra3_1mm$qval.y.x)+mean(nitrogen_ra1_ra2_ra3_1mm$qval.y.y))/3
n7_2_1c = (mean(nitrogen_ra1_ra2_ra3_2mm$qval.x.x) + mean(nitrogen_ra1_ra2_ra3_2mm$qval.y.x)+mean(nitrogen_ra1_ra2_ra3_2mm$qval.y.y))/3
R8_1_1c = (mean(RNA_synthetase_ra1_ra2_ra3_1mm$qval.x.x) + mean(RNA_synthetase_ra1_ra2_ra3_1mm$qval.y.x)+mean(RNA_synthetase_ra1_ra2_ra3_1mm$qval.y.y))/3
R8_2_1c = (mean(RNA_synthetase_ra1_ra2_ra3_2mm$qval.x.x) + mean(RNA_synthetase_ra1_ra2_ra3_2mm$qval.y.x)+mean(RNA_synthetase_ra1_ra2_ra3_2mm$qval.y.y))/3
G9_1_1c = (mean(G_protein_ra1_ra2_ra3_1mm$qval.x.x) + mean(G_protein_ra1_ra2_ra3_1mm$qval.y.x)+mean(G_protein_ra1_ra2_ra3_1mm$qval.y.y))/3
G9_2_1c = (mean(G_protein_ra1_ra2_ra3_2mm$qval.x.x) + mean(G_protein_ra1_ra2_ra3_2mm$qval.y.x)+mean(G_protein_ra1_ra2_ra3_2mm$qval.y.y))/3
l10_1_1c = (mean(light_ra1_ra2_ra3_1mm$qval.x.x) + mean(light_ra1_ra2_ra3_1mm$qval.y.x)+mean(light_ra1_ra2_ra3_1mm$qval.y.y))/3
l10_2_1c = (mean(light_ra1_ra2_ra3_2mm$qval.x.x) + mean(light_ra1_ra2_ra3_2mm$qval.y.x)+mean(light_ra1_ra2_ra3_2mm$qval.y.y))/3
r11_1_1c = (mean(redox_ra1_ra2_ra3_1mm$qval.x.x) + mean(redox_ra1_ra2_ra3_1mm$qval.y.x)+mean(redox_ra1_ra2_ra3_1mm$qval.y.y))/3
r11_2_1c = (mean(redox_ra1_ra2_ra3_2mm$qval.x.x) + mean(redox_ra1_ra2_ra3_2mm$qval.y.x)+mean(redox_ra1_ra2_ra3_2mm$qval.y.y))/3
m12_1_1c = (mean(monosaccharide_ra1_ra2_ra3_1mm$qval.x.x) + mean(monosaccharide_ra1_ra2_ra3_1mm$qval.y.x)+mean(monosaccharide_ra1_ra2_ra3_1mm$qval.y.y))/3
m12_2_1c = (mean(monosaccharide_ra1_ra2_ra3_2mm$qval.x.x) + mean(monosaccharide_ra1_ra2_ra3_2mm$qval.y.x)+mean(monosaccharide_ra1_ra2_ra3_2mm$qval.y.y))/3
c13_1_1c = (mean(cell_wall_ra1_ra2_ra3_1mm$qval.x.x) + mean(cell_wall_ra1_ra2_ra3_1mm$qval.y.x)+mean(cell_wall_ra1_ra2_ra3_1mm$qval.y.y))/3
c13_2_1c = (mean(cell_wall_ra1_ra2_ra3_2mm$qval.x.x) + mean(cell_wall_ra1_ra2_ra3_2mm$qval.y.x)+mean(cell_wall_ra1_ra2_ra3_2mm$qval.y.y))/3
s14_1_1c = (mean(sucrose_ra1_ra2_ra3_1mm$qval.x.x) + mean(sucrose_ra1_ra2_ra3_1mm$qval.y.x)+mean(sucrose_ra1_ra2_ra3_1mm$qval.y.y))/3
s14_2_1c = (mean(sucrose_ra1_ra2_ra3_2mm$qval.x.x) + mean(sucrose_ra1_ra2_ra3_2mm$qval.y.x)+mean(sucrose_ra1_ra2_ra3_2mm$qval.y.y))/3
t15_1_1c = (mean(trehalose_ra1_ra2_ra3_1mm$qval.x.x) + mean(trehalose_ra1_ra2_ra3_1mm$qval.y.x)+mean(trehalose_ra1_ra2_ra3_1mm$qval.y.y))/3
t15_2_1c = (mean(trehalose_ra1_ra2_ra3_2mm$qval.x.x) + mean(trehalose_ra1_ra2_ra3_2mm$qval.y.x)+mean(trehalose_ra1_ra2_ra3_2mm$qval.y.y))/3
i16_1_1c = (mean(intracellular_protein_ra1_ra2_ra3_1mm$qval.x.x) + mean(intracellular_protein_ra1_ra2_ra3_1mm$qval.y.x)+mean(intracellular_protein_ra1_ra2_ra3_1mm$qval.y.y))/3
i16_2_1c = (mean(intracellular_protein_ra1_ra2_ra3_2mm$qval.x.x) + mean(intracellular_protein_ra1_ra2_ra3_2mm$qval.y.x)+mean(intracellular_protein_ra1_ra2_ra3_2mm$qval.y.y))/3
s17_1_1c = (mean(small_GTPase_ra1_ra2_ra3_1mm$qval.x.x) + mean(small_GTPase_ra1_ra2_ra3_1mm$qval.y.x)+mean(small_GTPase_ra1_ra2_ra3_1mm$qval.y.y))/3
s17_2_1c = (mean(small_GTPase_ra1_ra2_ra3_2mm$qval.x.x) + mean(small_GTPase_ra1_ra2_ra3_2mm$qval.y.x)+mean(small_GTPase_ra1_ra2_ra3_2mm$qval.y.y))/3
i18_1_1c = (mean(intracellular_signaling_ra1_ra2_ra3_1mm$qval.x.x) + mean(intracellular_signaling_ra1_ra2_ra3_1mm$qval.y.x)+mean(intracellular_signaling_ra1_ra2_ra3_1mm$qval.y.y))/3
i18_2_1c = (mean(intracellular_signaling_ra1_ra2_ra3_2mm$qval.x.x) + mean(intracellular_signaling_ra1_ra2_ra3_2mm$qval.y.x)+mean(intracellular_signaling_ra1_ra2_ra3_2mm$qval.y.y))/3
p19_1_1c = (mean(proteasome_ra1_ra2_ra3_1mm$qval.x.x) + mean(proteasome_ra1_ra2_ra3_1mm$qval.y.x)+mean(proteasome_ra1_ra2_ra3_1mm$qval.y.y))/3
p19_2_1c = (mean(proteasome_ra1_ra2_ra3_2mm$qval.x.x) + mean(proteasome_ra1_ra2_ra3_2mm$qval.y.x)+mean(proteasome_ra1_ra2_ra3_2mm$qval.y.y))/3
n20_1_1c = (mean(nucleocytoplasmic_ra1_ra2_ra3_1mm$qval.x.x) + mean(nucleocytoplasmic_ra1_ra2_ra3_1mm$qval.y.x)+mean(nucleocytoplasmic_ra1_ra2_ra3_1mm$qval.y.y))/3
n20_2_1c = (mean(nucleocytoplasmic_ra1_ra2_ra3_2mm$qval.x.x) + mean(nucleocytoplasmic_ra1_ra2_ra3_2mm$qval.y.x)+mean(nucleocytoplasmic_ra1_ra2_ra3_2mm$qval.y.y))/3
m21_1_1c = (mean(mitochondrial_transport_ra1_ra2_ra3_1mm$qval.x.x) + mean(mitochondrial_transport_ra1_ra2_ra3_1mm$qval.y.x)+mean(mitochondrial_transport_ra1_ra2_ra3_1mm$qval.y.y))/3
m21_2_1c = (mean(mitochondrial_transport_ra1_ra2_ra3_2mm$qval.x.x) + mean(mitochondrial_transport_ra1_ra2_ra3_2mm$qval.y.x)+mean(mitochondrial_transport_ra1_ra2_ra3_2mm$qval.y.y))/3
i22_1_1c = (mean(ion_trans_ra1_ra2_ra3_1mm$qval.x.x) + mean(ion_trans_ra1_ra2_ra3_1mm$qval.y.x)+mean(ion_trans_ra1_ra2_ra3_1mm$qval.y.y))/3
i22_2_1c = (mean(ion_trans_ra1_ra2_ra3_2mm$qval.x.x) + mean(ion_trans_ra1_ra2_ra3_2mm$qval.y.x)+mean(ion_trans_ra1_ra2_ra3_2mm$qval.y.y))/3
A23_1_1c = (mean(ATP_synthesis_coupled_ra1_ra2_ra3_1mm$qval.x.x) + mean(ATP_synthesis_coupled_ra1_ra2_ra3_1mm$qval.y.x)+ mean(ATP_synthesis_coupled_ra1_ra2_ra3_1mm$qval.y.y))/3
A23_2_1c = (mean(ATP_synthesis_coupled_ra1_ra2_ra3_2mm$qval.x.x) + mean(ATP_synthesis_coupled_ra1_ra2_ra3_2mm$qval.y.x)+mean(ATP_synthesis_coupled_ra1_ra2_ra3_2mm$qval.y.y))/3
# column 2
n1_1_2c = (mean(nucleosome_ra1_ra2_1mm$qval.x) + mean(nucleosome_ra1_ra2_1mm$qval.y))/2
n1_2_2c = (mean(nucleosome_ra1_ra2_2mm$qval.x) + mean(nucleosome_ra1_ra2_2mm$qval.y))/2
c2_1_2c = (mean(chromatin_ra1_ra2_1mm$qval.x) + mean(chromatin_ra1_ra2_1mm$qval.y))/2
c2_2_2c = (mean(chromatin_ra1_ra2_2mm$qval.x) + mean(chromatin_ra1_ra2_2mm$qval.y))/2
v3_1_2c = (mean(vesicle_ra1_ra2_1mm$qval.x) + mean(vesicle_ra1_ra2_1mm$qval.y))/2
v3_2_2c = (mean(vesicle_ra1_ra2_2mm$qval.x) + mean(vesicle_ra1_ra2_2mm$qval.y))/2
t4_1_2c = (mean(temperature_ra1_ra2_1mm$qval.x) + mean(temperature_ra1_ra2_1mm$qval.y))/2
t4_2_2c = (mean(temperature_ra1_ra2_2mm$qval.x) + mean(temperature_ra1_ra2_2mm$qval.y))/2
s5_1_2c = (mean(stress_ra1_ra2_1mm$qval.x) + mean(stress_ra1_ra2_1mm$qval.y))/2
s5_2_2c = (mean(stress_ra1_ra2_2mm$qval.x) + mean(stress_ra1_ra2_2mm$qval.y))/2
t6_1_2c = (mean(transcription_ra1_ra2_1mm$qval.x) + mean(transcription_ra1_ra2_1mm$qval.y))/2
t6_2_2c = (mean(transcription_ra1_ra2_2mm$qval.x) + mean(transcription_ra1_ra2_2mm$qval.y))/2
n7_1_2c = (mean(nitrogen_ra1_ra2_1mm$qval.x) + mean(nitrogen_ra1_ra2_1mm$qval.y))/2
n7_2_2c = (mean(nitrogen_ra1_ra2_2mm$qval.x) + mean(nitrogen_ra1_ra2_2mm$qval.y))/2
R8_1_2c = (mean(RNA_synthetase_ra1_ra2_1mm$qval.x) + mean(RNA_synthetase_ra1_ra2_1mm$qval.y))/2
R8_2_2c = (mean(RNA_synthetase_ra1_ra2_2mm$qval.x) + mean(RNA_synthetase_ra1_ra2_2mm$qval.y))/2
G9_1_2c = (mean(G_protein_ra1_ra2_1mm$qval.x) + mean(G_protein_ra1_ra2_1mm$qval.y))/2
G9_2_2c = (mean(G_protein_ra1_ra2_2mm$qval.x) + mean(G_protein_ra1_ra2_2mm$qval.y))/2
l10_1_2c = (mean(light_ra1_ra2_1mm$qval.x) + mean(light_ra1_ra2_1mm$qval.y))/2
l10_2_2c = (mean(light_ra1_ra2_2mm$qval.x) + mean(light_ra1_ra2_2mm$qval.y))/2
r11_1_2c = (mean(redox_ra1_ra2_1mm$qval.x) + mean(redox_ra1_ra2_1mm$qval.y))/2
r11_2_2c = (mean(redox_ra1_ra2_2mm$qval.x) + mean(redox_ra1_ra2_2mm$qval.y))/2
m12_1_2c = (mean(monosaccharide_ra1_ra2_1mm$qval.x) + mean(monosaccharide_ra1_ra2_1mm$qval.y))/2
m12_2_2c = (mean(monosaccharide_ra1_ra2_2mm$qval.x) + mean(monosaccharide_ra1_ra2_2mm$qval.y))/2
c13_1_2c = (mean(cell_wall_ra1_ra2_1mm$qval.x) + mean(cell_wall_ra1_ra2_1mm$qval.y))/2
c13_2_2c = (mean(cell_wall_ra1_ra2_2mm$qval.x) + mean(cell_wall_ra1_ra2_2mm$qval.y))/2
s14_1_2c = (mean(sucrose_ra1_ra2_1mm$qval.x) + mean(sucrose_ra1_ra2_1mm$qval.y))/2
s14_2_2c = (mean(sucrose_ra1_ra2_2mm$qval.x) + mean(sucrose_ra1_ra2_2mm$qval.y))/2
t15_1_2c = (mean(trehalose_ra1_ra2_1mm$qval.x) + mean(trehalose_ra1_ra2_1mm$qval.y))/2
t15_2_2c = (mean(trehalose_ra1_ra2_2mm$qval.x) + mean(trehalose_ra1_ra2_2mm$qval.y))/2
i16_1_2c = (mean(intracellular_protein_ra1_ra2_1mm$qval.x) + mean(intracellular_protein_ra1_ra2_1mm$qval.y))/2
i16_2_2c = (mean(intracellular_protein_ra1_ra2_2mm$qval.x) + mean(intracellular_protein_ra1_ra2_2mm$qval.y))/2
s17_1_2c = (mean(small_GTPase_ra1_ra2_1mm$qval.x) + mean(small_GTPase_ra1_ra2_1mm$qval.y))/2
s17_2_2c = (mean(small_GTPase_ra1_ra2_2mm$qval.x) + mean(small_GTPase_ra1_ra2_2mm$qval.y))/2
i18_1_2c = (mean(intracellular_signaling_ra1_ra2_1mm$qval.x) + mean(intracellular_signaling_ra1_ra2_1mm$qval.y))/2
i18_2_2c = (mean(intracellular_signaling_ra1_ra2_2mm$qval.x) + mean(intracellular_signaling_ra1_ra2_2mm$qval.y))/2
p19_1_2c = (mean(proteasome_ra1_ra2_1mm$qval.x) + mean(proteasome_ra1_ra2_1mm$qval.y))/2
p19_2_2c = (mean(proteasome_ra1_ra2_2mm$qval.x) + mean(proteasome_ra1_ra2_2mm$qval.y))/2
n20_1_2c = (mean(nucleocytoplasmic_ra1_ra2_1mm$qval.x) + mean(nucleocytoplasmic_ra1_ra2_1mm$qval.y))/2
n20_2_2c = (mean(nucleocytoplasmic_ra1_ra2_2mm$qval.x) + mean(nucleocytoplasmic_ra1_ra2_2mm$qval.y))/2
m21_1_2c = (mean(mitochondrial_transport_ra1_ra2_1mm$qval.x) + mean(mitochondrial_transport_ra1_ra2_1mm$qval.y))/2
m21_2_2c = (mean(mitochondrial_transport_ra1_ra2_2mm$qval.x) + mean(mitochondrial_transport_ra1_ra2_2mm$qval.y))/2
i22_1_2c = (mean(ion_trans_ra1_ra2_1mm$qval.x) + mean(ion_trans_ra1_ra2_1mm$qval.y))/2
i22_2_2c = (mean(ion_trans_ra1_ra2_2mm$qval.x) + mean(ion_trans_ra1_ra2_2mm$qval.y))/2
A23_1_2c = (mean(ATP_synthesis_coupled_ra1_ra2_1mm$qval.x) + mean(ATP_synthesis_coupled_ra1_ra2_1mm$qval.y))/2
A23_2_2c = (mean(ATP_synthesis_coupled_ra1_ra2_2mm$qval.x) + mean(ATP_synthesis_coupled_ra1_ra2_2mm$qval.y))/2
# column 3
n1_1_3c = (mean(nucleosome_ra2_ra3_1mm$qval.x) + mean(nucleosome_ra2_ra3_1mm$qval.y))/2
n1_2_3c = (mean(nucleosome_ra2_ra3_2mm$qval.x) + mean(nucleosome_ra2_ra3_2mm$qval.y))/2
c2_1_3c = (mean(chromatin_ra2_ra3_1mm$qval.x) + mean(chromatin_ra2_ra3_1mm$qval.y))/2
c2_2_3c = (mean(chromatin_ra2_ra3_2mm$qval.x) + mean(chromatin_ra2_ra3_2mm$qval.y))/2
v3_1_3c = (mean(vesicle_ra2_ra3_1mm$qval.x) + mean(vesicle_ra2_ra3_1mm$qval.y))/2
v3_2_3c = (mean(vesicle_ra2_ra3_2mm$qval.x) + mean(vesicle_ra2_ra3_2mm$qval.y))/2
t4_1_3c = (mean(temperature_ra2_ra3_1mm$qval.x) + mean(temperature_ra2_ra3_1mm$qval.y))/2
t4_2_3c = (mean(temperature_ra2_ra3_2mm$qval.x) + mean(temperature_ra2_ra3_2mm$qval.y))/2
s5_1_3c = (mean(stress_ra2_ra3_1mm$qval.x) + mean(stress_ra2_ra3_1mm$qval.y))/2
s5_2_3c = (mean(stress_ra2_ra3_2mm$qval.x) + mean(stress_ra2_ra3_2mm$qval.y))/2
t6_1_3c = (mean(transcription_ra2_ra3_1mm$qval.x) + mean(transcription_ra2_ra3_1mm$qval.y))/2
t6_2_3c = (mean(transcription_ra2_ra3_2mm$qval.x) + mean(transcription_ra2_ra3_2mm$qval.y))/2
n7_1_3c = (mean(nitrogen_ra2_ra3_1mm$qval.x) + mean(nitrogen_ra2_ra3_1mm$qval.y))/2
n7_2_3c = (mean(nitrogen_ra2_ra3_2mm$qval.x) + mean(nitrogen_ra2_ra3_2mm$qval.y))/2
R8_1_3c = (mean(RNA_synthetase_ra2_ra3_1mm$qval.x) + mean(RNA_synthetase_ra2_ra3_1mm$qval.y))/2
R8_2_3c = (mean(RNA_synthetase_ra2_ra3_2mm$qval.x) + mean(RNA_synthetase_ra2_ra3_2mm$qval.y))/2
G9_1_3c = (mean(G_protein_ra2_ra3_1mm$qval.x) + mean(G_protein_ra2_ra3_1mm$qval.y))/2
G9_2_3c = (mean(G_protein_ra2_ra3_2mm$qval.x) + mean(G_protein_ra2_ra3_2mm$qval.y))/2
l10_1_3c = (mean(light_ra2_ra3_1mm$qval.x) + mean(light_ra2_ra3_1mm$qval.y))/2
l10_2_3c = (mean(light_ra2_ra3_2mm$qval.x) + mean(light_ra2_ra3_2mm$qval.y))/2
r11_1_3c = (mean(redox_ra2_ra3_1mm$qval.x) + mean(redox_ra2_ra3_1mm$qval.y))/2
r11_2_3c = (mean(redox_ra2_ra3_2mm$qval.x) + mean(redox_ra2_ra3_2mm$qval.y))/2
m12_1_3c = (mean(monosaccharide_ra2_ra3_1mm$qval.x) + mean(monosaccharide_ra2_ra3_1mm$qval.y))/2
m12_2_3c = (mean(monosaccharide_ra2_ra3_2mm$qval.x) + mean(monosaccharide_ra2_ra3_2mm$qval.y))/2
c13_1_3c = (mean(cell_wall_ra2_ra3_1mm$qval.x) + mean(cell_wall_ra2_ra3_1mm$qval.y))/2
c13_2_3c = (mean(cell_wall_ra2_ra3_2mm$qval.x) + mean(cell_wall_ra2_ra3_2mm$qval.y))/2
s14_1_3c = (mean(sucrose_ra2_ra3_1mm$qval.x) + mean(sucrose_ra2_ra3_1mm$qval.y))/2
s14_2_3c = (mean(sucrose_ra2_ra3_2mm$qval.x) + mean(sucrose_ra2_ra3_2mm$qval.y))/2
t15_1_3c = (mean(trehalose_ra2_ra3_1mm$qval.x) + mean(trehalose_ra2_ra3_1mm$qval.y))/2
t15_2_3c = (mean(trehalose_ra2_ra3_2mm$qval.x) + mean(trehalose_ra2_ra3_2mm$qval.y))/2
i16_1_3c = (mean(intracellular_protein_ra2_ra3_1mm$qval.x) + mean(intracellular_protein_ra2_ra3_1mm$qval.y))/2
i16_2_3c = (mean(intracellular_protein_ra2_ra3_2mm$qval.x) + mean(intracellular_protein_ra2_ra3_2mm$qval.y))/2
s17_1_3c = (mean(small_GTPase_ra2_ra3_1mm$qval.x) + mean(small_GTPase_ra2_ra3_1mm$qval.y))/2
s17_2_3c = (mean(small_GTPase_ra2_ra3_2mm$qval.x) + mean(small_GTPase_ra2_ra3_2mm$qval.y))/2
i18_1_3c = (mean(intracellular_signaling_ra2_ra3_1mm$qval.x) + mean(intracellular_signaling_ra2_ra3_1mm$qval.y))/2
i18_2_3c = (mean(intracellular_signaling_ra2_ra3_2mm$qval.x) + mean(intracellular_signaling_ra2_ra3_2mm$qval.y))/2
p19_1_3c = (mean(proteasome_ra2_ra3_1mm$qval.x) + mean(proteasome_ra2_ra3_1mm$qval.y))/2
p19_2_3c = (mean(proteasome_ra2_ra3_2mm$qval.x) + mean(proteasome_ra2_ra3_2mm$qval.y))/2
n20_1_3c = (mean(nucleocytoplasmic_ra2_ra3_1mm$qval.x) + mean(nucleocytoplasmic_ra2_ra3_1mm$qval.y))/2
n20_2_3c = (mean(nucleocytoplasmic_ra2_ra3_2mm$qval.x) + mean(nucleocytoplasmic_ra2_ra3_2mm$qval.y))/2
m21_1_3c = (mean(mitochondrial_transport_ra2_ra3_1mm$qval.x) + mean(mitochondrial_transport_ra2_ra3_1mm$qval.y))/2
m21_2_3c = (mean(mitochondrial_transport_ra2_ra3_2mm$qval.x) + mean(mitochondrial_transport_ra2_ra3_2mm$qval.y))/2
i22_1_3c = (mean(ion_trans_ra2_ra3_1mm$qval.x) + mean(ion_trans_ra2_ra3_1mm$qval.y))/2
i22_2_3c = (mean(ion_trans_ra2_ra3_2mm$qval.x) + mean(ion_trans_ra2_ra3_2mm$qval.y))/2
A23_1_3c = (mean(ATP_synthesis_coupled_ra2_ra3_1mm$qval.x) + mean(ATP_synthesis_coupled_ra2_ra3_1mm$qval.y))/2
A23_2_3c = (mean(ATP_synthesis_coupled_ra2_ra3_2mm$qval.x) + mean(ATP_synthesis_coupled_ra2_ra3_2mm$qval.y))/2
# column 4
n1_1_4c = mean(nucleosome_GO_ra1_1$qval)
n1_2_4c = mean(nucleosome_GO_ra1_2$qval)
c2_1_4c = mean(chromatin_GO_ra1_1$qval)
c2_2_4c = mean(chromatin_GO_ra1_2$qval)
v3_1_4c = mean(vesicle_GO_ra1_1$qval)
v3_2_4c = mean(vesicle_GO_ra1_2$qval)
t4_1_4c = mean(temperature_GO_ra1_1$qval)
t4_2_4c = mean(temperature_GO_ra1_2$qval)
s5_1_4c = mean(stress_GO_ra1_1$qval)
s5_2_4c = mean(stress_GO_ra1_2$qval)
t6_1_4c = mean(transcription_GO_ra1_1$qval)
t6_2_4c = mean(transcription_GO_ra1_2$qval)
n7_1_4c = mean(nitrogen_GO_ra1_1$qval)
n7_2_4c = mean(nitrogen_GO_ra1_2$qval)
R8_1_4c = mean(RNA_synthetase_GO_ra1_1$qval)
R8_2_4c = mean(RNA_synthetase_GO_ra1_2$qval)
G9_1_4c = mean(G_protein_GO_ra1_1$qval)
G9_2_4c = mean(G_protein_GO_ra1_2$qval)
l10_1_4c = mean(light_GO_ra1_1$qval)
l10_2_4c = mean(light_GO_ra1_2$qval)
r11_1_4c = mean(redox_GO_ra1_1$qval)
r11_2_4c = mean(redox_GO_ra1_2$qval)
m12_1_4c = mean(monosaccharide_GO_ra1_1$qval)
m12_2_4c = mean(monosaccharide_GO_ra1_2$qval)
c13_1_4c = mean(cell_wall_GO_ra1_1$qval)
c13_2_4c = mean(cell_wall_GO_ra1_2$qval)
s14_1_4c = mean(sucrose_GO_ra1_1$qval)
s14_2_4c = mean(sucrose_GO_ra1_2$qval)
t15_1_4c = mean(trehalose_GO_ra1_1$qval)
t15_2_4c = mean(trehalose_GO_ra1_2$qval)
i16_1_4c = mean(intracellular_protein_GO_ra1_1$qval)
i16_2_4c = mean(intracellular_protein_GO_ra1_2$qval)
s17_1_4c = mean(small_GTPase_GO_ra1_1$qval)
s17_2_4c = mean(small_GTPase_GO_ra1_2$qval)
i18_1_4c = mean(intracellular_signaling_GO_ra1_1$qval)
i18_2_4c = mean(intracellular_signaling_GO_ra1_2$qval)
p19_1_4c = mean(proteasome_GO_ra1_1$qval)
p19_2_4c = mean(proteasome_GO_ra1_2$qval)
n20_1_4c = mean(nucleocytoplasmic_GO_ra1_1$qval)
n20_2_4c = mean(nucleocytoplasmic_GO_ra1_2$qval)
m21_1_4c = mean(mitochondrial_transport_GO_ra1_1$qval)
m21_2_4c = mean(mitochondrial_transport_GO_ra1_2$qval)
i22_1_4c = mean(ion_trans_GO_ra1_1$qval)
i22_2_4c = mean(ion_trans_GO_ra1_2$qval)
A23_1_4c = mean(ATP_synthesis_coupled_GO_ra1_1$qval)
A23_2_4c = mean(ATP_synthesis_coupled_GO_ra1_2$qval)
# column 5
n1_1_5c = mean(nucleosome_GO_ra3_1$qval)
n1_2_5c = mean(nucleosome_GO_ra3_2$qval)
c2_1_5c = mean(chromatin_GO_ra3_1$qval)
c2_2_5c = mean(chromatin_GO_ra3_2$qval)
v3_1_5c = mean(vesicle_GO_ra3_1$qval)
v3_2_5c = mean(vesicle_GO_ra3_2$qval)
t4_1_5c = mean(temperature_GO_ra3_1$qval)
t4_2_5c = mean(temperature_GO_ra3_2$qval)
s5_1_5c = mean(stress_GO_ra3_1$qval)
s5_2_5c = mean(stress_GO_ra3_2$qval)
t6_1_5c = mean(transcription_GO_ra3_1$qval)
t6_2_5c = mean(transcription_GO_ra3_2$qval)
n7_1_5c = mean(nitrogen_GO_ra3_1$qval)
n7_2_5c = mean(nitrogen_GO_ra3_2$qval)
R8_1_5c = mean(RNA_synthetase_GO_ra3_1$qval)
R8_2_5c = mean(RNA_synthetase_GO_ra3_2$qval)
G9_1_5c = mean(G_protein_GO_ra3_1$qval)
G9_2_5c = mean(G_protein_GO_ra3_2$qval)
l10_1_5c = mean(light_GO_ra3_1$qval)
l10_2_5c = mean(light_GO_ra3_2$qval)
r11_1_5c = mean(redox_GO_ra3_1$qval)
r11_2_5c = mean(redox_GO_ra3_2$qval)
m12_1_5c = mean(monosaccharide_GO_ra3_1$qval)
m12_2_5c = mean(monosaccharide_GO_ra3_2$qval)
c13_1_5c = mean(cell_wall_GO_ra3_1$qval)
c13_2_5c = mean(cell_wall_GO_ra3_2$qval)
s14_1_5c = mean(sucrose_GO_ra3_1$qval)
s14_2_5c = mean(sucrose_GO_ra3_2$qval)
t15_1_5c = mean(trehalose_GO_ra3_1$qval)
t15_2_5c = mean(trehalose_GO_ra3_2$qval)
i16_1_5c = mean(intracellular_protein_GO_ra3_1$qval)
i16_2_5c = mean(intracellular_protein_GO_ra3_2$qval)
s17_1_5c = mean(small_GTPase_GO_ra3_1$qval)
s17_2_5c = mean(small_GTPase_GO_ra3_2$qval)
i18_1_5c = mean(intracellular_signaling_GO_ra3_1$qval)
i18_2_5c = mean(intracellular_signaling_GO_ra3_2$qval)
p19_1_5c = mean(proteasome_GO_ra3_1$qval)
p19_2_5c = mean(proteasome_GO_ra3_2$qval)
n20_1_5c = mean(nucleocytoplasmic_GO_ra3_1$qval)
n20_2_5c = mean(nucleocytoplasmic_GO_ra3_2$qval)
m21_1_5c = mean(mitochondrial_transport_GO_ra3_1$qval)
m21_2_5c = mean(mitochondrial_transport_GO_ra3_2$qval)
i22_1_5c = mean(ion_trans_GO_ra3_1$qval)
i22_2_5c = mean(ion_trans_GO_ra3_2$qval)
A23_1_5c = mean(ATP_synthesis_coupled_GO_ra3_1$qval)
A23_2_5c = mean(ATP_synthesis_coupled_GO_ra3_2$qval)

# data.table for p value of each GO category
C1_1mm = c(n1_1_1c, c2_1_1c, v3_1_1c, t4_1_1c, s5_1_1c, t6_1_1c, n7_1_1c, R8_1_1c, G9_1_1c, l10_1_1c, r11_1_1c, m12_1_1c, c13_1_1c, s14_1_1c, t15_1_1c, i16_1_1c, s17_1_1c, i18_1_1c, p19_1_1c, n20_1_1c, m21_1_1c, i22_1_1c, A23_1_1c)
C1_2mm = c(n1_2_1c, c2_2_1c, v3_2_1c, t4_2_1c, s5_2_1c, t6_2_1c, n7_2_1c, R8_2_1c, G9_2_1c, l10_2_1c, r11_2_1c, m12_2_1c, c13_2_1c, s14_2_1c, t15_2_1c, i16_2_1c, s17_2_1c, i18_2_1c, p19_2_1c, n20_2_1c, m21_2_1c, i22_2_1c, A23_2_1c)
C2_1mm = c(n1_1_2c, c2_1_2c, v3_1_2c, t4_1_2c, s5_1_2c, t6_1_2c, n7_1_2c, R8_1_2c, G9_1_2c, l10_1_2c, r11_1_2c, m12_1_2c, c13_1_2c, s14_1_2c, t15_1_2c, i16_1_2c, s17_1_2c, i18_1_2c, p19_1_2c, n20_1_2c, m21_1_2c, i22_1_2c, A23_1_2c)
C2_2mm = c(n1_2_2c, c2_2_2c, v3_2_2c, t4_2_2c, s5_2_2c, t6_2_2c, n7_2_2c, R8_2_2c, G9_2_2c, l10_2_2c, r11_2_2c, m12_2_2c, c13_2_2c, s14_2_2c, t15_2_2c, i16_2_2c, s17_2_2c, i18_2_2c, p19_2_2c, n20_2_2c, m21_2_2c, i22_2_2c, A23_2_2c)
C3_1mm = c(n1_1_3c, c2_1_3c, v3_1_3c, t4_1_3c, s5_1_3c, t6_1_3c, n7_1_3c, R8_1_3c, G9_1_3c, l10_1_3c, r11_1_3c, m12_1_3c, c13_1_3c, s14_1_3c, t15_1_3c, i16_1_3c, s17_1_3c, i18_1_3c, p19_1_3c, n20_1_3c, m21_1_3c, i22_1_3c, A23_1_3c)
C3_2mm = c(n1_2_3c, c2_2_3c, v3_2_3c, t4_2_3c, s5_2_3c, t6_2_3c, n7_2_3c, R8_2_3c, G9_2_3c, l10_2_3c, r11_2_3c, m12_2_3c, c13_2_3c, s14_2_3c, t15_2_3c, i16_2_3c, s17_2_3c, i18_2_3c, p19_2_3c, n20_2_3c, m21_2_3c, i22_2_3c, A23_2_3c)
C4_1mm = c(n1_1_4c, c2_1_4c, v3_1_4c, t4_1_4c, s5_1_4c, t6_1_4c, n7_1_4c, R8_1_4c, G9_1_4c, l10_1_4c, r11_1_4c, m12_1_4c, c13_1_4c, s14_1_4c, t15_1_4c, i16_1_4c, s17_1_4c, i18_1_4c, p19_1_4c, n20_1_4c, m21_1_4c, i22_1_4c, A23_1_4c)
C4_2mm = c(n1_2_4c, c2_2_4c, v3_2_4c, t4_2_4c, s5_2_4c, t6_2_4c, n7_2_4c, R8_2_4c, G9_2_4c, l10_2_4c, r11_2_4c, m12_2_4c, c13_2_4c, s14_2_4c, t15_2_4c, i16_2_4c, s17_2_4c, i18_2_4c, p19_2_4c, n20_2_4c, m21_2_4c, i22_2_4c, A23_2_4c)
C5_1mm = c(n1_1_5c, c2_1_5c, v3_1_5c, t4_1_5c, s5_1_5c, t6_1_5c, n7_1_5c, R8_1_5c, G9_1_5c, l10_1_5c, r11_1_5c, m12_1_5c, c13_1_5c, s14_1_5c, t15_1_5c, i16_1_5c, s17_1_5c, i18_1_5c, p19_1_5c, n20_1_5c, m21_1_5c, i22_1_5c, A23_1_5c)
C5_2mm = c(n1_2_5c, c2_2_5c, v3_2_5c, t4_2_5c, s5_2_5c, t6_2_5c, n7_2_5c, R8_2_5c, G9_2_5c, l10_2_5c, r11_2_5c, m12_2_5c, c13_2_5c, s14_2_5c, t15_2_5c, i16_2_5c, s17_2_5c, i18_2_5c, p19_2_5c, n20_2_5c, m21_2_5c, i22_2_5c, A23_2_5c)
GO_enrich = data.frame(C1_1mm, C1_2mm, C2_1mm, C2_2mm, C3_1mm, C3_2mm, C4_1mm, C4_2mm, C5_1mm, C5_2mm)

#Naming the raws
GO_enrichment = c("nucleosome", "chromatin", "vesicle", "temperature", "stress", "transcription", "nitrogen", "RNA_synthetase", "G_protein", "light", "redox", "monosaccharide", "cell_wall", "sucrose", "trehalose", "intracellular_protein", "small_GTPase", "intracellular_signaling", "proteasome", "nucleocytoplasmic", "mitochondrial_transport", "ion_trans", "ATP_synthesis_coupled")
GO_enrich = cbind(GO_enrichment, GO_enrich)
write.csv(GO_enrich, file = "GO_enrich.csv") #leave a file just in case

# Necessary preparation for Heatmap
install.packages("reshape2")
install.packages("gplots")
install.packages("RColorBrewer")
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(gplots)
row.names(GO_enrich) = GO_enrich$GO_enrichment
GO_enrich = GO_enrich[,2:11]
GO_enr_mat = data.matrix(GO_enrich)

# Plotting the Heatmap
breaks = seq(1e-192,1e-2,length.out=1000)
gradient1 = colorpanel( sum( breaks[-1]<=0.001 ), "red", "orange" )
gradient2 = colorpanel( sum( breaks[-1]>0.001 ), "orange", "grey" )
hm.colors = c(gradient1,gradient2)
GO_enr_heat = heatmap.2(GO_enr_mat, Rowv=NA, Colv=NA, dendrogram = "none", breaks = breaks, scale="none", density.info="none", col=hm.colors, trace = "none", na.color="white", margins=c(5,10))
