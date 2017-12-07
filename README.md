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

##########################################
# for ra1@2mm 
f2_ra1_2 = filter(f2_df, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
ra1_2 = as.character(f2_ra1_2$gene_id)

# for ra2@2mm 
f2_ra2_2 = filter(f2_df, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
ra2_2 = as.character(f2_ra2_2$gene_id)

# for ra3@2mm 
f2_ra3_2 = filter(f2_df, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")
ra3_2 = as.character(f2_ra3_2$gene_id)

# Venn Diagram for 2mm
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
## Fig 2C
# ra1
grep "EntrezGene" S_Table6.csv > nucleosome.csv
grep "nucleosome" S_Table6.csv >> nucleosome.csv

grep "EntrezGene" S_Table6.csv > chromatin.csv
grep "chromatin" S_Table6.csv >> chromatin.csv

grep "EntrezGene" S_Table6.csv > vesicle.csv
grep "vesicle" S_Table6.csv >> vesicle.csv

grep "EntrezGene" S_Table6.csv > temperature.csv
grep "temperature" S_Table6.csv >> temperature.csv

grep "EntrezGene" S_Table6.csv > stress.csv
grep "stress" S_Table6.csv >> stress.csv

grep "EntrezGene" S_Table6.csv > transcription.csv
grep "transcription" S_Table6.csv >> transcription.csv

grep "EntrezGene" S_Table6.csv > nitrogen.csv
grep "nitrogen" S_Table6.csv >> nitrogen.csv

grep "EntrezGene" S_Table6.csv > RNA_synthetase.csv
grep "RNA synthetase" S_Table6.csv >> RNA_synthetase.csv

grep "EntrezGene" S_Table6.csv > G_protein.csv
grep "G-protein" S_Table6.csv >> G_protein.csv

grep "EntrezGene" S_Table6.csv > light.csv
grep "light" S_Table6.csv >> light.csv

grep "EntrezGene" S_Table6.csv > redox.csv
grep "redox" S_Table6.csv >> redox.csv

grep "EntrezGene" S_Table6.csv > monosaccharide.csv
grep "monosaccharide" S_Table6.csv >> monosaccharide.csv

grep "EntrezGene" S_Table6.csv > cell_wall.csv
grep "cell wall" S_Table6.csv >> cell_wall.csv

grep "EntrezGene" S_Table6.csv > sucrose.csv
grep "sucrose" S_Table6.csv >> sucrose.csv

grep "EntrezGene" S_Table6.csv > trehalose.csv
grep "trehalose" S_Table6.csv >> trehalose.csv

grep "EntrezGene" S_Table6.csv > intracellular_protein.csv
grep "intracellular protein" S_Table6.csv >> intracellular_protein.csv

grep "EntrezGene" S_Table6.csv > small_GTPase.csv
grep "small GTPase" S_Table6.csv >> small_GTPase.csv

grep "EntrezGene" S_Table6.csv > intracellular_signaling.csv
grep "intracellular signaling" S_Table6.csv >> intracellular_signaling.csv

grep "EntrezGene" S_Table6.csv > proteasome.csv
grep "proteasome" S_Table6.csv >> proteasome.csv

grep "EntrezGene" S_Table6.csv > nucleocytoplasmic.csv
grep "nucleocytoplasmic" S_Table6.csv >> nucleocytoplasmic.csv

grep "EntrezGene" S_Table6.csv > mitochondrial_transport.csv
grep "mitochondrial transport" S_Table6.csv >> mitochondrial_transport.csv

grep "EntrezGene" S_Table6.csv > ion_trans.csv
grep "ion trans" S_Table6.csv >> ion_trans.csv

grep "EntrezGene" S_Table6.csv > ATP_synthesis_coupled.csv
grep "ATP synthesis coupled" S_Table6.csv >> ATP_synthesis_coupled.csv
