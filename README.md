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

# Back to R
library(dplyr)
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
# ra1+ra2+ra3@1mm
nucleosome_ra1_ra2_ra3_1mm = join_all (list (nucleosome_GO_ra1_1, nucleosome_GO_ra2_1, nucleosome_GO_ra3_1), by = "gene_id")
# ra1+ra2@1mm
nucleosome_ra1_ra2_1mm = merge (nucleosome_GO_ra1_1, nucleosome_GO_ra2_1, by = "gene_id")
# ra2+ra3@1mm
nucleosome_ra2_ra3_1mm = merge (nucleosome_GO_ra2_1, nucleosome_GO_ra3_1, by = "gene_id")
# ra1+ra2+ra3@2mm
nucleosome_ra1_ra2_ra3_2mm = join_all (list (nucleosome_GO_ra1_2, nucleosome_GO_ra2_2, nucleosome_GO_ra3_1), by = "gene_id")
# ra1+ra2@2mm
nucleosome_ra1_ra2_2mm = merge (nucleosome_GO_ra1_2, nucleosome_GO_ra2_2, by = "gene_id")
# ra2+ra3@2mm
nucleosome_ra2_ra3_2mm = merge (nucleosome_GO_ra2_2, nucleosome_GO_ra3_2, by = "gene_id")
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
# ra1+ra2+ra3@1mm
chromatin_ra1_ra2_ra3_1mm = join_all (list (chromatin_GO_ra1_1, chromatin_GO_ra2_1, chromatin_GO_ra3_1), by = "gene_id")
# ra1+ra2@1mm
chromatin_ra1_ra2_1mm = merge (chromatin_GO_ra1_1, chromatin_GO_ra2_1, by = "gene_id")
# ra2+ra3@1mm
chromatin_ra2_ra3_1mm = merge (chromatin_GO_ra2_1, chromatin_GO_ra3_1, by = "gene_id")
# ra1+ra2+ra3@2mm
chromatin_ra1_ra2_ra3_2mm = join_all (list (chromatin_GO_ra1_2, chromatin_GO_ra2_2, chromatin_GO_ra3_1), by = "gene_id")
# ra1+ra2@2mm
chromatin_ra1_ra2_2mm = merge (chromatin_GO_ra1_2, chromatin_GO_ra2_2, by = "gene_id")
# ra2+ra3@2mm
chromatin_ra2_ra3_2mm = merge (chromatin_GO_ra2_2, chromatin_GO_ra3_2, by = "gene_id")
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
# ra1+ra2+ra3@1mm
vesicle_ra1_ra2_ra3_1mm = join_all (list (vesicle_GO_ra1_1, vesicle_GO_ra2_1, vesicle_GO_ra3_1), by = "gene_id")
# ra1+ra2@1mm
vesicle_ra1_ra2_1mm = merge (vesicle_GO_ra1_1, vesicle_GO_ra2_1, by = "gene_id")
# ra2+ra3@1mm
vesicle_ra2_ra3_1mm = merge (vesicle_GO_ra2_1, vesicle_GO_ra3_1, by = "gene_id")
# ra1+ra2+ra3@2mm
vesicle_ra1_ra2_ra3_2mm = join_all (list (vesicle_GO_ra1_2, vesicle_GO_ra2_2, vesicle_GO_ra3_1), by = "gene_id")
# ra1+ra2@2mm
vesicle_ra1_ra2_2mm = merge (vesicle_GO_ra1_2, vesicle_GO_ra2_2, by = "gene_id")
# ra2+ra3@2mm
vesicle_ra2_ra3_2mm = merge (vesicle_GO_ra2_2, vesicle_GO_ra3_2, by = "gene_id")
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
# ra1+ra2+ra3@1mm
temperature_ra1_ra2_ra3_1mm = join_all (list (temperature_GO_ra1_1, temperature_GO_ra2_1, temperature_GO_ra3_1), by = "gene_id")
# ra1+ra2@1mm
temperature_ra1_ra2_1mm = merge (temperature_GO_ra1_1, temperature_GO_ra2_1, by = "gene_id")
# ra2+ra3@1mm
temperature_ra2_ra3_1mm = merge (temperature_GO_ra2_1, temperature_GO_ra3_1, by = "gene_id")
# ra1+ra2+ra3@2mm
temperature_ra1_ra2_ra3_2mm = join_all (list (temperature_GO_ra1_2, temperature_GO_ra2_2, temperature_GO_ra3_1), by = "gene_id")
# ra1+ra2@2mm
temperature_ra1_ra2_2mm = merge (temperature_GO_ra1_2, temperature_GO_ra2_2, by = "gene_id")
# ra2+ra3@2mm
temperature_ra2_ra3_2mm = merge (temperature_GO_ra2_2, temperature_GO_ra3_2, by = "gene_id")
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
# ra1+ra2+ra3@1mm
stress_ra1_ra2_ra3_1mm = join_all (list (stress_GO_ra1_1, stress_GO_ra2_1, stress_GO_ra3_1), by = "gene_id")
# ra1+ra2@1mm
stress_ra1_ra2_1mm = merge (stress_GO_ra1_1, stress_GO_ra2_1, by = "gene_id")
# ra2+ra3@1mm
stress_ra2_ra3_1mm = merge (stress_GO_ra2_1, stress_GO_ra3_1, by = "gene_id")
# ra1+ra2+ra3@2mm
stress_ra1_ra2_ra3_2mm = join_all (list (stress_GO_ra1_2, stress_GO_ra2_2, stress_GO_ra3_1), by = "gene_id")
# ra1+ra2@2mm
stress_ra1_ra2_2mm = merge (stress_GO_ra1_2, stress_GO_ra2_2, by = "gene_id")
# ra2+ra3@2mm
stress_ra2_ra3_2mm = merge (stress_GO_ra2_2, stress_GO_ra3_2, by = "gene_id")
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
# ra1+ra2+ra3@1mm
transcription_ra1_ra2_ra3_1mm = join_all (list (transcription_GO_ra1_1, transcription_GO_ra2_1, transcription_GO_ra3_1), by = "gene_id")
# ra1+ra2@1mm
transcription_ra1_ra2_1mm = merge (transcription_GO_ra1_1, transcription_GO_ra2_1, by = "gene_id")
# ra2+ra3@1mm
transcription_ra2_ra3_1mm = merge (transcription_GO_ra2_1, transcription_GO_ra3_1, by = "gene_id")
# ra1+ra2+ra3@2mm
transcription_ra1_ra2_ra3_2mm = join_all (list (transcription_GO_ra1_2, transcription_GO_ra2_2, transcription_GO_ra3_1), by = "gene_id")
# ra1+ra2@2mm
transcription_ra1_ra2_2mm = merge (transcription_GO_ra1_2, transcription_GO_ra2_2, by = "gene_id")
# ra2+ra3@2mm
transcription_ra2_ra3_2mm = merge (transcription_GO_ra2_2, transcription_GO_ra3_2, by = "gene_id")
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
# ra1+ra2+ra3@1mm
nitrogen_ra1_ra2_ra3_1mm = join_all (list (nitrogen_GO_ra1_1, nitrogen_GO_ra2_1, nitrogen_GO_ra3_1), by = "gene_id")
# ra1+ra2@1mm
nitrogen_ra1_ra2_1mm = merge (nitrogen_GO_ra1_1, nitrogen_GO_ra2_1, by = "gene_id")
# ra2+ra3@1mm
nitrogen_ra2_ra3_1mm = merge (nitrogen_GO_ra2_1, nitrogen_GO_ra3_1, by = "gene_id")
# ra1+ra2+ra3@2mm
nitrogen_ra1_ra2_ra3_2mm = join_all (list (nitrogen_GO_ra1_2, nitrogen_GO_ra2_2, nitrogen_GO_ra3_1), by = "gene_id")
# ra1+ra2@2mm
nitrogen_ra1_ra2_2mm = merge (nitrogen_GO_ra1_2, nitrogen_GO_ra2_2, by = "gene_id")
# ra2+ra3@2mm
nitrogen_ra2_ra3_2mm = merge (nitrogen_GO_ra2_2, nitrogen_GO_ra3_2, by = "gene_id")
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
# ra1+ra2+ra3@1mm
RNA_synthetase_ra1_ra2_ra3_1mm = join_all (list (RNA_synthetase_GO_ra1_1, RNA_synthetase_GO_ra2_1, RNA_synthetase_GO_ra3_1), by = "gene_id")
# ra1+ra2@1mm
RNA_synthetase_ra1_ra2_1mm = merge (RNA_synthetase_GO_ra1_1, RNA_synthetase_GO_ra2_1, by = "gene_id")
# ra2+ra3@1mm
RNA_synthetase_ra2_ra3_1mm = merge (RNA_synthetase_GO_ra2_1, RNA_synthetase_GO_ra3_1, by = "gene_id")
# ra1+ra2+ra3@2mm
RNA_synthetase_ra1_ra2_ra3_2mm = join_all (list (RNA_synthetase_GO_ra1_2, RNA_synthetase_GO_ra2_2, RNA_synthetase_GO_ra3_1), by = "gene_id")
# ra1+ra2@2mm
RNA_synthetase_ra1_ra2_2mm = merge (RNA_synthetase_GO_ra1_2, RNA_synthetase_GO_ra2_2, by = "gene_id")
# ra2+ra3@2mm
RNA_synthetase_ra2_ra3_2mm = merge (RNA_synthetase_GO_ra2_2, RNA_synthetase_GO_ra3_2, by = "gene_id")
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
# ra1+ra2+ra3@1mm
G_protein_ra1_ra2_ra3_1mm = join_all (list (G_protein_GO_ra1_1, G_protein_GO_ra2_1, G_protein_GO_ra3_1), by = "gene_id")
# ra1+ra2@1mm
G_protein_ra1_ra2_1mm = merge (G_protein_GO_ra1_1, G_protein_GO_ra2_1, by = "gene_id")
# ra2+ra3@1mm
G_protein_ra2_ra3_1mm = merge (G_protein_GO_ra2_1, G_protein_GO_ra3_1, by = "gene_id")
# ra1+ra2+ra3@2mm
G_protein_ra1_ra2_ra3_2mm = join_all (list (G_protein_GO_ra1_2, G_protein_GO_ra2_2, G_protein_GO_ra3_1), by = "gene_id")
# ra1+ra2@2mm
G_protein_ra1_ra2_2mm = merge (G_protein_GO_ra1_2, G_protein_GO_ra2_2, by = "gene_id")
# ra2+ra3@2mm
G_protein_ra2_ra3_2mm = merge (G_protein_GO_ra2_2, G_protein_GO_ra3_2, by = "gene_id")
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
# ra1+ra2+ra3@1mm
light_ra1_ra2_ra3_1mm = join_all (list (light_GO_ra1_1, light_GO_ra2_1, light_GO_ra3_1), by = "gene_id")
# ra1+ra2@1mm
light_ra1_ra2_1mm = merge (light_GO_ra1_1, light_GO_ra2_1, by = "gene_id")
# ra2+ra3@1mm
light_ra2_ra3_1mm = merge (light_GO_ra2_1, light_GO_ra3_1, by = "gene_id")
# ra1+ra2+ra3@2mm
light_ra1_ra2_ra3_2mm = join_all (list (light_GO_ra1_2, light_GO_ra2_2, light_GO_ra3_1), by = "gene_id")
# ra1+ra2@2mm
light_ra1_ra2_2mm = merge (light_GO_ra1_2, light_GO_ra2_2, by = "gene_id")
# ra2+ra3@2mm
light_ra2_ra3_2mm = merge (light_GO_ra2_2, light_GO_ra3_2, by = "gene_id")
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
# ra1+ra2+ra3@1mm
redox_ra1_ra2_ra3_1mm = join_all (list (redox_GO_ra1_1, redox_GO_ra2_1, redox_GO_ra3_1), by = "gene_id")
# ra1+ra2@1mm
redox_ra1_ra2_1mm = merge (redox_GO_ra1_1, redox_GO_ra2_1, by = "gene_id")
# ra2+ra3@1mm
redox_ra2_ra3_1mm = merge (redox_GO_ra2_1, redox_GO_ra3_1, by = "gene_id")
# ra1+ra2+ra3@2mm
redox_ra1_ra2_ra3_2mm = join_all (list (redox_GO_ra1_2, redox_GO_ra2_2, redox_GO_ra3_1), by = "gene_id")
# ra1+ra2@2mm
redox_ra1_ra2_2mm = merge (redox_GO_ra1_2, redox_GO_ra2_2, by = "gene_id")
# ra2+ra3@2mm
redox_ra2_ra3_2mm = merge (redox_GO_ra2_2, redox_GO_ra3_2, by = "gene_id")
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
# ra1+ra2+ra3@1mm
monosaccharide_ra1_ra2_ra3_1mm = join_all (list (monosaccharide_GO_ra1_1, monosaccharide_GO_ra2_1, monosaccharide_GO_ra3_1), by = "gene_id")
# ra1+ra2@1mm
monosaccharide_ra1_ra2_1mm = merge (monosaccharide_GO_ra1_1, monosaccharide_GO_ra2_1, by = "gene_id")
# ra2+ra3@1mm
monosaccharide_ra2_ra3_1mm = merge (monosaccharide_GO_ra2_1, monosaccharide_GO_ra3_1, by = "gene_id")
# ra1+ra2+ra3@2mm
monosaccharide_ra1_ra2_ra3_2mm = join_all (list (monosaccharide_GO_ra1_2, monosaccharide_GO_ra2_2, monosaccharide_GO_ra3_1), by = "gene_id")
# ra1+ra2@2mm
monosaccharide_ra1_ra2_2mm = merge (monosaccharide_GO_ra1_2, monosaccharide_GO_ra2_2, by = "gene_id")
# ra2+ra3@2mm
monosaccharide_ra2_ra3_2mm = merge (monosaccharide_GO_ra2_2, monosaccharide_GO_ra3_2, by = "gene_id")
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
# ra1+ra2+ra3@1mm
cell_wall_ra1_ra2_ra3_1mm = join_all (list (cell_wall_GO_ra1_1, cell_wall_GO_ra2_1, cell_wall_GO_ra3_1), by = "gene_id")
# ra1+ra2@1mm
cell_wall_ra1_ra2_1mm = merge (cell_wall_GO_ra1_1, cell_wall_GO_ra2_1, by = "gene_id")
# ra2+ra3@1mm
cell_wall_ra2_ra3_1mm = merge (cell_wall_GO_ra2_1, cell_wall_GO_ra3_1, by = "gene_id")
# ra1+ra2+ra3@2mm
cell_wall_ra1_ra2_ra3_2mm = join_all (list (cell_wall_GO_ra1_2, cell_wall_GO_ra2_2, cell_wall_GO_ra3_1), by = "gene_id")
# ra1+ra2@2mm
cell_wall_ra1_ra2_2mm = merge (cell_wall_GO_ra1_2, cell_wall_GO_ra2_2, by = "gene_id")
# ra2+ra3@2mm
cell_wall_ra2_ra3_2mm = merge (cell_wall_GO_ra2_2, cell_wall_GO_ra3_2, by = "gene_id")
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
# ra1+ra2+ra3@1mm
sucrose_ra1_ra2_ra3_1mm = join_all (list (sucrose_GO_ra1_1, sucrose_GO_ra2_1, sucrose_GO_ra3_1), by = "gene_id")
# ra1+ra2@1mm
sucrose_ra1_ra2_1mm = merge (sucrose_GO_ra1_1, sucrose_GO_ra2_1, by = "gene_id")
# ra2+ra3@1mm
sucrose_ra2_ra3_1mm = merge (sucrose_GO_ra2_1, sucrose_GO_ra3_1, by = "gene_id")
# ra1+ra2+ra3@2mm
sucrose_ra1_ra2_ra3_2mm = join_all (list (sucrose_GO_ra1_2, sucrose_GO_ra2_2, sucrose_GO_ra3_1), by = "gene_id")
# ra1+ra2@2mm
sucrose_ra1_ra2_2mm = merge (sucrose_GO_ra1_2, sucrose_GO_ra2_2, by = "gene_id")
# ra2+ra3@2mm
sucrose_ra2_ra3_2mm = merge (sucrose_GO_ra2_2, sucrose_GO_ra3_2, by = "gene_id")
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
# ra1+ra2+ra3@1mm
trehalose_ra1_ra2_ra3_1mm = join_all (list (trehalose_GO_ra1_1, trehalose_GO_ra2_1, trehalose_GO_ra3_1), by = "gene_id")
# ra1+ra2@1mm
trehalose_ra1_ra2_1mm = merge (trehalose_GO_ra1_1, trehalose_GO_ra2_1, by = "gene_id")
# ra2+ra3@1mm
trehalose_ra2_ra3_1mm = merge (trehalose_GO_ra2_1, trehalose_GO_ra3_1, by = "gene_id")
# ra1+ra2+ra3@2mm
trehalose_ra1_ra2_ra3_2mm = join_all (list (trehalose_GO_ra1_2, trehalose_GO_ra2_2, trehalose_GO_ra3_1), by = "gene_id")
# ra1+ra2@2mm
trehalose_ra1_ra2_2mm = merge (trehalose_GO_ra1_2, trehalose_GO_ra2_2, by = "gene_id")
# ra2+ra3@2mm
trehalose_ra2_ra3_2mm = merge (trehalose_GO_ra2_2, trehalose_GO_ra3_2, by = "gene_id")
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
# ra1+ra2+ra3@1mm
intracellular_protein_ra1_ra2_ra3_1mm = join_all (list (intracellular_protein_GO_ra1_1, intracellular_protein_GO_ra2_1, intracellular_protein_GO_ra3_1), by = "gene_id")
# ra1+ra2@1mm
intracellular_protein_ra1_ra2_1mm = merge (intracellular_protein_GO_ra1_1, intracellular_protein_GO_ra2_1, by = "gene_id")
# ra2+ra3@1mm
intracellular_protein_ra2_ra3_1mm = merge (intracellular_protein_GO_ra2_1, intracellular_protein_GO_ra3_1, by = "gene_id")
# ra1+ra2+ra3@2mm
intracellular_protein_ra1_ra2_ra3_2mm = join_all (list (intracellular_protein_GO_ra1_2, intracellular_protein_GO_ra2_2, intracellular_protein_GO_ra3_1), by = "gene_id")
# ra1+ra2@2mm
intracellular_protein_ra1_ra2_2mm = merge (intracellular_protein_GO_ra1_2, intracellular_protein_GO_ra2_2, by = "gene_id")
# ra2+ra3@2mm
intracellular_protein_ra2_ra3_2mm = merge (intracellular_protein_GO_ra2_2, intracellular_protein_GO_ra3_2, by = "gene_id")
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
# ra1+ra2+ra3@1mm
small_GTPase_ra1_ra2_ra3_1mm = join_all (list (small_GTPase_GO_ra1_1, small_GTPase_GO_ra2_1, small_GTPase_GO_ra3_1), by = "gene_id")
# ra1+ra2@1mm
small_GTPase_ra1_ra2_1mm = merge (small_GTPase_GO_ra1_1, small_GTPase_GO_ra2_1, by = "gene_id")
# ra2+ra3@1mm
small_GTPase_ra2_ra3_1mm = merge (small_GTPase_GO_ra2_1, small_GTPase_GO_ra3_1, by = "gene_id")
# ra1+ra2+ra3@2mm
small_GTPase_ra1_ra2_ra3_2mm = join_all (list (small_GTPase_GO_ra1_2, small_GTPase_GO_ra2_2, small_GTPase_GO_ra3_1), by = "gene_id")
# ra1+ra2@2mm
small_GTPase_ra1_ra2_2mm = merge (small_GTPase_GO_ra1_2, small_GTPase_GO_ra2_2, by = "gene_id")
# ra2+ra3@2mm
small_GTPase_ra2_ra3_2mm = merge (small_GTPase_GO_ra2_2, small_GTPase_GO_ra3_2, by = "gene_id")
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
# ra1+ra2+ra3@1mm
intracellular_signaling_ra1_ra2_ra3_1mm = join_all (list (intracellular_signaling_GO_ra1_1, intracellular_signaling_GO_ra2_1, intracellular_signaling_GO_ra3_1), by = "gene_id")
# ra1+ra2@1mm
intracellular_signaling_ra1_ra2_1mm = merge (intracellular_signaling_GO_ra1_1, intracellular_signaling_GO_ra2_1, by = "gene_id")
# ra2+ra3@1mm
intracellular_signaling_ra2_ra3_1mm = merge (intracellular_signaling_GO_ra2_1, intracellular_signaling_GO_ra3_1, by = "gene_id")
# ra1+ra2+ra3@2mm
intracellular_signaling_ra1_ra2_ra3_2mm = join_all (list (intracellular_signaling_GO_ra1_2, intracellular_signaling_GO_ra2_2, intracellular_signaling_GO_ra3_1), by = "gene_id")
# ra1+ra2@2mm
intracellular_signaling_ra1_ra2_2mm = merge (intracellular_signaling_GO_ra1_2, intracellular_signaling_GO_ra2_2, by = "gene_id")
# ra2+ra3@2mm
intracellular_signaling_ra2_ra3_2mm = merge (intracellular_signaling_GO_ra2_2, intracellular_signaling_GO_ra3_2, by = "gene_id")
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
# ra1+ra2+ra3@1mm
proteasome_ra1_ra2_ra3_1mm = join_all (list (proteasome_GO_ra1_1, proteasome_GO_ra2_1, proteasome_GO_ra3_1), by = "gene_id")
# ra1+ra2@1mm
proteasome_ra1_ra2_1mm = merge (proteasome_GO_ra1_1, proteasome_GO_ra2_1, by = "gene_id")
# ra2+ra3@1mm
proteasome_ra2_ra3_1mm = merge (proteasome_GO_ra2_1, proteasome_GO_ra3_1, by = "gene_id")
# ra1+ra2+ra3@2mm
proteasome_ra1_ra2_ra3_2mm = join_all (list (proteasome_GO_ra1_2, proteasome_GO_ra2_2, proteasome_GO_ra3_1), by = "gene_id")
# ra1+ra2@2mm
proteasome_ra1_ra2_2mm = merge (proteasome_GO_ra1_2, proteasome_GO_ra2_2, by = "gene_id")
# ra2+ra3@2mm
proteasome_ra2_ra3_2mm = merge (proteasome_GO_ra2_2, proteasome_GO_ra3_2, by = "gene_id")
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
# ra1+ra2+ra3@1mm
nucleocytoplasmic_ra1_ra2_ra3_1mm = join_all (list (nucleocytoplasmic_GO_ra1_1, nucleocytoplasmic_GO_ra2_1, nucleocytoplasmic_GO_ra3_1), by = "gene_id")
# ra1+ra2@1mm
nucleocytoplasmic_ra1_ra2_1mm = merge (nucleocytoplasmic_GO_ra1_1, nucleocytoplasmic_GO_ra2_1, by = "gene_id")
# ra2+ra3@1mm
nucleocytoplasmic_ra2_ra3_1mm = merge (nucleocytoplasmic_GO_ra2_1, nucleocytoplasmic_GO_ra3_1, by = "gene_id")
# ra1+ra2+ra3@2mm
nucleocytoplasmic_ra1_ra2_ra3_2mm = join_all (list (nucleocytoplasmic_GO_ra1_2, nucleocytoplasmic_GO_ra2_2, nucleocytoplasmic_GO_ra3_1), by = "gene_id")
# ra1+ra2@2mm
nucleocytoplasmic_ra1_ra2_2mm = merge (nucleocytoplasmic_GO_ra1_2, nucleocytoplasmic_GO_ra2_2, by = "gene_id")
# ra2+ra3@2mm
nucleocytoplasmic_ra2_ra3_2mm = merge (nucleocytoplasmic_GO_ra2_2, nucleocytoplasmic_GO_ra3_2, by = "gene_id")
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
# ra1+ra2+ra3@1mm
mitochondrial_transport_ra1_ra2_ra3_1mm = join_all (list (mitochondrial_transport_GO_ra1_1, mitochondrial_transport_GO_ra2_1, mitochondrial_transport_GO_ra3_1), by = "gene_id")
# ra1+ra2@1mm
mitochondrial_transport_ra1_ra2_1mm = merge (mitochondrial_transport_GO_ra1_1, mitochondrial_transport_GO_ra2_1, by = "gene_id")
# ra2+ra3@1mm
mitochondrial_transport_ra2_ra3_1mm = merge (mitochondrial_transport_GO_ra2_1, mitochondrial_transport_GO_ra3_1, by = "gene_id")
# ra1+ra2+ra3@2mm
mitochondrial_transport_ra1_ra2_ra3_2mm = join_all (list (mitochondrial_transport_GO_ra1_2, mitochondrial_transport_GO_ra2_2, mitochondrial_transport_GO_ra3_1), by = "gene_id")
# ra1+ra2@2mm
mitochondrial_transport_ra1_ra2_2mm = merge (mitochondrial_transport_GO_ra1_2, mitochondrial_transport_GO_ra2_2, by = "gene_id")
# ra2+ra3@2mm
mitochondrial_transport_ra2_ra3_2mm = merge (mitochondrial_transport_GO_ra2_2, mitochondrial_transport_GO_ra3_2, by = "gene_id")
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
# ra1+ra2+ra3@1mm
ion_trans_ra1_ra2_ra3_1mm = join_all (list (ion_trans_GO_ra1_1, ion_trans_GO_ra2_1, ion_trans_GO_ra3_1), by = "gene_id")
# ra1+ra2@1mm
ion_trans_ra1_ra2_1mm = merge (ion_trans_GO_ra1_1, ion_trans_GO_ra2_1, by = "gene_id")
# ra2+ra3@1mm
ion_trans_ra2_ra3_1mm = merge (ion_trans_GO_ra2_1, ion_trans_GO_ra3_1, by = "gene_id")
# ra1+ra2+ra3@2mm
ion_trans_ra1_ra2_ra3_2mm = join_all (list (ion_trans_GO_ra1_2, ion_trans_GO_ra2_2, ion_trans_GO_ra3_1), by = "gene_id")
# ra1+ra2@2mm
ion_trans_ra1_ra2_2mm = merge (ion_trans_GO_ra1_2, ion_trans_GO_ra2_2, by = "gene_id")
# ra2+ra3@2mm
ion_trans_ra2_ra3_2mm = merge (ion_trans_GO_ra2_2, ion_trans_GO_ra3_2, by = "gene_id")
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
# ra1+ra2+ra3@1mm
ATP_synthesis_coupled_ra1_ra2_ra3_1mm = join_all (list (ATP_synthesis_coupled_GO_ra1_1, ATP_synthesis_coupled_GO_ra2_1, ATP_synthesis_coupled_GO_ra3_1), by = "gene_id")
# ra1+ra2@1mm
ATP_synthesis_coupled_ra1_ra2_1mm = merge (ATP_synthesis_coupled_GO_ra1_1, ATP_synthesis_coupled_GO_ra2_1, by = "gene_id")
# ra2+ra3@1mm
ATP_synthesis_coupled_ra2_ra3_1mm = merge (ATP_synthesis_coupled_GO_ra2_1, ATP_synthesis_coupled_GO_ra3_1, by = "gene_id")
# ra1+ra2+ra3@2mm
ATP_synthesis_coupled_ra1_ra2_ra3_2mm = join_all (list (ATP_synthesis_coupled_GO_ra1_2, ATP_synthesis_coupled_GO_ra2_2, ATP_synthesis_coupled_GO_ra3_1), by = "gene_id")
# ra1+ra2@2mm
ATP_synthesis_coupled_ra1_ra2_2mm = merge (ATP_synthesis_coupled_GO_ra1_2, ATP_synthesis_coupled_GO_ra2_2, by = "gene_id")
# ra2+ra3@2mm
ATP_synthesis_coupled_ra2_ra3_2mm = merge (ATP_synthesis_coupled_GO_ra2_2, ATP_synthesis_coupled_GO_ra3_2, by = "gene_id")


# Getting average P_value


1_1mm = c(
# column 1
(mean(GGG_ra2_ra3_2mm$qval.x) + mean(GGG_ra2_ra3_2mm$qval.y))/2
1n_1_1c
1n_2_1c
2c_1_1c
2c_2_1c
3v_1_1c
3v_2_1c
4t_1_1c
4t_2_1c
5s_1_1c
5s_2_1c
6t_1_1c
6t_2_1c
7n_1_1c
7n_2_1c
8R_1_1c
8R_2_1c
9G_1_1c
9G_2_1c
10l_1_1c
10l_2_1c
11r_1_1c
11r_2_1c
12m_1_1c
12m_2_1c
13c_1_1c
13c_2_1c
14s_1_1c
14s_2_1c
15t_1_1c
15t_2_1c
16i_1_1c
16i_2_1c
17s_1_1c
17s_2_1c
18i_1_1c
18i_2_1c
19p_1_1c
19p_2_1c
20n_1_1c
20n_2_1c
21m_1_1c
21m_2_1c
22i_1_1c
22i_2_1c
23A_1_1c
23A_2_1c







# for ra1@1mm 
TF_GO_ra1_1 = filter(f2_TF, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
TF_GO_ra1_1 = arrange (TF_GO_ra1_1, gene_id)
# for ra2@1mm 
TF_GO_ra2_1 = filter(f2_TF, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
TF_GO_ra2_1 = arrange (TF_GO_ra2_1, gene_id)
# for ra3@1mm 
TF_GO_ra3_1 = filter(f2_TF, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
TF_GO_ra3_1 = arrange (TF_GO_ra3_1, gene_id)
# for ra1@2mm 
TF_GO_ra1_2 = filter(f2_TF, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
TF_GO_ra1_2 = arrange (TF_GO_ra1_2, gene_id)
# for ra2@2mm 
TF_GO_ra2_2 = filter(f2_TF, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
TF_GO_ra2_2 = arrange (TF_GO_ra2_2, gene_id)
# for ra3@2mm 
TF_GO_ra3_2 = filter(f2_TF, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")
TF_GO_ra3_2 = arrange (TF_GO_ra3_2, gene_id)
# ra1+ra2+ra3@1mm
TF_ra1_ra2_ra3_1mm = join_all (list (TF_GO_ra1_1, TF_GO_ra2_1, TF_GO_ra3_1), by = "gene_id")
# ra1+ra2@1mm
TF_ra1_ra2_1mm = merge (TF_GO_ra1_1, TF_GO_ra2_1, by = "gene_id")
# ra2+ra3@1mm
TF_ra2_ra3_1mm = merge (TF_GO_ra2_1, TF_GO_ra3_1, by = "gene_id")
# ra1+ra2+ra3@2mm
TF_ra1_ra2_ra3_2mm = join_all (list (TF_GO_ra1_2, TF_GO_ra2_2, TF_GO_ra3_1), by = "gene_id")
# ra1+ra2@2mm
TF_ra1_ra2_2mm = merge (TF_GO_ra1_2, TF_GO_ra2_2, by = "gene_id")
# ra2+ra3@2mm
TF_ra2_ra3_2mm = merge (TF_GO_ra2_2, TF_GO_ra3_2, by = "gene_id")
#



# @2mm
# ra1+ra2+ra3
# ra1+ra2
# ra2+ra3
# ra1
# ra3


