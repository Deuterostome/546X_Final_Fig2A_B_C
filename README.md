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
GO_ra1_1 = filter(f2_nucleosome, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
# for ra2@1mm 
GO_ra2_1 = filter(f2_nucleosome, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
# for ra3@1mm 
GO_ra3_1 = filter(f2_nucleosome, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
# for ra1@2mm 
GO_ra1_2 = filter(f2_nucleosome, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
# for ra2@2mm 
GO_ra2_2 = filter(f2_nucleosome, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
# for ra3@2mm 
GO_ra3_2 = filter(f2_nucleosome, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")

# chromatin
chromatin = read.csv("chromatin.csv", header = TRUE)
chromatin = tbl_df(chromatin)
chromatin_df = data.frame(chromatin$maize.gene.id)
chromatin_1 = arrange(chromatin_df, chromatin.maize.gene.id)
f2_chromatin = merge(f2_df_1, chromatin_1, by.x = "gene_id", by.y = "chromatin.maize.gene.id")
# vesicle
vesicle = read.csv("vesicle.csv", header = TRUE)
vesicle = tbl_df(vesicle)
vesicle_df = data.frame(vesicle$maize.gene.id)
vesicle_1 = arrange(vesicle_df, vesicle.maize.gene.id)
f2_vesicle = merge(f2_df_1, vesicle_1, by.x = "gene_id", by.y = "vesicle.maize.gene.id")
# temperature
temperature = read.csv("temperature.csv", header = TRUE)
temperature = tbl_df(temperature)
temperature_df = data.frame(temperature$maize.gene.id)
temperature_1 = arrange(temperature_df, temperature.maize.gene.id)
f2_temperature = merge(f2_df_1, temperature_1, by.x = "gene_id", by.y = "temperature.maize.gene.id")
# stress
stress = read.csv("stress.csv", header = TRUE)
stress = tbl_df(stress)
stress_df = data.frame(stress$maize.gene.id)
stress_1 = arrange(stress_df, stress.maize.gene.id)
f2_stress = merge(f2_df_1, stress_1, by.x = "gene_id", by.y = "stress.maize.gene.id")
# transcription
transcription = read.csv("transcription.csv", header = TRUE)
transcription = tbl_df(transcription)
transcription_df = data.frame(transcription$maize.gene.id)
transcription_1 = arrange(transcription_df, transcription.maize.gene.id)
f2_transcription = merge(f2_df_1, transcription_1, by.x = "gene_id", by.y = "transcription.maize.gene.id")
# nitrogen
nitrogen = read.csv("nitrogen.csv", header = TRUE)
nitrogen = tbl_df(nitrogen)
nitrogen_df = data.frame(nitrogen$maize.gene.id)
nitrogen_1 = arrange(nitrogen_df, nitrogen.maize.gene.id)
f2_nitrogen = merge(f2_df_1, nitrogen_1, by.x = "gene_id", by.y = "nitrogen.maize.gene.id")
# RNA_synthetase
RNA_synthetase = read.csv("RNA_synthetase.csv", header = TRUE)
RNA_synthetase = tbl_df(RNA_synthetase)
RNA_synthetase_df = data.frame(RNA_synthetase$maize.gene.id)
RNA_synthetase_1 = arrange(RNA_synthetase_df, RNA_synthetase.maize.gene.id)
f2_RNA_synthetase = merge(f2_df_1, RNA_synthetase_1, by.x = "gene_id", by.y = "RNA_synthetase.maize.gene.id")
# G_protein
G_protein = read.csv("G_protein.csv", header = TRUE)
G_protein = tbl_df(G_protein)
G_protein_df = data.frame(G_protein$maize.gene.id)
G_protein_1 = arrange(G_protein_df, G_protein.maize.gene.id)
f2_G_protein = merge(f2_df_1, G_protein_1, by.x = "gene_id", by.y = "G_protein.maize.gene.id")
# light
light = read.csv("light.csv", header = TRUE)
light = tbl_df(light)
light_df = data.frame(light$maize.gene.id)
light_1 = arrange(light_df, light.maize.gene.id)
f2_light = merge(f2_df_1, light_1, by.x = "gene_id", by.y = "light.maize.gene.id")
# redox
redox = read.csv("redox.csv", header = TRUE)
redox = tbl_df(redox)
redox_df = data.frame(redox$maize.gene.id)
redox_1 = arrange(redox_df, redox.maize.gene.id)
f2_redox = merge(f2_df_1, redox_1, by.x = "gene_id", by.y = "redox.maize.gene.id")
# monosaccharide
monosaccharide = read.csv("monosaccharide.csv", header = TRUE)
monosaccharide = tbl_df(monosaccharide)
monosaccharide_df = data.frame(monosaccharide$maize.gene.id)
monosaccharide_1 = arrange(monosaccharide_df, monosaccharide.maize.gene.id)
f2_monosaccharide = merge(f2_df_1, monosaccharide_1, by.x = "gene_id", by.y = "monosaccharide.maize.gene.id")
# cell_wall
cell_wall = read.csv("cell_wall.csv", header = TRUE)
cell_wall = tbl_df(cell_wall)
cell_wall_df = data.frame(cell_wall$maize.gene.id)
cell_wall_1 = arrange(cell_wall_df, cell_wall.maize.gene.id)
f2_cell_wall = merge(f2_df_1, cell_wall_1, by.x = "gene_id", by.y = "cell_wall.maize.gene.id")
# sucrose
sucrose = read.csv("sucrose.csv", header = TRUE)
sucrose = tbl_df(sucrose)
sucrose_df = data.frame(sucrose$maize.gene.id)
sucrose_1 = arrange(sucrose_df, sucrose.maize.gene.id)
f2_sucrose = merge(f2_df_1, sucrose_1, by.x = "gene_id", by.y = "sucrose.maize.gene.id")
# trehalose
trehalose = read.csv("trehalose.csv", header = TRUE)
trehalose = tbl_df(trehalose)
trehalose_df = data.frame(trehalose$maize.gene.id)
trehalose_1 = arrange(trehalose_df, trehalose.maize.gene.id)
f2_trehalose = merge(f2_df_1, trehalose_1, by.x = "gene_id", by.y = "trehalose.maize.gene.id")
# intracellular_protein
intracellular_protein = read.csv("intracellular_protein.csv", header = TRUE)
intracellular_protein = tbl_df(intracellular_protein)
intracellular_protein_df = data.frame(intracellular_protein$maize.gene.id)
intracellular_protein_1 = arrange(intracellular_protein_df, intracellular_protein.maize.gene.id)
f2_intracellular_protein = merge(f2_df_1, intracellular_protein_1, by.x = "gene_id", by.y = "intracellular_protein.maize.gene.id")
# small_gtpase
small_gtpase = read.csv("small_gtpase.csv", header = TRUE)
small_gtpase = tbl_df(small_gtpase)
small_gtpase_df = data.frame(small_gtpase$maize.gene.id)
small_gtpase_1 = arrange(small_gtpase_df, small_gtpase.maize.gene.id)
f2_small_gtpase = merge(f2_df_1, small_gtpase_1, by.x = "gene_id", by.y = "small_gtpase.maize.gene.id")
# intracellular_signaling
intracellular_signaling = read.csv("intracellular_signaling.csv", header = TRUE)
intracellular_signaling = tbl_df(intracellular_signaling)
intracellular_signaling_df = data.frame(intracellular_signaling$maize.gene.id)
intracellular_signaling_1 = arrange(intracellular_signaling_df, intracellular_signaling.maize.gene.id)
f2_intracellular_signaling = merge(f2_df_1, intracellular_signaling_1, by.x = "gene_id", by.y = "intracellular_signaling.maize.gene.id")
# proteasome
proteasome = read.csv("proteasome.csv", header = TRUE)
proteasome = tbl_df(proteasome)
proteasome_df = data.frame(proteasome$maize.gene.id)
proteasome_1 = arrange(proteasome_df, proteasome.maize.gene.id)
f2_proteasome = merge(f2_df_1, proteasome_1, by.x = "gene_id", by.y = "proteasome.maize.gene.id")
# nucleocytoplasmic
nucleocytoplasmic = read.csv("nucleocytoplasmic.csv", header = TRUE)
nucleocytoplasmic = tbl_df(nucleocytoplasmic)
nucleocytoplasmic_df = data.frame(nucleocytoplasmic$maize.gene.id)
nucleocytoplasmic_1 = arrange(nucleocytoplasmic_df, nucleocytoplasmic.maize.gene.id)
f2_nucleocytoplasmic = merge(f2_df_1, nucleocytoplasmic_1, by.x = "gene_id", by.y = "nucleocytoplasmic.maize.gene.id")
# mitochondrial_transport
mitochondrial_transport = read.csv("mitochondrial_transport.csv", header = TRUE)
mitochondrial_transport = tbl_df(mitochondrial_transport)
mitochondrial_transport_df = data.frame(mitochondrial_transport$maize.gene.id)
mitochondrial_transport_1 = arrange(mitochondrial_transport_df, mitochondrial_transport.maize.gene.id)
f2_mitochondrial_transport = merge(f2_df_1, mitochondrial_transport_1, by.x = "gene_id", by.y = "mitochondrial_transport.maize.gene.id")
# ion_trans
ion_trans = read.csv("ion_trans.csv", header = TRUE)
ion_trans = tbl_df(ion_trans)
ion_trans_df = data.frame(ion_trans$maize.gene.id)
ion_trans_1 = arrange(ion_trans_df, ion_trans.maize.gene.id)
f2_ion_trans = merge(f2_df_1, ion_trans_1, by.x = "gene_id", by.y = "ion_trans.maize.gene.id")
# ATP_synthesis_coupled
ATP_synthesis_coupled = read.csv("ATP_synthesis_coupled.csv", header = TRUE)
ATP_synthesis_coupled = tbl_df(ATP_synthesis_coupled)
ATP_synthesis_coupled_df = data.frame(ATP_synthesis_coupled$maize.gene.id)
ATP_synthesis_coupled_1 = arrange(ATP_synthesis_coupled_df, ATP_synthesis_coupled.maize.gene.id)
f2_ATP_synthesis_coupled = merge(f2_df_1, ATP_synthesis_coupled_1, by.x = "gene_id", by.y = "ATP_synthesis_coupled.maize.gene.id")

# for ra1@1mm 
GO_ra1_1 = filter(f2_TF, series == "mut_series", size == "1mm", q2 == "ra1", significant == "yes")
# for ra2@1mm 
GO_ra2_1 = filter(f2_TF, series == "mut_series", size == "1mm", q2 == "ra2", significant == "yes")
# for ra3@1mm 
GO_ra3_1 = filter(f2_TF, series == "mut_series", size == "1mm", q2 == "ra3", significant == "yes")
# for ra1@2mm 
GO_ra1_2 = filter(f2_TF, series == "mut_series", size == "2mm", q2 == "ra1", significant == "yes")
# for ra2@2mm 
GO_ra2_2 = filter(f2_TF, series == "mut_series", size == "2mm", q2 == "ra2", significant == "yes")
# for ra3@2mm 
GO_ra3_2 = filter(f2_TF, series == "mut_series", size == "2mm", q2 == "ra3", significant == "yes")


