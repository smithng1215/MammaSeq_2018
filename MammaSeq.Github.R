#############################################################################################################
# set working directory to jamaOncology_2016
userDirectory <- "~/Desktop/MammaSeq_2018_Master"
setwd(userDirectory)

# prerequisite packages, be sure these are already installed using installRequiredPackages.R
requiredPackages.cran <- c(
  "reshape2",
  "tidyr",
  "ggplot2")

requiredPackages.bioconductor <- c(
  "ComplexHeatmap")

lapply(c(requiredPackages.bioconductor,requiredPackages.cran), require, character.only = TRUE)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

#############################################################################################################
###################################### Define Functions #####################################################

ColorScaleUP <- function(ColorValue) {
  if (ColorValue > 3 ) {
    return(1)
  }
}

ColorScaleDOWN <- function(ColorValue) {
  if (ColorValue < 1 ) {
    return(1)
  }
}

#############################################################################################################
###################################### Load Data Frame ######################################################
load("mamma_full_withChasm.RData")
load("cnv.data.Rdata")
load("clinical.data.Rdata")
load("oncokb.download.171006.Rdata")
load("cfdna.Rdata")

#############################################################################################################
######################################## snp Filtering ######################################################

mamma.filter.1=Mamma[Mamma$Position %in% names(which(table(Mamma$Position) < 10)), ]
mamma.filter.2=mamma.filter.1[-grep('SY',mamma.filter.1$Sequence.ontology),]
mamma.filter.3=mamma.filter.2[mamma.filter.2$STB >= 0.5 & mamma.filter.2$STB <= 0.6, ]
mamma.filter.4=mamma.filter.3[mamma.filter.3$ExAC.total.AF < 0.01 | mamma.filter.3$X1000.Genomes.AF < 0.01, ]
mamma.final=mamma.filter.4[mamma.filter.4$AF < 0.9, ]


##############################################################################################################
####################################### Figure 1 #############################################################

AMP=cnv[cnv$copy.number > 6 ,]
DEL=cnv[cnv$copy.number < 1,]

AMP$color.translation=sapply(AMP$copy.number,FUN=ColorScaleUP)
DEL$color.translation=sapply(DEL$copy.number,FUN=ColorScaleDOWN)
AMP <- dcast(AMP,gene~Sample.ID,value.var = "color.translation",fun.aggregate = sum)
DEL <- dcast(DEL,gene~Sample.ID,value.var = "color.translation",fun.aggregate = sum)
rownames(AMP)=AMP$gene
rownames(DEL)=DEL$gene
AMP=AMP[,-1]
DEL=DEL[,-1]
AMP=as.matrix(AMP)
DEL=as.matrix(DEL)

MUT=mamma.final[,c(5:7)]
MS=as.numeric(gsub("MS|CS|FI|ID|SS|SG|FD|II","1",MUT$Sequence.ontology))
MUT$Sequence.ontology=MS
MUT=MUT[-grep("KRAS|IDH1",MUT$HUGO.symbol),]
MUT$HUGO.symbol <- factor(MUT$HUGO.symbol)
MUT=dcast(MUT, HUGO.symbol~Sample.ID ,value.var = "Sequence.ontology",fun.aggregate = sum)
rownames(MUT)=MUT[,1]
MUT=MUT[,-1]
MUT=as.matrix(MUT)


alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  DEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "blue", col = NA))
  },
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "red", col = NA))
  },
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#008000", col = NA))
  }
)


mat.list <- list(MUT=MUT, AMP=AMP, DEL=DEL)
mat.list.uni <- unify_mat_list(mat.list)

col = c("MUT" = "#008000", "AMP" = "red", "DEL" = "blue")

tumors=clinical.data[,c(1,2,5)]
colnames(tumors)=c("Histology","Receptor Status","Tumor Type")
tumors$`Tumor Type`=gsub("PR","Primary",tumors$`Tumor Type`)
tumors$`Tumor Type`=gsub("MET","Metastasis",tumors$`Tumor Type`)
ha_column = HeatmapAnnotation(df = tumors, 
                              col = list( `Tumor Type` = c("Primary" =  "#000000", "Metastasis" = "#787878"),
                                          `Receptor Status` = c("ER+"= "#000080","ER+/HER2+"="#9370DB","HER2+"="#DB7093","Triple Negative"="#ADD8E6","N/A"="#CCCCCC"),
                                           Histology = c("IDC"="#FF6347","ILC"="#40E0D0","IDC/ILC"="#EE82EE","N/A"="#CCCCCC"))) 

tumor.ht=oncoPrint(mat.list.uni, get_type = function(x) strsplit(x, ";")[[1]],
          bottom_annotation = ha_column,
          alter_fun = alter_fun, col = col, 
          column_title = NULL,
          remove_empty_columns = T,
          heatmap_legend_param = list(title = "Alternation", at = c("AMP", "DEL", "MUT"), 
                                      labels = c("Amplification", "Deletion", "Mutation")))


pdf("figure1.oncoprint.pdf",height = 14,width=6)
draw(tumor.ht,heatmap_legend_side = "bottom")
dev.off()


##############################################################################################################
######################################## Overlap with oncoKB database ########################################

oncokb$combined=as.character(interaction(oncokb$Gene, oncokb$Alteration))
mamma.final$combined=as.character(interaction(mamma.final$HUGO.symbol,mamma.final$Protein.sequence.change))
overlap=mamma.final[complete.cases(match(mamma.final$combined,oncokb$combined)),]
overlap=merge(overlap,oncokb,by="combined")
oncokb.table=as.data.frame(cbind( as.character(overlap$Sample.ID),
                                  as.character(overlap$HUGO.symbol),
                                  as.character(overlap$Sample.ID.x),
                                  as.character(overlap$Gene),
                                  as.character(overlap$Protein.sequence.change),
                                  as.character(overlap$Level),
                                  as.character(overlap$AF)))

oncokb.table=oncokb.table[order(oncokb.table$V5),] # Mutations for table 3

#############################################################################################################
#################################### cfdna snp Filtering ####################################################

cfdna.filter.1=cfdna[-grep('SY',cfdna$Sequence.ontology),]
cfdna.filter.2=cfdna.filter.1[cfdna.filter.1$Position %in% names(which(table(cfdna.filter.1$Position) < 4)), ]
cfdna.filter.3=cfdna.filter.2[cfdna.filter.2$ExAC.total.allele.frequency < 0.01 | cfdna.filter.2$X1000.Genomes.allele.frequency < 0.01, ]
cfdna.final=cfdna.filter.3[cfdna.filter.3$STB >= 0.5 & cfdna.filter.3$STB <= 0.6, ]

############################################################################################################
#################################### cfdna overlap with oncoKB #############################################

cfdna$combined=as.character(interaction(cfdna$HUGO.symbol,cfdna$Protein.sequence.change))
cfdna.overlap=cfdna[complete.cases(match(cfdna$combined,oncokb$combined)),]
cfdna.overlap=merge(cfdna.overlap,oncokbMS,by="combined")
cfdna.overlap.final=as.data.frame(cbind(as.character(cfdna.overlap$HUGO.symbol),as.character(cfdna.overlap$Sample.ID),
                           as.character(cfdna.overlap$Protein.sequence.change),
                           as.character(cfdna.overlap$Level.x),
                           as.character(cfdna.overlap$AF)))

############################################################################################################
##################################### Figure 3 oncoprint ###################################################

Mutation=cfdna.final[,c(5:7)]
cfdna.MS=as.numeric(gsub("MS|CS|FI|ID|SS|SG|SL|FD|II","1",Mutation$Sequence.ontology))
Mutation$Sequence.ontology=cfdna.MS
Mutation=dcast(Mutation, HUGO.symbol~Sample.ID ,value.var = "Sequence.ontology",fun.aggregate = sum)
rownames(Mutation)=Mutation[,1]
Mutation=Mutation[,-1]

no_mut=cbind(rep(0,29),rep(0,29),rep(0,29),rep(0,29))
colnames(no_mut)=c("CF_28_Draw_4","CF_23_Draw_2","CF_20_Draw_2","CF_26")
Mutation=cbind(Mutation,no_mut)
Mutation=as.matrix(Mutation)

cfdna.overlap.final$V5=rep(1,4)
cfdna.oncoKB.mut=dcast(cfdna.overlap.final,V1~V2,value.var = 'V5',fun.aggregate = sum)
rownames(cfdna.oncoKB.mut)=cfdna.oncoKB.mut$V1
cfdna.oncoKB.mut=as.matrix(cfdna.oncoKB.mut[,2:5])

alter_fun_list=list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "navy", col = NA))
  },
  cfdna.oncoKB.mut = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "purple", col = NA))
  })

mat.list <- list(Mutation = Mutation,cfdna.oncoKB.mut=cfdna.oncoKB.mut)
mat.list.uni <- unify_mat_list(mat.list)
order=colnames(Mutation)
order=order[c(1,2,13,3,4,5,12,14,6,7,8,9,11,10)]
col = c(Mutation = "navy",cfdna.oncoKB.mut="purple")

cfdna.ht= oncoPrint(mat.list.uni, get_type = function(x) strsplit(x, ";")[[1]], alter_fun = alter_fun_list, col = col,
              show_column_names = TRUE,
              remove_empty_columns = FALSE,
              column_order = order,
              row_names_gp = gpar(fontsize = 11),
              heatmap_legend_param = list(title = "Alternation", at = c("Mutation", "cfdna.oncoKB.mut"), 
                                          labels = c("Non-actionable Mutation", "Level 3 Annotated Mutation"), labels_gp = gpar(fontsize = 11),title_gp = gpar(fontsize = 11)))



pdf("figure3.oncoprint.pdf",height=7.5,width=4.5)
draw(cfdna.ht, heatmap_legend_side = "bottom")
dev.off()





