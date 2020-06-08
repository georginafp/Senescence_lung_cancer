########### DADES DEL TCGA PROJECTE - LUAD ##########



#Adjudico el meu directori de treball
#setwd("C:/Users/Georgina/Documents")
          setwd("C:/Users/Georgina/Desktop/LUAD")
getwd()
#Crido la llibreria de TCGAbiolinks per tal de fer els an?lisis pertinents


#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
# BiocManager::install(version = "3.11")
# 
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("TCGAbiolinks")
# 


#install.packages("installr")
#library(installr)
#updateR()

library(TCGAbiolinks)
#Descarrego les dades cliniques gracies a la funci? GCDquery_clinic
#clinical0 <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")
#Carrego un fitxer R que cont? les dades cliniques, que est? al mateix directori
#load("clinical0_LUSC.Rdata")
load("clinical0_LUAD.Rdata")
#Aplico la funci? summary per a veure un resum del contingut de la variable clinical0 (Minims, Maxims, llargada etc...)
summary(clinical0)

#Converteixo la variable clinical a factor amb la funci? as.factor i la funci? apply per a aplicar-ho a la funci? clinical0.
#El n?mero 2 fa refer?ncia a que ho apliquem a les columnes
clinical <- apply(clinical0, 2, as.factor)
#Converteixo la variable que pr?viament he convertit a factor, en un dataframe
clinical<-data.frame(clinical)
#Aplico la funci? summary altre vegada per tal de veure un resum del contingut del nou dataframe
summary(clinical)

#Elimino la segona columna de cada individu ja que estan duplicats per a cada tractament (quimioter?pia i radioter?pia) perqu? la segona fila est? buida
idna<-which(apply(clinical, 2, function(x) sum(is.na(x)))==1044)
idnarows<-which(apply(clinical, 1, function(x) sum(is.na(x)))==56)
idnarows
#Quan aplico aquesta ja no els tinc duplicats, ja no tinc 1008 si no 522 (que s?n els individus totals en les BDD del TCGA) perque els he eliminat
clinical<-clinical[-idnarows,-idna]
dim(clinical)

clinical<-data.frame(clinical)

#Creem una iteraci? per a que per cada individu (que en tenim 522), em concateni els valors de la variable "treatment or therapy" per tal de:
#Tenir per cada individu si ha fet Quimioterapia i/o Radioterapia
#Exemple:    Individu 1 --> yes yes voldr? dir que ha fet tant quimioterapia com radioterapia
for(i in 1:522){
  clinical$treatment[i]<-paste(clinical0$treatment_or_therapy[(2*i-1)],clinical0$treatment_or_therapy[(2*i)])
}

#Convertim la variable treatment a factor
clinical$treatment<-factor(clinical$treatment)
#Aplico la funci? summary per tal de veure quants individus hi ha de cada condici?.
#?s a dir, quants han fet tant quimioter?pia com radioter?pia, quants nom?s quimioter?pia, quants nom?s radioter?pia i quants no han rebut tractament
summary(clinical$treatment)
#Faig un summary de la variable anterior per a veure si els n?meros totals quadren amb els nous adjudicats a la funci? treatment (ja concatenats)
summary(clinical$treatment_or_therapy)

#Elimino totes aquelles columnes de la variabel clinical que la seva dimensi? sigui 1
idunic<-which(apply(clinical,2,function(x) dim(table(x))==1))
clinical<-clinical[,-idunic]
summary(clinical)

#Estableixo que la variable clinical$chemo es la mateixa que clinical$treatment_or_therapy original
clinical$chemo<-clinical$treatment_or_therapy
#Visualitzo que es correcte
clinical$chemo
#Faig una taula per veure els valors de treatment i de la variable creada com a chemo
table(clinical$treatment)
table(clinical$chemo)

#Replico el valor 0 a tots els 522 casos
clinical$chemo_or_radio<-rep(0,522)
clinical$chemo_or_radio

#Creo la variable "id" per a que em retorni un vector amb les posicions dels pacients que s'assumeixen com a "no not reported"
id<-which(clinical$treatment=="no not reported")
id

#Ara assigno un 1 a aquells pacients que s'assumien com a "no not reported"
clinical$chemo_or_radio[id]<- 1
clinical$chemo_or_radio[id]

#Creo la variable "id" per a que em retorni un vector amb les posicions dels pacients que s'assumeixen com a "not reported no"
id<-which(clinical$treatment=="not reported no")
id
#Ara assigno un 1 a aquells pacients que s'assumien com a "not reported no"
#clinical$chemo_or_radio[id]<- 1

#Creo la variable "id" per a que em retorni un vector amb les posicions dels pacients que s'assumeixen com a "not reported not reported"
id<-which(clinical$treatment=="not reported not reported")
#Ara assigno un 1 a aquells pacients que s'assumien com a "not reported not reported"
clinical$chemo_or_radio[id]<- 1

#Creo la variable "id" per a que em retorni un vector amb les posicions dels pacients que s'assumeixen com nomes han fet radioterapia
id<-which(clinical$treatment=="no yes")
#Li assigno un 2 a aquells que nom?s han fet radioterapia
clinical$chemo_or_radio[id]<- 2

#Creo la variable "id" per a que em retorni un vector amb les posicions dels pacients que s'assumeixen com nomes han fet radioterapia i no s'ha reportat si tambe quimioterapia
id<-which(clinical$treatment=="not reported yes")
#Li assigno un 2 a aquells que nom?s han fet radioterapia i que no s'ha reportat si tambe quimioterapia
#clinical$chemo_or_radio[id]<- 2

#Creo la variable "id" per a que em retorni un vector amb les posicions dels pacients que s'assumeixen com nomes han fet quimioterapia i no s'ha reportat si tambe radioterapia
id<-which(clinical$treatment=="yes not reported")
#Li assigno un 2 a aquells que nom?s han fet quimioterapia i que no s'ha reportat si tambe radioterapia
#clinical$chemo_or_radio[id]<- 2

#Creo la variable "id" per a que em retorni un vector amb les posicions dels pacients que s'assumeixen com nomes han fet quimioterapia i no radioterapia
id<-which(clinical$treatment=="yes no")
#Li assigno un 2 a aquells que nom?s han fet quimioterapia i  no radioterapia
clinical$chemo_or_radio[id]<- 2

#Creo la variable "id" per a que em retorni un vector amb les posicions dels pacients que s'assumeixen com nomes han fet quimioterapia i radioterapia
id<-which(clinical$treatment=="yes yes")
#Li assigno un 2 a aquells que nom?s han fet quimioterapia i radioterapia
clinical$chemo_or_radio[id]<- 2

#Converteixo la variable "chemo_or_ratio" a factor amb els nivells 0,1,2 corresponents a les condicions "no", "not reported" i "yes" respectivament
clinical$chemo_or_radio<-factor(clinical$chemo_or_radio, levels=c(0,1,2), labels=c("no", "not reported", "yes"))
clinical$chemo_or_radio
#Aplico la funci? summary per a veure les dimensions de cada condicio tant de la variable "chemo_or_ratio" com de la variable "chemo"
summary(clinical$chemo_or_radio)
summary(clinical$chemo)

#Demano l'identificaci? d'aquells pacients que SI han fet quimioterapia
namechemoyes<-clinical$submitter_id[clinical$chemo=="yes"]
namechemoyes
#Aplico la funcio lenght per a saber quants n'han fet
length(namechemoyes)  # 109 han fet quimio

#Demano l'identificaci? d'aquells pacients que NO han fet quimioterapia
namechemono<-clinical$submitter_id[clinical$chemo=="no"]
namechemono
#Aplico la funcio lenght per a saber quants no n'han fet
length(namechemono)  # 367 no


#PROCEDIM A L'ANALISI D'EXPRESSI?
#Descarreguem les dades d'expressio filtrant primer el projecte de TCGA
#Que siguin d'expressio
#legacy significa que ja ha estat alineat a hg19
#Que siguin dades de quantificacio
#Que la plataforma utilitzada sigui Illumina
#I que les dades siguin de tumors primaris ja que l'analisi hauria de ser previ a un desenvolupament de met?stasis per exemple
query <- GDCquery(project = "TCGA-LUAD",
                  legacy=TRUE,
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  platform = "Illumina HiSeq",
                  file.type = "results")
#Descarregem les dades
GDCdownload(query)

#Necesstitem fer la preparacio per a poder treballar amb objectes d'R i aix? es el que ens fa la funcio GDCprepare i ho guardem en un fitxer en format .rda
#Preparem les dades en forma de SummarizedExperiment, com una agrupacio
luad.exp <- GDCprepare(query,
                       save = TRUE,
                       summarizedExperiment = TRUE,
                       save.filename = "LUADIllumina_HiSeq.rda")

#Carreguem el fitxer previament guardat en el proc?s de preparaci?
rse <- luad.exp
#Veig que sí que es un summarizedExperiment
class(rse)


#PRE PROCESSAMENT DE LES DADES
#Hem de pre-processar les dades per tal de posteriorment normalitzar-les i ho guardem en un fitxer en format .png
dataPrep_LUAD <- TCGAanalyze_Preprocessing(object = rse,
                                           cor.cut = 0.6,
                                           datatype = "raw_count",
                                           filename = "LUAD_IlluminaHiSeq_RNASeqV2.png")


#NORMALITZACIO DE LES DADES
#volem normalitzar-les per tal de poder-les portar a un nivell comparable
#ho fem a partir del m?tode gcContent
dataNorm <- TCGAanalyze_Normalization(tabDF = cbind(dataPrep_LUAD),
                                      geneInfo = geneInfo,
                                      method = "gcContent")

#Aplico la funci? boxplot per a veure els outliers que hi havia un cop hem preparat les dades per comparar-ho amb un boxplot despr?s de normalitzar-les
boxplot(dataPrep_LUAD, outline = FALSE)
boxplot(dataNorm, outline = FALSE)



#Necessitem filtrar les dades per tal d'obtenir uns resultats significatius
#Ho fem a partir del m?tode de quantils
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile",
                                  qnt.cut = 0.25)
#Ho guardem en un fitxer de format .rda
save(dataFilt, file = paste0("LUAD_Norm_IlluminaHiSeq.rda"))
class(dataFilt)

#carrego el fitxer que hem creat quan hem preparat les dades
load("LUADIllumina_HiSeq.rda")


#Assigno un valor nul a la variable idchemo per a posteriorment recorrer tots els noms de la variable namechemo que hem creat a la l?nia 122
#Per tal de a partir d'una iteraci? en loop recorre'ls dins els que ja estaven filtrats i processats

#ho fem per la variable idchemoNO --> i per tant els que no han fet quimioterapia
idchemono<-NULL
for (name in namechemono){
  idchemono<-c(idchemono,grep(name, colnames(dataFilt)))
}

#ho fem per la variable idchemoYES --> i per tant els que SI han fet quimioterapia
idchemoyes<-NULL
for (name in namechemoyes){
  idchemoyes<-c(idchemoyes,grep(name, colnames(dataFilt)))
}

#A partir de la funci? TCGAanalyze_DEA obtinc els gens diferencialment expressats entre les dues condicions:
#Els tractats amb quimioterapia
#Els no tractats
dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,idchemoyes],
                            mat2 = dataFilt[,idchemono],
                            method="glmLRT",
                            Cond1type = "Chemo",
                            Cond2type = "No-Chemo")


??TCGAanalyze_DEA

library(gridExtra)
pdf("DEGs.pdf", height=11, width=10)
grid.table(dataDEGs)
dev.off()



dataDEGs
colnames(dataDEGs)
dataDEGsFilt <- dataDEGs[abs(dataDEGs$logFC) >= 1,]
datano <- dataFilt[,idchemono]
datayes <- dataFilt[,idchemoyes]
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGsFilt,"chemo","no-chemo",
                                          datayes,datano)


#com me l'ha enviat la Malu
plot(dataDEGs$logFC,-log(dataDEGs$PValue))
?plot

plotdata <- as.matrix(dataDEGs)
#MA PLOT
library(limma)
plotMA(dataDEGs, ylim = c(-4, 3) , xlim= c(-1,3), main="MA plot from DEGs")




#HISTOGRAMA
hist(dataDEGs$PValue, breaks=20, col="red")


#ANEM A VEURE SI EL QUE ESTÀ SOBRE-EXPRESSAT ENS SURT EN ALGUN DELS PATHWAYS RELACIONATS AMB SENESCÈNCIA:
GenelistComplete<- rownames(assay(rse,1))
dataDEGsFiltLevel$GeneID<-0

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("clusterProfiler")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("org.Hs.eg.db")

library(org.Hs.eg.db)
library(clusterProfiler)
eg=as.data.frame(bitr(dataDEGsFiltLevel$mRNA,
                         fromType="SYMBOL",
                         toType="ENTREZID",
                         OrgDb="org.Hs.eg.db"))

eg<-eg[!duplicated(eg$SYMBOL),]
dataDEGsFiltLevel<-dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA %in% eg$SYMBOL,]
dataDEGsFiltLevel<-dataDEGsFiltLevel[order(dataDEGsFiltLevel$mRNA,decreasing=FALSE),]
eg<-eg[order(eg$SYMBOL,decreasing=FALSE),]
#Hauria de retornar-nos TRUE:
table(eg$SYMBOL == dataDEGsFiltLevel$mRNA)   
all(eg$SYMBOL == dataDEGsFiltLevel$mRNA)
dataDEGsFiltLevel$GeneID<-eg$ENTREZID
dataDEGsFiltLevel_sub<-subset(dataDEGsFiltLevel, select = c("GeneID", "logFC"))
genelistDEGs<-as.numeric(dataDEGsFiltLevel_sub$logFC)
names(genelistDEGs)<-dataDEGsFiltLevel_sub$GeneID



#### extreure la taula de DEGsFilt (198 gens)


library(ggplot2)
library(gridExtra)

df <- dataDEGsFilt
pvalue <- dataDEGsFilt$PValue
index10<-order(pvalue, decreasing = F)
#Ordeno la meva taula creada amb noms de gens i p-valor de manera ordenada

gens_ordenats<-df[index10,]


png("test.png", height=15000, width=21000)
p<-tableGrob(gens_ordenats)
grid.arrange(p)
dev.off()



#BiocManager::install("pathview")
require("pathview")
library(pathview)

#Com que volem saber si hi ha relació amb el TFG-B en podem mirar el seu pathway 
#Red defines genes that are up-regulated and green defines genes that are down-regulated.

#hsa04350 --> TGF-B signalling pathway id in KEGG
hsa04350<-pathview(gene.data  = genelistDEGs,
                   pathway.id = "hsa04350",
                   species    = "hsa",
                   limit      = list(gene=as.integer(max(abs(genelistDEGs)))))
                
#hsa04218  --> cellular senescence 
hsa04218<-pathview(gene.data  = genelistDEGs,
                    pathway.id = "hsa04218",
                    species    = "hsa",
                    limit      = list(gene=as.integer(max(abs(genelistDEGs)))))



#NO VEIG RES NI SOBRE-EXPRESSAT NI INFRA-EXPRESSAT --> surten les mateixes imatges de referencia que al KEGG 


#Procedim a l'analisi d'enriquiment per un conjunt de gens (GeneList) gracies a la funcio TCGAanalyze_EAcomplete
#Realitzem un barplot gracies a la funcio TCGAvisualize_EAbarplot pero nomes en mostro 20 per a fer-ho mes visual
#Representem el proces biologic, el component celular, la funci? molecular i el pathway
#Se'm genera automaticament un PDF al directori de treball amb els 4 gr?fics
ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes chemo Vs non-chemo",
                                RegulonList = rownames(dataDEGs))
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),
                        GOBPTab = ansEA$ResBP,
                        GOCCTab = ansEA$ResCC,
                        GOMFTab = ansEA$ResMF,
                        PathTab = ansEA$ResPat,
                        nRGTab = rownames(dataDEGs),
                        nBar = 20)

getwd()
patient <- clinical$bcr_patient_barcode
patient
class(dataDEGsFiltLevel_sub)



# VOLCANOPLOT
load("C:/Users/Georgina/Desktop/LUAD/resultatsLUAD.Rdata")
with(dataDEGs, plot(logFC, -log10(PValue), pch=20, main="Volcano plot DEGs LUAD", xlim=c(-7,7)))


#si el p-valor adj<0.05
with(subset(dataDEGs, FDR<.05 ), points(logFC, -log10(PValue), pch=20, col="red"))
#si logFC>1
with(subset(dataDEGs, abs(logFC)>1), points(logFC, -log10(PValue), pch=20, col="orange"))
#si p-valor adj<0.05 i logFC>1
with(subset(dataDEGs, FDR<.05 & abs(logFC)>1), points(logFC, -log10(PValue), pch=20, col="green"))
legend(-7, 33, legend=c("Adjusted p-value < 0.05", "logFC > 1", "Adjusted p-value < 0.05 and logFC > 1"),
       col=c("red", "orange", "green"), lty=1:2, cex=0.8)




####HEATMAP
dim(dataFilt)
index <- which(dataDEGsFilt$PValue<0.0000000000001)
dataHeatMap <- dataFilt[index,union(idchemono,idchemoyes)]
dim(dataHeatMap)
colnames(dataHeatMap)


library(RColorBrewer)
#coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)
#heatmap(dataHeatMap, scale="row", cexRow=1.5, col= colorRampPalette(brewer.pal(8, "Blues"))(25))

mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
library("gplots")
heatmap(as.matrix(dataHeatMap), scale="row")
heatmap.2(as.matrix(dataHeatMap),col=rev(morecols(50)),trace="none", main="HeatMap LUAD project",scale="row")






clust.rows <- hclust(as.dist(1-cor(t(dataHeatMap))),method="ward.D2")
clust.cols <- hclust(as.dist(1-cor(dataHeatMap)),method="ward.D2")  


#Estableixo una paleta de colors
heatcol<-colorRampPalette(c("green", "Black","red"), space = "rgb")

# #El fem per columnes
# #Crido la funció utilizant la matriu "data.clus" guardada (extreta de les dades normalitzades) 
# #De la paleta de color n'extrec 256 colors (del vermell al verd passaant pel negre)
heatm<-heatmap.2(as.matrix(dataHeatMap), col = heatcol(256),
                 dendrogram="column", Colv=as.dendrogram(clust.cols),    #dendrogrames 
                 Rowv=as.dendrogram(clust.rows),     #escalem per files 
                 scale="row",cexRow=0.3, cexCol=0.5,    
                 main="",key=TRUE,keysize=1,density.info="none",trace="none")
# #Això ho fem per guardar una imatge que dona R. 








####                      GSEA  --> SENESCENCE SIGNATURE
gene_name<-rownames(dataDEGs)
pvalue<- dataDEGs$PValue

#Faig una taula amb aquestes dues columnes
table1<-data.frame(gene_name,pvalue)
#Els ordeno en creixent per a tenir els més significatius (p-valor més baix) a dalt de tot en una variable de nom index
index<-order(pvalue, decreasing = F)
#Ordeno la meva taula creada amb noms de gens i p-valor de manera ordenada
ordered_table<-table1[index,]
#Afegeixo una columna més que es per la signatura de senescencia però nomes la omplo amb 0s
ordered_table$senescence_signature<-rep(0,14899)


#Vector amb els noms dels gens de la signatura de senescencia
senescence_vector=c("C9orf3","CTGF", "IFNG", "IGFB2", "IGFB3","ALDH1A3", "AOPEP", "CCN2", "CCND1", "CD44", "CDKN1A", "CDKN1C", "CDKN2A", "CDKN2B", "CDKN2D", "CITED2", "CLTB", "COL1A2", "CREG1", "CRYAB", "CXCL14", "CYP1B1", "EIF2S2", "ESM1", "F3", "FILIP1L", "FN1", "GSN", "GUK1", "HBS1L", "HPS5", "HSPA2", "HTATIP2", "IFI16", "IFNG", "IGFBP1", "IGFB2", "IGFB3", "IGFBP4", "IGFBP5", "IGFBP6", "IGFBP7", "IGSF3", "ING1", "IRF5", "IRF7", "ISG15", "MAP1LC3B", "MAP2K3", "MDM2", "MMP1", "NDN", "NME2", "NRG1", "OPTN", "PEA15", "RAB13", "RAB31", "RAB5B", "RABGGTA", "RAC1", "RBL2", "RGL2", "RHOB", "RRAS", "S100A11", "SERPINB2", "SERPINE1", "SMPD1", "SMURF2", "SOD1", "SPARC", "STAT1", "TES", "TFAP2A", "TGFB1I1", "THBS1", "TNFAIP2", "TNFAIP3", "TP53", "TSPYL5", "VIM")
length(senescence_vector)
#vector amb els gens diferencialment expressats meus 
#Que son els que vull veure si coincideixen (o estan dins -la majoria) amb els del vector de la signatura de senescència 
gene_name<- c(gene_name)

#Li poso 1s a tots aquells gens diferencialment expressats que coincideixen amb els del vector de la signatura de senescencia 
for (s in senescence_vector){
  gene_exists <- grep(s, gene_name)
  if (length(gene_exists) > 0){
    ordered_table[ordered_table$gene_name == s, "senescence_signature"] <- 1
  }
  
  
}


#Al denominar-se de diferent manera segons la base de dades, he hagut d'accedir un per un amb els diferents gens que no es denominem igual als dos vectors

# #AOPEP --> C9orf3
# ordered_table$gene_name[2045] 
# ordered_table$senescence_signature[2045] <- 1
# 
# #CCN2 --> CTGF
# ordered_table$gene_name[3267] 
# ordered_table$senescence_signature[3267] <- 1
# 
# #IFNG --> no hi és
# 
# ordered_table$gene_name[14174] 
# ordered_table$senescence_signature[14174] <- 1
# 
# #IGFB2 --> no té cap 1 però l'hauria de tenir 
# ordered_table$gene_name[2648] 
# ordered_table$senescence_signature[2648] <- 1
# 
# #IGFB3 --> no té cap 2 però l'hauria de tenir 
# ordered_table$gene_name[1364] 
# ordered_table$senescence_signature[1364] <- 1



#Hypergeometric test from an ordered table
pvalgsea<-NULL
x<-seq(1,1000,by=10)
m<-sum(ordered_table$senescence_signature) # m = nombre total de senesc
n<-nrow(ordered_table)-m # n= nombre total no de senesc
for (i in x){
  (q<-sum(ordered_table$senescence_signature[(1:i)]))  # q = nombre de gens senesc entre els top 100
  pval<-phyper(q-1, m, n, i, lower.tail = FALSE)
  pvalgsea<-c(pvalgsea,pval)
}

y<-pvalgsea  # p-values sense log
plot(x, y, type="o", col="red", main= "GSEA for senescence signature", xlab="ordered genes", ylab="p-value")
#no hi ha enriquiment perque ho hem graficat segons el p-valor
#Com més petit més significatiu --> els "top genes" (els 100 primers) no tenen els p-valors més petits
#per això no hi ha enriquiment





####                      GSEA  --> TFG SIGNATURE
gene_name<-rownames(dataDEGs)
pvalue<- dataDEGs$PValue

#Faig una taula amb aquestes dues columnes
table2<-data.frame(gene_name,pvalue)
#Els ordeno en creixent per a tenir els més significatius (p-valor més baix) a dalt de tot en una variable de nom index
index1<-order(pvalue, decreasing = F)
#Ordeno la meva taula creada amb noms de gens i p-valor de manera ordenada
ordered_table<-table2[index1,]
#Afegeixo una columna més que es per la signatura de senescencia però nomes la omplo amb 0s
ordered_table$tgf_signature<-rep(0,14899)




#Vector amb els noms dels gens de la signatura de senescencia
tgf_vector=c("ACVR1", "ACVR1B", "ACVR1C", "ACVR2A", "ACVR2B", "ACVRL1", "AMHR2", "ATF2", "BAMBI", "BMP1", "BMP10", "BMP15", "BMP2", "BMP3", "BMP4","BMP5", "BMP6", "BMP7", "BMP8A", "BMP8B", "BMPR1A", "BMPR1B", "BMPR2", "CITED1", "CITED2", "CREBBP", "DCP1B", "EP300", "FKBP1A", "FOSL1", "FOXH1", "GDF1", "GDF10", "GDF11", "GDF15", "GDF2", "GDF3", "GDF5", "GDF6", "GDF7", "GDF9", "GDNF", "HRAS", "INHBA", "INHBB", "INHBC", "INHBE", "JUN", "JUNB", "JUND", "LEFTY1", "LEFTY2", "MAP3K7", "MAP3K7CL", "MAPK1", "MAPK10", "MAPK11", "MAPK12", "MAPK13", "MAPK14", "MAPK3", "MAPK8", "MAPK9", "MSTN", "NODAL", "NRAS", "RRAS", "SKI", "SKIL", "SMAD1", "SMAD2", "SMAD3", "SMAD4", "SMAD5", "SMAD6", "SMAD7", "SMAD9", "SMURF1", "SMURF2", "SNIP1", "TAB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "TLL1", "TLL2", "ZFYVE9")
length(tgf_vector)
#vector amb els gens diferencialment expressats meus 
#Que son els que vull veure si coincideixen (o estan dins -la majoria) amb els del vector de la signatura de senescència 
#gene_name1<- c(gene_name1)



#Li poso 1s a tots aquells gens diferencialment expressats que coincideixen amb els del vector de la signatura de senescencia 
for (p in tgf_vector){
  gene_exist <- grep(p, gene_name)
  if (length(gene_exist) > 0){
    ordered_table[ordered_table$gene_name == p, "tgf_signature"] <- 1
  }
  
  
}





#Denominació diferent --> Hem de buscar noms alterantius (a GeneCards per exemple) i tornar-los a entrar com a 1s a ordered_table1:

#AMHR2 --> MRI1 no hi es 
#BMP10 --> no hi és 
#BMP15 --> no hi és 
#CITED1 --> no hi és 
#FOXH1 --> no hi és 
#GDF1 --> no hi és 
#GDF2 --> no hi és 
#GDF3 --> no hi és 
#GDF6 --> SGMS1 no hi és
#GDF7 --> no hi és 
#GDF9 --> no hi és 
#GDNF --> ATF2  ja està contat
#INHBC --> no hi és 
#INHBE --> no hi és 
#MAP3K7CL --> c21orf7 no hi és 
#NODAL --> no hi és 
#TLL2 --> no hi és 



#Hypergeometric test from an ordered table
pvalgsea<-NULL
x<-seq(1,1000,by=10)
m<-sum(ordered_table$tgf_signature) # m = nombre total de senesc
m
n<-nrow(ordered_table)-m # n= nombre total no de senesc
for (i in x){
  (q<-sum(ordered_table$tgf_signature[(1:i)]))  # q = nombre de gens senesc entre els top 100
  pval<-phyper(q-1, m, n, i, lower.tail = FALSE)
  pvalgsea<-c(pvalgsea,pval)
}


length(x)
length(y)
y<-pvalgsea  # p-values sense log
plot(x, y, type="o", col="blue", main= "GSEA 4 TGF signature", xlab="ordered genes", ylab="p-value")
#no hi ha enriquiment perque ho hem graficat segons el p-valor
#Com més petit més significatiu --> els "top genes" (els 100 primers) no tenen els p-valors més petits
#per això no hi ha enriquiment







rm(list=ls()) #will remove ALL objects 


#             ANALISIS DE SUPERVIVENCIA

#install.packages("survival")
#install.packages("KMsurv")
library(survival)
library(KMsurv)


#Descarrego i creo un data frame només amb la informacio que m'interessa pels analisis de supervivencia
#ara nomes volem veure diferencies entre chemo i non-chemo

#Descarrego nomes la informacio que minteresa

clin_df = clinical[,c("submitter_id",
                      "vital_status",
                      "days_to_last_follow_up",
                      "days_to_death",
                      "chemo")]
idDead<-which(clin_df$vital_status=="Dead")
idAlive<-which(clin_df$vital_status=="Alive")

clin_df$time<-rep(0,522)
for (i in idDead){
  clin_df$time[i]<-as.numeric(as.character(clinical$days_to_death[i]))
}
for (i in idAlive){
  clin_df$time[i]<-as.numeric(as.character(clinical$days_to_last_follow_up[i]))
}


#Discerneixo entre vius i morts (poso 1 a tots els morts)
clin_df$status<-rep(0,522)
clin_df[clin_df$vital_status == "Dead", "status"] <- 1

#Discerneixo entre tractats i no tractats (1 a tots els tractats i 2 als no-tractats)
clin_df$treat<-rep(0,522)
clin_df[clin_df$chemo == "yes", "treat"] <- 1
clin_df[clin_df$chemo == "no", "treat"] <- 2
clin_df[clin_df$chemo == "not reported", "treat"] <- NA

#Elimino els valors NAs (entre ells els no-reported)
data <- clin_df
library(tidyr)
data%>% drop_na()

#Faig l'analisi de supervivencia i amb el seu grafic de Kaplan-Meier
kmbytreatment<-survfit(Surv(data$time, data$status)~ data$treat)
summary(kmbytreatment)
plot(kmbytreatment, col=2:3, main="Kaplan-Meier LUAD project", ylab= "Survival probability", xlab="timeline (days)")
legend("topright",col=2:3, legend=c("chemo","non-chemo"), lty=1)

#Test estadistic per determinar si els diferencies entre les dues corbes son significatives
survdiff(Surv(data$time, data$status)~ data$treat)




#ara volem veure les diferencies en la supervivencia dels pacients amb diferents estadis del tumor
clin_df1 = clinical[,c("submitter_id",
                      "vital_status",
                      "days_to_last_follow_up",
                      "days_to_death",
                      "tumor_stage")]
idDead<-which(clin_df1$vital_status=="Dead")
idAlive<-which(clin_df1$vital_status=="Alive")

clin_df1$time<-rep(0,522)
for (i in idDead){
  clin_df1$time[i]<-as.numeric(as.character(clinical$days_to_death[i]))
}
for (i in idAlive){
  clin_df1$time[i]<-as.numeric(as.character(clinical$days_to_last_follow_up[i]))
}





#Li poso 1s a tot el que sigui Alive
clin_df1$status<-rep(0,522)
clin_df1[clin_df1$vital_status == "Dead", "status"] <- 1


clin_df1$stage<-rep(0,522)
clin_df1$stage

clin_df1[clin_df1$tumor_stage == "not reported", "stage"] <- NA
clin_df1[clin_df1$tumor_stage == "stage ia", "stage"] <- 1
clin_df1[clin_df1$tumor_stage == "stage i", "stage"] <- 1
clin_df1[clin_df1$tumor_stage == "stage ib", "stage"] <- 1
clin_df1[clin_df1$tumor_stage == "stage iia", "stage"] <- 2
clin_df1[clin_df1$tumor_stage == "stage ii", "stage"] <- 2
clin_df1[clin_df1$tumor_stage == "stage iib", "stage"] <- 2
clin_df1[clin_df1$tumor_stage == "stage iiia", "stage"] <- 3
clin_df1[clin_df1$tumor_stage == "stage iiib", "stage"] <- 3
clin_df1[clin_df1$tumor_stage == "stage iv", "stage"] <- 4



data <- clin_df1
data
library(tidyr)
data%>% drop_na()
clin_df1



clin_df1$stage
clin_df1$days_to_last_follow_up <- as.numeric(as.character(clin_df1$days_to_last_follow_up))
str(clin_df1)

kmbytreatment<-survfit(Surv(clin_df1$time, clin_df1$status)~ clin_df1$stage)
summary(kmbytreatment)
plot(kmbytreatment, col=9:15, main="Kaplan-Meier LUAD project for stages", xlab= "timeline (days)", ylab="Survival probability")
legend("topright",col=9:15, legend=c("stage 1","stage 2", "stage 3", "stage 4"), lty=1)

#Per saber si les 5 corbes són diferents?
survdiff(Surv(clin_df1$time, clin_df1$status)~ clin_df1$stage)



## relació entre estadi del tumor i els dos grups d'interès
taula <- table(clin_df1$stage,clin_df$chemo)
taula[,-2]


clin_df3 = clinical[,c("submitter_id",
                      "vital_status",
                      "days_to_last_follow_up",
                      "days_to_death",
                      "chemo", 
                      "tumor_stage",
                      "gender",
                      "age_at_index")]
idDead<-which(clin_df3$vital_status=="Dead")
idAlive<-which(clin_df3$vital_status=="Alive")

clin_df3$time<-rep(0,522)
for (i in idDead){
  clin_df3$time[i]<-as.numeric(as.character(clinical$days_to_death[i]))
}
for (i in idAlive){
  clin_df3$time[i]<-as.numeric(as.character(clinical$days_to_last_follow_up[i]))
}




clin_df3$status<-rep(0,522)
clin_df3[clin_df3$vital_status == "Dead", "status"] <- 1

#Li poso 1s a tot el que sigui Chemo
clin_df3$treat<-rep(0,522)
clin_df3[clin_df3$chemo == "yes", "treat"] <- 2
clin_df3[clin_df3$chemo == "no", "treat"] <- 1
clin_df3[clin_df3$chemo == "not reported", "treat"] <- NA


clin_df3$stage<-rep(0,522)
clin_df3$stage

clin_df3[clin_df3$tumor_stage == "not reported", "stage"] <- NA
clin_df3[clin_df3$tumor_stage == "stage ia", "stage"] <- 1
clin_df3[clin_df3$tumor_stage == "stage i", "stage"] <- 1
clin_df3[clin_df3$tumor_stage == "stage ib", "stage"] <- 1
clin_df3[clin_df3$tumor_stage == "stage iia", "stage"] <- 2
clin_df3[clin_df3$tumor_stage == "stage ii", "stage"] <- 2
clin_df3[clin_df3$tumor_stage == "stage iib", "stage"] <- 2
clin_df3[clin_df3$tumor_stage == "stage iiia", "stage"] <- 3
clin_df3[clin_df3$tumor_stage == "stage iiib", "stage"] <- 3
clin_df3[clin_df3$tumor_stage == "stage iv", "stage"] <- 4






data <- clin_df3
data
library(tidyr)
data%>% drop_na()
clin_df3




clin_df3$days_to_last_follow_up <- as.numeric(as.character(clin_df3$days_to_last_follow_up))
str(clin_df3)


cox.model<-with(clin_df3, coxph(Surv(clin_df3$time, clin_df3$status)~ clin_df3$treat+clin_df3$stage))
cox.model

###### Tant el tractament com l'stage són significatius




##Afegim l'edat i el sexe

clin_df3$sexe<-rep(0,522)
clin_df3[clin_df3$gender == "male", "sexe"] <- 1

clin_df3$age_at_index <- as.numeric(as.character(clin_df3$age_at_index))


cox.model<-with(clin_df3, coxph(Surv(clin_df3$time, clin_df3$status)~ clin_df3$treat+clin_df3$stage+clin_df3$gender+clin_df3$age_at_index))
cox.model

data <- clin_df3
data
library(tidyr)
data%>% drop_na()
clin_df1



##Hi ha molts més pacients en els no-tractats i tambe dels no-tractats en estadis avançats, de manera que:
#com ja hem vist amb el gràfic de superviviència segons l'estadi del tumor, si l'estadi estigues afectant la supervivència obtinguda
# en el primer gràfic, veuriem molta més quantitat de tumors en estadis avançats en els tractats amb quimioterapia
#pero com no és aixi
#es conclou que l'estadi del tumor no afecta la supervivència entre els dos grups d'interes
#pero esta clar que si que afecta a la supervivència dels pacients (sense discernir entre els grups de tractament)


taula2 <- table(clinical$age_at_index)
taula2

### RELACIÓ SIGNATURA SENESCÈNCIA I SUPERVIVÈNCIA

o<-as.numeric(rownames(ordered_table))
o
index3<-order(o)
table2<-ordered_table[index3,]
table2
index4<-which(table2$senescence_signature==1)
index4
table2$gene_name[index4]
dataSenes<-dataFilt[index4,]
dataSenes

#Després transforma els gens (els escalo a mitjana 0 i desviació 1 per tal que tots els gens siguin comparables)

dataSenes2<-apply(dataSenes,2,scale)

#Trasposo la matriu:
  
dataSenes2<-t(dataSenes2)

#Calculo la mitjana d'expressió dels 71 gens:

meanSenes<-apply(dataSenes2,1,mean)
plot(meanSenes)

#Ara hem d'agregar aquesta informació dels gens a les dades clíniques (el problema és que tenim 600 mostres i 522 pacients perquè hi ha mostres que són rèpliques o controls). Miro quins són els índexos que coincideixen:
  
index5<-NULL
for (i in (1:nrow(clinical))){
  posicio<-grep(clinical$submitter_id[i],colnames(dataSenes))[1]
  index5<-c(index5,posicio)
}

#Afegeixo la mitjana dels gens a la taula clinical_patient:
clinical$patient<-rep(0,522)
clinical_patient <- rep(0,522)
clinical_patient$meanSen <- rep(0,522)
clinical_patient$meanSen<-meanSenes[index5]

#Calculo la mitjana d'aquesta mitjana:

m<-mean(na.omit(clinical_patient$meanSen))
m

#Miro quins individus tenen una expressió de senescència per sobre la mitjana:

clinical_patient$sen<-ifelse(clinical_patient$meanSen>m,1,0)
clinical_patient$sen

#Ara hauries de fer l'anàlisi de supervivència amb aquesta variable "sen"
#Afegeixo la variable sen a clin_Df amb la informació que necessito per fer l'analisi de supervivenica
clin_df <- data.frame(clin_df,clinical_patient$sen)

#elimino els NAs
data <- clin_df
data
library(tidyr)
data%>% drop_na()
clin_df

#Passo a numerica la variable "temps"
str(data)
data$days_to_last_follow_up <- as.numeric(as.character(data$days_to_last_follow_up))
str(data)



library(survival)
library(KMsurv)


#Faig l'anàlisi
kmbytreatment<-survfit(Surv(data$time, data$status)~ data$clinical_patient.sen)
summary(kmbytreatment)
plot(kmbytreatment, col=2:3, main="Kaplan-Meier LUAD project", xlab= "timeline (days)", ylab="Survival probability")
legend("topright",col=2:3, legend=c("High expression","Low expression"), lty=1)

#saber si les dues corbes són diferents
survdiff(Surv(data$days_to_last_follow_up, data$status)~ data$clinical_patient.sen)


## el mateix però per TFG-B

o<-as.numeric(rownames(ordered_table))
o
index3<-order(o)
table2<-ordered_table[index3,]
table2
index4<-which(table2$tgf_signature==1)
index4
table2$gene_name[index4]
dataSenes<-dataFilt[index4,]
dataSenes

#Després transforma els gens (els escalo a mitjana 0 i desviació 1 per tal que tots els gens siguin comparables)

dataSenes2<-apply(dataSenes,2,scale)

#Trasposo la matriu:

dataSenes2<-t(dataSenes2)

#Calculo la mitjana d'expressió dels 71 gens:

meanSenes<-apply(dataSenes2,1,mean)
plot(meanSenes)

#Ara hem d'agregar aquesta informació dels gens a les dades clíniques (el problema és que tenim 600 mostres i 522 pacients perquè hi ha mostres que són rèpliques o controls). Miro quins són els índexos que coincideixen:

index5<-NULL
for (i in (1:nrow(clinical))){
  posicio<-grep(clinical$submitter_id[i],colnames(dataSenes))[1]
  index5<-c(index5,posicio)
}

#Afegeixo la mitjana dels gens a la taula clinical_patient:
clinical$patient<-rep(0,522)
clinical_patient <- rep(0,522)
clinical_patient$meanSen <- rep(0,522)
clinical_patient$meanSen<-meanSenes[index5]

#Calculo la mitjana d'aquesta mitjana:

m<-mean(na.omit(clinical_patient$meanSen))
m

#Miro quins individus tenen una expressió de senescència per sobre la mitjana:

clinical_patient$sen<-ifelse(clinical_patient$meanSen>m,1,0)
clinical_patient$sen

#Ara hauries de fer l'anàlisi de supervivència amb aquesta variable "sen"
#Afegeixo la variable sen a clin_Df amb la informació que necessito per fer l'analisi de supervivenica
clin_df <- data.frame(clin_df3,clinical_patient$sen)

#elimino els NAs
data <- clin_df
data
library(tidyr)
data%>% drop_na()
clin_df

#Passo a numerica la variable "temps"
str(data)
data$days_to_last_follow_up <- as.numeric(as.character(data$days_to_last_follow_up))
str(data)



library(survival)
library(KMsurv)


#Faig l'anàlisi
kmbytreatment<-survfit(Surv(data$time, data$status)~ data$clinical_patient.sen)
summary(kmbytreatment)
plot(kmbytreatment, col=2:3, main="Kaplan-Meier LUAD project", xlab= "timeline (days)", ylab="Survival probability")
legend("topright",col=2:3, legend=c("High expression","Low expression"), lty=1)

#saber si les dues corbes són diferents
survdiff(Surv(data$days_to_last_follow_up, data$status)~ data$clinical_patient.sen)



library(tidyr)
library(survival)
library(KMsurv)



## MIREM LA TAXA DE SUPERVIVÈNCIA FILTRANT PER ESTADIS: 

clin_df = clinical[,c("submitter_id",
                      "vital_status",
                      "days_to_last_follow_up",
                      "days_to_death",
                      "chemo")]
idDead<-which(clin_df$vital_status=="Dead")
idAlive<-which(clin_df$vital_status=="Alive")

clin_df$time<-rep(0,522)
for (i in idDead){
  clin_df$time[i]<-as.numeric(as.character(clinical$days_to_death[i]))
}
for (i in idAlive){
  clin_df$time[i]<-as.numeric(as.character(clinical$days_to_last_follow_up[i]))
}

#Li poso 1s a tot el que sigui Alive
clin_df$status<-rep(0,522)

#clin_df$status<-rep(0,395)

clin_df[clin_df$vital_status == "Dead", "status"] <- 1

#Li poso 1s a tot el que sigui Chemo
clin_df$treat<-rep(0,522)
clin_df[clin_df$chemo == "yes", "treat"] <- 1
clin_df[clin_df$chemo == "no", "treat"] <- 2
clin_df[clin_df$chemo == "not reported", "treat"] <- NA




##estadi 4
clin_df2<-clin_df[clin_df1$stage==4,]
data<-clin_df2
# library(tidyr)
data%>% drop_na()

kmbytreatment<-survfit(Surv(data$time, data$status)~ data$treat)
summary(kmbytreatment)
plot(kmbytreatment, col=2:3, main="Kaplan-Meier LUAD project STAGE 4", ylab=" Survival probability", xlab="timeline (days)")
legend("topright",col=2:3, legend=c("chemo","non-chemo"), lty=1)

#saber si les dues corbes s?n diferents
survdiff(Surv(data$time, data$status)~ data$treat)





### estadi 3:
clin_df3<-clin_df[clin_df1$stage==3,]
data<-clin_df3
# library(tidyr)
data%>% drop_na()

kmbytreatment<-survfit(Surv(data$time, data$status)~ data$treat)
summary(kmbytreatment)
plot(kmbytreatment, col=2:3, main="Kaplan-Meier LUAD project STAGE 3", ylab=" Survival probability", xlab="timeline (days)")
legend("topright",col=2:3, legend=c("chemo","non-chemo"), lty=1)

#saber si les dues corbes s?n diferents
survdiff(Surv(data$time, data$status)~ data$treat)



## estadi 2
clin_df4<-clin_df[clin_df1$stage==2,]
data<-clin_df4
# library(tidyr)
data%>% drop_na()

kmbytreatment<-survfit(Surv(data$time, data$status)~ data$treat)
summary(kmbytreatment)
plot(kmbytreatment, col=2:3, main="Kaplan-Meier LUAD project STAGE 2", ylab=" Survival probability", xlab="timeline (days)")
legend("topright",col=2:3, legend=c("chemo","non-chemo"), lty=1)

#saber si les dues corbes s?n diferents
survdiff(Surv(data$time, data$status)~ data$treat)


## estadi 1
clin_df5<-clin_df[clin_df1$stage==1,]
data<-clin_df5
# library(tidyr)
data%>% drop_na()

kmbytreatment<-survfit(Surv(data$time, data$status)~ data$treat)
summary(kmbytreatment)
plot(kmbytreatment, col=2:3, main="Kaplan-Meier LUAD project STAGE 1", ylab=" Survival probability", xlab="timeline (days)")
legend("topright",col=2:3, legend=c("chemo","non-chemo"), lty=1)

#saber si les dues corbes s?n diferents
survdiff(Surv(data$time, data$status)~ data$treat)






#guarda tots les objectes en forma d'imatge
#save.image("resultats.RData")
#save(rse, file="rse.RData")
#browseVignettes("TCGAbiolinks")
#packageVersion("TCGAbiolinks")
#Versió 2.14.1
save.image("resultatsLUAD_08062020.Rdata")
