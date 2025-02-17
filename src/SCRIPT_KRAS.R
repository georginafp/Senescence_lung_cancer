#### MICROARRAY KRAS ANALYSIS ###


load("C:/Users/Georgina/Desktop/KRAS/GSE65258/resultatsKRAS.Rdata")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("GEOquery")
library(GEOquery) 
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("oligo")
library(oligo)
#Assignem el nostre directori de treball: 
workingDir<-"C:/Users/Georgina/Desktop/KRAS"
setwd(workingDir)
#Comprovem que estem al directori correcte:
getwd()
#Mirem quins arxius/fitxers hi ha en aquest directori per acabar de confirmar que estem al correcte:
list.files()
#Ens baixem les dades de GEO per a convertir-les a una matriu i poder-les processar:
#Entrem el codi de la s�rie de dades dins de GEO:
gseid<-"GSE65258"
#Ens baixem les dades NORMALS
gse <- getGEO(gseid,GSEMatrix=TRUE) 
#MIREM EL PRIMER ITEM
gse[[1]]
#Obtenim el phenoData
pD <- pData(gse[[1]]) 
pD
#Ens baixem les DADES CRUES, per tal de poder fer un control de qualitat d'aquestes
GEO<-getGEOSuppFiles(gseid,makeDirectory=TRUE) 
setwd(paste(workingDir,gseid,sep="/"))
#Utilitzem la funci� untar per tal de descomprimir el fitxer ".tar" directament dins del directori 
untar("GSE65258_RAW.tar", exdir = getwd())
######Anem direcrament a la carpeta del nostre directori per a comprovar que s'ha descomprimit correctament:
######i que per tant ens trobem amb uns fitxers "GZ" de format ".cel"
#Comprovem amb la funci� list.files() que tenim el fitxers descomprimit i que tenim un fitxer per cada mostra: 
list.files()
#Carreguem/llistem les dades pr�viament carregades seleccionant-les a partir de :
#"pattern" per tal de que nom�s em carregui aquells fitxerts que siguin "CEL.gz" per tal de:
#Assegurar-me que nom�s seleccionem els fitxers que volem
celFiles = list.files( pattern = "CEL.gz")
#Proseguim doncs, carregant les dades a partir de la funcio "read.celfiles()" del paquet Oligo:
GEOFS <- read.celfiles(celFiles)


#1. Comen�em a EXPLORAR les dades:
GEOFS 
#Veiem que estem davant d'un GeneFeatureSet
class(GEOFS) 
slotNames(GEOFS)
annotation(GEOFS)
phenoData(GEOFS)
pData(GEOFS)
#dimensio
dim(exprs(GEOFS))

#Analitzem els primers valors/elements de la matriu:
head(exprs(GEOFS))
#Veiem que els noms s�n molt llargs, i aix� ho podem solucionar amb la funcio  "gsub()":
#El que estem fent amb aquesta funci� �s substituir una part del nom de les variables per " " per tal d'aixi tenir noms m�s curts
sampleNames(GEOFS)
sampleNames(GEOFS)<-gsub(".CEL.gz","",sampleNames(GEOFS))
#Comprovem que la escursaci� dels noms s'ha fet correctament:
sampleNames(GEOFS)
#Mirem ara una altra vegada els primers elements de la matriu
#Veiem que ja cont� els noms escursats
head(exprs(GEOFS))
#A la matriu d'intensitats hem de tenir les sondes com a files i les mostres com a columnes


#2. Un cop tenim les dades i veiem com �s l'objecte, fem el CONTROL DE QUALITAT:
#Per a fer el control de qualitat seguirem uns passos: 
#2.1. Mirar les imatges 
#2.2. Gr�fic de densitats (Histograma)
#2.3. Gr�fic de tipus Boxplot
#2.4. Gr�fic de tipus MAplot

#Primer pas del control de qualitat:
image(GEOFS)
#No hi ha taques ni anomalies i per tant, podem donar les dades com a bones i continuar amb l'an�lisis.
#petita taca en alguna imatge que creiem que no ens interferir� molt en els resultats 

#Ara hem de generar un histograma i un boxplot de les mostres de l'experiment
#L'histograma et permet veure com es distribueixen les intensitats en quant a valors:
#El boxplot et representa la distribuci� d'intensitats en base a la mediana (caixes de 25, mediana i 75)
#El MAplot on mira la intensitat de les diferents mostres respecte de la resta
#El que esperariem veure s�n uns nivells d'intensitats iguals

#Assignem una paleta de colors (rainbow) 
colos<-rainbow(6)
#Esperem que hi hagi homogeneitat en totes les distribucions i veiem que els colors no segueixen un patro igual
hist(GEOFS,target="core",col=colos)
#Fem el gr�fic de tipus boxplot:
boxplot(GEOFS,target="core",col=colos,las=3,cex.axis=0.5)
#No segueixen una distribuci� semblant 
#Fem el gr�fic de tipus MAplot :
MAplot(GEOFS)
#Veiem un n�vol de punts que son les diferents sondes en el microarray en cada mesura
#Sabem que com: 
#M�s llunyana a linia vermella del nostre 0, m�s diferencia hi haura d'aquell microarray en comparacio als altres
#Per tant podem veure moltes difer�ncies entre ells.
#Fem el mateix pero de dos en dos:
MAplot(GEOFS,pairs=TRUE)



#3. Ara que ja hem decidit que el control de qualitat �s correcte, procedim a la NORMALITZACI� DE LES DADES:
#Aquest apartat el que far� ser� eliminar petits biaixos t�cnics en quant a distribuci� de les mostres
#I tamb� eliminar� problemes t�cnics en el processament de les mostres al laboratori
#Per tant el que estarem fent amb aquest pas �s corregir el soroll de background 

#3.1. Per QUANTILS:
# --> Equipara cadascun dels quantils, fent com com una mena de mitjana dels gr�fics que hi ha en aquell quantil. 
# --> �s una mesura de centralitat dels quantils, per tant arrosega les corbes per a que siguin equiparables. 
# --> Arrossega les corbes de densitat de les nostres mostres a una corba de mitjana que ell mateix ha generat. 

#3.2. SUMARITZACI�:
# --> Tenim moltes dades de les diferents sondes que s'han unit als exons i el que volem �s sumaritzar-ho tot en un
# --> Volem obtenir un �nic valor pel meu transcrit. 

#3.3. TRANSFORMACI� A BASE LOGAR�TMICA:
# --> Aplicar-ho tot a base log ja que:
# per aplicar un t-test necessitem que les dades estiguin distribuides de manera normal per m�todes param�trics 

library(oligo)

#Per a fer-ho utilitzem la funcio RMA del paquet de oligo 
GEOFS.rma<-rma(GEOFS,background=TRUE, normalize=TRUE, subset= NULL)
#Tornem a mirar la dimensi� ja que al corregir el background i fer el pas de sumaritzaci� l'hauriem d'haver redu�t:
dim(exprs(GEOFS.rma))
#Efectivament hem reduit la dimensi� de la matriu
#Observo les meves dades normalitzades 
GEOFS.rma
#Veiem que estem treballant amb un ExpressionSet (ja que les dades ja estan normalitzades)
class(GEOFS.rma)


#4. Ara passem al pas de AGREGACI� DE LES DADES
#Aquest pas ens permet detectar si hi ha outliers. 
#Estarem generant clusters jer�rquics per a veure com s'agrupen els gr�fics segons un dandograma
#Per� per a fer aquest pas necessitem: 
#4.1. Les dist�ncies entre mostres --> definir dist�ncia entre mostres o dist�ncies basades en la correlaci�
#4.2. Un m�tode de lincatge --> per a decidir si una mostra est� "a prop" o "lluny".
#Guardo la matriu de intensitats de les dades normalitzades a la variable "x"
x<-exprs(GEOFS.rma)
#Per a fer el dandrograma:
#Utilitzo la funcio hclust que est� dins del paquet stats
#Creo la dist�ncia basada en la correlaci� 
#Utilitzo el m�tode de: "ward.d2"
#genero l'objecte clus.cor.ward i en faig un plot
clust.cor.ward<- hclust(as.dist(1-cor(x)),method="ward.D2")
plot(clust.cor.ward, main="hierarchical clustering", hang=-1,cex=0.6)
#Probem el mateix per� amb diferent m�todes de lincatge:
clust.cor.complete<- hclust(as.dist(1-cor(x)),method="complete")
plot(clust.cor.complete, main="hierarchical clustering", hang=-1,cex=0.6)
#Altre m�tode de lincatge:
clust.cor.complete<- hclust(as.dist(1-cor(x)),method="average")
plot(clust.cor.complete, main="hierarchical clustering", hang=-1,cex=0.6)

#Farem el mateix an�lisi per� amb una reducci� de la dimensi� de l'espai:
#La t�cnica de la reducci� de la dimensi� el que fara �s:
#Dels milers de variables que tenim representades les reduir� a 3 per tal de poder-les representar
#Fent una combinaci� de gens i explorar al m�xim la variabilitat
#Per a continuar veient com s'agrupen les mostres entre si 
#Crido la llibreria necess�ria:
library(scatterplot3d)
#Per aplicar el m�tode aplico la funci� "prcom" (de les meves dades) i ho he de transposar.
#Veiem que tenim 6 components principals:
summary(pca.filt <- prcomp(t(x), cor = TRUE )) 
#Ho fem per saber quanta variabilitat teniem
#Un cop tenim la informaci� grafiquem:
pca3d<-scatterplot3d(x=pca.filt$x[,1],y=pca.filt$x[,2],z=pca.filt$x[,3],  xlab='PC1', ylab='PC2', zlab='PC3', main='PCA', pch=16,col.grid="lightblue")
#Afegim el nom de les mostres al gr�fic:
text(pca3d$xyz.convert(pca.filt$x+0.5), labels=rownames(pca.filt$x), cex=0.6)


#5. Ara prosseguirem amb l'ANALISI D'EXPRESSIO DIFERENCIAL
#Quan tenim les variables continues per comparar-les fem un t-test que ens permet comparar les dues condicions d'una mateixa variable o dues. 
#Aquest an�lisi afegeix un terme de moderaci� al denominador per corregir la limitaci� en el n�mero de mostres 
#El m�tode que utilizem �s "limma" - linear models for microarray analysis
#A m�s tenim la informaci� relacionada amb el fold-change (mitjana d'un grup respecte de la mitjana de l'altre)
#Farem el adjusted p-value

#Per tant, per cada gen tenim 3 mesures: 
#1. foldchange diferencies entre les mitjanes dels grups
#2. p-valor resultat del t-test moderat tenint en compte mitjanes, desviaci� est�ndard i mostres
#3. pvalor adjustat que es el p-valor ajustat tenint en compte el FDR

#On nom�s escollirem els que tenen un adj p-value < 0.05
#Cridem la llibreria necess�ria per l'an�lisis:
library(limma)
#Necessitem un vector de condicions i el guararem dins l'objecte "cond"
cond<-as.factor(c(rep("normal",8),rep("hyperplastic",7), rep("normal",8), (rep("FullBlown",7))))
#Necesitem definir el diseny de la matriu:
#Creem per cada categoria una columna, i a cada element de la columna li assignem la condici� de la categoria que est� representant
#Com en codis binaris m'assignar� 1 i 0 depenent de la condici� en la que estiguem.
design<-model.matrix(~0+cond)
colnames(design)<-gsub("cond","",colnames(design))
rownames(design)<-sampleNames(GEOFS)
design

#La funci� "lmfit" �s l'adjust del model. Per tant estem aplicant el model:
fit<-lmFit(GEOFS.rma,design)
#Apliquem els contrastos: 
#Estem comparant hyperplastic i Full-Blown  en el nivells de la matriu del disseny:
contrast.matrix<-makeContrasts(hyperplastic-FullBlown,levels=design)
contrast.matrix
#Aplico la matriu de contrastos:
fit2<-contrasts.fit(fit,contrast.matrix)
#Aplico l'adjust baiesi� i tinc el resultat fite:
fite<-eBayes(fit2)
fite


#Per extreure els resultats
#Posarem coeficient 1 perqu� tenim una comparaci�
#Number infinit perqu� vull extreure les comparacions de tots els gens
#i ladjust �s l'adjust per comparacions multiples:
#Adjust per FDR:
top.table<-topTable(fite,coef=1,number=Inf,adjust="BH")
results<-decideTests(fite)
table(results) 


#Anem a veure com es distribueixen els p-valors:
#Els "breaks" s�n el n�mero de talls que vull en el meu histograma
#El "Main" �s el t�tol que li adjudico al gr�fic
hist(top.table$P.Value,breaks=100,main="results P-value prostate adenocarcinoma")
#Podem considerar una bona distribuci� dels p-valors 

#Com que nom�s volem agafar els valors inferiors a 0.05 de la topable:
#Selecciono la columna que siguin m�s petits que 0.05
results.p0.05<-top.table[top.table$P.Value<0.05,]
dim(results.p0.05)
#A m�s vull que tinguin un logfoldchange superiors a 1
#guardo els valors com a result.p005:
results.p0.05.logFC1<-top.table[top.table$P.Value<0.05 & abs(top.table$logFC)>1,]
dim(results.p0.05.logFC1)


#6. Ara prosseguim al pas de ANOTACI� perqu� volem passar del ID del gen al nom del gen
#Cridem les llibreries necess�ries per aquest proc�s:
library(annotate)
??annotate

##### Hem d'utilitzar una db de mouse per a poder procedir amb l'anotaci� 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("mouse4302.db")
library(mouse4302.db)
mouse4302()
#Aquest pas, si estigu�ssim tractant amb ua altre esp�cie com "human" ho haur�em de canviar i utilitzar una de homo sapiens



#Construim un fitxer "html" que ens permetr� guardar o sumaritzar la informaci� dels gens diferencialment expressats
#El dat ens fa una matriu d'expressi� de les dades normalitazdes 
#Extracci� dels valors crus de les meves 30 mostres de la matriu normalitzada
dat <- exprs(GEOFS.rma)[rownames(results.p0.05.logFC1),]
#Guardem cada variable seva en un vector 
logFC <- results.p0.05.logFC1$logFC
pval <- results.p0.05.logFC1$P.Value
adj.pval<-results.p0.05.logFC1$adj.P.Val

#Accedeixo a les informacions que necessito:
#El s�mbol
symbol<-mget(rownames(results.p0.05.logFC1), env=mouse4302SYMBOL)
#El nom
name<-mget(rownames(results.p0.05.logFC1), env=mouse4302GENENAME)

#Generaci� del fitxer "html" 
affyids<-rownames(results.p0.05.logFC1)
genelist <- list(affyids)
filename <- "ResultatsKRAS1.html"
title <- "Gens diferencialment expressats quan el p<0.05 i  el logFC>1"
othernames <- list(symbol,name,round(logFC, 1), round(pval, 4), round(adj.pval, 4)) 
head <- c("Probe ID", "Symbol", "Gene Name", "logFC", "p-value","adj.p-value",sampleNames(GEOFS.rma))
repository <- list("affy")
htmlpage(genelist, filename, title, othernames, head, repository = repository)
#PER A SABER ON L'HA GUARDAT
getwd()


#Generem un heat map. Aquest �s una agrupaci� de dos cl�sters jer�rquics posats en files i columnes de la meva matriu de dades normalitzades
#El verd significa que �s negatiu i el vermell significa que �s positiu
#Cada gen es representa en una fila. La condici� que tenim representada per exemple, la verda significar� infra-expressats 
library(gplots)
#El fem per files:
data.clus<-exprs(GEOFS.rma[rownames(results.p0.05.logFC1),])
#Els rownames s�n els simbols obtinguts en el paquet d'anotacions 
rownames(data.clus)<-symbol
colnames(data.clus)


#Faig els clusters:

        #### AQUI �S ON TINC ERROR ### 
memory.size(max = FALSE)
memory.limit(size = 20000)

# Quan vull fer el cl�ster per a les files em diu que 'R' no pot ubicar un vector de mida 1.5Gb. 
clust.rows <- hclust(as.dist(1-cor(t(data.clus))),method="ward.D2")
clust.cols <- hclust(as.dist(1-cor(data.clus)),method="ward.D2")  


#Estableixo una paleta de colors
 heatcol<-colorRampPalette(c("green", "Black","red"), space = "rgb")
 pdf("heatmap_KRAS.pdf") 
# 
# #El fem per columnes
# #Crido la funci� utilizant la matriu "data.clus" guardada (extreta de les dades normalitzades) 
# #De la paleta de color n'extrec 256 colors (del vermell al verd passaant pel negre)
heatm<-heatmap.2(as.matrix(data.clus), col = heatcol(256),
                  dendrogram="column", Colv=as.dendrogram(clust.cols),    #dendrogrames 
                  Rowv=as.dendrogram(clust.rows),     #escalem per files 
                  scale="row",cexRow=0.3, cexCol=0.5,    
                  main="",key=TRUE,keysize=1,density.info="none",trace="none")
# #Aix� ho fem per guardar una imatge que dona R. 
# #I nom�s he donar la instrucci� del nom ja que ja li estic dient que la extensi� ha de ser ".pdf"
# #Per a guardar-lo al directori que he establert incialemnt:
 dev.off() 


#Ho he provat d'una altra manera per�:
library("gplots")
#No me'l genera
heatmap.2(as.matrix(data.clus), col=topo.colors(75), scale="none", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)
heatmap(as.matrix(data.clus))


x <- scale(as.matrix(data.clus))
heatmap.2(x,
          scale     = "none",
          trace     = "none",
          col       = heatcol(256),
          distfun   = function(x) as.dist(1-cor(t(x))), 
          hclustfun = function(x) hclust(x, method="ave"),
          margins=c(4,7)
)



# VOLCANOPLOT
with(top.table, plot(logFC, -log10(P.Value), pch=20, main="Volcano plot DEGs KRAS", xlim=c(-10,10)))

#si el p-valor adj<0.05
with(subset(top.table, adj.P.Val<.05 ), points(logFC, -log10(P.Value), pch=20, col="red"))
#si logFC>1
with(subset(top.table, abs(logFC)>1), points(logFC, -log10(P.Value), pch=20, col="orange"))
#si p-valor adj<0.05 i logFC>1
with(subset(top.table, adj.P.Val<.05 & abs(logFC)>1), points(logFC, -log10(P.Value), pch=20, col="green"))
legend(-10, -3, legend=c("Adjusted p-value < 0.05", "logFC > 1", "Adjusted p-value < 0.05 and logFC > 1"),
      col=c("red", "orange", "green"), lty=1:2, cex=0.8)





####                      GSEA  --> SENESCENCE SIGNATURE
#results.p0.05.logFC1 = genelist
gene_name<-rownames(results.p0.05.logFC1)
symbol<-mget(rownames(results.p0.05.logFC1), env=mouse4302SYMBOL)
length(symbol)
results1 <- results.p0.05.logFC1

results1$gene_name <- rep (0,20122)
results1$gene_name <- symbol
colnames(results1)


results1$logFC <- NULL
results1$AveExpr <- NULL
results1$t <- NULL
results1$adj.P.Val <- NULL
results1$B <- NULL
tablex <- as.matrix (results1)
table2 <- as.data.frame(tablex)
gene_name <- table2$gene_name
gene_name
pvalue<- table2$P.Value


#Afegeixo una columna m�s que es per la signatura de senescencia per� nomes la omplo amb 0s
table2$senescence_signature<-rep(0,20122)


#Vector amb els noms dels gens de la signatura de senescencia
senescence_vector=c("Acvr1", "Acvr1b", "Acvr1c", "Acvr2a", "Acvr2b", "Acvrl1", "Amhr2", "Atf2", "Bambi", "Bmp1", "Bmp10", "Bmp15", "Bmp2", "Bmp3", "Bmp4","Bmp5", "Bmp6", "Bmp7", "Bmp8a", "Bmp8b", "Bmpr1a", "Bmpr1b", "Bmpr2", "Cited1", "Cited2", "Crebbp", "Dcp1b", "Ep300", "Fkbp1a", "Fosl1", "Foxh1", "Gdf1", "Gdf10", "Gdf11", "Gdf15", "Gdf2", "Gdf3", "Gdf5", "Gdf6", "Gdf7", "Gdf9", "Gdnf", "Hras", "Inhba", "Inhbb", "Inhbc", "Inhbe", "Jun", "Junb", "Jund", "Lefty1", "Lefty2", "Map3k7", "Map3k7cl", "Mapk1", "Mapk10", "Mapk11", "Mapk12", "Mapk13", "Mapk14", "Mapk3", "Mapk8", "Mapk9", "Mstn", "Nodal", "Nras", "Rras", "Ski", "Skil", "Smad1", "Smad2", "Smad3", "Smad4", "Smad5", "Smad6", "Smad7", "Smad9", "Smurf1", "Smurf2", "Snip1", "Tab1", "Tgfb2", "Tgfb3", "Tgfbr1", "Tgfbr2", "Tll1", "Tll2", "Zfyv9")
length(senescence_vector)
#vector amb els gens diferencialment expressats meus 
#Que son els que vull veure si coincideixen (o estan dins -la majoria) amb els del vector de la signatura de senesc�ncia 

#Li poso 1s a tots aquells gens diferencialment expressats que coincideixen amb els del vector de la signatura de senescencia 
for (s in senescence_vector){
  gene_exists <- grep(s, gene_name)
  if (length(gene_exists) > 0){
    table2[table2$gene_name == s, "senescence_signature"] <- 1
  }
  
  
}


#Surten m�s de 77 perqu� n'hi ha alguns repetits pero no n'esborrar� cap perqu� no crec que influeixi molt en els resultats



#Hypergeometric test from an ordered table
pvalgsea<-NULL
x<-seq(1,1000,by=10)
m<-sum(table2$senescence_signature) # m = nombre total de senesc
m
n<-nrow(table2)-m # n= nombre total no de senesc
n
for (i in x){
  (q<-sum(table2$senescence_signature[(1:i)]))  # q = nombre de gens senesc entre els top 100
  pval<-phyper(q-1, m, n, i, lower.tail = FALSE)
  pvalgsea<-c(pvalgsea,pval)
}

y<-pvalgsea  # p-values sense log
length(y)
length(x)
plot(x, y, type="o", col="orange", main= "GSEA for senescence signature in KRAS project", xlab="ordered genes", ylab="p-value")
#no hi ha enriquiment perque ho hem graficat segons el p-valor
#Com m�s petit m�s significatiu --> els "top genes" (els 100 primers) no tenen els p-valors m�s petits
#per aix� no hi ha enriquiment




library(ggplot2)
library(gridExtra)
# 
#



## extraiem els 30 primers resultats dels DEGs per la mem�ria
df <- results.p0.05.logFC1
df[1:30,]
png("DEGS_KRAS.png", height=2000, width=2000)
p<-tableGrob(df)
grid.arrange(p)
dev.off()




length(results)
length(results.p0.05.logFC1$P.Value)






####                      GSEA  --> TGF SIGNATURE
#Encara que posi senescence es de TFG!!! 

#results.p0.05.logFC1 = genelist
gene_name<-rownames(results.p0.05.logFC1)
symbol<-mget(rownames(results.p0.05.logFC1), env=mouse4302SYMBOL)
length(symbol)
results1 <- results.p0.05.logFC1

results1$gene_name <- rep (0,20122)
results1$gene_name <- symbol
colnames(results1)


results1$logFC <- NULL
results1$AveExpr <- NULL
results1$t <- NULL
results1$adj.P.Val <- NULL
results1$B <- NULL
tablex <- as.matrix (results1)
table2 <- as.data.frame(tablex)
gene_name <- table2$gene_name
gene_name
pvalue<- table2$P.Value



#Afegeixo una columna m�s que es per la signatura de senescencia per� nomes la omplo amb 0s
table2$senescence_signature<-rep(0,20122)


#Vector amb els noms dels gens de la signatura de senescencia
senescence_vector=c("Acvr1", "Acvr1b", "Acvr1c", "Acvr2a", "Acvr2b", "Acvrl1", "Amhr2", "Atf2", "Bambi", "Bmp1", "Bmp10", "Bmp15", "Bmp2", "Bmp3", "Bmp4","Bmp5", "Bmp6", "Bmp7", "Bmp8a", "Bmp8b", "Bmpr1a", "Bmpr1b", "Bmpr2", "Cited1", "Cited2", "Crebbp", "Dcp1b", "Ep300", "Fkbp1a", "Fosl1", "Foxh1", "Gdf1", "Gdf10", "Gdf11", "Gdf15", "Gdf2", "Gdf3", "Gdf5", "Gdf6", "Gdf7", "Gdf9", "Gdnf", "Hras", "Inhba", "Inhbb", "Inhbc", "Inhbe", "Jun", "Junb", "Jund", "Lefty1", "Lefty2", "Map3k7", "Map3k7cl", "Mapk1", "Mapk10", "Mapk11", "Mapk12", "Mapk13", "Mapk14", "Mapk3", "Mapk8", "Mapk9", "Mstn", "Nodal", "Nras", "Rras", "Ski", "Skil", "Smad1", "Smad2", "Smad3", "Smad4", "Smad5", "Smad6", "Smad7", "Smad9", "Smurf1", "Smurf2", "Snip1", "Tab1", "Tgfb2", "Tgfb3", "Tgfbr1", "Tgfbr2", "Tll1", "Tll2", "Zfyve9", "Tsr1")
length(senescence_vector)
#vector amb els gens diferencialment expressats meus 
#Que son els que vull veure si coincideixen (o estan dins -la majoria) amb els del vector de la signatura de senesc�ncia 

#Li poso 1s a tots aquells gens diferencialment expressats que coincideixen amb els del vector de la signatura de senescencia 
for (s in senescence_vector){
  gene_exists <- grep(s, gene_name)
  if (length(gene_exists) > 0){
    table2[table2$gene_name == s, "senescence_signature"] <- 1
  }
  
  
}


#acvrl, bmp10, bmp2, bmp3, bmp6, bmp7, bmp8a, bmpr1a, bmpr2, cited1, ep300, gdf1, gdf10, gdf2, gdf3, gdf5, gdf6, gdf7, gdf9, gdnf, inhba, inhbc, jun, Mapk38, mapk3CL, mapk10, mapk11, mstn, NRas, Smad9, Tll1, Tll2,



#Hypergeometric test from an ordered table
pvalgsea<-NULL
x<-seq(1,1000,by=10)
m<-sum(table2$senescence_signature) # m = nombre total de senesc
m
n<-nrow(table2)-m # n= nombre total no de senesc
n
for (i in x){
  (q<-sum(table2$senescence_signature[(1:i)]))  # q = nombre de gens senesc entre els top 100
  pval<-phyper(q-1, m, n, i, lower.tail = FALSE)
  pvalgsea<-c(pvalgsea,pval)
}

y<-pvalgsea  # p-values sense log
length(y)
length(x)
plot(x, y, type="o", col="orange", main= "GSEA for TGF signature in KRAS project", xlab="ordered genes", ylab="p-value")
#no hi ha enriquiment perque ho hem graficat segons el p-valor
#Com m�s petit m�s significatiu --> els "top genes" (els 100 primers) no tenen els p-valors m�s petits
#per aix� no hi ha enriquiment




### n'hi ha 73 de 88 





save.image("resultatsKRAS.Rdata")
