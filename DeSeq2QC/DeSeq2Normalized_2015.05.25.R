

# FUNCTIONS THAT WORK ON DDS OBJECT
# PURPOSE TO EXPLORE NORMALIZED DATASET
# COULD BE USED DOWNSTREM THE FIRST SCRIPT TO PERFORM NORMALIZED EXPLORATION


library(RColorBrewer)
library(gplots)
library(DESeq2)


# Generation of objects useful for function testing


# creation of dds object 
SampleTable <- read.table('/Users/cheronicristina/Documents/Work/RNASeq_SBDS/7.DeSeq2/PoliTot/A.DifferentialExpression/Input/SampleTable.txt',header=TRUE)
Directory <- '/Users/cheronicristina/Documents/Work/RNASeq_SBDS/6.HTSeq/HTSeq_Count'
DDS <- DESeqDataSetFromHTSeqCount(sampleTable = SampleTable, directory = Directory, design= ~ Condition)
colData(DDS)$Condition<-factor(colData(DDS)$Condition, levels=c("CSIG_Poli", "PCIG_Poli", "CSIG_Tot","PCIG_Tot"))
DDS <- estimateSizeFactors(DDS)
DDS <- estimateDispersions(DDS)
DDS <- nbinomWaldTest(DDS)

NormalizedCount<-as.data.frame(counts(DDS, normalized=TRUE))
RLogMatrix <- assay(rlog(DDS,blind=T))
VstMatrix <- assay(varianceStabilizingTransformation(DDS,blind=T))

ResPoli <- results(DDS, alpha=0.05, contrast=c("Condition","PCIG_Poli","CSIG_Poli"), theta=c(0.65,0.65))


# ====================================
# ===  1.  EUCLIDEAN HEATMAP FUNCTION
# ====================================

# === The purpose is to generate an euclidean heatmapt starting from a matrix of read counts values (or other values) 
# === I would like the function to be autonomously applied both to save directly the graph in PDF and to produce a graph that can then be saved autonomously

# Arguments of the function:
# 1. normalizedcounts (matrix or dataframe): values of normalized and transformed read counts on which euclidean distance will be calculated. Values will not be further transformed
# 2. sidecolors (vector): color vector to identify in the heatmap the samples belonging to the same experimental group. Length must be the same of the number of samples in normalizedcount; order must be the same of samples in normalizedcounts
# 4. outdirectory (string): path to the output directory. Default is current directory 
# 5. pdf: specifies if a pdf file should be produced
# 6. title: title for the plot. Default is the name of the normalizedcount matrix



EuclideanHeatmap <- function(normalizedcounts, outdirectory=getwd(), pdf=TRUE, sidecolors=NULL, title=NULL) {
	# Generates an heatmap representing euclidean distance across samples
	# Relys on heatmap.2 function og gplots2 library
	
# 1. Definition of color palette for the heatmap
	hmcol<-colorRampPalette(brewer.pal(9,"GnBu"))(100)
	
# 2. distance function on the transposed matrix so to calculate sample-sample distance 
	DistMatrix <- as.matrix(dist(t(normalizedcounts)))
	
# 3. Setting of plot title
	if (is.null(title)) {
		title = deparse(substitute(normalizedcounts))
	}

# 4. For pdf, generation of directory and file	
	if (pdf == TRUE) {
		if (file.exists(paste(outdirectory, '/', 'Plots', sep='')) == FALSE) {
		dir.create(paste(outdirectory, '/', 'Plots', sep=''))
		}
		pdf(paste(outdirectory, '/', 'Plots/', title, '_EuclideanHeatmap.pdf', sep=''), width=12, height=12, pointsize=12, bg='white') 
	}
	
# 5. heatmap generation 
	if (is.null(sidecolors)) {
	heatmap.2(DistMatrix, Colv=TRUE, symm=TRUE, scale="none", margins=c(12,12), col = hmcol, main=paste(title, "\n Euclidean distance"),  trace="none", keysize=1.2, revC=T, cex.main=2.5, tracecol=F, cexRow=1.5, cexCol=1.5)
	} else {
		heatmap.2(DistMatrix, Colv=TRUE, symm=TRUE, scale="none", margins=c(12,12), col = hmcol, main=paste(title, "\n Euclidean distance"), ColSideColors=sidecolors,  trace="none", keysize=1.2, revC=T, cex.main=3, tracecol=F, cexRow=1.5, cexCol=1.5)
	}

if (pdf == TRUE) {
dev.off() }

}

	

# ===Test

# Test on the complete dataset 
EuclideanHeatmap(normalizedcounts=RLogMatrix, pdf=FALSE)

# Add color and a title 
colors <- c(rep("#1A9850",3),rep("#D73027",3),rep("#66BD63",3),rep("#F46D43",3))
EuclideanHeatmap(normalizedcounts=RLogMatrix, pdf=FALSE, sidecolors = colors, title = 'RNASeq SBDS Read Counts')

# Save as pdf
OutDirectory = '/Users/cheronicristina/Documents/Work/R/RPackages/FinalFunctions/DeSeq2QC'
EuclideanHeatmap(normalizedcounts=RLogMatrix, pdf=TRUE, sidecolors = colors, title = 'RNASeq SBDS Read Counts')


# Ok funziona - dovrei vedere per adattare le dimensioni del pdf in caso di un campione molto numeroso


# ====================================
# ===  2.  PEARSON HEATMAP FUNCTION
# ====================================

# === The purpose is to generate an euclidean heatmapt starting from a matrix of read counts values (or other values) 
# === I would like the function to be autonomously applied both to save directly the graph in PDF and to produce a graph that can then be saved autonomously

# Arguments of the function:
# 1. normalizedcounts (matrix or dataframe): values of normalized and transformed read counts on which euclidean distance will be calculated. Values will not be further transformed
# 2. sidecolors (vector): color vector to identify in the heatmap the samples belonging to the same experimental group. Length must be the same of the number of samples in normalizedcount; order must be the same of samples in normalizedcounts
# 4. outdirectory (string): path to the output directory. Default is current directory 
# 5. pdf: specifies if a pdf file should be produced
# 6. title: title for the plot. Default is the name of the normalizedcount matrix



PearsonHeatmap <- function(normalizedcounts, outdirectory=getwd(), pdf=TRUE, sidecolors=NULL, title=NULL) {
	# Generates an heatmap representing euclidean distance across samples
	# Relys on heatmap.2 function og gplots2 library
	
	
# 1. Correlation function on the transposed matrix so to calculate sample-sample distance 
	correlation <- cor(normalizedcounts)

# 2. Definition of color palette for the heatmap
	my.breaks = seq(min(correlation), 1, by = (1-min(correlation))/100)
	hmcol<-colorRampPalette(brewer.pal(9,"GnBu"))(100)
	
# 3. Setting of plot title
	if (is.null(title)) {
		title = deparse(substitute(normalizedcounts))
	}

# 4. For pdf, generation of directory and file	
	if (pdf == TRUE) {
		if (file.exists(paste(outdirectory, '/', 'Plots', sep='')) == FALSE) {
		dir.create(paste(outdirectory, '/', 'Plots', sep=''))
		}
		pdf(paste(outdirectory, '/', 'Plots/', title, '_PearsonHeatmap.pdf', sep=''), width=12, height=12, pointsize=12, bg='white') 
	}
	
# 5. heatmap generation 
	if (is.null(sidecolors)) {
	heatmap.2(correlation, Colv=TRUE, symm=TRUE, scale="none", margins=c(12,12), col = hmcol, main=paste(title, "\n Pearson Correlation"),  trace="none", keysize=1.2, revC=T, cex.main=2.5, tracecol=F, cexRow=1.5, cexCol=1.5, cellnote=round(correlation,2), notecol="black",notecex=0.9, breaks=my.breaks)
	} else {
	heatmap.2(correlation, Colv=TRUE, symm=TRUE, scale="none", margins=c(12,12), col = hmcol, main=paste(title, "\n Pearson Correlation"),  trace="none", keysize=1.2, revC=T, cex.main=2.5, tracecol=F, cexRow=1.5, cexCol=1.5, cellnote=round(correlation,2), notecol="black",notecex=0.9, breaks=my.breaks, ColSideColors=sidecolors)
	}

if (pdf == TRUE) {
dev.off() }

}


# ===Test

# Test on the complete dataset 
PearsonHeatmap(normalizedcounts=NormalizedCount, pdf=FALSE)
PearsonHeatmap(normalizedcounts=RLogMatrix, pdf=FALSE)

# Add color and a title 
colors <- c(rep("#1A9850",3),rep("#D73027",3),rep("#66BD63",3),rep("#F46D43",3))
PearsonHeatmap(normalizedcounts=RLogMatrix, pdf=FALSE, sidecolors = colors, title = 'RNASeq SBDS Read Counts')

# Save as pdf
OutDirectory = '/Users/cheronicristina/Documents/Work/R/RPackages/FinalFunctions/DeSeq2QC'
PearsonHeatmap(normalizedcounts=RLogMatrix, pdf=TRUE, sidecolors = colors, title = 'RNASeq SBDS Read Counts')



# ====================================
# ===  3.  HEATMAPNORMALIZED FUNCTION
# ====================================


NormalizedExploration <- function(dds, res, transformation='rlt', pdf=TRUE, PCA=TRUE, outdirectory=getwd(), savetable=TRUE, title=NULL, tested=TRUE) {

#1. Data organization
	NormalizedCount <- as.data.frame(counts(dds, normalized=TRUE))	
	
	if (transformation == 'rlt') {
		Trans <- rlog(dds,blind=T)  # devo lasciare questo step intermedio perchè viene utilizzato per la PCA
		Transformed <- assay(Trans)
	}

	if (transformation == 'vst') {
		Trans <- varianceStabilizingTransformation(dds,blind=T)
		Transformed <- assay(varianceStabilizingTransformation(dds,blind=T))	
	}
	
	if (tested == TRUE) {
	Tested <- is.na(res$padj) == FALSE
	Trans <- Trans[Tested,]
	Transformed <- Transformed[Tested,]	
	}

#3. Table saving
	
	if (is.null(title)) {
		title = deparse(substitute(dds))
	}
	
	if (savetable == TRUE) {
		if (file.exists(paste(outdirectory, '/', 'Tables', sep='')) == FALSE) {
		dir.create(paste(outdirectory, '/', 'Tables', sep=''))
		}
		write.table(NormalizedCount, file=paste(outdirectory, '/', 'Tables/', title, '_NormalizedReadCounts.txt', sep=''), sep='\t') 
		write.table(Transformed, file=paste(outdirectory, '/', 'Tables/', title, '_TransformedReadCounts.txt', sep=''), sep='\t') 
	}

#2. Heatmap generation
	EuclideanHeatmap(normalizedcounts=Transformed, pdf=pdf, outdirectory=outdirectory, title=title)
	
	PearsonHeatmap(normalizedcounts=Transformed, pdf=pdf, outdirectory=outdirectory, title=title)

# PCA generation
if (PCA == TRUE) {
	if (file.exists(paste(outdirectory, '/', 'Plots', sep='')) == FALSE) {
	dir.create(paste(outdirectory, '/', 'Plots', sep=''))
		}
	pdf(paste(outdirectory, '/', 'Plots/', title, '_PCA.pdf', sep=''), width=8, height=7, pointsize=12, bg='white')
print(plotPCA(Trans, names(colData(dds))[1]))
dev.off()
}
	
}



# =========== TEST

OutDirectory = '/Users/cheronicristina/Documents/Work/R/RPackages/FinalFunctions/DeSeq2QC'
NormalizedExploration(dds=DDS, outdirectory=OutDirectory, tested=FALSE, pdf=TRUE, PCA=FALSE)

NormalizedExploration(dds=DDS, res=ResPoli, transformation='rlt', outdirectory=OutDirectory, savetable=TRUE, title='RNASeq SBDS RLT', tested=TRUE, pdf=TRUE, PCA=TRUE)

NormalizedExploration(dds=DDS, res=ResPoli, transformation='vst', outdirectory=OutDirectory, savetable=TRUE, title='RNASeq SBDS VST', tested=TRUE, pdf=TRUE, PCA=TRUE)

# Everything seems to work properly 




