# Update: 14 May 2015


# ====================================
# ===  1.  STARLOG PLOT FUNCTION
# ====================================

# === The purpose is to generate a function that can be used fro plotting reads relative to different STAR diagnostic parameters
# === I would like the function to be autonomously applied both to save directly the graph in PDF and to produce a graph that can then be saved autonomously


# Arguments of the function:
# 1. InputReads: number of reads to be plottes
# 2. SampleVector (vector): names of the samples - the order must be the same of InputReads
# 3. GroupVector (vector?): Experimental Group for each sample. If not specified, each sample is plotted independently- The order must be the same of InputReads
# 4. OutDirectory (string): path to the output directory. Default is current directory 
# 5. Pdf: specifies if a pdf file should be produced
# 6. Main: title for the plot. Default is the name of the vector for input reads

ReadStripPlot <- function(InputReads, SampleVector, GroupVector=NA, OutDirectory=getwd(), Pdf=TRUE, Main = NULL) {

	# I perform log10 transformation that will be used in the representation
	LogReads = log10(InputReads)
	# I calculate the z-score for each value, It will be used for the conditional color of the data points
	Zscore = abs((LogReads-mean(LogReads))/sd(LogReads))

	if (is.na(GroupVector[1])==TRUE) {
		GroupVector = SampleVector
	}

# Setting of plot title
	if (is.null(Main)) {
		Main = deparse(substitute(InputReads))
	}

# For pdf, generation of directory and file	
	if (Pdf == TRUE) {
		if (file.exists(paste(OutDirectory, '/', 'Plots', sep='')) == FALSE) {
		dir.create(paste(OutDirectory, '/', 'Plots', sep=''))
		}
		pdf(paste(OutDirectory, '/', 'Plots/', Main, '.pdf', sep=''), width=length(unique(GroupVector)), height=length(unique(GroupVector))/3*2.2, pointsize=(8 + length(unique(GroupVector))/100*80), bg='white') 
		par(mar=c(7,5,3,3))
	}	

# # Generation of the plot	
	stripchart(LogReads~GroupVector, method='jitter', jitter=0.1, vertical=T,pch=21, cex=1.5, bg='cornflowerblue', main=paste(Main, '\n','Number of Reads -log10-'), col.main='darkorchid', ylab='Number of Paired Reads', las=2, cex.lab=0.8, cex.axis=0.75, xlim = c(0,length(unique(GroupVector))+1))
	text(LogReads~factor(GroupVector), label=SampleVector, cex=ifelse(Zscore > 2,0.7,0.5), col=ifelse(Zscore > 2,"darkorchid","grey"), pos=4, font=2) 
	abline(h= seq(min(LogReads), max(LogReads), by=0.1), col='grey', lty=2) 

# For pdf, the connection is closed	
	if (Pdf == TRUE) {
dev.off() 
}

}


# ====== TESTS

# Load of data
ReadDataFrame <- read.table('/Users/cheronicristina/Documents/Work/R/RPackages/FinalFunctions/AlignmentQC/RElaborated/Tables/StarNumbers.txt', sep='\t', header=T) 


# Test of plot function without production of plot
# The setting of OutDirectory in this case is not important since no plot is ploduced
Prova <- ReadStripPlot(ReadDataFrame$Number.of.input.reads, SampleVector=ReadDataFrame$Parameter, OutDirectory=getwd(), Pdf=FALSE)

# Test of plot function with saving of the plot
OutputDir <- '/Users/cheronicristina/Documents/Work/R/RPackages/FinalFunctions/AlignmentQC'
ReadStripPlot(ReadDataFrame$Number.of.input.reads, SampleVector=ReadDataFrame$Parameter,  GroupVector=NA, OutDirectory=OutputDir, Pdf=TRUE)
# Specifying title:
ReadStripPlot(ReadDataFrame$Number.of.input.reads, SampleVector=ReadDataFrame$Parameter, Main='STAR_Input_Reads', GroupVector=NA, OutDirectory=OutputDir, Pdf=TRUE)


# Test of plot function giving a group vector
Groups <- c(rep('A',8), rep('B',8), rep('C',8), rep('D',8), rep('E',8), rep('F',8))
ReadStripPlot(ReadDataFrame$Number.of.input.reads, SampleVector=ReadDataFrame$Parameter, Main='STAR_Input_Reads', GroupVector=Groups, OutDirectory=OutputDir, Pdf=TRUE)



# ====================================
# ===  2.  STARLOG FUNCTION
# ====================================

# === STAR OUTPUT

# Purpose of the function is to retrieve the parameters deriving from STAR log final out file

# Arguments
# 1. InDirectory (string): path to the directory in which are stored the STAR output files that contain .final.out Log files. Default is the current directory
# 2. OutDirectory (string): path to the directory in which the results will be stored; default is the input directory
# 2. par: numeric vector that reports the columns (??) to select for the output; default are all the columns produced by star, once eliminated header sections

STARLog <- function(InDirectory = getwd(), OutDirectory = InDirectory, par=c(1:27), Plot=TRUE) {

# 1. Directories for tables 

# I create the directory in which store the results
if (file.exists(paste(OutDirectory, '/', 'RElaborated', sep='')) == FALSE) {
	dir.create(paste(OutDirectory, '/', 'RElaborated', sep=''))
	}
# and then the subdirectory for txt files
if (file.exists(paste(OutDirectory, '/', 'RElaborated/Tables', sep='')) == FALSE) {
	dir.create(paste(OutDirectory, '/', 'RElaborated/Tables', sep=''))
	}

# 2. List file and creation of dataframe

# Using list.files I identify all the files that are in the specified directory
AllFiles <- list.files(InDirectory)
# I then use -grep to select specifically final log file produced by Star
UsefulFiles <- AllFiles[grep('.final.out', AllFiles)]

# Cicle on input files to generate the dataframe 
for (i in 1:length(UsefulFiles)) {
	# I generate the correct file name
	fileName <- paste(InDirectory, '/', UsefulFiles[i], sep='')	
	# I generate a Lines object by reading file lines
	con <- file(fileName,open="r")
	Lines <- readLines(con)
	close(con)
	
	# I eliminate the header reads that separate each section
	SelLines <- Lines[-c(5,8,23,28)]
	SelLines <- Lines[grep('\t', Lines)]	
	# I use strsplit to split the line based on the space 
	SplitLines <- strsplit(SelLines,'\t')
	# The command retrieve a dataframe from the list; I select only the number columns
	Numbers <- do.call('rbind', lapply(SplitLines, rbind))[par,] 
	colnames(Numbers) <- c('Parameter', UsefulFiles[i])

	# To create the final dataframe:	
	if (i == 1) {
		StarNumbers <- Numbers
	} else {
		StarNumbers <- merge(StarNumbers, Numbers)
	}			
}

# 3. Adjust the obtained data frame

# I eliminate, in the Parameter column, leading spaces and | symbol 
trim.leading <- function (x)  sub("^\\s+", "", x)
StarNumbers[,1] <- trim.leading(StarNumbers[,1])
StarNumbers[,1]  <- gsub(' |', '', StarNumbers[,1], fixed=T) 
# I transpose the matrix so to have the samples as rows and the parameters as columns
StarNumbersFinal <- t(StarNumbers)

# 4. Write the dataframe on output directory
# Write the table reporting the value of the parameters for each sample in the directory:
write.table(StarNumbersFinal, file=paste(OutDirectory, '/', 'RElaborated/Tables/StarNumbers.txt', sep=''), sep='\t', quote=F, col.names=F)


# 4. Generation of plots 
if (Plot == TRUE) {

DataFrame <- as.data.frame(StarNumbersFinal[-1,])
colnames(DataFrame) <- make.names(StarNumbersFinal[1,])
Samples = rownames(DataFrame)

# 4.1 Plot for Input reads
ReadStripPlot(InputReads=as.numeric(as.character(DataFrame$Number.of.input.reads)), SampleVector=rownames(DataFrame), GroupVector=NA, OutDirectory=OutDirectory, Pdf=TRUE, Main='STAR_Input_Reads') 

# 4.2 Plot for Uniquely Mapped Reads
ReadStripPlot(InputReads=as.numeric(as.character(DataFrame$Uniquely.mapped.reads.number)), SampleVector=rownames(DataFrame), GroupVector=NA, OutDirectory=OutDirectory, Pdf=TRUE, Main='STAR_Uniquely_Mapped_Reads') 

# 4.3 Plot for Percentage of Uniquely Mapped Reads
InputPercentages = 10^(as.numeric(gsub('%', '', as.character(DataFrame$Uniquely.mapped.reads..))))
ReadStripPlot(InputReads=InputPercentages, SampleVector=rownames(DataFrame), GroupVector=NA, OutDirectory=OutDirectory, Pdf=TRUE, Main='STAR_Uniquely_Mapped_Reads_Percentage') 

}
}	



# == Example: 
InputDir <- '/Users/cheronicristina/Documents/Work/RNASeq_Interferon/RNASeq_Ens77/1.STAR_Aligned'
OutputDir <- '/Users/cheronicristina/Documents/Work/R/RPackages/FinalFunctions/AlignmentQC'
STARLog(InDirectory=InputDir, OutDirectory=OutputDir)

#Example with specified plot title
STARLog(InDirectory=InputDir, OutDirectory=OutputDir)



# ============= FINO A QUI FUNZIONA TUTTO 



# Provo basandomi sulla stessa funzione per il grafico a fare anche per gli altri tipi di parametri relativi a STAR

# MINOR PROBLEM: dovrei modificare la griglia per il grafico in percentuale




































# === 
# In generale OK, però dovrei sistemare i grafici in modo da inserirli direttamente nelle funzioni corrispondenti. 
# Posso poi lasciare anche una funzione separata che possa agire su un input specificato
# Quindi il modo migliore sarebbe utilizzare la funzione grafica creata esternamente all'interno della funzione di elaborazione del testo

