library(Rsubread)
library(edgeR)
library(limma)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(gplots)
library(ggplot2)

#works for mouse and human only
species <- readline(prompt = "Enter species name: ")

#load metadata from working directory
print("Select the metadata")
sampleinfo <- read.delim(file.choose())

#explore the metadata
head(sampleinfo)

#downloading required files
if (species == "mouse" | species == "Mouse") {link <- "ftp://ftp.ensembl.org/pub/release-101/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz"
download.file(link,"mousegenome.fa.gz")
} else if (species == "human" | species == "Human")
{link <- "ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
download.file(link,"humangenome.fa.gz")
} else {
  print("some other species provided")
}

#loads all fastq files in working directory
fastq.files <- list.files(pattern = ".fastq.gz$", full.names = TRUE)


#building index
if (species == "mouse" | species == "Mouse") {buildindex(basename = "mm10", reference = "mousegenome.fa.gz")
} else if (species == "human" | species == "Human")
{buildindex(basename = "hg19", reference = "humangenome.fa.gz")
} else {
  print("some error")
}


#aligning reads to reference genome
if (species == "mouse" | species == "Mouse") {align(index="mm10",readfile1=fastq.files)
} else if (species == "human" | species == "Human")
{align(index="hg19",readfile1=fastq.files)
} else {
  print("some error")
}

#BAM files saved in directory
bam.files <- list.files(pattern = ".BAM$", full.names = TRUE)

#generating count matrix
if (species == "mouse" | species == "Mouse") {fc <- featureCounts(bam.files, annot.inbuilt="mm10")
} else if (species == "human" | species == "Human")
{fc <- featureCounts(bam.files, annot.inbuilt="hg19")
} else {
  print("some error")
}

countdata <- fc$counts
colnames(countdata) <- substr(colnames(countdata),start=1,stop=10)

#making DGEList object
y <- DGEList(countdata)

#grouping data (specific to metadata; needs to be manually done)
group <- paste(sampleinfo$CellType,sampleinfo$Status,sep=".")
group <- factor(group)

#adding group info to DGElist
y$samples$group <- group


#adding gene annotation
if (species == "mouse" | species == "Mouse") {ann <- select(org.Mm.eg.db,keys=rownames(y$counts),columns=c("ENTREZID","SYMBOL","GENENAME"))
y$genes <- ann
} else if (species == "human" | species == "Human")
{ann <- select(org.Hs.eg.db,keys=rownames(y$counts),columns=c("ENTREZID","SYMBOL","GENENAME"))
y$genes <- ann
} else {
  print("some error")
}


#obtaining counts per million (cpm)
CPM <- cpm(countdata)

#cpm greater than 0.5 are kept
thresh <- CPM > 0.5
keep <- rowSums(thresh) >= 2

#cpm plot
par(mfrow=c(1,2))
plot(CPM[,1],countdata[,1])
abline(v=0.5)
# limiting the x and y-axis values to see what is happening at smaller counts
plot(CPM[,1],countdata[,1],ylim=c(0,50),xlim=c(0,3))
abline(v=0.5)


#bar graph of library sizes
par(mfrow=c(1,1))
barplot(y$samples$lib.size/1e06, names=colnames(y), las=2, ann=FALSE, cex.names=0.75)
mtext(side = 1, text = "Samples", line = 4)
mtext(side = 2, text = "Library size (millions)", line = 3)
title("Barplot of library sizes")


# log2 counts per million
logcounts <- cpm(y,log=TRUE)

# Checking distributions of raw counts
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)

# adding a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

#mds plot
plotMDS(y)


#variance for each row in logcount matrix
var_genes <- apply(logcounts, 1, var)

#getting gene names for top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]

#subset logcount matrix
highly_variable_lcpm <- logcounts[select_var,]

#heatmap
colnames(highly_variable_lcpm)<- group
heatmap(highly_variable_lcpm, col = topo.colors(50), margin=c(10,6))

# Calculating mean for each gene
mean_readCounts <- apply(countdata[,1:3], 1, mean)

# Calculating variance for each gene
var_readCounts <- apply(countdata[,1:3], 1, var)

#dispersion plot
df <- data.frame(mean_readCounts, var_readCounts)
ggplot(df) +
  geom_point(aes(x=mean_readCounts, y= var_readCounts)) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("Mean counts per gene") +
  ylab("Variance per gene") +
  labs(title = "Dispersion")

#TMM normalisation
y <- calcNormFactors(y)

#creating the design matrix without intercept
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

#voom transform
v <- voom(y,design,plot = TRUE)

#logcpm plots (unnormalised & voom transformed)
par(mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
#adding a blue line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
#adding a blue line that corresponds to the median logCPM
abline(h=median(v$E),col="blue")


#fitting to a linear model
fit <- lmFit(v)

#comparing between groups (the comparison groups needs to be manually selected)
cont.matrix <- makeContrasts(diffexp = luminal.pregnant - luminal.lactate,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)

#estimation of t-statistics and p values
fit.cont <- eBayes(fit.cont)
summa.fit <- decideTests(fit.cont)

par(mfrow=c(1,1))
#MA plot
plotMD(fit.cont,coef=1,status=summa.fit[,"diffexp"], values = c(-1, 1), hl.col=c("blue","red"))

#volcano plot for top 100 DE genes
volcanoplot(fit.cont,coef=1,highlight=100,names=fit.cont$genes$SYMBOL, main="diffexp")

