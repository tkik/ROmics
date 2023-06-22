##### preliminaries
library(RnBeads)


#GEO_ACCESSION<-"GSE29290"
#GEO_URL<-"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29290"

###### load data from GEO

rnb.options(disk.dump.big.matrices = FALSE)

#rnb.set<-rnb.read.geo("GSE29290")


##if we start with idat files ($this is a different one!)

#data <- getGEO('GSE175758')
# metadata <-  data$GSE175758_series_matrix.txt.gz@phenoData@data
# metadata$title <- gsub("(sample )([[:digit:]]{5}) (genomic DNA from .+)", "\\2", metadata$title)
# metadata <- metadata[c("GSM5345770", "GSM5345771", "GSM5345772", "GSM5345830", "GSM5345831", "GSM5345832"),]
# metadata  <- metadata[,c("title", "geo_accession", "day:ch1", "gender:ch1")]
# metadata$barcode <- c("8622007099_R01C01", "8622007099_R03C01", "3998998208_R02C02", "8622007081_R02C01", "8622007082_R06C01", "3999078058_R05C01")
# write.csv(metadata, "data/sample_annotation.csv", row.names=F)
#idat.dir <- "data/idat"
#sample.annotation <- "data/sample_annotation.csv"

#data.source <- c(idat.dir, sample.annotation)

#result <- rnb.execute.import(data.source=data.source,
#                         data.type="infinium.idat.dir")

#saveRDS(result, "data/raw_rnbset_GSE175758.RDS")

rnb.options(analysis.name = "Re-analysis of HCT116 450k data from Dedeurwaerder et al.")
rnb.options(identifiers.column = "Sample_ID")

non_hct116<-which(!pheno(rnb.set)[["source_name_ch1"]] %in% "Colon")
rnb.set<-remove.samples(rnb.set, non_hct116)

## clean up the phenotypic table

pheno(rnb.set)
rnb.set@pheno<-pheno(rnb.set)[,sapply(pheno(rnb.set), function(x) length(unique(x)))>1]
rnb.set@pheno<-pheno(rnb.set)[,c(1:3)]
colnames(rnb.set@pheno)<-c("Sample_Title", "GEO_accession", "Group")

rnb.set@pheno$Replicate<-gsub("(genomic DNA of HCT116 )(WT|DKO)( cells \\()(r[[:digit:]])\\)", "\\4", pheno(rnb.set)$Sample_Title)

rnb.set@pheno$Sample_Group<-rep("WT", nrow(pheno(rnb.set)))
rnb.set@pheno$Sample_Group[grep("DKO", pheno(rnb.set)[["Sample_Title"]])]<-"DKO"
rnb.set@pheno$Sample_ID<-paste(sep="_",  pheno(rnb.set)[["Sample_Group"]], pheno(rnb.set)[["Replicate"]])
rnb.set@pheno<-pheno(rnb.set)[c("GEO_accession", "Sample_ID", "Sample_Group", "Replicate")]
pheno(rnb.set)

report.dir <- "output/diffmeth"

num.cores <- 6
parallel.setup(num.cores)

saveRDS(rnb.set, "data/raw_rnb.set.RDS")

rnb.run.qc(rnb.set, dir.reports = report.dir)

rnb.options(normalization.method = "bmiq")
rnb.options(normalization.background.method = "none")
rnb.options(filtering.coverage.threshold = 3)


result <- rnb.run.preprocessing(rnb.set, dir.reports=report.dir)
rnb.run.exploratory(rnb.set, dir.reports = report.dir)

rnb.set <- result$rnb.set
saveRDS(rnb.set, "data/normalized_rnb.set.RDS")


rnb.options(differential.comparison.columns = c("Sample_Group"))
rnb.options(columns.pairing = c("Sample_Group"="Replicate"))


 rnb.run.differential(rnb.set, dir.reports=report.dir)


rnb.set <- rnb.execute.age.prediction(rnb.set)
immune.content <- rnb.execute.lump(rnb.set)


rnb.options(identifiers.column = "geo_accession")
library("minfiData")
library(conumee)
data(MsetEx)
MsetEx
minfi.data <- CNV.load(MsetEx)
minfi.controls <- pData(MsetEx)$status == "normal"

rnb.set <- readRDS("data/raw_rnbset_GSE175758.RDS")
intensity <- M(rnb.set, row.names = T) + U(rnb.set, row.names = T)
intensity_control <- getMeth(MsetEx) + getUnmeth(MsetEx)
overlap <- intersect(rownames(intensity), rownames(intensity_control))
intensity <- intensity[overlap,]
intensity_control <- intensity_control[overlap, minfi.controls]


intensity <- apply(intensity, 2, as.numeric)
intensity <- as.data.frame(cbind(intensity, intensity_control))
cnv.data <- CNV.load(intensity)


data(exclude_regions)
data(detail_regions)
head(detail_regions, n = 2)

anno <- CNV.create_anno(array_type = "450k", exclude_regions = exclude_regions, detail_regions = detail_regions)

x <- CNV.fit( query = cnv.data[1], ref=cnv.data[7:9], anno)

x <- CNV.bin(x)
x <- CNV.detail(x)
x <- CNV.segment(x)

CNV.genomeplot(x)

CNV.genomeplot(x, chr = "chr6")
CNV.detailplot(x, name = "PTEN")
CNV.detailplot_wrap(x)

x <- CNV.fit( query = cnv.data[2], ref=cnv.data[7:9], anno)

x <- CNV.bin(x)
x <- CNV.detail(x)
x <- CNV.segment(x)

CNV.genomeplot(x)

CNV.genomeplot(x, chr = "chr6")
CNV.detailplot(x, name = "PTEN")
CNV.detailplot_wrap(x)
