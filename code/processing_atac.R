rm(list = ls(all.names = T)); gc()
library(data.table)
library(tidyverse)

df <- fread("~/SraRunTable.txt") %>%
  dplyr::filter(gender_identity == "cis-female")



write.table(df$Run, file = "/drive-pool/data/peter_data/sc_data/brca/GSE168837/gse168837_sra.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
cat(df$Run,sep=",", file = "/drive-pool/data/peter_data/sc_data/brca/GSE168837/cr_atac_sra.txt")

## Gonna focus on this sample for nowSRR13956915 



path <- "/drive-pool/data/peter_data/sc_data/brca/GSE168837/fastq/"
files <- list.files(path = path, pattern = "SRR13956915")

## After cellranger atac count... checking if all these barcode are present and for the correct sample
barcodes <- fread(file = "~/sample345/outs/filtered_peak_bc_matrix/barcodes.tsv", header = FALSE)

## Filter data for cis female single cell data
sample_type <- c("CF", "Cis female") 

## ATAC data
atac.se <- readRDS("/drive-pool/data/peter_data/sc_data/brca/GSE168837/GSE168837_peakMatrix_update.rds")
atac.se <- atac.se[, atac.se$SampleType == sample_type]

atac.data <- as.data.frame(atac.se@colData) %>% dplyr::filter(SampleName == "CF-318-813") %>% 
  rownames_to_column(var = "sample") %>% dplyr::select(sample) #%>%
  
atac.df <- atac.data %>% separate(sample, c("id", "barcode"), sep = "#")

knott.barcode <- unique(atac.df$barcode)


table(knott.barcode %in% unique(barcodes$V1))

####################
## ArchR analysis ##
library(ArchR)
library(tidyverse)
set.seed(514)

addArchRGenome("hg38")

outdir <- "/drive-pool/data/peter_data/sc_data/brca/GSE168837/archr/"

outdir <- "~"

inputFiles <- c("/drive-pool/data/peter_data/sc_data/brca/GSE168837/cellranger_count_output/CF-0404/outs/fragments.tsv.gz",
                "/drive-pool/data/peter_data/sc_data/brca/GSE168837/cellranger_count_output/CF-1380/outs/fragments.tsv.gz",
                "/drive-pool/data/peter_data/sc_data/brca/GSE168837/cellranger_count_output/CF-2099/outs/fragments.tsv.gz",
                "/drive-pool/data/peter_data/sc_data/brca/GSE168837/cellranger_count_output/CF-2797/outs/fragments.tsv.gz",
                "/drive-pool/data/peter_data/sc_data/brca/GSE168837/cellranger_count_output/CF-4014/outs/fragments.tsv.gz",
                "/drive-pool/data/peter_data/sc_data/brca/GSE168837/cellranger_count_output/CF-428-112/outs/fragments.tsv.gz",
                "/drive-pool/data/peter_data/sc_data/brca/GSE168837/cellranger_count_output/CF-7780/outs/fragments.tsv.gz",
                "/drive-pool/data/peter_data/sc_data/brca/GSE168837/cellranger_count_output/CF-318-813/outs/fragments.tsv.gz")

names(inputFiles) <- c("CF-0404", "CF-1380", "CF-2099", "CF-2797", "CF-4014", "CF-428-112", "CF-7780", "CF-318-813")

inputFiles

ArrowFiles = createArrowFiles(
  inputFiles=inputFiles,
  sampleNames=names(inputFiles),
  minTSS=2,
  minFrags=500,
  maxFrags=1e8,
  addTileMat=TRUE, 
  addGeneScoreMat = TRUE)

projCis <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = outdir,
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

saveArchRProject(ArchRProj = projCis, outputDirectory = "Save-projCis", load = FALSE)
loadArchRProject("~/Save-projCis/")

doubScores <- addDoubletScores(
  input = projCis,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
)

## Investigate the cells with different cutoffs
peakobj = doubScores

require(reshape2)
newdf = melt(as.data.frame(peakobj@cellColData[,c("TSSEnrichment", "NucleosomeRatio", "BlacklistRatio", "DoubletEnrichment")]))

P = ggplot(newdf, aes(x=value, fill=variable), colour="black") +
  geom_histogram() +
  facet_wrap(~variable, scales="free", nrow=1) +
  scale_fill_brewer(name="", palette="Set1") +
  xlab("") +
  theme_bw(base_size=14) +
  theme(legend.position="bottom")
P

filtered.obj <- filterDoublets(doubScores)

newdf = melt(as.data.frame(filtered.obj@cellColData[,c("PromoterRatio", "NucleosomeRatio", "BlacklistRatio", "DoubletEnrichment")]))

P = ggplot(newdf, aes(x=value, fill=variable), colour="black") +
  geom_histogram() +
  facet_wrap(~variable, scales="free", nrow=1) +
  scale_fill_brewer(name="", palette="Set1") +
  xlab("") +
  theme_bw(base_size=14) +
  theme(legend.position="bottom")
P


idx_cells = which(filtered.obj@cellColData$PromoterRatio >= 0.1 &
                    filtered.obj@cellColData$PromoterRatio <= 0.8 &
                    filtered.obj@cellColData$nFrags >= 1e3 &
                    filtered.obj@cellColData$nFrags <= 5e4 &
                    filtered.obj@cellColData$NucleosomeRatio < 4)
# filtered.obj@cellColData$DoubletEnrichment < 3.5)

## Checking barcode with knott's published paper.
orig.paper <- as.data.frame(atac.se@colData) %>% rownames_to_column(var = "id") %>%
  dplyr::select(id, predictedCellType) %>%
  separate(id, c("id", "data", "barcode"), sep = "_|#")

cellid <- as.data.frame(projCis@cellColData) %>% rownames_to_column(var = "id") %>%
  dplyr::select(id) %>%
  separate(id, c("id", "barcode"), sep = "#") %>% 
  rownames_to_column(var = "row_id") %>%
  mutate(status = 1)

match.cells <- left_join(x = orig.paper, y = cellid, by = c("id", "barcode"))
table(match.cells$status, match.cells$id, useNA = c("always"))
keep.cells <- match.cells %>% dplyr::filter(!is.na(status)) %>% pull(row_id) %>% as.numeric()
## We are keeping cells that are found in the original paper

filtered.obj <- filtered.obj[idx_cells, ]

saveArchRProject(ArchRProj = filtered.obj, outputDirectory = "Save-projCis", load = FALSE)

# LSI
# What if we use peak matrix instead of TileMatrix
filtered.obj <- addIterativeLSI(
  ArchRProj = filtered.obj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 5,
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2),
    sampleCells = 10000,
    n.start = 10
  ),
  varFeatures = 25000,
  dimsToUse = 1:30,
  force=TRUE
)


batchpath = "/drive-pool/data/peter_data/sc_data/brca/GSE168837/batchinfo_v2.txt"
batchdf = fread(batchpath, header=T, sep="\t") %>% column_to_rownames(var = "Sample")
filtered.obj@cellColData$SampleName = gsub("_ATAC", "", filtered.obj@cellColData$Sample)
filtered.obj@cellColData$BatchID = batchdf[as.character(filtered.obj@cellColData$SampleName), 1]

filtered.obj <- addHarmony(
  ArchRProj = filtered.obj,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "BatchID",
  force=TRUE
)

saveArchRProject(ArchRProj = filtered.obj, outputDirectory = "Save-projCis", load = FALSE)

minDist = 0.2
nNeighbors = 60
filtered.obj <- addUMAP(
  ArchRProj = filtered.obj,
  reducedDims = "Harmony",
  name = "AfterHarmony",
  nNeighbors = nNeighbors,
  minDist = minDist,
  metric = "cosine", force=TRUE
)

metadf = as.data.frame(filtered.obj@cellColData)

saveRDS(metadf, file="metadata_archr.RDS", compress=TRUE)

## Plot umap
metadf = as.data.frame(filtered.obj@cellColData)
umapdf = as.data.frame(filtered.obj@embeddings@listData$AfterHarmony@listData$df)
colnames(umapdf) = c("UMAP1", "UMAP2")
#umapdf$CellType = metadf[rownames(umapdf), "predictedCellType"]
# textdf = do.call("rbind", lapply(unique(umapdf$CellType), function(x){
#   tempdf = umapdf[umapdf$CellType==x, ]
#   outdf = data.frame(UMAP1=mean(tempdf[,1]), UMAP2=mean(tempdf[,2]), CellType=x)
#   return(outdf)}))
# umapdf$CellType = factor(umapdf$CellType, levels=cols[,1])
# umapdf$SampleType = ifelse(grepl("CF", metadf[rownames(umapdf), "SampleName"]), "Cis female", "Trans male")
# textdf$SampleType = "Cis female"

umapdf <- umapdf %>% mutate(id = rownames(.)) %>% separate(id, c("id", "barcode"), sep = "#") %>%
  left_join(x = ., y = match.cells[, c("id", "barcode", "predictedCellType")], by = c("id", "barcode"))
test <- umapdf %>% dplyr::filter(!is.na(predictedCellType))
# cell.color <- data.frame(predictedCellType = unique(test$predictedCellType),
#                          color = c("#58986b", "#58c8c9", "grey", "#055820", 
#                                    "#758609", "#5c1505", "#bac95a", "#08305d",
#                                    "#6b3a71", "#d5626f"))
group.colors <- c("LUM_HR-pos"="#58986b", "LUM_HR-neg"="#58c8c9", "Lymphoid"="grey", "Basal"="#055820", 
                  "Myeloid"="#758609", "Blood_EC"="#5c1505", "Adipocyte"="#bac95a", "Fibroblast"="#08305d", 
                  "Lymph_EC"="#6b3a71", "Vasc.Acc."="#d5626f")

colordf <- test %>% dplyr::select(predictedCellType) %>%
  left_join(x = ., y = cell.color, by = c("predictedCellType"))
#outdir_temp = paste(outdir, "Plots/temp", sep="/")
P = ggplot(test, aes(x=UMAP1, y=UMAP2, colour = predictedCellType)) +
  geom_point(alpha=0.8) +
  scale_colour_manual(values=group.colors) +
  #geom_label_repel(data=textdf, aes(label=CellType)) +
  theme_bw(base_size=14)
P


filtered.obj <- addClusters(
  input = filtered.obj,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)

## For each cluster, re-perform LSI/UMAP and save it
saveArchRProject(ArchRProj = filtered.obj, outputDirectory = "Save-projCis", load = FALSE)

################
## Call peaks ##
################

## First we're going to filter out the single cells found in the paper. We're going to lose some samples (~1k). 
## We want to keep things consisent with the paper
cellid <- as.data.frame(filtered.obj@cellColData) %>% rownames_to_column(var = "id") %>%
  dplyr::select(id) %>%
  separate(id, c("id", "barcode"), sep = "#")
orig.paper <- as.data.frame(atac.se@colData) %>% rownames_to_column(var = "id") %>%
  dplyr::select(id, predictedCellType) %>%
  separate(id, c("id", "data", "barcode"), sep = "_|#") %>%
  mutate(status = 1) %>%
  dplyr::select(-data)
match.cells <- left_join(x = cellid, y = orig.paper, by = c("id", "barcode"))

filtered.obj@cellColData$Status <- match.cells$status
filtered.obj@cellColData$predictedCellType <- match.cells$predictedCellType
## Keep only samples with Status == 1; 1 means barcode and id matches the orig. paper
idx_cells <- which(filtered.obj@cellColData$Status == 1)
filtered.obj <- filtered.obj[idx_cells,]

## Prepare data for atac
filtered.obj@cellColData$SampleType = ifelse(grepl("CF", filtered.obj@cellColData$SampleName), "CF", "TM")
filtered.obj@cellColData$CellTypeAndSampleType = paste(filtered.obj@cellColData$predictedCellType, filtered.obj@cellColData$SampleType)
saveArchRProject(ArchRProj = filtered.obj, outputDirectory = "Save-projCis", load = FALSE)

filtered.obj <- addGroupCoverages(ArchRProj = filtered.obj, groupBy = "CellTypeAndSampleType", force=TRUE)

pathToMacs2 = "/home/p3nguyen/miniconda3/bin/macs2"

filtered.obj <- addReproduciblePeakSet(
  ArchRProj = filtered.obj,
  peakMethod="Macs2",
  additionalParams="--nomodel --nolambda --bdg --SPMR --call-summits",
  groupBy = "CellTypeAndSampleType",
  pathToMacs2 = pathToMacs2
)

peakset = getPeakSet(filtered.obj)
filtered.obj = addPeakMatrix(filtered.obj)

saveArchRProject(ArchRProj = filtered.obj, outputDirectory = "Save-projCis", load = FALSE)


peakmat <- getMatrixFromProject(filtered.obj, "PeakMatrix")
saveRDS(peakmat, file="peakMatrix.RDS", compress=TRUE)


## Perform umap and etc. with peak matrix
filtered.obj <- addIterativeLSI(
  ArchRProj = filtered.obj,
  useMatrix = "PeakMatrix",
  name = "IterativeLSIPeak",
  iterations = 5,
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2),
    sampleCells = 10000,
    n.start = 10
  ),
  varFeatures = 25000,
  dimsToUse = 1:30,
  force=TRUE
)

filtered.obj <- addHarmony(
  ArchRProj = filtered.obj,
  reducedDims = "IterativeLSIPeak",
  name = "HarmonyPeak",
  groupBy = "BatchID",
  force=TRUE
)


minDist = 0.1
nNeighbors = 60
filtered.obj <- addUMAP(
  ArchRProj = filtered.obj,
  reducedDims = "HarmonyPeak",
  name = "AfterHarmonyPeak",
  nNeighbors = nNeighbors,
  minDist = minDist,
  metric = "cosine", force=TRUE
)

filtered.obj <- addClusters(
  input = filtered.obj,
  reducedDims = "HarmonyPeak",
  method = "Seurat",
  name = "ClustersPeaks",
  resolution = 0.8, force=TRUE
)

saveArchRProject(ArchRProj = filtered.obj, outputDirectory = "Save-projCis", load = FALSE)

## Add co-accessibility and gene-links


markersPeaks <- getMarkerFeatures(
  ArchRProj = filtered.obj,
  useMatrix = "PeakMatrix",
  groupBy = "CellTypeAndSampleType",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

cutoff_condition = "FDR <= 0.05 & Log2FC >= 0.25"

markerList <- getMarkers(markersPeaks, cutOff = cutoff_condition)
#save(markerList, file = "/drive-pool/data/peter_data/sc_data/brca/GSE168837/arrowFiles/PeakCalls/markerList.rda")

dir.create(sprintf("/PeakCalls/peakTsvFiles"))
outdir <- "/drive-pool/data/peter_data/sc_data/brca/GSE168837/arrowFiles/"

for(clustername in names(markerList)){
  peakdf = as.data.frame(markerList[[clustername]])
  if(nrow(peakdf) > 0){
    outpath = sprintf("%s/PeakCalls/peakTsvFiles/%s.tsv",
                      outdir, gsub(" ", "_", clustername))
    cat(sprintf("Saving to %s\n", outpath))
    write.table(peakdf, file=outpath, sep="\t", quote=F, row.names=F)
  }
}

heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks,
  cutOff = cutoff_condition,
  transpose = TRUE
)

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 12, height = 8, ArchRProj = filtered.obj, addDOC = FALSE)

pma <- markerPlot(seMarker = markersPeaks, name = "LUM_HR-pos CF", cutOff = cutoff_condition, plotAs = "MA")

pv <- markerPlot(seMarker = markersPeaks, name = "LUM_HR-pos CF", cutOff = cutoff_condition, plotAs = "Volcano")

plotPDF(pma, pv, name = "LUM_HR-pos-CisFemale-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = filtered.obj, addDOC = FALSE)


## Linking atac marker peaks to snps
load("/drive-pool/data/peter_data/genetic_programming/code/brca/gwas_finemap.rda")
load("/drive-pool/data/peter_data/genetic_programming/code/brca/snpinLD.rda")

## Get all snps (lead snp +  ld)
all.snps <- gwas_finemap[,c("CHR_ID", "CHR_POS", "SNP_ID_CURRENT")] %>%
  dplyr::rename("SNPS"="SNP_ID_CURRENT") %>%
  rbind(.,
        ld[,c("CHR_ID", "CHR_POS", "SNPS")]) %>%
  distinct()

## Make gwas snps into granges
gwas <- GRanges(seqnames = Rle(paste0("chr", all.snps$CHR_ID)),
                range = IRanges(start = as.numeric(all.snps$CHR_POS) - 500,
                                end = as.numeric(all.snps$CHR_POS) + 499),
                locus = all.snps$SNPS,
                chr_pos = as.numeric(all.snps$CHR_POS))

load("/drive-pool/data/peter_data/sc_data/brca/GSE168837/arrowFiles/PeakCalls/markerList.rda")


celltype <- names(markerList@listData)
atac.marker <- list()

for(n in celltype) {
  marker.df <- as.data.frame(markerList@listData[[n]])
  if(nrow(marker.df) <= 1) { ## If there are no markers for celltype then skip
    next
  } else {
    p <- unlist(strsplit(x = n, split = " "))[1] ## Get only celltype name. Don't care about CF at the end
    ## Get atac-seq probes and make into granges
    marker.granges <- GRanges(seqnames = Rle(marker.df$seqnames),
                              range = IRanges(start = marker.df$start,
                                              end = marker.df$end),
                              idx = marker.df$idx)
    marker.granges <- unique(marker.granges)
    ## Find overlaps between snp and atac region
    hits <- findOverlaps(query = gwas, subject = marker.granges)
    ## Extract results
    # marker.hit <- as.data.frame(gwas) %>% .[hits@from,] %>% remove_rownames() %>% 
    #   dplyr::select(locus) %>% distinct() %>%
    #   mutate(celltype = p)
    marker.hit <- cbind(as.data.frame(gwas) %>% .[hits@from,],
                        as.data.frame(marker.granges) %>% .[hits@to,]) %>%
      .[,6:13] %>% remove_rownames() %>%
      mutate(celltype = p)
    atac.marker[[p]] <- marker.hit
  }
}
atac.marker <- as.data.frame(do.call(rbind, atac.marker)) %>% remove_rownames() %>%
  mutate(atac = 1) %>% 
  dplyr::select(locus, celltype, atac) %>%
  distinct() %>%
  left_join(x = ., 
            y = distinct(ld[,c("SNPS", "LEAD_SNP")]), 
            by = c("locus"="SNPS"))

## Anything that is NA in LEAD_SNP is column the locus is the lead snp
atac.marker[which(is.na(atac.marker$LEAD_SNP)), "LEAD_SNP"] <- atac.marker[which(is.na(atac.marker$LEAD_SNP)), "locus"] 

save(atac.marker, file = "/drive-pool/data/peter_data/genetic_programming/code/brca/ld_atac.marker.rda")

#######################################
## Finding Common Peaks in ATAC data ##
loadArchRProject("/drive-pool/data/peter_data/sc_data/brca/GSE168837/arrowFiles/Save-projCis/")
projCis <- readRDS(file = "/drive-pool/data/peter_data/sc_data/brca/GSE168837/arrowFiles/Save-projCis/Save-ArchR-Project.rds")

## Checking to see available matrix in data
getAvailableMatrices(projCis)

## Getting peaks by cell type
peakset <- getPeakSet(projCis)
names(peakset) <- NULL

## Make gwas snps into granges
gwas <- GRanges(seqnames = Rle(paste0("chr", all.snps$CHR_ID)),
                range = IRanges(start = as.numeric(all.snps$CHR_POS) - 500,
                                end = as.numeric(all.snps$CHR_POS) + 499),
                locus = all.snps$SNPS,
                chr_pos = as.numeric(all.snps$CHR_POS))

hits <- findOverlaps(query = gwas, subject = peakset)
## Extract results
# marker.hit <- as.data.frame(gwas) %>% .[hits@from,] %>% remove_rownames() %>% 
#   dplyr::select(locus) %>% distinct() %>%
#   mutate(celltype = p)
common.atac.hit <- cbind(as.data.frame(gwas) %>% .[hits@from,],
                    as.data.frame(peakset) %>% .[hits@to,]) %>% 
  remove_rownames() %>%
  dplyr::select(locus, GroupReplicate, peakType, nearestGene, distToGeneStart, distToTSS) %>% 
  mutate(GroupReplicate = gsub(pattern = "\\s.*", replacement = "", x = GroupReplicate),
         common.atac = 1) %>%
  distinct() %>%
  left_join(x = ., 
            y = distinct(ld[,c("SNPS", "LEAD_SNP")]), 
            by = c("locus"="SNPS"))

common.atac.hit[which(is.na(common.atac.hit$LEAD_SNP)), "LEAD_SNP"] <- common.atac.hit[which(is.na(common.atac.hit$LEAD_SNP)), "locus"] 

save(common.atac.hit, file = "/drive-pool/data/peter_data/genetic_programming/code/brca/ld_common_atac_peak.rda")
load("/drive-pool/data/peter_data/genetic_programming/code/brca/ld_common_atac_peak.rda")

## Test stuff to see peaks
peakset <- getPeakSet(projCis)
peakset <- peakset[names(peakset) %in% c("Fibroblast CF", "LUM_HR-neg CF"),]

rs45631563.region <- GRanges(seqnames = Rle("chr10"),
                             range = IRanges(start = 81500000,
                                             end = 86500000))


p <- plotBrowserTrack(
  ArchRProj = projCis, 
  plotSummary = c("bulkTrack", "featureTrack", "geneTrack"),
  sizes = c(10, 3, 4),
  groupBy = "CellTypeAndSampleType", 
  #geneSymbol = c("FGFR2"),
  region = rs45631563.region,
  features =  peakset,
  #upstream = 10000,
  #downstream = 10000
)

dev.off()
png(filename = "~/ld_atactest.png", height = 150, width = 150, unit = "mm", res = 300)
grid::grid.draw(p)
dev.off()

atac.peaks <- data.frame(sample = names(peakset),
                         seqnames = peakset@seqnames,
                         start = peakset@ranges@start,
                         end = peakset@ranges@start + peakset@ranges@width - 1,
                         width = peakset@ranges@width)



