####################
## ArchR analysis ##
library(ArchR)
library(data.table)
library(tidyverse)
set.seed(514)

## Specify genome 
addArchRGenome("hg38")

outdir <- "/drive-pool/data/peter_data/sc_data/brca/GSE168837/archr/" ## output dir

## Fragment files obtained cell ranger. Check cellranger_count.py
inputFiles <- c("/drive-pool/data/peter_data/sc_data/brca/GSE168837/cellranger_count_output/CF-0404/outs/fragments.tsv.gz",
                "/drive-pool/data/peter_data/sc_data/brca/GSE168837/cellranger_count_output/CF-1380/outs/fragments.tsv.gz",
                "/drive-pool/data/peter_data/sc_data/brca/GSE168837/cellranger_count_output/CF-2099/outs/fragments.tsv.gz",
                "/drive-pool/data/peter_data/sc_data/brca/GSE168837/cellranger_count_output/CF-2797/outs/fragments.tsv.gz",
                "/drive-pool/data/peter_data/sc_data/brca/GSE168837/cellranger_count_output/CF-4014/outs/fragments.tsv.gz",
                "/drive-pool/data/peter_data/sc_data/brca/GSE168837/cellranger_count_output/CF-428-112/outs/fragments.tsv.gz",
                "/drive-pool/data/peter_data/sc_data/brca/GSE168837/cellranger_count_output/CF-7780/outs/fragments.tsv.gz",
                "/drive-pool/data/peter_data/sc_data/brca/GSE168837/cellranger_count_output/CF-318-813/outs/fragments.tsv.gz")

## Sample name
names(inputFiles) <- c("CF-0404", "CF-1380", "CF-2099", "CF-2797", "CF-4014", "CF-428-112", "CF-7780", "CF-318-813")

## Create arrowfiles required for ArchR
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
  copyArrows = TRUE #This is recommended so that if you modify the Arrow files you have an original copy for later usage.
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

## Before filtering plot
# require(reshape2)
# newdf = melt(as.data.frame(peakobj@cellColData[,c("TSSEnrichment", "NucleosomeRatio", "BlacklistRatio", "DoubletEnrichment")]))
# 
# P = ggplot(newdf, aes(x=value, fill=variable), colour="black") +
#   geom_histogram() +
#   facet_wrap(~variable, scales="free", nrow=1) +
#   scale_fill_brewer(name="", palette="Set1") +
#   xlab("") +
#   theme_bw(base_size=14) +
#   theme(legend.position="bottom")
# P

filtered.obj <- filterDoublets(doubScores)

## After filtering plot 
# newdf = melt(as.data.frame(filtered.obj@cellColData[,c("PromoterRatio", "NucleosomeRatio", "BlacklistRatio", "DoubletEnrichment")]))
# 
# P = ggplot(newdf, aes(x=value, fill=variable), colour="black") +
#   geom_histogram() +
#   facet_wrap(~variable, scales="free", nrow=1) +
#   scale_fill_brewer(name="", palette="Set1") +
#   xlab("") +
#   theme_bw(base_size=14) +
#   theme(legend.position="bottom")
# P

## Additional filters
idx_cells = which(filtered.obj@cellColData$PromoterRatio >= 0.1 &
                    filtered.obj@cellColData$PromoterRatio <= 0.8 &
                    filtered.obj@cellColData$nFrags >= 1e3 &
                    filtered.obj@cellColData$nFrags <= 5e4 &
                    filtered.obj@cellColData$NucleosomeRatio < 4)


## Checking barcode with knott's published paper.
## ATAC data
## Acquired peakMatrix from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168837
atac.se <- readRDS("/drive-pool/data/peter_data/sc_data/brca/GSE168837/GSE168837_peakMatrix_update.rds")
atac.se <- atac.se[, atac.se$SampleType == sample_type]

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

## Plot umap to compare between ours and knott's data (checking reproducibility)
metadf = as.data.frame(filtered.obj@cellColData)
umapdf = as.data.frame(filtered.obj@embeddings@listData$AfterHarmony@listData$df)
colnames(umapdf) = c("UMAP1", "UMAP2")


umapdf <- umapdf %>% mutate(id = rownames(.)) %>% separate(id, c("id", "barcode"), sep = "#") %>%
  left_join(x = ., y = match.cells[, c("id", "barcode", "predictedCellType")], by = c("id", "barcode"))
umapdf <- umapdf %>% dplyr::filter(!is.na(predictedCellType))

group.colors <- c("LUM_HR-pos"="#58986b", "LUM_HR-neg"="#58c8c9", "Lymphoid"="grey", "Basal"="#055820", 
                  "Myeloid"="#758609", "Blood_EC"="#5c1505", "Adipocyte"="#bac95a", "Fibroblast"="#08305d", 
                  "Lymph_EC"="#6b3a71", "Vasc.Acc."="#d5626f")

colordf <- umapdf %>% dplyr::select(predictedCellType) %>%
  left_join(x = ., y = cell.color, by = c("predictedCellType"))

P = ggplot(umapdf, aes(x=UMAP1, y=UMAP2, colour = predictedCellType)) +
  geom_point(alpha=0.8) +
  scale_colour_manual(values=group.colors) +
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

## Used macs2 to find peaks
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

#########################
## Get MarkerATAC list ## Used for isMarkerATAC
#########################
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

## Linking atac marker peaks to snps
load("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/gwas_snps.rda")
load("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/snpinLD.rda")

## Get all snps (lead snp +  ld)
all.snps <- gwas_snps[,c("CHR_ID", "CHR_POS", "SNP_ID_CURRENT")] %>%
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

## Get list of cell types
celltype <- names(markerList@listData)
atac.marker <- list() ## Store result from for loop below

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

save(atac.marker, file = "/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/ld_atac.marker.rda")

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
## Find overlaps between snp and atac region
hits <- findOverlaps(query = gwas, subject = peakset)
## Extract results
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

save(common.atac.hit, file = "/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/ld_common_atac_peak.rda")
load("/drive-pool/data/peter_data/repo/genetic_algorithm_gwas_sc/gwas_gp_manu_data/ld_common_atac_peak.rda")

