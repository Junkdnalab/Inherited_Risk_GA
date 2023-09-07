rm(list = ls(all.names = T)); gc()

library(data.table)
library(tidyverse)
library(quantsmooth)

## Filtering for fine map snps
gwas_catalog <- fread(file = "/drive-pool/data/peter_data/genetic_programming/gwas/gwas_catalog_v1.0.2-associations_e109_r2023-06-03.tsv",
                      select = c("REGION", "CHR_ID", "CHR_POS", "SNPS", "CURRENT_SNP_ID","P-VALUE", "PUBMEDID", "DATE", "MAPPED_TRAIT")) %>%
  dplyr::filter(PUBMEDID %in% c(32424353, 29059683, 29058716, 23535733, 25751625, 25038754, 23544012, 20453838)) %>% ## These are all the pubmed id needed for finemap
  dplyr::rename("PVALUE" = "P-VALUE") %>%
  mutate(PVALUE = as.numeric(PVALUE))

## susceptibility variants from genome-wide analyses (obj names is finemap_pmid)
finemap_32424353 <- c("rs78378222", "rs9712235", "rs141526427", "rs6065254", "rs7760611", "rs7924772", "rs206435",
                      "rs17215231", "rs4742903", "chr1:145126177", "rs79518236", "rs138044103", "rs12962334", "rs188092014",
                      "rs1375631", "rs13039563", "rs1061657", "rs495367", "rs2464195", "rs5776993", "chr12:29140260", "rs13277568",
                      "rs9808759", "rs2886671", "rs17743054", "rs4602255", "rs142890050", "rs11652463", "rs10838267", "rs13256025",
                      "rs11065822", "rs34052812"
                      )

finemap_29059683 <- c("rs2992756", "rs4233486", "rs79724016", "rs1707302", "rs140850326", "rs17426269", "rs7529522",
                      "rs4971059", "rs35383942", "rs11117758", "rs113577745", "rs6725517", "rs71801447", "rs12479355",
                      "rs6805189", "rs13066793", "rs9833888", "rs34207738", "rs58058861", "rs6815814", "chr4:84370124",
                      "rs10022462", "rs77528541", "rs116095464", "rs72749841", "rs35951924", "rs6882649", "rs6596100",
                      "rs4562056", "rs3819405", "rs2223621", "rs71557345", "rs12207986", "rs6569648", "rs7971", "rs17156577",
                      "rs17268829", "rs71559437", "rs514192", "rs12546444", 'rs58847541', "rs1895062", "rs10760444", 
                      "rs8176636", "rs67958007", "rs140936696", "rs6597981", "rs202049448", "rs206966", "rs10623258",
                      "rs28539243", "rs2432539", "rs4496150", "rs72826962", "rs2532263", "rs117618124", "rs78269692",
                      "rs2594714", "rs2965183", "rs71338792", "rs16991615", "rs6122906", "rs738321", "rs73161324",
                      "rs28512361"
                      )

finemap_29058716 <- c("rs200648189", "rs6569648", "rs66823261", "rs17350191", "rs11374964", "rs74911261", "rs11076805",
                      "rs36194942", "rs322144", "rs113701136"
                      )

finemap_23535733 <- c("rs4245739", "rs6678914", "rs12710696", "rs11075995"
                      )
## Known finemap from 29058716 supplemental table 10. Finemap from 31 diff source
known.finemap <- c("rs616488", "rs11552449", "rs11249433", "rs12405132", "rs12048493", "rs4951011", "rs72755295",
                   "rs4849887", "rs2016394", "rs1550623", "rs1830298", "rs34005590", "rs4442975", "rs16857609",
                   "rs6762644", "rs4973768", "rs12493607", "rs6796502", "rs1053338", "rs9790517", "rs6828523", 
                   "rs3215401", "rs13162653", "rs2012709", "rs10941679", "rs62355902", "rs10472076", "rs1353747",
                   "rs7707921", "rs10474352", "rs1432679", "rs11242675", "rs9348512", "rs204247", "rs9257408",
                   "rs17529111", "rs9485372", "rs9397437", "rs6964587", "rs4593472", "rs11977670", "rs720475",
                   "rs9693444", "rs13365225", "rs6472903", "rs2943559", "rs13267382", "rs13281615", "rs11780156",
                   "rs1011970", "rs10759243", "rs676256", "rs10816625", "rs13294895", "rs2380205", "rs7072776",
                   "rs11814448", "rs10995201", "rs704010", "rs7904519", "rs11199914", "rs35054928", "rs45631563",
                   "rs2981578", "rs3817198", "rs3903072", "rs554219", "rs75915166", "rs11820646", "rs12422552",
                   "rs7297051", "rs17356907", "rs1292011", "rs11571833", "rs2236007", "rs2588809", "rs999737",
                   "rs941764", "rs11627032", "rs2290203", "rs4784227", "rs17817449", "rs13329835", "chr17:29230520",
                   "rs2787486", "rs745570", "rs527616", "rs1436904", "rs6507583", "rs4808801", "rs3760982", "rs2823093",
                   "rs17879961", "rs132390", "chr22:39359355", "rs6001930"
                   )
## Get total finemap 
total_finemap <- unique(c(finemap_23535733, finemap_29058716, finemap_29059683, finemap_32424353, known.finemap))

gwas_finemap <- gwas_catalog %>% dplyr::filter(SNPS %in% total_finemap) %>%
  rownames_to_column(var = "rowid")

## Clean up obj no longer needed
rm(finemap_23535733,finemap_29058716, finemap_29059683, finemap_32424353, known.finemap, total_finemap); gc()

## Get rid of dups from diff study. Will choose SNP from most recent study with the highest pval
dup.snps <- gwas_finemap$SNPS[which(duplicated(gwas_finemap$SNPS))] %>% unique()

for(n in dup.snps) {
  df.dup <- gwas_finemap %>% dplyr::filter(SNPS %in% n) %>% arrange(desc(DATE), PVALUE) %>% .[-1,] %>% .$rowid
  gwas_finemap <- gwas_finemap %>% dplyr::filter(!rowid %in% df.dup)
}

rm(df.dup, dup.snps, n); gc()

## There are 7 snps without chr_id chr_pos and region
gwas_finemap[which(gwas_finemap$SNPS == "rs2380205"), "REGION"] <- quantsmooth::position2Cytoband(10,5926740,units="hg38")
gwas_finemap[which(gwas_finemap$SNPS == "rs2380205"), "CHR_ID"] <- "10"
gwas_finemap[which(gwas_finemap$SNPS == "rs2380205"), "CHR_POS"] <- "5926740"

gwas_finemap[which(gwas_finemap$SNPS == "chr12:29140260"), "REGION"] <- quantsmooth::position2Cytoband(12,29140260,units="hg38")
gwas_finemap[which(gwas_finemap$SNPS == "chr12:29140260"), "CHR_ID"] <- "12"
gwas_finemap[which(gwas_finemap$SNPS == "chr12:29140260"), "CHR_POS"] <- "29140260"

gwas_finemap[which(gwas_finemap$SNPS == "chr1:145126177"), "REGION"] <- quantsmooth::position2Cytoband(1,145126177,units="hg38")
gwas_finemap[which(gwas_finemap$SNPS == "chr1:145126177"), "CHR_ID"] <- "1"
gwas_finemap[which(gwas_finemap$SNPS == "chr1:145126177"), "CHR_POS"] <- "145126177"

gwas_finemap[which(gwas_finemap$SNPS == "chr4:84370124"), "REGION"] <- quantsmooth::position2Cytoband(4,84370124,units="hg38")
gwas_finemap[which(gwas_finemap$SNPS == "chr4:84370124"), "CHR_ID"] <- "4"
gwas_finemap[which(gwas_finemap$SNPS == "chr4:84370124"), "CHR_POS"] <- "84370124"

gwas_finemap[which(gwas_finemap$SNPS == "chr22:39359355"), "REGION"] <- quantsmooth::position2Cytoband(22,39359355,units="hg38")
gwas_finemap[which(gwas_finemap$SNPS == "chr22:39359355"), "CHR_ID"] <- "22"
gwas_finemap[which(gwas_finemap$SNPS == "chr22:39359355"), "CHR_POS"] <- "39359355"

gwas_finemap[which(gwas_finemap$SNPS == "rs8176636"), "REGION"] <- "9q34.2"
gwas_finemap[which(gwas_finemap$SNPS == "rs8176636"), "CHR_ID"] <- "9"
gwas_finemap[which(gwas_finemap$SNPS == "rs8176636"), "CHR_POS"] <- "136151579"

gwas_finemap[which(gwas_finemap$SNPS == "chr17:29230520"), "REGION"] <- quantsmooth::position2Cytoband(17,29230520,units="hg38")
gwas_finemap[which(gwas_finemap$SNPS == "chr17:29230520"), "CHR_ID"] <- "17"
gwas_finemap[which(gwas_finemap$SNPS == "chr17:29230520"), "CHR_POS"] <- "29230520"

gwas_finemap <- gwas_finemap %>% dplyr::select(REGION, CHR_ID, CHR_POS, SNPS)
gwas_finemap$CHR_POS <- as.numeric(gwas_finemap$CHR_POS)
#save(gwas_catalog, gwas_finemap, file = "/drive-pool/data/peter_data/genetic_programming/code/brca/gwas_finemap.rda")

rm(list = ls(all.names = T));

library(RIdeogram)
load("/drive-pool/data/peter_data/genetic_programming/code/brca/gwas_finemap.rda")

data(human_karyotype, package="RIdeogram")
data(gene_density, package="RIdeogram")

## SNPS
#>    Type    Shape Chr    Start      End  color
#> 1  tRNA   circle   6 69204486 69204568 6a3d9a
#> 2  rRNA      box   3 68882967 68883091 33a02c
#> 3  rRNA      box   5 55777469 55777587 33a02c
#> 4  rRNA      box  21 25202207 25202315 33a02c
#> 5 miRNA triangle   1 86357632 86357687 ff7f00
#> 6 miRNA triangle  11 74399237 74399333 ff7f00
#> 

snp_marks <- data.frame(Type = "SNP",
                        Shape = "circle",
                        Chr = gwas_finemap$CHR_ID,
                        Start = gwas_finemap$CHR_POS,
                        End = gwas_finemap$CHR_POS,
                        color = "6a3d9a")


ideogram(karyotype = human_karyotype, label = snp_marks, label_type = "marker", overlaid = gene_density)
convertSVG("chromosome.svg", device = "png")

## SNPclip

#devtools::install_github("CBIIT/LDlinkR")
library(LDlinkR)

gwas_catalog[which(gwas_catalog$SNPS == "4:84370124"), "SNPS"] <- "chr4:84370124" ## Manual fix. chr missing.

## Need to fill in missing CHR_ID for filtering in SNPclip. It likes chr id to be separated
chr_row <- grep(pattern = "chr", x = gwas_catalog$SNPS)

for(n in chr_row) {
  chr_id <- gsub(pattern = ":.*", replacement = "", x = gwas_catalog[n, "SNPS"]) %>%
    gsub(pattern = ".*chr", replacement = "", x = .)
  gwas_catalog[n, "CHR_ID"] <- chr_id
}

## Manually fix these snps
snps.fix <- gwas_catalog[which(gwas_catalog$CHR_ID == ""), "SNPS"]

gwas_catalog[which(gwas_catalog$SNPS == "rs2380205"), "CHR_ID"] <- "10"
## The other 8 are unmappable to grch38. wont use in analysis
rm(snps.fix, n, chr_id); gc()


## get unique snp id
chr_id <- unique(gwas_catalog$CHR_ID) %>% .[!. %in% ""]
snpclip_res <- list()

for(n in chr_id) {
  message("Working on chr", n)
  snps <- gwas_catalog %>% dplyr::filter(CHR_ID %in% n) %>% .$SNPS
  ## Run SNPclip
  res <- SNPclip(snps = snps,
                token = "9dfb615d7a15" ## Token acquired from registering on ldlink.nih.gov
  )
  snpclip_res[[n]] <- res
}
snpclip_res <- as.data.frame(do.call(rbind, snpclip_res))

## Filter for variants that were removed due to being in LD with another snp
snpLD <- snpclip_res[grep(pattern = "Variant in LD with", x = snpclip_res$Details),] 
## String split the detail info
snpLDdetail <- as.data.frame(str_split_fixed(string = snpLD$Details, pattern = " ", n = 8)) %>% .[,c(5,6)] %>% 
  'colnames<-' (c("snp_in_ld", "r2")) %>%
  mutate(r2 = gsub(pattern = ".*=", replacement = "", x = r2)) %>% ## Fix R2 text
  mutate(r2 = gsub(pattern = ").*", replacement = "", x = r2)) %>% 
  mutate(r2 = as.numeric(r2)) ## Convert R2 character to numeric

## Joined string split details and snpclip res together
snpLD <- cbind(snpLD, snpLDdetail); rm(snpLDdetail, snpclip_res, n, res, snps); gc()
#save(snpLD, gwas_catalog, file = "/drive-pool/data/peter_data/genetic_programming/code/brca/snpLD.rda")

rm(list = ls(all.names = T)); gc()
load("/drive-pool/data/peter_data/genetic_programming/code/brca/gwas_finemap.rda")
load("/drive-pool/data/peter_data/genetic_programming/code/brca/snpLD.rda") ## Loading this will overwrite gwas_catalog (this has a more updated version of gwas_catalog)

## Map snps in ld to lead snp (206 snps)


snp_in_ld <- unique(snpLD$snp_in_ld)
table(gwas_finemap$SNPS %in% unique(snpLD$snp_in_ld)) ## 47 lead snps found not removed
table(gwas_finemap$SNPS %in% unique(snpLD$RS_Number)) ## 62 lead snps found that were removed

## Create a map file to 
test <- snpLD %>% dplyr::filter(snp_in_ld %in% gwas_finemap$SNPS)


snps <- gwas_finemap$SNPS
counter <- 1

ldproxy_res <- list()

for(snp in snps) {
  message("Working on snp",counter, ": ", snp)
  
  res <- LDproxy(
    snp = snp,
    pop = "EUR",
    r2d = "r2",
    token = "9dfb615d7a15",
    file = FALSE,
    genome_build = "grch38",
    api_root = "https://ldlink.nih.gov/LDlinkRest"
  )
  res$snp <- snp
  ## Update counter: counter is for knowing where it is.
  counter <- counter + 1
  ldproxy_res[[snp]] <- res
}
ldproxy_res <- as.data.frame(do.call(rbind, ldproxy_res))

## Identify snps that did not have any snp in ld with it
rm.res <- data.frame(ncol = sapply(ldproxy_res, ncol)) %>% dplyr::filter(ncol == 2) %>% rownames()
ldproxy_res <- ldproxy_res[!names(ldproxy_res) %in% rm.res]

ldproxy_res <- as.data.frame(do.call(rbind, ldproxy_res)) %>% remove_rownames()
#save(ldproxy_res, file = "/drive-pool/data/peter_data/genetic_programming/code/brca/snpLDproxy.rda")

load("/drive-pool/data/peter_data/genetic_programming/code/brca/snpLDproxy.rda")
## Filter ld proxy for snps with an R2 of 0.9 and MAF of 0.05 or greater
ld <- ldproxy_res %>% dplyr::filter(R2 >= 0.9,
                                    MAF >= 0.05) %>%
  dplyr::select(RS_Number, Coord, snp)

## There are cases where the RS_number and snp (lead snp) are the same. Remove those cases.
ld <- ld[-which(ld$RS_Number == ld$snp),] 

## Are there any lead snps in ld with another snp
table(unique(ld$snp) %in% ld$RS_Number) ## No so we're good.

length(unique(ld$snp)) ## There's 139 lead snps that have ld with other snps

## Some snps dont have rsid, just coord. Lets add corrd to rsid for those cases
ld[which(ld$RS_Number == "."), "RS_Number"] <- ld[which(ld$RS_Number == "."), "Coord"]

ld <- ld %>% separate(Coord, c("CHR_ID", "CHR_POS"), sep = ":") %>%
  mutate(CHR_ID = gsub(pattern = ".*chr", replacement = "", x = CHR_ID)) %>%
  dplyr::rename("SNPS"="RS_Number", "LEAD_SNP"="snp")

# ld$snp.cytoband <- apply(X = ld, MARGIN = 1, FUN = function(x) {
#   quantsmooth::position2Cytoband(x["CHR_ID"], as.numeric(x["CHR_POS"]), units="hg38")
# })

ld.gwas <- GRanges(seqnames = Rle(paste0("chr", ld$CHR_ID)),
                range = IRanges(start = as.numeric(ld$CHR_POS) - 500,
                                end = as.numeric(ld$CHR_POS) + 499),
                locus = ld$SNPS,
                chr_pos = as.numeric(ld$CHR_POS))

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
    hits <- findOverlaps(query = ld.gwas, subject = marker.granges)
    ## Extract results
    # marker.hit <- as.data.frame(gwas) %>% .[hits@from,] %>% remove_rownames() %>% 
    #   dplyr::select(locus) %>% distinct() %>%
    #   mutate(celltype = p)
    marker.hit <- cbind(as.data.frame(ld.gwas) %>% .[hits@from,],
                        as.data.frame(marker.granges) %>% .[hits@to,]) %>%
      .[,c(6:13)] %>% remove_rownames() %>%
      mutate(celltype = p)
    atac.marker[[p]] <- marker.hit
  }
}
## convert list into data.frame and select relevant columns
atac.marker <- as.data.frame(do.call(rbind, atac.marker)) %>% remove_rownames() %>%
  mutate(atac = 1) %>% 
  dplyr::select(locus, celltype, atac) %>%
  distinct()

## add lead snp info so we can check if any new atac markers found that wasn't found using lead snp
ld.atac.marker <- atac.marker %>% left_join(x = ., y = ld[, c("SNPS", "LEAD_SNP")], 
                                            by = c("locus" = "SNPS"),
                                            relationship = "many-to-many")

load("/drive-pool/data/peter_data/genetic_programming/code/brca/atac.marker.rda")
test <- left_join(x = ld.atac.marker, y = atac.marker, by = c("LEAD_SNP"="locus"))

length(unique(test$LEAD_SNP))
length(unique(atac.marker$locus))


#########################################
## Prep data for pathogenicity calling ##

## Load lead snp info
## We need to get the ref and atl allele for the finemap snp. ldproxy gave has that info for its res. 
rm(list = ls(all.names = T)); gc()
gwas.full <- fread(file = "/drive-pool/data/peter_data/genetic_programming/gwas/gwas_catalog_v1.0.2-associations_e109_r2023-06-03.tsv")
load("/drive-pool/data/peter_data/genetic_programming/code/brca/gwas_finemap.rda"); rm(gwas_catalog)

## Sometimes the finemap snp has merged with another rsid and so we need to get the updated snp
snps <- gwas.full %>% dplyr::filter(SNPS %in% gwas_finemap$SNPS,
                                    MAPPED_TRAIT %in% c("breast carcinoma", 
                                                        "estrogen-receptor negative breast cancer",
                                                        "triple-negative breast cancer")) %>% ## Filter gwas.full for finemap snps
  dplyr::select(SNPS, SNP_ID_CURRENT, CONTEXT, MAPPED_TRAIT, DATE, PUBMEDID, STUDY, LINK) %>%  ## Keep relevant column
  distinct() ## Keep unique observations

unique(snps$MAPPED_TRAIT)

unique.snps <- unique(snps$SNPS)

snp.recent.pubmed <- list()
for(n in unique.snps) {
  df <- snps %>% dplyr::filter(SNPS %in% n) %>%
    .[order(as.Date(.$DATE, format="%Y/%m/%D"), decreasing = TRUE),]
  snp.recent.pubmed[[n]] <- df[1,]
}
snp.recent.pubmed <- as.data.frame(do.call(rbind, snp.recent.pubmed))

snps.rsid <- grep(pattern = "rs", x = snps$SNPS)
snps.chr <- grep(pattern = "chr", x = snps$SNPS)

## Add rs in front of the current snp id
for(rsid in snps.rsid) {
  snps[rsid, "SNP_ID_CURRENT"] <- paste0("rs", snps[rsid,"SNP_ID_CURRENT"])
}
## Copy original snp (chr?:????) to SNP_ID_CURRENT. It's missing
for(schr in snps.chr) {
  snps[schr, "SNP_ID_CURRENT"] <- snps[schr,"SNPS"]
}
rm(snps.chr, schr, rsid, snps.rsid); gc()

rsid <- snps[grep(pattern = "rs", x = snps$SNP_ID_CURRENT),]$SNP_ID_CURRENT

library(BSgenome.Hsapiens.UCSC.hg38)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)

snps.mb <- snps.from.rsid(rsid,
                          dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh38,
                          search.genome = BSgenome.Hsapiens.UCSC.hg38)


dbsnp <- SNPlocs.Hsapiens.dbSNP155.GRCh38 ## Load dbsnp data

df.snps <- unique.snps[grep(pattern = "rs", x = unique.snps)]
dbsnp.res <- snpsById(x = dbsnp, ids = df.snps, ifnotfound = "drop")
  
  ## have to convert to UCSC style chromosomes
seqlevelsStyle(dbsnp.res) <- "UCSC"
z <- inferRefAndAltAlleles(dbsnp.res, BSgenome.Hsapiens.UCSC.hg38)
mcols(dbsnp.res) <- cbind(mcols(dbsnp.res), z)
dbsnp.res

missing.snp <- dbsnp.res$RefSNP_id

snps.mb <- snps.from.rsid(rsid,
                          dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh38,
                          search.genome = BSgenome.Hsapiens.UCSC.hg38)

## Create bed file for chr cases
snps.chr <- snps[grep(pattern = "chr", x = snps$SNPS),]$SNPS
test.df <- gwas.full %>% dplyr::filter(SNPS %in% snps.chr)
snps.bed.file <- data.frame(chrom = gsub(pattern = ":.*", replacement = "", x = snps.chr),
                            chromStart = as.integer(gsub(pattern = ".*:", replacement = "", x = snps.chr)) - 1,
                            chromEnd = as.integer(gsub(pattern = ".*:", replacement = "", x = snps.chr)),
                            name = paste0(snps.chr),
                            score = 0,
                            strand = "+")

df <- GRanges(snps.bed.file[,1], IRanges(snps.bed.file[,2], snps.bed.file[,3]), snps.bed.file$strand,
              name = snps.bed.file[,4],
              score = snps.bed.file[,5])
genome(df) <- "hg38"

snps.bed.file <- df
export(snps.bed.file, "~/file.bed", "bed")

#write.table(snps.bed.file, "~/file.bed", quote = F)
snps.bed.file <- "~/file.bed"

read.table("~/file.bed")
rtracklayer::import("~/file.bed", format = "bed")

rtracklayer::import(snps.bed.file)

snps.mb <- snps.from.file(snps.bed.file,
                          search.genome = BSgenome.Hsapiens.UCSC.hg38,
                          format = "bed")


## Load snps in ld with lead snp
load("/drive-pool/data/peter_data/genetic_programming/code/brca/snpLDproxy.rda")
## Filter ld proxy for snps with an R2 of 0.9 and MAF of 0.05 or greater
ld <- ldproxy_res %>% dplyr::filter(R2 >= 0.9,
                                    MAF >= 0.05) 
ld <- ld[-which(ld$RS_Number == ld$snp),] ## Get rid of observations where RS_Number and snp are the same. 
rm(ldproxy_res); gc()



ld_current_snp <- all_gwas %>% dplyr::filter(SNPS %in% unique(ld$RS_Number))

## Fixing deletions ##
deletions <- which(for_caddcapice$Tumor_Seq_Allele2 == "-") ## Which observations have deletions?

for(n in deletions) {
  print(n)
  seq.length <- nchar(for_caddcapice[n, "Reference_Allele"]) - 1 ## Get number of dna sequence
  variant.coord <- paste0("chr", for_caddcapice[n,"Chromosome"],":", for_caddcapice[n,"Start_Position"]-1,"-",for_caddcapice[n,"Start_Position"] + seq.length) ## Create coordinate for browser
  variant.grange <- GRanges(variant.coord)
  variant.refseq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, variant.grange) ## Get dna seq
  
  ## Update for_caddcapice df
  for_caddcapice[n,"Start_Position"] <- for_caddcapice[n,"Start_Position"]-1 ## Update pos with extra nucleotide prior to deleterion
  for_caddcapice[n, "Reference_Allele"] <- as.character(variant.refseq[[1]]) ## Nucleotide seq with 1 nucleo before
  for_caddcapice[n, "Tumor_Seq_Allele2"] <- as.character(variant.refseq[[1]][1]) ## Nucleotide seq after deletion
}


