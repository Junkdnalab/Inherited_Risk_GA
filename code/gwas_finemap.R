rm(list = ls(all.names = T)); gc()

library(data.table)
library(tidyverse)
library(quantsmooth)

## Filtering for fine map snps
gwas_catalog <- fread(file = "/drive-pool/data/peter_data/genetic_programming/gwas/gwas_catalog_v1.0.2-associations_e109_r2023-06-03.tsv",
                      select = c("REGION", "CHR_ID", "CHR_POS", "SNPS", "SNP_ID_CURRENT", "P-VALUE", "PUBMEDID", "DATE", "MAPPED_TRAIT")) %>%
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

## Update SNP_ID_CURRENT. Add rs to front of rsid
snps.rsid <- grep(pattern = "rs", x = gwas_finemap$SNPS)
snps.chr <- grep(pattern = "chr", x = gwas_finemap$SNPS)

## Add rs in front of the current snp id
for(rsid in snps.rsid) {
  gwas_finemap[rsid, "SNP_ID_CURRENT"] <- paste0("rs", gwas_finemap[rsid,"SNP_ID_CURRENT"])
}
## Copy original snp (chr?:????) to SNP_ID_CURRENT. It's missing
for(schr in snps.chr) {
  gwas_finemap[schr, "SNP_ID_CURRENT"] <- gwas_finemap[schr,"SNPS"]
}
rm(snps.chr, schr, rsid, snps.rsid); gc()

## How many snps name were changed b/c of the SNP_ID_CURRENT?
table(gwas_finemap$SNPS == gwas_finemap$SNP_ID_CURRENT) ## 6 snps were changed

gwas_finemap$CHR_POS <- as.numeric(gwas_finemap$CHR_POS)
#save(gwas_catalog, gwas_finemap, file = "/drive-pool/data/peter_data/genetic_programming/code/brca/gwas_finemap.rda")

rm(list = ls(all.names = T)); gc()

###################################
## Fix snps that don't have rsid ##
###################################
load("/drive-pool/data/peter_data/genetic_programming/code/brca/gwas_finemap.rda")
## Get all current snp id
snps <- gwas_finemap$SNP_ID_CURRENT
no.rsid <- snps[grep(pattern = "chr", x = snps)] ## There are 5 snps that don't have rsid. lets see why.

## Used ncbi variation viewer to track the rsid 
## Prob not found because these are hg19 coord and gwas catalog never liftover to hg38
#chr12:29140260 rs1027113
#chr1:145126177 rs587712509
#chr4:84370124 rs1721204963
#chr22:39359355 rs868638441
#chr17:29230520 rs62070644 (Used PMID 25751625 supplemental to get this rsid - correlated 56k bp away from each other)

updated.rsid <- data.frame(SNPS = c("chr12:29140260", "chr1:145126177", "chr4:84370124", "chr22:39359355", "chr17:29230520"),
                           SNPS_updated = c("rs1027113", "rs587712509", "rs1721204963", "rs868638441", "rs62070644"))
## Using snps.from.rsid to get coord to update gwas_finemap
snps.mb <- snps.from.rsid(rsid = updated.rsid$SNPS_updated,
                          dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh38,
                          search.genome = BSgenome.Hsapiens.UCSC.hg38)
names(snps.mb) <- NULL
snps.mb <- as.data.frame(snps.mb)
snps.mb <- snps.mb %>% dplyr::select(seqnames, start, SNP_id) %>%
  distinct() %>%
  dplyr::rename("SNPS_updated"="SNP_id") %>%
  left_join(x = ., y = updated.rsid, by = "SNPS_updated")
rm(snps, no.rsid); gc()

gwas_finemap <- gwas_finemap %>% left_join(x = ., y = snps.mb, by = "SNPS")

#save(gwas_catalog, gwas_finemap, file = "/drive-pool/data/peter_data/genetic_programming/code/brca/gwas_finemap.rda")

## Find which rows need to have the SNP_ID_CURRENT and CHR_POS updated
chr_to_rsid <- which(!is.na(gwas_finemap$SNPS_updated)) ## Should be only 5 b/c we had 5 chr ids
## for loop to apply update
for(n in chr_to_rsid) {
  ## update snp_id_current
  gwas_finemap[n, "SNP_ID_CURRENT"] <- gwas_finemap[n, "SNPS_updated"]
  ## update chr_pos
  gwas_finemap[n, "CHR_POS"] <- gwas_finemap[n, "start"]
}
gwas_finemap <- gwas_finemap %>% dplyr::select(-seqnames, -start, -SNPS_updated)

library(RIdeogram)
load("/drive-pool/data/peter_data/genetic_programming/code/brca/gwas_finemap.rda")

## Load data for drawing ideogram
data(human_karyotype, package="RIdeogram")
data(gene_density, package="RIdeogram")

## Create df for marking snps on the ideogram
snp_marks <- data.frame(Type = "SNP",
                        Shape = "circle",
                        Chr = gwas_finemap$CHR_ID,
                        Start = gwas_finemap$CHR_POS,
                        End = gwas_finemap$CHR_POS,
                        color = "6a3d9a")


ideogram(karyotype = human_karyotype, label = snp_marks, label_type = "marker", overlaid = gene_density)
convertSVG("chromosome.svg", device = "png")

rm(gene_density, gwas_catalog, human_karyotype, snp_marks); gc()

#devtools::install_github("CBIIT/LDlinkR")
library(LDlinkR)

snps <- gwas_finemap$SNP_ID_CURRENT
counter <- 1
## Creating list to store res from for loop below. LDproxy does not like multi queries
ldproxy_res <- list()

## For loop to find snps in proxy to the lead snp
for(snp in snps) {
  message("Working on snp",counter, ": ", snp)
  
  res <- LDproxy(
    snp = snp,
    pop = c("CEU", "TSI", "GBR", "IBS"), ## All european pop
    r2d = "r2",
    token = "9dfb615d7a15", ## Unique token. Need to register on ldlink to get your own token
    file = FALSE,
    genome_build = "grch38",
    api_root = "https://ldlink.nih.gov/LDlinkRest"
  )
  res$snp <- snp
  ## Update counter: counter is for knowing where it is.
  counter <- counter + 1
  ldproxy_res[[snp]] <- res
}

## Identify snps that did not have any snp in ld with it
rm.res <- data.frame(ncol = sapply(ldproxy_res, ncol)) %>% dplyr::filter(ncol == 2) %>% rownames()
ldproxy_res <- ldproxy_res[!names(ldproxy_res) %in% rm.res]

ldproxy_res <- as.data.frame(do.call(rbind, ldproxy_res)) %>% remove_rownames()
#save(ldproxy_res, file = "/drive-pool/data/peter_data/genetic_programming/code/brca/snpLDproxy.rda")

rm(list = ls(all.names = T)); gc()
load("/drive-pool/data/peter_data/genetic_programming/code/brca/snpLDproxy.rda")
## Filter ld proxy for snps with an R2 of 0.9 and MAF of 0.05 or greater
ld <- ldproxy_res %>% dplyr::filter(R2 >= 0.6,
                                    Dprime >= 0.9) %>%
  dplyr::select(RS_Number, Coord, snp)

## There are cases where the RS_number and snp (lead snp) are the same. Dont know why, but remove those cases.
ld <- ld[-which(ld$RS_Number == ld$snp),] 

## Are there any lead snps in ld with another snp
table(unique(ld$snp) %in% ld$RS_Number) ## No so we're good.

plot.n.snps <- as.data.frame(table(ld$snp)) ## There's 174 lead snps that have ld with other snps
ggplot(data = plot.n.snps, aes(x = Freq)) +
  geom_histogram() +
  theme_minimal() +
  xlab("SNPs in LD") + ylab("Count") +
  ggtitle("Distribution of SNPs in LD for 174 breast cancer risk SNPs")


## Some snps dont have rsid, just coord. Lets add corrd to rsid for those cases
ld[which(ld$RS_Number == "."), "RS_Number"] <- ld[which(ld$RS_Number == "."), "Coord"]

ld <- ld %>% separate(Coord, c("CHR_ID", "CHR_POS"), sep = ":") %>%
  mutate(CHR_ID = gsub(pattern = ".*chr", replacement = "", x = CHR_ID)) %>%
  dplyr::rename("SNPS"="RS_Number", "LEAD_SNP"="snp")


